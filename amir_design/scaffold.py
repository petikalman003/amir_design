"""A published miR-E insert and explicit alternative vector templates.

The default is the 119-nt miR-E fragment between the XhoI/EcoRI sites in the
published pCF806-shRen.713 vector: Amen et al. (2022), Table S1, Vectors C3.
It preserves the optimized downstream motif and the source's fixed flanks.
Expression-vector sequence and user-specific cloning junctions are excluded.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


MIRE_PUBLICATION = "https://doi.org/10.1016/j.celrep.2013.11.020"
MIRE_TEMPLATE_PUBLICATION = "https://doi.org/10.1016/j.crmeth.2022.100239"
MIRE_TEMPLATE_SOURCE = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9308131/supplementaryFiles"
MIRE_TEMPLATE_WORKBOOK_SHA256 = "e8a0906bad08e713c052d085976514adc15c793a2f81ccb7875915af4c6887f0"
MIRE_PROVENANCE = {
    "backbone": "miR-E (optimized human miR-30a)",
    "publication": MIRE_PUBLICATION,
    "template_publication": MIRE_TEMPLATE_PUBLICATION,
    "template_source": MIRE_TEMPLATE_SOURCE,
    "source_table": "Table S1 (mmc2.xlsx), Vectors sheet, cell C3",
    "source_workbook_sha256": MIRE_TEMPLATE_WORKBOOK_SHA256,
    "reference_vector": "pCF806-shRen.713 (UT4GEPIR); Addgene #186712",
    "reference_vector_length_nt": 10109,
    "reference_fragment_1based_inclusive": [4424, 4542],
    "default_template_status": "verified published 119-nt miR-E fragment",
    "insert_length_nt": 119,
    "boundaries": "between complete XhoI/EcoRI sites; both restriction sites excluded",
    "downstream_processing_motif": "ACTTC (DCNNC), suffix positions 16-20",
    "cloning_junctions": "not included; vector-specific",
}


class ScaffoldUnavailableError(ValueError):
    """A complete miR-E template is required before exporting an insert."""


def _validate_dna(sequence: str, label: str, *, expected_length: int | None = None) -> None:
    if not isinstance(sequence, str) or not sequence or set(sequence) - set("ACGT"):
        raise ValueError(f"{label} must contain uppercase A/C/G/T only")
    if expected_length is not None and len(sequence) != expected_length:
        raise ValueError(f"{label} must be {expected_length} nt")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@dataclass(frozen=True)
class MireTemplate:
    """A declared complete insert around two variable 22-nt stem arms.

    ``prefix`` and ``suffix`` must include all fixed hairpin context that will
    be synthesized. They do not automatically include vector cloning junctions.
    ``passenger_rule`` is explicit because a passenger mismatch changes the
    insert; neither the guide nor its strand identity is changed.
    """

    name: str
    prefix: str
    loop: str
    suffix: str
    source: str
    passenger_rule: str
    passenger_first_base: str | None = None

    def __post_init__(self) -> None:
        if not self.name or not self.source:
            raise ValueError("A synthesis template needs a name and source specification")
        for label in ("prefix", "loop", "suffix"):
            _validate_dna(getattr(self, label), label)
        allowed_rules = {"reverse_complement", "shrnai_first_base_mismatch", "fixed_first_base"}
        if self.passenger_rule not in allowed_rules:
            raise ValueError(f"passenger_rule must be one of {sorted(allowed_rules)}")
        if self.passenger_rule == "fixed_first_base":
            _validate_dna(self.passenger_first_base, "passenger_first_base", expected_length=1)
        elif self.passenger_first_base is not None:
            raise ValueError("passenger_first_base is only used by fixed_first_base")

    @property
    def provenance(self) -> dict[str, object]:
        provenance = {
            "backbone": self.name,
            "source": self.source,
            "template_origin": "explicitly supplied full template",
            "passenger_rule": self.passenger_rule,
            "passenger_first_base": self.passenger_first_base,
            "template_sha256": hashlib.sha256(
                json.dumps(self.to_dict(), sort_keys=True).encode("utf-8")
            ).hexdigest(),
            "cloning_junctions": "include only junctions explicitly supplied in prefix/suffix",
        }
        if self == DEFAULT_MIRE_TEMPLATE:
            provenance.update(MIRE_PROVENANCE)
            provenance["template_origin"] = "verified published miR-E fragment"
        return provenance

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


# Source vector (Table S1, Vectors C3) has the following exact 131-nt fragment:
# XhoI + prefix27 + passenger22 + loop19 + guide22 + suffix29 + EcoRI.
# Preserve all fixed inter-site bases, removing only both 6-nt enzyme sites.
# The published guide TAGATAAGCATTATAATTCCTA and passenger
# CAGGAATTATAATGCTTATCTA validate the basal mismatch rule. The same rule was
# cross-checked against all 35 single-hairpin 97-mer examples in Table S1.
DEFAULT_MIRE_TEMPLATE: MireTemplate | None = MireTemplate(
    name="published miR-E fragment (119 nt; restriction sites excluded)",
    prefix="AAGGTATATTGCTGTTGACAGTGAGCG",
    loop="TAGTGAAGCCACAGATGTA",
    suffix="TGCCTACTGCCTCGGACTTCAAGGGGCTA",
    source=f"{MIRE_TEMPLATE_PUBLICATION}; Table S1 (mmc2.xlsx), Vectors C3, pCF806-shRen.713",
    passenger_rule="shrnai_first_base_mismatch",
)


def load_mire_template(path: Path) -> MireTemplate:
    """Read an explicitly supplied full-template JSON without executing content."""
    with Path(path).open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("miR-E template JSON must be an object")
    fields = {"name", "prefix", "loop", "suffix", "source", "passenger_rule", "passenger_first_base"}
    unknown = set(data) - fields
    if unknown:
        raise ValueError(f"Unknown miR-E template fields: {sorted(unknown)}")
    try:
        return MireTemplate(**data)
    except TypeError as exc:
        raise ValueError("miR-E template JSON is missing required fields") from exc


def load_template(path: Path) -> MireTemplate:
    """Load and structurally validate an explicitly supplied miR-E template."""
    return load_mire_template(path)


def build_mire_insert(guide_dna: str, *, template: MireTemplate | None = None) -> str:
    """Build the declared insert, preserving the mature guide exactly.

    The default preserves the published miR-E inter-site fragment. Alternative
    vector templates must be explicit; no vector-specific junctions are added.
    """
    _validate_dna(guide_dna, "Mature guide", expected_length=22)
    selected = template if template is not None else DEFAULT_MIRE_TEMPLATE
    if selected is None:
        raise ScaffoldUnavailableError(
            "A complete miR-E synthesis template has not been verified. Supply a full "
            "template JSON from the intended vector/source; Supplementary Figure S1 "
            "and shRNAI's scoring scaffold do not establish the complete insert."
        )
    passenger = _reverse_complement(guide_dna)
    if selected.passenger_rule == "shrnai_first_base_mismatch":
        first = "C" if guide_dna[-1] in "AT" else "A"
        passenger = first + passenger[1:]
    elif selected.passenger_rule == "fixed_first_base":
        passenger = selected.passenger_first_base + passenger[1:]
    return selected.prefix + passenger + selected.loop + guide_dna + selected.suffix
