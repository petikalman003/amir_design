"""Regression checks for the published default and explicit assembly fixtures."""

import json

import pytest

import amir_design.scaffold as scaffold_module

from amir_design.scaffold import (
    MireTemplate,
    DEFAULT_MIRE_TEMPLATE,
    ScaffoldUnavailableError,
    build_mire_insert,
    load_mire_template,
)


GUIDE = "ATAGTTTCAAACATCATCTTGT"
TARGET = "ACAAGATGATGTTTGAAACTAT"


def fixture_template(**overrides):
    values = {
        "name": "explicit assembly test fixture, not a biological template",
        "prefix": "AACCGG",
        "loop": "TTAA",
        "suffix": "CCGGTT",
        "source": "unit-test assembly fixture",
        "passenger_rule": "reverse_complement",
    }
    values.update(overrides)
    return MireTemplate(**values)


def test_missing_template_still_fails_clearly(monkeypatch):
    monkeypatch.setattr(scaffold_module, "DEFAULT_MIRE_TEMPLATE", None)
    with pytest.raises(ScaffoldUnavailableError, match="complete miR-E synthesis template"):
        build_mire_insert(GUIDE)


def test_published_default_matches_source_vector_fragment():
    source_guide = "TAGATAAGCATTATAATTCCTA"
    source_fragment = (
        "CTCGAGAAGGTATATTGCTGTTGACAGTGAGCGCAGGAATTATAATGCTTATCTATAGTGAAGCCACAGATGTA"
        "TAGATAAGCATTATAATTCCTATGCCTACTGCCTCGGACTTCAAGGGGCTAGAATTC"
    )
    insert = build_mire_insert(source_guide)
    assert len(source_fragment) == 131
    assert len(insert) == 119
    assert "CTCGAG" + insert + "GAATTC" == source_fragment
    assert DEFAULT_MIRE_TEMPLATE.suffix[15:20] == "ACTTC"
    assert DEFAULT_MIRE_TEMPLATE.provenance["source_table"].endswith("cell C3")
    assert DEFAULT_MIRE_TEMPLATE.provenance["reference_fragment_1based_inclusive"] == [4424, 4542]


def test_published_default_preserves_mature_guide_for_another_candidate():
    insert = build_mire_insert(GUIDE)
    assert len(insert) == 119
    guide_start = len(DEFAULT_MIRE_TEMPLATE.prefix) + 22 + len(DEFAULT_MIRE_TEMPLATE.loop)
    assert insert[guide_start:guide_start + 22] == GUIDE


def test_explicit_template_preserves_guide_and_orientation():
    template = fixture_template()
    insert = build_mire_insert(GUIDE, template=template)
    assert insert == "AACCGG" + TARGET + "TTAA" + GUIDE + "CCGGTT"
    assert len(insert) == 60


def test_declared_passenger_mismatch_leaves_guide_unchanged():
    template = fixture_template(passenger_rule="shrnai_first_base_mismatch")
    insert = build_mire_insert(GUIDE, template=template)
    assert insert == "AACCGG" + "C" + TARGET[1:] + "TTAA" + GUIDE + "CCGGTT"
    cg_terminal = GUIDE[:-1] + "G"
    alternate = build_mire_insert(cg_terminal, template=template)
    assert alternate[6] == "A"
    assert alternate[-28:-6] == cg_terminal


def test_fixed_first_base_is_explicit():
    template = fixture_template(passenger_rule="fixed_first_base", passenger_first_base="G")
    assert build_mire_insert(GUIDE, template=template)[6:28] == "G" + TARGET[1:]


def test_template_json_round_trip_and_unknown_fields(tmp_path):
    template = fixture_template()
    path = tmp_path / "template.json"
    path.write_text(json.dumps(template.to_dict()))
    assert load_mire_template(path) == template
    invalid = template.to_dict() | {"shell_command": "ignored content must not be executed"}
    path.write_text(json.dumps(invalid))
    with pytest.raises(ValueError, match="Unknown"):
        load_mire_template(path)


@pytest.mark.parametrize("guide", ["A" * 21, "N" * 22, "a" * 22, "U" * 22])
def test_invalid_guide_rejected(guide):
    with pytest.raises(ValueError):
        build_mire_insert(guide, template=fixture_template())


def test_ambiguous_or_missing_template_information_rejected():
    with pytest.raises(ValueError, match="uppercase"):
        fixture_template(prefix="NNNN")
    with pytest.raises(ValueError, match="source"):
        fixture_template(source="")
    with pytest.raises(ValueError, match="passenger_first_base"):
        fixture_template(passenger_rule="fixed_first_base")
