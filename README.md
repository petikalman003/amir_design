# amiR Design

A Streamlit workflow connecting the uploaded Design.py selection rules,
[shRNAI](https://github.com/ParkSJ-91/shRNAI) activity prediction, and
[SeedMatchR](https://github.com/tacazares/SeedMatchR) full-transcript seed screening.
Submit one CDS, select Human or Mouse, and download one Excel workbook containing
candidate guides, predicted scores, off-target sites, target isoform coverage, and
the reference/model provenance.

## Run

From this repository, on the cloud machine's Debian 13 environment:

```bash
bash scripts/setup.sh
export PATH=/workspace/.runtime/bin:$PATH
.venv/bin/python -m amir_design.reference_cli --species both --root data/references
bash scripts/start.sh
```

`setup.sh` installs Python 3.11.12, frozen Python dependencies, R/Biostrings, and
BLAST+ into workspace-local directories. It fetches the upstream tools at pinned
commits, checks the shRNAI source/model hashes, and preserves existing dependency
checkouts rather than resetting edits. `bash scripts/setup.sh --references` also
prepares both official references. A verified published miR-E insert template is
included in the application.

The cloud task already has an isolated checkout. Use it directly; do not create
a worktree for routine development. Dependencies and reference data are ignored
local files; processes must be restarted in a new task.

Runtime installation and both official reference downloads are validated in the
current machine. The required download hostnames are saved in the environment
draft. Review/save the configuration and publish the prepared environment to
retain its dependencies and references for later tasks; live server processes
must still be restarted.

## References and isoforms

The defaults are GENCODE human **v47 / GRCh38** and mouse **M36 / GRCm39** full
transcript FASTAs. Preparation verifies the archive against the publisher's
MD5SUMS, records SHA-256 hashes, and atomically builds a normalized FASTA, SQLite
gene/transcript mapping, and BLAST index. Both downstream tools use the same
reference. No canonical, protein-coding, or longest-transcript filter is applied.
Noncoding, NMD, and retained-intron transcript records remain included.

CDS mapping requires an exact full-length, plus-strand match to one gene.
Ambiguous gene assignments and unmatched sequences stop with an explicit error.
Every retained candidate must have its complete 22-base target site in **every
annotated transcript of that gene**. If no candidate qualifies, the report says
`no_eligible_candidate`; no fallback guide set is generated.

Custom full-transcript references are supported:

```bash
export PATH=/workspace/.runtime/bin:$PATH
.venv/bin/python -m amir_design.reference_cli --species human \
  --root /path/to/custom-references --fasta /path/to/transcripts.fa.gz \
  --release my-release --assembly my-assembly
```

Headers must be official GENCODE headers or `gene_id|transcript_id|gene_symbol`.
An optional `--sha256` checks a source-provided checksum. Existing bundles are
verified and reused; a changed source requires a different root.

## Sequence and score definitions

- Design outputs are 22-base **target-sense** DNA sites. Mature guides are their
  exact reverse complements, written 5′→3′. Their RNA sequence is exported too.
- The original GC, homopolymer, and position filters are retained. The original
  target-position-3 exclusion differs from a literal guide-strand reading of
  Supplementary Figure S1. The form exposes both profiles and records the
  selected rule. Default: uploaded Design.py rules (target position 3 is not A).
- Every eligible unique sequence is scored; there is no early top-200/top-five
  cutoff. Duplicate sequences retain their CDS positions.
- The real shRNAI models run on CPU using the upstream encoder and unchanged
  weights. Legacy integral float pooling settings are normalized in memory.
  The raw score is exploratory activity ranking, **not a knockdown percentage**.
  It uses shRNAI's supplied scaffold; transfer to a supplied miR-E construct is
  not experimentally validated by this application.
- The actual pinned SeedMatchR `get_seed()` and `SeedMatchR()` sequence functions
  are sourced in R. Unrelated RNA-seq report/annotation dependencies are not
  installed. No DESeq2 or experimental expression data is required.
- Exact `mer8`, `mer7m8`, and `mer7A1` seed matches are searched across complete
  spliced transcripts. Mismatches, indels, and wobble matches are disabled.
  Six-base searches are supported by the adapter but excluded from the default
  form/report; unsearched classes are not reported as zero hits.
- All transcripts of the intended gene are excluded from off-target summaries.
  Gene counts, transcript counts, and transcript-site counts are separate.
  Overlapping classes at one binding location keep the strongest label:
  8mer > 7mer-m8 > 7mer-A1. Coordinates are 1-based inclusive transcript positions,
  not genomic positions. Sites indicate potential binding, not demonstrated
  off-target silencing.

The report sorts by predicted activity descending, then off-target gene count
ascending. These are separate metrics, with no arbitrary combined score.

## miR-E synthesis output

The [Kopasz paper](https://doi.org/10.1186/s12915-022-01262-x) specifies **miR-E**,
the [Fellmann optimized miR-30a backbone](https://doi.org/10.1016/j.celrep.2013.11.020).
Its Figure S1 explicitly omits complete amiR structures. The application's
default is therefore verified against author-supplied vector sequences in
[Amen et al. (2022), Table S1](https://doi.org/10.1016/j.crmeth.2022.100239),
Vectors sheet cell C3, pCF806-shRen.713 / Addgene #186712. The official EBI
supplement was checked and its SHA-256 recorded in scaffold provenance.

The default exports the **119-nt miR-E fragment between the XhoI/EcoRI sites**,
preserving the fixed source flanks and optimized downstream processing motif.
Both restriction sites are excluded. Source-fragment reconstruction and mature
guide preservation are regression-tested. This is a synthesizable hairpin insert;
it does not include the complete expression vector or your vector-specific cloning
junctions. The model score still comes from shRNAI's supplied 113-nt scoring scaffold.

For a different vector insertion format, supply an explicit verified template
JSON and set `AMIR_MIRE_TEMPLATE` to its path:

| Field | Required content |
| --- | --- |
| `name` | Exact backbone/vector identification |
| `prefix` | Verified fixed DNA before the 22-base passenger strand |
| `loop` | Verified intervening loop DNA |
| `suffix` | Verified fixed DNA after the 22-base mature guide |
| `source` | Publication, vector accession, or laboratory sequence specification |
| `passenger_rule` | `reverse_complement`, `shrnai_first_base_mismatch`, or `fixed_first_base` |
| `passenger_first_base` | Required single DNA base only for `fixed_first_base` |

All sequence fields require uppercase A/C/G/T. Assembly preserves the mature
guide exactly. Cloning junctions appear only if explicitly included in the
template's prefix/suffix; otherwise the output is an unflanked insert rather
than vector-ready cloning oligos. A declared source and template hash are
recorded, but structural validation cannot independently establish that a
user-supplied sequence is the intended biological backbone. If a default template
is deliberately disabled, analysis and mature guide exports remain available,
with empty insert fields and explicit `template_required` status.

## Configuration and validation

Optional variables: `AMIR_REFERENCE_ROOT` (parent containing human/mouse bundles),
`AMIR_TOOLS_ROOT`, `AMIR_RSCRIPT`, `AMIR_MIRE_TEMPLATE`, `AMIR_RUNTIME_DIR`, and
`AMIR_PORT` (default 8501). No API credentials are required.

```bash
PATH=/workspace/.runtime/bin:$PATH TF_CPP_MIN_LOG_LEVEL=2 \
  TF_NUM_INTRAOP_THREADS=2 TF_NUM_INTEROP_THREADS=1 \
  .venv/bin/python -m pytest -q
.venv/bin/ruff check .
```

Tests execute real BLAST, real shRNAI model inference, real R/SeedMatchR searches,
and a Streamlit submission/download using explicitly synthetic transcript
fixtures. Fixtures are not installed as human/mouse production references.
Scoring regression tests reproduce the upstream notebook's example scores.
Missing-template and no-shared-isoform outcomes are tested separately. Real CDS
examples are in `examples/`: human LCE3A and mouse Krtap2-21, extracted using CDS
coordinates in the checksum-verified official GENCODE protein-coding transcript
archives. Both examples passed scoring and full-reference screening. Their short
seed matches are plentiful; these examples demonstrate operation rather than
endorsing either guide as a suitable experimental choice.

Large screens can produce millions of transcript sites. Workbook generation
streams rows and splits sheets at Excel's row limit. The R/Python screening
stage still materializes hit tables. Large candidate batches need sufficient RAM;
reference preparation and representative screens do not benchmark every gene or
the largest possible candidate set.
