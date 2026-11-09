# Bioutils Collection

[![PyPI version](https://img.shields.io/pypi/v/bioutils-collection.svg)](https://pypi.org/project/bioutils-collection/)
[![Python versions](https://img.shields.io/pypi/pyversions/bioutils-collection.svg)](https://pypi.org/project/bioutils-collection/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Bioinformatics utilities** — 80+ functions across 12 specialized modules for
> sequence analysis, alignment, annotation, FASTA processing, and molecular biology.

Current version: **0.3.3**. See the [release notes](CHANGELOG.md#033---2026-11-08)
for correctness fixes, documentation updates, and validation details.
The PyPI badges above report the version and Python support published on PyPI.

## Installation

```bash
pip install bioutils-collection
```

Requires Python 3.10+ and NumPy, SciPy, scikit-learn, and Numba. Runtime dependencies
are declared in `pyproject.toml`. To install a local checkout:

```bash
pip install -r requirements.txt
```

## Quick start

```python
from bioutils_collection import (
    gc_content,
    needleman_wunsch,
    parse_fasta,
    reverse_complement,
    translate_dna_to_protein,
)

assert reverse_complement("ATCGATCG") == "CGATCGAT"
assert gc_content("ATCGATCG") == 50.0  # Percentage, from 0 to 100

# Alignment returns the score first, followed by the aligned sequences.
score, aligned1, aligned2 = needleman_wunsch("GATTACA", "GCATGCU")
assert translate_dna_to_protein("ATGGCCTAA") == "MA*"

# parse_fasta accepts FASTA text.
for header, sequence in parse_fasta(">example\nATGGCCTAA\n"):
    print(f"{header}: {len(sequence)} bp")
```

## Modules

- `alignment_functions`: global/local alignment, distances, identity, score ratios.
- `sequence_operations`: complements, ORFs, CpG islands, k-mers, minimizers, syncmers.
- `translation_functions`: transcription, genetic codes, batch/streaming translation.
- `gc_functions`: GC percentages, GC skew, and windowed GC analysis.
- `sequence_statistics`: codon usage, entropy, conservation, protein properties.
- `data_validation`: DNA and protein validation.
- `fasta_misc`: FASTA parsing, writing, filtering, splitting, and format conversion.
- `annotation_functions`: BED/GFF/GTF/VCF parsing, conversion, and statistics.
- `motif_functions`: motif search, consensus generation, and pattern matching.
- `repeat_functions`: tandem repeats and palindrome detection.
- `restriction_functions`: restriction enzyme site identification.
- `clustering_functions`: motif clustering.

Functions are available from their modules; common utilities are also exported
from `bioutils_collection`. Docstrings describe inputs, outputs, and validation.

BED parsing preserves zero-based, half-open coordinates and skips track/browser
metadata. Records on chromosomes named `track` or `browser` are retained when
their tab-delimited coordinate fields are integers. BLAST score ratios require
finite scores: the reference must be positive and the target non-negative.

## Sequence analysis example

```python
from bioutils_collection import find_orfs, gc_content_windows
from bioutils_collection.sequence_operations.find_cpg_islands import find_cpg_islands
from bioutils_collection.sequence_statistics.amino_acid_composition import (
    amino_acid_composition,
)
from bioutils_collection.sequence_statistics.calculate_isoelectric_point import (
    calculate_isoelectric_point,
)

dna_sequence = "ATG" + "GCC" * 100 + "TAA"
orfs = [orf for orf in find_orfs(dna_sequence) if len(orf[2]) >= 300]
gc_profile = gc_content_windows(dna_sequence, window_size=100, step_size=25)
cpg_islands = find_cpg_islands(dna_sequence)
protein_seq = "MAAAAA"
pi = calculate_isoelectric_point(protein_seq)
composition = amino_acid_composition(protein_seq)
```

GC windows return `(start, end, percentage)` tuples with zero-based, exclusive
end positions. Incomplete trailing windows are skipped.

## Development and testing

```bash
git clone https://github.com/MForofontov/bioutils-collection.git
cd bioutils-collection
python -m venv .venv
# Activate .venv using the command for your shell, then:
python -m pip install -e ".[dev]"
python -m pytest
python -m pytest -m alignment
python -m pytest --cov=bioutils_collection --cov-report=html
python -m mypy bioutils_collection
python -m ruff check .
python -m ruff format --check .
python scripts/check_documentation.py --report dist/documentation-report.json
```

See [the major check plan](MAJOR_CHECK_PLAN.md) for the staged correctness,
packaging, file safety, and performance audit.

To reproduce performance measurements separately from functional tests:

```bash
python benchmarks/major_audit.py --output benchmarks/major_audit_results.json
python benchmarks/paired_audit.py --output benchmarks/paired_audit_results.json
python benchmarks/alignment_recheck.py --output benchmarks/alignment_recheck_results.json
python scripts/check_distribution.py --artifact-dir dist --report dist/distribution-report.json
```

The benchmark uses deterministic inputs, checks reference equivalence, and records
repeated timings and peak Python allocations. Native allocations and worker memory
are not included in that memory measurement. Timings depend on the machine.
The paired runner additionally records immutable baseline/current source manifests,
five counterbalanced timing samples, and separately measured process peak RSS.
It labels corrected-reference comparisons and does not claim speedups when old
outputs differ. Runtime/memory retention criteria are in [AGENTS.md](AGENTS.md).

The alignment recheck records actual square/unequal sequence dimensions and uses
calibrated batches of at least 100 ms across five counterbalanced rounds. It
includes an unchanged edit-distance control; drift above 10% makes timing evidence
inconclusive. Repeat once without concurrent heavy checks and retain both reports.
Persistent control drift does not support a new speed claim.

The documentation checker parses every package docstring before executing the
examples in a separate process with a temporary working directory and a timeout.
File examples create temporary inputs, stochastic examples assert invariants,
and large-job examples use small inputs. The unit suite runs this check as well.

Use a fresh artifact directory for distribution checks. The verifier derives the
version from package metadata and installs the wheel and sdist into separate
temporary virtual environments, checks imports/examples outside the checkout,
and rejects stale package sources. Default checks resolve dependencies freshly.
For offline local checks, `--reuse-dependencies --no-build-isolation` can reuse
available dependencies/build tools; the report explicitly records this limitation.

## Conventions and scientific models

- BED parsers retain zero-based, half-open coordinates; GFF/GTF parsers retain
  one-based, inclusive coordinates. `annotation_to_bed` expects one-based,
  inclusive input, including its legacy `chrom` alias. Convert a parsed BED start
  by adding one before passing it to this converter or `annotation_statistics`.
  Annotation statistics sum record lengths, including overlaps repeatedly.
- GFF/GTF writers accept their corresponding parser records as well as existing
  writer field names. Existing writer fields take precedence. GFF attribute
  dictionaries retain percent-encoded values; writers do not decode them.
- Closed syncmers select a k-mer when its leftmost minimum s-mer is first or last.
  Open syncmers retain offset zero. Both pattern-matching modes include overlaps.
- Translation uses the NCBI internal-codon AAs rows. Alternative initiation and
  context-dependent termination are not modeled; table 32 is not implemented.
  Empty records are skipped by file translation. Dictionary translation keeps the
  last duplicate ID; streaming translation preserves each record. Output replaces
  its destination only after successful processing; input aliases are rejected.
  Streaming memory grows with the largest record, rather than total file size.
- Protein pI solves zero net charge using the fixed pKa values in its docstring.
  Approximate mass tables are rounded to 0.1 Da and include linear-chain bond and
  terminal corrections. DNA/RNA masses assume a single strand, 5-prime phosphate,
  and 3-prime hydroxyl. These models omit modifications and sequence context.
- ENC uses population homozygosity without a finite-sample correction; sparse
  sequences are unreliable unbiased estimates. Sequence complexity reports a
  unique-substring fraction. Melting temperature uses composition approximations,
  rather than nearest-neighbor thermodynamics. DNA validation computes GC over
  standard A/T/G/C bases only; ambiguous symbols are excluded from its denominator.
- FASTA statistics return the mean of per-record GC fractions (0 to 1), count
  G/C case-insensitively, and include ambiguous characters in each denominator.
  Empty records contribute zero. GenBank conversion accepts LF/CRLF/CR line
  endings and removes whitespace between sequence groups.
- Restriction sites must be nonempty strings; matches are literal,
  case-insensitive, overlapping, and reported at zero-based start positions.
- Low-complexity masking requires a finite threshold in [0, 2] and masks the
  union of windows strictly below it. N counts as a fifth entropy symbol;
  repeats above the entropy threshold remain unchanged.

## Contributing

Fork the repository, create a branch, add regression tests for fixes, and run
relevant checks before submitting a pull request. Include type hints and
regression tests, and update documentation when behavior changes.
Follow [AGENTS.md](AGENTS.md): always update documentation to reflect the change,
add release notes, and update the package version once per completed changeset.
Fixes, documentation, and tooling use a patch increment unless a larger change
requires a minor or breaking-change increment. Version updates do not publish.

## Contact and related projects

- Author: [Mykyta Forofontov](https://github.com/MForofontov)
- [Report bugs or request features](https://github.com/MForofontov/bioutils-collection/issues)
- [PyPI package](https://pypi.org/project/bioutils-collection/)
- [Biopython](https://biopython.org/) and [scikit-bio](https://scikit-bio.org/)
- [pyutils-collection](https://github.com/MForofontov/pyutils-collection)

## License

MIT © 2026 Mykyta Forofontov. See [LICENSE](LICENSE).
