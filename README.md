# Bioutils Collection

Bioinformatics utilities for sequence analysis, alignment, annotation, FASTA
processing, and molecular biology workflows.

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
```

See [the major check plan](MAJOR_CHECK_PLAN.md) for the staged correctness,
packaging, file safety, and performance audit.

## Contributing

Fork the repository, create a branch, add regression tests for fixes, and run
relevant checks before submitting a pull request. Include type hints and
regression tests, and update documentation when behavior changes.

## Contact and related projects

- Author: [Mykyta Forofontov](https://github.com/MForofontov)
- [Report bugs or request features](https://github.com/MForofontov/bioutils-collection/issues)
- [PyPI package](https://pypi.org/project/bioutils-collection/)
- [Biopython](https://biopython.org/) and [scikit-bio](https://scikit-bio.org/)
- [pyutils-collection](https://github.com/MForofontov/pyutils-collection)

## License

MIT © 2026 Mykyta Forofontov. See [LICENSE](LICENSE).
