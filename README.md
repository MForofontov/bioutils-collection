# 🧬 Bioutils Collection

[![PyPI version](https://badge.fury.io/py/bioutils-collection.svg)](https://pypi.org/project/bioutils-collection/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Production-ready bioinformatics toolkit** - 77+ optimized functions for sequence analysis, alignment, annotation, and molecular biology workflows.

## ✨ Highlights

- 🚀 **77+ specialized functions** across 13 bioinformatics domains
- 🔒 **Fully typed** with complete type hints (mypy strict)
- 📊 **Research-grade algorithms** - Needleman-Wunsch, Smith-Waterman, and more
- ⚡ **Performance optimized** for large-scale genomic data
- ✅ **Extensively tested** with comprehensive test coverage
- 📝 **Self-documenting** - NumPy-style docstrings with examples

## 📦 Installation

```bash
pip install bioutils-collection
```

**Requirements:** Python 3.10+ with numpy, scipy, and scikit-learn

## 🎯 Quick Start

```python
from bioutils_collection import (
    reverse_complement,
    gc_content,
    needleman_wunsch,
    parse_fasta,
    translate_dna_to_protein,
)

# Sequence manipulation
seq = "ATCGATCG"
rev_comp = reverse_complement(seq)  # "CGATCGAT"

# Calculate GC content
gc = gc_content("ATCGATCG")  # 0.5

# Global sequence alignment
seq1, seq2 = "GATTACA", "GCATGCU"
aligned1, aligned2, score = needleman_wunsch(seq1, seq2)

# Parse FASTA files
for header, sequence in parse_fasta("genome.fasta"):
    print(f"{header}: {len(sequence)} bp")

# Translate DNA to protein
protein = translate_dna_to_protein("ATGGCCTAA")  # "MA*"
```

## 🧬 Modules

### Core Sequence Operations
- **`alignment_functions`** - Pairwise & multiple sequence alignment (Needleman-Wunsch, Smith-Waterman, BLAST score ratio)
- **`sequence_operations`** - Reverse complement, ORF finding, CpG islands, low-complexity filtering
- **`translation_functions`** - DNA↔RNA transcription, translation with custom codon tables

### Sequence Analysis & Statistics
- **`gc_functions`** - GC content, GC skew, windowed GC profiling
- **`sequence_statistics`** - Codon usage (CAI, ENC, RSCU), melting temp, isoelectric point, amino acid composition
- **`data_validation`** - DNA/RNA/protein sequence validation

### File I/O & Parsing
- **`fasta_misc`** - FASTA parsing, writing, filtering, splitting, concatenation, primer generation
- **`annotation_functions`** - BED/GFF/GTF/VCF parsing and conversion, annotation statistics

### Pattern & Motif Discovery
- **`motif_functions`** - Motif search, consensus generation, pattern matching
- **`repeat_functions`** - Tandem repeat finder, palindrome detection
- **`restriction_functions`** - Restriction enzyme site identification
- **`clustering_functions`** - Motif clustering and grouping

## 🔬 Use Cases

**Genomic Analysis**
```python
from bioutils_collection import find_orfs, gc_content_windows, find_cpg_islands

# Find all ORFs in a sequence
orfs = find_orfs(dna_sequence, min_length=300)

# Sliding window GC analysis
gc_proDevelopment

```bash
# Clone repository
git clone https://github.com/MForofontov/bioutils-collection.git
cd bioutils-collection

# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run specific test categories
pytest -m alignment
pytest -m fasta
pytest -m translation

# Type checking
mypy bioutils_collection

# Linting
ruff check .

# Coverage report
pytest --cov=bioutils_collection --cov-report=html
```

## 📚 API Documentation

All functions include:
- ✅ **Complete type hints** for static analysis
- 📖 **NumPy-style docstrings** with parameter descriptions
- 💡 **Usage examples** in docstrings
- ⚠️ **Complexity notes** for performance-critical code
- 📎 **Algorithm references** where applicable

Example:
```python
from bioutils_collection import needleman_wunsch
help(needleman_wunsch)  # Comprehensive documentation
```

## 🤝 Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Add comprehensive tests
4. Ensure all tests pass (`pytest`)
5. Add type hints and docstrings
6. Submit a pull request

**Development Guidelines:**
- Follow existing code style (ruff formatting)
- Add tests for all new functions
- Update documentation
- Keep functions focused and single-purpose

## 🔗 Related Projects

- **[BioPython](https://biopython.org/)** - Comprehensive biological computation library
- **[scikit-bio](http://scikit-bio.org/)** - Scientific Python library for bioinformatics
- **[pyutils-collection](https://github.com/MForofontov/pyutils-collection)** - General Python utilities (sister project)

# Protein properties
pi = calculate_isoelectric_point(protein_seq)
composition = amino_acid_composition(protein_seq)

# Primer design
tm = melting_temperature("ATCGATCGATCG")
```

## 🧪 Testing

```bash
# Run all tests
pytest

# Run specific module tests
pytest -m alignment
pytest -m fasta
pytest -m translation

# Run with coverage
pytest --cov=bioutils_collection --cov-report=html
```

## 📖 Documentation

Each function includes:
- **� License

MIT License - see [LICENSE](LICENSE) for details.

## 📊 Project Stats

- **77+ Functions** across 13 specialized modules
- **670+ Tests** with comprehensive coverage
- **Type-safe** with mypy strict mode
- **Python 3.10+** with modern type hints

## 📮 Contact & Support

- **Author:** Mykyta Forofontov
- **Repository:** [github.com/MForofontov/bioutils-collection](https://github.com/MForofontov/bioutils-collection)
- **Issues:** [Report bugs or request features](https://github.com/MForofontov/bioutils-collection/issues)
- **PyPI:** [pypi.org/project/bioutils-collection](https://pypi.org/project/bioutils-collection/)

---

Star this repo if you find it useful!

## Changelog

### v0.2.0 (January 28, 2026)

#### Removed - 5 Redundant Functions

**Alignment Functions:**
- **phylogenetic_distance** - Redundant with `hamming_distance(s1, s2) / len(s1)`
  - Just a simple division wrapper, no added value

**Sequence Statistics:**
- **dinucleotide_frequency** - Redundant with `kmer_frequency(seq, 2)`
  - Exact same functionality, just calls kmer_frequency internally
- **nucleotide_frequency** - Redundant with `kmer_frequency(seq, 1)`
  - Functionally identical output to kmer_frequency with k=1

**Sequence Operations:**
- **sequence_masking** - Incorrect implementation
  - Only masked simple runs (AAAA→NNNN), didn't detect true low-complexity
  - Users should use `remove_low_complexity_regions` (proper entropy-based)
- **reverse_sequence** - Trivial one-liner
  - Just `seq[::-1]`, no need for dedicated function

#### Fixed - 2 Functions

**blast_score_ratio:**
- Added validation for negative target_score values
- Added ValueError in docstring Raises section

**restriction_site_finder:**
- Added type checking for sequence (must be str)
- Added type checking for sites (must be Sequence)
- Added empty validation for both parameters
- Updated docstring with TypeError and ValueError details

#### Verified
- All remaining ~70 functions verified as unique and non-redundant
- 692 tests passing, 6 skipped
- No broken dependencies or import errors

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🔗 Related Projects

- [BioPython](https://biopython.org/) - Comprehensive bioinformatics toolkit
- [scikit-bio](http://scikit-bio.org/) - Scientific Python library for bioinformatics

## 📮 Contact

**Author:** Mykyta Forofontov  
**Repository:** https://github.com/MForofontov/bioutils-collection  
**Issues:** https://github.com/MForofontov/bioutils-collection/issues
