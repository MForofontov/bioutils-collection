# Bioutils Collection

> **Comprehensive bioinformatics utilities** - 77+ production-ready functions for sequence analysis, alignment, annotation, and molecular biology workflows.

## 🎯 What is This?

A curated collection of **77+ bioinformatics functions** across **13 specialized modules** - designed for genomics research, sequence analysis, and molecular biology applications. Each function is self-contained with type hints, comprehensive docstrings, and handles edge cases gracefully.

**Philosophy:**
- 🧬 **Research-grade** - Production-ready bioinformatics algorithms
- 🔒 **Type-safe** - Complete type hints (Python 3.10+)
- 📝 **Self-documenting** - NumPy-style docstrings with examples
- ✅ **Well-tested** - Comprehensive test coverage
- ⚡ **Efficient** - Optimized for large-scale sequence analysis

## 📦 Quick Start

```bash
# Install from PyPI
pip install bioutils-collection

# Or install for development
git clone https://github.com/MForofontov/bioutils-collection.git
cd bioutils-collection
pip install -e ".[dev]"
```

## 🧬 Modules Overview

| Module | Functions | Purpose |
|--------|-----------|---------|
| 🧬 **alignment_functions** | Sequence alignment algorithms (Needleman-Wunsch, Smith-Waterman) |
| 📝 **annotation_functions** | Genomic annotation parsing and statistics (BED, VCF, GFF) |
| 🔬 **clustering_functions** | Sequence clustering and similarity analysis |
| 🔬 **data_validation** | Sequence validation (DNA, RNA, protein) |
| 📄 **fasta_misc** | FASTA file operations (parsing, writing, filtering) |
| 🧮 **gc_functions** | GC content analysis and skew calculations |
| 🔍 **motif_functions** | Motif finding and pattern matching |
| 🔁 **repeat_functions** | Repeat detection and analysis |
| ✂️ **restriction_functions** | Restriction enzyme site identification |
| 🧬 **sequence_operations** | Sequence manipulation (reverse complement, shuffling, filtering) |
| 📊 **sequence_statistics** | Sequence statistical analysis |
| 🔄 **translation_functions** | DNA/RNA/protein translation and transcription |

## 🚀 Example Usage

```python
from bioutils_collection.bioinformatics_functions import (
    reverse_complement,
    calculate_gc_content,
    needleman_wunsch_align,
    parse_fasta,
)

# Reverse complement
seq = "ATCGATCG"
rev_comp = reverse_complement(seq)  # "CGATCGAT"

# GC content
gc = calculate_gc_content("ATCGATCG")  # 0.5

# Sequence alignment
seq1 = "GATTACA"
seq2 = "GCATGCU"
aligned1, aligned2, score = needleman_wunsch_align(seq1, seq2)

# Parse FASTA
sequences = parse_fasta("sequences.fasta")
for header, seq in sequences:
    print(f"{header}: {len(seq)} bp")
```

## 📚 Key Features

### Sequence Analysis
- DNA/RNA/protein validation
- GC content and AT/GC skew analysis
- Codon usage analysis
- Sequence complexity metrics

### Alignment & Comparison
- Needleman-Wunsch (global alignment)
- Smith-Waterman (local alignment)
- Sequence similarity metrics
- Multiple sequence alignment utilities

### Annotation & Parsing
- FASTA file I/O with quality control
- BED file parsing and statistics
- VCF file parsing
- GFF annotation handling

### Pattern & Motif Analysis
- Motif discovery algorithms
- Restriction enzyme site finding
- Repeat sequence detection
- Pattern matching utilities

### Translation & Transcription
- DNA to RNA transcription
- RNA to DNA conversion
- DNA/RNA to protein translation
- Codon table support

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
- **Type hints** for all parameters and returns
- **NumPy-style docstrings** with detailed descriptions
- **Usage examples** with expected outputs
- **Complexity notes** for performance-critical operations
- **References** to algorithms and methods

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🔗 Related Projects

- [pyutils-collection](https://github.com/MForofontov/pyutils-collection) - General Python utilities library
- [BioPython](https://biopython.org/) - Comprehensive bioinformatics toolkit
- [scikit-bio](http://scikit-bio.org/) - Scientific Python library for bioinformatics

## 📮 Contact

**Author:** Mykyta Forofontov  
**Repository:** https://github.com/MForofontov/bioutils-collection  
**Issues:** https://github.com/MForofontov/bioutils-collection/issues
