# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2026-01-28

### Added
- **translate_dna_fast** - Numba JIT-compiled translation for 2.7x speedup on large sequences
- **translate_batch** - Parallel translation of multiple sequences using multiprocessing
- **translate_large_sequence** - Chunk-based parallel translation for sequences >10MB
- **translate_fasta** - Translate sequences from FASTA files with in-memory processing
- **translate_fasta_streaming** - Memory-efficient streaming translation for huge FASTA files
- **genetic_code_tables** - Support for all 27 NCBI genetic code tables (1-6, 9-16, 21-31, 33)

### Removed
- **phylogenetic_distance** - Redundant with `hamming_distance(s1, s2) / len(s1)` (simple division wrapper)
- **dinucleotide_frequency** - Redundant with `kmer_frequency(seq, 2)` (identical functionality)
- **nucleotide_frequency** - Redundant with `kmer_frequency(seq, 1)` (identical functionality)
- **sequence_masking** - Incorrect implementation (only masked simple runs, use `remove_low_complexity_regions` instead)
- **reverse_sequence** - Trivial one-liner (`seq[::-1]`)

### Fixed
- **blast_score_ratio** - Added validation for negative target_score values
- **restriction_site_finder** - Added type checking and empty validation for parameters

### Changed
- Updated all `__init__.py` files to remove deleted function exports
- Removed 5 test files for deleted functions

### Verified
- All remaining ~70 functions verified as unique and non-redundant
- 692 tests passing, 6 skipped
- No broken dependencies or import errors

## [0.1.1] - Previous Release

Initial stable release with 77+ bioinformatics functions.
