# Changelog

All notable changes to this package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-08-04

### Added
- **generate_minimizers** - Core minimizer generation using sliding window (Roberts et al. 2004)
- **canonical_minimizers** - Strand-independent minimizers using reverse complement
- **minimizer_sketch** - Compact sequence representation with position mapping
- **minimizer_positions** - Extract all minimizer positions from sequence
- **unique_minimizer_positions** - Deduplicated minimizer positions
- **minimizer_density** - Calculate minimizer density/compression ratio
- **canonical_minimizer_sketch** - Canonical minimizers with position mapping
- **minimizer_sketch_from_generator** - Memory-efficient streaming minimizer sketch
- **generate_syncmers** - Syncmer variant generation with open/closed methods (Edgar 2021)
- **syncmer_positions** - Extract all syncmer positions from sequence
- **syncmer_density** - Calculate syncmer density/compression ratio
- **sequence_to_kmers_with_positions** - Split sequence into (kmer, position) tuples
- 240 new unit tests for all minimizer/syncmer/utility functions (20 tests per function)

### Changed
- **generate_minimizers** - Refactored to use `sequence_to_kmers` instead of inline list comprehension
- **kmer_frequency** - Refactored to use `sequence_to_kmers` instead of inline generator
- **sequence_statistics** - Refactored to use `gc_content` instead of inline GC calculation
- **minimizer_positions** - Refactored to use `sequence_to_kmers_with_positions`
- **unique_minimizer_positions** - Refactored to use `sequence_to_kmers_with_positions`
- **minimizer_sketch** - Refactored to use `sequence_to_kmers_with_positions`
- **find_cpg_islands** - O(n) sliding-window detection with merged overlapping islands
- **levenshtein_distance** - Two-row DP for O(min(n,m)) space complexity
- **codon_adaptation_index** - `reference_weights` is now required (breaking)
- Added academic References sections to 9 algorithmic functions:
  - **needleman_wunsch** - Needleman & Wunsch (1970)
  - **smith_waterman** - Smith & Waterman (1981)
  - **generate_minimizers** - Roberts et al. (2004)
  - **generate_syncmers** - Edgar (2021)
  - **calculate_cai** - Sharp & Li (1987)
  - **calculate_enc** - Wright (1990)
  - **levenshtein_distance** - Levenshtein (1966)
  - **melting_temperature** - Wallace (1979), SantaLucia (1998)
  - **find_cpg_islands** - Gardiner-Garden & Frommer (1987)

### Fixed
- **translate_large_sequence** - Raise `ValueError` when `chunk_size` rounds below 3 (was crashing with `range(..., 0)`)
- **translate_dna_fast** - Reject invalid DNA bases instead of silently mapping to `X`
- **find_orfs** - Truncate ORF sequences to last complete codon boundary
- **blast_score_ratio** - Clamp result to documented 0–1 range
- **needleman_wunsch** / **smith_waterman** - Reject `bool` values for score parameters
- **pairwise_identity** - Added strict=True to zip() calls (B905 linting)
- **canonical_minimizers** - Replaced inline reverse complement dict with `reverse_complement` function
- **canonical_minimizer_sketch** - Replaced inline reverse complement dict with `reverse_complement` function
- **palindromic_sequence_finder** - Replaced private `_reverse_complement` helper with `reverse_complement` function
- **fasta_reverse_complement** - Replaced inline `str.maketrans` with `reverse_complement` function
- Code formatting improvements across 22 files using ruff

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

## [0.1.1] - Previous Release

Initial stable release with 77+ bioinformatics functions.
