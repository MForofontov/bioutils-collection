# Changelog

All notable changes to this package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.3] - 2026-11-08

Release prepared on 2026-11-08; publication is tracked by the GitHub release
and PyPI publishing workflow.

### Fixed
- Preserve BED records whose chromosome is named `track` or `browser` while
  continuing to skip metadata lines.
- Reject non-finite BLAST reference and target scores instead of returning NaN
  or treating an infinite target score as a perfect match.
- Count lowercase G/C bases in FASTA statistics while preserving the mean of
  per-record GC fractions and the existing empty/ambiguous-record denominators.
- Normalize GenBank LF/CRLF/CR line endings and remove tabs between sequence
  groups rather than emitting whitespace in FASTA sequences.
- Reject empty restriction recognition sites and non-string site elements.
- Reject NaN low-complexity masking thresholds, including on short/empty input.
- Correct executable alignment, FASTA statistics, and entropy-masking docstring
  examples.
- Correct ten malformed FASTA docstrings and 38 incorrect example results across
  sequence selection, annotation, motifs, repeats, validation and statistics.
- Make stochastic, file-based and large-job examples executable with small,
  temporary inputs and invariant assertions.

### Changed
- Reuse scoring-row references in global alignment while preserving full outputs
  and tie-breaking. Reject the local-alignment candidate because confirmation
  did not meet the retention gate on every medium/large fixture. Timings are
  environment-specific.
- Prepare version 0.3.3 with PyPI version, supported Python versions, and MIT
  license badges in the README, plus a summary of the package's 80+ functions
  across 12 module families.
- Explicitly identify version 0.3.3 in the README, link its release notes, and
  distinguish the source version from the live PyPI badges.

### Added
- Fresh main-sync audit evidence, BED/score-ratio regressions, a five-round paired
  benchmark with 103 equivalent outputs, and Linux symlink/file-safety checks.
- Regression fixtures for format whitespace, lowercase GC counts, invalid
  recognition sites, and masking thresholds, plus fresh paired audit evidence.
- A bounded documentation checker, exercised by the unit suite, and stronger
  independent masking, IUPAC, selector and alignment-tie fixtures.
- Calibrated square/unequal alignment benchmarks with unchanged controls,
  complete-output comparisons, and explicit inconclusive timing decisions.

## [0.3.2] - 2026-11-02

### Fixed
- Validate table identifiers, process counts, and codon-adjusted chunk sizes even
  for empty fast, batch, and chunked translation inputs.
- Corrected closed syncmers to select boundary minima with leftmost tie-breaking.
- Corrected NCBI tables 24/33 (AGG to lysine) and 31 (TGA to tryptophan).
- Made GFF/GTF writers accept native parser records; GFF records now work with
  BED conversion, filtering, and feature statistics without renaming fields.
- Handled BED track/browser headers, blank annotation lines, and embedded GFF
  FASTA sections. Malformed GTF/VCF records now raise errors instead of disappearing.
- Prevented both FASTA translators from overwriting input aliases; output is now
  written transactionally and preserved after processing, write, or commit errors.
- Rejected empty FASTA/FASTQ identifiers, mismatched FASTQ separator identifiers,
  and invalid Phred+33 characters. Both file translators share FASTA parsing.
- Included overlapping IUPAC matches, consistent with literal matching.
- Corrected CAI zero-weight handling and rejected non-finite/out-of-range weights.
- Replaced pKa averaging with a Henderson-Hasselbalch zero-charge pI solver,
  retaining the existing fixed pKa set and two-decimal output.
- Accounted for peptide-bond water and nucleotide termini in approximate masses,
  including the molecular weight reported by protein validation.
- Removed spurious trailing windows from streaming minimizer sketches and made
  zero-density sketches empty rather than raising division-by-zero errors.

### Changed
- Keep the previous k-mer counting implementation: paired measurements rejected
  the generator candidate's medium-input runtime penalty under the balanced policy.
- Keep hardware speed thresholds in explicit benchmarks; slow translation tests
  now check backend/process equivalence without machine-dependent speed assertions.
- Shared a monotonic-deque syncmer selector across generation, positions, and
  density; density no longer retains every selected position.
- Used deque selection for normal/canonical minimizer sketches and streaming
  sketches while preserving complete-window ordering and leftmost ties.
- Documented coordinate conventions, approximate scientific models, and corrected
  misleading examples. Public function signatures and return types are unchanged.

### Added
- Tracked wheel/sdist verification in separate virtual environments, wired into
  packaging CI, with source/version checks and installed README/public-import checks.
- Paired baseline/current benchmarks with complete source snapshots, five
  counterbalanced timing rounds, separate Python/process memory measurements,
  and explicit correctness and runtime/memory retention decisions.
- File-safety checks for partial writes, close failures, cancellation, missing
  destinations, permission preservation, and mandatory Linux symlink coverage.
- Canonical repository instructions, shared audit/documentation skills, and an
  always-applied rule requiring documentation and version updates with changes.
- Exhaustive small-input and analytical regression fixtures, including an
  independent NCBI snapshot covering all supported genetic-code identifiers.
- Reproducible correctness-first benchmarks for selectors, GC windows, k-mer
  counting, alignment, translation, and FASTA processing, with separate memory
  measurements and process/Numba startup measurements.

## [0.3.1] - 2026-10-30

### Fixed
- Fixed invalid TOML in project URLs and aligned the setuptools requirement with
  the license metadata syntax used by the package.
- Centralized runtime/development requirements and allowed NumPy versions
  compatible with Python 3.10.
- Prevented streaming FASTA translation from overwriting its input, including
  hard-link aliases, and validated translation tables before opening output.
- Validated in-memory FASTA translation tables before processing records or
  overwriting output, including for empty input files.
- Normalized FASTQ line endings before quality filtering and FASTA conversion.
- Replaced the private Numba kernel's reflected list with an immutable tuple.
- Corrected README examples for GC percentages, alignment return ordering,
  FASTA text input, and ORF length filtering.
- Resolved lint/formatting issues and updated Numba decoration for typed stubs.

### Changed
- Reduced overlapping GC window analysis to linear time while preserving output
  boundaries, rounding, and validation.
- Deferred scikit-learn loading until motif clustering is requested, avoiding
  its import overhead for unrelated sequence utilities and translation workers.
- Made mypy use the running interpreter version to match installed dependency
  stubs; Python compatibility remains checked through the test matrix.

### Added
- Regression tests for GC window equivalence and streaming translation file safety.
- Translation backend equivalence checks across all registered genetic code
  identifiers, custom codon tables, and FASTQ newline regression tests.
- Packaged the PEP 561 typing marker.
- Added test, static analysis, and packaging workflows required before publication.
- Added a staged major check plan with explicit acceptance criteria.

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
