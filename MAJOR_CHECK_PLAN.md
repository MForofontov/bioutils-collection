# Major check plan

This is a staged engineering check. Completion requires evidence from the commands
below; a passing unit suite alone does not establish scientific correctness.

## 1. Installation and packaging (highest priority)

- Parse `pyproject.toml`, build wheel and source distributions, and run
  `python -m twine check dist/*`.
- Install with both `pip install -r requirements.txt` and `pip install -e ".[dev]"`.
- Run `python -m pip check` and import the installed wheel outside the checkout.
- Check Python 3.10–3.13 on Linux and Windows, including Python 3.10 dependency
  resolution and minimum supported dependency versions.
- Confirm the wheel includes the typing marker promised by the package metadata.

Exit criteria: clean installations, valid distributions, all public imports work.

## 2. Functional and scientific correctness

- Run `python -m pytest -q --cov=bioutils_collection --cov-report=term-missing`.
- Add regression tests for confirmed bugs before changing behavior.
- Audit coordinates and boundary conventions in BED/GFF/GTF/VCF conversions.
- Compare alignment, translation tables, codon statistics, and protein properties
  against independently checked small fixtures and established reference tools.
- Cover empty inputs, ambiguous bases, lowercase input, malformed files,
  overlapping motifs, duplicate headers, and extreme parameter values.
- Exercise serial and multiprocessing translation, including Windows process
  startup, custom codon tables, errors, and ordering.

Exit criteria: passing regressions and documented coordinate/biological conventions.

## 3. File safety and streaming

- Check same-path input/output, symbolic links, interrupted writes, invalid FASTA
  headers, multiline records, and malformed FASTQ.
- Ensure an invalid parameter cannot truncate an existing output file.
- Measure peak memory for streaming APIs with many records and one huge record.

Exit criteria: no input data loss; streaming memory limits are documented.

## 4. Performance

- Benchmark representative small, medium, and large sequences; separate Numba
  compilation time from warmed execution and multiprocessing startup costs.
- Check GC windows at overlapping, disjoint, and gapped step sizes against a
  straightforward reference implementation before measuring speed.
- Profile minimizers, syncmers, k-mer counting, alignment, and FASTA processing.
- Record Python/dependency versions, inputs, repeated timings, and peak memory.
- Preserve public return values, ordering, rounding, and error behavior.

Exit criteria: measured improvements with equivalent outputs; no timing-based
assertions in ordinary CI.

## 5. Static checks, documentation, and release gates

- Run `python -m ruff check .`, `python -m ruff format --check .`, and
  `python -m mypy bioutils_collection`; triage existing failures by impact.
- Execute README examples and verify documented imports and return ordering.
- Require automated tests on pull requests and before the publication build.
- Review changelog, package version, license, and distribution contents.

Exit criteria: checks pass or remaining issues are explicitly recorded with scope.

## Current pass

- Fixed invalid TOML in the author profile URL key.
- Centralized runtime requirements in package metadata, removed conflicting pins,
  and allowed NumPy 2.0+ so the declared Python 3.10 support can resolve.
- Raised the setuptools floor to support the license metadata syntax in use.
- Optimized GC windows by counting only incoming/outgoing bases for overlaps.
- Added reference-equivalence tests for GC windows and file-preservation tests
  for streaming translation, including hard-link aliases and invalid tables.
- Corrected and simplified README examples without changing public APIs.
- Cleared Ruff lint/format checks and strict mypy checks for 103 source files.
- Added a Python 3.10–3.13/Linux and Python 3.12/Windows functional test matrix,
  static checks, package checks, and a publication dependency on these checks.
- Built valid wheel/source distributions and verified the wheel's typing marker.
- Validated an installed wheel import and executed every README Python example.
- Functional validation on Windows/Python 3.12: 974 passed, 9 slow tests excluded,
  with 94% package coverage. BLAS threads were limited to one for the final run
  after concurrent checks exhausted the host's paging-file capacity.
- The nine slow performance tests also passed in a separate run (983 total tests
  passed across the two runs). The CI platform matrix must pass for each
  published release.
- Replaced the private Numba kernel's reflected list with a tuple; targeted
  backend-equivalence tests treat this deprecation warning as an error.
  The clustering warning for identical motifs reflects duplicate input points
  and does not indicate a failed test.
- A local benchmark on 100,000 bases (window 1,000, step 10; best of three runs)
  matched the direct reference output and improved from 1.6113 s to 0.0290 s.
  This measurement is specific to this input and environment.
- Keep API names, return types, ordering, rounding, and coordinate conventions.
  Record scientific corrections that require behavior changes for a later release.
- Remaining phases above are the broader audit backlog, not claims of completion.

## File handling and translation follow-up

- Normalize FASTQ line endings before record validation and quality filtering;
  cover LF, CRLF, and CR inputs with and without quality filtering.
- Validate in-memory FASTA translation tables before reading records or opening
  output, including empty input; preserve existing output on configuration errors.
- Compare fast and Python translation for all 64 codons under every registered
  genetic-code identifier, plus a sparse custom table with multi-character values.
- Load scikit-learn on demand inside motif clustering. A fresh-process smoke
  check confirms basic package imports do not load scikit-learn.
- Targeted validation: 101 tests passed; Ruff, formatting, and mypy passed.
- Final full run on Windows/Python 3.12: 1,047 passed, including the slow tests,
  with 96% package coverage. Numba pending-deprecation warnings were treated as
  errors; only the expected identical-motif clustering warning remains.
- Python 3.10 dependency resolution succeeded in a binary-wheel dry run with
  NumPy 2.2.6, SciPy 1.15.3, scikit-learn 1.7.2, and Numba 0.68.1. This checks
  dependency compatibility, not execution on Python 3.10.
- Rebuilt wheel/source distributions and passed both metadata checks.

## Release 0.3.1 validation (2026-10-30)

- Repeated the full suite for the release version: 1,047 passed, 96% coverage,
  with only the expected identical-motif clustering warning.
- Ruff lint/formatting, strict mypy (103 source files), dependency consistency,
  distribution metadata, and isolated installed-wheel smoke checks passed.
- Publishing uses GitHub's trusted-publisher workflow after cross-platform CI.

## Major audit implementation (2026-10-30)

### Baseline and method

- Pulled `origin/main` with `--ff-only` and created `codex/major-audit` from
  `9e1cc06420f14940f9e640af554539be29886fb5`; the starting tree was clean.
- Baseline functional run: 1,038 passed, nine slow tests excluded, one expected
  identical-motif clustering warning. Baseline Ruff lint/format, strict mypy
  (103 source files), and dependency consistency passed. The earlier release's
  measured coverage was 96%; this baseline run did not remeasure coverage.
- Local interpreter: Windows/Python 3.12.14; NumPy 2.5.4, SciPy 1.18.1,
  scikit-learn 1.9.1, Numba 0.68.1. BLAS/OpenMP thread counts were limited to one
  for tests and benchmarks.
- Added regressions before fixes. Initial focused runs reproduced 14 selector,
  annotation, and file-safety failures; subsequent scientific/format checks
  reproduced additional defects rather than relying on baseline test coverage.

### Correctness and compatibility

- Alignment: enumerated all alignment-path scores for short binary strings and
  checked global/local optima, traceback scores, and sequence reconstruction.
- Translation: checked every supported numeric table against all 64 codons in an
  independent [NCBI AAs snapshot](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
  (source update 2024-09-23, retrieved 2026-10-30). Fixed tables 24/33 AGG and
  table 31 TGA. Existing cross-backend checks cover aliases and sparse custom
  tables; additional fixtures exercise serial/process ordering and uneven chunks.
- Sequence operations: fixed closed-syncmer boundary selection and leftmost ties
  against the [original syncmer definition](https://peerj.com/articles/10805/).
  Exhaustively checked short inputs. Fixed zero-density sketches and spurious
  streaming suffix windows; checked normal/canonical sketch results against
  direct windows, including repeats, short inputs, and nonconsecutive positions.
- Annotation: fixed GFF/GTF parser-to-writer round trips, GFF-to-BED field aliases,
  GFF filtering/statistics, BED metadata/blank lines, and GFF embedded FASTA.
  GTF/VCF malformed columns now raise line-numbered errors instead of silently
  dropping records. GFF encoded attribute values and list separators round trip.
- Coordinates: preserved native parser coordinates. Documented one-based,
  inclusive input for BED conversion and annotation statistics; BED parser
  records need start incremented before those APIs. Checked single-base features.
  Statistics sum lengths rather than calculating an interval union. See
  [UCSC BED conventions](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and the
  [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
- Files: shared FASTA parsing, rejected missing identifiers/sequence-before-header,
  and retained documented empty/duplicate-record behavior. Both file translators
  reject input aliases and use same-directory temporary output with flush/fsync,
  close, then replacement. Error tests cover processing, fsync, and replacement,
  and check output preservation and temporary-file cleanup. Existing destination
  modes are preserved; output symlinks continue to refer to the updated target.
- FASTQ: retained newline equivalence; rejected missing/mismatched identifiers,
  nonprintable Phred+33 quality characters, and trailing quality whitespace.
- Motifs/repeats/restriction: checked independent overlap and palindrome fixtures;
  corrected IUPAC matching to include overlaps, consistent with literal matching.
- Statistics/protein properties: fixed CAI zero-weight behavior and finite [0, 1]
  validation. Replaced pKa averaging with a zero-charge root solver using the
  existing pKa values, including C/Y and termini. Corrected bond/terminal water
  in approximate masses and protein validation. Added charge-balance, analytical
  mass, codon-frequency/RSCU/ENC extreme, entropy, and conservation fixtures.
  Charge balance and polymer condensation follow the methods illustrated by
  [Biopython's pI implementation](https://biopython.org/docs/latest/api/Bio.SeqUtils.IsoelectricPoint.html)
  and [sequence utilities](https://biopython.org/docs/latest/api/Bio.SeqUtils.html);
  different pKa sets and rounded mass tables intentionally produce different values.
- Clustering: verified an unambiguous duplicate-motif partition without depending
  on arbitrary cluster-label numbering. Data validation's existing GC denominator
  convention was preserved and its misleading examples corrected.
- No public signatures, return types, dependencies, or export
  names were changed. Corrected scientific outputs and validation are documented
  in the changelog; exact previous results are not promised for those bugs.

### Measured performance

- Reproduce with `python benchmarks/major_audit.py --output
  benchmarks/major_audit_results.json`. The checked-in JSON records 129 cases,
  deterministic seed 20261030, three timing repetitions, and separate peak Python
  allocation measurements. Inputs cover 1k/10k/100k bases, random/repetitive DNA,
  capped quadratic alignments, uneven batches, and FASTA files up to 1M bases.
- References use corrected biological/window semantics; optimizing an incorrect
  result is not counted as an improvement. Timings are environment-specific and
  shared-host load is not controlled. No speed thresholds were added to CI.
- At 100k bases (k=31, s=9), deque syncmer selection was about 11x faster on
  random input and 6–8x faster on repetitive input than direct s-mer windows.
  Counting density without positions reduced peak Python allocations from about
  498 KB to 103 KB (random) and 1.12 MB to 103 KB (repetitive), with no claim of
  a density-count runtime improvement on every input.
- At 100k bases (k=31, w=32), normal/canonical sketches improved about 1.7–3x
  against direct-window references, preserving ordering and ties.
- The initial k-mer candidate consumed a generator. On repetitive 100k input it
  reduced intermediate peak allocations from about 8 MB to 1.2 KB. Random input
  had similar peak allocations and variable/slower timings; the benefit is memory
  for low-diversity inputs, not a universal speed improvement.
  The paired hardening pass below subsequently rejected this candidate and
  restored the original counter under the agreed runtime/memory policy.
- Import startup, fresh-process Numba initialization, isolated empty-cache first
  calls, warmed translation, and process-pool startup are recorded separately.
  The latter measurements are not recommendations to parallelize tiny batches.
- Allocation measurements exclude preallocated inputs, native arrays, and worker
  RSS. Streaming memory is bounded by the largest FASTA record plus translation
  buffers; a single huge record is not processed in constant memory.

### Explicit model limits and remaining release gates

- ENC remains a population-homozygosity approximation without finite-sample
  correction; sparse sequences are not reliable unbiased estimates. Sequence
  complexity is a unique-substring fraction, not alphabet-normalized linguistic
  complexity. Tm remains a composition approximation. These limits are now explicit
  rather than presented as more sophisticated algorithms.
- Translation does not model alternative initiators or context-dependent stops;
  NCBI table 32 is unsupported. Protein pI uses fixed pKa values, and molecular
  masses use 0.1-Da-rounded tables without modifications. Minimal annotation
  parsers do not validate every format semantic or expose VCF sample columns.
- Local validation results are recorded below after the final run. The remote
  Linux/Python 3.10–3.13 and Windows/Python 3.12 matrix must still run on this
  branch before merge; local Windows checks cannot establish those results.
- No push, merge, or publication is part of this audit. On 2026-11-01, the user
  added documentation/version maintenance to the scope. Version is now 0.3.2;
  final distributions must be rebuilt and checked with that version.

### Final local validation

- Runtime audit on 2026-10-30: 1,100 functional tests passed, four symlink tests
  skipped because Windows symlink privileges were unavailable, and nine slow
  tests excluded from that run. The separate slow run passed all nine tests:
  1,109 passed in total. Coverage was 96.58% (97% rounded). The expected
  identical-motif clustering warning was the only functional-suite warning;
  Numba pending-deprecation warnings were treated as errors.
- Final 2026-11-01 checks after the version/instruction update: Ruff lint and
  formatting passed, strict mypy passed for 106 source files, and `pip check`
  reported no broken requirements. Package modules also parsed with Python
  3.10 syntax rules; this does not establish execution on that interpreter.
- Built 0.3.2 wheel and sdist and passed Twine metadata checks. The isolated
  build-dependency installation stalled; the final build used cached setuptools
  with `--no-isolation`, without changing runtime dependencies.
- Installed the final wheel into a temporary directory outside the checkout.
  Verified version 0.3.2 in both runtime and distribution metadata, imported all
  106 modules, resolved every declared public export, and executed both README
  Python examples. Checked wheel source equality, typing marker, and sdist license
  and typing-marker contents.
- Validated all four new/updated project skills across bioutils/pyutils, checked
  local links in changed guidance/documentation, and passed `git diff --check`.
  Both repositories now require documentation and one version update per completed
  changeset. The default chosen for this task bumps bioutils to 0.3.2; pyutils's
  package version remains 0.2.0 while its guidance is updated.
- Remote matrix execution remains pending. No changes were committed, pushed,
  merged, tagged, or published by this audit.

## Validation and hardening follow-up (2026-11-01)

### Configuration and transactional files

- Added empty-input regressions before the fix: 17 failures demonstrated that
  fast/batch/chunked translation could ignore invalid genetic-code identifiers,
  zero/invalid process counts, or chunk sizes below one codon. Valid empty results,
  sparse dictionaries, and multi-character custom translations are retained.
- Validate those settings before empty-input returns. Public signatures, exports,
  return types, and the documented approximate scientific models are unchanged.
- Expanded transactional-file tests to partial write and close failures,
  KeyboardInterrupt, missing destinations, preserved permission bits, and cleanup.
  Linux symlink errors now fail rather than silently skipping the tests.
- Windows functional run: **1,146 passed**, four symlink privilege skips, nine
  slow tests excluded; the separate slow run passed all nine, **1,155 total**.
  Coverage: **96.72%**. Only the expected identical-motif warning remains.
- A disposable Linux Python 3.12 container passed **57 targeted checks with no
  skips**, including all four symlink cases. The checkout was mounted read-only;
  code and test inputs were copied into the container. Image digest:
  `sha256:aa943ce411423c7582288f18707ff66eddcbaeb510502d747a72aeb34f30d159`.
  This is local Linux evidence, not a completed remote compatibility matrix.
- Slow translation tests retain their large-input/backend/process equivalence
  assertions. Hardware-dependent speed assertions now belong to explicit benchmarks.

### Paired performance evidence

- [Final paired results](benchmarks/paired_audit_results.json) record 97
  comparisons against `9e1cc06420f14940f9e640af554539be29886fb5`, complete tracked
  and nonignored-untracked source manifests, five counterbalanced timing rounds,
  separate allocation/native process peaks, and import/JIT startup measurements.
  The runtime package hash matches the final source files. Full-tree manifests
  describe measurement-time state; later documentation/tooling edits are separate.
- Of these comparisons, 85 have equivalent outputs. The other 12 compare old
  incorrect closed-syncmer or streaming-suffix results and explicitly carry no
  speedup claim. Corrected-reference comparisons are separately labeled.
- On 100k random/repetitive fixtures, equivalent baseline comparisons measured
  about **6x** for open syncmers, **3.6–3.9x** for minimizer sketches, and **1.6x**
  for canonical sketches. Corrected-reference comparisons measured **6.6–7.7x**
  for closed syncmers and **3.4–4.3x** for streaming sketches. These are local
  measurements, not portable speed guarantees.
- [Density supplement](benchmarks/paired_density_results.json) adds four
  equivalent medium/large open-density comparisons with separate native/Python
  memory measurements. All 24 equivalent medium/large retained selector cases
  across the final report and supplement satisfy the balanced retention gate.
- [Rejected candidate evidence](benchmarks/paired_candidate_results.json) records
  the generator counter experiment. Its medium-input median runtime ratios were
  1.21–1.24 despite allocation savings. It did not demonstrate a pass under the
  allowed 10% runtime cost, so the original counter was restored.
- Shared-host load is uncontrolled. Some unchanged controls vary substantially
  in timing, so a threshold label alone is not evidence of a code regression.
  Original counter results are controls, not a newly retained optimization.
  Native RSS includes imports, prepared inputs and JIT peaks; pre-operation peaks
  are reported. Python allocations exclude prepared inputs and are measured in
  separate untimed calls. File-output digests normalize CRLF to LF; source
  manifests hash raw bytes.
- The measuring program, references, and package sources are frozen in temporary
  snapshots for new runs. The density supplement also exercises this control
  isolation. No benchmark resets or checks out the working tree.

### Reproducible distribution gates

- Added tracked [distribution verification](scripts/check_distribution.py),
  replacing the ignored scratch verifier. It derives version from `_version.py`,
  rejects stale versions, inventories and source contents before installation,
  checks both distribution names/metadata, typing and licenses, and writes JSON status
  reports including failures. Negative artifact regressions are in the unit suite.
- The default verifier creates separate fresh virtual environments for wheel and
  sdist, clears inherited Python-path overrides, resolves dependencies, runs
  `pip check`, verifies installed provenance and version, imports all 106 modules,
  and executes both README Python examples outside the checkout.
- Rebuilt **0.3.2** wheel/sdist with installed setuptools 84.2.1 and
  `--no-isolation`; Twine passed both. The tracked verifier passed both fresh
  installations with fresh dependency resolution and isolated sdist build
  dependencies. No offline/dependency-reuse mode was used for this final gate.
- Final artifact SHA-256 values: wheel
  `d062ea658a278896fb778e060ecd0056aeb8a1fd6026af8c4360bdf0dc1b21b2`;
  sdist `44b6dbcf7a1167dcfd683535bbeedb230b9322d47867f287ef6e12c96688947e`.
  Local XML, coverage, build and installation reports are under ignored `dist/`.
- Packaging CI now calls the tracked verifier and uploads its report. Existing
  Linux Python 3.10–3.13 and Windows Python 3.12 jobs remain in place. They must
  run on the final pushed revision before merge; they have not run remotely here.
- Ruff lint/format, strict mypy (106 source files), and dependency consistency
  passed. Documentation and shared workflows now describe these commands and the
  balanced gate. Version remains 0.3.2 for this same uncommitted changeset; no
  pyutils files were changed in this continuation. No commit/push/publication.


### Release validation — 2026-11-02

Publication was explicitly requested for version **0.3.2**. Fresh Windows
Python 3.12.14 validation passed **1,146 functional tests and 9 slow tests**,
with four Windows symlink privilege skips and 96.72%
coverage. Ruff lint/format (217 files), strict mypy (106 source files), and
`pip check` passed. Package source hashes and dependency versions still match
the retained paired benchmark evidence. The prior Linux symlink run remains
57 passed with no skips; fresh remote platform results are a separate gate.

Fresh isolated wheel/sdist builds, Twine metadata checks, and separate fresh
wheel/source installations passed, including all public imports, provenance,
typing markers, licenses, version consistency, and both README examples.
Evidence is recorded locally under `dist/release-2026-11-02/`.

The release procedure pushes `codex/major-audit`, requires the exact pushed
commit to pass the complete remote CI matrix, and then publishes tag `v0.3.2`
through the existing PyPI trusted-publishing workflow. Remote CI and actual
PyPI artifact verification must succeed before publication is reported complete.
