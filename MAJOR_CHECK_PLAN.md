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

## Maintenance audit 0.3.3 — 2026-11-02

### Baseline and preservation

- Started on local `main` at `186a85fbbd70ac31cb3b013281cb7691d4365c59`.
  Saved the existing README/changelog/version diff and dependency inventory under
  ignored `dist/maintenance-audit/`. Fetch and fast-forward-only synchronization
  found no newer remote `main`; the local 0.3.2 release commit was preserved.
  Work continues on `codex/maintenance-audit`, including the existing 0.3.3 bump.
- Used the explicit `.venv/Scripts/python.exe` interpreter (Windows/Python
  3.12.14), editable `.[dev]` installation, and one BLAS/OpenMP thread. Baseline
  functional tests: **1,146 passed, four symlink privilege skips, nine slow tests
  deselected**, with the expected duplicate-motif clustering warning.

### Review and reproduced defects

- Reviewed the 12 module families: annotation coordinates/round trips, file
  translation and atomic replacement, genetic-code/backend/process equivalence,
  alignment scoring/traceback, selectors/GC windows, codon/protein statistics,
  motif/repeat/restriction matches, clustering, validation, and FASTA formats.
  Existing independent NCBI, analytical, exhaustive-short-input, and file-failure
  fixtures remain part of final validation; existing model limits still apply.
- Added regressions before fixes: **12 failures and two passing controls**.
  Corrected GenBank CR line endings and sequence tabs; FASTA statistics now count
  lowercase G/C while retaining mean-per-record fractions, ambiguous-base
  denominators and empty-record behavior. Restriction sites reject empty and
  non-string values. Masking rejects NaN even for empty/short input.
- Focused format/validation checks: **186 passed**. An initial cross-module run
  had a Numba worker import `MemoryError` on the shared host; its thread-limited
  rerun passed **72 checks**. Neither tests nor process-equivalence assertions
  were weakened. Local logs retain the unsuccessful attempt.
- Executing affected docstrings exposed incorrect alignment and entropy-masking
  examples. Corrected those examples and clarified case handling, GC fractions,
  literal restriction matching, and the union/strict-threshold masking rule.
  Public signatures, exports, coordinate conventions and return shapes remain.

### Profiling and retained candidates

- Profiled alignment, minimizer sketches, k-mer counting, translation and FASTA
  parsing on deterministic 1k/10k/100k input. Alignment loops repeatedly indexed
  scoring rows and characters; tried caching these references and using
  independent zero-filled rows. Only global alignment is retained after
  confirmation below. Scoring and traceback selection rules are unchanged.
- Compared **1,176 complete alignment tuples per algorithm** against the starting
  implementation across short binary sequences and six scoring configurations,
  including zero scores, positive gaps, and negative match scores: all matched.
  Existing exhaustive alignment-path fixtures independently check optimum scores.
- [Focused paired measurements](benchmarks/maintenance_alignment_results.json)
  use five counterbalanced rounds against the explicit starting commit. All 12
  outputs match; all eight medium/large cases meet the runtime retention gate.
  Initial global runtime ratios were 0.56–0.78 and local ratios 0.61–0.86.
  One tiny local case was neutral (0.95), so no universal speedup is claimed.
  Alignment matrices are capped: medium/large labels correspond to 100-by-50
  and 300-by-150 comparisons rather than quadratic 10k/100k-base matrices.
- Each report includes raw source manifests, dependency versions, timing samples,
  Python allocation measurements and separate process RSS. RSS includes imports
  and prepared inputs and is not described as Python allocation savings.
  The first broad benchmark attempt failed inside a worker; its rerun captures
  worker diagnostics. Unchanged controls are not new optimization claims.
- [Source equivalence](benchmarks/maintenance_source_equivalence.json) records
  the final raw package manifest and both measured package hashes. FASTA
  statistics docstring escaping changed before confirmation. After confirmation,
  only the rejected local-alignment candidate was restored: its executable AST
  matches the starting implementation. Every other package file, including the
  retained global alignment, matches confirmation byte-for-byte. Historical
  candidate snapshots are not presented as the complete final package.
- [Broad paired results](benchmarks/maintenance_paired_results.json) contain
  **103 equivalent comparisons**, including separate import/JIT/process startup
  measurements. Broad timing ratios conflict with the focused run: even unchanged
  large k-mer counters varied to 1.80–1.89 and unchanged large edit distance to
  1.15–2.21. These are uncontrolled shared-host measurements, not evidence of
  changes to those implementations. A focused confirmation with an unchanged
  edit-distance control determines final retention below.
- [Confirmation results](benchmarks/maintenance_confirmation_results.json) use
  **ten counterbalanced rounds** with an unchanged edit-distance control and
  separate memory probes. All 12 outputs match. All four medium/large global
  cases pass: runtime ratios **0.745–0.805** (about **19–25% less time**).
  The local candidate misses the gate on the large repetitive case (runtime
  0.986, Python allocation ratio 0.909), so it is **rejected** and the original
  implementation restored. No local-alignment speedup is retained or claimed.
  The unchanged control still varies (ratios 0.99–1.44); measurements remain
  specific to this shared host, with no portable speed guarantee.

### Final validation

- Final functional run: **1,161 passed, four Windows symlink privilege skips,
  nine slow tests deselected**, with **96.74%** package statement coverage.
  The sole warning is the expected duplicate-motif clustering warning.
  The separate slow run passed **all nine tests** (1,170 passes across both runs).
- Ruff lint and formatting passed (218 Python files); strict mypy passed for
  106 source files, and `pip check` passed. All **15 affected doctests** and local
  Markdown link-target checks passed. Modified public signatures are unchanged.
- Fresh isolated **0.3.3** wheel/sdist builds and both Twine checks passed. The
  distribution verifier passed separate fresh dependency-resolution environments
  and isolated sdist build dependencies: each installed distribution imported all
  106 modules and executed both README examples outside the checkout, with typing,
  licensing, source provenance, metadata and version checks.
- Wheel SHA-256: `61c246fe0f110ad0da2b778ccc02d32aefe767063139340f5f7f43c8c7aa8cf7`.
  Sdist SHA-256: `2aa26b33e291bbcd67bc5a09a44b1369095bcb29b8a66181f97f40e6a3eaf679`.
  Artifacts and the verifier report are under ignored
  `dist/maintenance-audit/artifacts-final/` and
  `dist/maintenance-audit/distribution-report.json`, respectively. Baseline,
  reproduction, unsuccessful-attempt and final test logs are preserved alongside.
- Runtime dependencies: NumPy 2.5.4, SciPy 1.18.1, scikit-learn 1.9.1, and
  Numba 0.68.1. Results apply to Windows/Python 3.12.14 and the explicit source
  snapshots, not unexecuted supported Python/platform combinations.

Fresh Linux symlink checks remain pending:
Docker did not respond within 12 seconds, and the only WSL distribution is
`docker-desktop`. Previous Linux results are historical. No push, tag, merge or
publication is authorized by this maintenance audit; remote CI must validate the
final pushed revision separately.

## Second audit continuation 0.3.3 — 2026-11-02

### Preserved baseline and documentation defects

- Continued the same unpublished changeset on `codex/maintenance-audit`; version
  remains 0.3.3. The complete tracked/nonignored-untracked starting tree is frozen
  under ignored `dist/second-audit/baseline/`, with a raw source manifest, recovery
  patch and dependency inventory alongside. Starting tree SHA-256:
  `6d219e3e0345b7f71f79040404d20ea83bc1cb09e639e4eaf4087e027d83ba6a`.
- The pre-implementation recheck passed 1,161 functional and nine slow tests,
  with four Windows symlink privilege skips. Static checks and distribution-source
  comparisons passed. These counts describe the starting tree, not final results.
- Added [bounded documentation verification](scripts/check_documentation.py)
  and a failing end-to-end regression before fixes. It reproduced **10 malformed
  FASTA docstrings**. After escaping repairs, execution reproduced **38 incorrect
  example outputs** across 208 examples; reports retain both discovery stages.
- Corrected examples against direct windows, literal/IUPAC expansions, known
  frequencies, documented coordinate conventions and existing independent
  scientific fixtures. Examples now reflect gap-excluding identity, thresholded
  consensus, leftmost complete minimizer windows, accepted protein X codes,
  population-approximation ENC and sorted logo alphabets. Runtime conventions
  and public signatures were preserved rather than changed to fit examples.
- Replaced missing-file examples with temporary inputs, stochastic outputs with
  invariant assertions, and oversized translation jobs with small executable
  examples. All **209 examples** then passed. The checker first parses every
  docstring, then executes examples in an isolated child process and temporary
  working directory, with a finite positive timeout and one numerical thread.
  Unit tests verify both real examples and rejection of malformed escaping.

### Independent correctness and performance review

- Replaced permissive masking assertions with exact expected strings. Added
  exhaustive short ACN entropy/interval-union fixtures, case and fifth-symbol N
  cases, nonfinite parameters even on empty input, all short DNA IUPAC expansion
  matches, direct canonical/noncanonical minimizer windows, and traceback ties.
  Existing independent fixtures cover all 12 module families, including NCBI
  codes, serial/process equivalence, annotation conversions, scientific model
  limits, malformed records and transactional file failure/cleanup scenarios.
- An AST comparison against the frozen starting package confirmed that all
  executable package code was unchanged by this continuation's documentation
  corrections (46 files changed only in docstrings). No new runtime defect was
  confirmed by these independent checks; no public API or dependency was added.
- Added [calibrated alignment verification](benchmarks/alignment_recheck.py):
  five counterbalanced rounds, batches of at least 100 ms after separate warmup
  and calibration, square/unequal 100x50, 300x300, 1000x100 and 1000x1000 inputs,
  random/repetitive/tie-heavy cases, complete output digests, frozen measuring
  code/source manifests, and a byte-identical edit-distance control. Python
  allocations, native process peak RSS and import startup are measured separately.
  JIT and process-pool timings are outside this alignment-only experiment.
  Report fields `rows`/`columns` are sequence lengths; the global-alignment
  dynamic-programming matrices include boundary cells and are respectively
  **101x51, 301x301, 1001x101 and 1001x1001**.
- Control drift above 10% makes the run inconclusive. Repeat at most once with
  no concurrent heavy checks and retain both reports; persistent instability
  cannot support a new speed claim. Final timing/retention outcomes follow below.
- [First calibrated run](benchmarks/second_audit_alignment_results.json) and
  [single repeat](benchmarks/second_audit_alignment_repeat_results.json) each
  produced **24 equivalent complete outputs**. Every reported timing batch was
  at least 100 ms. Unchanged control ratios were 0.818–1.190 initially and
  0.636–1.239 on repetition: **both timing experiments are inconclusive**.
  Raw global-alignment gates passed 12/12 and 11/12 cases respectively, but
  uncontrolled timing does not establish a new speed improvement or regression.
  No further timing repetitions or new optimization candidates were introduced.
  The previously retained global implementation is unchanged; its earlier
  retention evidence remains historical rather than a new quantitative claim.
- [Current source manifest](benchmarks/second_audit_source_manifest.json)
  records the 46 docstring-only file changes, runtime-AST equivalence to the
  starting working tree, and exact final package SHA-256
  `9ea1cfe38c7468c63b28b0bd5f2231f9d53cfd1c33b3ab41e6dd97d5f7edac71`.
  Both calibrated runs measured this exact package, including documentation.

### Final evidence

- The completed functional coverage run passed **1,177 tests**, with **four
  Windows symlink privilege skips**, nine slow tests deselected, and **96.74%**
  statement coverage (2,236 statements, 73 missed). The sole warning is the
  expected duplicate-motif clustering warning. Logs, JUnit and coverage JSON
  are preserved in ignored `dist/second-audit/`.
- On 2026-11-03, the separate slow run passed **all nine tests** (1,181
  deselected), and fresh executable documentation verification passed **209
  examples across 106 modules**. Per-file package hashes still match the source
  manifest above. All 20 local Markdown link targets in the root documentation
  were checked successfully.
- Validation attempts are distinguished from completed passes. An initial
  functional attempt stalled and was interrupted; the completed diagnostic
  functional run above passed. A subsequent fresh-environment functional attempt
  ended without a final result while waiting on the documentation subprocess;
  it is not counted as a pass. A separate earlier slow-test attempt encountered
  a Windows native access violation during NumPy import before tests ran.
  These logs are retained rather than overwritten or attributed to a confirmed
  package defect. The successful slow/documentation retries use a fresh
  `dist/second-audit/check-env/` with Python **3.12.14**, editable development
  dependencies, and numerical-library thread limits of one. The original `.venv`
  is preserved. Both dependency inventories are recorded; scientific versions
  remain NumPy 2.5.4, SciPy 1.18.1, scikit-learn 1.9.1 and Numba 0.68.1.

- Fresh-environment Ruff lint passed, formatting passed for **222 Python files**,
  strict mypy passed for **106 source files**, and `pip check` reported no broken
  requirements. The isolated wheel and sdist builds succeeded with setuptools
  84.2.1. Both artifacts passed Twine checks. Final artifacts are in the fresh
  ignored directory `dist/second-audit/artifacts/`:
  wheel SHA-256 `321daadd24e46b4b0add8c9b957b7040d7afd2e667e021d48a85c0a9f94b428f`;
  sdist SHA-256 `4b2022828288a7986aad6406dfd65df253cae60ff9246143e237bedd7e4ec3ed`.
- The final distribution verifier **passed** with fresh dependency resolution,
  separate wheel/sdist virtual environments and isolated sdist build dependencies.
  Each installation imported **106 modules** and executed **both README
  examples outside the checkout**. Metadata, version, public exports, `py.typed`,
  licensing and package-source provenance checks passed. The detailed report is
  preserved at `dist/second-audit/distribution-report.json`; its artifact hashes
  match the files above. This completes local validation for this continuation;
  unsuccessful attempts remain recorded separately.

The bounded Linux retry timed out waiting for Docker after eight seconds; WSL
listed only `docker-desktop`. Linux symlink/file-safety execution and the remote
CI platform matrix remain **pending**. Windows privilege skips do not certify
Linux behavior. No push, merge, tag or publication is part of this continuation.

Continuation recheck on 2026-11-03 again timed out after eight seconds for Docker;
WSL still listed only `docker-desktop`. Package-source and verified artifact
hashes are unchanged, and `git diff --check` passed. The bounded probe is retained
in ignored `dist/second-audit/platform-followup.json`. No platform test pass is
inferred from this availability check; the pending gates above remain open.

The remaining fresh-environment functional reproducibility check completed on
2026-11-03: **1,177 passed, four Windows symlink skips, nine slow tests deselected,
one expected clustering warning**, in 83.49 seconds. Statement coverage again
matched **96.74%** (2,236 statements, 73 missed). This successful rerun includes
the documentation subprocess and replaces no earlier failed/incomplete log;
`fresh-functional.log`, `fresh-functional.xml` and
`fresh-functional-coverage.json` are retained under `dist/second-audit/`.
Both the original and fresh environments have now completed the functional
suite against the same recorded package-source hashes.

## Main synchronization and audit continuation (2026-11-03)

### Preservation, baseline and fixes

- Backed up and SHA-256-verified all 244 tracked/nonignored untracked files in
  `D:\Projects\bioutils-audit-backup-20261103-141033` before stashing. Restored
  every file after synchronization; five had only Git line-ending conversions.
  The backup and named preservation stash remain available.
- Fetched origin, switched to `main`, pulled with `--ff-only`, returned to
  `codex/maintenance-audit` and merged `main`. Both operations were already up
  to date. Local `main` and the audit HEAD remain
  `186a85fbbd70ac31cb3b013281cb7691d4365c59`; remote `main` remains
  `9e1cc06420f14940f9e640af554539be29886fb5`. The local commit was preserved.
- Used the explicit Python 3.12.14 interpreter in
  `dist/second-audit/check-env/`, refreshed its editable development installation,
  and recorded dependencies under `dist/sync-audit/`. The default shell Python
  3.8 is unsupported and was not used. Numerical-library thread limits were one.
- Starting-tree functional coverage: **1,177 passed, four Windows symlink
  privilege skips, nine slow tests deselected**, one expected clustering warning.
  This is a fresh baseline run, distinct from the earlier historical results.
- Reviewed all 12 module families and the existing independent alignment,
  genetic-code, coordinate, selector, biological-statistics and file-safety
  fixtures. New regressions reproduced **six failures** before fixes: valid BED
  records on chromosomes named `track`/`browser` were discarded, and NaN/positive
  infinity bypassed BLAST score validation. Negative infinity already raised.
- Retain keyword-named BED chromosomes when tab-delimited coordinates are
  integers; continue skipping metadata. Require finite BLAST reference/target
  scores while preserving finite ratios and clamping. Focused checks passed
  **137 tests**. No public signatures, exports or coordinate conventions changed.
  Version remains **0.3.3**, the existing bump for this continuing changeset.

### Performance evidence

- [Fresh paired results](benchmarks/sync_audit_paired_results.json) explicitly
  select baseline `186a85fbbd70ac31cb3b013281cb7691d4365c59`, five counterbalanced
  timing rounds, all default workloads at 1k/10k/100k bases, and separate Python
  allocations, process RSS, import/JIT and process-startup measurements.
  **All 103 comparisons have equivalent complete outputs**, including the
  separately labeled corrected-reference cases.
- The existing global-alignment row-reuse implementation passes the numerical
  retention gate in all six comparisons; medium/large runtime ratios are
  **0.601–0.819**. Quadratic alignment inputs are capped: these medium/large
  cases use 100x50 and 300x150 sequence lengths, not 10k/100k alignment matrices.
- Unchanged controls show substantial drift: local-alignment ratios span
  0.476–1.664 and other controls also vary. Raw gate labels are not proof of
  regressions in unchanged code or a reliable new percentage speed claim.
  Retain the existing optimization with its historical evidence and this raw
  gate confirmation; introduce no additional optimization candidate.
- [Final package manifest](benchmarks/sync_audit_source_manifest.json) matches
  the measured package SHA-256
  `8e349fcdc744127fe5cf204c10cc3763b2b725785ae85bce4b4dee5727f2c931`.
  Full benchmark manifests describe measurement-time state; subsequent test
  formatting and documentation/evidence edits do not change package sources.

### Current validation

- Final Windows/Python 3.12.14 functional run: **1,188 passed, four symlink
  privilege skips, nine slow tests deselected**, one expected clustering warning;
  **96.79%** coverage (2,242 statements, 72 missed). The separate slow run passed
  **all nine tests**. Logs, XML and coverage JSON are under `dist/sync-audit/`.
- Docker became available during this continuation. On Linux/Python **3.12.16**,
  **41 targeted file-safety and new regression tests passed without skips**,
  including symlink alias rejection and output-symlink preservation. Tests ran
  on the container's `/tmp` filesystem after copying from a read-only checkout
  mount and installing editable development dependencies. Image `python:3.12-slim`
  ID: `sha256:aa943ce411423c7582288f18707ff66eddcbaeb510502d747a72aeb34f30d159`.
  This closes the previously pending targeted Linux symlink check; it does not
  establish a full Linux suite or the remote interpreter matrix.
- Fresh executable documentation verification passed **209 examples** across
  106 modules; all **22** final local root-documentation links resolved. Package modules
  also parse under Python 3.10 syntax rules, without claiming 3.10 execution.
- Ruff lint, formatting (**223 files**), strict mypy (**106 source files**) and
  dependency consistency passed. Corrected the new test's initial formatting
  failure before the final formatting pass. Reproduction failures remain in
  `dist/sync-audit/reproduced.log`; validation details are recorded in
  `dist/sync-audit/validation-summary.json`.
- Fresh isolated **0.3.3** wheel/sdist builds and both Twine metadata checks
  passed. Final artifacts are in `dist/sync-audit/artifacts/`, containing only
  the wheel and sdist. Wheel SHA-256:
  `e37f3e96fb86566ad3688885207233db8d2e5e5341432e371b7949c5d8f08168`;
  sdist SHA-256:
  `b2dc9ad4ffb250aa1a38918c697025cd5e19b4b8729d95c95fb352d42967a195`.
- The distribution verifier **passed** with separate fresh wheel/sdist virtual
  environments, fresh dependency resolution and isolated sdist build dependencies.
  Each installation imported **106 modules** and executed **both README examples
  outside the checkout**. Version, metadata, source inventory/provenance,
  public exports, typing marker and license checks passed. Dependency versions
  remained NumPy 2.5.4, SciPy 1.18.1, scikit-learn 1.9.1 and Numba 0.68.1.
  Report: `dist/sync-audit/distribution-report.json`; reported artifact hashes
  match the files above. This audit document is not included in the sdist, so
  recording these results does not change the verified artifact contents.
- Final preservation verification found no missing original files and changes
  to only the five intended existing files in this continuation (the two fixes,
  README, changelog and this evidence); existing audit edits remain intact.
  Remote CI remains pending on the final pushed revision. No commit, push,
  remote merge, tag or publication was performed.

## Documentation and branch verification (2026-11-08)

- Refreshed origin successfully using Git's Windows Schannel TLS backend after
  the default OpenSSL connection was reset. Remote `main` is still `9e1cc06`
  (version 0.3.1); local `main` and `codex/maintenance-audit` point to `186a85f`
  (version 0.3.2). Local `main` is one commit ahead of remote `main`.
- The continuing 0.3.3 changeset remains uncommitted: 59 tracked files modified
  and 15 untracked files. These changes are not committed on either main branch.
  No branch switch, commit, push or merge was performed during this check.
- README now explicitly identifies the unpublished checkout as 0.3.3 and links
  its release notes, distinguishing it from live published-package badges.
  Package source, installed distribution metadata and the latest changelog
  section agree on 0.3.3. Keep the existing bump for this continuing changeset.
- Fresh Windows/Python 3.12.14 verification: both documentation tests passed
  (including executable package docstrings), all 23 local root-documentation
  links resolved, dependency consistency and `git diff --check` passed.
- Fresh isolated wheel/sdist builds and both Twine checks passed. Artifacts are
  under `dist/version-check-20261108/artifacts/`. An initial nonisolated build
  failed because its interpreter lacked setuptools; the isolated build passed.
  Earlier functional/performance results above retain their original dates;
  this documentation check does not establish a new full-suite or remote CI run.
- Fresh distribution verification passed for both artifacts with separate fresh
  virtual environments, fresh dependency resolution and isolated sdist build
  dependencies. Each installation imported all 106 modules and executed both
  README examples outside the checkout; source and metadata checks passed for
  0.3.3. Report: `dist/version-check-20261108/distribution-report.json`.

## Version 0.3.3 release validation (2026-11-08)

- The user authorized updating GitHub `main` and publishing 0.3.3 to PyPI,
  with validation first. Preserve the existing version bump for this changeset.
  README identifies 0.3.3 and the release notes use the release preparation date.
- AST comparison against `186a85f` confirms unchanged function signatures in
  all modified package files. Seven functions change implementation; the other
  package changes correct docstrings or the version. Independent fixtures cover
  alignment ties, format whitespace, GC denominators, invalid/non-finite values,
  chromosome metadata keywords, entropy masking and selector/IUPAC results.
- Fresh Windows/Python 3.12.14 functional coverage run: **1,188 passed, four
  Windows symlink privilege skips, nine slow tests deselected**, with one expected
  identical-motif clustering warning; **96.79%** coverage (2,242 statements,
  72 missed). Separate slow run: **nine passed, 1,192 deselected**.
  Logs: `dist/publish-functional-serial-20261108.log` and
  `dist/publish-slow-serial-20261108.log`.
- Ruff lint, formatting (223 files), strict mypy (106 source files), dependency
  consistency and `git diff --check` passed. Initial concurrent checks encountered
  host memory pressure; interrupted functional/slow/build attempts are not passes.
  The complete functional and slow results above were obtained sequentially.
- Remote platform validation must run against the final pushed commit. Earlier
  Linux symlink evidence above remains historical, rather than a new local run.
- Final artifact verification completed on **2026-11-09** after continuation:
  isolated wheel/sdist builds and both Twine checks passed. Separate fresh
  installations with fresh dependency resolution and isolated sdist build
  dependencies each imported 106 modules and executed both README examples
  outside the checkout. Source inventory, version 0.3.3, exports, license and
  typing-marker checks passed. Report:
  `dist/publish-final-20261108/distribution-report.json`. The directory retains
  the preparation date; the verification completed the following day.
