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
