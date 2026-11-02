---
name: major-check
description: Audit bioutils-collection correctness, scientific conventions, file safety, packaging, portability, and measured performance. Use for major maintenance and broad optimization work; focused fixes need their relevant checks.
---

# Major check

Read root `AGENTS.md` and `MAJOR_CHECK_PLAN.md` first.

1. Record Git state, baseline revision, Python/platform, dependency versions, and
   actual baseline results. Preserve unrelated edits; synchronize Git only when
   requested. Never reset existing work.
2. Check all module families. Prioritize file preservation, annotation coordinate
   conventions, independent biological fixtures, malformed records, and backend
   or serial/process equivalence. Add regressions before fixing confirmed defects.
3. Profile before optimizing. Use `benchmarks/paired_audit.py` for immutable
   baseline/working-tree snapshots, source manifests, five counterbalanced timing
   rounds, separate Python allocation/process RSS checks, JIT initialization,
   warmed execution, and process startup. Use `benchmarks/major_audit.py` for
   additional corrected-reference cases. Distinguish references from old code;
   apply root AGENTS.md's balanced retention gate on equivalent medium/large cases.
4. Run full functional coverage and slow tests separately, then Ruff lint/format,
   strict mypy, dependency consistency, wheel/sdist builds, metadata checks, and
   installed-package/README smoke checks outside the checkout. Run
   `scripts/check_distribution.py --artifact-dir <fresh-artifact-directory>` for
   separate wheel/sdist virtual environments. Execute symlink checks on Linux;
   Windows privilege skips must not stand in for Linux evidence.
5. Always update documentation based on the changes and update the package version
   once per completed changeset as required by root `AGENTS.md`. Use the
   documentation-maintenance skill for synchronization and final-version checks.
6. Record current results, warnings, skips, model limitations, benchmark evidence,
   and pending remote platform checks. Earlier green CI cannot certify local edits.

Do not push, merge, tag, or publish unless requested.
