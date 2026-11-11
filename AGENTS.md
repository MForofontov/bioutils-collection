# Repository guidance

## Project and compatibility

- Package source is `bioutils_collection/`; unit tests are in `pytest/unit/`.
- `pyproject.toml` defines supported Python versions, dependencies, and tooling.
  `bioutils_collection/_version.py` is the version source of truth.
- Preserve public signatures, exports, return shapes, ordering, and coordinate
  conventions unless a documented bug correction requires a behavior change.
- Reproduce defects before fixing them. Validate scientific results against
  independent fixtures; passing the existing suite alone is insufficient.

## Documentation and version are part of every change

- Always update documentation based on the actual changes: affected docstrings,
  README examples, setup instructions, maintenance evidence, and release notes.
  Every completed changeset must include an accurate changelog entry; document
  internal maintenance even when no public usage changes. Do not make cosmetic
  edits solely to satisfy this requirement.
- Always update the package version once per completed changeset, including
  documentation/tooling changes. Check the baseline and current version first;
  do not bump again for individual files, retries, or continuations of the same
  changeset. Use a patch increment for fixes/docs/tooling, a minor increment for
  new public capabilities, and an appropriate breaking-change increment.
- Keep `_version.py`, package metadata, release notes, and documented version
  references consistent. Preserve historical release sections. A version bump
  does not authorize a tag, push, merge, or publication.
- Tie test counts and performance claims to the actual revision/environment.
  Distinguish current passes, historical results, skips, and pending remote CI.

## Verification

Install with `python -m pip install -e ".[dev]"`, selecting the interpreter
explicitly. Limit BLAS/OpenMP threads to one on memory-constrained hosts.

```bash
python -m pytest -q -m "not slow" --strict-markers --cov=bioutils_collection
python -m pytest -q -m slow --strict-markers
python -m ruff check .
python -m ruff format --check .
python -m mypy bioutils_collection
python -m pip check
python -m build
python -m twine check dist/*
python scripts/check_distribution.py --artifact-dir dist --report dist/distribution-report.json
```

Run focused regressions during implementation. Build and smoke-test the final
version outside the checkout. Remote CI must match the final pushed revision.
Do not weaken tests or claim unrun platforms passed.
Use a fresh artifact directory containing only the final wheel and sdist. The
distribution verifier installs them into separate fresh virtual environments.
Its explicit dependency-reuse option is for local offline checks, not proof of
clean dependency resolution.

For paired optimization evidence, run `python benchmarks/paired_audit.py --output
benchmarks/paired_audit_results.json`. It snapshots tracked and nonignored
untracked files and uses five counterbalanced timing rounds. Retain an optimization
only with equivalent outputs and at least 10% median runtime improvement, or 25%
Python allocation savings with no more than 10% runtime cost on representative
medium/large workloads. Report process RSS separately; never call incorrect old
outputs an equivalent baseline. Slow unit tests assert correctness, not speed.

## Shared workflows

- Broad audits: [.agents/skills/major-check/SKILL.md](.agents/skills/major-check/SKILL.md).
- Documentation/version maintenance:
  [.agents/skills/documentation-maintenance/SKILL.md](.agents/skills/documentation-maintenance/SKILL.md).
- Audit evidence and model limits: [MAJOR_CHECK_PLAN.md](MAJOR_CHECK_PLAN.md).
