---
name: documentation-maintenance
description: Keep bioutils-collection documentation, examples, release notes, package version, and repository instructions synchronized with implementation changes. Use for documentation/version maintenance or contributor workflow changes.
---

# Documentation and version maintenance

Read root `AGENTS.md`, `pyproject.toml`, and `bioutils_collection/_version.py`.

- Always update documentation based on the actual diff: affected public
  docstrings, usage examples, setup commands, scientific assumptions, and audit
  evidence. Every completed changeset needs an accurate changelog entry.
- Always update the package version once per completed changeset under the root
  versioning policy, including docs/tooling changes. Compare with the baseline
  first; an existing bump for the same changeset must not be incremented again.
- Synchronize metadata and version references. Preserve published history; record
  the new version and actual date without implying that it has been published.
  Version maintenance does not authorize tagging or publication.
- Execute changed examples using temporary destinations. Verify local Markdown
  links, public imports, and installed distribution metadata for the final version.
  Use `scripts/check_distribution.py --artifact-dir <artifact-directory>` to
  reject stale artifacts and verify independent wheel/sdist installations.
- Keep `AGENTS.md` canonical and tool rules as short routing instructions.
- Record exact validation outcomes, scientific model limits, skips, and unrun
  platforms. Do not copy historical results as fresh evidence.
- Label paired implementation benchmarks separately from corrected-reference
  benchmarks. Tie each to snapshot hashes and state both timing and memory costs;
  document rejected optimizations as well as accepted ones.
