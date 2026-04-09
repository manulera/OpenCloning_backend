# Follow-ups

## uv workspace: dependency isolation guardrails (not implemented yet)

Python and uv do not enforce per-package import isolation in a workspace: all members share one resolution and one environment, so a package can import a dependency that was only pulled in by another member or that was never declared in its own `pyproject.toml`. Uv documents this limitation under [When (not) to use workspaces](https://docs.astral.sh/uv/concepts/projects/workspaces/#when-not-to-use-workspaces).

There is no uv setting that turns on stricter checks; mitigations are process and tooling.

### Possible mitigations (pick what fits the repo later)

| Approach | What it helps with |
| -------- | ------------------ |
| [deptry](https://github.com/fpgmaas/deptry) or [FawltyDeps](https://github.com/tweag/FawltyDeps) in CI | Compare imports under each member’s `src` to that member’s declared dependencies (including “using another member’s transitive dependency” once multiple members exist). |
| `uv build` then install the wheel in a clean virtualenv and run smoke tests | Validates that the **published artifact** runs with only declared runtime dependencies. |
| [import-linter](https://github.com/seddonym/import-linter) (especially with 2+ members) | Enforce forbidden imports between packages or layers (e.g. CLI must not depend on another package’s internals). |
| Per-member CI commands (`uv run --package <name> pytest …`) | Keeps test entrypoints and selection aligned with a member; does not replace undeclared-import detection. |
| Team convention + review | Prefer importing only your package’s public API; extra checklist when adding new workspace members. |

### Single-member workspace

Until a second member exists, the main risk is the usual one: code imports a third-party name that “happens” to be installed transitively but is not listed in `opencloning`’s dependencies. **deptry** (or a wheel-in-clean-venv check) is the highest-value automated guardrail to add before or when the workspace grows.
