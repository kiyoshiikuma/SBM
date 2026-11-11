# Agent.md

Shared guidelines for agents working on SBM. Detailed specifications and validation criteria are maintained in [Handbook.md](Handbook.md). Code paths and commands in this document are relative to the repository root.

## Read first

- [Handbook.md](Handbook.md): local setup, implementation specifications, and the two layers of testing.
- [Specification template](Handbook.md#spec-template): design and acceptance criteria for a new effect.
- [CONTRIBUTING.md](../CONTRIBUTING.md): the upstream contribution workflow.
- [conventions.rst](../docs/source/conventions.rst): angles, spins, units, and conventions for external data.
- [adding_systematics.rst](../docs/source/adding_systematics.rst): registration and testing procedures.

If descriptions of physical conventions or APIs disagree, check `conventions.rst` and the current code, identify the discrepancy, and update the Handbook. Respect the user's edits and change only what is needed.

## Workspace files and design

- Maintain the shared guides in `agent_workspace/` in version control. Keep `.venv/` and machine-specific artifacts local; use `.git/info/exclude` for local exclusions. See Handbook Section 1 for setup instructions.
- For each new systematic effect, use the template to save a specification in `agent_workspace/`. Define the ideal TOD, TOD with the error, perturbation, spin decomposition, units and signs, approximation order, omitted cross terms, and independent validation before implementation.
- Record the reference equations, implementation module, public API, templates, and components that cannot be absorbed. Adding effect models does not automatically generate cross terms between them.

## Implementation guidelines

- Ordinary systematic functions return only the perturbation as a `SignalModel`. Keep scan, `HMaps`, and map-maker processing outside these functions.
- Add the function to the appropriate module in `sbm/systematics/`. Update `@systematic(data=..., templates=...)`, the public imports and `__all__` in `sbm/systematics/__init__.py`, and the function list in that module's docstring.
- Start helper function names with `_` and leave them out of `__all__`. The set of public callables, excluding `systematic`, must match `ALL`. `ALL` is a registry, not a mechanism for applying effects in sequence.
- Follow `P = Q + iU`: the polarization component is `(2, 0)` without a HWP and `(2, -4)` with a HWP, with coefficient `P/2`. Do not register both members of a conjugate pair.
- Build components from the spin-carrying maps of `Sky`, with spin weight equal to `n`. For external maps, check the transformation convention and make the spin weight explicit with `SpinMap(arr, s)`. Do not relabel a spin weight merely to pass a check.
- Avoid stripping spin information through `np.asarray` or by separating Q and U. Accepting plain ndarrays does not imply that a registered function passes the spin checks. The current `UNCHECKED` entries are `hwp_mueller` and `beam`. Add an exemption only with a reason and independent validation.
- End angle argument names with `_rad`. Document units, signs, data type, required maps and derivatives, approximation order, templates, and handling of invalid inputs in the docstring.
- Do not apply the external convention conversion performed by `HMaps.load/save` twice. Use `HMaps.load_legacy`, which has a different conversion, for SBM 0.4 files.
- **Only when adding or changing beam implementations**, apply the [beam-specific requirements in Handbook.md](Handbook.md#beam-spec). A beam model returns the full signal; do not add `sky_signal` again.
- Write code, comments, docstrings, documentation in `docs/source/`, and these guides and the Handbook in English.

## Separate validation into two layers

1. **Agreement at the same approximation order**: Add a case to `cases()` in `tests/test_systematics.py`, with the full model including the baseline signal and real-space TOD written independently of the spin decomposition. Truncate the reference TOD to the same order as an approximate model and retain the `assert_matches_tod` defaults, `rtol=atol=1e-9`. Here, `atol` is scaled by the maximum absolute TOD value. Check the main argument branches and the zero-error limit.
2. **Approximation truncation error**: For an approximate effect, compare with the exact signal in a separate test. Check the residual hierarchy as the retained order increases, or the convergence order as the error amplitude is varied. Do not loosen the tolerances in layer 1 to accommodate the approximation.

`RealSky.detector` uses SBM angles; `lbs_polarizer_tod` and `lbs_mueller_tod` use litebird_sim angles. State the reference TOD's angle convention in the specification.

For an effect that declares templates, also register `absorbed=lambda: (base, full)` using parameters whose perturbation can be absorbed exactly. Follow Handbook Section 7 for the absorption test's Basis and comparison criteria. Consider adding new templates to `HWP_BASES` or `NO_HWP_BASES` if they should be exercised against all existing cases. See the beam-specific requirements for additional tests when changing beams.

## Run and verify

Use the Python interpreter in `.venv` from the repository root. See Handbook Section 1 for development dependencies and the [baseline](Handbook.md#baseline) for existing results.

```sh
.venv/bin/python -m pytest tests -q
.venv/bin/python -m ruff --version
.venv/bin/python -m ruff check sbm tests --statistics
.venv/bin/python -m sphinx -b html -W docs/source docs/build/html
```

For implementation changes, run the relevant tests first, then the full suite. Use Ruff `0.17.4`, compare total and per-rule diagnostic counts with the baseline, and check changed code for new diagnostics. An unchanged count alone is insufficient. Distinguish a nonzero exit code caused by existing diagnostics from new problems in the report.

For documentation-only changes, check links, references, and consistency with the code. Also run Sphinx when changing `docs/source/`.

## Completion criteria

- The specification, implementation, and independent TOD tests describe the same physical model.
- Validate the zero-error limit, main argument branches, template absorption, and truncation error where an approximation is used.
- When adding a systematic, update its public registration, module docstring, `docs/source/systematics.rst`, `docs/source/reference.rst`, and `CHANGELOG.md`.
- Record the checks performed, their results, differences from the baseline, and untested scope in the effect's specification and the completion report.

Do not make bulk changes to unrelated equations, existing data, or notebook outputs. Keep generated artifacts, virtual environments, and large outputs out of commits.
