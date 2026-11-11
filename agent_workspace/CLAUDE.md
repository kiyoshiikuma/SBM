# CLAUDE.md

Development guide for Claude. Shared procedures are maintained in [Agent.md](Agent.md); detailed physical models and validation criteria are in [Handbook.md](Handbook.md). Read both before starting work.

Code paths and commands are relative to the repository root. Maintain the shared guides in `agent_workspace/` in version control. Keep `.venv/` and machine-specific artifacts local.

When working on a new systematic effect, follow this sequence:

1. Use the [specification template](Handbook.md#spec-template) to create an effect specification in `agent_workspace/`. Define the TOD, spin components, angle conventions, approximation order, public API, and independent validation first.
2. Return only the perturbation for ordinary effects and preserve the spin information carried by `Sky` and `SpinMap`. Keep registration, exports, and module docstring lists consistent; leave helper functions out of the public list.
3. Check agreement with real-space TOD at the same approximation order using the default tolerances. Test approximation truncation error against the exact signal separately; do not absorb it into the TOD agreement tolerance. Add absorption cases when templates are declared.
4. **Only when implementing or changing beams**, check the optical chain in the [beam-specific requirements](Handbook.md#beam-spec). A beam model is the full signal that replaces `sky_signal`.
5. Validate with `.venv/bin/python`. Use Ruff `0.17.4` and distinguish new diagnostics by comparing with the [baseline](Handbook.md#baseline). Add the validation results and limitations to the specification.

If conventions disagree, check [conventions.rst](../docs/source/conventions.rst) and the code, then update the Handbook as needed. Follow Agent.md for the language of code, docstrings, and upstream documentation, and for completion criteria. Avoid duplicating equations, detailed validation criteria, or baseline numbers in this file.

Maintain this guide in `agent_workspace/CLAUDE.md`, alongside `Agent.md` and `Handbook.md`.
