# Handbook.md — Implementation Specifications for New Systematic Effects

This handbook defines how to specify a new systematic effect in SBM, from its physical model to its public API and tests. Copy the [specification template](#spec-template) for each effect, fill it in, and maintain it alongside the implementation.

It applies to yusuke-takase/SBM version 0.4.2 at the installation commit 780b047da592ffe3f4e65e0275d925668d9ee65b. The tables in Sections 2–4 summarize the upstream documentation and code at that revision. If they disagree, treat [conventions.rst](../docs/source/conventions.rst) and the code as authoritative, and update this handbook to match the implementation.

## 1. Setup and Execution

Keep this document in agent_workspace/. Unless stated otherwise, code paths and commands in this document are relative to the repository root.

The shared guides in agent_workspace/ are versioned repository documentation. Keep .venv/ and machine-specific artifacts out of commits. Python 3.12's venv does not create its own .gitignore, so add the virtual environment to the local exclusion file once per clone.

```sh
printf '%s\n' '/.venv/' >> .git/info/exclude
git check-ignore -v .venv   # Should point to .git/info/exclude
```

SBM requires Python 3.10 or later. The following commands use Python 3.12 to create .venv at the repository root and install SBM in editable mode. If the system Python is older than 3.10, select a supported interpreter explicitly.

To create the environment, run the following from the repository root.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -e '.[test,plot,docs]' 'ruff==0.17.4'
.venv/bin/python -m pip check
.venv/bin/python -c 'import sbm; print(sbm.__version__)'
```

Replace the executable in the first command with another Python 3.10 or later interpreter if needed. Use `pip install -e .` for normal use, the test extra for tests, plot for plotting, and docs for Sphinx. Ruff is not included in these extras; install it separately with its version pinned. pyproject.toml does not select Ruff rules, so diagnostic counts depend on each version's default rules (see the [baseline](#baseline)).

```sh
source .venv/bin/activate
python -m pytest tests -q
```

The main dependencies are NumPy, healpy, h5py, and ducc0. Treat pyproject.toml as authoritative for dependency requirements.

## 2. Model Responsibilities and Implementation Locations

A normal systematic-effect function returns only the perturbation relative to an ideal detector, `delta d = d_with_error - d_ideal`, as a `SignalModel`. Keep it independent of the observing scan and map-making.

```python
from sbm import Basis, map_make
from sbm import systematics as sy

model = sy.sky_signal(sky, hwp=True) + sy.abs_pointing(
    sky, rho_rad=1e-4, chi_rad=0.3
)
basis = Basis.iqu(hwp=True).extend(sy.abs_pointing.templates)
result = map_make((hmaps, model), basis)
```

The caller prepares `sky` and `hmaps`. When adding several effects, specify the reference model assumed by each effect and the omitted cross terms. Addition alone does not generate cross terms automatically. For nonlinear responses, use products of `SignalModel` objects as needed, for example `d + beta * (d * d)`.

The beam is the exception: it returns the complete signal. Its usage is defined in the [beam-specific specifications](#beam-spec).

| Implementation location | Main purpose |
| --- | --- |
| sbm/systematics/sky.py | Ideal signal, `sky_signal` |
| sbm/systematics/pointing.py | Pointing offsets, HWP wedge, displacement expansions |
| sbm/systematics/optics.py | Public beam model |
| sbm/systematics/hwp.py | Non-ideal HWP Mueller matrices and instrumental polarization |
| sbm/systematics/detector.py | Gain, bandpass, and polarization-angle errors |
| sbm/beam.py, sbm/harmonics.py | Beam convolution and spherical-harmonic coefficients |

If an existing effect can express the model with an additional argument, avoid implementing the same model twice. Write upstream code and documentation in docs/source/ in English, following [CONTRIBUTING.md](../CONTRIBUTING.md).

## 3. Mathematical and Data Contracts

### 3.1 Angles and Spin Components

```text
d(psi, phi) = sum S[n,m] exp(i (n psi + m phi))
S[-n,-m] = conj(S[n,m])
h[n,m] = <exp(-i (n psi + m phi))>
P = Q + iU
```

`psi` is the crossing angle. `phi` is the HWP angle relative to the telescope, counted in the opposite direction to the crossing angle.

| Data | Ideal-signal components |
| --- | --- |
| Single detector, without HWP | (0, 0): I; (2, 0): P/2 |
| Single detector, with HWP | (0, 0): I; (2, -4): P/2 |
| Pair difference, without HWP | (2, 0): P/2; no intensity component |

Without HWP, the ideal TOD is `I + Q cos(2 psi) - U sin(2 psi)`. A pair difference is defined as `(d_T - d_B)/2` for two orthogonal detectors at the same position.

The relations to litebird_sim / Falcons.jl are `psi_SBM = -psi_external`, `phi_SBM = chi_external`, and `h_SBM[n,m] = h_external[n,-m]`. `HMaps.load` / `HMaps.save` convert the h-map file convention; do not apply the conversion again after loading. Read SBM 0.4 cross-link files with `HMaps.load_legacy`, which uses a different conversion: the n = 0 couples are the conjugates of the stored (0, -m) values.

### 3.2 Arrays and Spin Weights

- Component maps must be one-dimensional arrays of shape `(npix,)`, all with the same npix. Do not pass a scalar, such as an angle parameter, as a component map. Components with incompatible shapes raise an exception.
- The `(0, 0)` component must be real. A nonzero imaginary part raises `ValueError`.
- Register only one member of each conjugate pair. Supplying the same physical component at both `(n, m)` and `(-n, -m)` adds the contributions and counts it twice.
- The spin weight of `S[n,m]` must equal n. The HWP harmonic m does not contribute to the spin weight.
- Build components from `Sky` maps (`SpinMap`) and preserve their spin weights when multiplying by coefficients. `SignalModel` itself accepts plain ndarrays without checking their spin weights, but tests/test_spin.py fails registered functions that build components without spin weights. The only exemptions are `hwp_mueller` and `beam`, listed with reasons in `UNCHECKED` in that file. `sky.q` and `sky.u` also have no spin weight individually. For maps supplied from outside `Sky`, assign their spin weight with `SpinMap(arr, s)`. For example, `bandpass_mismatch` wraps temperature maps with `SpinMap(comp, 0)`.

| Sky map | Spin weight |
| --- | --- |
| `sky.i` | 0 |
| `sky.p` | 2 |
| `sky.eth_i`, `sky.eth_bar_i` | 1, -1 |
| `sky.eth_p`, `sky.eth_bar_p` | 3, 1 |
| `second_order['eth2_i']`, `second_order['lap_i']` | 2, 0 |
| `second_order['eth2_p']`, `second_order['eth_bar2_p']`, `second_order['lap_p']` | 4, 0, 2 |

Here, `second_order` means `sky.second_order`. Complex conjugation reverses the sign of the spin weight. Use covariant derivatives for polarization; do not replace them with derivatives that treat Q and U as independent scalar maps.

### 3.3 Units and Approximations

Angle argument names must end in `_rad`. Specify dimensionless quantities such as gains, temperature maps, derivatives, and beam normalization, and keep the final TOD units consistent. For a small-angle expansion, state the retained orders and omitted terms. NET is in muK sqrt(s), and the sampling rate is in Hz.

See [conventions.rst](../docs/source/conventions.rst) and [spin.py](../sbm/spin.py) for details.

## 4. Public API and Registration

Decorate each function with `@systematic(data=..., templates=...)`.

| data | Meaning |
| --- | --- |
| "pair" | Difference of two orthogonal detectors without HWP |
| "single" | Single detector without HWP |
| "single_hwp" | Single detector with HWP |
| "any" | Several data types, selected by the function arguments |

`templates` lists the couples to add to a Basis to absorb the first-order perturbation. Use it as in `Basis.iqu(hwp=True).extend(effect.templates)`. Also document components that templates cannot separate, such as an error at the same component as the original polarization. If the list cannot be fixed, for example because it depends on arguments, explain why it is empty and which Basis the caller should select.

An existing declaration:

```python
@systematic(data="any")
def polarization_angle(sky: Sky, delta_rad: float, hwp: bool = False) -> SignalModel:
    ...
```

This API is already implemented in detector.py; do not add another function with the same name. A new API's docstring must state its arguments, return value, required Sky maps, data type, error sign convention, perturbation order, and templates. Define and verify the behavior for unsupported argument combinations and invalid input shapes or ranges, including exceptions.

After adding a function, update both the imports and `__all__` in sbm/systematics/\_\_init\_\_.py, and add it to the function list in that module's docstring. The decorator registers a function in `ALL` when its module is imported. Therefore, functions in a module not imported by \_\_init\_\_.py are not registered by an ordinary `import sbm.systematics` alone (importing that module directly does register them). `ALL` is a registry, not a sequence that applies effects in order.

tests/test_systematics.py checks that the set of callables in `__all__` (excluding `systematic`) equals `ALL`. Keep helper functions out of `__all__`, give their names a leading `_`, and place them inside the module; `_offset_components` in pointing.py is an example. Constants (`IDEAL_HWP`, `DATA`) may be included in `__all__`.

<a id="spec-template"></a>

## 5. Specification Template for Each Effect

Copy the following into a design document for each effect. Fill in the physical model and independent validation equations before implementation, then append the results after implementation.

```markdown
# <Effect name> Implementation Specification

## Purpose and Scope
- Physical effect and reference equations or literature:
- Implementation module and public function name:
- Data type (pair / single / single_hwp / any):
- Assumptions and conditions outside the supported scope:

## Parameters
| Name | Meaning, sign, and coordinate frame | Unit | Shape, allowed range, and default |
| --- | --- | --- | --- |
| ... | ... | ... | ... |

## Physical Model
- Ideal TOD d_ideal:
- TOD with the error d_with_error:
- delta d = d_with_error - d_ideal:
- Retained perturbation orders, omitted cross terms, and validity range:
- Expected result at zero error:

## Spin Decomposition
| (n, m) | Expression for S[n,m] | Required Sky maps | Spin-weight check |
| --- | --- | --- | --- |
| ... | ... | ... | Must equal n |

## API and Composition
- Signature, return value, and required derivative maps:
- Composition with sky_signal:
- Templates, absorbed orders, and components that cannot be absorbed:
- Required h-map couples and Basis:
- Behavior for invalid inputs and unsupported conditions:

## Independent Validation and Acceptance Criteria
- Real-space TOD at the same approximation order (tests/test_systematics.py):
  - Equations independent of the spin decomposition and their angle convention (RealSky.detector / litebird_sim angles):
  - Case name, hwp flag, and main argument branches:
  - Absorbed case: parameters whose perturbation is absorbed exactly by the templates:
- Approximation truncation error (a separate test, following test_pointing_geometry.py):
  - Exact reference signal, error magnitudes to sample, and expected residual order:
- Zero-error limit, signs, and normalization:
- Tolerances: test_systematics uses the default rtol = atol = 1e-9. Values used in separate tests and their justification:
- Tests, public exports, docs, and CHANGELOG to update:
- Commands run, results (compare Ruff counts with the baseline), and remaining limitations:
```

For a beam specification, define the return value as the complete signal and include the additional items in the next section.

<a id="beam-spec"></a>

## 6. Beam-Specific Specifications: Optical Chain Review

The optical chain review in this section is an additional requirement only when implementing or modifying a beam.

### 6.1 Response Definitions and Ordering

When designing a beam, document the relative placement of the sky, HWP, optical responses included in the beam, and detector, along with the coordinate frame used to define each. In particular, include the following in the specification.

- Optical responses already included in the input `beam_alm`, and the new effect being added.
- Whether the new effect acts before or after the HWP.
- The order of Mueller matrices and rotations, beam normalization, and detector response.
- Confirmation that the same response is not applied through both `beam_alm` and `mueller`.

The current implementation's pencil-beam limit with HWP matches the following Stokes measurement.

```text
d = b_det^T R(-phi) M R(phi) R(-psi) S
S = (I, Q, U, V)^T
R(a): (Q, U) -> (Q cos(2a) + U sin(2a), -Q sin(2a) + U cos(2a))
```

Matrices act from right to left. M is the Mueller matrix in the HWP's own coordinate frame. Do not arbitrarily reorder matrices or rotations. This equation is a validation reference for the pencil-beam limit; it does not mean that convolution with an arbitrary asymmetric beam can be replaced by this point measurement. The reference is `point_response` in [test_hwp_convolution.py](../tests/test_hwp_convolution.py).

The current API does not automatically construct arbitrary optical chains. If a new layout is needed, determine whether the existing assumptions can represent it, and specify any API changes and independent reference computation required.

### 6.2 Complete-Signal Composition and Input Contracts

```python
from sbm import systematics as sy

model = sy.beam(sky_alm, beam_alm, nside, spin_k=range(5), hwp=True, mueller=M)
```

`sy.beam` returns the complete signal after the beam response, so do not add `sy.sky_signal(...)` to this `model`. If only a beam perturbation is needed, define it as the difference between complete signals from the erroneous and reference beams under the same conditions.

- `sky_alm` and `beam_alm` must be TEB or TEBV coefficients with the same lmax. Prefer `SphericalHarmonics`; litebird_sim coefficients and arrays in healpy order (mmax = lmax) are also accepted. Do not confuse IQUV Mueller-matrix components with TEB coefficients.
- `mueller` must be a finite, constant, real 4 x 4 matrix, specified with `hwp=True`. The default is the ideal HWP, diag(1, 1, -1, -1).
- `sy.hwp_mueller(sky, harmonics)` is a separate API for angle-dependent harmonic components and returns a perturbation relative to an ideal HWP. Do not treat it as the same API as the beam's `mueller` argument.
- Asymmetric beams may require odd moments. Do not assume that the default `spin_k=(0, 2, 4)` is sufficient; select the required moments and upper limit. The `range(5)` above is an example covering 0–4, not a guarantee of a sufficient upper limit in general.
- Supply the h-maps in `model.required_couples(basis.full)` and the differences between Basis couples needed for the orientation matrix. Include the HWP moments produced by a general Mueller matrix. `sbm.testing.assert_matches_tod` provides an example of constructing this set of couples.

### 6.3 Approximations and Additional Beam Validation

Distinguish convolution without HWP from the default shifted narrow-beam approximation with HWP. In the HWP case, if the HWP harmonic is m, the input beam order is n + m. Do not add an extra complex conjugation to the internal `C = T M T†`.

`low_ell="beam_sum"` assumes a Mueller matrix uniform in the beam-centered Ludwig-3 frame and replaces the I-input ell=0,1 2f terms arising from M_QI / M_UI with a beam multipole sum. It requires `hwp=True`. It does not make all channels or broad beams exact. See [hwp_beam_convolution.rst](../docs/source/hwp_beam_convolution.rst) for details.

When modifying a beam, validate against the pencil-beam Stokes matrix product, ideal-HWP limit, identity Mueller matrix, and an independent Wigner sum for asymmetric beams. For changes involving low ell or offsets, also check the approximation's validity range and convergence.

```sh
.venv/bin/python -m pytest -q tests/test_hwp_convolution.py tests/test_beam_low_ell.py tests/test_beam_pointing.py
```

## 7. Tests and Completion Criteria

Testing has two layers. Provide both for a new effect; an effect without an approximation needs only the first.

### 7.1 Agreement at the Same Approximation Order

Add a case for the new function to `cases()` in tests/test_systematics.py using `add(fn, name, hwp, model, tod, absorbed=None)`. `model` is the complete signal, normally `sky_signal + perturbation`; a beam model supplies the complete signal itself. `tod(pixels, psi, phi)` directly expresses the same signal in real space. `test_model_matches_tod` checks that `map_make` matches a least-squares fit of that TOD in each Basis in `NO_HWP_BASES` or `HWP_BASES`, according to the `hwp` flag. The tolerances are the `assert_matches_tod` defaults, rtol = atol = 1e-9, where atol is relative to the maximum TOD amplitude. Even for an approximate model, truncating the reference TOD at the same order makes agreement possible within numerical roundoff. Keep the default tolerances for this layer and assess the approximation's validity in Section 7.2.

- Write the expected TOD directly from the physical model in real space, without reconstructing it from the implementation's spin components. `sbm.testing.RealSky` provides random values and derivatives in each pixel, and `detector()` samples an offset pointing through a Taylor expansion using SBM angles. The polarization-angle and Mueller-matrix cases use `lbs_polarizer_tod` and `lbs_mueller_tod`, written in litebird_sim angles (psi' = -psi, chi = phi), to keep the reference independent of SBM's angle convention as well. Follow one of these approaches for new cases.
- A case with `hwp=False` is also checked for the absence of HWP harmonics (all m = 0).
- Give each main public branch its own case: single detector, pair, HWP on/off, `second_order`, and so on. If a case needs additional random numbers, use `np.random.default_rng` with a fixed seed; the `bandpass_mismatch` case is an example.
- Verify that an ordinary perturbation vanishes at zero error and that the composed signal returns to the reference signal.
- If templates are declared, provide `absorbed=lambda: (base, base + perturbation)`. `test_templates_absorb_the_first_order_perturbation` map-makes both models using `Basis.iqu(hwp=True)` for a case with `hwp=True`, or `Basis.qu(hwp=False)`, extended by the templates. It requires the polarization (`iqu[1:]`) to agree with atol = 1e-12. Select parameters whose perturbation the templates absorb exactly, such as a first-order perturbation with `second_order=False`.
- For new templates, consider adding a Basis containing them to `HWP_BASES` or `NO_HWP_BASES` so that TOD agreement is checked in that Basis for all cases as well.
- Check spin weights (Section 3.2). If an entry in `UNCHECKED` in tests/test_spin.py is needed, document the reason and the alternative independent validation.

The existing tests also check that every registered function has a case and that every function with templates has an absorbed case (`test_every_systematic_has_a_case`), that it is listed in docs/source/systematics.rst (`test_every_systematic_is_documented`), and that public exports agree with `ALL`.

### 7.2 Approximation Truncation Error

For an approximation such as a small-angle expansion, compare against the exact signal in a separate test and verify that the residual decreases at the next order beyond those retained as the error magnitude is varied. [test_pointing_geometry.py](../tests/test_pointing_geometry.py) demonstrates how to construct an exact reference signal and check the residual hierarchy. It samples a band-limited sky exactly at the displaced position and, at fixed RHO = 2e-3, checks that the maximum residual falls below 0.2 times the previous value as the model is improved from sky-only to first and then second order. The expected residual orders are rho, rho^2, and rho^3, but the current test does not vary rho to estimate those orders directly. When validating the convergence order of a new effect, also compare several error magnitudes. Detect incorrect directions, signs, weights, or missing transport of the polarization frame through the residual hierarchy or convergence order.

### 7.3 Completion Checklist

- [ ] Record the effect's specification, equations, units, signs, and approximation range.
- [ ] Update the function, decorator, public imports, `__all__`, and the docstring in \_\_init\_\_.py.
- [ ] Add real-space TOD comparisons at the same order, the required limit and branch cases, and an absorbed case if templates are declared. For approximations, add a convergence test for the truncation error.
- [ ] Update the table in docs/source/systematics.rst, docs/source/reference.rst, and CHANGELOG.md.
- [ ] Run the relevant tests and the full test suite, and record the results.
- [ ] Verify that the number of Ruff 0.17.4 diagnostics has not increased from the baseline and that Sphinx passes with `-W`.

For Ruff, check counts by rule and diagnostics in the changed code as well as the total count. A decrease in existing diagnostics can offset an increase in new ones, so equal counts alone do not establish that no problems were introduced.

```sh
.venv/bin/python -m pytest tests -q
.venv/bin/python -m ruff --version   # Must be 0.17.4
.venv/bin/python -m ruff check sbm tests --statistics
.venv/bin/python -m sphinx -b html -W docs/source docs/build/html
```

<a id="baseline"></a>

### 7.4 Baseline (2026-11-11, commit 780b047)

Python 3.12.2, SBM 0.4.2, NumPy 2.5.4. These results were obtained with no changes to the code in sbm/ or tests/.

| Check | Result |
| --- | --- |
| Editable install, `import sbm` | Successful; imports the source from this working directory |
| `pip check` | No dependency inconsistencies |
| `pytest tests -q` | 152 passed |
| Sphinx HTML build with `-W` | Successful |
| Ruff 0.17.4, `ruff check sbm tests` | 47 diagnostics (UP037 21, I001 6, RUF059 5, C408 4, and others) |

The 47 Ruff diagnostics are in the existing upstream code. Because pyproject.toml does not select rules, they depend on Ruff's default rules; Ruff 0.11.9 reports zero diagnostics. For future implementations, append the current results to the specification for each effect.

## 8. Reference Material

- [CONTRIBUTING.md](../CONTRIBUTING.md)
- [adding_systematics.rst](../docs/source/adding_systematics.rst)
- [conventions.rst](../docs/source/conventions.rst)
- [testing.rst](../docs/source/testing.rst)
- [hwp_beam_convolution.rst](../docs/source/hwp_beam_convolution.rst)
- [thesis_systematics.rst](../docs/source/thesis_systematics.rst) (maps the thesis equations to the code; use it for the specification's reference equations and literature)
- [adding_a_systematic.ipynb](../notebooks/adding_a_systematic.ipynb)
