# Model-handoff and normalization audit — 2026-11-13

**Scope.** Read-only inspection of the existing response, background, and
40 GeV/c² WIMP-control products; reconciliation against the declared Table I
input and signal-builder conventions. No new WIMP masses were generated. No
full-grid production, likelihood/FC fit, region optimization, large ensemble,
or declared-output materialization was performed. All executable checks ran
under `lc run local`; this report and its JSON are review artifacts, **not**
declared `results/` outputs.

## Verdict

- **PASS — common conventions and physical normalization.** All 20 explicit
  archive checks pass. Background and control-WIMP data use the same ROI,
  binning, gains, and 3.3-tonne-year WS2024 exposure. Full/display accounting
  closes. All nine focused `test_wimp_signal.py` tests pass in 8.257 s.
- **FAIL — background MC-variance metadata semantics.**
  `hist_mc_weight_variance` contains the raw per-bin weighted second moment
  `sum(w²)`, not a variance. This does **not** affect Table I normalization,
  but the array must not be consumed as a covariance or total-yield variance.
- **LIMITATION — exposure label and scope.** The actual products are both
  3.3 tonne-year WS2024. The active `combined_single` text still describes
  4.2 tonne-years. Do not multiply either product by 4.2/3.3. A combined
  exposure requires a separately validated WS2022 addition.
- **BLOCKER for a mass curve — the declared signal archive is absent.**
  `results/baseline/wimp_signal_models.hdf5` and its manifest have never been
  materialized. The existing 40 GeV/c² control is a temporary diagnostic pilot,
  not a substitute for the declared 26-mass archive.
- **Full-grid readiness is not fully established.** Interfaces, normalization,
  and the 40 GeV/c² pilot pass; the other 25 masses' response/tail paths,
  adaptive minimum-ROI budgets, production runtime, and full production-budget
  annual comparison have not been tested by this audit. No mass-grid run was
  launched. A production run needs its own scope and approval.

## Artifacts and provenance

| Artifact | Role | SHA-256 |
|---|---|---|
| `results/baseline/background_models.hdf5` | Active 12-source background | `87bbbcc6a2dfdce2af9a21fd35bca916b0c669f33f9364f00f895ceeceddaf26` |
| `results/baseline/response_model_calibration.hdf5` | Shared response / ER band projection | `7a1a9f3d85035a93e45e42b8ea4a242730f8895a24179a2a9decfbae9e89cce8` |
| `data/lz_ws2024_table1_counts.csv` | Declared Table I input | `7f9eaab6ccdc934f8f6d1101125db59a0e499922939ed96f745b0b3959371510` |
| `data/lz_ws2024_si_sensitivity.csv` | Declared mass-grid / reference-limit input | `10a18b02ed136e110fc69cfe862384dff375f1ce19945a437c5aceb639a5bc04` |
| `/tmp/lz-wimp-40-extended-control.hdf5` | Temporary 40 GeV/c² native control | `2dafc21df4743cb7eb25c6b67bcea29eef38e697e742e4b3e79953156449ab43` |

The background manifest was generated at git commit
`454688e5c0bd04236b8ef58db135c128ca94c035`. Its recipe uses
`input_pre_roi_curve`, `staged_conditional`, `ws2024_nest240`, `uniform_1d`,
and `factorized_contour_constrained`. `lc status` reports background/response
as **behind**, specifically `"made under an earlier environment"`, not stale
or physically failed. No output was materialized here.

The signal-builder source inspected is the current
`src/build_wimp_signal_models.py`; the temporary control predates some of its
optional metadata, including group-level `exposure_structure`. The control
still records `campaign=ws2024`, exposure 3.3, `input_pre_roi_curve`,
`pre_roi_ss_data_cuts`, `ws2024_nest240`, `uniform_1d`, and `shm_standard`.
Its native preparation/response seeds are 29010114/29010341, with source seed
29000107—the seed chosen by the current external-grid rule for 40 GeV/c².

## Common corrected-observable conventions

Both models reproduce the same explicit selection:

```text
3 < S1c < 80 phd
S2c >= 493 phd
log10(S2c) < 4.5
```

The shared background and WIMP arrays have 77 × 91 bins. Their edges match
`HIST_S1C_EDGES` and `HIST_Y_EDGES` exactly: 1 phd in S1c, with the lower Y
edge equal to `log10(493) = 2.69284691927723`, followed by 0.02-dex bins
through 4.5. No rounded 2.69 lower edge is used.

The response/background gains agree numerically: `g1 = 0.1122` and
`g2 = 34.02239534018503` phd/electron; the response records `SE =
46.85019810842267` phd. Background NR selection is the named pre-ROI
`pre_roi_ss_data_cuts` stage; WIMP uses the same stage and no Table-I or
published-limit normalization. No raw-S2 > 645 cut, extra coincidence factor,
total-ROI efficiency, or second corrected-space ROI was applied here.

## Background normalization and thinning

`hist_counts` uses the full unthinned generation; `weight` belongs to a capped
display sample and is independently normalized to Table I after thinning.

| Source | Table I mean | Full histogram sum | Display N | Full ROI N | Thinned? |
|---|---:|---:|---:|---:|---|
| 214Pb_betas | 743.0 | 743.0 | 2,034 | 2,034 | No |
| kr85_ar39_det_gammas | 162.0 | 162.0 | 1,980 | 1,980 | No |
| solar_nu_er | 102.0 | 102.0 | 2,970 | 2,970 | No |
| pb212_po218_betas | 62.7 | 62.7 | 2,413 | 2,413 | No |
| tritium_c14_betas | 58.3 | 58.3 | 30,342 | 30,342 | No |
| xe136_2nubb | 55.6 | 55.6 | 145 | 145 | No |
| xe124_dec | 19.4 | 19.4 | 37,105 | 37,105 | No |
| xe127_xe125_ec | 3.2 | 3.2 | 37,705 | 37,705 | No |
| accidental_coincidences | 2.8 | 2.8 | 40,000 | 40,000 | No |
| atm_nu_nr | 0.12 | 0.12 | 50,000 | 73,980 | Yes |
| b8_hep_nu_nr | 0.06 | 0.06 | 65 | 65 | No |
| detector_neutrons | 0.0 | 0.0 | 50,000 | 55,738 | Yes |

For every uncapped source, direct display-sample histograms reproduce both
the integer full-ROI histogram and its physical-weight normalization. For the
capped sources, the full integer histogram sums to the retained full-ROI
multiplicity and the display weights separately close to their Table I means.
All display events pass the common ROI, with positive physical weights except
the intentionally zero-weight detector-neutron central shape. Detector
neutrons remain useful as a unit-area morphology but have no nonzero physical
central rate to which a normalization ratio can be taken.

The twelve unrounded component means and histograms sum to **1209.18**. The
published rounded Table I `total` row is **1210**, a difference of **−0.82
events (−0.0678%)**; the component total should not be mistaken for a fit. The
naive quadrature of symmetric stand-in per-row sigmas is 91.6899149307 versus
the table's total 91; this is not a claim that all rows are independent or
that their uncertainties encode a covariance.

### MC-variance finding

`build_background_models.py:1215-1217` writes `Hw2`, the per-bin weighted sum
of squares, as `hist_mc_weight_variance`. There is no multinomial subtraction,
and no per-source total-yield MC-variance attribute in the existing archive.
The WIMP builder uses a different, explicit
`sum(w²) - sum(w)²/N` construction and its 40 GeV/c² one-block array reproduces
that formula numerically in this audit.

Because each source's constant-weight background N is fixed, summing bins—even
after the missing term—is not a legitimate total-yield variance; per-bin
counts are multinomially correlated. `table1_sigma` is a supplied physical
background uncertainty, not a generated-MC uncertainty. Before using MC errors
for region optimization, derive the required conditional-shape uncertainty from
`hist_mc_counts` and generation multiplicities, or explicitly correct and
record the array's semantics. Do not use the current name as a ready-made
covariance. This is not a demonstrated background-shape or rate error.

## Absolute WIMP normalization at 40 GeV/c²

The reconstructed scaling is

```text
Y_pre = R_native [events/kg/day/keV integrated over sampled native support]
        × 1000 kg/tonne × 365.25 day/year × 3.3 tonne-year
        × sigma_n / 1e-36 cm².
```

At `sigma_n = 1e-36 cm²` the existing representative-day control gives:

| Quantity | Value |
|---|---:|
| Native integrated rate | 1,091,198.2563265532 events/kg/day |
| Absolute preselection yield | 1,315,248,538,306.8025 events |
| Absolute full-ROI yield | 656,002,129,289.4882 events |
| Native generated | 120,000 |
| Efficiency-selected | 70,399 |
| Full-ROI retained | 59,852 |
| Display retained | 50,000 |
| Full-histogram MC relative sigma | 0.00408530424 |
| Effective conditional ROI after efficiency | 59,852 / 70,399 = 0.8501825310 |
| Efficiency × ROI acceptance | 59,852 / 120,000 = 0.4987666667 |

The native rate was re-derived from its matching preparation integrals and
per-block sample sizes, not from a fitted scalar or display weights. The
reconstructed weighted ROI yield, histogram sum, display-weight sum, and stored
direct yield agree at floating-point tolerance. The ordinary
`Y_pre × N_ROI/N_native` closure differs by only 2.0e-15 relative. Display
thinning is 59,852 → 50,000 with re-normalized weights; it is not the source
of the absolute-yield normalization. No added >100 keV tail is present; its
stored efficiency-folded fraction at 100–500 keV is 2.497e-9 and the 500–5000
keV unadded guard fraction is zero. The control's stored prepared-minus-native-
quadrature relative integral is +0.0002667916, with an independent-preparation
integral SD of 49.720054 events/kg/day. These stored numerical controls were
left unchanged, not tuned or reinterpreted as a new physics discrepancy.

The independently generated 12-phase annual control gives 655,872,837,040.3193
full-ROI events at the same reference cross section. Its relative difference
from the representative-day result is **+0.019713%**, versus 0.577755%
one-sigma independent-MC comparison precision; the archived comparison passes.
This is a 40 GeV/c² time-convention control, not phase-resolved sensitivity.

## Full-grid and downstream handoff

All nine existing WIMP tests pass, covering exact exposure/cross-section
linearity, native sampler passthrough, analytic tail weights, representative-day
semantics, and both small full-selection round trips. The current builder's
full recipe asks for exactly 26 masses with 120,000 native events minimum,
2,000 ROI events minimum, a maximum of 20 blocks per mass, and 30,000 native
grid-control events, with a separately controlled 40 GeV/c² annual comparison.
Those are builder defaults, not a runtime benchmark or a completed-mass-grid
validation.

A 40 GeV/c² optimized-FC pilot can be discussed using the audited control's
full histogram and native normalization, explicitly labeled a diagnostic. It
cannot stand in for `wimp_signal_models`, `wimp_fig4_validation`, a full FC
curve, or the eventually required PLR/FC cross-check. `lc materialize --check`
reports exactly the WIMP signal and its two validation outputs as planned but
never materialized; all other runnable outputs are behind only for the earlier
environment.

Before statistics implementation, separately specify (1) how the intended
WS2024-only 3.3-tonne-year scope relates to the active 4.2-tonne-year combined
label/reference target and (2) how optimized FC relates to the declared
median-toy, two-sided, power-constrained statistics conventions. Existing
`astra.yaml` was not changed by this audit. The optimized region should consume
full physical histograms, not capped display weights. Its cross-section scaling
at fixed mass must be `sigma_ratio = sigma / 1e-36 cm²`, applied once to the
absolute WIMP yield; Table I normalization applies only to backgrounds.

## Audit trail

- Read-only schema probe: `/tmp/lz-model-handoff-schema.py` and
  `/tmp/lz-model-handoff-schema.json`.
- Main numerical audit: `/tmp/lz-model-handoff-audit.py`; exact outputs copied
  to `model-handoff-normalization-audit.json` in the project root.
- Focused suite: `lc run local -- python -I -m unittest discover -s src -p test_wimp_signal.py -v` — **9 tests, OK**.
- Local allocation was released after the audit. No commits, `uv add`, output
  materializations, or changes to existing review/probe artifacts were made.
