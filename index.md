---
title: Spin-independent WIMP sensitivity of an LZ-like xenon TPC — a forward model validated against LZ WS2024
short_title: LZ-lightcone
---

# Summary

This project builds, from public information only, a forward model of a
dual-phase xenon dark-matter experiment modeled on LUX-ZEPLIN (LZ), and uses it
to compute the median sensitivity to spin-independent WIMP–nucleon scattering.
No observed events are used. NEST microphysics maps a WIMP signal and each of
LZ's twelve WS2024 background sources into corrected observables
(S1c, log10 S2c). A binned two-dimensional profile likelihood then turns these
templates into a median expected limit for the WS2024 exposure of 3.3
tonne-years.

The validation target is LZ's **WS2024-only** median expected limit,
5.76 × 10⁻⁴⁸ cm² at 40 GeV/c², from slide 24 of the collaboration's TeVPA talk.
It is like-for-like with our exposure. LZ's combined 4.2 tonne-year result
(5.1 × 10⁻⁴⁸ cm²) also includes the WS2022 sample, which we do not model.

At 40 GeV/c², our median power-constrained limit is
{astra:value col=median_pc_cm2 where="mass_gev=40" sig=3}`outputs.plr_sensitivity_curve` cm²
(bootstrap uncertainty
{astra:value col=median_bootstrap_sigma_cm2 where="mass_gev=40" sig=2}`outputs.plr_sensitivity_curve` cm²).
That is within the agreed 25% of LZ's value, and the curve tracks LZ's across
the published mass range ({astra:ref}`outputs.sensitivity_vs_lz_validation`).
An independent Feldman–Cousins cut-and-count analysis gives
{astra:value col=median_pc_cm2 where="mass_gev=40" sig=3}`outputs.fc_sensitivity_curve` cm²
at the same mass. All seven validations along the chain pass (see the table
below).

:::{astra} outputs.sensitivity_vs_lz_validation
:caption: **Headline result.** Median expected 90% CL upper limit on the spin-independent WIMP–nucleon cross section from the 2D profile likelihood (3.3 tonne-years, power-constrained), with ±1σ/±2σ bands, compared with LZ's WS2024-only median and bands. The lower panel shows the ratio to LZ.
:::

:::{astra} findings.plr_sensitivity_vs_lz_ws2024_only
:::

## Projection to 1000 live days

The same machinery, with signal and every background scaled to 1000 live days
(fractional background uncertainties held fixed), gives a 40 GeV/c² median of
{astra:value col=median_pc_cm2 where="mass_gev=40" sig=3}`outputs.plr_sensitivity_curve_1000d` cm²
from the profile likelihood, versus
{astra:value col=median_pc_cm2 where="mass_gev=40" sig=3}`outputs.fc_sensitivity_curve_1000d` cm²
from cut-and-count. At the WS2024 exposure, the two methods agree within about
10%. Both counting regions optimize to nearly background-free (about 0.7
expected events), where the extra shape information adds little. With more
exposure, the likelihood's use of the whole (S1c, log10 S2c) plane starts to
pay off above about 20 GeV/c².

:::{astra} outputs.plr_fc_exposure_comparison
:caption: **Projection.** Median sensitivity of the 2D profile likelihood and of the FC cut-and-count analysis at the WS2024 exposure (220 live days) and at 1000 live days, with the PLR/FC ratio below.
:::

:::{astra} findings.plr_gain_over_fc_at_1000d
:::

## Validation summary

The analysis was built in stages: microphysics, then backgrounds, then signal,
then statistics. Each stage was validated before the next was run.

| # | Validation | Output | Acceptance | Result |
|---|---|---|---|---|
| 1 | ER band vs LZ tritium band | {astra}`outputs.er_band_validation_ws2024` | median ≤ 5%; width ≤ 5% or within the reference's digitization floor | PASS |
| 2 | NR band vs LZ DD band | {astra}`outputs.nr_band_validation_ws2024` | as above | PASS |
| 3 | Efficiency passthrough | {astra}`outputs.efficiency_passthrough` | black vs HEPData ≤ 0.01; E50 = 5.4 ± 0.1 keVnr | PASS |
| 4 | Background spectra vs LZ Fig. 4 / Table I | {astra}`outputs.background_projection_validation` | six gates (counts, panels, shape RMS, tail) | PASS 6/6 |
| 5 | Median sensitivity vs LZ WS2024-only | {astra}`outputs.sensitivity_vs_lz_validation` | within 25% at 40 GeV/c²; shape across masses | PASS |
| 6 | FC cut-and-count vs PLR | {astra}`outputs.fc_plr_crosscheck` | consistent across all masses | PASS 26/26 |
| 7 | 40 GeV/c² signal region vs LZ Fig. 3 | {astra}`outputs.wimp_fig3_validation` | 68% region 0.68 ± 0.05; 95% region 0.95 ± 0.03 | PASS |

Each figure prints its own verdict. Validations 1 and 2 were also run in the
WS2022 configuration ({astra}`outputs.er_band_validation`,
{astra}`outputs.nr_band_validation`), and both pass there too.

# Model

## Configuration and decisions

The analysis is declared in `astra.yaml`; every choice below is a named
decision, with alternatives and rationale. The consequential ones:

- **Response:** NEST itself, through the nestpy bindings, with LZ's published
  WS2024 detector header ({astra}`decisions.response_model`,
  {astra}`decisions.nest_parameter_set`). Active:
  {astra:value}`decisions.response_model`, {astra:value}`decisions.nest_parameter_set`.
- **Detector effects:** {astra:value}`decisions.detector_effects`. A uniform
  response means LZ-style position corrections are the identity, so everything
  lives directly in corrected space (S1c, log10 S2c)
  ({astra}`decisions.detector_effects`).
- **Backgrounds:** every WS2024 Table I source simulated individually
  ({astra}`decisions.background_components`); accidentals built from LZ's
  factorized construction, constrained by the WS2024 Fig. 3 contours
  ({astra}`decisions.accidentals_model`).
- **Efficiencies:** public curves as inputs, applied in stages with the ROI
  applied exactly once ({astra}`decisions.efficiency_treatment`,
  {astra}`decisions.er_efficiency_treatment`; see [Efficiency methodology](#efficiency-methodology)).
- **Statistics:** binned 2D profile likelihood as the primary method, with
  Feldman–Cousins cut-and-count as a cross-check ({astra}`decisions.statistics`).
  Median of background-only toy limits
  ({astra}`decisions.median_sensitivity_definition`), two-sided test statistic
  ({astra}`decisions.test_statistic_sidedness`), power constraint at −1σ
  ({astra}`decisions.power_constraint`), and standard halo parameters, fixed
  rather than profiled ({astra}`decisions.halo_model`).
- **Exposure:** a single combined sample ({astra}`decisions.exposure_structure`),
  realized as the WS2024 3.3 tonne-years.

The full decision registry, with every option, is in
[the appendix](#decision-registry).

## Stage 1 — microphysics

{astra}`outputs.response_model_calibration` drives NEST's own
FullCalculation/GetS1/GetS2 chain with the WS2024 detector header. On its own,
NEST reproduces LZ's published g2, extraction efficiency and single-electron
size from the header parameters. The ER calibration band uses the true tritium
β spectrum ({astra}`decisions.calibration_er_spectrum`).

:::{astra} outputs.er_band_validation_ws2024
:caption: **Validation 1.** Simulated ER (tritium) band median and ±1σ quantiles in log10 S2c vs S1c, compared with LZ's WS2024 band (digitized from the TeVPA talk).
:::

:::{astra} outputs.nr_band_validation_ws2024
:caption: **Validation 2.** Simulated NR (DD-neutron) band compared with LZ's WS2024 band.
:::

## Stage 2 — efficiencies

The published NR efficiency is an input, not something we recreate. Validation 3
checks that it is read in faithfully. See
[Efficiency methodology](#efficiency-methodology) for how the stages are
applied.

:::{astra} outputs.efficiency_passthrough
:caption: **Validation 3.** The public efficiency stages: the pre-ROI Fig. 2 stage applied before NEST, and the HEPData cumulative total (which already includes the ROI), with the independent black-curve/HEPData consistency check.
:::

## Stage 3 — backgrounds

{astra}`outputs.background_models` holds each Table I source pushed through the
response and ROI, then normalized to its Table I count. Physical means are
imposed after the response, not inferred from the number of Monte Carlo events.
Validation 4 projects the summed model onto LZ's Fig. 4 axis: the distance from
the ER band median in units of the band width, in two S1c panels.

:::{astra} outputs.background_projection_validation
:caption: **Validation 4.** Sum of all twelve background sources (red, ±1σ finite-MC band) against LZ WS2024 Fig. 4's best-fit total and its uncertainty band, in the 3–40 and 40–80 phd S1c panels. The footer lists the six acceptance gates. 124Xe is drawn at the Fig.-4 shape-optimal L-shell charge yield q = 0.78; the stored model and gates use q = 0.65.
:::

Each source is also compared individually with LZ's drawn outline, where one
exists (see [Per-source background figures](#per-source-background-figures)).

## Stage 4 — WIMP signal

{astra}`outputs.wimp_signal_models` propagates NEST's native WIMP recoil
spectrum for the 26 published masses, through the pre-ROI efficiency, the same
response and a single ROI. It is absolutely normalized for 3.3 tonne-years,
with no fit to LZ. Two checkpoints at 40 GeV/c² compare it with LZ's figures.

:::{astra} outputs.wimp_fig4_validation
:caption: **WIMP checkpoint (Fig. 4).** Absolute 40 GeV/c² signal yield, displayed at LZ's observed limit, against LZ's drawn WIMP outline in both S1c panels, with residuals.
:::

:::{astra} outputs.wimp_fig3_validation
:caption: **Validation 7.** The 40 GeV/c² signal's 2D shape against LZ WS2024 Fig. 3's purple 68%/95% regions (extracted exactly from the paper's vector paths). The box reports the share of our ROI events inside each region and the gate verdict.
:::

## Stage 5 — statistics

The primary result, {astra}`outputs.plr_sensitivity_curve`, is a binned Poisson
likelihood over the (S1c, log10 S2c) plane, with all twelve background rates as
Gaussian-constrained nuisance parameters. Limits come from a toy-MC Neyman
construction, using 2000 signal toys per tested signal strength and 1000
background-only toys per mass. The likelihood core reproduces the exact
one-bin Feldman–Cousins intervals.

The cross-check, {astra}`outputs.fc_sensitivity_curve`, optimizes a counting
region per mass by signal-to-background ranking. It evaluates the median FC
limit exactly from the Poisson median count, so no toys are needed.

:::{astra} outputs.fc_plr_crosscheck
:caption: **Validation 6.** FC cut-and-count median sensitivity against the 2D PLR, with their ratio.
:::

:::{astra} findings.fc_plr_crosscheck_agreement
:::

# Limitations

- **Exposure:** WS2024's 3.3 tonne-years only. LZ's combined result adds WS2022,
  which needs its own response, efficiency and backgrounds.
- **Likelihood:** binned and single-sample, where LZ's is unbinned over six
  samples. The power-constrained −2σ band edge repeats −1σ by construction;
  the unconstrained bands are stored for comparison with LZ's.
- **Detector:** uniform response. Position-dependent broadening is the likely
  cause of the residual NR band width in the WS2024 configuration, and of the
  deep low-charge tail in Validation 4. `simplified_position` is the
  escalation option.
- **Efficiencies:** applied as energy-conditional averages. The full
  unpublished dependence on S1c, S2c and position is not recoverable.
- **Fixed inputs:** background constraints use Table I pre-fit uncertainties.
  Halo parameters are fixed.

# Efficiency methodology

These sections document how the public efficiency curves are applied, and the
before/after comparisons that justified each choice. The comparisons are
cross-checks, not acceptance gates. The model that results passes
Validations 1–7.

## Staged NR efficiency

The public Fig. 2 efficiency is a sequence of selections, not a single
pre-response factor. The cumulative total already includes the S2 and S1 ROI.
The choice of energy-space selection is explicit in
{astra}`decisions.efficiency_treatment`; the active treatment is
{astra:value}`decisions.efficiency_treatment`.

For the corrected chain, only the pre-ROI trigger plus SS/data-analysis stage
is applied in recoil-energy space. NEST and the explicit corrected-observable
ROI follow once, before Table I normalization. Consequently this correction
changes NR generation efficiencies and shapes, but not the imposed central
physical Table I counts. It is also the chain used for the absolute WIMP
signal yields.

:::{astra} decisions.efficiency_treatment
:::

{astra}`outputs.nr_efficiency_comparison` reruns complete, matched-budget
background models under the legacy and corrected selections. It retains full
Monte Carlo histograms, unit-area NR shape comparisons, and both evaluations
of the exclusive-total Validation 4 gates. The independent-seed,
fivefold-budget B8+hep study distinguishes sparse-bin fluctuations from the
effect of the efficiency change.

:::{astra} outputs.nr_efficiency_validation
:caption: Matched NR efficiency comparison. Shape differences conditional on the ROI are not changes to the imposed physical rates.
:::

## Staged ER efficiency

Electron-recoil selection follows the public cumulative stages of
{astra}`inputs.lz_er_efficiency_paper` rather than treating the published
total efficiency as another independent weight. The active choice is
{astra:value}`decisions.er_efficiency_treatment`.

For each sampled true recoil energy, the public trigger curve B(E) is applied
before NEST. NEST's own threefold coincidence and the explicit S1c > 3 phd
selection follow. Events that pass S1 then receive one conditional
single-scatter/data-analysis acceptance C(E) = G(E)/Y(E), evaluated at the
retained sampled true energy. Y(E) is the cumulative trigger, coincidence,
and lower-S1 efficiency. G(E) is supplied by the public 1D-ROI HEPData
curve t3 through 15 keV, where the user identified its overlap with the green
SS/data-analysis curve; above 15 keV, G(E) is held at G(15) = 0.959. This is
an input convention based on that visual identification, not a separately
fitted green curve.

B, Y, and G are **not** multiplied as three independent efficiencies. The
only energy-space selections are B and the one conditional G/Y rejection.
The remaining corrected-observable ROI is S1c < 80 phd, S2c ≥ 493 phd,
and log10(S2c) < 4.5, applied exactly once. There is no additional raw-S2
cut, threefold coincidence, 1D-analysis energy cut, or total-efficiency
factor. NEST does not reconstruct or unfold the true energy used by this
selection.

This construction is an energy-dependent average acceptance. It does not
recover the unpublished dependence of reconstruction and quality cuts on
S1c, S2c, position, or their correlations; that limitation also applies to
the NR energy-space treatment. Table I still fixes each background's central
physical mean after selection. Accidentals remain an observable-space model
and never pass through this ER energy chain.

{astra}`outputs.er_efficiency_comparison` archives complete old/new
regeneration of all twelve sources at matched budgets, all eight ER energy-
space source shapes, exact NR/accidentals controls, and both evaluations of
the six unchanged Validation 4 gates. Cumulative-Y and the independent
HEPData total-2D curve are forward diagnostics only: neither is multiplied
into the production selection. Both modes pass the same six gates.

:::{astra} decisions.er_efficiency_treatment
:::

:::{astra} outputs.er_efficiency_validation
:caption: Complete staged ER selection comparison. The conditional G/Y stage uses retained sampled true energy. Unit-shape differences and Monte Carlo multiplicities are distinct from the unchanged Table-I physical means; the independent total-2D curve is never another selection factor.
:::

# Appendix

## Per-source background figures

Each figure shows one Table I source in both Fig. 4 S1c panels. The red dashed
line is this work, with a ±1σ finite-MC band. The solid black line is LZ's
drawn outline, where Fig. 4 draws one, with a gray band for the Table I
normalization uncertainty. LZ's "Other ER betas" outline holds exactly
85Kr/39Ar/detector gammas plus 136Xe, so that figure plots the sum.

:::{astra} outputs.bkg_source_pb214_betas
:::

:::{astra} outputs.bkg_source_kr85_ar39_det_gammas
:::

:::{astra} outputs.bkg_source_xe136_2nubb
:::

:::{astra} outputs.bkg_source_xe124_dec
:caption: 124Xe double electron capture. Red: L-shell charge yield q = 0.78 (Fig.-4 shape-optimal); green: LZ's fitted q = 0.70 ± 0.04. Both are analytic rescalings of the stored q = 0.65 model.
:::

:::{astra} outputs.bkg_source_xe127_xe125_ec
:::

:::{astra} outputs.bkg_source_accidental_coincidences
:::

:::{astra} outputs.bkg_source_b8_hep_nu_nr
:::

:::{astra} outputs.bkg_source_pb212_po218_betas
:::

:::{astra} outputs.bkg_source_solar_nu_er
:::

:::{astra} outputs.bkg_source_tritium_c14_betas
:::

:::{astra} outputs.bkg_source_atm_nu_nr
:::

:::{astra} outputs.bkg_source_detector_neutrons
:::

## Efficiency findings

:::{astra} findings.nr_pre_roi_chain
:::

:::{astra} findings.er_staged_conditional_chain
:::

## Decision registry

:::{astra} decisions
:::

## Outputs

:::{astra} outputs
:::
