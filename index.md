# LZ-lightcone

This project builds a forward model of a dual-phase xenon experiment modeled
on LUX-ZEPLIN (LZ). It propagates WIMP signals and the background components
through corrected observables for a future spin-independent sensitivity
calculation; it does not use observed events as input.

## Background component models

The construction uses the published {astra}`inputs.lz_table1_counts_reference`
and the calibrated NEST response in {astra}`outputs.response_model_calibration`.
The energy-to-observable response is checked separately in the
{astra}`outputs.er_band_validation_ws2024` and
{astra}`outputs.nr_band_validation_ws2024` comparisons.

{astra}`outputs.background_models` stores the per-source physical expectations
and generation diagnostics. The physical source means are imposed after the
response and ROI, not inferred from the number of retained Monte Carlo events.
Its comparison with LZ's published background projection is
{astra}`outputs.background_projection_validation`.

:::{astra} outputs.background_projection_validation
:caption: Background-projection validation. The numeric sidecar records all six acceptance gates and the exclusive-source total.
:::

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
physical Table I counts. It is also the chain required for future absolute
WIMP signal yields.

:::{astra} outputs.efficiency_passthrough
:caption: The public efficiency stages and the independent black-curve/HEPData consistency check.
:::

:::{astra} decisions.efficiency_treatment
:::

## Quantifying the correction

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
The remaining corrected-observable ROI is S1c < 80 phd, S2c >= 493 phd,
and log10(S2c) < 4.5, applied exactly once. There is no additional raw-S2
cut, threefold coincidence, 1D-analysis energy cut, or total-efficiency
factor. NEST does not reconstruct or unfold the true energy used by this
selection.

This construction is an energy-dependent average acceptance. It does not
recover the unpublished dependence of reconstruction and quality cuts on
S1c, S2c, position, or their correlations; that limitation also applies to
the NR energy-space treatment. Table I still fixes each background's central
physical mean after selection. Monte Carlo multiplicities, effective sample
sizes, and observable shapes are simulation outputs, not changes to those
physical rates. Accidentals remain an observable-space model and never pass
through this ER energy chain.

{astra}`outputs.er_efficiency_comparison` archives complete old/new
regeneration of all twelve sources at matched budgets, all eight ER energy-
space source shapes, exact NR/accidentals controls, and both evaluations of
the six unchanged Validation 4 gates. Cumulative-Y and the independent
HEPData total-2D curve are forward diagnostics only: neither is multiplied
into the production selection. Residual disagreement is retained rather than
fitted away, with central-curve extraction and finite-MC precision stated
separately.

At full production budgets, the active model and both matched comparison
modes pass the same six gates. The active-model reproduction, full/display
normalization, and NR/accidental controls pass as well. Measured stage
acceptances and shape changes are recorded in
{astra}`findings.er_staged_source_acceptances_and_shapes`; the independent
forward-check and summed-background results are in
{astra}`findings.er_efficiency_forward_and_validation4`. Residual differences
from cumulative Y and independent total-2D t4 remain diagnostics of the
central-curve, energy-conditional approximation, not additional production
weights.

:::{astra} outputs.er_efficiency_validation
:caption: Complete staged ER selection comparison. The conditional G/Y stage uses retained sampled true energy. Unit-shape differences and Monte Carlo multiplicities are distinct from the unchanged Table-I physical means; the independent total-2D curve is never another selection factor.
:::
