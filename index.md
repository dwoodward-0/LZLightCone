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
