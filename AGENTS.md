# LZ-lightcone

A from-scratch model of a dual-phase xenon dark-matter experiment, modeled on
LUX-ZEPLIN (LZ), that computes sensitivity to spin-independent WIMP-nucleon
scattering. The spec is `astra.yaml` — the source of truth for structure,
decisions, and evidence. The analysis is a forward model: no observed data is
used; the validation target is LZ's published **median sensitivity** (expected
limit) from the WS2024 result (arXiv 2410.17036): **5.1e-48 cm2 at 40 GeV/c2**.

The user's working mode, agreed during scoping: build and validate each
component (microphysics → backgrounds → statistics) *before* running
end-to-end, with a validation at every stage. The explicit validation list
(six items, with acceptance criteria) is declared as outputs in the first
post-scoping session — see Project Notes.

## Commands

```bash
uvx astra-tools@0.2.18 validate                        # validate spec + universes (run after every astra.yaml edit)
uvx astra-tools@0.2.18 validate astra.yaml --verify-evidence   # + verify all quotes against cached PDFs
uvx astra-tools@0.2.18 paper list                      # cached papers (5 analysis papers, version-pinned)
lc status                                              # state of declared outputs
```

When implementation starts: `uv add <pkg>` (never pip), probes via
`lc run local -- python src/...`, real outputs via `lc materialize local` —
see the lightcone skill; never write into `results/` by hand.

## Project Notes

### Scoping outcome (2026-11-09)

- **Eleven decisions**, all defaults confirmed by the user. The consequential
  ones and why:
  - `response_model: nest_library` — **use the NEST library itself via the
    nestpy Python bindings** (user's call: "no need to write new code if we
    don't need to"), not our own reimplementation. The released nestpy wheel
    ships `nestpy.detectors.LZ_WS2024` as a built-in; the downloaded HEPData
    header (`lz_nest_header_ws2024`) is the cross-check that the packaged
    parameters match LZ's published ones, with a documented source-build
    fallback if they ever differ. The original NEST-review and original-NEST
    paper insights still back the physics understanding; the
    `nest_implementation` option (own reimplementation) remains in the
    decision space as a non-default alternative.
  - `nest_parameter_set: ws2024_nest240` — LZ's WS2024 used NEST 2.4.0
    (g1 = 0.112, g2 = 34.0, SES 46.9 phd, extraction 72.5%, drift 97 V/cm)
    and published a reproducing header on HEPData (record 155182); that header
    becomes a declared input. The WS2022 set (NEST 2.3.7, g1 = 0.114,
    g2 = 47.1) is kept as a comparison option.
  - `efficiency_treatment: input_published_curve` — the published NR
    efficiency curve is an **input**, not recreated; event-level cut details
    are not public (user's explicit call).
  - `statistics: both_with_crosscheck` — 2D binned profile likelihood in
    (S1c, log10 S2c) as primary + cut-and-count with Feldman-Cousins 90% CL
    as cross-check. LZ itself is **unbinned** over six samples; our binned 2D
    likelihood is a documented variation.
  - `background_components: all_table1` — every Table I source simulated
    individually (12 background sources, 214Pb dominant at 743 ± 88 events;
    total 1210 ± 91 pre-fit for the 3.3 t-y WS2024 exposure).
  - `detector_effects: uniform_1d` — no position-dependent effects; **all
    outputs, the ROI, and the likelihood live in corrected-observable space
    (S1c, log10 S2c)** — with a uniform response the LZ-style corrections are
    the identity, so the generated response is already corrected-space
    (user probed this explicitly; do not generate raw S1/S2 in this option).
    Electron-lifetime / field effects are absorbed into effective g1/g2 and
    LZ's tuned (reduced) band widths. If ER/NR bands validate too narrow,
    `simplified_position` is the escalation universe.
  - `median_sensitivity_definition: median_of_toy_limits` — median of
    background-only toy-MC limits (LZ's Fig. 5 construction); toy MC, not
    asymptotic, per the conventions paper (XENON1T saw ~50% bias).
  - `exposure_structure: combined_single` — pool 4.2 tonne-years into one
    sample; the six-sample split only buys nuisance decoupling.
  - `test_statistic_sidedness: two_sided`, `power_constraint: applied_1sigma`
    (both = LZ), `halo_model: shm_standard` with parameters **fixed, not
    profiled** (conventions paper).

- **Agreed validation list** (to be declared as outputs with acceptance
  criteria in the first post-scoping session, before implementation):
  1. ER yield/band vs LZ tritium band — median and width better than ~5% in
     log10(S2c) (LZ's own tuning bar).
  2. NR yield/band vs DD-neutron band, same style.
  3. Efficiency passthrough — published curve plotted vs LZ Fig. 2 (sanity
     only; it is an input).
  4. Background spectra — per-component projections vs LZ. Note LZ Fig. 4 is
     a **1D projection along the ER-band median**, not a stacked 2D spectrum;
     compare projections + total counts vs Table I pre-fit values.
  5. Median sensitivity curve vs LZ Fig. 5 dashed line — minimum ≈
     5.1e-48 cm2 at 40 GeV/c2, shape agreement.
  6. Statistical cross-check — FC cut-and-count curve vs 2D-PLR curve
     consistency.

- **Data quirks / paper facts that matter downstream**:
  - LZ's published median uses **post-fit** background estimates from
    background-only pseudo-experiments; a background-only forward model
    matches this naturally (fit returns true rates).
  - Table I counts are for the 3.3 t-y WS2024 exposure only; the combined
    analysis adds the 0.9 t-y WS2022 sample.
  - Detector neutrons are derived **in situ** (veto-tagged sideband), and
    accidental coincidences are **data-driven** — both need a modeling
    decision when we get there (published-description simulation).
  - Electron-capture decays have enhanced recombination:
    Qy^L/Qy^beta = 0.87 ± 0.03; needed for the 124Xe/127Xe backgrounds, which
    sit off the ER band. 124Xe LL charge yield is a nuisance bounded
    [0.65, 0.87] in LZ.
  - ROI: 3 < S1c < 80 phd, log10(S2c) < 4.5, uncorrected S2 > 645 phd
    (14.5 extracted electrons) — WS2024 values (WS2022 used S2 > 600 phd,
    S2c < 1e5 phd).

- **Literature**: 5 papers cached version-pinned (see `astra paper list`);
  141 prior insights, all 209 evidence quotes verified against the PDFs.
  Extraction was done per-paper by subagents against the versioned cache
  entries; dropped quotes were never re-added by hand.

- **nestpy ↔ HEPData header cross-check (2026-11-09, probe
  `src/probe_nest_crosscheck.py`)**: the PyPI wheel's
  `nestpy.detectors.LZ_WS2024` matches the HEPData header on 24/25 scalars
  and on the full NR/ER yield vectors. Two deltas, both resolved in favor of
  the nestpy/NEST-master values:
  1. `s2_thr`: header 0.0 vs nestpy 783.0 (phe). 783 phe = 645 phd x
     (1 + P_dphe = 1.214) — it is the ROI uncorrected-S2 bound expressed in
     phe; upstream NEST bakes the ROI threshold in, the HEPData header
     applies it externally. **Our driver sets s2_thr back to 0** and applies
     the ROI at analysis level.
  2. Width-vector element order: the HEPData header *comment* lists
     `0.15505, 0.46894` at indices 8/9, but NEST.cpp consumes [8] = width,
     [9] = center of the omega distribution, and both NEST master's
     LZ_WS2024.hh and the nestpy wheel have `0.46894, 0.15505` — consistent
     with NEST's own default vector structure (0.205, 0.45). The header
     comment is the transposed one. nestpy's packaged
     `ws_nr_er_width_parameters` (13 elements, two trailing zeros) matches
     NEST master exactly.
  Operational rule: **use nestpy's built-in detector and its `ws_*` vectors
  as-is** (self-consistent with its vendored NEST), pass them explicitly to
  GetYields/GetQuanta per the header's instructions, and let validations 1/2
  (band width at the 5% level) be the empirical arbiter.

- **Post-scoping progress** (2026-11-09): inputs/outputs skeleton declared.
  - Inputs (all in `data/`): `lz_nest_header_ws2024` (HEPData DOI
    10.17182/hepdata.155182.v2/r1), `lz_efficiency_curve` and
    `lz_si_sensitivity_reference` (HEPData ins2841863 tables t4/t1; the SI csv
    was rebuilt from the table YAML to carry all eight curves — the plain csv
    download only has the observed limit), `lz_table1_counts_reference`
    (authored from verified Table I quotes).
  - Outputs (11, no recipes yet): `response_model_calibration` (hdf5),
    `background_models` (hdf5), `wimp_signal_models` (hdf5),
    `plr_sensitivity_curve` (csv), `fc_sensitivity_curve` (csv), and the six
    validations (`er_band_validation`, `nr_band_validation`,
    `efficiency_passthrough`, `background_projection_validation`,
    `sensitivity_vs_lz_validation`, `fc_plr_crosscheck`) with acceptance
    criteria in their descriptions.
  - Still to declare when wiring validations 1/2/4: digitized LZ band and
    Fig. 4 projection references (from WS2024 figures; not in HEPData).
  - Remaining agenda: implement with staged validation; then
    background-variation studies.
