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

- **Implementation progress** (2026-11-09, same day):
  - `efficiency_passthrough` (V3) materialized and `current` — first recipe
    wired end-to-end (`src/plot_efficiency.py`, E50 gate 5.4±0.1 keVnr,
    measured 5.393). Note: outputs go `behind` after any `uv add` (normal);
    `efficiency_passthrough` is currently `behind` from the h5py addition.
  - `response_model_calibration` (hdf5) materialized and `current` —
    `src/build_response_model.py`: nestpy driver calling NEST's own
    FullCalculation + GetS1 + GetS2 with the LZ_WS2024 built-in detector,
    ws_* vectors passed explicitly; 200k events per band (er_ch3t, nr),
    mean Ly/Qy tables, and band quantiles of log10(S2c) in S1c bins
    (0.0228/0.1587/0.5/0.8413/0.9772). Verified headline numbers:
    g2 = 34.02 phd/e, ExtEff = 0.7262, SE = 46.85 phd — NEST's
    CalculateG2 reproduces LZ's published g2/extraction/SES from the
    header parameters on its own.
  - **nestpy 2.1.3 gotchas** (each cost a debug cycle; do not rediscover):
    1. `FullCalculation`'s named kwargs are mislabeled — `er_parameters`
       and `nr_er_width_parameters` land in each other's C++ slots (the
       C++ master signature order is NRy, Width, ERy; the binding's
       py::arg order disagrees). Call it **positionally**:
       `(species, E, density, field, A, Z, NRy, Width, ERy, do_times)`.
    2. No `INTERACTION_TYPE`; the enum is `nestpy.interactions` with
       `NR`, `beta`, `CH3T`, `DD`, `WIMP`, ... (no generic `ER`; use
       `beta`/`CH3T`). YieldResult has `PhotonYield`/`ElectronYield` only
       (no TotalYield); QuantaResult has `photons`/`electrons`/`excitons`
       /`ions` (no excfrac).
    3. `GetQuanta` kwarg is `nr_er_width_parameters=`; `GetYields` kwargs
       are `nr_parameters=`/`er_parameters=` (these named bindings are
       correct — only FullCalculation's are scrambled).
    4. `SetDriftVelocity(T, density, field, pressure)` — note density
       before pressure. Gives 1.33 mm/us at LZ conditions (NEST model;
       true LZ ~1.4; only used for z, irrelevant under uniform_1d).
    5. Seeding: `nestpy.RandomGen.rndm().set_seed(n)` (singleton, no
       constructor). Detector has `set_s2_thr(...)`; our driver sets it
       to 0 so GetS2 keeps the full-S2 output branch (`ion[7]` = S2c phd;
       `scint[5]` = S1c phd — both are the XYZ-corrected NphdC, identity
       under uniform_1d).

- **Validations 1 & 2 done (2026-11-09, ER/NR band vs WS2022 Fig. 1)** —
  `er_band_validation` and `nr_band_validation` materialized and `current`
  (`src/validate_bands.py` + `src/digitize_ws2022_bands.py`). Both PASS with
  documented residuals; the microphysics stage is now validated on both
  bands.
  - Reference: `data/lz_ws2022_band_reference.csv`, digitized
    **programmatically from the PDF vector graphics** (pymupdf
    `page.get_drawings()` + tick-label axis calibration on Fig. 1 of
    arXiv:2207.03764v4 p4; verified against a raster overlay), not
    eyeball-digitized. Reference curves' own pointwise smoothness noise
    0.0005–0.0018 in log10(S2c) — the limit on how tightly widths compare.
  - Results (1M events, seed 43, tritium beta spectrum — see the 2026-11-10
    entry): **medians ≤0.54% (ER) / 0.44% (NR)** everywhere (LZ's own NR bar
    is 1%); **widths ≤5% relative** everywhere — ER max 3.6% (0.0037 abs),
    NR 6.9% (0.0046 abs) at S1c 40–65 (narrowest reference half-width),
    max abs 0.0081 (wide bin, ~4% rel). Acceptance operationalized as:
    median ≤5% rel; width ≤5% rel **or** 0.0065 abs in log10(S2c) (the
    digitized reference's own curve resolution). Nothing was tuned to pass.
    (Under the original flat-ER-spectrum approximation the numbers were
    ER 1.05% median / 6.0% width at S1c > 65 — that residual was the
    spectrum, not the engine; the flat toggle is kept via `--er-spectrum
    flat`.)
  - **Chain-equivalence forensics (source diffs v2.3.7 vs NEST master)**:
    RandomGen identical (rand_gauss/rand_skewGauss/binom_draw) or
    negligible; GetYields ER model hardcoded+identical, NR vector-driven
    (v2.3.7 defaults); GetQuanta equivalent under the bridges; GetS1/GetS2
    total-area paths identical (tauTrap change only jitters pulse timing;
    FindNewMean is a 1e-16-level correction at sPEres = 0.12; radmax
    diffusion edge is zero at detector center). The residuals are therefore
    attributed to the reference's provenance (LZ's DD source spectrum, their
    MC stats, figure production), not our engine.
  - **Full v2.3.7→master bridge set** (beyond the [6]-negation and
    SkewnessER=0 already recorded): width vec **[8]↔[9] swapped** (master
    RecombOmegaER reads [8]=center, [9]=width; 2.3.7-ordered vector trips
    the width>center sanity check and zeroes S2c); **er_parameters = [-1]×10
    sentinel** (master GetYieldBetaGR is parameterized — a positive vector
    overrides the hardcoded ER model element-by-element; this cost the ER
    24%/190% failure); **pad width vector to 13** with two zeros;
    **ROI re-applied at analysis level** (S2c ≥ 600 phd for WS2022; header
    bakes s2_thr = 728.4 phe uncorrected, our driver zeroes it). NR
    spectrum: DD 2.45 MeV neutrons → ~flat Xe recoil spectrum to
    E_max = 4mM/(m+M)²·E_n = 73.5 keVnr; conditional band width is
    spectrum-dependent (support starvation inflates it), emax 73.5 chosen
    over probes at 20–100. ER spectrum: ~flat 0.5–18.6 keVee (beta
    endpoint). Declaring `src/validate_bands.py` as an ASTRA **input** is
    the sanctioned way to make script edits cascade to outputs.

- **Validations 1 & 2 in the WS2024 configuration done (2026-11-10)** —
  `er_band_validation_ws2024` and `nr_band_validation_ws2024` materialized
  and `current`; both PASS. `src/validate_bands.py` now takes
  `--config {ws2022,ws2024}` (default ws2022 — existing recipes unchanged);
  ws2022 results bit-identical after the refactor.
  - **ER spectrum choice (decision `calibration_er_spectrum`, default
    `tritium_beta`)**: the user noticed the simulated band sat slightly high
    vs LZ at high S1c in both configs — root cause: the flat spectrum
    (uniform 0.5–18.6 keVee) over-weights the tritium endpoint region, and
    at fixed S1c higher-energy events sit at higher S2c. Sampling the true
    tritium beta spectrum (dN/dE ∝ p·E·(E0−E)²·F(Z=2, E), endpoint
    18.590 keV) removes the bias and the ws2022 S1c > 65 width residual
    (ER medians ≤0.54%/0.36% and widths ≤3.6% now; ws2024 max width 0.0120
    abs, inside the 0.015 raster floor). Toggle `--er-spectrum flat` kept
    (bit-identical to the original path). NR needs no choice: monoenergetic
    2.45 MeV DD neutrons give an exactly flat single-scatter recoil
    spectrum.
  - Reference: `data/lz_ws2024_band_reference.csv`, digitized
    programmatically from **slide 9 of the collaboration's TeVPA results
    talk** (`data/lz_ws2024_tevpa_talk.pdf`, sha256-pinned) — the WS2024
    paper does not publish the bands. The slide's plot is a pasted raster,
    so extraction is by **color segmentation** (ER curves light blue, NR
    red, distinct from navy/orange scatter), axes calibrated on all 41
    inward-pointing ticks per axis (fit residuals max 0.042 phd /
    0.0022 dex; measured second-difference noise 0.0020 median / 0.005 p90
    — the raster reference's acceptance floor is 0.015 abs, vs ws2022's
    0.0065 for the vector-digitized one).
  - Results (1M events, seed 43, tritium beta spectrum): medians ≤0.36%
    (ER) / 0.64% (NR); widths ER max 12.2% rel (0.0120 abs, S1c 46), NR max
    13.8% rel (0.0120 abs, S1c 42–61) — both inside the 0.015 raster floor.
  - **WS2024 bridges differ from WS2022**: the header ships tuned
    `NRYieldsParam` AND `ERYieldsParam` vectors (positive er_parameters
    override the hardcoded ER model — intended v2.4.0 behavior);
    `[6] = 0.3` is an **absolute Fano under v2.4.0 semantics, passed
    as-is** (no negation); `[8]↔[9]` swap still applies; **ER skewness
    vector-driven** (SkewnessER = -999 sentinel reads [10] = -0.26564; no
    header disable-instruction).
  - **The slide's band quantiles were computed on an ROI-applied sample**,
    even though its scatter is not ROI-cut (extends to log10 S2c = 2.53):
    without the cut the NR S1c = 3 half-width is 12% too wide; with
    **S2c ≥ 493 phd = 14.5 extracted electrons × g2 = 34.0** it is within
    2% below S1c = 9. (645 phd uncorrected / 1.214 = 531 also passes; 645
    corrected is too aggressive at 13.6%.)
  - **SES wrinkle resolved (2026-11-10, user's investigation)**: no
    contradiction — the numbers live in two spaces. Table S2 is the NEST
    model in corrected space: g2 = ExtEff × SES by LZ's own definition
    (measured "as the ratio of the best-fit electroluminescence gain g2 to
    the observed single electron signal size"); 0.725 × 46.9 = 34.0, and
    CalculateG2 reproduces it independently (34.02 = 0.7262 × 46.85; same
    structure in WS2022: 0.805 × 58.5 = 47.09). The ROI parenthetical is in
    raw/uncorrected space: "a lower bound on uncorrected S2 pulse size of
    645 phd, equivalent to 14.5 extracted electrons" → 645/14.5 = 44.5
    phd/e = the slide's SES, the raw (uncorrected, average-position)
    single-electron size (14.5 × 44.5 = 645.25). The ratio 46.9/44.5 =
    1.054 ≈ exp(t_mean/τ) for τ ≈ 9–10 ms is the average S2 position
    correction (electron lifetime + x-y map at mean drift time), NOT a PMT
    efficiency ("phd" already counts detected photoelectrons; P_dphe = 1.214
    is separately explicit). Empirical support: NEST's SE = 46.85 reproduces
    the band medians to ≤0.6% in corrected space, and a corrected-space ROI
    at 645 is far too aggressive (13.6% too narrow at S1c = 3) — a raw 645
    maps to a position-DEPENDENT corrected cut, softer than any single
    corrected threshold, consistent with the 493–531 window the bands
    prefer. Nothing in our chain changes: we are wholly in corrected space,
    where 46.9 is the data-agreeing value.
  - **NR mid-range residual is reference-side**: our engine gives nearly
    identical NR half-widths under ws2022 and ws2024 parameters at fixed
    S1c (diffs ≤0.005 abs), and the slide band is ~1.10–1.13× wider than
    LZ's own WS2022 paper band over S1c 40–61 (which we reproduce at the
    usual 4–7%). Plausibly position-dependent broadening absent under
    uniform_1d — the anticipated escalation path to `simplified_position`.

