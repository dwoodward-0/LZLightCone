# LZ-lightcone

A from-scratch forward model of a dual-phase xenon dark-matter experiment,
modeled on LUX-ZEPLIN (LZ), that computes the median sensitivity to
spin-independent WIMP–nucleon scattering from public information only. No
observed data is used. NEST microphysics (through nestpy) maps a WIMP signal
and each of LZ's twelve WS2024 background sources into corrected observables
(S1c, log10 S2c). A binned 2D profile likelihood, cross-checked by
Feldman–Cousins cut-and-count, turns these into expected limits for the WS2024
3.3 tonne-year exposure. The validation target is LZ's WS2024-only median
expected limit.

**📄 Read the full report: <https://dwoodward-0.github.io/LZLightCone/>**

It covers the headline results, all seven validations (which all pass), the
modeling decisions, the efficiency methodology and the per-source background
comparisons.

## Headline results

Median expected 90% CL limit vs LZ's WS2024-only result:

![Sensitivity vs LZ WS2024-only](https://dwoodward-0.github.io/LZLightCone/figures/sensitivity_vs_lz_validation.png)

Projection to 1000 live days, profile likelihood vs cut-and-count:

![1000 live-day projection](https://dwoodward-0.github.io/LZLightCone/figures/plr_fc_exposure_comparison.png)

## How the project is organized

- **`astra.yaml`** is the spec and source of truth. It declares every input,
  output, modeling decision (with its alternatives) and finding, with evidence
  quotes from the literature. Validate it with
  `uvx astra-tools@0.2.18 validate`.
- **`src/`** holds one script per output. Each output's recipe in
  `astra.yaml` says exactly how it was made.
- **`results/baseline/`** holds the outputs. Each sits beside a
  `.<id>.manifest.json` run record (recipe, decisions, git commit, input
  hashes).
- **`data/`** holds the declared inputs: HEPData tables, the LZ NEST detector
  header, and curves digitized from the public papers by `src/digitize_*.py`.
- **`index.md`** is the MyST source of the report. Every number in it is a
  live reference to a declared output.
- **`AGENTS.md`** holds the working notes kept during the analysis: decision
  rationale, debugging findings and provenance details.

Outputs are produced and tracked with the [Lightcone](https://github.com/LightconeResearch)
`lc` CLI (`lc status`, `lc materialize`).

## Getting the data

`results/` and `data/` are stored with [git-annex](https://git-annex.branchable.com/),
with the content in this repository's GitHub LFS storage. Browsing them on
GitHub shows small pointer files rather than the figures and archives. To get
the content:

```bash
git clone git@github.com:dwoodward-0/LZLightCone.git
cd LZLightCone
git annex init
git annex enableremote origin-lfs
git annex get .
```

Not distributed:

- **The five LZ papers and slides** used as inputs (public, but not ours to
  redistribute). [`PAPERS.md`](PAPERS.md) gives download links and SHA256
  hashes.
- **`data/colleague_feldman_cousins/FeldmanCousins.py`**, a colleague's
  Feldman–Cousins implementation used by the cut-and-count cross-check.

`git annex get .` reports these as unavailable, which is expected.

## Rebuilding the report site

```bash
tools/publish_pages.sh   # builds with MyST, strips HDF5 archives, pushes gh-pages
```
