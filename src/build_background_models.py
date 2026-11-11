"""Background models (output: background_models, hdf5).

Builds the (S1c, log10 S2c) distribution of each of the twelve WS2024
Table I background sources (decision background_components.all_table1),
pooled into one 4.2 tonne-year sample (decision
exposure_structure.combined_single), by sampling each source's recoil
energy from its published spectral shape and pushing it through the same
nestpy response chain as the validated band work.

Spectral shapes (all public, LZ's own):
  - ER sources: the digitized per-source curves of the LZ 2018 projection
    paper (Fig. 4 right panel, input lz2018_component_spectra) — 214Pb
    (222Rn curve), 212Pb (220Rn curve), 85Kr, 136Xe 2nubb, solar-nu ER
    (RRPA), Det.+Sur.+Env. ER.
  - NR-neutrino sources: the digitized Fig. 5 curves 8B, hep, Atm.
  - Detector neutrons: lz_neutron_spectrum_2018 (already digitized).
  - Tritium + 14C: analytic beta spectra (endpoints 18.590 / 156.476 keV,
    Fermi functions for the Z=2 / Z=7 daughters), activity ratio 8:1
    (tritium:14C, WS2024 text).
  - 124Xe DEC: monoenergetic lines 10.00 keV (BR 1.4%) and 5.98 keV
    (BR 0.8%) with the published EC charge suppression
    Qy_L/Qy_beta = 0.87 +/- 0.03 applied at the quanta boundary (see
    ResponseEngine.generate: S2 from the suppressed electron count, the
    recombined charge added to S1 as a Gaussian).
  - 127Xe/125Xe EC: 5.2 keV L-shell line with the same suppression
    (REVIEW ITEM: coarse stand-in for the full EC cascade).
  - Accidentals: reconstructed in observable space (no NEST) as a
    factorized P(S1c)*P(log10 S2c) PDF (LZ's own construction,
    arXiv:2211.17120 Sec. V F) fitted to the WS2024 Fig. 3 68%/95%
    quantile contours (decision accidentals_model.factorized_contour_constrained);
    the alternative sr1_marginals_rescaled option rescales the digitized
    SR1 marginals instead.

Response: nestpy LZ_WS2024 built-in detector + the ws2024 header's tuned
NR/ER yield vectors and width vector (same bridges as the ws2024 band
validations: [8]/[9] swap, pad to 13, s2_thr = 0, drift field from
FitEF = 97 V/cm), uniform_1d (decision detector_effects.uniform_1d —
generated observables are already corrected-space). ER events use the
CH3T species (identical to the validated ER chain).

Efficiency: the published NR data-selection efficiency curve
(lz_efficiency_curve) is applied by rejection to the three NR components
(neutrons, atm nu, 8B+hep nu), whose Table I counts include it. ER
components get no efficiency folding (LZ's ER event-level cuts are not
public); their normalization absorbs the efficiency.

ROI (WS2024, corrected space, applied at analysis level): 3 < S1c < 80
phd, log10(S2c) < 4.5, S2c >= 493 phd (= 14.5 extracted electrons x
g2 = 34.0; the convention established by the ws2024 band validations).

Normalization: each component is normalized to its Table I pre-fit mean
count (lz_table1_counts_reference); the stored per-event weights are the
Table I counts, so component totals match Table I by construction and
the fit stage later floats the rates. Weighted samples + a fine 2D
histogram are stored per source, plus the accidentals fit diagnostics
and an SR1-marginals morphology cross-check.
"""

import argparse
import csv
import re
import sys
import time
from pathlib import Path

import h5py
import nestpy
import numpy as np

DENSITY = 2.90  # g/cm3 (same as the production response model)

# WS2024 ROI in corrected space (see module docstring)
ROI_S1C_MIN, ROI_S1C_MAX = 3.0, 80.0
ROI_LOG10S2C_MAX = 4.5
ROI_S2C_MIN = 493.0

# observable-space histogram grid stored per source
HIST_S1C_EDGES = np.arange(3.0, 80.5, 1.0)
HIST_Y_EDGES = np.arange(2.69, 4.51, 0.02)

# output sample cap per source (weighted samples are thinned to this)
MAX_SAMPLES = 50_000

# EC charge suppression (pre-fit): Qy_L/Qy_beta = 0.87 +/- 0.03
Q_EC = 0.87

# tritium / 14C analytic spectra
TRITIUM_ENDPOINT = 18.590  # keV
C14_ENDPOINT = 156.476  # keV
TRITIUM_C14_RATIO = 8.0 / 1.0  # activity ratio (WS2024 text)
ALPHA = 1.0 / 137.036

# 124Xe DEC lines (keVee, branch fraction within the decay)
XE124_LINES = [(10.00, 1.4), (5.98, 0.8)]
# 127Xe/125Xe EC: single L-shell stand-in (REVIEW ITEM, see docstring)
XE127_EC_LINE = 5.2  # keVee

# Table I row -> spectral shape. 'kind': curve (2018 ER close-up),
# nr_curve (2018 Fig. 5 + efficiency), analytic, accidentals.
SOURCE_SPECS = {
    "214Pb_betas": ("curve", "er_pb214"),
    "kr85_ar39_det_gammas": ("curve", "er_kr85"),
    "solar_nu_er": ("curve", "er_solar_nu_rrpa"),
    "pb212_po218_betas": ("curve", "er_pb212"),
    "tritium_c14_betas": ("analytic", "tritium_c14"),
    "xe136_2nubb": ("curve", "er_xe136"),
    "xe124_dec": ("analytic", "xe124_lines"),
    "xe127_xe125_ec": ("analytic", "xe127_ec"),
    "accidental_coincidences": ("accidentals", None),
    "atm_nu_nr": ("nr_curve", "nr_atm"),
    "b8_hep_nu_nr": ("nr_curve", "nr_b8+nr_hep"),
    "detector_neutrons": ("nr_curve", "neutron_2018"),
}

# documented modeling notes stored as hdf5 attrs per source
SOURCE_NOTES = {
    "214Pb_betas": "222Rn curve of LZ-2018 Fig. 4 (214Pb beta spectrum as "
                   "modeled by LZ, incl. all branches); WS2024's 743-count "
                   "row exceeds the naked-beta 11.0% BR count, consistent "
                   "with the summed spectrum",
    "kr85_ar39_det_gammas": "85Kr curve stands in for the whole row: the "
                            "WS2024 text states these ER sources produce "
                            "flat spectra essentially indistinguishable in "
                            "the ROI; the det-gamma high-E tail (>155 keV) "
                            "lies outside the S1c ROI regardless",
    "solar_nu_er": "Solar nu (free electron) (RRPA) curve; RRPA suppression "
                   "below ~25 keV included (the published solid curve)",
    "pb212_po218_betas": "220Rn curve (= 212Pb beta spectrum; the "
                         "backgrounds paper considered only 212Pb in the "
                         "220Rn chain; 218Po beta is subdominant)",
    "tritium_c14_betas": "analytic beta spectra, endpoints 18.590/156.476 "
                         "keV, Fermi functions F(Z=2)/F(Z=7), activity "
                         "ratio 8:1 tritium:14C",
    "xe136_2nubb": "136Xe curve of LZ-2018 Fig. 4 (2nu-beta-beta spectrum)",
    "xe124_dec": "monoenergetic lines 10.00 keV (BR 1.4%) + 5.98 keV "
                 "(BR 0.8%), EC charge suppression q=0.87 applied at the "
                 "quanta boundary (S2 from 0.87 Ne; recombined charge "
                 "added to S1 as Gaussian(dNph*g1, sqrt(dNph*g1)))",
    "xe127_xe125_ec": "5.2 keV L-shell EC line stand-in for the full "
                      "127Xe/125Xe cascade, q=0.87 (REVIEW ITEM: coarse)",
    "accidental_coincidences": "factorized P(S1c)*P(log10 S2c) PDF "
                               "(LZ's construction) fitted to the WS2024 "
                               "Fig. 3 68%/95% quantile contours; "
                               "observable-space, no NEST response",
    "atm_nu_nr": "Atm curve of LZ-2018 Fig. 5; published NR efficiency "
                 "applied by rejection",
    "b8_hep_nu_nr": "8B + hep curves of LZ-2018 Fig. 5, weighted by their "
                    "efficiency-folded ROI integrals from the curves' own "
                    "absolute normalization; published NR efficiency "
                    "applied by rejection",
    "detector_neutrons": "Det.+Sur.+Env. neutron spectrum "
                         "(lz_neutron_spectrum_2018); published NR "
                         "efficiency applied by rejection",
}


def parse_header_vectors(header_path):
    """Same extraction as the band validations (ws2024 header)."""
    text = Path(header_path).read_text()

    def braced_vector(name):
        m = re.search(rf"{name}\s*=\s*\{{([^}}]*)\}}", text)
        if not m:
            return None
        return [float(x) for x in re.findall(r"-?\d+\.?\d*", m.group(1))]

    width = braced_vector("NRERWidthsParam")
    if width is None:
        raise ValueError("width vector not found in header")
    if len(width) != 11:
        raise ValueError(f"width vector: expected 11 elements, got {len(width)}")
    return width, braced_vector("NRYieldsParam"), braced_vector("ERYieldsParam")


def read_csv_data(path):
    """Data rows of a provenance-headered CSV (skips # comments + header)."""
    rows = []
    with open(path) as f:
        for row in csv.reader(f):
            if not row or not row[0] or row[0].startswith("#"):
                continue
            if any(ch.isalpha() for ch in row[0]):
                continue  # header line
            rows.append(row)
    return rows


def load_component_spectra(path):
    """{component: (kind, array)} from lz2018_component_spectra.csv.

    kind 'bin': array (n, 3) = [e_lo, e_hi, rate]; kind 'sample':
    (n, 2) = [e, rate]."""
    out = {}
    rows = list(csv.reader(Path(path).read_text().splitlines()))
    for row in rows:
        if not row or not row[0] or row[0].startswith("#"):
            continue
        if row[0] == "figure":
            continue  # header (first fields are alphabetic: fig4, ...)
        fig, comp, kind, elo, ehi, lgr = row[:6]
        lo, hi, r = float(elo), float(ehi), 10.0 ** float(lgr)
        entry = out.setdefault(comp, [kind, []])
        assert entry[0] == kind
        entry[1].append((lo, hi, r))
    # kind 'sample': collapse to (n, 2) = [e, rate] (e_hi == e_lo there;
    # keeping the 3-column form would invite reading column 1 as rate)
    return {k: (v[0], np.array(v[1])[:, [0, 2]] if v[0] == "sample"
                else np.array(v[1])) for k, v in out.items()}


def spectrum_pdf_sampler(kind, arr, rng):
    """Return sampler(n) -> energies ~ the spectrum (inverse CDF).

    bin rows: P(bin) ~ rate * width. sample rows: piecewise-constant
    density = rate between consecutive vertex energies (trapezoid
    weighting)."""
    if kind == "bin":
        edges = np.append(arr[:, 0], arr[-1, 1])
        rate = arr[:, 2]
        w = rate * np.diff(edges)
        cdf = np.cumsum(w)
        cdf /= cdf[-1]

        def sampler(n):
            u = rng.uniform(0.0, cdf[-1] * (1 - 1e-12), n)
            idx = np.searchsorted(cdf, u)
            return rng.uniform(edges[idx], edges[idx + 1])
        return sampler, (edges[0], edges[-1])

    e = arr[:, 0]
    r = arr[:, 1]
    w = 0.5 * (r[1:] + r[:-1]) * np.diff(e)
    cdf = np.concatenate([[0.0], np.cumsum(w)])
    cdf /= cdf[-1]

    mean_curve = float(np.sum(w * 0.5 * (e[1:] + e[:-1])) / np.sum(w))
    _checked = False

    def sampler(n):
        nonlocal _checked
        u = rng.uniform(0.0, cdf[-1] * (1 - 1e-12), n)
        idx = np.clip(np.searchsorted(cdf, u) - 1, 0, len(w) - 1)
        # sample uniformly within the segment (segment-wise flat is a
        # good approximation at the figure's vertex density)
        out = e[idx] + rng.uniform(0.0, 1.0, n) * (e[idx + 1] - e[idx])
        # self-check on the first call: the sampled mean must track the
        # curve's mean (catches rate-column mixups, which otherwise
        # bias silently)
        if not _checked:
            assert abs(out.mean() - mean_curve) < 0.1 * mean_curve + 0.5, \
                f"sampler mean {out.mean():.2f} != curve mean " \
                f"{mean_curve:.2f} keV"
            _checked = True
        return out
    return sampler, (e[0], e[-1])


def sample_tritium_c14(rng, n):
    """Mixed tritium + 14C beta spectrum (8:1 by activity)."""
    def pdf(E, E0, Z):
        p = np.sqrt(E * (E + 2.0 * 511.0))
        eta = Z * ALPHA * E / p
        fermi = 2.0 * np.pi * eta / (1.0 - np.exp(-2.0 * np.pi * eta))
        return np.where((E > 0) & (E <= E0), p * E * (E0 - E) ** 2 * fermi, 0.0)

    grids = []
    for E0, Z in [(TRITIUM_ENDPOINT, 2.0), (C14_ENDPOINT, 7.0)]:
        grid = np.linspace(0.05, E0, 3000)
        f = pdf(grid, E0, Z)
        cdf = np.concatenate([[0.0], np.cumsum(f)])
        cdf /= cdf[-1]
        grids.append((grid, cdf))
    # mixture: draw the isotope by activity weight, then invert its CDF
    out = np.empty(n)
    iso = rng.choice(2, size=n, p=[8.0 / 9.0, 1.0 / 9.0])
    for k, (grid, cdf) in enumerate(grids):
        sel = iso == k
        if not sel.any():
            continue
        u = rng.uniform(0.0, 1.0, sel.sum())
        # cdf[j] = P(E <= grid[j-1]) for j >= 1; invert against cdf[1:]
        out[sel] = np.interp(u, cdf[1:], grid)
    return out


class ResponseEngine:
    """The validated ws2024 GetYields -> GetQuanta -> GetS1/GetS2 chain."""

    def __init__(self, header_path, seed):
        self.det = nestpy.detectors.LZ_WS2024()
        self.det.set_s2_thr(0.0)  # ROI applied at analysis level
        self.nc = nestpy.NESTcalc(self.det)
        nestpy.RandomGen.rndm().set_seed(seed)
        self.rng = np.random.default_rng(seed)
        self.drift_field = self.det.FitEF(0, 0, 0)  # 97 V/cm (ws2024)

        width, nr_y, er_y = parse_header_vectors(header_path)
        width_vec = list(width)
        width_vec[8], width_vec[9] = width_vec[9], width_vec[8]  # [8]<->[9]
        self.width_vec = width_vec + [0.0, 0.0]  # pad to 13
        if nr_y is None or er_y is None:
            raise ValueError("ws2024 header must carry tuned NR/ER vectors")
        self.nr_params, self.er_params = nr_y, er_y

        self.g2_params = self.nc.CalculateG2(0)
        self.g2 = self.g2_params[3]
        if abs(self.g2 - 34.0) > 0.2:
            raise ValueError(f"g2 {self.g2:.3f} != published 34.0")
        self.drift_velocity = self.nc.SetDriftVelocity(
            self.det.get_T_Kelvin(), DENSITY, self.drift_field,
            self.det.get_p_bar(),
        )

    def generate(self, energies, nr, n_events=None, q_ec=None):
        """Push energies through the response; returns (S1c, S2c).

        nr: True for the NR chain (NR species, A=131.3, Z=54), False for
        ER (CH3T species, as in the validated ER band). q_ec: optional EC
        charge suppression (Qy_L/Qy_beta), applied at the quanta boundary:
        S2 is generated from the suppressed electron count (GetS2 takes
        the count directly), and the recombined charge's photon
        contribution is added to S1c as a Gaussian with the correct mean
        (dNph * g1) and counting variance (sqrt(dNph * g1)).
        nestpy's QuantaResult is immutable (def_readonly), so the extra
        photons cannot be pushed through GetS1 itself; the approximation
        is subdominant to the published +/-0.03 uncertainty on q."""
        species = nestpy.interactions.NR if nr else nestpy.interactions.CH3T
        A, Z = (131.293, 54.0) if nr else (0.00054858, 1.0)
        n = len(energies)
        det, nc = self.det, self.nc
        dt_min, dt_max = det.get_dt_min(), det.get_dt_max()
        top_drift = det.get_TopDrift()
        g1 = det.get_g1()
        s1c = np.empty(n)
        s2c = np.empty(n)
        t0 = time.time()
        for i in range(n):
            energy = float(energies[i])
            yields = nc.GetYields(
                species, energy, density=DENSITY,
                drift_field=self.drift_field, A=A, Z=Z,
                nr_parameters=self.nr_params,
                er_parameters=self.er_params,
            )
            quanta = nc.GetQuanta(
                yields, DENSITY,
                nr_er_width_parameters=self.width_vec,
                SkewnessER=-999.0,  # vector-driven (ws2024 ER: [10])
            )
            dphot = 0.0
            if q_ec is not None:
                ne = quanta.electrons
                dphot = (1.0 - q_ec) * ne  # charge re-emitted as light
            dt = self.rng.uniform(dt_min, dt_max)
            z = top_drift - dt * self.drift_velocity
            pos = (0.0, 0.0, z)  # uniform_1d
            wf_time, wf_amp = [], []
            scint = nc.GetS1(
                quanta, *pos, *pos, self.drift_velocity, self.drift_velocity,
                species, i, self.drift_field, energy,
                nestpy.S1CalculationMode.Full, 0, wf_time, wf_amp,
            )
            ne_for_s2 = int(round(q_ec * quanta.electrons)) \
                if q_ec is not None else quanta.electrons
            ion = nc.GetS2(
                ne_for_s2, *pos, *pos, dt, self.drift_velocity, i,
                self.drift_field, nestpy.S2CalculationMode.Full, 0,
                wf_time, wf_amp, self.g2_params,
            )
            s1c[i] = scint[5]
            if dphot > 0.0:
                mean_add = dphot * g1
                s1c[i] += self.rng.normal(mean_add, np.sqrt(mean_add))
            s2c[i] = ion[7]
        rate = n / (time.time() - t0)
        print(f"    {rate:.0f} events/s", file=sys.stderr)
        return s1c, s2c


def roi_mask(s1c, s2c):
    y = np.log10(np.maximum(s2c, 1e-6))
    return ((s1c > ROI_S1C_MIN) & (s1c < ROI_S1C_MAX)
            & (y < ROI_LOG10S2C_MAX) & (s2c >= ROI_S2C_MIN))


def load_efficiency(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    arr = np.array(rows)
    assert np.all(np.diff(arr[:, 0]) > 0)

    def eff(e):
        return np.interp(e, arr[:, 0], arr[:, 1])
    return eff


def read_table1(path):
    """{source: (mean, sigma)} from the Table I reference CSV. Not
    read_csv_data: source names are alphabetic ('total' included)."""
    out = {}
    rows = list(csv.reader(Path(path).read_text().splitlines()))
    for row in rows:
        if not row or not row[0] or row[0].startswith("#"):
            continue
        if row[0] == "source":
            continue  # header
        out[row[0]] = (float(row[1]), float(row[2]))
    assert "total" in out
    return out


# ---------------- accidentals reconstruction ----------------

def fit_factorized_accidentals(contour_path, n_s1_knots=13, n_y_knots=12,
                               lam=0.5):
    """Fit log f1(S1c) + log f2(y) = c_q to the 68%/95% quantile contours.

    The quantile contour of a unimodal PDF is a level set of the PDF, so
    on each contour log f1 + log f2 is constant. Solve on piecewise-
    linear knot parameterizations with second-difference smoothness and
    the gauge u(S1c=3) = 0. The level-set equations alone are degenerate
    (the constant solution u+v = c68 = c95 satisfies every pointwise
    equation), so the fit anchors the level GAP, c68 - c95 = 1, which
    excludes constants and picks one member of the remaining scale
    family; the physical scale is recovered afterwards by requiring the
    PDF p ~ exp(alpha*(u+v)) to enclose 68% of its ROI mass inside the
    fitted q68 contour (the enclosed 95% volume is the consistency
    check).
    Contour points lying ON the ROI edges are dropped: there the contour
    is clipped by the plotted window (it continues outside it), so it
    does not constrain the PDF's level sets — keeping them would force
    the marginals flat along the edges. S1c knots stop at 60 phd (the
    q95 contour's reach); beyond that f1 is extrapolated with the local
    log-slope (see below). Returns (s1c_grid, f1, y_grid, f2, rms,
    (s1_knots, u_knots, y_knots, v_knots, c68, c95, vol68, vol95))."""
    pts = {"q68": [], "q95": []}
    n_dropped = 0
    for row in csv.reader(Path(contour_path).read_text().splitlines()):
        if not row or not row[0] or row[0].startswith("#"):
            continue
        if row[0] not in ("q68", "q95"):
            continue  # header
        s1, y = float(row[1]), float(row[2])
        on_edge = (abs(s1 - ROI_S1C_MIN) < 0.1 or abs(y - np.log10(
            ROI_S2C_MIN)) < 0.005 or abs(y - ROI_LOG10S2C_MAX) < 0.005)
        if on_edge:
            n_dropped += 1
        else:
            pts[row[0]].append((s1, y))
    assert set(pts) == {"q68", "q95"}
    print(f"  accidentals: dropped {n_dropped} edge-clipped contour "
          f"points, kept "
          f"{len(pts['q68']) + len(pts['q95'])}")

    # knots
    s1_knots = np.linspace(ROI_S1C_MIN, 60.0, n_s1_knots)
    y_knots = np.linspace(np.log10(ROI_S2C_MIN), ROI_LOG10S2C_MAX, n_y_knots)
    K, M = n_s1_knots, n_y_knots

    def interp_row(x, knots):
        """Piecewise-linear basis row: value at x over the knot grid."""
        row = np.zeros(len(knots))
        idx = np.clip(np.searchsorted(knots, x) - 1, 0, len(knots) - 2)
        t = (x - knots[idx]) / (knots[idx + 1] - knots[idx])
        row[idx] = 1.0 - t
        row[idx + 1] = t
        return row

    A_rows, b_rows = [], []
    for q in ("q68", "q95"):
        for s1, y in pts[q]:
            row = np.zeros(K + M + 2)
            row[:K] = interp_row(s1, s1_knots)
            row[K:K + M] = interp_row(y, y_knots)
            row[K + M + ("q68", "q95").index(q)] = -1.0
            A_rows.append(row)
            b_rows.append(0.0)
    # gauge: u(S1c=3) = 0 (fixes the additive translation; the overall
    # normalization of the PDF cancels in the sampling)
    row = np.zeros(K + M + 2)
    row[0] = 1.0
    A_rows.append(row)
    b_rows.append(0.0)
    # level-gap anchor: c68 - c95 = 1. The level-set equations alone are
    # degenerate — the constant solution u+v = c68 = c95 satisfies every
    # pointwise equation, and so does every near-flat member that the
    # min-norm criterion prefers. Fixing the level GAP excludes
    # constants and picks one member of the remaining scale family; the
    # physical scale is recovered below from the 68% quantile volume.
    row = np.zeros(K + M + 2)
    row[K + M] = 1.0
    row[K + M + 1] = -1.0
    A_rows.append(row)
    b_rows.append(1.0)
    # smoothness: second differences on both knot vectors
    for vec_off, n_knots in ((0, K), (K, M)):
        for i in range(1, n_knots - 1):
            row = np.zeros(K + M + 2)
            row[vec_off + i - 1] = lam
            row[vec_off + i] = -2.0 * lam
            row[vec_off + i + 1] = lam
            A_rows.append(row)
            b_rows.append(0.0)
    Amat = np.array(A_rows)
    bvec = np.concatenate([
        np.zeros(len(pts["q68"]) + len(pts["q95"])),  # level equations
        [0.0],  # gauge
        [1.0],  # gap anchor
        np.zeros((K - 2) + (M - 2)),  # smoothness
    ])
    theta, *_ = np.linalg.lstsq(Amat, bvec, rcond=None)
    u_knots = theta[:K]
    v_knots = theta[K:K + M]
    c68_fit, c95_fit = theta[K + M], theta[K + M + 1]

    # recover the physical scale: find alpha such that the PDF
    # p ~ exp(alpha*(u+v)) puts 68% of its ROI mass inside the fitted
    # q68 level set {u+v > c68_fit}. The contour SETS are invariant
    # under the scale family; only the integrand weighting depends on
    # alpha. The 95% volume is the consistency check.
    s1_g = np.arange(ROI_S1C_MIN + 0.05, ROI_S1C_MAX, 0.25)
    y_g = np.arange(np.log10(ROI_S2C_MIN) + 0.01, ROI_LOG10S2C_MAX, 0.02)
    w = (np.interp(s1_g, s1_knots, u_knots)[:, None]
         + np.interp(y_g, y_knots, v_knots)[None, :])

    def volume_above(level, alpha):
        p = np.exp(alpha * (w - w.max()))  # avoid overflow
        return p[w > level].sum() / p.sum()

    lo, hi = 1e-4, 1e4
    if volume_above(c68_fit, hi) < 0.68:
        print("  WARNING: no alpha puts 68% inside the q68 contour "
              "(volume saturates below target); keeping alpha = "
              f"{lo:.3g}", file=sys.stderr)
        alpha = lo
    else:
        while hi / lo > 1.001:
            mid = np.sqrt(lo * hi)
            if volume_above(c68_fit, mid) > 0.68:
                hi = mid
            else:
                lo = mid
        alpha = np.sqrt(lo * hi)
    u_knots = alpha * u_knots
    v_knots = alpha * v_knots
    c68 = alpha * c68_fit
    c95 = alpha * c95_fit
    vol68 = volume_above(c68_fit, alpha)
    vol95 = volume_above(c95_fit, alpha)

    # data residual in log space (fit quality of the level-set
    # equations), evaluated on the final scaled solution — the residuals
    # scale with alpha, so this is the physically meaningful number
    res = []
    for q in ("q68", "q95"):
        cq = c68 if q == "q68" else c95
        for s1, y in pts[q]:
            model = (np.interp(s1, s1_knots, u_knots)
                     + np.interp(y, y_knots, v_knots) - cq)
            res.append(model)
    rms = float(np.sqrt(np.mean(np.square(res))))
    print(f"  accidentals: level-fit RMS {rms:.4f} dex; scale alpha = "
          f"{alpha:.3f}; enclosed volumes {vol68:.4f} (target 0.68), "
          f"{vol95:.4f} (target 0.95)")

    # evaluate on the ROI grid. Beyond the last S1c knot (60 phd, the q95
    # contour's reach) the shape is publicly unconstrained: extrapolate
    # f1 with the local log-slope, floored at zero slope (decaying tail,
    # never a rise) — a flat hold would overestimate the tail.
    s1_grid = np.arange(ROI_S1C_MIN + 0.05, ROI_S1C_MAX, 0.1)
    ds1 = s1_knots[-1] - s1_knots[-2]
    tail_slope = min((u_knots[-1] - u_knots[-2]) / ds1, 0.0)
    u_right = u_knots[-1] + tail_slope * (s1_grid[-1] - s1_knots[-1])
    f1 = np.exp(np.interp(s1_grid, s1_knots, u_knots, right=u_right))
    y_grid = np.arange(np.log10(ROI_S2C_MIN) + 0.002, ROI_LOG10S2C_MAX, 0.005)
    f2 = np.exp(np.interp(y_grid, y_knots, v_knots))
    return s1_grid, f1, y_grid, f2, rms, (s1_knots, u_knots, y_knots,
                                          v_knots, c68, c95, vol68, vol95)


def load_sr1_marginals(path):
    """The two SR1 marginal functions as (x, rate) arrays."""
    s1x, s1r, yx, yr = [], [], [], []
    mode = None
    for row in csv.reader(Path(path).read_text().splitlines()):
        if not row or not row[0] or row[0].startswith("#"):
            continue
        if row[0] == "s1c_phd":
            mode = "s1"
            continue
        if row[0] == "log10_s2c":
            mode = "y"
            continue
        try:
            x, v = float(row[0]), float(row[1])
        except ValueError:
            continue
        if mode == "s1":
            s1x.append(x), s1r.append(v)
        else:
            yx.append(x), yr.append(v)
    return (np.array(s1x), np.array(s1r)), (np.array(yx), np.array(yr))


def sample_from_marginals(rng, x, pdf_x, y, pdf_y, n):
    """Independent draws from two normalized-on-grid marginals."""
    p1 = pdf_x / pdf_x.sum()
    p2 = pdf_y / pdf_y.sum()
    c1 = np.cumsum(p1)
    c2 = np.cumsum(p2)
    s1 = np.interp(rng.uniform(0, 1, n), c1, x)
    yy = np.interp(rng.uniform(0, 1, n), c2, y)
    return s1, yy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--spectra", required=True,
                    help="lz2018_component_spectra.csv")
    ap.add_argument("--neutron-spectrum", required=True)
    ap.add_argument("--efficiency", required=True)
    ap.add_argument("--table1", required=True)
    ap.add_argument("--contours", required=True)
    ap.add_argument("--sr1-marginals", required=True)
    ap.add_argument("--nest-header", required=True)
    ap.add_argument("--seed", type=int, default=44)
    ap.add_argument("--n-er", type=int, default=40_000,
                    help="generated events per ER source (pre-ROI)")
    ap.add_argument("--n-nr", type=int, default=120_000,
                    help="generated events per NR source (pre-ROI)")
    ap.add_argument("--n-b8", type=int, default=600_000,
                    help="generated events for the combined 8B+hep source "
                         "(its efficiency-folded ROI pass rate is ~1e-3, "
                         "so it needs far more generated events than the "
                         "other NR sources to fill its ROI histogram)")
    ap.add_argument("--background-components", default="all_table1")
    ap.add_argument("--exposure-structure", default="combined_single")
    ap.add_argument("--nest-parameter-set", default="ws2024_nest240")
    ap.add_argument("--detector-effects", default="uniform_1d")
    ap.add_argument("--accidentals-model",
                    default="factorized_contour_constrained")
    args = ap.parse_args()
    for name, val in (("background_components", args.background_components),
                      ("exposure_structure", args.exposure_structure),
                      ("nest_parameter_set", args.nest_parameter_set),
                      ("detector_effects", args.detector_effects)):
        allowed = {"background_components": {"all_table1"},
                   "exposure_structure": {"combined_single"},
                   "nest_parameter_set": {"ws2024_nest240"},
                   "detector_effects": {"uniform_1d"}}[name]
        if val not in allowed:
            print(f"FAIL: unsupported {name} option {val!r} "
                  f"(supported: {sorted(allowed)})", file=sys.stderr)
            return 1
    if args.accidentals_model not in {"factorized_contour_constrained",
                                      "sr1_marginals_rescaled"}:
        print(f"FAIL: unsupported accidentals_model "
              f"{args.accidentals_model!r}", file=sys.stderr)
        return 1

    table1 = read_table1(args.table1)
    spectra = load_component_spectra(args.spectra)
    eff = load_efficiency(args.efficiency)
    # neutron spectrum CSV: (e_lo, e_hi, log10_rate) bins
    neutron_bins = np.array(
        [[float(r[0]), float(r[1]), 10.0 ** float(r[2])]
         for r in read_csv_data(args.neutron_spectrum)]
    )
    spectra["neutron_2018"] = ("bin", neutron_bins)
    # combined 8B+hep entry handled by splitting its mixture weight below

    engine = ResponseEngine(args.nest_header, args.seed)
    print(f"response: g2 = {engine.g2:.3f} phd/e, field = "
          f"{engine.drift_field:.1f} V/cm")

    with h5py.File(args.output, "w") as f:
        f.attrs["seed"] = args.seed
        f.attrs["nest_version"] = nestpy.__nest_version__
        f.attrs["g1"] = engine.det.get_g1()
        f.attrs["g2"] = engine.g2
        f.attrs["drift_field_Vcm"] = engine.drift_field
        f.attrs["roi_s1c"] = (ROI_S1C_MIN, ROI_S1C_MAX)
        f.attrs["roi_log10s2c_max"] = ROI_LOG10S2C_MAX
        f.attrs["roi_s2c_min_phd"] = ROI_S2C_MIN
        f.attrs["q_ec"] = Q_EC
        f.attrs["decisions"] = (
            f"background_components={args.background_components}, "
            f"exposure_structure={args.exposure_structure}, "
            f"nest_parameter_set={args.nest_parameter_set}, "
            f"detector_effects={args.detector_effects}, "
            f"accidentals_model={args.accidentals_model}"
        )
        f.attrs["n_generated"] = f"er={args.n_er}, nr={args.n_nr}"

        fit_diagnostics = f.create_group("accidentals_fit")

        for source, (kind, shape_key) in SOURCE_SPECS.items():
            if source not in table1:
                print(f"FAIL: Table I reference lacks {source}",
                      file=sys.stderr)
                return 1
            mean_count, sigma = table1[source]
            print(f"\n{source}: Table I {mean_count} +/- {sigma}")

            q_ec = None
            if kind == "accidentals":
                # observable-space PDF; no response chain
                if args.accidentals_model == "factorized_contour_constrained":
                    (s1_grid, f1, y_grid, f2, rms,
                     knots) = fit_factorized_accidentals(args.contours)
                    s1, y = sample_from_marginals(
                        engine.rng, s1_grid, f1, y_grid, f2, args.n_er)
                    # SR1 morphology cross-check (shape only), over the
                    # S1c range where the SR1 marginal has support (it
                    # starts at its own ROI edge ~3.4 phd); beyond an
                    # interpolated curve's support the held-endpoint
                    # values would inflate the difference spuriously
                    (sr1x, sr1r), (sr1yx, sr1yr) = load_sr1_marginals(
                        args.sr1_marginals)
                    fit_s1_norm = f1 / np.trapezoid(f1, s1_grid)
                    overlap = (s1_grid >= sr1x.min()) & (s1_grid <= sr1x.max())
                    sr1_norm = np.interp(s1_grid[overlap], sr1x, sr1r)
                    sr1_norm = sr1_norm / np.trapezoid(sr1_norm,
                                                       s1_grid[overlap])
                    log_ratio = np.log10(fit_s1_norm[overlap] / sr1_norm)
                    s1_xcheck_maxlog = float(np.abs(log_ratio).max())
                    s1_xcheck_medlog = float(np.median(np.abs(log_ratio)))
                    fit_diagnostics.attrs["contour_rms_log10"] = rms
                    fit_diagnostics.attrs["s1_marginal_vs_sr1_max_logdiff"] = \
                        s1_xcheck_maxlog
                    fit_diagnostics.attrs["s1_marginal_vs_sr1_median_logdiff"] \
                        = s1_xcheck_medlog
                    fit_diagnostics.attrs["s1_xcheck_overlap_phd"] = (
                        float(s1_grid[overlap].min()),
                        float(s1_grid[overlap].max()))
                    fit_diagnostics.attrs["s1_knots"] = knots[0]
                    fit_diagnostics.attrs["s1_knot_values"] = knots[1]
                    fit_diagnostics.attrs["y_knots"] = knots[2]
                    fit_diagnostics.attrs["y_knot_values"] = knots[3]
                    fit_diagnostics.attrs["c68"] = knots[4]
                    fit_diagnostics.attrs["c95"] = knots[5]
                    print(f"  accidentals fit: contour RMS {rms:.3f} dex; "
                          f"vs SR1 max|dlog| {s1_xcheck_maxlog:.2f} "
                          f"(median {s1_xcheck_medlog:.2f})")
                else:  # sr1_marginals_rescaled
                    (sr1x, sr1r), (sr1yx, sr1yr) = load_sr1_marginals(
                        args.sr1_marginals)
                    s1_grid = np.arange(ROI_S1C_MIN + 0.05, ROI_S1C_MAX, 0.1)
                    f1 = np.interp(s1_grid, sr1x, sr1r)
                    y_grid = np.arange(np.log10(ROI_S2C_MIN) + 0.002,
                                       ROI_LOG10S2C_MAX, 0.005)
                    f2 = np.interp(y_grid, sr1yx, sr1yr)
                    s1, y = sample_from_marginals(
                        engine.rng, s1_grid, f1, y_grid, f2, args.n_er)
                    fit_diagnostics.attrs["contour_rms_log10"] = float("nan")
                in_roi = np.ones(len(s1), dtype=bool)
                s2c = 10.0 ** y
                base_w = np.ones(len(s1))
            else:
                # energy-space sources
                if shape_key == "nr_b8+nr_hep":
                    # weight 8B vs hep by their efficiency-folded integrals
                    # from the curves' own absolute normalization
                    # (both are 'bin' step histograms)
                    ints = {}
                    for comp in ("nr_b8", "nr_hep"):
                        _, arr2 = spectra[comp]
                        mids = 0.5 * (arr2[:, 0] + arr2[:, 1])
                        ints[comp] = float(np.sum(
                            arr2[:, 2] * eff(mids)
                            * (arr2[:, 1] - arr2[:, 0])))
                    w_b8 = ints["nr_b8"] / (ints["nr_b8"] + ints["nr_hep"])
                    print(f"  8B:hep efficiency-folded weight = "
                          f"{w_b8:.3f}:{1 - w_b8:.3f}")
                    n_b8 = int(round(args.n_b8 * w_b8))
                    parts = []
                    for comp, n_part in (("nr_b8", n_b8),
                                         ("nr_hep", args.n_b8 - n_b8)):
                        k2, arr2 = spectra[comp]
                        sampler2, _ = spectrum_pdf_sampler(k2, arr2,
                                                           engine.rng)
                        parts.append(sampler2(n_part))
                    energies = np.concatenate(parts)
                elif kind in ("curve", "nr_curve"):
                    k2, arr2 = spectra[shape_key]
                    sampler, _ = spectrum_pdf_sampler(k2, arr2, engine.rng)
                    energies = sampler(args.n_er)
                else:
                    energies = None  # analytic: sampled below

                nr = kind == "nr_curve"
                if nr:
                    keep = engine.rng.uniform(0.0, 1.0, len(energies)) \
                        < eff(energies)
                    energies = energies[keep]
                    print(f"  efficiency: {keep.mean():.4%} of events pass")
                if energies is None:  # analytic spectra
                    if shape_key in ("xe124_lines", "xe127_ec"):
                        q_ec = Q_EC
                        energies = np.full(args.n_er, float(
                            XE124_LINES[0][0] if shape_key == "xe124_lines"
                            else XE127_EC_LINE))
                        if shape_key == "xe124_lines":
                            # mixture of the two lines by branch ratio
                            w_lo = XE124_LINES[0][1] / (XE124_LINES[0][1]
                                                        + XE124_LINES[1][1])
                            is_hi = engine.rng.uniform(size=args.n_er) > w_lo
                            energies[is_hi] = XE124_LINES[1][0]
                    elif shape_key == "tritium_c14":
                        energies = sample_tritium_c14(engine.rng, args.n_er)
                    else:
                        raise ValueError(f"unhandled analytic {shape_key}")

                s1c_all, s2c_all = engine.generate(energies, nr=nr,
                                                   q_ec=q_ec)
                m = roi_mask(s1c_all, s2c_all)
                s1, s2c = s1c_all[m], s2c_all[m]
                y = np.log10(s2c)
                base_w = np.ones(len(s1))
                in_roi = m
                print(f"  ROI acceptance {m.mean():.4f} "
                      f"({m.sum()} of {len(m)} events)")

            n_kept = len(s1)
            if n_kept == 0:
                print(f"FAIL: {source} has no ROI events", file=sys.stderr)
                return 1
            w = base_w * (mean_count / base_w.sum())

            # thin weighted samples for storage
            if n_kept > MAX_SAMPLES:
                sel = engine.rng.choice(n_kept, MAX_SAMPLES, replace=False)
            else:
                sel = np.arange(n_kept)

            grp = f.create_group(source)
            grp.attrs["table1_mean"] = mean_count
            grp.attrs["table1_sigma"] = sigma
            grp.attrs["n_generated"] = int(
                args.n_er if kind in ("curve", "analytic", "accidentals")
                else (args.n_b8 if shape_key == "nr_b8+nr_hep"
                      else args.n_nr))
            grp.attrs["n_in_roi"] = n_kept
            grp.attrs["spectrum"] = str(shape_key)
            grp.attrs["notes"] = SOURCE_NOTES[source]
            grp.create_dataset("s1c", data=s1[sel])
            grp.create_dataset("log10_s2c", data=y[sel])
            grp.create_dataset("weight", data=w[sel])
            H, _, _ = np.histogram2d(
                s1, y, bins=[HIST_S1C_EDGES, HIST_Y_EDGES],
                weights=w)
            grp.create_dataset("hist_counts", data=H)
            grp.create_dataset("hist_s1c_edges", data=HIST_S1C_EDGES)
            grp.create_dataset("hist_y_edges", data=HIST_Y_EDGES)
            print(f"  stored {len(sel)} samples; total weight "
                  f"{w.sum():.2f}")

    print(f"\nPASS: wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
