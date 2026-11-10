"""Response model calibration (output: response_model_calibration, hdf5).

Runs NEST itself (via nestpy) with LZ's WS2024 parameter set to produce
the microphysics tables everything downstream consumes:

  - mean ER and NR yields (Ly, Qy) on an energy grid,
  - toy events in corrected-observable space (S1c, S2c) for the tritium
    (CH3T) ER band and the DD-style NR band, with NEST's own detector
    model (g1, gas gain, double-PHE, noise) applied via GetS1/GetS2,
  - banded summaries: quantiles of log10(S2c) in S1c bins.

Detector: nestpy's built-in LZ_WS2024, used as-is per the operational rule
recorded in AGENTS.md (validated against the HEPData header by
src/probe_nest_crosscheck.py). One documented deviation: s2_thr is set to
0 so GetS2 keeps the full-S2 (phe/phd) output branch instead of flagging
pulses below the baked-in ROI bound (783 phe); the ROI is applied at
analysis level.

Position handling: uniform_1d (decision detector_effects) — events are
generated at fixed (x, y) = (0, 0) with drift time sampled uniformly over
the detector's dt range, so the lifetime correction cancels in
expectation and the generated observables are already corrected-space
( FitS1 = FitS2 = 1.0, PosRes = 0 in the LZ header mean exactly this).
"""

import argparse
import sys
import time

import h5py
import nestpy
import numpy as np

DENSITY = 2.90  # g/cm3, liquid xenon at LZ operating conditions
DRIFT_FIELD = 96.5  # V/cm

# quantile levels for the band summaries: median +/- 1sigma, +/- 2sigma
QUANTILE_LEVELS = [0.0228, 0.1587, 0.5, 0.8413, 0.9772]

# S1c binning for the band summaries (phd)
S1C_BIN_EDGES = np.arange(0.0, 82.0, 2.0)

MIN_EVENTS_PER_BIN = 50  # quantiles reported only for bins with enough events

# interaction per band: CH3T is LZ's tritium ER calibration source,
# NR is the generic nuclear recoil (DD-neutron-like)
BANDS = {
    "er_ch3t": {
        "interaction": nestpy.interactions.CH3T,
        "A": 0.00054858,  # electron mass in amu
        "Z": 1.0,
        "emin": 0.5,  # keVee
        "emax": 40.0,
        "yield_emin": 0.5,
        "yield_emax": 40.0,
    },
    "nr": {
        "interaction": nestpy.interactions.NR,
        "A": 131.293,
        "Z": 54.0,
        "emin": 1.0,  # keVnr
        "emax": 100.0,
        "yield_emin": 1.0,
        "yield_emax": 100.0,
    },
}


def generate_band_events(nc, det, band, n_events, g2_params, rng,
                         drift_velocity):
    """Simulate n_events uniformly in energy; return (S1c, S2c) in phd."""
    species = band["interaction"]
    s1c = np.empty(n_events)
    s2c = np.empty(n_events)
    dt_min, dt_max = det.get_dt_min(), det.get_dt_max()
    top_drift = det.get_TopDrift()

    t0 = time.time()
    for i in range(n_events):
        energy = rng.uniform(band["emin"], band["emax"])
        # positional call, C++ order (nestpy 2.1.3 mislabels these kwargs:
        # the named er_parameters/nr_er_width_parameters land in each
        # other's slots -- see probe in AGENTS.md notes)
        result = nc.FullCalculation(
            species,
            energy,
            DENSITY,
            DRIFT_FIELD,
            band["A"],
            band["Z"],
            list(det.ws_nr_parameters),
            list(det.ws_nr_er_width_parameters),
            list(det.ws_er_parameters),
            False,
        )
        dt = rng.uniform(dt_min, dt_max)
        z = top_drift - dt * drift_velocity
        pos = (0.0, 0.0, z)  # uniform_1d: detector center line
        wf_time = []  # not requested (outputTiming = 0)
        wf_amp = []
        scint = nc.GetS1(
            result.quanta,
            *pos,
            *pos,
            drift_velocity,
            drift_velocity,
            species,
            i,
            DRIFT_FIELD,
            energy,
            nestpy.S1CalculationMode.Full,
            0,
            wf_time,
            wf_amp,
        )
        ion = nc.GetS2(
            result.quanta.electrons,
            *pos,
            *pos,
            dt,
            drift_velocity,
            i,
            DRIFT_FIELD,
            nestpy.S2CalculationMode.Full,
            0,
            wf_time,
            wf_amp,
            g2_params,
        )
        # scint[5] = NphdC (S1c, phd), ion[7] = NphdC (S2c, phd)
        s1c[i] = scint[5]
        s2c[i] = ion[7]
    rate = n_events / (time.time() - t0)
    print(f"  {rate:.0f} events/s", file=sys.stderr)
    return s1c, s2c


def mean_yield_table(nc, band, n_points=200):
    """Mean Ly/Qy on an energy grid (no fluctuations)."""
    species = band["interaction"]
    energies = np.linspace(band["yield_emin"], band["yield_emax"], n_points)
    ly = np.empty(n_points)
    qy = np.empty(n_points)
    for i, energy in enumerate(energies):
        y = nc.GetYields(
            species,
            float(energy),
            density=DENSITY,
            drift_field=DRIFT_FIELD,
            A=band["A"],
            Z=band["Z"],
        )
        ly[i] = y.PhotonYield / energy
        qy[i] = y.ElectronYield / energy
    return energies, ly, qy


def band_summary(s1c, s2c):
    """Quantiles of log10(S2c) in S1c bins (NaN where too few events)."""
    log10_s2c = np.log10(s2c)
    n_bins = len(S1C_BIN_EDGES) - 1
    centers = np.empty(n_bins)
    counts = np.empty(n_bins, dtype=int)
    quantiles = np.empty((len(QUANTILE_LEVELS), n_bins))
    for b in range(n_bins):
        lo, hi = S1C_BIN_EDGES[b], S1C_BIN_EDGES[b + 1]
        sel = (s1c >= lo) & (s1c < hi)
        counts[b] = sel.sum()
        centers[b] = 0.5 * (lo + hi)
        if counts[b] >= MIN_EVENTS_PER_BIN:
            quantiles[:, b] = np.quantile(
                log10_s2c[sel], QUANTILE_LEVELS
            )
        else:
            quantiles[:, b] = np.nan
    return centers, counts, quantiles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True, help="hdf5 path to write")
    ap.add_argument("--n-events", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--er-emin", type=float, default=BANDS["er_ch3t"]["emin"])
    ap.add_argument("--er-emax", type=float, default=BANDS["er_ch3t"]["emax"])
    ap.add_argument("--nr-emin", type=float, default=BANDS["nr"]["emin"])
    ap.add_argument("--nr-emax", type=float, default=BANDS["nr"]["emax"])
    args = ap.parse_args()
    BANDS["er_ch3t"]["emin"], BANDS["er_ch3t"]["emax"] = (
        args.er_emin,
        args.er_emax,
    )
    BANDS["nr"]["emin"], BANDS["nr"]["emax"] = args.nr_emin, args.nr_emax

    det = nestpy.detectors.LZ_WS2024()
    det.set_s2_thr(0.0)  # documented deviation: ROI applied at analysis level
    nc = nestpy.NESTcalc(det)
    nestpy.RandomGen.rndm().set_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    g2_params = nc.CalculateG2(0)
    g2 = g2_params[3]
    print(f"g2 (ExtEff * SE) = {g2:.3f} phd/electron")
    print(f"ExtEff = {g2_params[1]:.4f}, SE = {g2_params[2]:.3f} phd")
    drift_velocity = nc.SetDriftVelocity(
        det.get_T_Kelvin(), DENSITY, DRIFT_FIELD, det.get_p_bar()
    )
    print(f"drift velocity = {drift_velocity:.4f} mm/us")

    with h5py.File(args.output, "w") as f:
        f.attrs["g1"] = det.get_g1()
        f.attrs["g2"] = g2
        f.attrs["g1_gas"] = det.get_g1_gas()
        f.attrs["ExtEff"] = g2_params[1]
        f.attrs["SE"] = g2_params[2]
        f.attrs["P_dphe"] = det.get_P_dphe()
        f.attrs["sPEres"] = det.get_sPEres()
        f.attrs["noiseLinear_S1"] = det.get_noiseLinear()
        f.attrs["s2Fano"] = det.get_s2Fano()
        f.attrs["E_gas_kVcm"] = det.get_E_gas()
        f.attrs["drift_field_Vcm"] = DRIFT_FIELD
        f.attrs["density_gcm3"] = DENSITY
        f.attrs["eLife_us"] = det.get_eLife_us()
        f.attrs["dt_min_us"] = det.get_dt_min()
        f.attrs["dt_max_us"] = det.get_dt_max()
        f.attrs["nest_version"] = nestpy.__nest_version__
        f.attrs["seed"] = args.seed
        f.attrs["n_events_per_band"] = args.n_events
        f.attrs["quantile_levels"] = np.array(QUANTILE_LEVELS)
        f.attrs["s1c_bin_edges"] = S1C_BIN_EDGES

        for name, band in BANDS.items():
            print(f"band {name}: {args.n_events} events, "
                  f"E in [{band['emin']}, {band['emax']}] keV")
            s1c, s2c = generate_band_events(
                nc, det, band, args.n_events, g2_params, rng, drift_velocity
            )
            centers, counts, quant = band_summary(s1c, s2c)
            energies, ly, qy = mean_yield_table(nc, band)

            grp = f.create_group(name)
            grp.attrs["interaction"] = band["interaction"].name
            grp.attrs["emin_keV"] = band["emin"]
            grp.attrs["emax_keV"] = band["emax"]
            grp.create_dataset("s1c_phd", data=s1c)
            grp.create_dataset("s2c_phd", data=s2c)
            grp.create_dataset("band_s1c_centers", data=centers)
            grp.create_dataset("band_counts", data=counts)
            grp.create_dataset(
                "band_log10s2c_quantiles", data=quant
            )
            grp.create_dataset("yield_energy_keV", data=energies)
            grp.create_dataset("yield_ly_per_keV", data=ly)
            grp.create_dataset("yield_qy_per_keV", data=qy)

            # sanity gate: band median must rise with S1c, tolerating
            # small quantile wobble well below the 5% validation bar
            med = quant[2]
            finite = np.isfinite(med)
            steps = np.diff(med[finite])
            if np.any(steps < -0.05):
                print(
                    f"FAIL: {name} band median not monotonic in S1c "
                    f"(largest step {steps.min():.3f})",
                    file=sys.stderr,
                )
                return 1
            if np.nanmean(quant[2][:5]) <= 0:
                print(f"FAIL: {name} band median non-positive at low S1c",
                      file=sys.stderr)
                return 1
            med5 = med[~np.isnan(med)]
            print(f"  median log10(S2c) at S1c={centers[np.isfinite(med)][0]:.0f}"
                  f" -> {med5[0]:.3f}; last finite bin "
                  f"S1c={centers[np.isfinite(med)][-1]:.0f} -> {med5[-1]:.3f}")

    print(f"PASS: wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
