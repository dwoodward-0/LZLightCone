"""Focused tests and native-spectrum probes for absolute WIMP signal models.

Run as a unittest suite, or use ``--native-probe output.json`` for the detailed
40 GeV/c2 pre-response numerical study. Probe files belong in /tmp, not in the
materialized results tree.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import nestpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_background_models import ResponseEngine, load_efficiency, roi_mask
from build_wimp_signal_models import generate_mass, write_model
from validate_backgrounds import load_band
from wimp_native import (
    REFERENCE_SIGMA_CM2,
    REPRESENTATIVE_DAY,
    earth_speed_km_s,
    expected_yield,
    integrate_native,
    native_provenance,
    native_tail_weight,
    phase_midpoints,
    sample_native,
    sample_native_tail,
)


def selection_diagnostics(sample, efficiency, seed):
    """True-energy efficiency only; response and the one ROI come later."""
    e = sample["energies"]
    valid = e > 0.0
    eps = efficiency(np.maximum(e, 0.0))
    eps[~valid] = 0.0  # failed native attempts cannot pass selection
    rng = np.random.default_rng(seed)
    selected = rng.uniform(size=len(e)) < eps
    w_pre = sample["metadata"]["integral"] / len(e)
    w_selected = w_pre * selected
    p = float(np.mean(selected))
    sigma_p = float(np.sqrt(max(p * (1.0 - p), 0.0) / len(e)))
    selected_mean = float(np.mean(e[selected])) if np.any(selected) else None
    mean_all = float(np.mean(e))
    return {
        "n_generated": len(e),
        "n_efficiency_selected": int(np.sum(selected)),
        "efficiency_acceptance": p,
        "efficiency_acceptance_binomial_sigma": sigma_p,
        "efficiency_folded_rate": sample["metadata"]["integral"] * p,
        "efficiency_folded_rate_mc_sigma": sample["metadata"]["integral"] * sigma_p,
        "energy_mean_kevnr": mean_all,
        "energy_sem_kevnr": float(np.std(e) / np.sqrt(len(e))),
        "selected_energy_mean_kevnr": selected_mean,
        "n_above_efficiency_support": int(np.sum(e > efficiency.energy_kevnr[-1])),
        "n_selected_above_efficiency_support": int(
            np.sum(selected & (e > efficiency.energy_kevnr[-1]))),
    }, selected, w_selected


def native_probe(args):
    mass = args.mass
    seed = args.seed
    eff = load_efficiency(args.efficiency, "pre_roi_ss_data_cuts")
    e_grid = np.linspace(0.0, 100.0, 20001)
    tail_grid = np.linspace(100.0, 500.0, 20001)
    base = sample_native(mass, args.n_events, seed)
    base_stats, base_selected, base_w = selection_diagnostics(base, eff, seed + 1)
    rate_quad = integrate_native(mass, REPRESENTATIVE_DAY, e_grid, seed + 2, eff)
    tail_quad = integrate_native(mass, REPRESENTATIVE_DAY, tail_grid, seed + 3, eff)
    report = {
        "purpose": "Pre-response WIMP numerical checks; not likelihood production.",
        "provenance": native_provenance(),
        "representative_day": REPRESENTATIVE_DAY,
        "representative_earth_speed_km_s": earth_speed_km_s(REPRESENTATIVE_DAY),
        "baseline": {
            "metadata": base["metadata"], "selection": base_stats,
            "preselection_yield_3p3ty_at_reference_sigma": expected_yield(
                base["metadata"]["integral"], 3.3),
        },
        "native_integral_vs_independent_quadrature": {
            "prepared_integral": base["metadata"]["integral"],
            "independent_integral_0_100": rate_quad["integral"],
            "isotope_grid_standard_error": rate_quad["isotope_grid_standard_error"],
            "efficiency_folded_quadrature_0_100":
                rate_quad["efficiency_folded_integral"],
            "independent_grid_step_keVnr": float(e_grid[1] - e_grid[0]),
        },
        "above_100_keVnr": {
            "integral": tail_quad["integral"],
            "isotope_grid_standard_error": tail_quad["isotope_grid_standard_error"],
            "relative_to_0_100": tail_quad["integral"] / rate_quad["integral"],
            "efficiency_folded_integral": tail_quad["efficiency_folded_integral"],
            "efficiency_folded_relative_to_0_100":
                tail_quad["efficiency_folded_integral"] /
                rate_quad["efficiency_folded_integral"],
        },
        "efficiency_support_kevnr": [
            float(eff.energy_kevnr[0]), float(eff.energy_kevnr[-1])],
        "efficiency_endpoint_convention": "linear interpolation with held endpoints",
    }

    # Check that the added bulk path is an exact passthrough, not a new sampler.
    control_n = min(args.n_events, 5000)
    control_seed = seed + 10
    nestpy.RandomGen.rndm().set_seed(control_seed)
    original = nestpy.spectra.WIMP(control_n, mass=mass, eStep=5.0,
                                   day=REPRESENTATIVE_DAY)
    control = sample_native(mass, control_n, control_seed)
    report["same_seed_original_sampler_identical"] = bool(
        np.array_equal(original, control["energies"]))
    if not report["same_seed_original_sampler_identical"]:
        raise AssertionError("WIMP metadata wrapper changed the native sample stream")

    grids = []
    for i, e_step in enumerate(args.e_steps):
        record = {"e_step": e_step}
        try:
            item = sample_native(mass, args.grid_events, seed + 100 + i,
                                 e_step=e_step)
            stats, selected, _ = selection_diagnostics(item, eff, seed + 1000 + i)
            record.update({"metadata": item["metadata"], "selection": stats})
            record["efficiency_acceptance_difference_from_default"] = (
                stats["efficiency_acceptance"] - base_stats["efficiency_acceptance"])
            record["native_rate_relative_difference_from_default"] = (
                item["metadata"]["integral"] / base["metadata"]["integral"] - 1.0)
        except (RuntimeError, ValueError) as exc:
            record["error"] = str(exc)
        grids.append(record)
    report["grid_checks"] = grids
    report["independent_seed_preparations"] = [
        sample_native(mass, 0, seed + 2000 + i)["metadata"] for i in range(3)
    ]

    if args.skip_annual:
        out = Path(args.native_probe)
        out.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"baseline": base["metadata"],
                          "selection": base_stats,
                          "native_integral_vs_independent_quadrature":
                              report["native_integral_vs_independent_quadrature"],
                          "above_100_keVnr": report["above_100_keVnr"],
                          "grid_checks": grids,
                          "independent_seed_preparations":
                              report["independent_seed_preparations"]}, indent=2))
        print(f"wrote {out}")
        return 0

    phase_records = []
    phase_hist = np.zeros(50)
    energy_edges = np.linspace(0.0, 100.0, 51)
    phase_rate_sum = 0.0
    phase_eff_rate_sum = 0.0
    phase_eff_rate_var_sum = 0.0
    for i, day in enumerate(phase_midpoints(12)):
        item = sample_native(mass, args.phase_events, seed + 3000 + i,
                             day=float(day))
        stats, selected, _ = selection_diagnostics(item, eff, seed + 4000 + i)
        n = args.phase_events
        phase_rate_sum += item["metadata"]["integral"]
        phase_eff_rate_sum += stats["efficiency_folded_rate"]
        phase_eff_rate_var_sum += stats["efficiency_folded_rate_mc_sigma"] ** 2
        phase_hist += np.histogram(
            item["energies"], bins=energy_edges,
            weights=np.full(n, item["metadata"]["integral"] / n))[0]
        phase_records.append({"day": float(day),
                              "earth_speed_km_s": earth_speed_km_s(float(day)),
                              "metadata": item["metadata"], "selection": stats})
    annual_rate = phase_rate_sum / len(phase_records)
    annual_eff_rate = phase_eff_rate_sum / len(phase_records)
    annual_eff_sigma = np.sqrt(phase_eff_rate_var_sum) / len(phase_records)
    annual_rate_sem = float(np.std([
        x["metadata"]["integral"] for x in phase_records], ddof=1) /
                            np.sqrt(len(phase_records)))
    unit_hist = phase_hist / phase_hist.sum()
    base_hist = np.histogram(base["energies"], bins=energy_edges)[0].astype(float)
    base_hist /= base_hist.sum()
    finite = unit_hist > 0
    expected = len(base["energies"]) * unit_hist[finite]
    observed = len(base["energies"]) * base_hist[finite]
    energy_bin_z_max = float(np.max(np.abs(observed - expected) /
                                    np.sqrt(expected)))
    annual = {
        "construction": "12 equally weighted midpoint phases; 40 GeV/c2 diagnostic only",
        "native_rate_mean": annual_rate,
        "native_rate_phase_sem_including_isotope_noise": annual_rate_sem,
        "efficiency_folded_rate_mean": annual_eff_rate,
        "efficiency_folded_rate_mc_sigma": annual_eff_sigma,
        "phases": phase_records,
        "representative_minus_annual_relative_rate":
            base["metadata"]["integral"] / annual_rate - 1.0,
        "representative_minus_annual_efficiency_folded_relative_rate":
            base_stats["efficiency_folded_rate"] / annual_eff_rate - 1.0,
        "energy_histogram_max_bin_z_approx": energy_bin_z_max,
        "energy_histogram_bin_z_is_diagnostic_not_a_shape_fit": True,
    }
    # Exact representative-annual energy-mean residual, using the phase-specific
    # native-rate weights (not an arbitrary arithmetic mean of event means).
    annual_energy_mean = 0.0
    for item in phase_records:
        annual_energy_mean += (item["metadata"]["integral"] *
                               item["selection"]["energy_mean_kevnr"])
    annual_energy_mean /= phase_rate_sum
    annual["representative_minus_annual_energy_mean_kevnr"] = (
        base_stats["energy_mean_kevnr"] - annual_energy_mean)
    annual["annual_energy_mean_kevnr"] = annual_energy_mean
    rate_err = float(np.sqrt(annual_rate_sem ** 2 + np.std([
        x["integral"] for x in report["independent_seed_preparations"]], ddof=1) ** 2)
                    / annual_rate)
    eff_err = float(np.sqrt((annual_eff_sigma / annual_eff_rate) ** 2 +
                            (base_stats["efficiency_folded_rate_mc_sigma"] /
                             base_stats["efficiency_folded_rate"]) ** 2))
    annual["relative_rate_comparison_uncertainty"] = rate_err
    annual["relative_efficiency_folded_rate_comparison_uncertainty"] = eff_err
    annual["representative_approximation_within_3sigma"] = bool(
        abs(annual["representative_minus_annual_relative_rate"]) <=
        max(3.0 * rate_err, 1e-4) and
        abs(annual["representative_minus_annual_efficiency_folded_relative_rate"])
        <= 3.0 * eff_err)
    report["annual_average_check"] = annual
    report["interpretation"] = (
        "Annual modulation is not modeled as a search observable. The direct annual "
        "mixture is checked only at the pilot mass; full response/Fig.4 shape is "
        "checked separately before selecting the representative-day approximation.")

    out = Path(args.native_probe)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in (
        "representative_day", "baseline", "native_integral_vs_independent_quadrature",
        "above_100_keVnr", "same_seed_original_sampler_identical")}, indent=2))
    print("grid summaries:", json.dumps(grids, indent=2))
    print("annual summaries:", json.dumps(
        {k: v for k, v in annual.items() if k != "phases"}, indent=2))
    print(f"wrote {out}")
    return 0


class TestWIMPNormalization(unittest.TestCase):
    def test_representative_day_zeroes_periodic_term(self):
        self.assertAlmostEqual(earth_speed_km_s(REPRESENTATIVE_DAY), 250.2, places=12)
        self.assertNotAlmostEqual(earth_speed_km_s(0.0), 250.2, places=3)

    def test_exact_absolute_scaling(self):
        ref = expected_yield(12345.0, 3.3)
        self.assertEqual(expected_yield(12345.0, 3.3 * 2.0), 2.0 * ref)
        self.assertEqual(expected_yield(12345.0, 3.3, 2.0 * REFERENCE_SIGMA_CM2),
                         2.0 * ref)
        self.assertEqual(expected_yield(12345.0, 3.3, 0.0), 0.0)

    def test_phase_midpoints(self):
        days = phase_midpoints(12)
        self.assertAlmostEqual(days[0], 365.24 / 24)
        self.assertAlmostEqual(days[-1] - days[0], 365.24 * 11 / 12)

    def test_metadata_wrapper_preserves_sampler(self):
        seed, n, mass = 17000912, 100, 40.0
        nestpy.RandomGen.rndm().set_seed(seed)
        original = nestpy.spectra.WIMP(n, mass=mass, day=REPRESENTATIVE_DAY)
        wrapped = sample_native(mass, n, seed)
        self.assertTrue(np.array_equal(original, wrapped["energies"]))
        self.assertEqual(wrapped["metadata"]["n_samples"], n)
        self.assertGreater(wrapped["metadata"]["integral"], 0)
        with self.assertRaises(ValueError):
            sample_native(mass, -1, seed)

    def test_tail_sampler_and_weights_are_consistent(self):
        # Exact uniform analytic control, separately from NEST's physics.
        grid = np.linspace(100, 101, 101)
        diagnostic = {"tail_extension_used": True, "tail_energy_kevnr": grid,
                      "tail_native_rate_boundary_kevnr": [100.0, 101.0],
                      "tail_differential_rate": np.full_like(grid, 2.0)}
        e = sample_native_tail(diagnostic, 20000, 81828)
        self.assertTrue(np.all(e >= 100) and np.all(e < 101))
        self.assertAlmostEqual(np.mean(e), 100.5, delta=0.005)
        weights = native_tail_weight(diagnostic, e, 9.0 / len(e))
        # Mixture probabilities encode rates; do not weight dR/dE twice.
        self.assertAlmostEqual(float(weights.sum()), 9.0, places=10)
        with self.assertRaises(ValueError):
            sample_native_tail({"tail_extension_used": True,
                                "tail_energy_kevnr": grid,
                                "tail_differential_rate": -grid}, 1, 1)

    def test_linear_tail_mixture_closure(self):
        # Integral = (1+3)/2*2 + (3+5)/2*2 = 12; mean offset 22/9.
        grid = np.array([100.0, 102.0, 104.0])
        density = np.array([1.0, 3.0, 5.0])
        n = 60000
        diagnostic = {"tail_extension_used": True, "tail_energy_kevnr": grid,
                      "tail_native_rate_boundary_kevnr": [100.0, 104.0],
                      "tail_differential_rate": density}
        prepared_rate, tail_rate = 4.0, 12.0
        n_tail = np.random.default_rng(9182).binomial(
            n, tail_rate / (prepared_rate + tail_rate))
        energies = sample_native_tail(diagnostic, n_tail, 81829)
        weight = (prepared_rate + tail_rate) / n
        weights = native_tail_weight(diagnostic, energies, weight)
        self.assertAlmostEqual(np.mean(energies), 100 + 22 / 9, delta=0.02)
        # Exact fixed-total mixture closure for every tail-count realization.
        total = (n - n_tail) * weight + float(weights.sum())
        self.assertAlmostEqual(total, prepared_rate + tail_rate, places=10)

    def test_native_drate_is_exposed(self):
        nestpy.RandomGen.rndm().set_seed(17000913)
        value = nestpy.spectra.WIMP_dRate(5.0, 40.0, REPRESENTATIVE_DAY)
        self.assertTrue(np.isfinite(value))
        self.assertGreater(value, 0.0)
        quad = integrate_native(40.0, REPRESENTATIVE_DAY, np.linspace(0, 1, 101),
                                17000914)
        self.assertGreater(quad["integral"], 0.0)
class TestFullWIMPSelection(unittest.TestCase):
    """Small controls; the separate 40 GeV pilot supplies physical precision."""

    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parent.parent
        cls.engine = ResponseEngine(
            cls.root / "data/lz_ws2024_nest_header.hh", 182773)

    def setUp(self):
        # These helpers are functions: storing them as class attributes makes
        # Python bind them as methods receiving the test instance. Keep them
        # on the test instance instead.
        self.eff = load_efficiency(
            self.root / "data/lz_ws2024_fig2_efficiency.csv",
            "pre_roi_ss_data_cuts")
        self.band = load_band(
            self.root / "results/baseline/response_model_calibration.hdf5")

    def args(self, mode, n=1200):
        return argparse.Namespace(
            time_mode=mode, n_events=n, min_native_events=n,
            min_roi_events=1, max_blocks=1, e_step=5.0,
            n_grid_events=100, support_treatment="native_drate_extension",
            exposure=3.3, band=self.band, u_edges=np.linspace(-6.0, 4.0, 19),
            nest_parameter_set="ws2024_nest240", detector_effects="uniform_1d",
            response_model="nest_library",
            efficiency_treatment="input_pre_roi_curve", halo_model="shm_standard",
            exposure_structure="combined_single")

    def check_full_selection(self, m):
        self.assertTrue(np.all(m["s1c"] > 3))
        self.assertTrue(np.all(m["s1c"] < 80))
        self.assertTrue(np.all(m["s2c"] >= 493))
        self.assertTrue(np.all(m["s2c"] < 10**4.5))
        self.assertTrue(np.all(roi_mask(m["s1c"], m["s2c"])))
        self.assertEqual(m["n_native"], sum(len(r["energies"]) for r in m["records"]))
        self.assertEqual(int(m["hist_mc"].sum()), m["n_roi"])
        self.assertEqual(len(m["s1c"]), m["n_roi"])
        self.assertAlmostEqual(m["hist"].sum(), m["absolute_yield"])
        self.assertAlmostEqual(m["hist_var"].sum(), m["mc_variance"])
        edges = np.linspace(-6.0, 4.0, 19)
        mu, sigma = self.band(m["s1c"])
        u = (m["y"] - mu) / sigma
        for panel in m["u_hist"]:
            lo, hi = {"top": (3, 40), "bottom": (40, 80)}[panel]
            selected = (m["s1c"] > lo) & (m["s1c"] <= hi)
            np.testing.assert_array_equal(
                m["u_mc"][panel], np.histogram(u[selected], bins=edges)[0])
            np.testing.assert_allclose(
                m["u_hist"][panel],
                np.histogram(u[selected], bins=edges, weights=m["weights"][selected])[0],
                rtol=1e-11, atol=0)
        np.testing.assert_allclose(m["absolute_yield"], m["accepted_yield_closure"],
                                   rtol=1e-11, atol=0.0)
        for r in m["records"]:
            self.assertGreater(r["weight"], 0)
            self.assertTrue(np.all(r["eff_energies"] > 0))
            self.assertEqual(len(r["native_weights"]), len(r["eff_energies"]))
            self.assertEqual(len(r["native_weights"]), int(r["selected"].sum()))
            self.assertEqual(len(r["selected"]), len(r["energies"]))
        return m

    def test_annual_phase_weights_sum_to_one_year(self):
        m = self.check_full_selection(generate_mass(
            40.0, self.args("annual_12phase"), 17373001, self.engine, self.eff))
        self.assertEqual(len(m["records"]), 12)
        self.assertEqual(m["n_native"], 1200)
        self.assertTrue(all(r["phase_share"] == 1.0 / 12 for r in m["records"]))
        self.assertAlmostEqual(sum(r["phase_yield"] for r in m["records"]),
                               m["reference_yield"], places=3)
        self.assertLess(0.2, m["roi_acceptance"])
        self.assertLess(m["roi_acceptance"], 0.8)

    def test_representative_day_full_selection_roundtrip(self):
        args = self.args("representative_day")
        m = self.check_full_selection(generate_mass(
            40.0, args, 17373002, self.engine, self.eff))
        # The complete own-seed chain must reproduce, not merely reclose.
        self.assertFalse(m["support_diagnostics"]["tail_extension_used"])
        self.assertEqual(m["n_tail_native_generated"], 0)
        again = generate_mass(40.0, args, 17373002, self.engine, self.eff)
        self.assertTrue(np.array_equal(m["hist_mc"], again["hist_mc"]))
        self.assertTrue(np.array_equal(m["u_mc"]["top"], again["u_mc"]["top"]))
        self.assertTrue(np.array_equal(m["s1c"], again["s1c"]))
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "wimp.hdf5"
            with h5py.File(path, "w") as f:
                g = write_model(f, None, m, args)
                self.assertAlmostEqual(g["hist_counts"][:].sum(),
                                       g.attrs["accepted_yield_direct_closure"])
                self.assertEqual(int(g["hist_mc_counts"][:].sum()), m["n_roi"])
                self.assertAlmostEqual(g["weight"][:].sum(),
                                       g.attrs["absolute_roi_yield_at_reference_sigma"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-probe", help="Write detailed native JSON here")
    parser.add_argument("--skip-annual", action="store_true",
                        help="Support/grid probes at non-pilot masses need no modulation work")
    parser.add_argument("--mass", type=float, default=40.0)
    parser.add_argument("--seed", type=int, default=17000071)
    parser.add_argument("--n-events", type=int, default=120000)
    parser.add_argument("--grid-events", type=int, default=30000)
    parser.add_argument("--phase-events", type=int, default=10000)
    parser.add_argument("--e-steps", type=float, nargs="+",
                        default=[2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0])
    parser.add_argument("--efficiency",
                        default="data/lz_ws2024_fig2_efficiency.csv")
    args = parser.parse_args()
    if args.native_probe:
        sys.exit(native_probe(args))
    unittest.main(argv=[sys.argv[0]], verbosity=2)
