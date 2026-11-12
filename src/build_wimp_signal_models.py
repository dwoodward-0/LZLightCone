"""Build absolute WS2024 WIMP templates from NEST's native WIMP spectrum.

Chain: native true recoil energies -> public pre-ROI NR efficiency rejection ->
validated ResponseEngine(nr=True) -> exactly one corrected-observable ROI ->
absolute native-rate yield. The archive is at NEST's native 1e-36 cm2 reference
cross section and 3.3 tonne-years WS2024 exposure. There is no Table-I or
figure normalization. A direct 12-phase annual comparison is made only at
40 GeV/c2, never as a modulation measurement or a per-mass production grid.

The later 4.2-tonne-year combined model must add a separately validated
0.9-tonne-year WS2022 campaign, with its response/selection/ROI. It must not
be made by scaling this archive or its backgrounds by 4.2/3.3.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_background_models import (
    HIST_S1C_EDGES, HIST_Y_EDGES, MAX_SAMPLES, ROI_S1C_MIN,
    ResponseEngine, load_efficiency, roi_mask,
)
from validate_backgrounds import PANEL_S1C, load_band, load_reference
from wimp_native import (
    NATIVE_INTEGRAL_EMAX_KEV, REFERENCE_SIGMA_CM2, REPRESENTATIVE_DAY,
    expected_yield, integrate_native, native_provenance, native_support_diagnostics,
    native_tail_weight, phase_midpoints, sample_native, sample_native_tail,
)


def read_mass_grid(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows or "mass_gev" not in rows[0] or "limit" not in rows[0]:
        raise ValueError("SI reference must contain mass_gev and observed limit")
    if len(rows) != 26 or len({float(r["mass_gev"]) for r in rows}) != 26:
        raise ValueError("The pinned SI reference must define exactly 26 unique masses")
    if not all(float(r["mass_gev"]) > 0 and float(r["limit"]) > 0 for r in rows):
        raise ValueError("Published mass/reference-limit table must be positive")
    return rows


def source_seed(base_seed, mass_rows, mass):
    """Stable in the external published grid, independent of requested masses."""
    index = next(i for i, row in enumerate(mass_rows)
                 if float(row["mass_gev"]) == float(mass))
    return int(base_seed + 1000003 * index)


def selected_counts(s1c, y):
    mc, _, _ = np.histogram2d(s1c, y, bins=[HIST_S1C_EDGES, HIST_Y_EDGES])
    bin_index = np.nonzero(mc)
    counts = mc[bin_index]
    fractions = counts.astype(float)
    return bin_index, counts, fractions


def generate_mass(mass, args, mass_seed, engine, efficiency):
    """One independent native/response realization; never table-normalized."""
    support_cache = {}
    kind = args.time_mode
    days = [REPRESENTATIVE_DAY] if kind == "representative_day" else list(
        phase_midpoints(12))
    n_blocks = 0
    n_native_total = 0
    records = []
    while n_native_total < args.min_native_events or (
            sum(int(np.sum(roi_mask(r["s1c"], r["s2c"]))) for r in records)
            < args.min_roi_events) and n_blocks < args.max_blocks:
        for phase_index, day in enumerate(days):
            block_seed = mass_seed + 10007 * (phase_index + 1) + (
                n_blocks * 1_000_003 if kind == "annual_12phase" else n_blocks)
            # Use the same total initial budget for one- and 12-phase models.
            start = n_native_total if kind == "representative_day" else (
                n_blocks * args.n_events + phase_index *
                (args.n_events // len(days)))
            stop = n_native_total + args.n_events if kind == "representative_day" else (
                (n_blocks + 1) * args.n_events if phase_index == len(days) - 1 else
                n_blocks * args.n_events + (phase_index + 1) *
                (args.n_events // len(days)))
            n = stop - start
            prep_probe = sample_native(mass, 0, block_seed, day=day, e_step=args.e_step)
            if day not in support_cache:
                support_cache[day] = native_support_diagnostics(
                    mass, day, block_seed, efficiency,
                    prep_probe["metadata"]["integral"],
                    n_grid_events=getattr(args, "n_grid_events", 30000))
            support = support_cache[day]
            n_tail = (int(np.random.default_rng(block_seed + 119).binomial(
                n, support["native_integral_100_500"] /
                (prep_probe["metadata"]["integral"] +
                 support["native_integral_100_500"])))
                if support["tail_extension_used"] else 0)
            # No tail: retain the validated native sample stream byte-for-byte.
            sample = sample_native(mass, n - n_tail, block_seed,
                                   day=day, e_step=args.e_step)
            e = np.concatenate([sample["energies"], sample_native_tail(
                support, n_tail, block_seed + 121)])
            meta = sample["metadata"]
            tail_rate = (support["native_integral_100_500"]
                         if support["tail_extension_used"] else 0.0)
            meta["native_added_tail_integral"] = float(tail_rate)
            meta["total_integral_including_sampled_tail"] = meta["integral"] + tail_rate
            meta["n_native_grid_samples"] = len(sample["energies"])
            meta["n_tail_samples"] = n_tail
            meta["n_samples"] = n
            meta["n_zero_samples"] = int(np.sum(e == 0))
            meta["n_negative_samples"] = int(np.sum(e < 0))
            meta["n_nonfinite_samples"] = int(np.sum(~np.isfinite(e)))
            meta["n_above_native_integral_support"] = int(np.sum(e >= 100))
            meta["n_above_xmax"] = int(np.sum(e > meta["xMax"]))
            meta["energy_min_kevnr"] = float(np.min(e))
            meta["energy_max_kevnr"] = float(np.max(e))
            meta["support_upper_kevnr"] = 500.0 if n_tail else 100.0
            weight = expected_yield(
                meta["total_integral_including_sampled_tail"], args.exposure) / len(days) / n
            native_weights = np.full(n, weight)
            if n_tail:
                native_weights[-n_tail:] = native_tail_weight(
                    support, e[-n_tail:], weight)
            sample["metadata"] = meta
            good = e > 0.0  # native sampler's failed attempts cannot generate light
            rng = np.random.default_rng(block_seed + 113)
            epsilon = efficiency(np.maximum(e, 0.0))
            epsilon[~good] = 0.0
            selected = rng.uniform(size=n) < epsilon
            eff_energies = e[selected]
            engine.reseed(block_seed + 227)
            s1c, s2c = engine.generate(eff_energies, nr=True)
            if np.any(eff_energies == 0):
                raise AssertionError("Failed native sample entered NEST")
            y = np.log10(np.maximum(s2c, 1e-300))
            meta = sample["metadata"]
            meta["n_above_efficiency_support"] = int(np.sum(
                e > efficiency.energy_kevnr[-1]))
            meta["n_efficiency_selected_above_efficiency_support"] = int(np.sum(
                eff_energies > efficiency.energy_kevnr[-1]))
            meta["efficiency_stage"] = efficiency.stage
            meta["efficiency_acceptance"] = float(np.sum(selected) / n)
            meta["efficiency_binomial_sigma"] = float(np.sqrt(
                meta["efficiency_acceptance"] *
                (1 - meta["efficiency_acceptance"]) / n))
            phase_yield = expected_yield(
                meta["total_integral_including_sampled_tail"], args.exposure) / len(days)
            records.append({
                "meta": meta, "phase_yield": phase_yield, "phase_share": 1.0 / len(days),
                "energies": e, "eff_energies": eff_energies,
                "native_weights": native_weights[selected], "selected": selected,
                "s1c": s1c, "s2c": s2c, "y": y, "weight": weight, "support": support,
                "phase": phase_index, "seed": block_seed,
                "outside_support": eff_energies > efficiency.energy_kevnr[-1],
            })
            n_native_total += n
        n_blocks += 1

    if not records or not sum(len(r["eff_energies"]) for r in records):
        raise RuntimeError(f"{mass} GeV/c2 has no post-efficiency events")
    n_roi = sum(int(np.sum(roi_mask(r["s1c"], r["s2c"]))) for r in records)
    if n_roi < args.min_roi_events:
        raise RuntimeError(f"{mass}: only {n_roi} ROI MC events after {n_blocks} blocks")

    hist = np.zeros((len(HIST_S1C_EDGES) - 1, len(HIST_Y_EDGES) - 1))
    hist_mc = np.zeros_like(hist)
    # Sum native multinomial variances of every independently seeded sample.
    # p is the probability of landing in an observable bin through all stages;
    # using n_selected rather than the phase's original n would inflate it.
    hist_var = np.zeros_like(hist)
    u_hist = {p: np.zeros(len(args.u_edges) - 1) for p in PANEL_S1C}
    u_mc = {p: np.zeros(len(args.u_edges) - 1) for p in PANEL_S1C}
    u_var = {p: np.zeros(len(args.u_edges) - 1) for p in PANEL_S1C}
    native_rate_numerator = 0.0
    for r in records:
        n = len(r["energies"])
        native_rate_numerator += n * r["meta"]["total_integral_including_sampled_tail"]
        accepted = roi_mask(r["s1c"], r["s2c"])
        hw = r["native_weights"][roi_mask(r["s1c"], r["s2c"])]
        hist_mc += np.histogram2d(r["s1c"][accepted], r["y"][accepted],
                                   bins=[HIST_S1C_EDGES, HIST_Y_EDGES])[0]
        hist_block = np.histogram2d(
            r["s1c"][accepted], r["y"][accepted],
            bins=[HIST_S1C_EDGES, HIST_Y_EDGES], weights=hw)[0]
        hist += hist_block
        p_all, n_all, _ = selected_counts(r["s1c"], r["y"])
        all_w = np.where(r["selected"], r["weight"], 0.0)
        all_w[r["selected"]] = r["native_weights"]
        selected_w = all_w[r["selected"]]
        sumw2 = np.histogram2d(
            r["s1c"], r["y"], bins=[HIST_S1C_EDGES, HIST_Y_EDGES],
            weights=selected_w ** 2)[0]
        sumw = np.histogram2d(
            r["s1c"], r["y"], bins=[HIST_S1C_EDGES, HIST_Y_EDGES],
            weights=selected_w)[0]
        # General importance-weighted multinomial variance. With one constant
        # weight this reduces exactly to the original w^2 n p (1-p).
        hist_var += sumw2 - sumw ** 2 / n
        if np.any(hist_var < -1e-8 * max(hist.max(), 1) ** 2):
            raise AssertionError("Negative full histogram MC variance")

        mu, sig = args.band(r["s1c"])
        u = (r["y"] - mu) / sig
        for panel, (lo, hi) in PANEL_S1C.items():
            on_panel = (r["s1c"] > lo) & (r["s1c"] <= hi)
            panel_w = selected_w[on_panel]
            sumw_panel = np.histogram(u[on_panel], bins=args.u_edges,
                                      weights=panel_w)[0]
            sumw2_panel = np.histogram(u[on_panel], bins=args.u_edges,
                                       weights=panel_w ** 2)[0]
            u_var[panel] += sumw2_panel - sumw_panel ** 2 / n
            accepted_panel = accepted & on_panel
            u_hist[panel] += np.histogram(
                u[accepted_panel], bins=args.u_edges, weights=hw[on_panel[accepted]])[0]
            u_mc[panel] += np.histogram(u[accepted_panel], bins=args.u_edges)[0]

    s1c = np.concatenate([r["s1c"][roi_mask(r["s1c"], r["s2c"])] for r in records])
    s2c = np.concatenate([r["s2c"][roi_mask(r["s1c"], r["s2c"])] for r in records])
    y = np.log10(s2c)
    energies = np.concatenate([
        r["eff_energies"][roi_mask(r["s1c"], r["s2c"])] for r in records])
    weights = np.concatenate([
        r["native_weights"][roi_mask(r["s1c"], r["s2c"])] for r in records])
    phases = np.concatenate([
        np.full(np.sum(roi_mask(r["s1c"], r["s2c"])), r["phase"], dtype=int)
        for r in records])
    if len(s1c) != int(hist_mc.sum()) or not np.isclose(hist.sum(), weights.sum()):
        raise AssertionError("Full observable histogram and retained events disagree")
    display_rng = np.random.default_rng(mass_seed + 300_000_007)
    if len(s1c) > MAX_SAMPLES:
        display = np.sort(display_rng.choice(len(s1c), MAX_SAMPLES, replace=False))
    else:
        display = np.arange(len(s1c))
    display_weights = weights[display].copy()
    display_weights *= weights.sum() / display_weights.sum()
    absolute_yield = float(hist.sum())
    native_rate = native_rate_numerator / n_native_total
    reference_yield = expected_yield(native_rate, args.exposure)
    mc_variance = float(hist_var.sum())
    closure_sigma = float(np.sqrt(mc_variance))
    accepted_yield_closure = float(weights.sum())
    # The histogram sums many large floating-point weights; a few ULPs of
    # relative cancellation must not be mistaken for a normalization change.
    if not np.isclose(accepted_yield_closure, absolute_yield,
                      rtol=1e-11, atol=0.0):
        raise AssertionError(
            f"Absolute yield does not close: {absolute_yield} vs "
            f"{accepted_yield_closure} (rel "
            f"{absolute_yield / accepted_yield_closure - 1:.3g})")
    if not np.isclose(expected_yield(native_rate * 2, args.exposure),
                      reference_yield * 2):
        raise AssertionError("Exposure-scaled native rate is not exactly linear")
    return {
        "mass": float(mass), "time_mode": kind, "days": days,
        "native_rate": native_rate, "reference_yield": reference_yield,
        "absolute_yield": absolute_yield, "mc_variance": mc_variance,
        "n_native": n_native_total, "n_roi": n_roi, "n_blocks": n_blocks,
        "source_seed": mass_seed, "records": records,
        "hist": hist, "hist_mc": hist_mc, "hist_var": hist_var,
        "u_hist": u_hist, "u_mc": u_mc, "u_var": u_var,
        "s1c": s1c, "s2c": s2c, "y": y, "energies": energies,
        "weights": weights, "display": display, "display_weights": display_weights,
        "phases": phases, "reference_yield_mc_closure": float(
            reference_yield * n_roi / n_native_total),
        "accepted_yield_closure": accepted_yield_closure,
        "roi_acceptance": float(n_roi / n_native_total),
        "roi_acceptance_relative_mc_sigma": closure_sigma / absolute_yield
        if absolute_yield else None,
        "support_diagnostics": next(iter(support_cache.values())),
        "support_diagnostics_by_phase": {str(k): v for k, v in support_cache.items()},
        "n_tail_native_generated": sum(r["meta"]["n_tail_samples"] for r in records),
        "n_tail_efficiency_selected": sum(int(np.sum(
            r["eff_energies"] >= NATIVE_INTEGRAL_EMAX_KEV)) for r in records),
        "n_tail_roi": sum(int(np.sum(
            r["eff_energies"][roi_mask(r["s1c"], r["s2c"])] >=
            NATIVE_INTEGRAL_EMAX_KEV)) for r in records),
        "tail_absolute_yield": float(sum(
            r["native_weights"][roi_mask(r["s1c"], r["s2c"])][
                r["eff_energies"][roi_mask(r["s1c"], r["s2c"])] >=
                NATIVE_INTEGRAL_EMAX_KEV].sum() for r in records)),
        "native_rate_preparation_relative_sd": float(np.std(
            [r["meta"]["total_integral_including_sampled_tail"]
             for r in records]) / native_rate)
        if len(records) > 1 else 0.0,
    }


def write_model(f, key, m, args, store_diagnostic=True):
    g = f.create_group(f"ws2024/mass_{m['mass']:g}gev" if key is None else key)
    g.attrs.update({
        "campaign": "ws2024", "exposure_tonne_years": args.exposure,
        "reference_sigma_cm2": REFERENCE_SIGMA_CM2,
        "mass_gev": m["mass"], "time_mode": m["time_mode"],
        "native_rate_events_kg_day": m["native_rate"],
        "preselection_yield_at_reference_sigma": m["reference_yield"],
        "absolute_roi_yield_at_reference_sigma": m["absolute_yield"],
        "roi_yield_mc_variance": m["mc_variance"],
        "roi_yield_mc_closure": m["reference_yield_mc_closure"],
        "accepted_yield_direct_closure": m["accepted_yield_closure"],
        "n_native_generated": m["n_native"], "n_full_roi_mc": m["n_roi"],
        "n_seeded_blocks": m["n_blocks"], "source_seed": m["source_seed"],
        "support_treatment": args.support_treatment,
        "n_tail_native_generated": m["n_tail_native_generated"],
        "n_tail_efficiency_selected": m["n_tail_efficiency_selected"],
        "n_tail_roi": m["n_tail_roi"],
        "tail_absolute_roi_yield_at_reference_sigma": m["tail_absolute_yield"],
        "native_tail_fraction_of_roi_yield":
            m["tail_absolute_yield"] / m["absolute_yield"],
        "efficiency_treatment": args.efficiency_treatment,
        "response_model": args.response_model,
        "exposure_structure": args.exposure_structure,
        "efficiency_stage": "pre_roi_ss_data_cuts",
        "nest_parameter_set": args.nest_parameter_set,
        "detector_effects": args.detector_effects,
        "halo_model": args.halo_model,
        "native_spectrum": "NEST TestSpectra WIMP_with_metadata; isotope and Helm native",
        "roi": "S1c>3; S1c<80; S2c>=493 phd; log10(S2c)<4.5; applied once",
        "normalization": "native integral *1000*365.25*3.3; no Table-I/figure normalization",
        "variance": "full native multinomial sum(w2)-sum(w)^2/N per block before thinning",
        "combined_exposure": "not 4.2 tonne-years; explicit future WS2022 component required",
    })
    g.create_dataset("s1c", data=m["s1c"][m["display"]])
    g.create_dataset("log10_s2c", data=m["y"][m["display"]])
    g.create_dataset("recoil_energy_kevnr", data=m["energies"][m["display"]])
    g.create_dataset("weight", data=m["display_weights"])
    if m["time_mode"] != "representative_day" or store_diagnostic:
        g.create_dataset("phase_index", data=m["phases"][m["display"]])
    for name, data in (("hist_counts", m["hist"]), ("hist_mc_counts", m["hist_mc"]),
                       ("hist_mc_weight_variance", m["hist_var"]),
                       ("hist_s1c_edges", HIST_S1C_EDGES),
                       ("hist_y_edges", HIST_Y_EDGES),
                       ("u_edges", args.u_edges)):
        g.create_dataset(name, data=data)
    for panel, (lo, hi) in PANEL_S1C.items():
        pg = g.create_group(f"fig4_{panel}")
        pg.attrs["s1c_lower_exclusive_phd"] = lo
        pg.attrs["s1c_upper_inclusive_phd"] = hi
        pg.create_dataset("counts", data=m["u_hist"][panel])
        pg.create_dataset("mc_counts", data=m["u_mc"][panel])
        pg.create_dataset("mc_variance", data=m["u_var"][panel])
    sg = g.create_group("native_support_diagnostics")
    for day, diagnostics in m["support_diagnostics_by_phase"].items():
        dg = sg.create_group(f"day_{float(day):.6f}")
        arrays = {k: diagnostics[k] for k in (
            "tail_energy_kevnr", "tail_differential_rate")}
        compact = {k: v for k, v in diagnostics.items() if k not in arrays}
        dg.create_dataset("summary_json", data=json.dumps(compact),
                          dtype=h5py.string_dtype("utf-8"))
        dg.create_dataset("tail_energy_kevnr", data=arrays["tail_energy_kevnr"])
        dg.create_dataset("tail_differential_rate", data=arrays["tail_differential_rate"])
    ng = g.create_group("native_preparations")
    for i, r in enumerate(m["records"]):
        rg = ng.create_group(f"block_{i:03d}_phase_{r['phase']:02d}")
        rg.attrs.update({k: v for k, v in r["meta"].items() if v is not None})
        rg.attrs["phase_fraction_of_annual_native_yield"] = r["phase_share"]
        rg.attrs["phase_preselection_native_reference_yield"] = r["phase_yield"]
        rg.attrs["response_seed"] = r["seed"] + 227
        rg.create_dataset("generated_energies", data=r["energies"])
        rg.create_dataset("efficiency_selected_energies", data=r["eff_energies"])
        rg.create_dataset("roi_recoil_energies", data=r["eff_energies"][
            roi_mask(r["s1c"], r["s2c"])])
        rg.create_dataset("outside_efficiency_support", data=r["outside_support"])
    return g


def compare_phases(primary, annual):
    """Independent full-MC comparison; no fitting to the LZ reference."""
    relative = primary["absolute_yield"] / annual["absolute_yield"] - 1
    unc = np.sqrt(primary["mc_variance"] / primary["absolute_yield"] ** 2 +
                  annual["mc_variance"] / annual["absolute_yield"] ** 2)
    shapes = {}
    for panel in PANEL_S1C:
        pp = primary["u_hist"][panel] / primary["absolute_yield"]
        ap = annual["u_hist"][panel] / annual["absolute_yield"]
        pv = primary["u_var"][panel] / primary["absolute_yield"] ** 2
        av = annual["u_var"][panel] / annual["absolute_yield"] ** 2
        shapes[panel] = {
            "bins": float(np.sum(pp > 0)),
            "max_abs_unit_shape_difference": float(np.max(np.abs(pp - ap))),
            "max_abs_difference_mc_sigma": float(np.max(np.divide(
                np.abs(pp - ap), np.sqrt(pv + av),
                out=np.zeros_like(pp), where=(pv + av) > 0))),
        }
    return {
        "relative_absolute_yield_difference_representative_minus_annual": relative,
        "relative_yield_difference_mc_sigma": float(unc),
        "representative_within_3_mc_sigma_of_annual": bool(abs(relative) <= 3 * unc),
        "shape_differences": shapes,
        "note": "Independent native/response realizations; no common-RNG cancellation or LZ tuning",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("header", "efficiency", "response", "sensitivity", "fig4_reference", "output"):
        ap.add_argument(f"--{name.replace('_', '-')}", required=True)
    ap.add_argument("--masses", type=float, nargs="+", help="Default: published 26-mass grid")
    ap.add_argument("--exposure", type=float, default=3.3)
    ap.add_argument("--seed", type=int, default=17000071)
    ap.add_argument("--n-events", type=int, default=120000)
    ap.add_argument("--min-native-events", type=int, default=120000)
    ap.add_argument("--min-roi-events", type=int, default=2000)
    ap.add_argument("--max-blocks", type=int, default=20)
    ap.add_argument("--e-step", type=float, default=5.0)
    ap.add_argument("--support-treatment", required=True,
                    choices=("native_drate_extension",))
    ap.add_argument("--n-grid-events", type=int, default=30000,
                    help="Independent sampler/grid controls per mass")
    ap.add_argument("--time-mode", choices=("representative_day", "annual_12phase"),
                    default="representative_day")
    ap.add_argument("--annual-comparison", action="store_true",
                    help="At 40 GeV/c2 only, independently generate the annual control")
    ap.add_argument("--efficiency-treatment", required=True,
                    choices=("input_pre_roi_curve",))
    ap.add_argument("--response-model", required=True, choices=("nest_library",))
    ap.add_argument("--nest-parameter-set", required=True, choices=("ws2024_nest240",))
    ap.add_argument("--detector-effects", required=True, choices=("uniform_1d",))
    ap.add_argument("--exposure-structure", required=True, choices=("combined_single",))
    ap.add_argument("--halo-model", required=True, choices=("shm_standard",))
    args = ap.parse_args()
    if args.exposure != 3.3:
        raise ValueError("This archive is explicitly the validated 3.3-t-y WS2024 chain")
    if args.time_mode == "annual_12phase" and args.annual_comparison:
        raise ValueError("Annual comparison should diagnose the representative-day model")
    if args.annual_comparison and (args.masses is not None and
                                 list(map(float, args.masses)) != [40.0]):
        raise ValueError("Annual response comparison is restricted to the 40 GeV pilot")
    rows = read_mass_grid(args.sensitivity)
    masses = args.masses if args.masses is not None else [float(r["mass_gev"]) for r in rows]
    for mass in masses:
        if not any(float(r["mass_gev"]) == mass for r in rows):
            raise ValueError(f"Mass {mass} is not from the published SI sensitivity grid")
    ref = load_reference(args.fig4_reference)
    tot = ref[("top", "total_bkg")]
    args.u_edges = float(tot[0, 1]) + float(np.median(np.diff(tot[:, 1]))) * np.arange(len(tot) + 1)
    args.band = load_band(args.response)
    eff = load_efficiency(args.efficiency, "pre_roi_ss_data_cuts")
    engine = ResponseEngine(args.header, args.seed)
    limits = {float(r["mass_gev"]): float(r["limit"]) for r in rows}
    selected = {}
    output = Path(args.output)
    with h5py.File(output, "w") as f:
        f.attrs.update({
            "campaign_scope": "WS2024 only, 3.3 tonne-years central exposure",
            "reference_sigma_cm2": REFERENCE_SIGMA_CM2,
            "reference_sigma_meaning": "NEST-native SI per-nucleon reference, not LZ limit",
            "comparison_limit_40gev_cm2": limits[40.0],
            "comparison_median_40gev_cm2": float(next(
                r["median"] for r in rows if float(r["mass_gev"]) == 40.0)),
            "package_and_patch_json": json.dumps(native_provenance()),
            "native_support_treatment": args.support_treatment,
            "native_support_definition": (
                "native prepared 0-100 keV integral + native-dRate linearly-interpolated "
                "100-500 keV tail when efficiency-folded fraction >1e-4; 500-5000 guard "
                "must be <1e-4; no yield rescaling to figures"),
            "annual_convention": (
                "NEST scalar v_e; representative day zeroes periodic term; not day=0; "
                "direct 12-phase annual diagnostic at 40 GeV/c2 only"),
            "historical_decision_guard_dex_metrics": "not used for sparse WIMP reference",
        })
        for mass in masses:
            seed = source_seed(args.seed, rows, mass)
            print(f"Generating {mass:g} GeV/c2: seed {seed}, mode {args.time_mode}")
            model = generate_mass(mass, args, seed, engine, eff)
            write_model(f, None, model, args)
            selected[mass] = model
            print(json.dumps({k: model[k] for k in (
                "native_rate", "reference_yield", "absolute_yield", "n_native", "n_roi",
                "roi_acceptance", "roi_acceptance_relative_mc_sigma")}, indent=2))
            if mass == 40.0 and args.annual_comparison:
                original_mode = args.time_mode
                args.time_mode = "annual_12phase"
                annual = generate_mass(mass, args, seed + 500_000_000, engine, eff)
                args.time_mode = original_mode
                comparison = compare_phases(model, annual)
                write_model(f, "diagnostics/annual_12phase_40gev", annual, args)
                dg = f.create_group("diagnostics/representative_vs_annual_40gev")
                dg.create_dataset("summary", data=json.dumps(comparison),
                                  dtype=h5py.string_dtype("utf-8"))
                print("Annual response comparison:", json.dumps(comparison, indent=2))
                if not comparison["representative_within_3_mc_sigma_of_annual"]:
                    print("FAIL: representative-day absolute yield differs from direct annual average")
                    return 2
        f.attrs["native_metadata_json"] = json.dumps(
            {str(mass): [r["meta"] for r in model["records"]]
             for mass, model in selected.items()})
        f.attrs["n_generated_masses"] = len(selected)
        f.attrs["completed"] = True
    print(f"PASS: wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
