"""Complete matched-budget comparison of the eight ER background sources.

The staged and legacy models regenerate all twelve sources with independent
per-source random streams. NR sources and observable-space accidentals must
be identical between ER treatments. Table-I central physical counts are fixed;
retained MC multiplicity, effective sample size and shapes can change.

The archive contains both full models, all stage counts, full/display weights,
MC-aware unit-area S1c/log10-S2c/u comparisons, and both models' six-gate
Validation 4 metrics. Independent forward checks against cumulative Y and
HEPData total-2D t4 are diagnostic only, never another selection factor.
"""
import argparse
import json
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import numpy as np

from build_background_models import SOURCE_SPECS, load_efficiency, make_parser
from check_er_efficiency_forward import forward_checks
from compare_nr_efficiency import (
    NR_SOURCES, copy_source, dump_metadata, generate, json_value,
    record_shapes, same_model, validate)
from er_efficiency import load_staged_er
from validate_backgrounds import load_band

MODES = {
    "legacy_roi_only": "legacy_no_energy_efficiency",
    "staged_conditional": "staged_conditional",
}
ER_SOURCES = tuple(source for source, (kind, _) in SOURCE_SPECS.items()
                   if kind in ("curve", "analytic"))
UNAFFECTED_SOURCES = NR_SOURCES + ("accidental_coincidences",)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for action in make_parser()._actions:
        if action.dest == "help":
            continue
        options = {key: getattr(action, key) for key in (
            "default", "type", "choices", "required", "nargs", "const",
            "metavar", "help")}
        p.add_argument(*action.option_strings, **options)
    p.add_argument("--background-models", required=True,
                   help="Full model generated under the active ER decision")
    p.add_argument("--response", required=True, help="ER-band calibration HDF5")
    p.add_argument("--validation-script", required=True)
    p.add_argument("--reference", required=True, help="Digitized Fig. 4 CSV")
    p.add_argument("--er-2d-efficiency", required=True,
                   help="Independent HEPData t4 comparison only")
    p.add_argument("--forward-check-statistics", type=int, default=5000,
                   help="Independent NEST fixed-energy MC samples per 0.25-keV point")
    p.add_argument("--dark", action="store_true", help=argparse.SUPPRESS)
    p.add_argument(
        "--allow-reduced-budget-validation-failure", action="store_true",
        help="Probe only: archive finite-MC Validation-4 failures below the "
             "production budgets. The declared recipe never uses this flag.")
    return p


def safe_acceptance(numerator, denominator):
    p = numerator / denominator if denominator else None
    sigma = np.sqrt(p * (1 - p) / denominator) if p is not None else None
    return None if p is None else float(p), None if sigma is None else float(sigma)


def source_metrics(group, mode):
    a = group.attrs
    full_mc = int(group["hist_mc_counts"][:].sum())
    weights = group["hist_counts"][:].sum()
    display = group["weight"][:].sum()
    assert full_mc == a["n_in_roi"]
    assert np.isclose(weights, a["table1_mean"], atol=1e-9)
    assert np.isclose(display, a["table1_mean"], atol=1e-9)
    variance = group["hist_mc_weight_variance"][:].sum()
    n_conditional = (int(a["n_after_er_conditional"])
                     if mode == "staged_conditional" else None)
    trigger_count = (int(a["n_er_after_trigger"])
                     if mode == "staged_conditional" else None)
    s1_mc = int(a["n_after_s1_stage"])
    expected_y = float(a["er_y_true_energy_mean"])
    y_ratio = s1_mc / a["n_er_after_trigger"] if a["n_er_after_trigger"] else None
    return {
        "source_seed": int(a["source_seed"]),
        "er_stage_applicability": str(a["er_efficiency_stage"]),
        "table1_central_physical_count": float(a["table1_mean"]),
        "table1_count_uncertainty": float(a["table1_sigma"]),
        "n_true_energy_samples": int(a["n_generated"]),
        "n_after_trigger": trigger_count,
        "n_negative_s1_nest_coincidence_flags": int(a["n_negative_s1_nest_coincidence_flags"]),
        "n_after_native_coincidence_and_s1": s1_mc,
        "n_after_conditional_g_over_y": n_conditional,
        "full_roi_mc_count": full_mc,
        "retained_display_count": len(group["s1c"]),
        "full_histogram_physical_count": float(weights),
        "display_physical_count": float(display),
        "full_histogram_effective_sample_size": float(weights**2 / variance) if variance > 0 else None,
        "trigger_acceptance": (safe_acceptance(trigger_count, a["n_generated"])
                               if trigger_count is not None else None),
        "s1_stage_acceptance_conditional_on_trigger": (
            safe_acceptance(s1_mc, trigger_count)
            if trigger_count is not None else None),
        "g_over_y_stage_acceptance_conditional_on_s1": (
            safe_acceptance(n_conditional, s1_mc)
            if n_conditional is not None else None),
        "remaining_roi_acceptance_after_energy_selection": safe_acceptance(
            full_mc, n_conditional if n_conditional is not None else s1_mc),
        "overall_simulation_acceptance": safe_acceptance(full_mc, a["n_generated"]),
        "cumulative_y_expected_in_trigger_sample": expected_y if np.isfinite(expected_y) else None,
        "cumulative_y_expected_binomial_sigma": float(a["er_y_true_energy_binomial_sigma"]),
        "cumulative_y_mc_measured": y_ratio,
        "cumulative_y_mc_minus_true_energy_prediction": y_ratio - expected_y if y_ratio is not None and np.isfinite(expected_y) else None,
        "generated_outside_trigger_support": int(a["n_er_generated_outside_trigger_support"]),
        "generated_outside_g_support": int(a["n_er_generated_outside_g_support"]),
        "g15_plateau_post_conditional_count": int(a["n_er_g15_plateau_extrapolated_in_sample"]),
        "g15_plateau_in_roi_count": int(a["n_er_g15_plateau_extrapolated_in_roi"]),
        "y_below_0p01_generated_count": int(a["n_er_y_below_0p01_in_sample"]),
        "y_below_0p01_post_s1_count": int(a["n_er_y_below_0p01_after_s1"]),
        "count_interpretation": "Generated, accepted, and retained counts are simulation multiplicities. The Table-I central physical mean and final full/display weights are fixed independently of those multiplicities.",
        "shape_interpretation": "Full unthinned unit-area observable-shape comparison; display samples are capped only for visualization and do not determine these MC estimates.",
    }


def main():
    args = parser().parse_args()
    mg_args = argparse.Namespace(**vars(args))
    output = Path(args.output)
    active_mode = args.er_efficiency_treatment
    b, y, ge, gv, curve_diagnostics = load_staged_er(
        args.er_staged_efficiency, args.er_1d_efficiency)
    g_input = load_efficiency(args.er_1d_efficiency)
    assert np.allclose(gv[ge <= 15], g_input.values[g_input.energy_kevnr <= 15])
    with h5py.File(args.background_models) as production:
        if str(production.attrs["er_efficiency_treatment"]) != active_mode:
            raise RuntimeError("Declared model and active ER decision disagree")
        if str(production.attrs["efficiency_treatment"]) != args.efficiency_treatment:
            raise RuntimeError("NR treatment must remain the chosen common mode")
        production_counts = {
            source: float(production[source]["hist_counts"][:].sum())
            for source in SOURCE_SPECS}

    columns, forward_values, forward_summary = forward_checks(
        args.er_staged_efficiency, args.er_1d_efficiency,
        args.er_2d_efficiency, args.nest_header,
        statistics=args.forward_check_statistics,
        seed=args.seed + 17_000_017)

    with tempfile.TemporaryDirectory(prefix="lz_er_efficiency_") as work:
        work = Path(work)
        generated = {}
        for label, mode in MODES.items():
            path = work / f"{label}.hdf5"
            generated[label] = generate(
                args, mg_args, mode, path, decision="er_efficiency_treatment")
        active_label = next(label for label, mode in MODES.items()
                            if mode == active_mode)
        identical, model_checks = same_model(args.background_models,
                                             generated[active_label])
        _, cross_mode_checks = same_model(
            generated["legacy_roi_only"], generated["staged_conditional"],
            # Accidentals retain their observable-space PDF but record this
            # intentional ER-decision label; NR marks ER not_applicable.
            ignore_attributes=("er_efficiency_treatment",))
        non_er_equal = all(check["datasets_equal"] and check["attributes_equal"]
                           for check in cross_mode_checks
                           if check["source"] in UNAFFECTED_SOURCES)
        if not identical or not non_er_equal:
            raise RuntimeError(f"ER model controls failed: active match={identical}, "
                               f"non-ER equality={non_er_equal}; {model_checks}; "
                               f"cross-mode={cross_mode_checks}")
        validations = {}
        for label, mode in MODES.items():
            metrics = validate(args, generated[label], args.efficiency_treatment,
                               work, mg_args, er_mode=mode)
            if metrics.get("er_efficiency_treatment") != mode:
                raise RuntimeError("Validation 4 ran with the wrong ER mode")
            validations[label] = metrics
        validations_pass = all(
            metrics["acceptance"]["pass"] for metrics in validations.values())
        reduced_budget = any(
            getattr(args, name) < make_parser().get_default(name)
            for name in ("n_er", "n_nr", "n_b8"))
        if not validations_pass and not (
                args.allow_reduced_budget_validation_failure and reduced_budget):
            raise RuntimeError(
                f"Both modes must pass the six existing gates: {validations}")
        if not validations_pass:
            print("WARNING: reduced-budget probe archives Validation-4 "
                  "finite-MC failures; this is not a production PASS.")

        metrics = {label: {} for label in MODES}
        with h5py.File(output, "w") as archive:
            archive.attrs["title"] = "ER staged conditional selection: complete matched twelve-source comparison"
            archive.attrs["schema_version"] = 1
            archive.attrs["active_er_efficiency_treatment"] = active_mode
            archive.attrs["comparison_modes"] = ",".join(MODES)
            archive.attrs["table1_counts_not_rates_note"] = "Physical central background counts are imposed from Table I. MC multiplicities and observable shapes are simulation outputs, not changes in physical rates."
            archive.attrs["forward_checks_not_factors"] = True
            archive.attrs["reduced_budget_diagnostic"] = reduced_budget
            archive.attrs["both_validation4_modes_pass"] = validations_pass
            for label, mode in MODES.items():
                with h5py.File(generated[label]) as model:
                    group = archive.create_group(label)
                    for key, value in model.attrs.items():
                        group.attrs[key] = value
                    for source in SOURCE_SPECS:
                        copy_source(model, source, group)
                        metrics[label][source] = source_metrics(group[source], mode)
                        if not np.isclose(metrics[label][source]["table1_central_physical_count"],
                                          production_counts[source], atol=1e-9):
                            raise RuntimeError(f"Physical Table-I mean changed for {source}")
                    for source in ER_SOURCES:
                        if not group[source].attrs["explicit_roi_applied"]:
                            raise RuntimeError(f"ER explicit ROI absent for {source}")
                        if mode == "staged_conditional":
                            if group[source].attrs["er_efficiency_stage"] != "trigger__nest_coincidence_s1__conditional_g_over_y":
                                raise RuntimeError("Unexpected staged selection label")
            shape_metrics = record_shapes(archive, load_band(args.response),
                                          sources=ER_SOURCES, modes=MODES)
            curves = archive.create_group("efficiency_inputs")
            curves.create_dataset("trigger_B", data=np.column_stack([b.energy_kevnr, b.values]))
            curves.create_dataset("cumulative_s1_Y", data=np.column_stack([y.energy_kevnr, y.values]))
            curves.create_dataset("G_overlap_then_G15", data=np.column_stack([ge, gv]))
            ratio = np.zeros_like(ge)
            positive = y(ge) > 1e-12
            ratio[positive] = gv[positive] / y(ge)[positive]
            curves.create_dataset("conditional_G_over_Y", data=np.column_stack([ge, ratio]))
            t4 = load_efficiency(args.er_2d_efficiency)
            curves.create_dataset("independent_total_2d_t4", data=np.column_stack([t4.energy_kevnr, t4.values]))
            for key in curves:
                curves[key].attrs["role"] = "conditional construction input" if key != "independent_total_2d_t4" else "diagnostic only; never apply in production"
            archive["efficiency_inputs"].attrs["diagnostics_json"] = dump_metadata(curve_diagnostics)
            forward_group = archive.create_group("forward_checks")
            dataset = forward_group.create_dataset("fixed_energy", data=forward_values)
            dataset.attrs["columns"] = columns
            forward_group.attrs["statistics_per_energy"] = args.forward_check_statistics
            forward_group.attrs["summary_json"] = dump_metadata(forward_summary)
            validation_group = archive.create_group("validation4")
            for label, values in validations.items():
                validation_group.create_dataset(label, data=dump_metadata(values))
            controls = {
                "active_model_exactly_reproduced": identical,
                "nr_and_observable_space_accidentals_identical_between_modes": non_er_equal,
                "all_twelve_source_physical_counts_unchanged": True,
                "full_histogram_matches_retained_full_mc": True,
                "full_and_display_samples_match_table1_means": True,
                "both_modes_report_and_pass_the_same_six_validation4_gates": validations_pass,
                "reduced_budget_diagnostic_not_production": reduced_budget,
                "no_improvement_in_fig4_agreement_required": True,
                "source_specific_seeding_changes_old_mc_realizations_not_physics": True,
                "matched_seed_is_not_common_random_numbers_after_rejection": True,
            }
            summary = {
                "active_er_efficiency_treatment": active_mode,
                "comparison_modes": MODES,
                "er_sources": list(ER_SOURCES),
                "unchanged_sources": list(UNAFFECTED_SOURCES),
                "source_metrics": metrics,
                "unit_shape_metrics": shape_metrics,
                "validation4": validations,
                "controls": controls,
                "input_diagnostics": curve_diagnostics,
                "forward_checks": forward_summary,
                "production_match_detail": model_checks,
                "non_er_match_detail": [check for check in cross_mode_checks
                                         if check["source"] in UNAFFECTED_SOURCES],
            }
            archive["summary"] = dump_metadata(summary)

    print("\nER controls PASS: active-model reproduction, unchanged NR/accidentals, "
          "Table-I normalization; both modes' six-gate Validation 4: "
          f"{'PASS' if validations_pass else 'FINITE-MC FAILURE (reduced-budget probe only)'}")
    for source in ER_SOURCES:
        for label in MODES:
            values = metrics[label][source]
            print(f"{source} [{label}]: {values['n_true_energy_samples']} -> "
                  f"trigger {values['n_after_trigger']} -> S1 {values['n_after_native_coincidence_and_s1']} "
                  f"-> G/Y {values['n_after_conditional_g_over_y']} -> ROI "
                  f"{values['full_roi_mc_count']}, display {values['retained_display_count']}; "
                  f"physical {values['table1_central_physical_count']}")
    for label, validation in validations.items():
        print(f"{label}: Validation 4 total RMS "
              f"{validation['residuals_dex']['total']['rms']:.9f} dex, tail "
              f"{validation['tail_abs_events']['total']:.9f} events, "
              f"PASS={validation['acceptance']['pass']}")
    print(json.dumps(forward_summary, indent=2))
    print(f"PASS: wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
