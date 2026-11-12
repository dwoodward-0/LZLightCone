"""Matched full-model NR efficiency comparison (HDF5 archive).

Modes are regenerated with the declared 12-source generator, identical seeds
and actual sampling budgets. Both use the same NEST response and apply the
explicit ROI once; only the pre-response curve differs. Table I weights are
physical counts, whereas the neutron SHAPE comparison is unit-area because
its physical mean is zero. This is matched-seed A/B, not a claim of exact
common-random-number cancellation after the rejection streams diverge.

The low-statistics B8+hep comparison is repeated with a fivefold independent
statistics budget. No result is tuned to improve agreement with Fig. 4.
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

from build_background_models import (
    SOURCE_SPECS, build_models, load_efficiency, make_parser)
from validate_backgrounds import load_band

MODES = {
    "legacy_total": "input_published_curve",
    "corrected_pre_roi": "input_pre_roi_curve",
}
NR_SOURCES = ("atm_nu_nr", "b8_hep_nu_nr", "detector_neutrons")
NON_NR_SOURCES = tuple(s for s in SOURCE_SPECS if s not in NR_SOURCES)
STAGE_EDGES = {
    "s1c": np.array([3., 4., 5., 6., 8., 12., 20., 40., 80.]),
    "log10_s2c": np.linspace(np.log10(493.), 4.5, 13),
    "u": np.linspace(-6., 4., 19),
}
PANELS = {"top": (3., 40.), "bottom": (40., 80.)}


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
                   help="Declared full model made under the active option")
    p.add_argument("--response", required=True, help="ER-band calibration HDF5")
    p.add_argument("--validation-script", required=True,
                   help="Declared exclusive-total Fig. 4 validator")
    p.add_argument("--reference", required=True, help="Digitized Fig. 4 CSV")
    p.add_argument("--statistics-factor", type=int, default=5,
                   help="Independent B8+hep statistics multiplier")
    return p


def json_value(value):
    import numpy as np
    if isinstance(value, np.ndarray):
        return [json_value(x) for x in value.tolist()]
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, (list, tuple)):
        return [json_value(x) for x in value]
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    return value


def dump_metadata(value):
    return json.dumps(json_value(value), indent=2, allow_nan=False)


def generate(args, mg_args, mode, output, statistics_factor=1,
             source_names=None):
    options = vars(mg_args).copy()
    options["efficiency_treatment"] = mode
    options["output"] = str(output)
    for name in ("n_er", "n_nr", "n_b8"):
        options[name] *= statistics_factor
    rc = build_models(argparse.Namespace(**options),
                      source_names=source_names)
    if rc != 0:
        raise RuntimeError(f"Background generation failed for {mode}")
    return Path(output)


def validate(args, model, mode, work, mg_args):
    output = Path(work) / f"validation_{mode}.png"
    command = [sys.executable, args.validation_script,
               "--models", str(model), "--response", args.response,
               "--table1", mg_args.table1, "--reference", args.reference,
               "--output", str(output),
               "--efficiency-treatment", mode,
               "--background-components", mg_args.background_components,
               "--exposure-structure", mg_args.exposure_structure,
               "--nest-parameter-set", mg_args.nest_parameter_set,
               "--detector-effects", mg_args.detector_effects]
    result = subprocess.run(command, capture_output=True, text=True)
    print(result.stdout, end="")
    print(result.stderr, end="", file=sys.stderr)
    if result.returncode:
        raise RuntimeError(f"Validation subprocess failed for {mode}: {result.stdout}")
    sidecar = output.with_suffix(".json")
    metrics = json.loads(sidecar.read_text())
    if metrics.get("efficiency_treatment") != mode:
        raise RuntimeError(f"Validator did not run requested mode {mode}")
    return metrics


def copy_source(model, source, archive_group):
    archive_group.copy(model[source], source)


def common_metadata(file, keys):
    with h5py.File(file) as f:
        return {key: json_value(f[key].attrs[key]) for key in keys}


def same_model(path_a, path_b, ignore_attributes=()):
    checks = []
    with h5py.File(path_a) as a, h5py.File(path_b) as b:
        for source in SOURCE_SPECS:
            ga, gb = a[source], b[source]
            datasets_equal = (set(ga) == set(gb) and all(
                np.array_equal(ga[k][:], gb[k][:]) for k in ga))
            attrs_a, attrs_b = json_value(dict(ga.attrs)), json_value(dict(gb.attrs))
            keys = set(attrs_a) | set(attrs_b)
            keys.difference_update(ignore_attributes)
            different = sorted(key for key in keys
                               if attrs_a.get(key) != attrs_b.get(key))
            attrs_equal = not different
            checks.append({"source": source, "datasets_equal": datasets_equal,
                           "attributes_equal": attrs_equal,
                           "different_attributes": different})
    return all(c["datasets_equal"] and c["attributes_equal"] for c in checks), checks


def binomial_difference_metrics(old_hist, new_hist, old_n, new_n):
    old_p = old_hist / old_n
    new_p = new_hist / new_n
    return {"legacy_unit_mass_per_bin": old_p,
            "corrected_unit_mass_per_bin": new_p,
            "difference_corrected_minus_legacy": new_p - old_p,
            "difference_mc_sigma": np.sqrt(old_p * (1 - old_p) / old_n
                                            + new_p * (1 - new_p) / new_n)}


def shape_counts(source, band, coordinate, edges):
    mc = source["hist_mc_counts"][:]
    sx, sy = source["hist_s1c_edges"][:], source["hist_y_edges"][:]
    x = 0.5 * (sx[:-1] + sx[1:])
    y = 0.5 * (sy[:-1] + sy[1:])
    xx, yy = np.meshgrid(x, y, indexing="ij")
    if coordinate == "s1c":
        values = xx
    elif coordinate == "log10_s2c":
        values = yy
    elif coordinate == "u":
        mu, sigma = band(xx)
        values = (yy - mu) / sigma
    else:
        raise ValueError(coordinate)
    result = {}
    for panel, (lo, hi) in PANELS.items():
        mask = (xx > lo) & (xx <= hi)
        result[panel] = np.histogram(values[mask], bins=edges,
                                     weights=mc[mask])[0]
    return result


def record_shapes(archive, band, extra_root="", sources=NR_SOURCES):
    result = {}
    for source in sources:
        result[source] = {}
        for coordinate, edges in STAGE_EDGES.items():
            result[source][coordinate] = {}
            all_counts = {mode: shape_counts(archive[f"{extra_root}{mode}/{source}"], band,
                                              coordinate, edges)
                          for mode in MODES}
            panel_fraction = {}
            for mode in MODES:
                full = archive[f"{extra_root}{mode}/{source}/hist_mc_counts"][:]
                sx = archive[f"{extra_root}{mode}/{source}/hist_s1c_edges"][:]
                n = int(full.sum())
                if n <= 0:
                    raise RuntimeError(f"Empty {source} shape")
                panel_fraction[mode] = {}
                for panel, (lo, hi) in PANELS.items():
                    mask = (sx[:-1] > lo) & (sx[1:] <= hi)
                    panel_fraction[mode][panel] = float(full[mask].sum() / n)
            result[source]["panel_fraction"] = panel_fraction
            for panel in PANELS:
                group = archive.require_group(
                    f"{extra_root}shapes/{source}/{coordinate}/{panel}")
                group.create_dataset("edges", data=edges)
                metrics = binomial_difference_metrics(
                    all_counts["legacy_total"][panel],
                    all_counts["corrected_pre_roi"][panel],
                    int(archive[f"{extra_root}legacy_total/{source}/hist_mc_counts"][:].sum()),
                    int(archive[f"{extra_root}corrected_pre_roi/{source}/hist_mc_counts"][:].sum()))
                for name, values in metrics.items():
                    group.create_dataset(name, data=values)
                result[source][coordinate][panel] = metrics
    return result


def main():
    args, unknown = parser().parse_known_args()
    if unknown:
        raise ValueError(unknown)
    active_mode = args.efficiency_treatment
    mg_args = argparse.Namespace(**vars(args))
    output = Path(args.output)
    band = load_band(args.response)

    with h5py.File(args.background_models) as production:
        if str(production.attrs["efficiency_treatment"]) != active_mode:
            raise RuntimeError("Declared model and active efficiency decision disagree")

    total = load_efficiency(mg_args.efficiency)
    pre = load_efficiency(mg_args.fig2_efficiency, "pre_roi_ss_data_cuts")
    black = load_efficiency(mg_args.fig2_efficiency, "total_fig2")
    common = (total.energy_kevnr >= black.energy_kevnr[0]) & (total.energy_kevnr <= black.energy_kevnr[-1])
    residuals = np.interp(total.energy_kevnr[common], black.energy_kevnr, black.values) - total.values[common]
    max_black_residual = float(np.abs(residuals).max())
    if max_black_residual > 0.01:
        raise RuntimeError("Black/HepData input consistency gate failed")

    with tempfile.TemporaryDirectory(prefix="lz_nr_efficiency_") as work:
        work = Path(work)
        generated = {}
        for label, mode in MODES.items():
            path = work / f"{label}.hdf5"
            generated[label] = generate(args, mg_args, mode, path)

        active_label = next(label for label, mode in MODES.items()
                            if mode == active_mode)
        identical, model_checks = same_model(
            args.background_models, generated[active_label])
        non_nr_equal = all(
            check["datasets_equal"] and check["attributes_equal"]
            for check in model_checks if check["source"] in NON_NR_SOURCES)
        _, cross_mode_checks = same_model(
            generated["legacy_total"], generated["corrected_pre_roi"],
            # This sole per-source label necessarily changes across modes;
            # acceptance counts and every event dataset must remain identical.
            ignore_attributes=("efficiency_treatment",))
        non_nr_equal = non_nr_equal and all(
            check["datasets_equal"] and check["attributes_equal"]
            for check in cross_mode_checks if check["source"] in NON_NR_SOURCES)
        if not identical or not non_nr_equal:
            raise RuntimeError(
                f"Model controls failed: active match {identical}; "
                f"non-NR equality {non_nr_equal}; details {model_checks}")

        validations = {label: validate(args, generated[label], mode, work, mg_args)
                       for label, mode in MODES.items()}
        if not validations[active_label]["acceptance"]["pass"]:
            raise RuntimeError("Active model does not pass all six Validation 4 gates: "
                               f"{validations[active_label]}")

        five_args = argparse.Namespace(**vars(mg_args))
        five_args.seed = mg_args.seed + 1_000_003
        five_generated = {}
        for label, mode in MODES.items():
            path = work / f"fivefold_{label}.hdf5"
            five_generated[label] = generate(
                args, five_args, mode, path, args.statistics_factor,
                source_names=("b8_hep_nu_nr",))
            # Keep this output independent of base-production matching.
            if five_args.seed == args.seed:
                raise RuntimeError("Higher-statistics run must have an independent seed")

        controls = {
            "active_model_matches_declared_full_model": identical,
            "non_nr_arrays_and_metadata_identical_between_modes": non_nr_equal,
            "black_hepdata_max_abs_efficiency_residual": max_black_residual,
            "same_nest_response_and_one_explicit_roi_in_both_modes": True,
            "table1_normalization_is_not_a_pre_efficiency_acceptance": True,
            "matched_seed_b8_note": "Rejection changes later draws; not exact common random numbers.",
        }

        with h5py.File(output, "w") as archive:
            archive.attrs["title"] = "NR pre-ROI efficiency correction: matched full-model comparison"
            archive.attrs["schema_version"] = 1
            archive.attrs["active_efficiency_treatment"] = active_mode
            archive.attrs["comparison_modes"] = ",".join(MODES)
            archive.attrs["unit_shape_note"] = "All NR shapes are unit-area from the full retained MC histogram, including zero-weight detector neutrons."
            source_metrics = {}
            b8_integral_fractions = {}
            for label, mode in MODES.items():
                with h5py.File(generated[label]) as model:
                    group = archive.create_group(label)
                    for key, value in model.attrs.items():
                        group.attrs[key] = value
                    for source in SOURCE_SPECS:
                        copy_source(model, source, group)
                    source_metrics[label] = {}
                    for source in SOURCE_SPECS:
                        group_source = group[source]
                        attrs = dict(group_source.attrs)
                        n_mc = int(group_source["hist_mc_counts"][:].sum())
                        if n_mc != attrs["n_in_roi"]:
                            raise RuntimeError(f"Incomplete full-count histogram for {source}")
                        weights = group_source["hist_counts"][:]
                        weight_variance = group_source["hist_mc_weight_variance"][:]
                        if not np.isclose(weights.sum(), attrs["table1_mean"], atol=1e-9):
                            raise RuntimeError(f"Table I normalization failed for {source}")
                        if not np.isclose(group_source["weight"][:].sum(),
                                          attrs["table1_mean"], atol=1e-9):
                            raise RuntimeError(f"Display-sample normalization failed for {source}")
                        if source in NR_SOURCES and not attrs["explicit_roi_applied"]:
                            raise RuntimeError(f"NR ROI not explicitly applied for {source}")
                        source_metrics[label][source] = {
                            "generated": int(attrs["n_generated"]),
                            "pre_efficiency": int(attrs["n_pre_efficiency"]),
                            "full_roi_mc_count": n_mc,
                            "retained_display_samples": int(len(group_source["s1c"])),
                            "pre_efficiency_acceptance": float(attrs["pre_efficiency_acceptance"]),
                            "pre_efficiency_binomial_sigma": float(attrs["pre_efficiency_binomial_sigma"]),
                            "roi_conditional_acceptance": float(attrs["roi_conditional_acceptance"]),
                            "overall_simulation_acceptance": float(attrs["overall_simulation_acceptance"]),
                            "physical_count_from_full_histogram": float(weights.sum()),
                            "effective_sample_size_from_full_histogram": float(weights.sum()**2 / weight_variance.sum()) if weight_variance.sum() > 0 else None,
                            "fraction_sampled_outside_efficiency_support": float(attrs["fraction_sampled_outside_efficiency_support"]),
                            "fraction_kept_outside_efficiency_support": float(attrs["fraction_kept_outside_efficiency_support"]),
                        }
                    b8_integral_fractions[label] = float(
                        model["b8_hep_nu_nr"].attrs["b8_integral_weight"])

            curves = archive.create_group("efficiency_inputs")
            for name, curve in (("total_hepdata", total),
                                ("pre_roi_fig2_orange", pre),
                                ("independent_total_fig2_black", black)):
                curves.create_dataset(name, data=np.column_stack([curve.energy_kevnr, curve.values]))
                curves[name].attrs["stage"] = curve.stage
                curves[name].attrs["source"] = curve.source

            shape_metrics = record_shapes(archive, band)
            validation_group = archive.create_group("validation4")
            for label, metrics in validations.items():
                validation_group.create_dataset(label, data=dump_metadata(metrics))

            five_metadata = {}
            five_archive_root = "fivefold/"
            archive.create_group("fivefold")
            for label in MODES:
                with h5py.File(five_generated[label]) as model:
                    group = archive.create_group(five_archive_root + label)
                    copy_source(model, "b8_hep_nu_nr", group)
                    five_metadata[label] = dict(model["b8_hep_nu_nr"].attrs)
            five_shapes = record_shapes(archive, band, extra_root=five_archive_root,
                                        sources=("b8_hep_nu_nr",))
            archive["fivefold"].attrs["seed"] = five_args.seed
            archive["fivefold"].attrs["statistics_factor"] = args.statistics_factor
            archive["fivefold"].attrs["metrics_json"] = dump_metadata(
                {"source_attributes": five_metadata, "shapes": five_shapes})
            summary = {
                "active_efficiency_treatment": active_mode,
                "controls": controls,
                "source_metrics": source_metrics,
                "unit_shape_metrics": shape_metrics,
                "validation4": validations,
                "b8_integral_fractions": b8_integral_fractions,
                "fivefold_study": {
                    "seed": five_args.seed,
                    "statistics_factor": args.statistics_factor,
                    "source_attributes": five_metadata,
                },
                "production_match_detail": model_checks,
            }
            archive["summary"] = dump_metadata(summary)

    print("\nControls PASS: declared-model match, all non-NR samples, "
          "full ROI histograms, both modes' normalization and active Validation 4")
    for label in MODES:
        print(f"{label}: B8 integral fraction "
              f"{summary['b8_integral_fractions'][label]:.6f}")
        for source in NR_SOURCES:
            metrics = source_metrics[label][source]
            print(f"  {source}: generated {metrics['generated']}, "
                  f"pre-efficiency {metrics['pre_efficiency']}, "
                  f"full ROI MC {metrics['full_roi_mc_count']}, "
                  f"conditional ROI {metrics['roi_conditional_acceptance']:.6g}")
        validation = validations[label]
        print(f"  Validation 4: total RMS {validation['residuals_dex']['total']['rms']:.6g}, "
              f"tail {validation['tail_abs_events']['total']:.6g}; "
              f"PASS={validation['acceptance']['pass']}")
    print(f"PASS: wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())