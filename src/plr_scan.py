"""2D-PLR toy-MC median sensitivity over the full WIMP mass grid.

Runs plr_pilot.run_mass (identical construction to the 40 GeV/c2 pilot) for every
mass in the WIMP archive with one shared worker pool; mass i is seeded with
SeedSequence([seed, i]). Writes the sensitivity CSV to --output-csv (the declared
plr_sensitivity_curve when run as a recipe) and per-mass JSON detail to --outdir
(a scratch directory; a private temporary one if omitted). Without --output-csv the
CSV goes to --outdir/plr_sensitivity_scan.csv as before.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path
import sys
import tempfile
import time

import h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from plr_pilot import _init, add_common_args, check_scratch, run_mass
from stat_decisions import add_decision_args


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    add_common_args(ap)
    ap.add_argument("--outdir", type=Path, help="scratch dir for per-mass JSON (default: private temp dir)")
    ap.add_argument("--output-csv", type=Path, help="sensitivity CSV path, opened exactly as given")
    add_decision_args(ap, ["statistics", "background_components", "exposure_structure", "median_sensitivity_definition", "power_constraint", "test_statistic_sidedness"], "plr")
    ap.add_argument("--masses", nargs="*", type=float, help="subset of masses (GeV); default all")
    args = ap.parse_args()
    if args.outdir is None:
        args.outdir = Path(tempfile.mkdtemp(prefix="plr_scan_"))
    check_scratch(args.outdir)
    if not args.outdir.is_dir():
        raise ValueError("--outdir must exist")
    print(f"per-mass JSON -> {args.outdir}", flush=True)
    with h5py.File(args.signal, "r") as f:
        keys = sorted(f["ws2024"].keys(), key=lambda k: float(f[f"ws2024/{k}"].attrs["mass_gev"]))
        masses = {k: float(f[f"ws2024/{k}"].attrs["mass_gev"]) for k in keys}
    if args.masses:
        keys = [k for k in keys if masses[k] in set(args.masses)]
    t0 = time.perf_counter()
    rows = []
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.workers, initializer=_init, initargs=(args.signal, args.background, args.exposure_scale)) as pool:
        for i, k in enumerate(sorted(masses, key=masses.get)):
            if k not in keys:
                continue
            print(f"=== {masses[k]:g} GeV/c2 ===", flush=True)
            r = run_mass(pool, args, k, mass_index=i)
            (args.outdir / f"{k}.json").write_text(json.dumps(r, indent=1) + "\n")
            pc, raw = r["pc_limit_quantiles_cm2"], r["limit_quantiles_cm2"]
            rows.append([r["mass_gev"], pc["median"], pc["minus2sigma"], pc["minus1sigma"], pc["plus1sigma"],
                         pc["plus2sigma"], raw["median"], r["median_bootstrap_sigma_cm2"],
                         r["bkg_limits_beyond_grid"], r["bkg_limits_at_or_below_grid_start"],
                         r["fit_failures"]["signal_toys"] + r["fit_failures"]["bkg_toys"],
                         raw["minus2sigma"], raw["minus1sigma"], raw["plus1sigma"], raw["plus2sigma"],
                         r["exposure_scale"]])
            print(f"    elapsed {time.perf_counter() - t0:.0f} s", flush=True)
    header = ("mass_gev,median_pc_cm2,minus2sigma_pc_cm2,minus1sigma_pc_cm2,plus1sigma_pc_cm2,plus2sigma_pc_cm2,"
              "median_unconstrained_cm2,median_bootstrap_sigma_cm2,bkg_limits_beyond_grid,"
              "bkg_limits_at_grid_start,fit_failures,minus2sigma_cm2,minus1sigma_cm2,plus1sigma_cm2,"
              "plus2sigma_cm2,exposure_scale")
    csv_path = args.output_csv or args.outdir / "plr_sensitivity_scan.csv"
    csv_path.write_text(
        header + "\n" + "\n".join(",".join(f"{x:.10g}" for x in row) for row in rows) + "\n")
    print(f"done: {len(rows)} masses, {time.perf_counter() - t0:.0f} s")


if __name__ == "__main__":
    main()
