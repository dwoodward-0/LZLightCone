"""Convert public HEPData ER efficiency tables into authored CSV inputs.

This authored-input utility verifies the complete 0..20-keV, 0.25-keV central
value tables. It does not invent errors or a covariance matrix. t3 (total 1D
ROI) serves as G(E) through 15 keV only because the user independently
verified the green/black overlap. t4 (total 2D ROI) is an independent future
forward-model check, not another selection factor.

The input YAML files are verbatim public HEPData table downloads. Outputs are
probe paths for inspection before authoring the CSVs into data/.
"""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
import yaml

TABLES = {
    "1d": {
        "table": "t3",
        "doi": "10.17182/hepdata.181350.v1/t3",
        "url": "https://www.hepdata.net/download/table/181350/WS2024%201D%20ROI%20efficiency/1/yaml",
        "name": "WS2024 1D ROI efficiency",
        "role": "G(E) through 15 keV, then hold at G(15); not applied as a second total efficiency",
    },
    "2d": {
        "table": "t4",
        "doi": "10.17182/hepdata.181350.v1/t4",
        "url": "https://www.hepdata.net/download/table/181350/WS2024%202D%20ROI%20efficiency/1/yaml",
        "name": "WS2024 2D ROI efficiency",
        "role": "independent forward-model consistency check; never a factor in the staged ER selection",
    },
}


def load_verified(path):
    table = yaml.safe_load(Path(path).read_text())
    independent = table["independent_variables"]
    dependent = table["dependent_variables"]
    if len(independent) != 1 or len(dependent) != 1:
        raise RuntimeError(f"{path}: expected one independent and one dependent variable")
    if independent[0]["header"]["name"] != "Energy [keV]":
        raise RuntimeError(f"{path}: unexpected energy label")
    if dependent[0]["header"]["name"] != "Efficiency":
        raise RuntimeError(f"{path}: unexpected efficiency label")
    energy = np.array([float(row["value"]) for row in independent[0]["values"]])
    efficiency = np.array([float(row["value"]) for row in dependent[0]["values"]])
    if not np.array_equal(energy, np.arange(81) * 0.25):
        raise RuntimeError(f"{path}: table is not the complete 0.25-keV grid")
    if not np.all(np.isfinite(efficiency)) or np.any(efficiency < 0) or np.any(efficiency > 1):
        raise RuntimeError(f"{path}: nonphysical central efficiency")
    serialized = Path(path).read_bytes()
    return energy, efficiency, hashlib.sha256(serialized).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kind", required=True, choices=TABLES)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    info = TABLES[args.kind]
    energy, efficiency, digest = load_verified(args.input)
    lines = [
        f"# Public HEPData {info['doi']}: {info['name']}, central values.",
        f"# Public record 181350, {info['table']}; verbatim YAML sha256={digest}",
        f"# Source: {info['url']}",
        f"# Analysis role: {info['role']}.",
        "# Energy units keVee, linear interpolation, held endpoints outside the tabulated support.",
        "# No numerical uncertainties or covariance matrix are included in the YAML dependent variable; none are invented here.",
        "# Green/black overlap through 15 keV is the user's visual observation, not an assertion that black and green are identical at all energies.",
        "energy_kevee,efficiency",
    ]
    with Path(args.output).open("w") as f:
        f.write("\n".join(lines) + "\n")
        writer = csv.writer(f, lineterminator="\n")
        for e, eff in zip(energy, efficiency):
            writer.writerow([f"{e:.2f}", f"{eff:.7f}"])
    print(f"{info['name']}: {len(energy)} central points, YAML sha256={digest}")
    if args.kind == "1d":
        print(f"G(15)={efficiency[energy == 15.0][0]:.7f}; black total-1D at 20 keV={efficiency[-1]:.7f}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
