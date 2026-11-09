"""Cross-check nestpy's built-in LZ_WS2024 detector against the HEPData header.

Probe (not a recipe): prints the detector's scalar parameters and the
WS2024 yield/width parameter vectors, side by side with the values in
data/lz_ws2024_nest_header.hh. The verdict (match/mismatch) is read by a human
and recorded in the project notes; nothing is written to results/.
"""

import nestpy

HEADER_SCALARS = {
    # from data/lz_ws2024_nest_header.hh (LZ_Result2024.hh)
    "g1": 0.1122,
    "sPEres": 0.338,
    "sPEthr": 0.10,
    "sPEeff": 1.0,
    "P_dphe": 0.214,
    "coinWind": 150,
    "coinLevel": 3,
    "numPMTs": 481,
    "g1_gas": 0.076404,
    "s2Fano": 4.0,
    "s2_thr": 0.0,
    "E_gas": 8.301178,
    "eLife_us": 9000.0,
    "T_Kelvin": 175.9,
    "p_bar": 1.86,
    "dtCntr": 525.0,
    "dt_min": 71.0,
    "dt_max": 1029.0,
    "radius": 688.0,
    "radmax": 728.0,
    "TopDrift": 1461.0,
    "anode": 1469.0,
    "gate": 1456.0,
    "cathode": 0.0,
    "drift_field_Vcm": 96.5,  # FitEF return value
}

HEADER_VECTORS = {
    # yield-model tuning vectors required by the header (lines 20-22)
    "NRYieldsParam": [10.19, 1.11, 0.0498, -0.0533, 12.46, 0.2942, 1.899,
                      0.3197, 2.066, 0.509, 0.996, 0.999],
    "ERYieldsParam": [12.4886, 85.0, 0.6050, 2.14687, 25.721, -1.0, 59.651,
                      3.6869, 0.2872, 0.1121],
    "NRERWidthsParam": [0.404, 0.393, 0.0383, 0.497, 0.1906, 2.220, 0.3,
                        0.04311, 0.15505, 0.46894, -0.26564],
}

d = nestpy.detectors.LZ_WS2024()
print("detector name:", d.get_name())
print()

print("=== scalar parameters: header vs nestpy built-in ===")
for key, expected in HEADER_SCALARS.items():
    if key == "drift_field_Vcm":
        actual = d.FitEF(0, 0, 0)
    elif key in ("sPEres", "sPEeff"):
        actual = d.get_sPEres() if key == "sPEres" else d.get_sPEeff()
    elif key == "sPEthr":
        actual = d.get_sPEthr()
    elif key in ("coinWind", "coinLevel"):
        actual = d.coinWind if key == "coinWind" else d.coinLevel
    elif key == "numPMTs":
        actual = d.get_numPMTs()
    elif key == "s2Fano":
        actual = d.get_s2Fano()
    elif key == "s2_thr":
        actual = d.get_s2_thr()
    elif key == "eLife_us":
        actual = d.get_eLife_us()
    elif key == "T_Kelvin":
        actual = d.get_T_Kelvin()
    elif key == "p_bar":
        actual = d.get_p_bar()
    elif key == "numPMTs":
        actual = d.get_numPMTs()
    else:
        actual = getattr(d, key, None)
    ok = "OK " if actual is not None and abs(float(actual) - expected) < 1e-6 * max(1, abs(expected)) else "DIFF"
    print(f"  {ok} {key:14s} header={expected!r:12} nestpy={actual!r}")

print()
print("=== yield-model parameter vectors ===")
for attr in ["nr_parameters", "er_parameters", "nr_er_width_parameters",
             "ws_nr_parameters", "ws_er_parameters", "ws_nr_er_width_parameters",
             "henr_nr_parameters", "henr_er_parameters",
             "henr_nr_er_width_parameters"]:
    vec = getattr(d, attr, None)
    if vec is None:
        print(f"  {attr}: (absent)")
        continue
    try:
        vals = list(vec)
    except TypeError:
        print(f"  {attr}: not iterable -> {vec!r}")
        continue
    print(f"  {attr} = {[round(float(v), 6) for v in vals]}")

print()
print("=== header vectors for direct comparison ===")
for name, vec in HEADER_VECTORS.items():
    print(f"  {name} = {vec}")
