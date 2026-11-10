"""Cross-check nestpy's built-in LZ_WS2022 detector against the HEPData header.

Probe (not a recipe): prints the detector's scalar parameters side by side
with the values in data/lz_ws2022_nest_header.hh (LZ_SR1.hh, HEPData record
10.17182/hepdata.144760.v1/r1). The verdict (match/mismatch) is read by a
human and recorded in the project notes; nothing is written to results/.

Note on the width vector: the header does not carry yield/width vectors as
scalars; its comment block instructs the user to call GetQuanta with the
edited vector {0.4, 0.4, 0.04, 0.5, 0.19, 2.25, 0.0015, 0.046452, 0.205,
0.45, -0.2} and SkewnessER = 0.0 for ER. The comparison band generator
(src/generate_ws2022_bands.py) passes exactly that. Also relevant: NEST
v2.4.x (what nestpy vendors) interprets NRERWidthsParam[6] > 0 as an
absolute Fano value, while v2.3.7 (what LZ SR1 ran) added [6]*sqrt(Nq)*
sqrt(E_field) unconditionally; v2.4.x's [6] < 0 branch reproduces the
v2.3.7 behavior, so the generator negates that one element (documented
there).
"""

import nestpy

HEADER_SCALARS = {
    # from data/lz_ws2022_nest_header.hh (LZ_SR1.hh)
    "g1": 0.113569,
    "sPEres": 0.31,
    "sPEthr": 0.10,
    "sPEeff": 1.0,
    "P_dphe": 0.214,
    "coinWind": 300,
    "coinLevel": 3,
    "numPMTs": 484,
    "g1_gas": 0.092103545,
    "s2Fano": 2.0,
    "s2_thr": 600.0 * 1.214,  # 728.4 phe = ROI uncorrected-S2 bound
    "E_gas": 8.42417,
    "eLife_us": 6500.0,
    "T_Kelvin": 174.1,
    "p_bar": 1.79,
    "dtCntr": 462.5,
    "dt_min": 86.0,
    "dt_max": 936.5,
    "radius": 688.0,
    "radmax": 728.0,
    "TopDrift": 1461.0,
    "anode": 1469.0,
    "gate": 1456.0,
    "cathode": 0.0,
    "drift_field_Vcm": 192.0,  # FitEF return value
}

d = nestpy.detectors.LZ_WS2022()
print("detector name:", d.get_name())
print("nest version vendored by nestpy:", nestpy.__nest_version__)
print()

print("=== scalar parameters: header vs nestpy built-in ===")
n_ok = n_diff = 0
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
    else:
        actual = getattr(d, key, None)
    ok = actual is not None and abs(float(actual) - expected) < 1e-6 * max(
        1, abs(expected)
    )
    n_ok += ok
    n_diff += not ok
    print(f"  {'OK ' if ok else 'DIFF'} {key:14s} header={expected!r:12} nestpy={actual!r}")

print(f"\n{n_ok} match, {n_diff} differ")
