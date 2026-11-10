"""Probe: calling convention for nestpy GetYields/GetQuanta with the LZ
WS2024 parameter vectors passed explicitly (per the HEPData header's
instructions), plus sanity values for human inspection.

Nothing written to results/ — verdict is read by a human.
"""

import nestpy

d = nestpy.detectors.LZ_WS2024()
nc = nestpy.NESTcalc(d)

print("GetYields doc:")
print(nc.GetYields.__doc__)
print("GetQuanta doc:")
print(nc.GetQuanta.__doc__)

nr_params = list(d.ws_nr_parameters)
er_params = list(d.ws_er_parameters)
width_params = list(d.ws_nr_er_width_parameters)
print()
print("ws_nr_parameters :", nr_params)
print("ws_er_parameters :", er_params)
print("ws_width_params  :", width_params)

DENSITY = 2.90  # g/cm3, LXe at LZ operating conditions
FIELD = 96.5  # V/cm (LZ drift field)

NR = nestpy.interactions.NR  # NR
ER = nestpy.interactions.beta  # ER (electron recoil, generic beta)

print()
print("=== mean yields (GetQuanta mean path, fluctuations off via quench?) ===")
for label, itype, energy, A, Z in [
    ("NR  5 keVnr ", NR, 5.0, 131.293, 54.0),
    ("NR 10 keVnr ", NR, 10.0, 131.293, 54.0),
    ("NR 30 keVnr ", NR, 30.0, 131.293, 54.0),
    ("ER  2.7 keV (tritium endpoint region)", ER, 2.7, 0.00054858, 1.0),
    ("ER  5.2 keV (tritium mean-ish)      ", ER, 5.2, 0.00054858, 1.0),
    ("ER 10.0 keV ", ER, 10.0, 0.00054858, 1.0),
]:
    y = nc.GetYields(
        itype,
        energy,
        density=DENSITY,
        drift_field=FIELD,
        A=A,
        Z=Z,
        nr_parameters=nr_params,
        er_parameters=er_params,
    )
    nq = y.PhotonYield + y.ElectronYield
    print(
        f"  {label}: Nq={nq:8.2f}  "
        f"Nph={y.PhotonYield:7.2f}  Ne={y.ElectronYield:7.2f}  "
        f"Ly={y.PhotonYield / energy:6.2f}/keV  Qy={y.ElectronYield / energy:6.2f}/keV"
    )

print()
print("=== fluctuated quanta (5 draws, ER 5.2 keV) ===")
y = nc.GetYields(
    ER,
    5.2,
    density=DENSITY,
    drift_field=FIELD,
    A=0.00054858,
    Z=1.0,
    nr_parameters=nr_params,
    er_parameters=er_params,
)
for i in range(5):
    q = nc.GetQuanta(y, density=DENSITY, nr_er_width_parameters=width_params)
    print(
        f"  draw {i}: photons={q.photons:3d}  electrons={q.electrons:3d}"
    )
