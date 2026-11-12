"""Public staged ER selection curves and guarded G/Y construction.

B: blue Fig. 1 trigger, applied before NEST.
Y: yellow cumulative trigger + threefold coincidence + S1c > 3 phd.
G: SS/data-analysis cumulative efficiency. The user verified green/black
   overlap through 15 keV, so HEPData t3 central values supply G on that
   range; G is held at G(15) above 15 keV rather than following the black
   total-1D ROI decline.
C: G/Y, the average conditional SS/data-analysis acceptance after NEST's
   built-in threefold coincidence and an explicit S1c > 3 selection.

Only C is multiplied into the chain; B, Y, G are not multiplied together.
"""
from pathlib import Path

import numpy as np

from build_background_models import load_efficiency

# One Fig.-1 central-path stroke is 0.00493 efficiency. Allow two full
# strokes (0.00987) plus HEPData's 0.0005 central-value rounding. Ordering
# beyond this explicit extraction/rounding tolerance is a real input problem,
# not something to clamp silently. Valid raw C may exceed one only within
# this same tolerance, after which the small correction is recorded.
ORDERING_TOLERANCE = 0.01
Y_ZERO = 1e-12


def load_staged_er(staged_path, one_d_path, plateau_start=15.0):
    b = load_efficiency(staged_path, stage="trigger")
    y = load_efficiency(staged_path, stage="trigger_coincidence_s1")
    g_input = load_efficiency(one_d_path)
    if g_input.energy_kevnr[0] != 0 or g_input.energy_kevnr[-1] != 20:
        raise ValueError("ER G(E) must span 0..20 keV in HEPData t3")
    if g_input.energy_kevnr[np.argmin(np.abs(g_input.energy_kevnr - plateau_start))] != plateau_start:
        raise ValueError("G(15) must be explicitly tabulated")
    g_energy = np.arange(0., 20.0001, .25)
    g_values = np.interp(np.minimum(g_energy, plateau_start),
                         g_input.energy_kevnr, g_input.values)
    if not np.allclose(g_values[g_energy <= plateau_start],
                       g_input.values[g_input.energy_kevnr <= plateau_start]):
        raise ValueError("G on the overlap range differs from direct HEPData values")
    e = np.arange(0., 20.0001, .001)
    yv = y(e)
    gv = np.interp(e, g_energy, g_values)
    positive = yv > Y_ZERO
    if np.any(yv <= Y_ZERO):
        # At a zero cumulative S1-stage denominator, C is defined as zero
        # only if the already-published subsequent total is also consistent
        # with zero. No regularized 0/0 ratio is sampled.
        inconsistent = (yv <= Y_ZERO) & (gv > ORDERING_TOLERANCE)
        if inconsistent.any():
            raise ValueError("G(E) is nonzero where the cumulative S1 stage vanishes")
    max_ordering = float(np.max(gv - yv))
    if max_ordering > ORDERING_TOLERANCE:
        raise ValueError(f"G(E) exceeds Y(E) by {max_ordering}, above "
                         f"the {ORDERING_TOLERANCE} extraction tolerance")
    raw = np.zeros_like(e)
    raw[positive] = gv[positive] / yv[positive]
    if not np.all(np.isfinite(raw)) or np.any(raw < -ORDERING_TOLERANCE):
        raise ValueError("Nonfinite or negative conditional ER acceptance")
    max_raw = float(raw.max())
    if max_raw > 1 + ORDERING_TOLERANCE:
        raise ValueError(f"G/Y reaches {max_raw}: exceeds one by more than "
                         f"the {ORDERING_TOLERANCE} extraction tolerance")
    # The production inputs currently never reach this correction. Keeping
    # it explicit guards future input updates without masking a broad error.
    values = np.minimum(raw, 1.0)
    diagnostics = {
        "max_g_minus_y": max_ordering,
        "min_g_minus_y": float(np.min(gv - yv)),
        "max_raw_conditional": max_raw,
        "max_conditional_clipping": float(np.max(raw - values)),
        "ordering_tolerance": ORDERING_TOLERANCE,
        "denominator_zero_criterion": Y_ZERO,
        "plateau_start_kevee": plateau_start,
        "plateau_efficiency": float(g_input.values[g_input.energy_kevnr == plateau_start][0]),
        "low_denominator_diagnostic": "Raw C is evaluated only for Y>1e-12. Central-curve precision is limited near threshold; no arbitrary ratio-one substitution is used where Y is small but positive.",
        "uncertainty_note": "Central curves only. HEPData t3 publishes no numerical covariance. Fig.-1 stroke, interpolation and green/black overlap limitations are not treated as a cancellation of uncertainties.",
    }
    return b, y, g_energy, g_values, diagnostics
