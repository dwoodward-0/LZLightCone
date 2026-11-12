"""Independent fixed-energy ER forward checks for the staged selection.

This module does not supply another selection factor. It samples NEST at
fixed true energies and compares (a) B(E) times native coincidence/S1
acceptance with the cumulative yellow curve Y(E), and (b) the staged total
acceptance with the independent HEPData total-2D curve. These are central-
curve diagnostics, not assumptions that covariance cancels in G/Y.
"""
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_background_models import ROI_S1C_MIN, ResponseEngine, load_efficiency, roi_mask
from er_efficiency import ORDERING_TOLERANCE, load_staged_er


def forward_checks(staged_path, one_d_path, two_d_path, header_path,
                   statistics=5000, seed=17000061):
    if statistics < 1:
        raise ValueError("Forward-check MC statistics must be positive")
    b, y, ge, gv, diagnostics = load_staged_er(staged_path, one_d_path)
    total2d = load_efficiency(two_d_path)
    energies = total2d.energy_kevnr.copy()
    trigger = b(energies)
    yellow = y(energies)
    green = np.interp(energies, ge, gv)
    conditional = np.zeros_like(energies)
    positive = yellow > 1e-12
    conditional[positive] = np.minimum(green[positive] / yellow[positive], 1.)
    s1_acceptance = np.zeros_like(energies)
    roi_acceptance = np.zeros_like(energies)
    s1_sigma = np.zeros_like(energies)
    roi_sigma = np.zeros_like(energies)
    negative_flags = np.zeros_like(energies, dtype=int)
    engine = ResponseEngine(header_path, int(seed))
    for index, energy in enumerate(energies):
        if energy == 0:
            continue
        s1c, s2c = engine.generate(np.full(statistics, energy), nr=False)
        negative_flags[index] = np.count_nonzero(s1c < 0)
        passes = roi_mask(s1c, s2c)
        assert np.all(s1c[passes] > ROI_S1C_MIN)
        s1_acceptance[index] = np.mean(s1c > ROI_S1C_MIN)
        roi_acceptance[index] = np.mean(passes)
        s1_sigma[index] = np.sqrt(s1_acceptance[index] * (1 - s1_acceptance[index]) / statistics)
        roi_sigma[index] = np.sqrt(roi_acceptance[index] * (1 - roi_acceptance[index]) / statistics)
    predicted_y = trigger * s1_acceptance
    predicted_total2d = trigger * conditional * roi_acceptance
    # The only uncertainty term here is binomial MC. The extraction/rounding
    # allowance is a precision allowance, NOT a statistically fitted error;
    # no unpublished source uncertainty or G/Y covariance is fabricated.
    central_precision = 2 * ORDERING_TOLERANCE
    diagnostic_tolerance_y = 2 * trigger * s1_sigma + central_precision
    diagnostic_tolerance_total = 2 * trigger * conditional * roi_sigma + central_precision
    values = np.column_stack([
        energies, trigger, yellow, green, conditional,
        s1_acceptance, s1_sigma, roi_acceptance, roi_sigma,
        predicted_y, predicted_total2d, total2d.values,
        predicted_y - yellow, predicted_total2d - total2d.values,
        diagnostic_tolerance_y, diagnostic_tolerance_total,
    ])
    columns = (
        "energy_kevee", "trigger_B", "published_cumulative_Y",
        "green_overlap_G_then_G15", "conditional_G_over_Y",
        "nest_s1_acceptance", "nest_s1_binomial_sigma",
        "nest_full_roi_acceptance", "nest_full_roi_binomial_sigma",
        "predicted_cumulative_Y", "predicted_total_2d", "published_total_2d",
        "s1_difference_model_minus_published", "total_2d_difference_model_minus_published",
        "s1_finite_mc_plus_precision_allowance", "total_2d_finite_mc_plus_precision_allowance",
    )
    assert values.shape[1] == len(columns)
    low_y = (yellow > 0) & (yellow < .01)
    summary = {
        "statistics_per_energy": int(statistics),
        "seed": int(seed),
        "energy_min_kevee": float(energies[energies > 0].min()),
        "energy_max_kevee": float(energies.max()),
        "input_diagnostics": diagnostics,
        "n_negative_s1_coincidence_flags": int(negative_flags.sum()),
        "n_low_y_grid_points": int(np.count_nonzero(low_y)),
        "low_y_interpretation": "Near threshold, interpolation at the scale of the digitized central-path stroke can dominate G/Y. Low-Y points are reported, not regularized to one or used as acceptance gates.",
        "central_precision_allowance": central_precision,
        "allowance_note": "Twice the 0.01 extraction/rounding tolerance, added to two finite-MC binomial sigmas. This is a declared diagnostic precision allowance, not a recovered HEPData uncertainty or covariance. Disagreement is reported rather than hidden or tuned.",
        "s1_max_abs_difference_for_y_at_least_0p01": float(np.max(np.abs(predicted_y[~low_y] - yellow[~low_y]))),
        "total_2d_max_abs_difference": float(np.max(np.abs(predicted_total2d - total2d.values))),
        "total_2d_rms_difference": float(np.sqrt(np.mean((predicted_total2d - total2d.values) ** 2))),
        "total_2d_difference_at_15_kevee": float(predicted_total2d[energies == 15][0] - total2d.values[energies == 15][0]),
        "total_2d_difference_at_20_kevee": float(predicted_total2d[energies == 20][0] - total2d.values[energies == 20][0]),
        "n_s1_within_diagnostic_tolerance_y_at_least_0p01": int(np.count_nonzero(np.abs(predicted_y[~low_y] - yellow[~low_y]) <= diagnostic_tolerance_y[~low_y])),
        "n_s1_compared_y_at_least_0p01": int(np.count_nonzero(~low_y)),
        "n_total_2d_within_diagnostic_tolerance": int(np.count_nonzero(np.abs(predicted_total2d - total2d.values) <= diagnostic_tolerance_total)),
        "n_total_2d_compared": int(len(energies)),
        "production_note": "Nothing in this module is multiplied into production. t4 is never another efficiency factor. The G(15) plateau intentionally does not follow t4's high-energy decline; the observed residual is a diagnostic of the energy-conditional approximation and unpublished observable-dependent selections, not an automatic model failure.",
    }
    return columns, values, summary
