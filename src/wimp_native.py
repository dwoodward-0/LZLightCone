"""Checked Python access to NEST's native SI WIMP spectrum and rate.

This module does not implement a competing scattering model. All differential
rates, isotope sampling, Helm form factors, halo parameters and recoil samples
come from the pinned nestpy source's TestSpectra implementation. The local
binding-only extension makes WIMP_dRate and the native preparation metadata
available to Python.

Native WIMP_dRate returns differential events / kg / day / keVnr at the
built-in SI per-nucleon reference cross section 1e-36 cm2. The preparation
integral is its native sum over [0, 100) keVnr in 1e-4 keV bins. Every such
differential-rate evaluation randomly selects a xenon isotope; independently
seeded native integrals therefore differ separately from response-event MC
noise. NEST's global RNG is used serially; callers reseed the response after
sampling when they need an independent response stream.
"""

from __future__ import annotations

from importlib.metadata import version
import math

import nestpy
import numpy as np

REFERENCE_SIGMA_CM2 = 1.0e-36
DAYS_PER_YEAR = 365.25
NEST_PERIOD_DAYS = 365.24
V_SUN_KM_S = 250.2
EARTH_MODULATION_AMPLITUDE_KM_S = 0.49 * 29.8
# cos(2*pi*(day/365.24 - 0.415)) = 0: the scalar Solar-speed setting.
REPRESENTATIVE_DAY = NEST_PERIOD_DAYS * (0.415 - 0.25)
NATIVE_INTEGRAL_EMAX_KEV = 100.0
TAIL_EMAX_KEV = 500.0
# Below this folded fraction the omitted tail is below a 0.01% yield effect;
# every larger tail is sampled, not absorbed into a normalization correction.
TAIL_FRACTION_NUMERICAL_FLOOR = 1.0e-4


def earth_speed_km_s(day: float) -> float:
    """The exact scalar Earth-speed convention in the locked NEST source."""
    phase = 2.0 * math.pi * (day / NEST_PERIOD_DAYS - 0.415)
    return V_SUN_KM_S + EARTH_MODULATION_AMPLITUDE_KM_S * math.cos(phase)


def phase_midpoints(n_phases: int = 12) -> np.ndarray:
    """Equally weighted phase midpoints over NEST's one-year period."""
    if n_phases < 1:
        raise ValueError("n_phases must be positive")
    return (np.arange(n_phases, dtype=float) + 0.5) * NEST_PERIOD_DAYS / n_phases


def native_provenance() -> dict[str, str]:
    return {
        "package": "nestpy",
        "version": version("nestpy"),
        "upstream_sdist_sha256":
            "258155894377e24c02a763e8e58319e01fbd42f12aa7a8c960084848b1a84e84",
        "patch": "vendor/NESTPY_WIMP_CONNECTOR.patch",
        "native_rate_units": "events/kg/day/keVnr",
        "native_integral_units": "events/kg/day",
        "reference_sigma_cm2": str(REFERENCE_SIGMA_CM2),
    }


def sample_native(
    mass_gev: float,
    n_events: int,
    seed: int,
    day: float = REPRESENTATIVE_DAY,
    e_step: float = 5.0,
) -> dict:
    """One NEST preparation, its samples and its matching normalization.

    These samples and the rate must stay together: a separately seeded
    preparation can have a different native isotope-sampling realization.
    Zero-length sampling probes support and normalization; production callers
    separately require nonzero samples and a positive retained-event count.
    """
    if not np.isfinite(mass_gev) or mass_gev <= 0:
        raise ValueError("mass_gev must be finite and positive")
    if n_events < 0:
        raise ValueError("n_events must be non-negative")
    if not np.isfinite(day) or not np.isfinite(e_step) or e_step <= 0:
        raise ValueError("day must be finite and e_step finite and positive")
    nestpy.RandomGen.rndm().set_seed(int(seed))
    out = nestpy.spectra.WIMP_with_metadata(
        n_events, mass=float(mass_gev), eStep=float(e_step), day=float(day))
    energies = np.asarray(out["samples"], dtype=float)
    metadata = {
        key: float(out[key]) for key in
        ("integral", "xMax", "divisor", "mass", "eStep", "day")
    }
    metadata.update({
        "seed": int(seed),
        "n_samples": int(n_events),
        "n_zero_samples": int(np.count_nonzero(energies == 0.0)),
        "n_negative_samples": int(np.count_nonzero(energies < 0.0)),
        "n_nonfinite_samples": int(np.count_nonzero(~np.isfinite(energies))),
        "n_above_native_integral_support": int(
            np.count_nonzero(energies >= NATIVE_INTEGRAL_EMAX_KEV)),
        "n_above_xmax": int(np.count_nonzero(energies > metadata["xMax"])),
        "energy_min_kevnr": float(np.min(energies)) if len(energies) else None,
        "energy_max_kevnr": float(np.max(energies)) if len(energies) else None,
    })
    if metadata["n_nonfinite_samples"] or metadata["n_negative_samples"]:
        raise ValueError("Native WIMP sampler returned invalid energies")
    if not math.isfinite(metadata["integral"]) or metadata["integral"] <= 0:
        raise ValueError("Native WIMP integral must be finite and positive")
    if not 0 < metadata["xMax"] <= NATIVE_INTEGRAL_EMAX_KEV:
        raise ValueError(f"Unexpected native xMax={metadata['xMax']}")
    return {"energies": energies, "metadata": metadata}


def integrate_native(
    mass_gev: float,
    day: float,
    energy_kevnr: np.ndarray,
    seed: int,
    efficiency=None,
) -> dict:
    """Quadrature of native WIMP_dRate on a supplied deterministic grid.

    This checks the native preparation integral on common support; it is not a
    replacement scattering calculation. The array of native rate evaluations
    is retained so callers can bin it on the same energy grid as the sampler.
    Reported standard errors quantify only the independent isotope draws at
    the grid points, not quadrature or downstream response-MC uncertainty.
    """
    e = np.asarray(energy_kevnr, dtype=float)
    if e.ndim != 1 or len(e) < 2 or np.any(np.diff(e) <= 0):
        raise ValueError("energy_kevnr must be a strictly increasing 1D grid")
    if not np.all(np.isfinite(e)) or e[0] < 0:
        raise ValueError("quadrature energies must be finite and non-negative")
    nestpy.RandomGen.rndm().set_seed(int(seed))
    density = np.array([
        nestpy.spectra.WIMP_dRate(float(x), float(mass_gev), float(day))
        for x in e
    ], dtype=float)
    if np.any(~np.isfinite(density)) or np.any(density < 0):
        raise ValueError("Native quadrature returned invalid rate")
    widths = np.diff(e)
    trap = 0.5 * (density[:-1] + density[1:])
    se = float(np.sqrt(np.sum((widths * 0.5 *
                              (density[1:] - density[:-1])) ** 2) / len(widths)))
    result = {
        "energy_kevnr": e,
        "density": density,
        "integral": float(np.sum(widths * trap)),
        "isotope_grid_standard_error": se,
        "seed": int(seed),
        "mass_gev": float(mass_gev),
        "day": float(day),
    }
    if efficiency is not None:
        epsilon = np.asarray(efficiency(e), dtype=float)
        if epsilon.shape != e.shape or np.any(~np.isfinite(epsilon)):
            raise ValueError("Invalid efficiency values on native grid")
        result["efficiency_folded_integral"] = float(
            np.sum(widths * 0.5 * (density[:-1] * epsilon[:-1] +
                                   density[1:] * epsilon[1:])))
    return result


def native_support_diagnostics(mass_gev, day, seed, efficiency,
                               prepared_integral, n_grid_events=30000,
                               e_step=5.0, e_steps=(5.0, 7.5, 10.0)):
    """Common-support quadrature, sampler controls, and native >100keV tail.

    The native dRate is the only differential spectrum used. The tail is an
    exactly sampled piecewise-linear histogram of that rate, not an exponential
    approximation or a fit to LZ. The upper [500, 5000] keV region is recorded
    as an unadded residual guard; at high masses its exponential tail is tiny.
    Quadrature/isotope uncertainty remains distinct from downstream event MC.
    """
    if e_step != 5.0:
        raise ValueError("Production grid is pinned to native eStep=5; controls may vary")
    fine = np.linspace(0.0, NATIVE_INTEGRAL_EMAX_KEV, 20001)
    tail_grid = np.linspace(NATIVE_INTEGRAL_EMAX_KEV, TAIL_EMAX_KEV, 20001)
    base_quad = integrate_native(mass_gev, day, fine, seed + 4100, efficiency)
    tail_quad = integrate_native(mass_gev, day, tail_grid, seed + 4200, efficiency)
    upper_grid = np.geomspace(TAIL_EMAX_KEV, 5000.0, 2001)
    upper_quad = integrate_native(mass_gev, day, upper_grid, seed + 4300, efficiency)
    tail_fraction = tail_quad["efficiency_folded_integral"] / (
        base_quad["efficiency_folded_integral"])
    if upper_quad["efficiency_folded_integral"] / (
            base_quad["efficiency_folded_integral"]) > TAIL_FRACTION_NUMERICAL_FLOOR:
        raise RuntimeError("Native tail beyond 500 keVnr is not negligible; do not truncate")
    use_tail = (tail_quad["efficiency_folded_integral"] > 0 and
                tail_fraction > TAIL_FRACTION_NUMERICAL_FLOOR)
    grid_records = []
    for i, step in enumerate(e_steps):
        # Invalid native grids are failures, not an opportunity to select a grid
        # by agreement with a reference figure. The pinned production default
        # must work; the independent controls document numerical stability.
        try:
            item = sample_native(mass_gev, n_grid_events,
                                 seed + 4500 + i, day=day, e_step=step)
        except (RuntimeError, ValueError) as exc:
            grid_records.append({"eStep": step, "numerical_error": str(exc)})
            continue
        e = item["energies"]
        eps = efficiency(np.maximum(e, 0.0))
        eps[e <= 0] = 0.0
        rng = np.random.default_rng(seed + 4600 + i)
        n_selected = int(np.sum(rng.uniform(size=len(e)) < eps))
        grid_records.append({
            "metadata": item["metadata"],
            "efficiency_acceptance": n_selected / len(e),
            "efficiency_acceptance_binomial_sigma": float(np.sqrt(
                (n_selected / len(e)) * (1 - n_selected / len(e)) / len(e))),
            "prepared_integral_relative_to_production":
                item["metadata"]["integral"] / prepared_integral - 1.0,
        })
    controls = []
    control_errors = []
    for i in range(3):
        try:
            controls.append(sample_native(
                mass_gev, 0, seed + 4700 + i, day=day)["metadata"])
        except (RuntimeError, ValueError) as exc:
            # A native control failure is evidence about this numerical path,
            # not permission to choose a preparation by agreement with LZ.
            control_errors.append({
                "seed": int(seed + 4700 + i), "e_step": float(e_step),
                "numerical_error": str(exc)})
    if len(controls) < 2:
        raise RuntimeError(
            f"Fewer than two successful independent native control preparations: "
            f"{control_errors}")
    control_spread = float(np.std([x["integral"] for x in controls], ddof=1))
    return {
        "native_integral_0_100": base_quad["integral"],
        "native_integral_0_100_isotope_grid_standard_error":
            base_quad["isotope_grid_standard_error"],
        "efficiency_folded_integral_0_100":
            base_quad["efficiency_folded_integral"],
        "native_integral_100_500": tail_quad["integral"],
        "native_integral_100_500_isotope_grid_standard_error":
            tail_quad["isotope_grid_standard_error"],
        "efficiency_folded_integral_100_500":
            tail_quad["efficiency_folded_integral"],
        "efficiency_folded_tail_fraction_relative_to_0_100": tail_fraction,
        "native_integral_500_5000_guard_not_added": upper_quad["integral"],
        "efficiency_folded_integral_500_5000_guard_not_added":
            upper_quad["efficiency_folded_integral"],
        "guard_fraction_of_prepared_integral":
            upper_quad["integral"] / prepared_integral,
        "tail_numerical_floor": TAIL_FRACTION_NUMERICAL_FLOOR,
        "tail_extension_used": bool(use_tail),
        "tail_native_rate_boundary_kevnr": [NATIVE_INTEGRAL_EMAX_KEV,
                                            TAIL_EMAX_KEV],
        "native_quadrature_is_independent_of_production_preparation": True,
        "quadrature_step_kevnr": float(fine[1] - fine[0]),
        "prepared_integral_minus_independent_0_100_relative":
            prepared_integral / base_quad["integral"] - 1.0,
        "grid_checks_numerical_only": grid_records,
        "independent_preparations": controls,
        "independent_preparation_numerical_failures": control_errors,
        "independent_preparation_n_successful": len(controls),
        "independent_preparation_integral_standard_deviation": control_spread,
        "tail_integral_is_added_not_rescaled": bool(use_tail),
        "tail_energy_kevnr": tail_grid if use_tail else np.empty(0),
        "tail_differential_rate": tail_quad["density"] if use_tail else np.empty(0),
    }


def sample_native_tail(diagnostics, n_events, seed):
    """Sample exactly the stored linearly-interpolated native-tail density."""
    e = diagnostics["tail_energy_kevnr"]
    density = diagnostics["tail_differential_rate"]
    if len(e) == 0 or not diagnostics["tail_extension_used"]:
        return np.empty(0)
    cell_rate = 0.5 * (density[:-1] + density[1:])
    cell_area = cell_rate * np.diff(e)
    if np.any(~np.isfinite(cell_area)) or np.any(cell_area < 0) \
            or not cell_area.sum() > 0:
        raise ValueError("Native tail has invalid sampling density")
    rng = np.random.default_rng(seed)
    index = rng.choice(len(cell_area), size=n_events, p=cell_area / cell_area.sum())
    left, width = e[index], np.diff(e)[index]
    a, b = density[index], np.diff(density)[index]
    area = cell_area[index]
    target = rng.uniform(size=n_events) * area
    if np.any(a <= 0):
        raise ValueError("Native tail interpolation needs strictly positive left rates")
    # Invert integral_0^x [a + b*t/width] dt = target. This is the
    # linear-density inverse CDF, not a flat-in-cell approximation.
    discriminant = a ** 2 + 2 * b / width * target
    x = 2 * target / (a + np.sqrt(np.maximum(discriminant, 0.0)))
    return left + x


def native_tail_weight(diagnostics, energies, weight_scale):
    """Weight of native-tail events in a fixed-total spectrum mixture.

    The tail count is Binomial(N, R_tail/(R_prepared+R_tail)), and each
    component is sampled from its normalized native density. Every event thus
    represents exactly the same absolute rate, R_total/N. A density/proposal
    ratio must NOT be applied again: that would double-count the tail rate.
    """
    if diagnostics["tail_extension_used"] and len(diagnostics["tail_energy_kevnr"]):
        low, high = diagnostics["tail_native_rate_boundary_kevnr"]
        if np.any(~np.isfinite(energies)) or np.any(energies < low) \
                or np.any(energies > high):
            raise ValueError("Tail event lies outside the recorded native support")
    return np.full(np.shape(energies), weight_scale, dtype=float)


def expected_yield(
    rate_events_per_kg_day: float,
    exposure_tonne_years: float,
    sigma_cm2: float = REFERENCE_SIGMA_CM2,
) -> float:
    """Absolute pre-selection yield with exact cross-section/exposure scaling."""
    for value in (rate_events_per_kg_day, exposure_tonne_years, sigma_cm2):
        if not np.isfinite(value) or value < 0:
            raise ValueError("Rate, exposure and cross section must be finite, >=0")
    return (rate_events_per_kg_day * 1000.0 * DAYS_PER_YEAR *
            exposure_tonne_years * sigma_cm2 / REFERENCE_SIGMA_CM2)
