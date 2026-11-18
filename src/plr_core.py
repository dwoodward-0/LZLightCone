"""Binned 2D profile-likelihood core for the (S1c, log10 S2c) sensitivity.

Model, per analysis bin i (bins with any signal or background support):
    nu_i(mu, theta) = mu * s_i + sum_j theta_j * f_ji + c_i
s_i: signal shape normalized to 1 (mu is in signal events in the ROI);
f_ji: unit-area shape of background source j; theta_j: its expected count
with a Gaussian constraint N(a_j | theta_j, sigma_j) and theta_j >= 0;
c_i: optional fixed (unconstrained) background, used only by the
one-bin tests. a_j are the global observables (Table I means for the
nominal model; fluctuated in toys).

NLL(mu, theta) = sum_i [nu_i - n_i ln nu_i] + sum_j (theta_j - a_j)^2 / (2 sigma_j^2)

Test statistic: two-sided q_mu = 2 [NLL(mu, theta_hat_hat) - NLL(mu_hat, theta_hat)]
with the physical bound mu_hat >= 0 (likelihood-ratio ordering with a bounded
best fit, the Feldman-Cousins construction in one bin without nuisances).

Fits: L-BFGS-B with analytic gradients in scaled parameters
z_0 = mu / mu_scale, z_j = theta_j / sigma_j.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

NU_FLOOR = 1e-30  # keeps n/nu finite when a bin's only source sits at theta = 0
THETA_START_FLOOR = 0.1  # starts theta >= 0.1 sigma: a negative global observable must not start a
                         # source at exactly 0 while the toy has counts in bins only it populates


@dataclass
class FitResult:
    nll: float
    mu: float
    theta: np.ndarray
    success: bool
    nfev: int


class BinnedLikelihood:
    def __init__(self, signal, shapes, means, sigmas, fixed_background=None, mu_scale=1.0):
        signal = np.asarray(signal, float).ravel()
        shapes = np.atleast_2d(np.asarray(shapes, float)).reshape(len(means), -1) if len(means) else \
            np.zeros((0, signal.size))
        fixed = np.zeros(signal.size) if fixed_background is None else np.asarray(fixed_background, float).ravel()
        support = (signal > 0) | (fixed > 0) | (shapes > 0).any(axis=0)
        if signal.sum() <= 0:
            raise ValueError("signal template is empty")
        self.support = support
        self.s = signal[support] / signal.sum()
        self.signal_roi_fraction = float(self.s.sum())  # 1 unless signal lies outside support (never)
        norms = shapes.sum(axis=1)
        if (norms <= 0).any():
            raise ValueError("every background shape needs positive area")
        self.f = shapes[:, support] / norms[:, None]
        self.c = fixed[support]
        self.means = np.asarray(means, float)
        self.sigmas = np.asarray(sigmas, float)
        if (self.sigmas <= 0).any():
            raise ValueError("constraint widths must be positive")
        self.J = len(self.means)
        self.mu_scale = float(mu_scale)

    # ---- model -----------------------------------------------------------
    def expected(self, mu, theta):
        return mu * self.s + theta @ self.f + self.c

    def project(self, hist):
        """Restrict a full-grid histogram (or count array) to the analysis bins."""
        return np.asarray(hist).ravel()[self.support]

    def _unpack(self, z):
        return z[0] * self.mu_scale, z[1:] * self.sigmas

    def _nll_grad(self, z, n, a, mu_fixed):
        if mu_fixed is None:
            mu, theta = self._unpack(z)
        else:
            mu, theta = mu_fixed, z * self.sigmas
        nu = np.maximum(self.expected(mu, theta), NU_FLOOR)
        nll = float(nu.sum() - n @ np.log(nu))
        r = 1.0 - n / nu
        pull = (theta - a) / self.sigmas
        nll += 0.5 * float(pull @ pull)
        g_theta = (self.f @ r + pull / self.sigmas) * self.sigmas
        if mu_fixed is not None:
            return nll, g_theta
        g_mu = float(self.s @ r) * self.mu_scale
        return nll, np.concatenate([[g_mu], g_theta])

    def nll(self, mu, theta, n, a=None):
        a = self.means if a is None else a
        z = np.concatenate([[mu / self.mu_scale], np.asarray(theta) / self.sigmas])
        return self._nll_grad(z, n, a, None)[0]

    # ---- fits ------------------------------------------------------------
    def fit(self, n, a=None, mu_fixed=None, start=None):
        a = self.means if a is None else np.asarray(a, float)
        n = np.asarray(n, float)
        theta0 = np.maximum(a if start is None else start[1], THETA_START_FLOOR * self.sigmas)
        if mu_fixed is None:
            mu0 = 0.0 if start is None else max(start[0], 0.0)
            if start is None:  # crude excess-over-background start
                mu0 = max(float(n.sum() - theta0.sum() - self.c.sum()), 0.0)
            z0 = np.concatenate([[mu0 / self.mu_scale], theta0 / self.sigmas])
            bounds = [(0.0, None)] * (1 + self.J)
        else:
            z0 = theta0 / self.sigmas
            bounds = [(0.0, None)] * self.J
        if len(z0) == 0:  # mu fixed and no nuisances
            nu = np.maximum(self.expected(mu_fixed, np.zeros(0)), NU_FLOOR)
            return FitResult(float(nu.sum() - n @ np.log(nu)), float(mu_fixed), np.zeros(0), True, 1)
        def run(z):
            r = minimize(self._nll_grad, z, args=(n, a, mu_fixed), jac=True, method="L-BFGS-B",
                         bounds=bounds, options={"maxiter": 2000, "ftol": 1e-13, "gtol": 1e-9})
            # On ABNORMAL line-search exits scipy's r.fun can be the last TRIAL value, not f(r.x).
            return r, self._nll_grad(r.x, n, a, mu_fixed)[0]

        res, fun = run(z0)
        nfev, success = res.nfev, res.success
        if not success:
            # Restarts: a mu = 0 start from the global observables (the crude excess start can sit
            # hundreds of events off when constraint widths are large), and a warm restart from
            # the stalled point. Keep the lowest NLL; success if any attempt converged.
            alts = [res.x]
            if mu_fixed is None:
                alts.append(np.concatenate([[0.0], np.maximum(a, THETA_START_FLOOR * self.sigmas) / self.sigmas]))
            for z in alts:
                r2, f2 = run(z)
                nfev += r2.nfev
                success = success or r2.success
                if f2 < fun:
                    res, fun = r2, f2
        if not np.isfinite(fun):
            raise FloatingPointError(f"non-finite NLL after fit (mu_fixed={mu_fixed})")
        if mu_fixed is None:
            mu, theta = self._unpack(res.x)
        else:
            mu, theta = float(mu_fixed), res.x * self.sigmas
        return FitResult(float(fun), float(mu), theta, bool(success), int(nfev))

    def q_mu(self, n, mu, a=None, free=None):
        """Two-sided profile-likelihood-ratio statistic (>= 0) and the free fit used."""
        free = self.fit(n, a) if free is None else free
        cond = self.fit(n, a, mu_fixed=mu, start=(mu, free.theta))
        return max(2.0 * (cond.nll - free.nll), 0.0), free

    # ---- toys ------------------------------------------------------------
    def toy(self, rng, mu, theta_true=None):
        """Poisson data in analysis bins + Gaussian-fluctuated global observables."""
        theta_true = self.means if theta_true is None else theta_true
        n = rng.poisson(self.expected(mu, theta_true)).astype(float)
        a = rng.normal(theta_true, self.sigmas)
        return n, a


def neyman_upper_limit(mu_grid, q_obs, q_crit):
    """Largest-mu crossing of q_obs(mu) - q_crit(mu) from below (linear interpolation).

    Returns the upper end of the region {q_obs <= q_crit} that contains its minimum.
    """
    d = np.asarray(q_obs) - np.asarray(q_crit)
    i0 = int(np.argmin(d))
    if d[i0] > 0:
        return float(mu_grid[0])  # nothing accepted on the grid: limit below grid start
    for i in range(i0, len(mu_grid) - 1):
        if d[i] <= 0 < d[i + 1]:
            return float(mu_grid[i] + (mu_grid[i + 1] - mu_grid[i]) * (-d[i]) / (d[i + 1] - d[i]))
    return float("inf")  # still accepted at the top of the grid
