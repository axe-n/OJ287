"""Run an orbital model over a time grid."""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from .model import BBHModel
# from .solver import solve_kepler
from .solver import solve_kepler_danby_array


def _cross(model: BBHModel, t: np.ndarray, p: np.ndarray) -> np.ndarray:
	x = []
	for j, (a, b) in enumerate(zip(p[:-1], p[1:])):
		if a == b:
			continue
		for i in range(math.ceil(min(a, b) / math.pi), math.floor(max(a, b) / math.pi) + 1):
			q = i * math.pi
			if a <= q <= b or b <= q <= a:
				if math.isclose(a, q, rel_tol=0.0, abs_tol=1e-12):
					x.append(t[j])
				elif math.isclose(b, q, rel_tol=0.0, abs_tol=1e-12):
					x.append(t[j + 1])
				else:
					x.append(brentq(
						lambda z: model.compute_phi(z) - q,
						t[j],
						t[j + 1],
						xtol=1e-12,
					))
	return np.unique(x)


def simulate_orbit(model: BBHModel, t_values: np.ndarray) -> dict:
	"""Compute l, u, v, phi, and times where phi = i*pi."""
	t = np.asarray(t_values, dtype=float)
	if t.ndim != 1:
		raise ValueError("t_values must be one-dimensional")
	if not np.all(np.isfinite(t)):
		raise ValueError("t_values must be finite")
	if t.size > 1 and np.any(np.diff(t) < 0.0):
		raise ValueError("t_values must be ordered")

	dt = t - model.t0
	la = (
		model.l0
		+ model.n * dt
		+ 0.5 * model.n_dot * dt**2
		+ (1.0 / 6.0) * model.n_ddot * dt**3
	)
	# u = np.array([solve_kepler(q, model.e_t) for q in la])
	u = solve_kepler_danby_array(la, model.e_t)
	v = u + 2.0 * np.arctan2(
		model.beta_phi * np.sin(u),
		1.0 - model.beta_phi * np.cos(u),
	)
	p = model.phi0 + (1.0 + model.k) * v
	return {
		"t": t,
		"l": la,
		"u": u,
		"v": v,
		"phi": p,
		"phase_crossings": _cross(model, t, p),
		"period": model.Per(),
	}


if __name__ == "__main__":
	m = BBHModel(0.0, 0.0, 0.3, 0.0, 1.0, 0.0, 0.0, 0.0, 0.3)
	r = simulate_orbit(m, np.linspace(0.0, 2.0 * math.pi, 1001))
	print(f"samples: {len(r['t'])}")
	print(f"phase crossings: {len(r['phase_crossings'])}")
