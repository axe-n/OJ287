"""Run an orbital model over a time grid."""

from __future__ import annotations

import math

import numpy as np

from .model import BBHModel
from .solver import solve_kepler


def _cross(t: np.ndarray, p: np.ndarray) -> np.ndarray:
	x = []
	for j, (a, b) in enumerate(zip(p[:-1], p[1:])):
		if a == b:
			continue
		for i in range(math.ceil(min(a, b) / math.pi), math.floor(max(a, b) / math.pi) + 1):
			q = i * math.pi
			if a <= q <= b or b <= q <= a:
				x.append(t[j] + (q - a) * (t[j + 1] - t[j]) / (b - a))
	return np.unique(x)


def simulate_orbit(model: BBHModel, t_values: np.ndarray) -> dict:
	"""Compute l, u, v, phi, and times where phi = i*pi."""
	t = np.asarray(t_values, dtype=float)
	if t.ndim != 1:
		raise ValueError("t_values must be one-dimensional")

	la = np.array([model.mean_anomaly(q) for q in t])
	u = np.array([solve_kepler(q, model.e_t) for q in la])
	v = np.array([model.true_anomaly(q) for q in u])
	p = np.array([model.phase_angle(a, b) for a, b in zip(u, v)])
	return {
		"t": t,
		"l": la,
		"u": u,
		"v": v,
		"phi": p,
		"phase_crossings": _cross(t, p),
	}


if __name__ == "__main__":
	m = BBHModel(0.0, 0.0, 0.3, 0.0, 1.0, 0.0, 0.0, 0.0, 0.3)
	r = simulate_orbit(m, np.linspace(0.0, 2.0 * math.pi, 1001))
	print(f"samples: {len(r['t'])}")
	print(f"phase crossings: {len(r['phase_crossings'])}")
