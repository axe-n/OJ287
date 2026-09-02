"""Compact model-independent binary orbital state."""

from __future__ import annotations

import math
from dataclasses import dataclass
from numbers import Real

from .physics import b_phi, comp_v, l
from .physics import phi as ph
from .solver import solve_kepler


def _fin(x: Real, name: str) -> float:
	if isinstance(x, bool) or not isinstance(x, Real):
		raise TypeError(f"{name} must be a real scalar")
	x = float(x)
	if not math.isfinite(x):
		raise ValueError(f"{name} must be finite")
	return x


def _ecc(x: Real, name: str) -> float:
	x = _fin(x, name)
	if not 0.0 <= x < 1.0:
		raise ValueError(f"{name} must satisfy 0 <= {name} < 1")
	return x


@dataclass
class BBHModel:
	"""Parameters for an elliptic binary orbital phase model.

	The mean anomaly is evolved as a cubic in time, Kepler's equation gives
	the eccentric anomaly ``u``, and ``phi = phi0 + (1 + k) * v`` gives phase.
	"""

	phi0: float
	k: float
	e_phi: float
	l0: float
	n: float
	n_dot: float
	n_ddot: float
	t0: float
	e_t: float
	beta_phi: float | None = None

	def __post_init__(self) -> None:
		self.phi0 = _fin(self.phi0, "phi0")
		self.k = _fin(self.k, "k")
		self.e_phi = _ecc(self.e_phi, "e_phi")
		self.l0 = _fin(self.l0, "l0")
		self.n = _fin(self.n, "n")
		self.n_dot = _fin(self.n_dot, "n_dot")
		self.n_ddot = _fin(self.n_ddot, "n_ddot")
		self.t0 = _fin(self.t0, "t0")
		self.e_t = _ecc(self.e_t, "e_t")
		if self.beta_phi is None:
			self.beta_phi = b_phi(self.e_phi)
		else:
			self.beta_phi = _fin(self.beta_phi, "beta_phi")
			if not 0.0 <= self.beta_phi < 1.0:
				raise ValueError("beta_phi must satisfy 0 <= beta_phi < 1")

	def mean_anomaly(self, t: float) -> float:
		"""Return ``l(t)`` from the cubic mean-anomaly expansion."""
		return l(t, self.l0, self.n, self.n_dot, self.n_ddot, self.t0)

	def true_anomaly(self, u: float) -> float:
		"""Return true anomaly ``v = (v - u) + u``."""
		return comp_v(u, self.e_phi)

	def phase_angle(self, u: float, v: float | None = None) -> float:
		"""Return phase from ``phi = phi0 + (1 + k) * v``."""
		if v is None:
			v = self.true_anomaly(u)
		return ph(self.phi0, self.k, v)

	def compute_phi(self, t: float, u: float | None = None) -> float:
		"""Solve for ``u`` at time ``t`` and return the orbital phase."""
		if u is None:
			u = solve_kepler(self.mean_anomaly(t), self.e_t)
		return self.phase_angle(u)
