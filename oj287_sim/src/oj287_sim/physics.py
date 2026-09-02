"""Simple, model-independent orbital utilities. Angles are in radians."""

from __future__ import annotations

import math
from numbers import Real


__all__ = [
	"b_phi",
	"v",
	"vu",
	"comp_v",
	"phi",
	"l",
]


def _fin(x: Real, name: str) -> float:
	if isinstance(x, bool) or not isinstance(x, Real):
		raise TypeError(f"{name} must be a real scalar")
	x = float(x)
	if not math.isfinite(x):
		raise ValueError(f"{name} must be finite")
	return x


def _ecc(e: Real, name: str = "e") -> float:
	e = _fin(e, name)
	if not 0.0 <= e < 1.0:
		raise ValueError(f"{name} must satisfy 0 <= {name} < 1")
	return e


def b_phi(e: float) -> float:
	"""Return beta_phi = (1 - sqrt(1 - e^2)) / e."""
	e = _ecc(e, "e_phi")
	if e == 0.0:
		return 0.0
	return e / (1.0 + math.sqrt(1.0 - e**2))


def v(u: float, e: float) -> float:
	"""Return v from u using a stable form of the true-anomaly equation."""
	u = _fin(u, "u")
	e = _ecc(e, "e_phi")
	h = 0.5 * u
	return 2.0 * math.atan2(
		math.sqrt(1.0 + e) * math.sin(h),
		math.sqrt(1.0 - e) * math.cos(h),
	)


def vu(u: float, b: float) -> float:
	"""Return v - u = 2 atan(b sin(u) / (1 - b cos(u)))."""
	u = _fin(u, "u")
	b = _fin(b, "beta_phi")
	if not 0.0 <= b < 1.0:
		raise ValueError("beta_phi must satisfy 0 <= beta_phi < 1")
	return 2.0 * math.atan2(b * math.sin(u), 1.0 - b * math.cos(u))


def comp_v(u: float, e: float) -> float:
	"""Return v = (v - u) + u."""
	u = _fin(u, "u")
	return u + vu(u, b_phi(_ecc(e, "e_phi")))


def phi(p0: float, k: float, x: float) -> float:
	"""Return phi = phi0 + (1 + k) v."""
	return _fin(p0, "phi0") + (1.0 + _fin(k, "k")) * _fin(x, "v")


def l(
	t: float,
	l0: float,
	n: float,
	nd: float,
	ndd: float,
	t0: float,
) -> float:
	"""Return l = l0 + n dt + 1/2 nd dt^2 + 1/6 ndd dt^3."""
	dt = _fin(t, "t") - _fin(t0, "t0")
	return (
		_fin(l0, "l0")
		+ _fin(n, "n") * dt
		+ 0.5 * _fin(nd, "n_dot") * dt**2
		+ (1.0 / 6.0) * _fin(ndd, "n_ddot") * dt**3
	)
