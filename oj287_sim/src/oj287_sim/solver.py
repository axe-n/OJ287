"""Stable scalar solvers for Kepler's equation."""

from __future__ import annotations

import math
from numbers import Real


def _fin(x: Real, name: str) -> float:
	if isinstance(x, bool) or not isinstance(x, Real):
		raise TypeError(f"{name} must be a real scalar")
	x = float(x)
	if not math.isfinite(x):
		raise ValueError(f"{name} must be finite")
	return x


def _ecc(e: Real) -> float:
	e = _fin(e, "e")
	if not 0.0 <= e < 1.0:
		raise ValueError("e must satisfy 0 <= e < 1")
	return e


def _args(l: Real, e: Real, tol: Real, max_iter: int) -> tuple[float, float, float, int]:
	l = _fin(l, "l")
	e = _ecc(e)
	tol = _fin(tol, "tol")
	if tol <= 0.0:
		raise ValueError("tol must be positive")
	if isinstance(max_iter, bool) or not isinstance(max_iter, int):
		raise TypeError("max_iter must be an integer")
	if max_iter < 1:
		raise ValueError("max_iter must be positive")
	return l, e, tol, max_iter


def _u0(l: float, e: float) -> float:
	"""Mikkola's cubic approximation and correction."""
	a = (1.0 - e) / (4.0 * e + 0.5)
	b = l / (2.0 * (4.0 * e + 0.5))
	d = math.sqrt(b * b + a * a * a)
	q = b + d
	if q == 0.0:
		q = b - d
	z = math.copysign(abs(q) ** (1.0 / 3.0), q)
	s = 0.0 if z == 0.0 else z - a / z
	w = s - 0.078 * s**5 / (1.0 + e)
	return l + e * (3.0 * w - 4.0 * w**3)


def _f(u: float, l: float, e: float) -> float:
	return u - e * math.sin(u) - l


def solve_kepler_mikkola(
	l: float,
	e: float,
	tol: float = 1e-12,
	max_iter: int = 50,
) -> float:
	"""Solve ``l = u - e*sin(u)`` from Mikkola's approximation.

	The cubic in ``s = sin(u/3)`` is normalized as
	``s**3 + 3*a*s - 2*b = 0``, where
	``a = (1-e)/(4*e+1/2)`` and ``b = l/(2*(4*e+1/2))``.
	A short Newton refinement makes the returned root meet ``tol``.
	"""
	l, e, tol, max_iter = _args(l, e, tol, max_iter)
	u = _u0(l, e)
	for _ in range(max_iter):
		f = _f(u, l, e)
		if abs(f) <= tol:
			return u
		u -= f / (1.0 - e * math.cos(u))
	if abs(_f(u, l, e)) <= tol:
		return u
	raise RuntimeError("Mikkola solver did not converge")


def solve_kepler_danby(
	l: float,
	e: float,
	tol: float = 1e-12,
	max_iter: int = 50,
) -> float:
	"""Solve Kepler's equation with Mikkola initialization and Danby steps."""
	l, e, tol, max_iter = _args(l, e, tol, max_iter)
	u = _u0(l, e)
	for _ in range(max_iter):
		f = _f(u, l, e)
		if abs(f) <= tol:
			return u
		fp = 1.0 - e * math.cos(u)
		fpp = e * math.sin(u)
		f3 = e * math.cos(u)
		f4 = -fpp
		u1 = -f / fp
		u2 = -f / (fp + 0.5 * fpp * u1)
		u3 = -f / (fp + 0.5 * fpp * u2 + f3 * u2**2 / 6.0)
		u4 = -f / (
			fp + 0.5 * fpp * u3 + f3 * u3**2 / 6.0 + f4 * u3**3 / 24.0
		)
		u += u4
	if abs(_f(u, l, e)) <= tol:
		return u
	raise RuntimeError("Danby solver did not converge")


def solve_kepler(
	l: float,
	e: float,
	tol: float = 1e-12,
	max_iter: int = 50,
) -> float:
	"""Return the default Mikkola-Danby solution of Kepler's equation."""
	return solve_kepler_danby(l, e, tol, max_iter)
