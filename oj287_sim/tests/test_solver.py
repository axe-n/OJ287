import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from oj287_sim.solver import solve_kepler, solve_kepler_danby, solve_kepler_mikkola


def test_kepler_circular_case():
	# For e=0, Kepler's equation has the exact solution u=l.
	l = 1.2
	assert solve_kepler(l, 0.0) == pytest.approx(l, abs=1e-14)


def test_kepler_eccentric_residual():
	# The numerical solution must make F(u) = u - e*sin(u) - l nearly zero.
	l, e = 1.0, 0.6
	for fn in (solve_kepler_mikkola, solve_kepler_danby, solve_kepler):
		u = fn(l, e)
		assert u - e * __import__("math").sin(u) == pytest.approx(l, abs=1e-12)
