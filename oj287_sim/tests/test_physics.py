import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from oj287_sim.main import simulate_orbit
from oj287_sim.model import BBHModel
from oj287_sim.physics import b_phi, phi


def test_phase_formula():
	# Check the model-independent phase relation directly.
	assert phi(0.4, 0.2, 1.5) == pytest.approx(0.4 + 1.2 * 1.5)


def test_beta_is_increasing():
	# beta_phi must increase from zero for increasing elliptic eccentricity.
	e = np.array([0.0, 0.2, 0.5, 0.8])
	b = np.array([b_phi(x) for x in e])
	assert b[0] == 0.0
	assert np.all(np.diff(b) > 0.0)


def test_phase_crossings():
	# With a circular orbit and k=0, phi(t)=t, so crossings are known.
	m = BBHModel(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
	t = np.linspace(0.0, 2.0 * np.pi, 5)
	r = simulate_orbit(m, t)
	assert r["phase_crossings"] == pytest.approx([0.0, np.pi, 2.0 * np.pi])
	assert r["period"] == pytest.approx(2.0 * np.pi)


def test_period_uses_model_n():
	m = BBHModel(0.0, 0.0, 0.2, 0.0, 0.8, 0.0, 0.0, 0.0, 0.3)
	r = simulate_orbit(m, np.array([0.0, 1.0]))
	assert r["period"] == pytest.approx(2.0 * np.pi / 0.8)


def test_time_values_must_be_ordered():
	m = BBHModel(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
	with pytest.raises(ValueError, match="ordered"):
		simulate_orbit(m, np.array([0.0, 2.0, 1.0]))


def test_beta_phi_is_derived():
	m = BBHModel(0.0, 0.0, 0.3, 0.0, 1.0, 0.0, 0.0, 0.0, 0.3)
	assert m.beta_phi == pytest.approx(b_phi(0.3))
