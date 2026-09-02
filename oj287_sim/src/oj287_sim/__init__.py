"""Model-independent tools for binary black hole orbital simulations."""

from .model import BBHModel
from .main import simulate_orbit
from .physics import phi as compute_phase
from .solver import solve_kepler


__version__ = "0.1.0"

__all__ = [
	"BBHModel",
	"compute_phase",
	"simulate_orbit",
	"solve_kepler",
	"__version__",
]
