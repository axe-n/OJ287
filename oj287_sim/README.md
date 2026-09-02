# oj287_sim

Small, model-independent Python tools for simulating binary black hole orbital
phase evolution. The code evolves the mean anomaly, solves Kepler's equation,
computes the true anomaly and phase, and estimates times where the phase is an
integer multiple of pi.

## Equations

For `dt = t - t0`, the mean anomaly is

```text
l = l0 + n*dt + 1/2*n_dot*dt^2 + 1/6*n_ddot*dt^3
```

The eccentric anomaly `u` satisfies Kepler's equation:

```text
l = u - e_t*sin(u)
```

The true anomaly is evaluated using

```text
v - u = 2*atan(beta_phi*sin(u)/(1 - beta_phi*cos(u)))
v = (v - u) + u
```

where

```text
beta_phi = (1 - sqrt(1 - e_phi^2))/e_phi
```

Finally,

```text
phi = phi0 + (1 + k)*v
```

All angles use radians.

## Install

From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .
```

## Run

Run the small built-in example:

```bash
python -m oj287_sim.main
```

The notebook demonstration is in `notebooks/demo.ipynb`.

## Solver

`solver.py` uses Mikkola's cubic approximation for an initial value of `u`.
Cardano's formula solves the normalized cubic in `sin(u/3)`, followed by
Mikkola's correction. Danby iterations then reduce the residual of
`u - e*sin(u) - l` to the requested tolerance.

## Workflow

1. Define the orbital parameters with `BBHModel`.
2. Generate a time grid.
3. Compute `l(t)` and solve for `u`.
4. Compute `v(t)` and `phi(t)`.
5. Estimate `phi = i*pi` event times by interpolation between samples.
6. Test numerical residuals before using results for analysis.

## Get the project

Clone the repository and enter its root directory:

```bash
git clone <repository-url>
cd oj287_sim
```

The root directory is the directory containing `pyproject.toml`.

## Package layout

```text
oj287_sim/
├── README.md
├── requirements.txt
├── pyproject.toml
├── notebooks/
│   └── demo.ipynb
├── src/
│   └── oj287_sim/
│       ├── __init__.py
│       ├── model.py
│       ├── physics.py
│       ├── solver.py
│       └── main.py
└── tests/
		├── test_solver.py
		└── test_physics.py
```

The `src/oj287_sim` directory is the importable package. The `src` layout
keeps package code separate from tests, notebooks, and project files.

## Use the Python API

```python
import numpy as np

from oj287_sim import BBHModel, simulate_orbit

m = BBHModel(
		phi0=0.0,
		k=0.02,
		e_phi=0.3,
		l0=0.0,
		n=1.0,
		n_dot=0.0,
		n_ddot=0.0,
		t0=0.0,
		e_t=0.3,
)

t = np.linspace(0.0, 20.0, 2001)
r = simulate_orbit(m, t)

print(r["phi"])
print(r["phase_crossings"])
```

The returned dictionary contains:

- `t`: input time values
- `l`: mean anomaly at each time
- `u`: eccentric anomaly from Kepler's equation
- `v`: true anomaly
- `phi`: orbital phase
- `phase_crossings`: interpolated times where `phi = i*pi`

For direct access to individual operations:

```python
from oj287_sim import compute_phase, solve_kepler

u = solve_kepler(l=1.0, e=0.3)
phase = compute_phase(phi0=0.0, k=0.02, v=1.5)
```

## Run the example

After installation, run the built-in example from the project root:

```bash
python -m oj287_sim.main
```

It prints the number of simulated samples and detected phase crossings.

## Run the tests

```bash
pytest -q
```

The tests check exact circular-orbit behavior, Kepler residuals, the phase
formula, monotonic behavior of `beta_phi`, and known circular phase crossings.

## Open the notebook

Install Jupyter if it is not already available, then start it from the project
root:

```bash
python -m pip install jupyter
jupyter notebook notebooks/demo.ipynb
```

The notebook shows model setup, a direct Kepler solve, a complete simulation,
crossing detection, and one plot of `phi(t)`.

## Parameter rules

- `e_t` and `e_phi` must satisfy `0 <= e < 1`.
- All scalar parameters must be finite real numbers.
- Angular parameters and outputs are in radians.
- `beta_phi` is computed automatically from `e_phi` unless supplied.
- The time array should be one-dimensional and ordered for meaningful event
	interpolation.

## Reproducibility

Record the parameter values, time-grid spacing, solver tolerance, and package
version with any scientific result. The default solver tolerance is `1e-12`.
The detected crossing times depend on the time-grid spacing because crossings
are initially estimated by linear interpolation between samples.

## Development install

For local development, install the package in editable mode:

```bash
python -m pip install -e .
```

Changes under `src/oj287_sim` are then available immediately to Python, the
tests, and the notebook.
