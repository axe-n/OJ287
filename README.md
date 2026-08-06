# OJ287
Major works on the OJ287 SMBH binary blazar system.

## Project structure

- `run_pipeline.py` — main runner script for the modeling pipeline.
- `oj287_modeling/` — package root.
  - `data/` — observed outburst data files, including `outbursts.csv`.
  - `notebooks/` — interactive notebooks for analysis and visualization.
  - `src/` — core pipeline modules:
    - `orbit.py` — analytical orbit functions and phase calculation.
    - `solver.py` — solver for disk-impact times where `phi(t) = i * pi`.
    - `fit_engine.py` — fit metrics such as chi-squared and likelihood.
- `test.py` — legacy scratch script. The production logic has been moved into `oj287_modeling/src/` and `run_pipeline.py`, so it can be removed if you no longer need it.

## How to use

1. Install dependencies:

```bash
pip install numpy scipy
```

2. Run the pipeline:

```bash
python run_pipeline.py
```

3. Use the package directly from Python:

```python
from oj287_modeling.src.solver import get_outburst_times
params = {
    't_0': 1886.0,
    'phi_0': 0.0,
    'k': 0.11,
    'e_t': 0.38,
    'delta_phi': 0.01,
    'n': 0.52,
    'n_dot': -1e-4,
    'n_ddot': -1e-7,
}
indices = [20, 21, 22]
results = get_outburst_times(indices, params)
```

## Notes

- `run_pipeline.py` is the current entrypoint for the pipeline.
- `oj287_modeling/src/` contains the core logic.
- `oj287_modeling/data/outbursts.csv` contains the observed outburst table.
