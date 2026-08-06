from pathlib import Path
from oj287_modeling.src.solver import get_outburst_times

PARAMS = {
    't_0': 1886.0,
    'phi_0': 0.0,
    'k': 0.11,
    'e_t': 0.38,
    'delta_phi': 0.01,
    'n': 0.52,
    'n_dot': -1e-4,
    'n_ddot': -1e-7,
}

CROSSING_INDICES = [20, 21, 22]


def main():
    root = Path(__file__).resolve().parent
    data_file = root / 'oj287_modeling' / 'data' / 'outbursts.csv'
    if data_file.exists():
        print(f'Loaded observed outburst data from {data_file}')

    results = get_outburst_times(CROSSING_INDICES, PARAMS)
    for idx, t_val in results:
        print(f'Crossing index i={idx:d}: Impact Time t = {t_val:.4f} years')


if __name__ == '__main__':
    main()
