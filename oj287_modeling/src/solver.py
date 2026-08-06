import numpy as np
from scipy.optimize import brentq
from .orbit import orbital_phase


def impact_time_for_index(i, params):
    t_0 = params['t_0']
    phi_0 = params['phi_0']
    k = params['k']
    e_t = params['e_t']
    delta_phi = params['delta_phi']
    n = params['n']
    n_dot = params['n_dot']
    n_ddot = params['n_ddot']

    def residual(t):
        return orbital_phase(t, t_0, phi_0, k, e_t, delta_phi, n, n_dot, n_ddot) - i * np.pi

    P = 2 * np.pi / n
    t_est = t_0 + (i * np.pi - phi_0) / (n * (1.0 + k))
    t_start = t_est - 0.25 * P
    t_end = t_est + 0.25 * P

    return brentq(residual, t_start, t_end)


def get_outburst_times(target_indices, params):
    outburst_times = []
    for i in target_indices:
        try:
            t_impact = impact_time_for_index(i, params)
            outburst_times.append((i, t_impact))
        except ValueError:
            print(f'Bracketing failed for index {i}. Adjusting search window...')
    return outburst_times
