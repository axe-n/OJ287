import numpy as np


def mean_anomaly(t, t_0, n, n_dot, n_ddot):
    dt = t - t_0
    return n * dt + 0.5 * n_dot * (dt**2) + (1.0 / 6.0) * n_ddot * (dt**3)


def solve_kepler(l, e_t, tol=1e-12, max_iter=100):
    u = l if e_t < 0.8 else l + e_t * np.sin(l)
    for _ in range(max_iter):
        f = u - e_t * np.sin(u) - l
        df = 1.0 - e_t * np.cos(u)
        u_next = u - f / df
        if np.abs(u_next - u) < tol:
            return u_next
        u = u_next
    return u


def true_anomaly(u, e_t, delta_phi):
    e_phi = e_t * (1.0 + delta_phi)
    return 2.0 * np.arctan2(
        np.sqrt(1.0 + e_phi) * np.sin(u / 2.0),
        np.sqrt(1.0 - e_phi) * np.cos(u / 2.0),
    )


def orbital_phase(t, t_0, phi_0, k, e_t, delta_phi, n, n_dot, n_ddot):
    l = mean_anomaly(t, t_0, n, n_dot, n_ddot)
    u = solve_kepler(l, e_t)
    f = true_anomaly(u, e_t, delta_phi)
    return phi_0 + (1.0 + k) * f
