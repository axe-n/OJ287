import numpy as np
from scipy.optimize import brentq


def get_outburst_times(target_indices, params):
    t_0 = params["t_0"]
    phi_0 = params["phi_0"]
    k = params["k"]
    e_t = params["e_t"]
    delta_phi = params["delta_phi"]
    n = params["n"]
    n_dot = params["n_dot"]
    n_ddot = params["n_ddot"]

    def orbital_phase(t):
        dt = t - t_0
        l = n * dt + 0.5 * n_dot * (dt**2) + (1.0 / 6.0) * n_ddot * (dt**3)

        u = l if e_t < 0.8 else l + e_t * np.sin(l)
        for _ in range(100):
            f_val = u - e_t * np.sin(u) - l
            f_deriv = 1.0 - e_t * np.cos(u)
            u_next = u - f_val / f_deriv
            if np.abs(u_next - u) < 1e-12:
                u = u_next
                break
            u = u_next

        e_phi = e_t * (1.0 + delta_phi)
        f = 2.0 * np.arctan2(
            np.sqrt(1.0 + e_phi) * np.sin(u / 2.0),
            np.sqrt(1.0 - e_phi) * np.cos(u / 2.0),
        )
        return phi_0 + (1.0 + k) * f

    outburst_times = []
    for i in target_indices:
        def residual(t, index=i):
            return orbital_phase(t) - index * np.pi

        P = 2 * np.pi / n
        t_est = t_0 + (i * np.pi - phi_0) / (n * (1.0 + k))
        t_start = t_est - 0.25 * P
        t_end = t_est + 0.25 * P

        try:
            t_impact = brentq(residual, t_start, t_end)
            outburst_times.append((i, t_impact))
        except ValueError:
            print(f"Bracketing failed for index {i}. Adjusting search window...")

    return outburst_times

# Example Parameter Dictionary (values are mock representations of a relativistic system)
params = {
    't_0': 1886.0,       # referepnce eoch (years) [17]
    'phi_0': 0.0,        # phase at t_0 (rad) [14]
    'k': 0.11,           # periastron precession rate per orbit [18]
    'e_t': 0.38,         # time eccentricity [19]
    'delta_phi': 0.01,   # post-Newtonian radial coordinate deformation [14]
    'n': 0.52,           # mean motion (rad/year), corresponding to ~12 yr period [18]
    'n_dot': -1e-4,      # first derivative of mean motion [4]
    'n_ddot': -1e-7      # second derivative of mean motion [4]
}

# Find the 1st, 2nd, and 3rd disk-crossings (i = 1, 2, 3)
crossings = [20, 21, 22]
results = get_outburst_times(crossings, params)

for idx, t_val in results:
    print(f"Crossing index i={idx:d}: Impact Time t = {t_val:.4f} years")
