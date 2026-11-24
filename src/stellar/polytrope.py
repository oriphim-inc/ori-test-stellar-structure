"""Lane-Emden polytropes and the Eddington standard (n = 3) model used as initial guess."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import physics as ph


def lane_emden(n, xi_max=20.0, rtol=1e-11, atol=1e-13):
    """Integrate theta'' + 2/xi theta' + theta^n = 0 to the first zero.

    Returns dict with xi1, dtheta at xi1, the dense solution and the derived constants.
    """
    xi0 = 1e-6
    th0 = 1 - xi0**2 / 6 + n * xi0**4 / 120
    dth0 = -xi0 / 3 + n * xi0**3 / 30

    def rhs(xi, y):
        th, dth = y
        return [dth, -np.sign(th) * abs(th) ** n - 2 * dth / xi]

    def surf(xi, y):
        return y[0]
    surf.terminal = True
    surf.direction = -1
    sol = solve_ivp(rhs, (xi0, xi_max), [th0, dth0], rtol=rtol, atol=atol,
                    events=surf, dense_output=True, method="DOP853")
    xi1 = sol.t_events[0][0]
    dth1 = sol.y_events[0][0][1]
    return {
        "n": n,
        "xi1": xi1,
        "dtheta1": dth1,
        "mass_const": -xi1**2 * dth1,             # -xi^2 theta' at xi1
        "rhoc_over_rhomean": -xi1 / (3 * dth1),
        "sol": sol,
    }


def eddington_beta(M, mu):
    """Eddington quartic 1-beta = 0.00298 (M/Msun)^2 mu^4 beta^4 (constant from n=3 constants)."""
    le = lane_emden(3)
    # exact coefficient: (1-beta)/(mu^4 beta^4) = (M/M3)^2  with
    # M3 = (4/pi) ... ; compute directly: M = 4 pi (K/(pi G))^(3/2) * mass_const, K from beta
    C = (3 * (ph.k_B / ph.m_u) ** 4 / ph.a_rad) ** 0.5
    M3 = 4 / np.sqrt(np.pi) * C / ph.G**1.5 * le["mass_const"]
    # K = [(k/(mu m_u))^4 3(1-beta)/(a beta^4)]^(1/3); M = 4 pi (K/(pi G))^(3/2) mass_const
    f = lambda b: (1 - b) / (mu**4 * b**4) - (M / M3) ** 2
    return brentq(f, 1e-6, 1 - 1e-14)


def polytrope_guess(M, R, comp, n=3):
    """Central P, rho, T of an n-polytrope of mass M and radius R (with Eddington beta for n=3)."""
    le = lane_emden(n)
    rho_mean = 3 * M / (4 * np.pi * R**3)
    rhoc = le["rhoc_over_rhomean"] * rho_mean
    alpha = R / le["xi1"]
    Pc = 4 * np.pi * ph.G * rhoc**2 * alpha**2 / (n + 1)
    b = eddington_beta(M, comp.mu) if n == 3 else 1.0
    Tc = b * Pc * comp.mu * ph.m_u / (ph.k_B * rhoc)
    return {"Pc": Pc, "rhoc": rhoc, "Tc": Tc, "beta": b, "lane_emden": le}


def polytrope_profile(M, R, comp, n=3, npts=400):
    """Full n=3 Eddington standard model profile (r, m, rho, P, T)."""
    g = polytrope_guess(M, R, comp, n)
    le = g["lane_emden"]
    xi = np.linspace(1e-6, le["xi1"], npts)
    th, dth = le["sol"].sol(xi)
    th = np.clip(th, 0, None)
    alpha = R / le["xi1"]
    rho = g["rhoc"] * th**n
    P = g["Pc"] * th ** (n + 1)
    T = g["Tc"] * th
    m = 4 * np.pi * alpha**3 * g["rhoc"] * (-xi**2 * dth)
    return {"r": alpha * xi, "m": m, "rho": rho, "P": P, "T": T, **g}
