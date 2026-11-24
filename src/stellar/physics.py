"""Microphysics: constants, EOS, opacity, nuclear energy generation, Gamow peak (cgs units).

Sources are given in PLAN.md (Kippenhahn, Weigert & Weiss 2012 = KWW; Carroll & Ostlie 2007 = C&O;
Hansen, Kawaler & Trimble 2004 = HKT).
"""
from dataclasses import dataclass
import numpy as np

# CODATA 2018 / IAU 2015 nominal values (cgs)
G = 6.67430e-8
c = 2.99792458e10
k_B = 1.380649e-16
m_u = 1.66053906660e-24
sigma_SB = 5.670374419e-5
a_rad = 4.0 * sigma_SB / c
h_planck = 6.62607015e-27
hbar = h_planck / (2 * np.pi)
e_charge = 4.80320471e-10          # esu
alpha_fs = e_charge**2 / (hbar * c)
keV = 1.602176634e-9               # erg
MSUN = 1.98847e33                  # IAU 2015 nominal GM / G
RSUN = 6.957e10
LSUN = 3.828e33


@dataclass(frozen=True)
class Composition:
    X: float = 0.70
    Z: float = 0.02
    f_cno: float = 0.7   # fraction of Z in CNO nuclei (assumption, PLAN.md)

    @property
    def Y(self):
        return 1.0 - self.X - self.Z

    @property
    def mu(self):
        """Mean molecular weight, fully ionized."""
        return 1.0 / (2 * self.X + 0.75 * self.Y + 0.5 * self.Z)

    @property
    def X_cno(self):
        return self.f_cno * self.Z


# ---------------------------------------------------------------- EOS
def density(P, T, comp):
    """rho from total P (gas + radiation). Returns negative/zero if P_gas <= 0."""
    Pgas = P - a_rad * T**4 / 3.0
    return Pgas * comp.mu * m_u / (k_B * T)


def pressure(rho, T, comp):
    return rho * k_B * T / (comp.mu * m_u) + a_rad * T**4 / 3.0


def beta(P, T):
    return 1.0 - a_rad * T**4 / (3.0 * P)


def nabla_ad(b):
    """Adiabatic gradient for ideal gas + radiation (KWW ch. 13)."""
    q = (1.0 - b) * (4.0 + b) / b**2
    return (1.0 + q) / (2.5 + 4.0 * q)


# ---------------------------------------------------------------- opacity
def kappa_kramers(rho, T, comp):
    """Bound-free (with guillotine factor) + free-free Kramers opacities (C&O eqs. 9.22-9.23)."""
    t_guil = 2.82 * (np.maximum(rho, 1e-30) * 1e3 * (1 + comp.X)) ** 0.2   # rho in kg/m^3 in the fit
    kbf = 4.34e25 / t_guil * comp.Z * (1 + comp.X) * rho * T**-3.5
    kff = 3.68e22 * (1 - comp.Z) * (1 + comp.X) * rho * T**-3.5
    return kbf + kff


def kappa_es(comp):
    return 0.2 * (1 + comp.X)


def kappa_hminus(rho, T, comp):
    """H^- opacity (HKT eq. 4.94; C&O eq. 9.24 converted to cgs)."""
    return 2.5e-31 * (comp.Z / 0.02) * np.sqrt(np.maximum(rho, 1e-30)) * T**9


def opacity(rho, T, comp, hminus=True):
    k = kappa_kramers(rho, T, comp) + kappa_es(comp)
    if hminus:
        kh = kappa_hminus(rho, T, comp)
        k = 1.0 / (1.0 / kh + 1.0 / k)
    return k


# ---------------------------------------------------------------- nuclear
def eps_pp(rho, T, comp):
    T9 = T * 1e-9
    g11 = 1 + 3.82 * T9 + 1.51 * T9**2 + 0.144 * T9**3 - 0.0114 * T9**4
    return 2.57e4 * g11 * rho * comp.X**2 * T9 ** (-2 / 3) * np.exp(-3.381 * T9 ** (-1 / 3))


def eps_cno(rho, T, comp):
    T9 = T * 1e-9
    g14 = 1 - 2.00 * T9 + 3.41 * T9**2 - 2.43 * T9**3
    return (8.24e25 * g14 * comp.X_cno * comp.X * rho * T9 ** (-2 / 3)
            * np.exp(-15.231 * T9 ** (-1 / 3) - (T9 / 0.8) ** 2))


def eps_nuc(rho, T, comp):
    return eps_pp(rho, T, comp) + eps_cno(rho, T, comp)


def temperature_exponent(f, T, rho=1.0, comp=Composition(), dlnT=1e-4):
    """Local nu = dln eps / dln T at fixed rho (numerical)."""
    up = f(rho, T * np.exp(dlnT), comp)
    dn = f(rho, T * np.exp(-dlnT), comp)
    return (np.log(up) - np.log(dn)) / (2 * dlnT)


# ---------------------------------------------------------------- Gamow peak
def gamow_peak(T, Z1, Z2, A1, A2):
    """Gamow peak energy E0 (erg) and 1/e width Delta (erg), from first principles.

    E0 = (b kT / 2)^(2/3), b = pi alpha Z1 Z2 sqrt(2 m_r c^2) = sqrt(E_G);  Delta = 4 sqrt(E0 kT / 3).
    """
    m_r = A1 * A2 / (A1 + A2) * m_u
    b = np.pi * alpha_fs * Z1 * Z2 * np.sqrt(2 * m_r * c**2)
    kT = k_B * T
    E0 = (b * kT / 2) ** (2 / 3)
    width = 4 * np.sqrt(E0 * kT / 3)
    return E0, width
