"""Static stellar structure by shooting to a fitting point (Schwarzschild method).

Unknowns x = (ln Pc, ln Tc, ln R, ln L). Outward integration from the centre and inward integration
from the photosphere meet at m_f; Newton iteration zeroes the mismatch in (ln r, ln P, ln T, L/L).
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import physics as ph


@dataclass
class Settings:
    comp: ph.Composition = field(default_factory=ph.Composition)
    hminus: bool = True
    fit_frac: float = 0.5           # fitting point m_f / M
    m0_frac: float = 1e-10          # starting mass of outward integration / M
    rtol: float = 1e-9
    atol: float = 1e-11
    method: str = "LSODA"
    newton_tol: float = 1e-7
    max_iter: int = 60
    fd_step: float = 1e-5


class IntegrationFailure(RuntimeError):
    pass


# ---------------------------------------------------------------- local quantities
def local(P, T, r, m, L, s):
    """Return rho, kappa, eps, nabla_rad, nabla_ad, nabla at one point."""
    comp = s.comp
    rho = ph.density(P, T, comp)
    if not np.isfinite(rho) or rho <= 0:
        raise IntegrationFailure("non-positive gas pressure")
    kap = ph.opacity(rho, T, comp, s.hminus)
    eps = ph.eps_nuc(rho, T, comp)
    nrad = 3 * kap * L * P / (16 * np.pi * ph.a_rad * ph.c * ph.G * m * T**4)
    nad = ph.nabla_ad(ph.beta(P, T))
    return rho, kap, eps, nrad, nad, min(nrad, nad)


def rhs(m, y, s, Lscale):
    lnr, lnP, lnT, l = y
    r, P, T, L = np.exp(lnr), np.exp(lnP), np.exp(lnT), l * Lscale
    rho, kap, eps, nrad, nad, nab = local(P, T, r, m, L, s)
    dlnr = 1.0 / (4 * np.pi * r**3 * rho)
    dlnP = -ph.G * m / (4 * np.pi * r**4 * P)
    dlnT = dlnP * nab
    dl = eps / Lscale
    return [dlnr, dlnP, dlnT, dl]


# ---------------------------------------------------------------- boundary conditions
def center_bc(Pc, Tc, M, s):
    """Series expansion at m0 (KWW 2012 §11.3)."""
    comp = s.comp
    m0 = s.m0_frac * M
    rhoc = ph.density(Pc, Tc, comp)
    if rhoc <= 0:
        raise IntegrationFailure("central P_gas <= 0")
    epsc = ph.eps_nuc(rhoc, Tc, comp)
    kapc = ph.opacity(rhoc, Tc, comp, s.hminus)
    r = (3 * m0 / (4 * np.pi * rhoc)) ** (1 / 3)
    P = Pc - 3 * ph.G / (8 * np.pi) * (4 * np.pi * rhoc / 3) ** (4 / 3) * m0 ** (2 / 3)
    L = epsc * m0
    nrad_c = 3 * kapc * epsc * Pc / (16 * np.pi * ph.a_rad * ph.c * ph.G * Tc**4)
    nad_c = ph.nabla_ad(ph.beta(Pc, Tc))
    if nrad_c > nad_c:   # convective centre
        lnT = np.log(Tc) - (np.pi / 6) ** (1 / 3) * ph.G * nad_c * rhoc ** (4 / 3) * m0 ** (2 / 3) / Pc
        T = np.exp(lnT)
    else:
        T4 = Tc**4 - 1 / (2 * ph.a_rad * ph.c) * (3 / (4 * np.pi)) ** (2 / 3) * kapc * epsc \
            * rhoc ** (4 / 3) * m0 ** (2 / 3)
        T = T4**0.25
    return m0, np.array([np.log(r), np.log(P), np.log(T), L])


def photosphere(R, L, M, s):
    """Grey Eddington photosphere: T = Teff, P_gas = 2 g / (3 kappa), P = P_gas + a T^4 / 3."""
    comp = s.comp
    Teff = (L / (4 * np.pi * R**2 * ph.sigma_SB)) ** 0.25
    g = ph.G * M / R**2

    def f(lnrho):
        rho = np.exp(lnrho)
        Pgas = rho * ph.k_B * Teff / (comp.mu * ph.m_u)
        return np.log(Pgas) - np.log(2 * g / (3 * ph.opacity(rho, Teff, comp, s.hminus)))
    lnrho = brentq(f, np.log(1e-20), np.log(1e2), xtol=1e-14)
    rho = np.exp(lnrho)
    P = ph.pressure(rho, Teff, comp)
    return Teff, P, rho


# ---------------------------------------------------------------- integration
def _integrate(m_start, y0, m_end, s, Lscale, dense=False, M=None):
    """Outward: independent variable m. Inward (M given): mass depth q = M - m, from 0."""
    if M is None:
        def f(m, y):
            return rhs(m, y, s, Lscale)
        span = (m_start, m_end)
    else:
        def f(q, y):
            return [-v for v in rhs(M - q, y, s, Lscale)]
        span = (0.0, M - m_end)
    try:
        sol = solve_ivp(f, span, y0, method=s.method, rtol=s.rtol,
                        atol=s.atol, dense_output=dense)
    except IntegrationFailure as e:
        raise IntegrationFailure(str(e))
    except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError) as e:
        raise IntegrationFailure(str(e))
    if sol.status != 0 or not np.all(np.isfinite(sol.y[:, -1])):
        raise IntegrationFailure(f"solve_ivp: {sol.message}")
    return sol


def shoot(x, M, s, dense=False):
    lnPc, lnTc, lnR, lnL = x
    Pc, Tc, R, L = np.exp(x)
    mf = s.fit_frac * M
    m0, y0 = center_bc(Pc, Tc, M, s)
    y0 = y0.copy(); y0[3] /= L
    out = _integrate(m0, y0, mf, s, L, dense)
    Teff, Ps, rhos = photosphere(R, L, M, s)
    ys = np.array([lnR, np.log(Ps), np.log(Teff), 1.0])
    inn = _integrate(M, ys, mf, s, L, dense, M=M)
    resid = out.y[:, -1] - inn.y[:, -1]
    return resid, out, inn


def solve(M, x0, s=None, verbose=False):
    """Newton iteration on x = (ln Pc, ln Tc, ln R, ln L). Returns (x, info)."""
    s = s or Settings()
    x = np.array(x0, float)
    F, _, _ = shoot(x, M, s)
    hist = [float(np.max(np.abs(F)))]
    for it in range(s.max_iter):
        if np.max(np.abs(F)) < s.newton_tol:
            break
        J = np.empty((4, 4))
        for j in range(4):
            dx = np.zeros(4); dx[j] = s.fd_step
            Fp, _, _ = shoot(x + dx, M, s)
            J[:, j] = (Fp - F) / s.fd_step
        step = -np.linalg.solve(J, F)
        # damping: limit to 0.3 in any log, then backtrack
        lam = min(1.0, 0.3 / np.max(np.abs(step)))
        for _ in range(12):
            try:
                Fn, _, _ = shoot(x + lam * step, M, s)
                if np.max(np.abs(Fn)) < np.max(np.abs(F)) or lam < 1e-2:
                    break
            except (IntegrationFailure, ValueError):
                pass
            lam *= 0.5
        else:
            raise IntegrationFailure("line search failed")
        x = x + lam * step
        F = Fn
        hist.append(float(np.max(np.abs(F))))
        if verbose:
            print(f"   it {it:2d}  lam={lam:.3f}  |F|={hist[-1]:.2e}  "
                  f"Pc={np.exp(x[0]):.3e} Tc={np.exp(x[1]):.3e} R={np.exp(x[2])/ph.RSUN:.4f} "
                  f"L={np.exp(x[3])/ph.LSUN:.4f}")
    converged = np.max(np.abs(F)) < s.newton_tol
    return x, {"converged": bool(converged), "resid": F, "history": hist, "iterations": len(hist) - 1}


# ---------------------------------------------------------------- profiles
def _refine_boundaries(t, sol, to_state, s, L, hmin=0.0):
    """Insert points bracketing each root of nabla_rad - nabla_ad along coordinate t (sorted ascending).

    to_state(t, y) -> (r, P, T, m, Lr). Roots are found with brentq on the dense ODE solution; points
    are inserted at a relative offset 1e-7 of the local cell either side, so the grid resolves the
    boundary independent of the base spacing.
    """
    def d(tt):
        y = sol(tt)
        r, P, T, m, Lr = to_state(tt, y)
        _, _, _, nrad, nad, _ = local(P, T, r, m, max(Lr, 1e-300), s)
        return nrad - nad
    vals = np.array([d(tt) for tt in t])
    extra = []
    for i in np.flatnonzero(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0):
        root = brentq(d, t[i], t[i + 1], xtol=1e-14 * max(abs(t[i]), abs(t[i + 1]), 1e-300), rtol=1e-14)
        cell = t[i + 1] - t[i]
        h = min(max(1e-7 * cell, hmin), 0.25 * cell)
        extra += [root - h, root + h]
    return np.unique(np.concatenate([t, extra]))


def profile(x, M, s, npts_each=600):
    """Evaluate a converged model on a grid (centre -> surface), with convective boundaries resolved."""
    _, out, inn = shoot(x, M, s, dense=True)
    L = np.exp(x[3])
    m0 = out.t[0]
    mf = s.fit_frac * M
    # inner segment in m: log-spaced near centre, then linear
    m_in = np.concatenate([np.geomspace(m0, mf * 1e-2, npts_each // 3, endpoint=False),
                           np.linspace(mf * 1e-2, mf, 2 * npts_each // 3)])
    # outer segment in mass depth q = M - m (ascending q for the root search)
    dq = np.concatenate([[0.0], np.geomspace(1e-14 * M, M - mf, npts_each)[:-1]])
    m_in = _refine_boundaries(m_in, out.sol, lambda t, y: (np.exp(y[0]), np.exp(y[1]), np.exp(y[2]), t, y[3] * L), s, L)
    dq = _refine_boundaries(dq, inn.sol, lambda t, y: (np.exp(y[0]), np.exp(y[1]), np.exp(y[2]), M - t, y[3] * L), s, L,
                            hmin=8 * np.spacing(M))   # keep m = M - q distinct in float
    dq = dq[::-1]
    m_out = M - dq
    y_in = out.sol(m_in)
    y_out = inn.sol(dq)
    m = np.concatenate([m_in, m_out])
    y = np.concatenate([y_in, y_out], axis=1)
    keep = np.concatenate([[True], np.diff(m) > 0])     # guard against float-coincident m near M
    m, y = m[keep], y[:, keep]
    r, P, T, Lr = np.exp(y[0]), np.exp(y[1]), np.exp(y[2]), y[3] * L
    rho = np.empty_like(m); kap = np.empty_like(m); eps = np.empty_like(m)
    nrad = np.empty_like(m); nad = np.empty_like(m)
    for i in range(len(m)):
        rho[i], kap[i], eps[i], nrad[i], nad[i], _ = local(P[i], T[i], r[i], m[i], max(Lr[i], 1e-300), s)
    comp = s.comp
    epp = ph.eps_pp(rho, T, comp)
    ecno = ph.eps_cno(rho, T, comp)
    return {"m": m, "r": r, "P": P, "T": T, "L": Lr, "rho": rho, "kappa": kap, "eps": eps,
            "eps_pp": epp, "eps_cno": ecno, "nabla_rad": nrad, "nabla_ad": nad,
            "convective": nrad > nad, "beta": ph.beta(P, T)}


def convective_zones(prof):
    """List of (m_in/M, m_out/M, r_in/R, r_out/R) for each convective zone."""
    conv = prof["convective"]
    M, R = prof["m"][-1], prof["r"][-1]
    zones = []
    i = 0
    n = len(conv)
    while i < n:
        if conv[i]:
            j = i
            while j + 1 < n and conv[j + 1]:
                j += 1
            zones.append((prof["m"][i] / M if i > 0 else 0.0, prof["m"][j] / M,
                          prof["r"][i] / R if i > 0 else 0.0, prof["r"][j] / R))
            i = j + 1
        else:
            i += 1
    return zones


def initial_guess(M, s):
    """n=3 Eddington standard model for Pc, Tc; homology for R, L (rough guesses only)."""
    from .polytrope import polytrope_guess
    Ms = M / ph.MSUN
    R = ph.RSUN * Ms**0.8
    L = ph.LSUN * Ms**3.5
    g = polytrope_guess(M, R, s.comp, 3)
    return np.log([g["Pc"], g["Tc"], R, L]), g
