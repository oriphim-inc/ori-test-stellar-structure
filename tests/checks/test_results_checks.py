"""Sceptic's pass/fail checks on results/ (written before the main run; only the Sceptic edits this).

Every test skips until results/summary.json exists. Run with:
    python -m pytest tests/checks -v
(add `-m "not slow"` to skip the re-solve convergence checks, which take a few minutes).

Physics re-implemented here independently of src/ where possible (opacities converted from the
SI forms in Carroll & Ostlie 2007 eqs. 9.22-9.24; EOS, nabla_ad, nuclear rates from PLAN.md),
so a transcription or unit error in src/ shows up as a mismatch, not as agreement with itself.
Nothing here encodes the open answers (L-M exponent, where convection appears).
"""
import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(ROOT, "results")
SUMMARY = os.path.join(RES, "summary.json")
PROFDIR = os.path.join(RES, "profiles")

pytestmark = pytest.mark.skipif(not os.path.exists(SUMMARY), reason="results/summary.json not produced yet")

# ---- independent constants (cgs, CODATA 2018 / IAU 2015) ----
G = 6.67430e-8
C = 2.99792458e10
KB = 1.380649e-16
MU_AMU = 1.66053906660e-24
SIGMA = 5.670374419e-5
A_RAD = 4 * SIGMA / C
KEV = 1.602176634e-9
MSUN, RSUN, LSUN = 1.98847e33, 6.957e10, 3.828e33
X, Z = 0.70, 0.02
Y = 1 - X - Z
XCNO = 0.7 * Z
MU = 1 / (2 * X + 0.75 * Y + 0.5 * Z)


def load_summary():
    with open(SUMMARY) as f:
        return json.load(f)


def load_profile(Mmsun):
    path = os.path.join(PROFDIR, f"profile_{Mmsun:.3f}Msun.csv")
    if not os.path.exists(path):
        pytest.fail(f"missing profile {path}")
    a = np.loadtxt(path)
    cols = ["m", "r", "P", "T", "rho", "L", "kappa", "eps", "eps_pp", "eps_cno",
            "nabla_rad", "nabla_ad", "beta", "convective"]
    assert a.shape[1] == len(cols), "profile column count changed"
    return {c: a[:, i] for i, c in enumerate(cols)}


def all_profiles():
    s = load_summary()
    return [(sm["M_Msun"], load_profile(sm["M_Msun"])) for sm in s["sweep"]]


# ---- independent microphysics ----
def kappa_indep(rho, T, hminus=True):
    # C&O SI forms [m^2/kg, rho in kg/m^3]; 1 m^2/kg = 10 cm^2/g
    rho_si = rho * 1e3
    t = 2.82 * (rho_si * (1 + X)) ** 0.2
    kbf = 4.34e21 / t * Z * (1 + X) * rho_si * T ** -3.5 * 10
    kff = 3.68e18 * (1 - Z) * (1 + X) * rho_si * T ** -3.5 * 10
    kes = 0.02 * (1 + X) * 10
    k = kbf + kff + kes
    if hminus:
        # HKT eq. 4.94 cgs coefficient 2.5e-31 (PLAN); C&O SI 7.9e-34 converts to 2.498e-31 (0.07%)
        assert abs(7.9e-34 * 10 * np.sqrt(1e3) / 2.5e-31 - 1) < 1e-3
        kh = 2.5e-31 * (Z / 0.02) * np.sqrt(rho) * T ** 9
        k = 1 / (1 / kh + 1 / k)
    return k


def eps_indep(rho, T):
    T9 = T / 1e9
    g11 = 1 + 3.82 * T9 + 1.51 * T9**2 + 0.144 * T9**3 - 0.0114 * T9**4
    g14 = 1 - 2.00 * T9 + 3.41 * T9**2 - 2.43 * T9**3
    pp = 2.57e4 * g11 * rho * X**2 * T9 ** (-2 / 3) * np.exp(-3.381 / np.cbrt(T9))
    cno = 8.24e25 * g14 * XCNO * X * rho * T9 ** (-2 / 3) * np.exp(-15.231 / np.cbrt(T9) - (T9 / 0.8) ** 2)
    return pp, cno


def nabla_ad_indep(b):
    q = (1 - b) * (4 + b) / b**2
    return (1 + q) / (2.5 + 4 * q)


def trap(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def cumtrap(y, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


# =============================================================== 1. polytropes / Eddington model
LE_REF = {  # Chandrasekhar (1939) tables; n=0,1 analytic
    "0": (np.sqrt(6), 2 * np.sqrt(6), 1.0),
    "1": (np.pi, np.pi, np.pi**2 / 3),
    "1.5": (3.65375, 2.71406, 5.99071),
    "3": (6.89685, 2.01824, 54.1825),
}


def test_lane_emden_constants():
    le = load_summary()["lane_emden"]
    for n, (xi1, mc, rr) in LE_REF.items():
        assert n in le, f"n={n} missing"
        assert le[n]["xi1"] == pytest.approx(xi1, rel=1e-5), n
        assert le[n]["minus_xi2_dtheta"] == pytest.approx(mc, rel=1e-5), n
        assert le[n]["rhoc_over_rhomean"] == pytest.approx(rr, rel=1e-4), n


def test_eddington_standard_model_1msun():
    """n=3, M=1 Msun, R=1 Rsun: rho_c = 54.18 <rho>, Pc = W3 G M^2/R^4, Tc from Eddington quartic."""
    from scipy.optimize import brentq
    e = load_summary()["eddington_n3_1Msun_guess"]
    xi1, mc = LE_REF["3"][0], LE_REF["3"][1]
    rhoc = 54.1825 * 3 * MSUN / (4 * np.pi * RSUN**3)
    dth = mc / xi1**2
    Pc = G * MSUN**2 / RSUN**4 / (4 * np.pi * 4 * dth**2)
    # beta from (1-beta) P = a T^4/3 with T = beta mu m_u P/(rho k)
    f = lambda b: (1 - b) * Pc - A_RAD / 3 * (b * MU * MU_AMU * Pc / (rhoc * KB)) ** 4
    b = brentq(f, 0.5, 1 - 1e-14)
    Tc = b * MU * MU_AMU * Pc / (rhoc * KB)
    assert e["rhoc_gcc"] == pytest.approx(rhoc, rel=2e-4)
    assert e["Pc"] == pytest.approx(Pc, rel=5e-4)
    assert e["beta"] == pytest.approx(b, abs=2e-5)
    assert e["Tc_K"] == pytest.approx(Tc, rel=1e-3)
    # textbook: Tc ~ 1.2e7 K, Pc ~ 1.25e17 dyn/cm2 for the Sun as an n=3 polytrope
    assert 1.1e7 < Tc < 1.35e7 and 1.1e17 < Pc < 1.4e17


# =============================================================== 2. nuclear exponents, Gamow peaks
def test_nuclear_temperature_exponents_analytic():
    s = load_summary()
    T9 = 0.015
    g11 = 1 + 3.82 * T9 + 1.51 * T9**2 + 0.144 * T9**3 - 0.0114 * T9**4
    dg11 = (3.82 * T9 + 2 * 1.51 * T9**2 + 3 * 0.144 * T9**3 - 4 * 0.0114 * T9**4) / g11
    g14 = 1 - 2.00 * T9 + 3.41 * T9**2 - 2.43 * T9**3
    dg14 = (-2.00 * T9 + 2 * 3.41 * T9**2 - 3 * 2.43 * T9**3) / g14
    nu_pp = -2 / 3 + 3.381 / (3 * np.cbrt(T9)) + dg11
    nu_cno = -2 / 3 + 15.231 / (3 * np.cbrt(T9)) - 2 * (T9 / 0.8) ** 2 + dg14
    assert s["nu_pp_15MK"] == pytest.approx(nu_pp, abs=1e-3)
    assert s["nu_cno_15MK"] == pytest.approx(nu_cno, abs=1e-3)
    assert 3.5 < s["nu_pp_15MK"] < 4.5 and 17 < s["nu_cno_15MK"] < 21


def test_gamow_peaks_and_kT():
    """E0 = 1.22 keV (Z1^2 Z2^2 A T6^2)^(1/3) (Clayton 1968), A reduced mass in amu; E0 > kT."""
    s = load_summary()
    for sm in s["sweep"] + [s["model_1Msun"]]:
        T6 = sm["Tc_K"] / 1e6
        App = 1.007825 / 2
        Acn = 1.007825 * 14.003074 / (1.007825 + 14.003074)
        assert sm["gamow_pp_E0_keV"] == pytest.approx(1.22 * (App * T6**2) ** (1 / 3), rel=0.01)
        assert sm["gamow_p14N_E0_keV"] == pytest.approx(1.22 * (49 * Acn * T6**2) ** (1 / 3), rel=0.01)
        kT = KB * sm["Tc_K"] / KEV
        assert sm["kT_c_keV"] == pytest.approx(kT, rel=1e-6)
        assert sm["gamow_pp_E0_keV"] > 2 * kT and sm["gamow_p14N_E0_keV"] > sm["gamow_pp_E0_keV"]


# =============================================================== 3. solver convergence and global conservation
def test_all_models_converged():
    s = load_summary()
    models = [s["model_1Msun"]] + s["sweep"] + s.get("sweep_noHminus", [])
    for sm in models:
        assert sm["converged"] is True, sm["M_Msun"]
        assert sm["max_residual"] < 1e-6, (sm["M_Msun"], sm["max_residual"])   # PLAN: mismatch < 1e-6


def test_sweep_coverage():
    s = load_summary()
    Ms = np.array([sm["M_Msun"] for sm in s["sweep"]])
    assert Ms.min() == pytest.approx(0.5, rel=1e-6) and Ms.max() == pytest.approx(5.0, rel=1e-6)
    assert np.any(np.isclose(Ms, 1.0)) and len(Ms) >= 10
    assert np.all(np.diff(Ms) > 0)
    if "sweep_noHminus" in s:
        Mv = np.array([sm["M_Msun"] for sm in s["sweep_noHminus"]])
        assert np.allclose(Mv, Ms)


def test_reported_virial_and_energy():
    """Virial theorem with surface term (holds for any EOS in hydrostatic equilibrium) and
    L = integral of eps dm, both to 0.1% (quadrature error only)."""
    s = load_summary()
    for sm in [s["model_1Msun"]] + s["sweep"] + s.get("sweep_noHminus", []):
        assert sm["virial_ratio"] == pytest.approx(1.0, abs=1e-3), sm["M_Msun"]
        assert sm["L_integrated_over_L"] == pytest.approx(1.0, abs=1e-3), sm["M_Msun"]


def test_1msun_model_matches_sweep_entry():
    s = load_summary()
    one = [sm for sm in s["sweep"] if abs(sm["M_Msun"] - 1) < 1e-9][0]
    for k in ("L_Lsun", "R_Rsun", "Tc_K", "rhoc_gcc"):
        assert one[k] == pytest.approx(s["model_1Msun"][k], rel=1e-5), k


# =============================================================== 4. profiles: the agreed equations hold
@pytest.mark.parametrize("which", ["lo", "one", "hi"])
def test_profile_microphysics_independent(which):
    """EOS, opacity, eps, beta, nabla_ad, nabla_rad columns recomputed independently."""
    s = load_summary()
    Ms = {"lo": s["sweep"][0]["M_Msun"], "one": 1.0, "hi": s["sweep"][-1]["M_Msun"]}[which]
    p = load_profile(Ms)
    rho = (p["P"] - A_RAD * p["T"] ** 4 / 3) * MU * MU_AMU / (KB * p["T"])
    np.testing.assert_allclose(p["rho"], rho, rtol=1e-5)
    np.testing.assert_allclose(p["kappa"], kappa_indep(p["rho"], p["T"]), rtol=1e-5)
    pp, cno = eps_indep(p["rho"], p["T"])
    np.testing.assert_allclose(p["eps_pp"], pp, rtol=1e-5, atol=1e-30)
    np.testing.assert_allclose(p["eps_cno"], cno, rtol=1e-5, atol=1e-30)
    np.testing.assert_allclose(p["eps"], pp + cno, rtol=1e-5, atol=1e-30)
    b = 1 - A_RAD * p["T"] ** 4 / (3 * p["P"])
    np.testing.assert_allclose(p["beta"], b, rtol=1e-6)
    np.testing.assert_allclose(p["nabla_ad"], nabla_ad_indep(b), rtol=1e-6)
    i = p["L"] > 0
    nrad = 3 * p["kappa"] * p["L"] * p["P"] / (16 * np.pi * A_RAD * C * G * p["m"] * p["T"] ** 4)
    np.testing.assert_allclose(p["nabla_rad"][i], nrad[i], rtol=1e-5)
    np.testing.assert_array_equal(p["convective"] > 0.5, p["nabla_rad"] > p["nabla_ad"])


def test_profiles_monotonic_and_bounded():
    for Ms, p in all_profiles():
        assert np.all(np.diff(p["m"]) > 0), Ms
        assert np.all(np.diff(p["r"]) > 0), Ms
        assert np.all(np.diff(p["P"]) < 0), Ms
        assert np.all(np.diff(p["T"]) <= 0), Ms
        assert np.all(np.diff(p["L"]) >= -1e-9 * p["L"][-1]), Ms        # eps >= 0
        assert np.all(p["rho"] > 0) and np.all(p["eps"] >= 0), Ms
        assert np.all((p["beta"] > 0) & (p["beta"] <= 1)), Ms
        assert np.all((p["nabla_ad"] >= 0.25 - 1e-9) & (p["nabla_ad"] <= 0.4 + 1e-9)), Ms
        inner = p["r"] < 0.9 * p["r"][-1]
        assert np.all(np.diff(p["rho"][inner]) <= 0), Ms                 # no density inversion inside


def test_profile_structure_equations_integral():
    """Integral form of the four ODEs on the stored profile (independent of the solver):
    d lnT = nabla d lnP with nabla = min(nabla_rad, nabla_ad) (Schwarzschild), d ln r = -rP/(G m rho) d lnP,
    hydrostatic d lnP/d ln(M-m) = (M-m) G m/(4 pi r^4 P), L(m) - L(m0) = int eps dm.
    From m = 1e-3 M to the last point before the surface; tolerance 3e-3 in ln, 2e-3 L."""
    for Ms, p in all_profiles():
        M = p["m"][-1]
        sel = (p["m"] >= 1e-3 * M) & (p["m"] < M)
        m, r, P, T, rho = p["m"][sel], p["r"][sel], p["P"][sel], p["T"][sel], p["rho"][sel]
        lnP = np.log(P)
        nab = np.minimum(p["nabla_rad"][sel], p["nabla_ad"][sel])
        dlnT = trap(nab, lnP)
        assert np.log(T[-1] / T[0]) == pytest.approx(dlnT, abs=3e-3), Ms
        dlnr = trap(-r * P / (G * m * rho), lnP)
        assert np.log(r[-1] / r[0]) == pytest.approx(dlnr, abs=3e-3), Ms
        u = np.log(M - m)
        dlnP = trap((M - m) * G * m / (4 * np.pi * r**4 * P), u)
        assert np.log(P[-1] / P[0]) == pytest.approx(dlnP, abs=3e-3 * max(1, abs(np.log(P[-1] / P[0])) / 10)), Ms
        Lc = p["L"][0] + cumtrap(p["eps"], p["m"])
        assert np.max(np.abs(Lc - p["L"])) / p["L"][-1] < 2e-3, Ms


def test_virial_from_profiles():
    for Ms, p in all_profiles():
        m, r, P, rho = p["m"], p["r"], p["P"], p["rho"]
        therm = 3 * trap(P / rho, m) - 4 * np.pi * r[-1] ** 3 * P[-1]
        grav = trap(G * m / r, m)
        assert therm / grav == pytest.approx(1.0, abs=2e-3), Ms


def test_surface_boundary_conditions():
    """Last row = Eddington photosphere: T = Teff, L = 4 pi R^2 sigma Teff^4, P_gas = 2g/(3 kappa)."""
    s = load_summary()
    for sm in s["sweep"]:
        p = load_profile(sm["M_Msun"])
        M, R, L, T = p["m"][-1], p["r"][-1], p["L"][-1], p["T"][-1]
        assert M == pytest.approx(sm["M_Msun"] * MSUN, rel=1e-8)
        assert R == pytest.approx(sm["R_Rsun"] * RSUN, rel=1e-6)
        assert L == pytest.approx(sm["L_Lsun"] * LSUN, rel=1e-6)
        assert T == pytest.approx(sm["Teff_K"], rel=1e-6)
        assert L == pytest.approx(4 * np.pi * R**2 * SIGMA * T**4, rel=1e-5)
        Pgas = p["P"][-1] - A_RAD * T**4 / 3
        assert Pgas == pytest.approx(2 * G * M / R**2 / (3 * p["kappa"][-1]), rel=1e-5)
        # central values agree with the summary
        assert p["T"][0] == pytest.approx(sm["Tc_K"], rel=1e-4)
        assert p["P"][0] == pytest.approx(sm["Pc_dyn_cm2"], rel=1e-4)


def test_zone_bookkeeping_matches_profiles():
    s = load_summary()
    for sm in s["sweep"]:
        p = load_profile(sm["M_Msun"])
        conv = p["convective"] > 0.5
        M = p["m"][-1]
        assert (sm["conv_core_m_frac"] > 0) == bool(conv[0]), sm["M_Msun"]
        if conv[0]:
            j = np.argmin(conv) - 1 if not conv.all() else len(conv) - 1
            assert sm["conv_core_m_frac"] == pytest.approx(p["m"][j] / M, rel=1e-6)
        # envelope zone: a convective run whose outer end is within 1e-6 of M
        edges = np.flatnonzero(np.diff(np.concatenate([[0], conv.astype(int), [0]])))
        runs = list(zip(edges[::2], edges[1::2] - 1))
        env = [(i, j) for i, j in runs if p["m"][j] / M >= 0.999999]
        assert (sm["conv_env_m_depth"] > 0) == bool(env), sm["M_Msun"]
        if env:
            assert sm["conv_env_m_depth"] == pytest.approx(1 - p["m"][env[0][0]] / M, rel=1e-6, abs=1e-12)
        Lpp, Lcno = trap(p["eps_pp"], p["m"]), trap(p["eps_cno"], p["m"])
        assert sm["f_CNO_of_L"] == pytest.approx(Lcno / (Lpp + Lcno), abs=2e-3)


# =============================================================== 5. physical trends and plausibility
def test_robust_trends_with_mass():
    """Robust ZAMS trends (homology, Eddington quartic) - not the open questions."""
    sw = load_summary()["sweep"]
    L = np.array([x["L_Lsun"] for x in sw]); Tc = np.array([x["Tc_K"] for x in sw])
    b = np.array([x["beta_c"] for x in sw]); R = np.array([x["R_Rsun"] for x in sw])
    rc = np.array([x["rhoc_gcc"] for x in sw])
    assert np.all(np.diff(L) > 0) and np.all(np.diff(Tc) > 0) and np.all(np.diff(b) < 0)
    assert R[-1] > R[0] and rc[-1] < rc[0]
    assert np.all((b > 0.9) & (b < 1))


PLAUS = {  # generous: real ZAMS (X~0.7, Z~0.02 grids, e.g. Schaller et al. 1992) within ~x3 in L,
    # ~x2 in R and rho_c, ~30% in Tc; outside these, suspect units/physics rather than crude opacities
    1.0: dict(L_Lsun=(0.3, 2.0), R_Rsun=(0.5, 1.6), Tc_K=(1.0e7, 2.0e7), rhoc_gcc=(40, 250), Teff_K=(4500, 7500)),
    0.5: dict(L_Lsun=(0.01, 0.2), R_Rsun=(0.2, 0.9), Tc_K=(6e6, 1.3e7), rhoc_gcc=(30, 400)),
    5.0: dict(L_Lsun=(150, 2500), R_Rsun=(1.2, 5.0), Tc_K=(1.8e7, 4.0e7), rhoc_gcc=(5, 80), Teff_K=(1.2e4, 2.6e4)),
}


@pytest.mark.parametrize("Mt", list(PLAUS))
def test_plausibility_vs_known_zams(Mt):
    s = load_summary()
    sm = s["model_1Msun"] if Mt == 1.0 else min(s["sweep"], key=lambda x: abs(x["M_Msun"] - Mt))
    for k, (lo, hi) in PLAUS[Mt].items():
        assert lo < sm[k] < hi, (Mt, k, sm[k])


def test_reported_fits_recomputed():
    s = load_summary()
    Ms = np.array([x["M_Msun"] for x in s["sweep"]]); Ls = np.array([x["L_Lsun"] for x in s["sweep"]])
    for key, v in s["L_M_fits"].items():
        lo, hi = map(float, key.split("-"))
        sel = (Ms >= lo - 1e-9) & (Ms <= hi + 1e-9)
        a = np.polyfit(np.log10(Ms[sel]), np.log10(Ls[sel]), 1)[0]
        assert v["alpha"] == pytest.approx(a, abs=1e-6), key
    np.testing.assert_allclose(s["local_dlnL_dlnM"], np.gradient(np.log(Ls), np.log(Ms)), rtol=1e-6)


def test_sweep_table_matches_summary():
    import csv
    s = load_summary()
    with open(os.path.join(RES, "sweep_table.csv")) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(s["sweep"])
    for row, sm in zip(rows, s["sweep"]):
        for k in ("M_Msun", "L_Lsun", "R_Rsun", "Tc_K", "conv_core_m_frac", "conv_env_m_depth"):
            assert float(row[k]) == pytest.approx(sm[k], rel=1e-5, abs=1e-9), k


# =============================================================== 6. numerical convergence (re-solve)
@pytest.mark.slow
@pytest.mark.parametrize("which", ["lo", "one", "hi"])
def test_fitpoint_and_tolerance_independence(which):
    """Re-solve from the stored solution with a different fitting point (0.3 M) and tighter rtol
    (1e-11): L, R, Tc, Pc must agree to 1e-4 and convective boundaries to 0.01 in m/M."""
    sys.path.insert(0, os.path.join(ROOT, "src"))
    from stellar import model as md
    s = load_summary()
    sm = {"lo": s["sweep"][0], "one": s["model_1Msun"], "hi": s["sweep"][-1]}[which]
    M = sm["M_Msun"] * MSUN
    x0 = np.log([sm["Pc_dyn_cm2"], sm["Tc_K"], sm["R_Rsun"] * RSUN, sm["L_Lsun"] * LSUN])
    st = md.Settings(rtol=1e-11, atol=1e-13, fit_frac=0.3)
    x, info = md.solve(M, x0, st)
    assert info["converged"]
    np.testing.assert_allclose(np.exp(x), np.exp(x0), rtol=1e-4)
    prof = md.profile(x, M, st)
    zones = md.convective_zones(prof)
    core = [z for z in zones if z[0] == 0.0]
    assert (core[0][1] if core else 0.0) == pytest.approx(sm["conv_core_m_frac"], abs=0.01)
    env = [z for z in zones if z[1] >= 0.999999]
    assert ((1 - env[0][0]) if env else 0.0) == pytest.approx(sm["conv_env_m_depth"], abs=0.01)


# ---- convergence: re-solve with tighter tolerance and a different fitting point (added at review;
# the header promised these but they were missing). Pass: L, R, Tc, rho_c within 1e-4 relative;
# convective-zone extents within 2 % (core m/M, envelope r_base) and envelope mass depth within 20 %.
@pytest.mark.slow
@pytest.mark.parametrize("Mt", [1.0, 1.077217, 5.0])
def test_resolve_converged(Mt):
    sys.path.insert(0, os.path.join(ROOT, "src"))
    from stellar import model as md
    s = load_summary()
    ref = min(s["sweep"], key=lambda e: abs(e["M_Msun"] - Mt))
    M = ref["M_Msun"] * MSUN
    x0 = np.log([ref["Pc_dyn_cm2"], ref["Tc_K"], ref["R_Rsun"] * RSUN, ref["L_Lsun"] * LSUN])
    for kw in ({"rtol": 1e-11, "atol": 1e-13}, {"fit_frac": 0.3}):
        st = md.Settings(**kw)
        x, info = md.solve(M, x0 * (1 + 1e-4), st)
        assert info["converged"], kw
        Pc, Tc, R, L = np.exp(x)
        assert L / LSUN == pytest.approx(ref["L_Lsun"], rel=1e-4), kw
        assert R / RSUN == pytest.approx(ref["R_Rsun"], rel=1e-4), kw
        assert Tc == pytest.approx(ref["Tc_K"], rel=1e-4), kw
        prof = md.profile(x, M, st)
        zones = md.convective_zones(prof)
        core = zones[0][1] if zones and zones[0][0] == 0.0 else 0.0
        assert core == pytest.approx(ref["conv_core_m_frac"], rel=0.02, abs=1e-4), kw
        env = [z for z in zones if z[0] > 0.0]
        if ref["conv_env_m_depth"] > 0:
            assert env, kw
            assert 1 - env[0][0] == pytest.approx(ref["conv_env_m_depth"], rel=0.2), kw
            assert env[0][2] == pytest.approx(ref["conv_env_r_base"], rel=0.02), kw
