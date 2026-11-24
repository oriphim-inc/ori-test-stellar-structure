import numpy as np
import pytest

from stellar import physics as ph
from stellar import polytrope as pt
from stellar import model as md


def test_lane_emden_n0_n1_analytic():
    assert pt.lane_emden(0)["xi1"] == pytest.approx(np.sqrt(6), rel=1e-7)
    le1 = pt.lane_emden(1)
    assert le1["xi1"] == pytest.approx(np.pi, rel=1e-7)
    assert le1["mass_const"] == pytest.approx(np.pi, rel=1e-6)   # -xi^2 theta' = pi for n=1


def test_lane_emden_n3_constants():
    le = pt.lane_emden(3)
    assert le["xi1"] == pytest.approx(6.89685, abs=2e-5)
    assert le["mass_const"] == pytest.approx(2.01824, abs=2e-5)
    assert le["rhoc_over_rhomean"] == pytest.approx(54.18, abs=0.01)


def test_eddington_quartic_sun():
    # 1 - beta ~ 0.003 mu^4 beta^4 for 1 Msun
    mu = 0.61
    b = pt.eddington_beta(ph.MSUN, mu)
    assert (1 - b) / (mu**4 * b**4) == pytest.approx(0.003, rel=0.02)


def test_mu_and_eos_roundtrip():
    comp = ph.Composition()
    assert comp.mu == pytest.approx(1 / (2 * 0.7 + 0.75 * 0.28 + 0.01))
    rho, T = 50.0, 1.2e7
    P = ph.pressure(rho, T, comp)
    assert ph.density(P, T, comp) == pytest.approx(rho, rel=1e-12)


def test_nabla_ad_limits():
    assert ph.nabla_ad(1.0) == pytest.approx(0.4)
    assert ph.nabla_ad(1e-6) == pytest.approx(0.25, rel=1e-4)


def test_nuclear_temperature_exponents():
    nu_pp = ph.temperature_exponent(ph.eps_pp, 1.5e7)
    nu_cno = ph.temperature_exponent(ph.eps_cno, 1.5e7)
    assert 3.5 < nu_pp < 4.5
    assert 17 < nu_cno < 21


def test_gamow_peak_formula():
    # E0 = 1.22 keV (Z1^2 Z2^2 A T6^2)^(1/3), A = reduced mass in amu (Clayton 1968)
    for (Z1, Z2, A1, A2) in [(1, 1, 1.0, 1.0), (1, 7, 1.0, 14.0)]:
        A = A1 * A2 / (A1 + A2)
        for T6 in (10, 15, 25):
            E0, _ = ph.gamow_peak(T6 * 1e6, Z1, Z2, A1, A2)
            ref = 1.22 * (Z1**2 * Z2**2 * A * T6**2) ** (1 / 3)
            assert E0 / ph.keV == pytest.approx(ref, rel=0.01)


def test_hminus_combination_limits():
    comp = ph.Composition()
    # hot: H- negligible -> Kramers + es ; cool: H- dominates
    k_hot = ph.opacity(1.0, 1e6, comp)
    assert k_hot == pytest.approx(ph.kappa_kramers(1.0, 1e6, comp) + ph.kappa_es(comp), rel=1e-6)
    k_cool = ph.opacity(1e-7, 5000.0, comp)
    assert k_cool == pytest.approx(ph.kappa_hminus(1e-7, 5000.0, comp), rel=1e-3)


def test_1msun_converges_and_is_consistent():
    s = md.Settings(rtol=1e-8, atol=1e-10)
    x0, _ = md.initial_guess(ph.MSUN, s)
    x, info = md.solve(ph.MSUN, x0, s)
    assert info["converged"]
    prof = md.profile(x, ph.MSUN, s)
    # monotonic structure
    assert np.all(np.diff(prof["r"]) > 0)
    assert np.all(np.diff(prof["P"]) < 0)
    assert np.all(np.diff(prof["T"]) <= 0)
    # fit-point independence: different m_f gives the same star
    s2 = md.Settings(rtol=1e-8, atol=1e-10, fit_frac=0.3)
    x2, info2 = md.solve(ph.MSUN, x, s2)
    assert info2["converged"]
    assert np.allclose(x2, x, atol=1e-4)
