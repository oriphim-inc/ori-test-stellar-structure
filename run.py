"""Reproduce every result: python run.py [--n-masses N] [--rtol R] [--out results]

1. Lane-Emden polytropes (n = 0, 1, 3) and the n = 3 Eddington standard model for 1 Msun.
2. 1 Msun ZAMS model by shooting (profiles, central conditions, Gamow peaks, checks).
3. Mass sweep 0.5-5 Msun (baseline opacity with H^-) and comparison sweep without H^-.
4. Figures and results/summary.json.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mplcache"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import ScalarFormatter  # noqa: E402


def plain_log_x(ax):
    ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.ticklabel_format(axis="x", style="plain"); ax.tick_params(axis="x", which="minor", labelsize=7)
from scipy.integrate import simpson  # noqa: E402

from stellar import physics as ph  # noqa: E402
from stellar import model as md  # noqa: E402
from stellar import polytrope as pt  # noqa: E402


def write_profile_csv(path, prof):
    cols = ["m", "r", "P", "T", "rho", "L", "kappa", "eps", "eps_pp", "eps_cno",
            "nabla_rad", "nabla_ad", "beta", "convective"]
    arr = np.column_stack([prof[c].astype(float) for c in cols])
    hdr = ("m[g] r[cm] P[dyn/cm2] T[K] rho[g/cm3] L[erg/s] kappa[cm2/g] eps[erg/g/s] eps_pp eps_cno "
           "nabla_rad nabla_ad beta convective(1=yes)")
    np.savetxt(path, arr, header=hdr, fmt="%.17g")   # lossless: m near the surface differs from M by ~1e-14 M


def diagnostics(prof, x, M):
    """Global consistency checks on a converged profile."""
    m, r, P, rho, T = prof["m"], prof["r"], prof["P"], prof["rho"], prof["T"]
    L = np.exp(x[3])
    Lint = simpson(prof["eps"], x=m) + prof["L"][0]
    # virial: Omega = -int G m / r dm ; 3 int P/rho dm + Omega = 4 pi R^3 P_surf (surface term)
    Omega = -simpson(ph.G * m / r, x=m)
    therm = 3 * simpson(P / rho, x=m)
    surf = 4 * np.pi * r[-1] ** 3 * P[-1]
    return {"L_integrated_over_L": float(Lint / L),
            "virial_ratio": float((therm - surf) / (-Omega))}


def model_summary(M, x, info, prof, s):
    comp = s.comp
    Pc, Tc, R, L = np.exp(x)
    rhoc = ph.density(Pc, Tc, comp)
    Teff = (L / (4 * np.pi * R**2 * ph.sigma_SB)) ** 0.25
    zones = md.convective_zones(prof)
    core = [z for z in zones if z[0] == 0.0]
    env = [z for z in zones if z[1] >= 0.999999 or z[3] >= 0.999999]
    Lpp = simpson(prof["eps_pp"], x=prof["m"]); Lcno = simpson(prof["eps_cno"], x=prof["m"])
    E0pp, Wpp = ph.gamow_peak(Tc, 1, 1, 1.007825, 1.007825)
    E0cn, Wcn = ph.gamow_peak(Tc, 1, 7, 1.007825, 14.003074)
    d = {
        "M_Msun": M / ph.MSUN, "L_Lsun": L / ph.LSUN, "R_Rsun": R / ph.RSUN, "Teff_K": Teff,
        "Tc_K": Tc, "rhoc_gcc": rhoc, "Pc_dyn_cm2": Pc, "beta_c": ph.beta(Pc, Tc),
        "converged": info["converged"], "max_residual": float(np.max(np.abs(info["resid"]))),
        "newton_iterations": info["iterations"],
        "f_CNO_of_L": float(Lcno / (Lpp + Lcno)),
        "conv_core_m_frac": core[0][1] if core else 0.0,
        "conv_core_r_frac": core[0][3] if core else 0.0,
        "conv_env_m_depth": (1 - env[0][0]) if env else 0.0,   # mass fraction in envelope CZ
        "conv_env_r_base": env[0][2] if env else 1.0,          # r/R of base of envelope CZ
        "zones_m_r": [list(map(float, z)) for z in zones],
        "gamow_pp_E0_keV": E0pp / ph.keV, "gamow_pp_width_keV": Wpp / ph.keV,
        "gamow_p14N_E0_keV": E0cn / ph.keV, "gamow_p14N_width_keV": Wcn / ph.keV,
        "kT_c_keV": ph.k_B * Tc / ph.keV,
    }
    d.update(diagnostics(prof, x, M))
    return d


def homology_guess(prev_masses, prev_x, M):
    """Extrapolate ln x linearly in ln M from the two nearest converged models."""
    lm = np.log(M)
    if len(prev_masses) == 1:
        # rough homology slopes d ln(Pc, Tc, R, L) / d ln M (guess only)
        return prev_x[0] + np.array([-1.0, 0.2, 0.8, 3.5]) * (lm - np.log(prev_masses[0]))
    l0, l1 = np.log(prev_masses[-2]), np.log(prev_masses[-1])
    return prev_x[-1] + (prev_x[-1] - prev_x[-2]) * (lm - l1) / (l1 - l0)


def solve_continued(pm, px, M, s, depth=0):
    """Solve at M from an extrapolated guess; on failure solve an intermediate mass first."""
    try:
        x, info = md.solve(M, homology_guess(pm, px, M), s)
        if info["converged"]:
            return x, info
    except (md.IntegrationFailure, ValueError, np.linalg.LinAlgError):
        if depth >= 4:
            raise
    if depth >= 4:
        return x, info
    Mmid = np.sqrt(pm[-1] * M)
    xm, im = solve_continued(pm, px, Mmid, s, depth + 1)
    if not im["converged"]:
        return xm, im
    return solve_continued(pm + [Mmid], px + [xm], M, s, depth + 1)


def sweep(masses_msun, s, x_1msun, label, verbose=True):
    """Continuation from 1 Msun upward and downward."""
    masses = np.array(sorted(masses_msun))
    results = {}
    i1 = int(np.argmin(abs(masses - 1.0)))
    for branch in (masses[i1:], masses[: i1 + 1][::-1]):
        pm, px = [ph.MSUN], [x_1msun]
        for Ms in branch:
            M = Ms * ph.MSUN
            t0 = time.time()
            if abs(Ms - 1.0) < 1e-12:
                x, info = md.solve(M, x_1msun, s)
            else:
                x, info = solve_continued(pm, px, M, s)
            prof = md.profile(x, M, s)
            summ = model_summary(M, x, info, prof, s)
            results[round(Ms, 6)] = (x, info, prof, summ)
            if verbose:
                print(f"  [{label}] M={Ms:5.3f} conv={info['converged']} it={info['iterations']:2d} "
                      f"L={summ['L_Lsun']:9.4f} R={summ['R_Rsun']:6.3f} Tc={summ['Tc_K']/1e6:6.2f}MK "
                      f"core={summ['conv_core_m_frac']:.3f} env={summ['conv_env_m_depth']:.2e} "
                      f"({time.time()-t0:.1f}s)", flush=True)
            if not info["converged"]:
                raise RuntimeError(f"model M={Ms} did not converge")
            if abs(Ms - 1.0) > 1e-12:
                pm.append(M); px.append(x)
    return [results[k] for k in sorted(results)]


def fit_powerlaw(Ms, Ls, lo, hi):
    sel = (Ms >= lo - 1e-9) & (Ms <= hi + 1e-9)
    p = np.polyfit(np.log10(Ms[sel]), np.log10(Ls[sel]), 1)
    return float(p[0]), float(10 ** p[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-masses", type=int, default=19, help="log-spaced masses in [0.5, 5] (plus 1.0)")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--skip-variant", action="store_true", help="skip the no-H^- comparison sweep")
    args = ap.parse_args()
    out = args.out
    os.makedirs(os.path.join(out, "profiles"), exist_ok=True)
    t_start = time.time()
    summary = {"settings": {"n_masses": args.n_masses, "rtol": args.rtol}}

    # ---------------------------------------------------------------- 1. polytropes
    print("1. Lane-Emden polytropes")
    poly = {}
    for n in (0, 1, 1.5, 3):
        le = pt.lane_emden(n)
        poly[str(n)] = {"xi1": le["xi1"], "minus_xi2_dtheta": le["mass_const"],
                        "rhoc_over_rhomean": le["rhoc_over_rhomean"]}
        print(f"   n={n}: xi1={le['xi1']:.6f}  -xi^2 theta'={le['mass_const']:.6f}  "
              f"rho_c/<rho>={le['rhoc_over_rhomean']:.4f}")
    summary["lane_emden"] = poly

    s = md.Settings(rtol=args.rtol, atol=args.rtol * 1e-2)
    M1 = ph.MSUN
    x0, g = md.initial_guess(M1, s)
    edd = {"Pc": g["Pc"], "Tc_K": g["Tc"], "rhoc_gcc": g["rhoc"], "beta": g["beta"],
           "R_assumed_Rsun": 1.0}
    summary["eddington_n3_1Msun_guess"] = edd
    print(f"   Eddington n=3, 1 Msun, R=1 Rsun: Tc={g['Tc']:.3e} K rho_c={g['rhoc']:.1f} "
          f"Pc={g['Pc']:.3e} beta={g['beta']:.5f}")

    # ---------------------------------------------------------------- 2. 1 Msun model
    print("2. 1 Msun model")
    x1, info1 = md.solve(M1, x0, s, verbose=True)
    if not info1["converged"]:
        raise RuntimeError("1 Msun model failed to converge")
    prof1 = md.profile(x1, M1, s)
    sum1 = model_summary(M1, x1, info1, prof1, s)
    summary["model_1Msun"] = sum1
    write_profile_csv(os.path.join(out, "profiles", "profile_1.000Msun.csv"), prof1)
    for k in ("L_Lsun", "R_Rsun", "Teff_K", "Tc_K", "rhoc_gcc", "Pc_dyn_cm2", "f_CNO_of_L",
              "conv_core_m_frac", "conv_env_m_depth", "conv_env_r_base",
              "gamow_pp_E0_keV", "gamow_p14N_E0_keV", "L_integrated_over_L", "virial_ratio"):
        print(f"   {k:22s} {sum1[k]:.5g}")

    # polytrope vs model comparison (same M, R)
    polyprof = pt.polytrope_profile(M1, np.exp(x1[2]), s.comp)
    # nuclear exponents at 15 MK
    summary["nu_pp_15MK"] = float(ph.temperature_exponent(ph.eps_pp, 1.5e7))
    summary["nu_cno_15MK"] = float(ph.temperature_exponent(ph.eps_cno, 1.5e7))
    print(f"   nu_pp(15MK)={summary['nu_pp_15MK']:.2f}  nu_CNO(15MK)={summary['nu_cno_15MK']:.2f}")

    # ---------------------------------------------------------------- 3. sweeps
    masses = np.unique(np.round(np.append(np.geomspace(0.5, 5.0, args.n_masses), 1.0), 6))
    print(f"3. Sweep over {len(masses)} masses (baseline: Kramers + es + H^-)")
    base = sweep(masses, s, x1, "H-")
    variant = None
    if not args.skip_variant:
        print("   comparison sweep: Kramers + es, no H^-")
        sv = md.Settings(rtol=args.rtol, atol=args.rtol * 1e-2, hminus=False)
        x1v, iv = md.solve(M1, x1, sv)
        variant = sweep(masses, sv, x1v, "noH-")

    for (_, _, prof, sm) in base:
        write_profile_csv(os.path.join(out, "profiles", f"profile_{sm['M_Msun']:.3f}Msun.csv"), prof)
    tab = [sm for (_, _, _, sm) in base]
    keys = ["M_Msun", "L_Lsun", "R_Rsun", "Teff_K", "Tc_K", "rhoc_gcc", "Pc_dyn_cm2", "beta_c",
            "f_CNO_of_L", "conv_core_m_frac", "conv_core_r_frac", "conv_env_m_depth",
            "conv_env_r_base", "gamow_pp_E0_keV", "gamow_p14N_E0_keV", "kT_c_keV",
            "virial_ratio", "L_integrated_over_L", "max_residual"]
    with open(os.path.join(out, "sweep_table.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for sm in tab:
            f.write(",".join(f"{sm[k]:.6g}" for k in keys) + "\n")
    summary["sweep"] = tab
    if variant:
        summary["sweep_noHminus"] = [sm for (_, _, _, sm) in variant]

    Ms = np.array([t["M_Msun"] for t in tab]); Ls = np.array([t["L_Lsun"] for t in tab])
    Rs = np.array([t["R_Rsun"] for t in tab])
    fits = {}
    for lo, hi in ((0.5, 5.0), (0.5, 1.0), (1.0, 2.0), (2.0, 5.0)):
        a, c0 = fit_powerlaw(Ms, Ls, lo, hi)
        fits[f"{lo}-{hi}"] = {"alpha": a, "L_at_1Msun_fit": c0}
    aR, _ = fit_powerlaw(Ms, Rs, 0.5, 5.0)
    summary["L_M_fits"] = fits
    summary["R_M_exponent_0.5-5"] = aR
    summary["local_dlnL_dlnM"] = list(map(float, np.gradient(np.log(Ls), np.log(Ms))))
    if variant:
        Lv = np.array([t["L_Lsun"] for t in summary["sweep_noHminus"]])
        summary["L_M_fit_noHminus_0.5-5"] = fit_powerlaw(Ms, Lv, 0.5, 5.0)[0]
    print("   L-M power-law fits:", {k: round(v["alpha"], 3) for k, v in fits.items()})

    # ---------------------------------------------------------------- 4. figures
    make_figures(out, prof1, polyprof, base, variant, sum1)
    summary["runtime_s"] = time.time() - t_start
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"done in {summary['runtime_s']:.1f} s -> {out}")


def shade_conv(ax, prof, x):
    conv = prof["convective"]
    for z in md.convective_zones(prof):
        ax.axvspan(z[2] if x == "r" else z[0], z[3] if x == "r" else z[1], color="orange", alpha=0.25, lw=0)


def make_figures(out, prof1, polyprof, base, variant, sum1):
    R = prof1["r"][-1]; x = prof1["r"] / R
    fig, axs = plt.subplots(3, 2, figsize=(11, 11), sharex=True)
    items = [("T", "T [K]", True), ("rho", r"$\rho$ [g cm$^{-3}$]", True), ("P", "P [dyn cm$^{-2}$]", True),
             ("L", r"$L_r / L_\odot$", False), ("m", r"$m / M_\odot$", False), ("nabla", r"$\nabla$", False)]
    for ax, (k, lab, lg) in zip(axs.flat, items):
        shade_conv(ax, prof1, "r")
        if k == "nabla":
            ax.plot(x, np.minimum(prof1["nabla_rad"], 5), label=r"$\nabla_{rad}$")
            ax.plot(x, prof1["nabla_ad"], label=r"$\nabla_{ad}$")
            ax.set_ylim(0, 1.5); ax.legend()
        elif k in ("L", "m"):
            ax.plot(x, prof1[k] / (ph.LSUN if k == "L" else ph.MSUN), "k")
        else:
            ax.plot(x, prof1[k], "k", label="model")
            ax.plot(polyprof["r"] / R, polyprof[k], "b--", label="n=3 polytrope (same M, R)")
            ax.legend(fontsize=8)
        if lg:
            ax.set_yscale("log")
        ax.set_ylabel(lab)
    for ax in axs[-1]:
        ax.set_xlabel("r / R")
    fig.suptitle(f"1 Msun ZAMS model: L={sum1['L_Lsun']:.3f} Lsun, R={sum1['R_Rsun']:.3f} Rsun, "
                 f"Tc={sum1['Tc_K']/1e6:.2f} MK (orange = convective)")
    fig.tight_layout(); fig.savefig(os.path.join(out, "profiles_1Msun.png"), dpi=130); plt.close(fig)

    # L-M relation
    Ms = np.array([b[3]["M_Msun"] for b in base]); Ls = np.array([b[3]["L_Lsun"] for b in base])
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    axs[0].loglog(Ms, Ls, "ko-", label="model (Kramers+es+H$^-$)")
    if variant:
        axs[0].loglog(Ms, [v[3]["L_Lsun"] for v in variant], "s--", color="gray", label="no H$^-$")
    for a in (3, 4, 5):
        axs[0].loglog(Ms, Ms**a, ":", lw=0.8, label=f"$M^{a}$")
    axs[0].set_xlabel("M / Msun"); axs[0].set_ylabel("L / Lsun"); axs[0].legend(fontsize=8)
    axs[1].semilogx(Ms, np.gradient(np.log(Ls), np.log(Ms)), "ko-")
    axs[1].set_xlabel("M / Msun"); axs[1].set_ylabel("d ln L / d ln M")
    for ax in axs:
        plain_log_x(ax)
    fig.tight_layout(); fig.savefig(os.path.join(out, "mass_luminosity.png"), dpi=130); plt.close(fig)

    # Kippenhahn-like convective zone diagram
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    Mlist = [b[3]["M_Msun"] for b in base]
    lm = np.log10(Mlist)
    edges = 10 ** np.concatenate([[lm[0] - (lm[1] - lm[0]) / 2], (lm[1:] + lm[:-1]) / 2,
                                  [lm[-1] + (lm[-1] - lm[-2]) / 2]])
    for ax, coord, lab in ((axs[0], (0, 1), "m / M"), (axs[1], (2, 3), "r / R")):
        for i, b in enumerate(base):
            for z in md.convective_zones(b[2]):
                ax.fill_between([edges[i], edges[i + 1]], z[coord[0]], z[coord[1]], color="tab:blue",
                                alpha=0.6, lw=0)
        ax.set_xscale("log"); ax.set_ylim(0, 1); ax.set_xlim(edges[0], edges[-1])
        ax.set_xlabel("M / Msun"); ax.set_ylabel(lab)
        ax.set_title(f"convective zones (blue) in {lab}, baseline opacity", fontsize=10)
    for sw, st, lab in ((base, "o-", "Kramers+es+H$^-$ (baseline)"), (variant, "s--", "Kramers+es (no H$^-$)")):
        if sw is None:
            continue
        Mv = np.array([b[3]["M_Msun"] for b in sw])
        env = np.array([b[3]["conv_env_m_depth"] for b in sw])
        core = np.array([b[3]["conv_core_m_frac"] for b in sw])
        axs[2].semilogy(Mv[env > 0], env[env > 0], st, color="tab:red", label="envelope CZ mass, " + lab)
        axs[2].semilogy(Mv[core > 0], core[core > 0], st, color="tab:blue", label="core CZ mass, " + lab)
    axs[2].set_xscale("log"); axs[2].set_xlabel("M / Msun"); axs[2].set_ylabel("convective mass fraction")
    axs[2].legend(fontsize=7); axs[2].set_title("log scale (envelope zones are thin)", fontsize=10)
    for ax in axs:
        plain_log_x(ax)
    fig.tight_layout(); fig.savefig(os.path.join(out, "convective_zones.png"), dpi=130); plt.close(fig)

    # profiles for selected masses
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    for b in base:
        Mm = b[3]["M_Msun"]
        if not any(abs(Mm - t) < 0.05 * t for t in (0.5, 1.0, 2.0, 5.0)):
            continue
        p = b[2]; xr = p["r"] / p["r"][-1]
        axs[0, 0].semilogy(xr, p["T"], label=f"{Mm:.2f}")
        axs[0, 1].semilogy(xr, p["rho"], label=f"{Mm:.2f}")
        axs[1, 0].plot(xr, p["L"] / p["L"][-1], label=f"{Mm:.2f}")
        axs[1, 1].plot(xr, np.minimum(p["nabla_rad"], 2) - p["nabla_ad"], label=f"{Mm:.2f}")
    axs[1, 1].axhline(0, color="k", lw=0.5)
    for ax, lab in zip(axs.flat, ["T [K]", r"$\rho$ [g/cc]", r"$L_r/L$", r"$\nabla_{rad}-\nabla_{ad}$ (>0 convective)"]):
        ax.set_ylabel(lab); ax.set_xlabel("r / R"); ax.legend(title="M/Msun", fontsize=8)
    axs[1, 1].set_ylim(-0.4, 1.0)
    fig.tight_layout(); fig.savefig(os.path.join(out, "profiles_vs_mass.png"), dpi=130); plt.close(fig)

    # central conditions
    Tc = [b[3]["Tc_K"] for b in base]; rc = [b[3]["rhoc_gcc"] for b in base]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(rc, Tc, "ko-")
    for b in base:
        ax.annotate(f"{b[3]['M_Msun']:.2g}", (b[3]["rhoc_gcc"], b[3]["Tc_K"]), fontsize=7)
    ax.set_xlabel(r"$\rho_c$ [g cm$^{-3}$]"); ax.set_ylabel(r"$T_c$ [K]")
    fig.tight_layout(); fig.savefig(os.path.join(out, "central_conditions.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
