> **A model star from first principles.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)
>
> How this repo was made: Ori built and ran the project and wrote the paper. Afterwards, Oriphim added the licence, edited the author line and wrote this note, with help from Claude Code.

# Stellar structure code: 1 Msun and a 0.5–5 Msun main-sequence sweep

By Oriphim, with Ori.

## What was asked
Build a stellar-structure code that solves the four structure equations (mass continuity,
hydrostatic equilibrium, energy generation and energy transport). Physics: ideal gas + radiation EOS,
Kramers-type opacity, pp + CNO burning, and Schwarzschild convection with an adiabatic gradient.
Solve by shooting with centre and surface boundary conditions, starting from an n = 3 Lane–Emden
polytrope. Then sweep 0.5–5 Msun to get the mass–luminosity relation and track the convective zones.
Deliverables: profiles, central conditions, L–M scaling and Gamow peaks. The full agreed set-up is in
[PLAN.md](PLAN.md).

## What was found
All results are for a homogeneous zero-age main-sequence (ZAMS) star with X = 0.70, Z = 0.02.

**1 Msun model.** It converges in 6 Newton iterations (fit mismatch 3e-8):

| quantity | model | real ZAMS Sun (approx.) |
|---|---|---|
| L | 1.106 Lsun | ~0.7 Lsun |
| R | 0.944 Rsun | ~0.89 Rsun |
| Teff | 6091 K | ~5600–5700 K |
| Tc | 14.55 MK | ~13.5 MK |
| rho_c | 91.4 g cm^-3 | ~80–90 |
| Pc | 1.79e17 dyn cm^-2 | ~1.3e17 |
| beta_c (= P_gas/P) | 0.9994 | |
| CNO share of L | 0.8 % | ~1 % |
| convective core | none (nabla_rad/nabla_ad = 0.97 at the centre) | none |
| convective envelope | 1.2e-8 of the mass, base at r/R = 0.991 | ~0.02–0.03 of the mass, base at ~0.7 R |

The model is too luminous by about 60 % and slightly too hot. Its convective envelope is far too thin:
this is a consequence of the agreed simple physics (no ionization zones, analytic opacities), not a
bug. See [LIMITS.md](LIMITS.md), item 1.

**Gamow peak** at the 1 Msun central temperature (kT_c = 1.25 keV): p+p gives E0 = 5.79 keV (1/e
width 6.2 keV); p+14N, the slowest CNO step, gives E0 = 26.1 keV (width 13.2 keV). Over the sweep the
p+p peak runs from 4.0 keV (0.5 Msun) to 8.7 keV (5 Msun). The rates' local temperature exponents at
15 MK are nu_pp = 3.96 and nu_CNO = 19.9.

**Mass–luminosity relation.** The overall power law over 0.5–5 Msun is **L ∝ M^4.44**. The local slope
d ln L / d ln M falls steadily, from 5.4 at 0.5 Msun (close to the Kramers homology limit of 5.5) to
3.4 at 5 Msun, heading towards the electron-scattering limit of 3. Piecewise fits give 5.30
(0.5–1), 4.59 (1–2) and 3.74 (2–5). L does not change when H^- is removed (only R does), so L is set
by the interior. Radius goes as R ∝ M^0.48.

**Convective zones** (results/convective_zones.png):
- **Convective core** appears between 1.00 and 1.08 Msun, as CNO takes over: 0.25 % of the mass at
  1.08 Msun, 4 % at 1.22 Msun, 13 % at 1.58 Msun and 28 % at 5 Msun. The CNO share of L rises from 2 %
  to 98 % over the same range.
- **Convective envelope** is deepest at low mass and shrinks steeply with M. It holds 1.7e-3 of the
  mass at 0.5 Msun (base at r/R = 0.87), 1e-8 at 1 Msun, and disappears above about 1.6–1.8 Msun.
  Without H^- opacity there is no envelope convection at any mass. The trend is right qualitatively,
  but the depths are orders of magnitude too small (LIMITS.md, items 1 and 3).

**Polytrope sanity check.** The numerical Lane–Emden n = 3 solution gives xi1 = 6.89685,
-xi1^2 theta' = 2.01824 and rho_c/<rho> = 54.18. The n = 0 and n = 1 cases match their analytic
solutions. The Eddington standard model for 1 Msun at R = 1 Rsun gives Tc = 12.1 MK, rho_c = 76 g/cc
and beta = 0.99957; this was the initial guess. The converged model is more centrally concentrated in
T and slightly less in rho than the n = 3 polytrope with the same M and R
(results/profiles_1Msun.png).

## Evidence
- `results/summary.json`: all numbers above. `results/sweep_table.csv`: one row per mass (central
  conditions, Teff, R, L, zone extents, Gamow peaks, virial and energy checks).
- `results/profiles/profile_<M>Msun.csv`: m, r, P, T, rho, L, kappa, eps (pp and CNO), nabla_rad,
  nabla_ad, beta and a convective flag, for 20 masses.
- Figures: `profiles_1Msun.png` (T, rho, P, L, m and nabla vs r/R, with convection shaded and the n=3
  polytrope overlaid), `profiles_vs_mass.png`, `mass_luminosity.png`, `convective_zones.png`,
  `central_conditions.png`.
- Checks: 9 unit tests (`tests/test_physics.py`) and 28 pass/fail checks written independently by the
  Sceptic (`tests/checks/`). They cover the polytrope constants, the virial theorem (to 1e-8), energy
  conservation (L = ∫eps dm to 2e-7), the structure equations re-integrated from the stored profiles,
  the boundary conditions, and re-solves with 100x tighter tolerance and a different fitting point
  (L, R and Tc agree to 1e-4). All 37 pass.

## Method in one paragraph
The four equations are integrated in mass coordinate with scipy LSODA (rtol 1e-9). Integration runs
outward from the centre, using series expansions, and inward from a grey Eddington photosphere
(T = Teff, P_gas = 2g/3kappa). The envelope is integrated in mass depth M − m so that the surface is
resolved. Newton iteration on (ln Pc, ln Tc, ln R, ln L) zeroes the mismatch at m = 0.5 M. Opacity is
Kramers bound-free (guillotine) + free-free + electron scattering, combined harmonically with an H^-
fit. The H^- term is an addition to the brief, because with pure Kramers opacity no envelope can
become convective (PLAN.md). Burning rates are from Kippenhahn, Weigert & Weiss (2012), without
screening. The sweep uses continuation in mass. Convective boundaries are root-found where
nabla_rad = nabla_ad.

## Limits
See **[LIMITS.md](LIMITS.md)** (written by the Sceptic). In short: this is a homogeneous ZAMS model
with analytic opacities, no partial ionization and no screening. L, Tc and the L–M slopes are
trustworthy only to tens of %. Envelope convection depths, and R and Teff at 0.5 Msun, are not
quantitative. The onset of the convective core (1.0–1.08 Msun) is marginal at the few-% level.

## Re-run
```
uv venv .venv && uv pip install --python .venv/bin/python numpy scipy matplotlib pytest
.venv/bin/python run.py                 # ~15 s; options: --n-masses N, --rtol R, --skip-variant
.venv/bin/python -m pytest tests -q     # unit tests + Sceptic checks
```
