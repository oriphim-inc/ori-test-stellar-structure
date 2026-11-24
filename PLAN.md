# PLAN — stellar structure code, 1 Msun and a 0.5–5 Msun sweep

## Question
Build a stellar-structure (static, chemically homogeneous, zero-age main sequence) code solving the four
structure equations, validate it against the n = 3 polytrope, compute a 1 Msun model, then sweep
M = 0.5–5 Msun to get the mass–luminosity relation and the location of convective zones.

Deliverables: profiles of T, rho, P, L and convective regions vs radius; central conditions;
L–M scaling; Gamow peak energies.

## Physics (what is agreed)
Equations in mass coordinate m (cgs):
- dr/dm = 1/(4 pi r^2 rho)
- dP/dm = -G m / (4 pi r^4)
- dL/dm = eps_pp + eps_CNO
- dT/dm = -G m T nabla / (4 pi r^4 P), nabla = nabla_rad if nabla_rad <= nabla_ad, else nabla_ad
  (Schwarzschild criterion; convection taken as exactly adiabatic, no mixing-length superadiabaticity).
- nabla_rad = 3 kappa L P / (16 pi a c G m T^4).

EOS: P = rho k T / (mu m_u) + a T^4 / 3, fully ionized, mu = 1/(2X + 3Y/4 + Z/2).
nabla_ad for ideal gas + radiation with beta = P_gas/P:
nabla_ad = [1 + (1-beta)(4+beta)/beta^2] / [5/2 + 4(1-beta)(4+beta)/beta^2]
(Kippenhahn, Weigert & Weiss 2012, *Stellar Structure and Evolution*, ch. 13).

Composition (ASSUMPTION, solar-like ZAMS): X = 0.70, Y = 0.28, Z = 0.02, X_CNO = 0.7 Z.

Opacity (cgs), Kramers-type plus the two terms needed for the limits of the mass range:
- bound-free Kramers: kappa_bf = 4.34e25 (1/t) Z (1+X) rho T^-3.5, guillotine t = 2.82 [rho(1+X)]^0.2
  with rho in kg m^-3 (Carroll & Ostlie 2007, *Intro. to Modern Astrophysics*, eqs. 9.22, and StatStar).
- free-free Kramers: kappa_ff = 3.68e22 (1-Z)(1+X) rho T^-3.5 (gaunt factor 1; Carroll & Ostlie eq. 9.23).
- electron scattering: kappa_es = 0.2 (1+X) (needed at high mass, where Kramers alone gives L ~ M^5.5).
- H^- : kappa_H- = 2.5e-31 (Z/0.02) rho^0.5 T^9 (Hansen, Kawaler & Trimble 2004, eq. 4.94;
  Carroll & Ostlie eq. 9.24). Combined as 1/kappa = 1/kappa_H- + 1/(kappa_bf + kappa_ff + kappa_es).
  REASON (an addition beyond "Kramers-type"): with pure Kramers + ionized ideal gas, the outer envelope
  approaches nabla -> 1/(n+1) = 0.235 < 0.4 and stays radiative, so no envelope convection can appear.
  A comparison sweep with pure Kramers + es (no H^-) is run to show this dependence.

Energy generation (no screening, psi = f = 1; KWW 2012 eqs. 18.63, 18.65):
- eps_pp  = 2.57e4 g11 rho X^2 T9^(-2/3) exp(-3.381 T9^(-1/3)),
  g11 = 1 + 3.82 T9 + 1.51 T9^2 + 0.144 T9^3 - 0.0114 T9^4
- eps_CNO = 8.24e25 g14 X_CNO X rho T9^(-2/3) exp(-15.231 T9^(-1/3) - (T9/0.8)^2),
  g14 = 1 - 2.00 T9 + 3.41 T9^2 - 2.43 T9^3
Their local exponents d ln eps/d ln T at 15 MK are computed and reported (expect ~4 and ~18-20).

Gamow peak: E0 = (pi alpha Z1 Z2 kT sqrt(m_r c^2 / 2))^(2/3), computed from constants at each model's
central temperature for p+p and p+14N; tested against E0 = 1.22 keV (Z1^2 Z2^2 A T6^2)^(1/3).

## Boundary conditions and method
- Centre: m -> 0 series expansions (KWW §11.3) for r, P, T, L from trial (Pc, Tc).
- Surface: grey Eddington photosphere: T = Teff with L = 4 pi R^2 sigma Teff^4,
  P_gas = (2/3) g / kappa (solved self-consistently for rho), P = P_gas + a T^4/3.
- Shooting to a fitting point (Schwarzschild method): outward from centre and inward from surface to
  m_f (default 0.5 M); Newton iteration on (ln Pc, ln Tc, ln R, ln L) to zero the mismatch in
  ln r, ln P, ln T, L. Finite-difference Jacobian, damped steps. ODEs integrated with scipy (LSODA).
- Polytrope: Lane–Emden solved numerically for n = 0, 1, 3 (n = 0, 1 checked analytically;
  n = 3: xi1 = 6.89685, -xi1^2 theta'(xi1) = 2.01824, rho_c/rho_mean = 54.18). The n = 3 Eddington
  standard model (Eddington quartic for beta) gives the initial guess of Pc, Tc; R and L from
  rough homology (R ~ M^0.8, L ~ M^3.5) as guesses only.
- Sweep: masses 0.5–5 Msun (default 19 log-spaced points plus 1.0), continuation from neighbouring
  converged models. Numerical settings (tolerances, fit point, number of masses) are tunable;
  run.py has --n-masses and --rtol.

## What counts as the answer
- Converged 1 Msun model (mismatch < 1e-6), with profiles and central conditions, compared with
  standard ZAMS values (order-of-magnitude agreement expected: Tc ~ 13-14 MK, rho_c ~ 80-90 g/cc,
  L ~ 0.7 Lsun, R ~ 0.9 Rsun for a real ZAMS Sun; this crude physics may differ by tens of %).
- L ∝ M^alpha: fitted alpha over the sweep and local slope; expected ~3-5 for this range.
- Convective core / envelope extents (in m/M and r/R) vs M; expected core convection to appear
  around 1.1-1.5 Msun when CNO takes over, envelope convection shrinking with M.
- Consistency checks: virial theorem, L_surface = ∫ eps dm, n=3 polytrope constants.

## Scope / limits agreed up front
Homogeneous ZAMS only (no evolution, no diffusion); no screening; no partial ionization in EOS;
adiabatic convection (no MLT); simple opacity fits. These limit quantitative accuracy, especially
at 0.5 Msun (cool, partially ionized envelope) and in the outer convective envelope.
