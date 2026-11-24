# LIMITS — what these models do NOT show

Written by the Sceptic after the main run (commit 51c8ff3), updated after the boundary-refinement
rerun (commit 0c2e342; 37/37 checks pass). The numerical solution is sound: every
model converges (mismatch < 1e-7), conserves energy (L = integral of eps dm to 2e-7) and satisfies the
virial theorem to 1e-8. Re-solving with 100x tighter ODE tolerance or a different fitting point
(0.3 M instead of 0.5 M) changes L, R and Tc by less than 1e-4 at 1.0, 1.08 and 5 Msun. The limits
below come from the **physics that was chosen** (PLAN.md), not from the numerics.

## Assumptions labelled in PLAN.md
- **Composition** X = 0.70, Y = 0.28, Z = 0.02, CNO = 0.7 Z. This is assumed, not fitted to the Sun.
  The results hold for this mixture only.
- **Zero-age and homogeneous.** There is no evolution, no diffusion and no 3He or CNO equilibrium
  history. The real present-day Sun (4.6 Gyr old, L = 1 Lsun, R = 1 Rsun) is not what was modelled.
- **No electron screening** (f = 1). Screening raises both rates by about 5-10 % at solar Tc, and
  raises CNO more than pp. Without it the model needs a slightly higher Tc.
- **Opacity is a set of analytic fits**: Kramers bound-free (guillotine) + free-free + electron
  scattering, combined harmonically with an H^- fit. There are no line opacities, no H/He bound-free
  opacity peak near 1e4-1e5 K, no molecules and no conductive opacity. H^- was **added beyond the
  brief** so that envelope convection can appear at all. The sweep without H^- shows that this
  choice alone decides whether any envelope convection appears.
- **EOS is fully ionized ideal gas + radiation.** There is no partial ionization (so nabla_ad = 0.4
  everywhere, even in the H/He ionization zones where it really falls to about 0.1), and no
  degeneracy or Coulomb corrections. The degeneracy and Coulomb corrections matter at the 1-2 % level
  for 0.5 Msun.
- **Convection is exactly adiabatic** (nabla = nabla_ad, with no mixing-length superadiabaticity).
  The Schwarzschild criterion is used, with no overshoot and no semiconvection.
- **Grey Eddington photosphere** (T = Teff at tau = 2/3, P_gas = 2g/(3 kappa)), placed at r = R,
  m = M, using the same opacity fits.

## What the numbers do not show
1. **The real Sun's convective envelope.** The 1 Msun model has an envelope convection zone of only
   1.2e-8 M, with its base at r/R = 0.991. The real ZAMS Sun's zone holds about 0.02-0.03 M, with its
   base near 0.7 R. Diagnosis: the zone is confined to where the H^- opacity (rising as T^9) meets
   the falling Kramers opacity, at T of about 1e4-5e4 K. Below that, kappa is about rho T^-3.5, so
   nabla_rad settles at about 0.22, the radiative-zero value 1/(n+1) for Kramers, which is well below
   nabla_ad = 0.4. A deep solar convection zone needs (a) the large real opacity at 1e5-2e6 K (H, He
   and metal bound-free and lines, which is 10-100x the Kramers fit there) and (b) the lowered
   nabla_ad in ionization zones. Neither is in this model. So the thin zone is an expected
   consequence of the physics as agreed, not a code bug (re-solves reproduce it to <20 %). It must
   not be quoted as a prediction of the solar convection-zone depth. The same applies to all
   envelope-zone depths in the sweep, which are underestimated by orders of magnitude for M <~ 1.3
   Msun, and to the mass (1.6-1.8 Msun) above which envelope convection "disappears".
2. **Luminosity and radius of the 1 Msun ZAMS star.** L = 1.11 Lsun against about 0.7 Lsun real, and
   Teff = 6090 K against about 5600-5700 K real. The low analytic opacity lets energy out too easily.
   These are the crude-physics differences of tens of % that PLAN.md anticipated, not validations.
3. **Low-mass end (0.5 Msun).** Teff = 2780 K and R = 0.72 Rsun, against about 3700-3800 K and
   0.45 Rsun real. This end lies outside the validity of the H^- fit and of a grey atmosphere with no
   molecules. A real 0.5 Msun star has a deep convective envelope (about 0.2-0.3 of the mass), not
   1.7e-3. R and Teff at M < 0.8 Msun are not trustworthy. L is much more robust (see 6).
4. **High-mass surface.** For M >~ 1.5 Msun the photosphere (Teff 1e4-2e4 K) uses the H^- fit far
   outside its range (it was made for about 3000-6000 K). Surface opacity and Teff there are
   unreliable. L and the interior are not affected (see 6).
5. **Onset mass of the convective core.** At the centre, nabla_rad/nabla_ad is 0.94-0.97 for every
   model from 0.5 to 1.0 Msun, crosses 1 between 1.0 (0.972) and 1.077 Msun (1.019), and reaches 5.9
   at 5 Msun. The onset mass (1.0-1.08 Msun here; about 1.1-1.2 Msun in published Z = 0.02 models) is
   therefore sensitive at the few-% level to central opacity, screening and composition. It is
   bracketed only by the sweep spacing. "No convective core at 1 Msun" is a 3 % margin, not a robust
   statement. Numerically, zone boundaries are now root-found on the ODE solution (commit 0c2e342):
   the core at 1.077 Msun (0.0024846 M) is the same to 8 digits for profile grids of 150-2400 points
   and fitting points 0.3 M and 0.5 M. But it is physically very sensitive near onset: a 0.03 %
   change in M shifts it by about 2 % (d ln M_core / d ln M of about 60). Small-core masses near onset
   should therefore be read as "a small core exists", not as precise values.
6. **L is insensitive to the surface treatment.** L is identical with and without H^- (to <6e-4),
   and only R and Teff change. This is the classical result for stars with radiative envelopes. It
   would **not** hold if the envelope were deeply convective (real M < ~0.8 Msun), so L at the
   low-mass end carries the same caveat as point 3, though more weakly.
7. **L-M slopes.** The local slope of 5.4 at 0.5 Msun is close to the Kramers homology limit
   (L ~ M^5.5 R^-0.5). It falls to 3.4 at 5 Msun as electron scattering takes over, but has not
   reached the pure-scattering limit of 3 within the sweep. Radiation pressure is negligible
   (beta_c >= 0.98), so the slope's approach to 3 and its turnover above 5 Msun were not tested. The
   overall exponent of 4.44 belongs to this opacity law. It is not a fit to observed stars.
8. **Not tested:** masses outside 0.5-5 Msun, other compositions, rotation, magnetic fields, mass
   loss, and any time dependence (no ages, lifetimes or evolutionary tracks).
