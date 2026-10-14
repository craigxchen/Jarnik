# fake.md working notes

## plan
- profiles: sdc (skew double core, round5 exact Thm 4.1), Paley twin profile (round5 construct §2; copies collapse at T=1),
  Sylvester H8/H16.
- T=1 box: A = {k' in {-1,0,1}^d : |k(P)| <= L for every relative pattern P (base point removed)}.
- tangent forms span = J_m^perp, dim H_A(m) (not dim V). f(u)^2 = sum_i Phi_{c_i}(u)^2 (ON basis of V) = ||G w(u)||^2,
  w(u)_a = (m!/a!) u^a over reduced (a_x<=2) degree-m monomials.
- kappa duality (nearstrata Prop 2.1): (D/2)^m kappa(u) = best uniform approx error of l_u^m on A by deg<m polys.

## key conceptual lemma (transfer)
- fake profile (patterns, weights) strictly admissible with norm-proportional margin sigma_lambda >= rho1 ||lambda||_1.
  Box A, V_H depend only on patterns and r (number of primes) and L; scaling tau = log p_min (primes in window) keeps A.
  So any nondegenerate u0 in TC_H with |<lambda,u0>| >= e^{-K||lambda||} gives gap-admissible fakes at tau large with
  direction u0 and kappa = 0.  At fixed W = O(n log n): need Gamma-type control.
- dim V small + regular nondegenerate common zero of Phi basis -> IFT -> Z_P contains (d - s)-dim nondegenerate manifold
  near t u0 -> fakes on Z_P (a.e. point gap-admissible by measure argument, tau large).
## data so far
- sdc7 L2: H=[1,7,28,77,161,251,188,43,1], zeta in [6,7] (LP full=7), m=8 dimV=1; cert zero sep>=0.249565 (sdc7_cert)
- sdc11 L2: H=[1,11,66,275,880,2243,3735,2157,5] (mod p 2 primes), zeta in [6,7]; V(A,8) exact invariant dim 5.
- pal7 L2: H=[1,7,28,77,161,251,241,95], rho=7, zeta<=7 (LP). pal7 L3 == Paley core box == M8' (perm (0,1,2,4,3,6,7,5)).
- pal11 L2: blocks s0=[1,3,14,55,176,449,932,703,88], s1=s2=[0,2,13,55,176,451,933,701,87] (Z5 mult 3) -> H(8)=436, rho=8, zeta<=7 (LP 22/3). faithful, overdetermined!
- tcnum sdc7 m7: kappa2*(s): .05:8.3e-8, .1:2.1e-6, .15:1.8e-5, .2:4.6e-5. pal7 m7: s=0 ~0 at collisions, .05:1.3e-7, .1:3.7e-6
- sdc11 m8: certified zeros sep>=0.005287 (Gamma 3.03), 0.010008 (Gamma 4.03), (F) ok -> Thm 2.2 -> fakes on Z_P.
  landscape: zeros up to s=0.01, positive beyond (4e-14 at .02 ... ); Newton from s=.02/.05 minimisers -> sep 7e-4, 3e-3.
- sdc11 m7 landscape: .02 6.5e-9, .05 3.9e-7, .08 6.6e-6, .12 8.9e-5, .16 2.1e-4 (progressions + far point).
- pal11 m8 landscape s=0: zero at collision (sorted [-0.65, 0 x7, .08,.09,.11, 1]).
