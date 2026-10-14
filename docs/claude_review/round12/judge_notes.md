# Judge notes, round 12

## Candidate finding J5 (uniform-Q0 collision hypothesis fails on Sylvester boxes)
- ZE-TC^coll_eps(Q0) (referee_t1rigour §6.2; also author's (Q1)+(Q2)): single radius sigma=(1/2)e^{-Q0 abar/2}
  for all collision subspaces with >= N2/Q0 pairs.
- Sylvester n=2^k, m'=ceil(4w), coordinate columns 1[x_i=0] widths <= 2j. If 8j < m', the 4-col span L4 (16 classes)
  is cheap -> in TC with kappa=0. If Q0 < N2/P(L4) (~16.06), generic point of L4 is >= 1/30 from all members of the family
  -> fails at all levels g >= 2ln15/Q0.
- If Q0 >= N2/P(L4), single-column span L1 (2 classes, o >= m'-2j) is in family, but points at distance ~sigma from L1
  have kappa <= 4^{m'} eta^{o}; violation iff o*Q0 > d roughly, i.e. 16(m'-2j) > d  (n < ~48w).
- Expected: fails for k=12..17. Per-subspace radii sigma_pi = (1/2)exp(-N2 abar/(2P(pi))) (Lemma 6.6 / Thm 6.8 part 1)
  avoid this: single column has N2/P ~ 2, o*2 < d.
- Hierarchy: ZE-TC^coll_sigma(per-subspace) => MIN(abar) => MIN^pair (Lemma 3.3 / gap region inside {Gamma<=abar}).

## Distance lemma for test point u = y + eta v (y in L generic, v distinct values within classes, grid spacing dv)
- partition pi' with P(pi') >= N2/Q0: either a block mixes L-classes (dist >= (gap - 2 eta)/2) or refines;
  b_max - 1 >= (n-1)/Q0; spread >= (b_max-1) dv eta; dist >= spread/2.

## J5 verified (judge_checks/uniformQ0_sylvester.py, persubspace_sylvester.py)
- uniform Q0: fails k=12..17 for every Q0>=1 at every g >= 5.4161 (< near-extremal level 5.52 at k=12). none from k>=18.
- per-subspace: fails k=12 (b=2,3) from g>=14.90 (W>=1.22e5 ~2.7 W_min); no failure k=13,14 from coordinate spans.

## Idea J6 (apolarity / fat points)
- J(A) contains J_w := (<k',v>^{h_A(v)+1} : v integer). Apolarity (Emsalem-Iarrobino): (J_w)^perp_D = degree-D forms
  vanishing to order >= D - h_A(v) at [v]. Hence TC(A,m) = {u: l_u^m in J(A)_m} contains Bs_R(L_m(fat points [v]^{m-h(v)})).
- Linear span of points in base locus iff sum h < m  == cheap. RNC through k+3 pts in P^k: sum h < 3m (heuristic recall).
- One cut: J(A_S(L)) = (x^3, sigma^{2L+1}) by SLP of k[x]/(x^3) (Stanley) + count => one-cut thm = SLP.
- TEST: is J(A) = J_w in degrees <= m for M8, one-cut, two-cut? (mod p cumulative equality proves equality over Q.)
- jw_test small: one-cut (3 boxes), two-cut (2 boxes), Syl8 L=1: J(A)=J_w up to m (exact via mod-p sandwich).
- M8 (L=2,m=8): J(A) != J_w from degree 6: dims J(A) [..,637,1541,2961] vs J_w [..,616,1485,2940]; 21 extra in deg 6 = C(7,2)?
  -> test J_proj2 = J_w + pullbacks of J(pi_F(A)) for pairs of columns.
- jproj_test: pairs of columns add nothing (mod p). jproj3_test: J_w + 3-col projections = J(M8) in all degrees <= 8 (EXACT by sandwich); 4-col same.
- jw_M8_p2: second prime same gap. jw_M8_exact6: EXACT over Q: dim(J_w)_6 = 616 < 637 = dim J(M8)_6 (gap 21).
- onecut_slp: J(A_S(L)) = (x^3, sigma^{2L+1}) in ALL degrees, 6 boxes, Hilbert fn = Lefschetz truncation (exact).
- README §0-§2 written.
