# notes (round 11 compute)

## N1. Reformulations (to be proved in compute.md)
- d = n-1, coordinates k' in Z^d, k_n = -sum k'. Positions: u_x = delta_x - delta_n (x<n); <k,delta> = <k',u>.
  Degenerate <=> u in a rational hyperplane of R^d (some nonzero k' in Z^d with <k',u>=0).
- Sign column beta in {+-1}^n: <k,beta> = 2 * sum_{x in S} k_x  where S = {x : beta_x = +1} (sum k = 0).
  So slab |<k,beta>| <= w  <=>  |k(S)| <= w/2 =: L. Width h_A(beta) = 2 * range_A(k(S)) <= 4L = 2w.
  Supercritical condition (README (ii)) 2w < m  <=>  h_A(beta_p) < m for each p.
- V(A,m) = annihilator of Poly_{<m}|_A. Hence base locus Bs = {w : (w^k)_{k in A} in Poly_{<m}|_A}:
  "the exponential k -> w^k agrees on A with a polynomial of degree < m".
- Let I(A) = vanishing ideal, J = gr I(A) = ideal of top-degree forms. Then V(A,m) != 0 iff J_m != S_m
  (iff (S/J)_m != 0). Tangent cone forms Phi_c span the Macaulay inverse system (J^perp)_m (= least space
  degree-m part, de Boor-Ron). TC := common zeros = {u : l_u^m in J_m}.
- If l_v takes s values on A (rational v), prod(l_v - val) in I(A), so l_v^s in J; v in TC if s <= m.
- A = -A (symmetric) => V = V_even + V_odd; on real torus the equations are sum c cos = 0 (even c),
  sum c sin = 0 (odd c).
- d=2: at most one cutting slab direction can satisfy 2w<m with V != 0 (two independent l with l^{2L+1} in J
  and m >= 4L+1 force J_m = S_m). n=3: slab on k_1 or k_2 gives a product box (Prop 5.3 / R1 territory);
  slab on k_3 = -(k1+k2) gives the hexagon {|k1|,|k2|<=T, |k1+k2|<=L}.

## N2. Tools
python-flint 0.9 available (fmpz_mat nullspace/lll, fmpz_mpoly gcd/factor, arb/acb). sympy 1.14, numpy, scipy, mpmath.

## N3. Log of runs

## N4. Key structural facts found (to be written as Props in compute.md)
- (Segment/parallelepiped criterion) If A contains a lattice parallelepiped P = k0 + sum t_i v_i (v_i lin. indep.,
  0<=t_i<=L_i) with sum L_i >= m, then c = Delta^alpha (alpha<=L, |alpha|=m), entries <= 2^m <= H, lies in V_H and
  sum c_k w^k = w^{k0} prod (w^{v_i}-1)^{alpha_i}: nonzero at every point with all w^{v_i} != 1. So no
  nondegenerate point of the WHOLE complex torus lies in Z_P / Bs. ZE trivially true then.
- Define pi(A) = max sum L_i over parallelepipeds in A; rho(A) = max{j : (S/J)_j != 0} (V(A,m)!=0 iff m <= rho(A)).
  pi <= rho. Regime-faithful model needs pi(A) < m <= rho(A).
- In the regime, for an actual cluster's own profile, pi(A) < m is forced by the pair/character inequalities
  (Prop 5.3(ii)-(iii) argument): P in A with sum L >= m gives sum_i L_i 2X_{v_i}/W <= m/2 => some X_{v_i} <= W/4.
  So models with pi(A) >= m (e.g. m <= 2T: coordinate segments longer than m) are NOT regime-faithful.
- Full-rank profiles with uniform L at d=2 are inadmissible (complete intersection socle 4L < m).
- explore3: n=3 admissible only for single-slab profiles; hexagon results: base locus curve part = (w0-w1)^j
  [+ (w0-1),(w1-1) when m <= 2T+1]; tangent cone gcd = (u0-u1)^j [+u0,u1]. All rational.
  BUT those instances have m <= 2T => coordinate segments of length 2T >= m... need re-check pi vs m.
- (Section principle) For a rank-j sublattice Lam' with basis v_1..v_j and k0, A' = A cap (k0 + Lam'), viewed in
  t-coords. w in Bs(A,m) => (w^{v_1},...,w^{v_j}) in Bs(A',m). If Bs(A',m) near e is contained in finitely many
  proper subtori of T^j, then every w in Bs(A,m) near e satisfies a character relation => degenerate.
  So ZE-real(A',m) for one section implies ZE-real(A,m). Parallelepipeds w/ sum L >= m are the base case.
  Regime: j-dim sections inside box have rho(A') <= 2Tj (CI of coordinate functionals) so need j >= m/2T ~ 0.37 n.
- Hexagon (n=3, slab on base point): rho = 2T+L (explore3), pi = 2T. Faithful window m in (2T, 2T+L].
- n=4 Hadamard (= skew double core q=3) profile: A ~ {(a,b,c) in [-L,L]^3, a+b+c even}, rho=3L < 4L = h(beta):
  condition (ii) infeasible, even (ii') range=2L<m<=3L feasible.
- ZONOTOPE certificate (generalises parallelepipeds): Z = k0 + sum_i [0,L_i] v_i (any number of generators, dependent
  allowed), lattice points k0+sum t_i v_i all in A, sum L_i >= m. c = prod_i Delta_{v_i}^{L_i} (truncate to total m)
  annihilates Poly_{<m}, entries <= 2^m, sum c_k w^k = w^{k0} prod (w^{v_i}-1)^{L_i}. Nondegenerate w excluded.
  Widths additive => regime argument identical to Prop 5.3(ii): excluded for actual nondegenerate clusters.
  zeta(A) := max sum L_i. Faithful needs zeta(A) < m <= rho(A).
- Hexagon IS a zonotope: (1,-1)x(2T-L), e0 x L, e1 x L: sum = 2T+L = rho. So hexagon not faithful. Every
  centrally symmetric lattice polygon is a lattice zonogon; conjecture rho = sum L_i => no faithful d=2 model.
- bs2d hexagon T=4 L=1 m=9: dimV=1, P = (w1-1)(w0-1)(w0-w1)^7 exactly = zonotope finite difference.

## N5. Results so far (runs)
- explore3/n3_hexagon: n=3 only single-slab profiles admissible; hexagon rho = zeta = 2T+L (Prop 2.4 proved).
- scan4+faithful4: ALL 355 distinct admissible n=4 slab-boxes (all full-rank pattern subsets, T<=4, L<=3,
  windows II or II') have zeta = rho  => no faithful model at n-1=3.
- scan5+faithful5 (12 random full-rank profiles, T<=3): running; so far zeta = rho everywhere.
- gap_search d=3 (60 random), d=4 T=1 (40), T=2 (25): zeta = rho always.
- KAPPA BOUND (proved): Hadamard profile order n: sum_{p>=2} |<v,h_p>| >= n for v != 0 in Z^n_0
  (Parseval + |<v,h_p>| <= |v|_1 <= |v|_2^2). So zeta(A) <= ((n-1)/n) max_p h_A(h_p) = 2((n-1)/n) rmax.
- had8 T=1: L=2: |A|=1037, rmax=4, rho_p=8 -> zeta <= 7 < 8 = rho: FAITHFUL model m=8 (II' only: h=8=m).
            L=3: |A|=1877, rmax=6, rho_p=11 -> zeta <= 10 < 11: faithful m=11.
- M8 = had8 T=1 L=2 m=8: dimV=42 exact, sup 6. TC forms: 42 octics, monomials need >=4 vars; coordinate 1,2,3-planes
  all in TC, 7 coordinate 4-planes (the S_p not containing base pt) in TC.
- tc_search (BFGS + mp refine) and bs_search: zeros converge badly (singular components); apparent "nondegenerate"
  zeros are unconverged points near 4-sparse planes (1/2,1/2,0,0,1/2,1/2,0) or paired (0,a,b,c,c,b,a). Not trustworthy.
- bs_search2 (tau=tan(u/2) polynomial, relative residual, LM): eps=0.1, 40 starts: 30 zeros, all with integer
  relations (|v|<=30, loose tol). Need exact verification of candidate subtori.
- Exact subtorus test: exp(H) subset Bs iff pushforward of every c in V under y(k)=(<k,g_j>) vanishes.
- Section principle useless for M8: sections dim<=3 have rho<=6; 4-dim sections with rho=8 must be full grid => zonotope.
- LP dual bound for zeta (lpzeta.py, exact rational certificate): had8 L=2 -> 7.
- M8' = had8 T=1 L=3 m=11: |A|=1877, dimV=7, sup 4, rho=11 exact. zeta <= 10 (kappa). FAITHFUL.
- FACTOR STRUCTURE (exact): every basis P_c = monomial * binomials * non-binomial factor(s). M8': 3 basis polys
  = (binomials of order 8) * cubic(12 terms); 4 = (binomials order 7) * quartic(24 terms). Every non-binomial
  factor is a FORM in x = w - 1 (homogeneous). => Bs_nondeg near e subset of the CONE V(G_1..G_7) in x-coords.
  M8: some non-binomial factor is NOT homogeneous in x (assertion failed) -> handle separately.
- !!! EXACT: M8' cone_sat: V(G_1..G_7) (x=w-1 cone) subset union{x_j=0} cup {x_i=x_j}  (unit ideal after
  saturating by x0..x6, x0-x1,...,x2-x3). Binomials of M8' basis: w_j-1, w0-w1, w3-w4. => Bs(M8') subset of the
  "collision" subtori {w_j=1} cup {w_i=w_j}: every point of Z_P has two coinciding positions. GLOBAL, complex.
- !!! EXACT: M8: 24 of 42 basis polys have form-type non-binomial parts; their cone saturates to unit ideal
  (237 s). Binomials: 27 types incl w0-w1*w6 etc. => Bs(M8) subset union of 27 binomial subtori + collision
  subtori. ZE holds globally for M8.
- explain.py numerics (superseded): 117 explained, 30 'unexplained' = poorly converged points on rational
  subspaces (relation detection failures).
- fakes_m8: 100/100 random U[0,1e-3]^8 positions certified off Z_P.
- model8.py T=1 L=2: sdc7 (skew double core q=7, 14 patterns): |A|=757 rho=8 LP zeta<=7: faithful m=8, dimV=1 (sup 42).
  pal7 (15 patterns) |A|=807 rho=7, LP 7: not faithful. rand1 rho 9 LP 9.5; rand2 rho 8 LP 8.5: undetermined (LP not < rho).
  sdc7 L=3: rho=10, LP 10.5: undetermined.
- sdc7 m=8: single P_c = mono * G(w-1), G irreducible octic FORM (687 terms in w). c even. Phi = -40320 * irreducible
  octic (196 terms). Exact sign scan on [-2,2]^7 integer points: Phi > 0 never, < 0 at 65604, = 0 at 12520.
  So Phi <= 0 apparently (negative semidefinite?). Float minimisation too noisy. Testing factorisation over Q(i).
- faithful5: 55 instances zeta = rho; instance 56 (T=3,L=3, |A|=877) hangs (LP bound 13.5 vs rho 13).
- !!! sdc7 m=8: Phi/8! has 196 terms, coefficients in {-2,-1,1,3}; EXACT witnesses: Phi/8!(u+) = 1.98e28 > 0 at
  u+ = [1020,-1592,-6221,-2363,5574,2741,-3682]; Phi/8!(-1,-1,-1,-1,0,0,0) = -1 < 0. Phi irreducible over Q.
  => ZE-real FAILS for sdc7 model (codim-1 linear system): F = sum c cos changes sign near 0; Baire + irreducibility
  argument gives nondegenerate zeros arbitrarily close to e. (Proof in compute.md 6bis.)
