# notes for split.md
## N1 setup (derived)
- split torus T=G_m^5/G_m over Q; points [u], u primitive integer 5-tuple (gcd 1), all u_j>0 near e.
- -K polytope R_5=conv(e_a-e_b). H_{-K}([u]) = (max u/min u) * lcm(u)/gcd(u) (judge J2 exact).
- dist([u],e) ~ D/U, D=max-min, U=min.
- A-point exponent gamma: dist^gamma H <= 1  <=>  (D/U)^gamma * (maxu/minu) * L <= 1.
- Divisor dictionary: d_j e_j = N with gcd(d)=gcd(e)=1  <=>  N=lcm(d). So [u] <-> N=lcm(u) has the 5 divisors u_j.
  Symmetric case u ~ sqrt(N): A-point <=> D <= N^{1/2 - 1/gamma}; gamma=4 <-> D <= N^{1/4}: Erdos #887 scale.
- pairwise Liouville: L >= prod u / prod_{a<b} gcd(u_a,u_b) >= U^5/D^10 -> gamma >= 10/3 (CC analogue).
- Inversion u -> L/u is an automorphism of (LM_5,-K,e).
- 1-PS orbits t^v with distinct coords: width_v = 2(max v-min v) >= 8. Useless.
- Bounded covering families of rational curves through e: alpha>=4 (Prop 4.3; check geometric).
## N2 literature (WebSearch/WebFetch 2026-10-28)
- erdosproblems.com does not resolve. Snippet: #887 (Erdos-Rosenfeld ErRo97): absolute K s.t. for every C, large n have <=K
  divisors in (n^{1/2}, n^{1/2}+C n^{1/4})? OPEN. ErRo: <= 1+C^2 (trivial: s=d+n/d integer in (2sqrt n, 2sqrt n + C^2)).
  ErRo: infinitely many n with 4 divisors in (sqrt n, sqrt n+16 n^{1/4}).
- github google-deepmind/formal-conjectures PR #7218: rosenfeld_4 variant false: m-czech forum post July 2026 (thread 887
  post 7374): n=(x^2-1)(x^2-4)(x^2-16)(x^2-49), z^2=56x^2+65 Pell; 4 divisors from balanced splits of {+-1,+-2,+-4,+-7},
  fifth = x^4-35x^2+56+2xz. C=31 (Lean proof uses C=40). Main #887 question still open per PR.
## N3 reduction split A_5(delta) => absolute-K #887, uniformly in window position y (derived, to write as Thm)
- divisors d_1..d_M of n in [y, y(1+C n^{-1/4})]; random 5-subset S; L_S=lcm(d_S)/gcd(d_S) = H/(max/min).
- per prime p^t||n: max-min = sum_{s=1}^t 1[S splits at threshold s]; P(no split) >= q^5+(1-q)^5 - O(1/M) >= 1/16 - O(1/M).
- E log L_S <= (15/16+O(1/M)) log n.  dist <= C n^{-1/4}.  => alpha <= 15/4 points always exist (gamma=15/4).
- with gamma0>15/4 positive proportion good; good subsets in Z (n large); Schwartz-Zippel => M bounded by M0(Z,delta), K absolute.
- formal-conjectures 886.lean: #886 (Ruzsa): eps>0, divisors in (n^{1/2}, n^{1/2}+n^{1/2-eps}) is O_eps(1)? open. ErRo97 4 divisors
  in window 16 n^{1/4}; ErRo bound 1+C^2. Snippet: known for eps>1/4 (trivial since window < n^{1/4}).
- Split CC analogues: Cilleruelo-Jimenez-Urroz "Divisors in a Dedekind domain" Acta Arith 85 (1998) 229-233 (x^2-dy^2=N, arcs
  N^alpha, alpha<=1/4-1/(8[k/2]+4) contain <=k pts); "The hyperbola xy=N" JTNB 12 (2000) 87-92; Granville-Jimenez-Urroz
  "The least common multiple and lattice points on hyperbolas" QJM 51 (2000) 343-352. Chan "Numbers with three close
  factorizations" (arXiv 1303.2249, Integers 18 (2018) A13); Chan-Holmes-Liu-Villarreal arXiv 2508.02818 (four close
  factorizations, Pell, 0.04742); arXiv 2507.07094 (three close factorizations, central lattice points). Letendre arXiv
  2503.12146 (Mar 2025) "Divisors of an integer in a short interval": Prop 1 bound for (n^theta,...) ; Conjecture 1.
  Hosts numdam, colgate, arxiv, themoonlight, althofer, erdosproblems DO NOT resolve; only snippets + github.
- #886 link: averaging gives H<=n^{15/16}, dist<=n^{-eps}: A-pt exponent 15/(16 eps). split A_5(delta) => #886 for
  eps > (15/64)(1+2delta) (range below 1/4!).  Liouville tau=10/3 gives only eps>9/32>1/4: nothing.
## N4 curve ratio LP (polynomial blocks): ratio = B.C/mu, B.C=sum over roots (max-min ord). 2|3 profile ideal: 20/6=10/3;
   balanced ideal: 30/8; with residual degree r and block degree s: 2|3 <4 iff s>r; balanced <4 iff s>2r. Overdetermined.
- TODO: cyclotomic-product curve families x^m-1 (unbounded degree rigid curves) -> ratio? could they be dense with ratio<4?
## N5 Lorentzian vs Euclidean window
- torus metric: dist = Delta d / X (multiplicative). Window at position X: Y = C X n^{-1/4} (Lorentz-invariant arc of xy=n:
  hyperbolic angle C n^{-1/4}; boosts (x,y)->(lam x, y/lam) = split rotations). lit18 §2.2 uses Euclidean Y~X^{1/2}; agree only at X=sqrt n.
- trivial ErRo bound at vertex: s = d + n/d integer. For lam^2 = p/q: s = q d + p e, bound ~ 1 + C^2 sqrt(pq). General lam:
  Dirichlet pq ~ n^{1/4}: useless. => hyperbolic phase log(X/sqrt n) is the split drifting phase; vertex = phase zero.
- S-CG (Lorentzian split CG): for each C, #{d|n: X<=d<=X(1+C n^{-1/4})} bounded uniformly in X. Known at X=sqrt n (1+C^2).
- Thm1: split A_5(delta), delta in (0,1/30) => S-CG with ABSOLUTE K (=> #887 parts.ii) and #886 for eps>(15/64)(1+2delta)+.
  proof: averaging (threshold indicators, hypergeometric), Markov on Y=1-logL_S/log n, Northcott, Schwartz-Zippel (deg F M^4).
## N6 KEY EXACT FACT (pell5.py, pell5_curve.py)
- m-czech tuple: D0=(x-1)(x-2)(x-4)(x+7) [off 0], f5 [off 2sqrt56-35=14.97], P_{2,4} [15], P_{1,4} [18], P_{1,2} [30]. C_actual=30.
- As a curve in split LM_5 over K=Q(conic)=Q(t): 16 rational places over x=+-1,+-2,+-4,+-7; vectors 0/1; exactly ONE place
  (7,53) constant; B.C=15, each of the 15 cuts hit exactly once (cut-balanced, 'kappa=1/2', impossible in twisted where
  E_A.C=E_{A^c}.C). Infinite places: 2 conjugate real-quadratic branches (z/x=+-sqrt56), differences O(x^2) vs x^4: mu0=2 each.
- alpha_C = B.C/(2 mu0) = 15/4 (Liouville quadratic irrational on P^1 sharp via cont. fractions/Pell). V-null (D_V.C=40=10*4 with half-width conv.)
- numerics gamma -> 3.9154 at x~1.4e14, log L/log n -> ~0.928 (15/16=0.9375 with constants). consistent.
- Two-branch covering bound: (2 + 3*2) mu0 <= B.C => alpha>=4 for covering families; so family can't cover.
- General (a,b,c), d=a+b+c: f5 = x^4-(S/2)x^2+m+2xz, m=abcd, z^2=Ax^2+B, A=(S^2/4+2m-e2)/4, B=(e3-Sm)/4 (e_i of squares).
  A r^2+B automatically square at roots => 16 rational places. 2-param family of curves -> 3-fold W (to check Jacobian rank).
- Density of 15/4 points requires infinitely many such W (unbounded degree patterns) -> K3-like density question; Vojta says no.
## N7 family results (pell_family.py)
- N(f5)=n identically; A = abc(a+b+c) = m; B = product of three quadrics /4; A r^2+B = square at each root (explicit).
- B.C=15, 15 cuts, 16 places at (1,2,4),(1,3,5),(2,3,7),(1,5,6),(1,3/2,5/3). Jacobian of (a,b,c,x)->torus: rank 3 at 4 random pts
  (one zero sv = scaling). => image W is a HYPERSURFACE (3-fold) swept by alpha=15/4 rational curves (2-param family).
- conic z^2 = m x^2 + B has rational points (at roots) -> rational curve over Q; sqrt(m) irrational generic -> Pell-quality pts.
- twisted: Galois pairs places 1_A <-> 1_{A^c} (real places have v=0), so per-cut weight even; balanced B.C>=30, needs mult 8.
- golden-8 (record) in torus language: 4 blocks = 4 cuts only, 5-subsets of the 8 odd words can't share a bit -> alpha=4 (critical).
## N8 misc
- scan.c to N<=2e8 (lcm=N, primitive, Lorentz window N^{-1/5}): no point with gamma<=4.0; min gamma 4.2772 (2^27). m-czech
  first member x=3656 has log n=65 >> log(2e8)=19: scan uninformative for 15/4 (low height, constants dominate). Label data only.
- hypergeom.py: min_K p(M,K) >= 1/16 - 2/M for M<=3000 exact; analytic: convexity of C(K,5) (2nd difference C(K,3)>=0)
  gives p >= (1/16) prod_{i=0}^4 (M-1-2i)/(M-i) >= (1/16)(1 - 15/(M-4)).
- Remark: CB(split LM_k at e) for all k => Ruzsa #886 in full: k-subsets, E log L_S <= (1-2^{1-k}) log n; eps>(1-2^{1-k})/(k-1).
- split cut-balanced Liouville: dist >= max_{ordered (a,b)} 1/prod_{A ni a, not b} n_A; sum over ordered pairs = sum_A |A||A^c| log n_A
  = 80w under cut balance; average 4w; log H = 15 w -> alpha>=15/4. (n_A, n_{A^c} separate heights in split: 30 ray heights.)
## N9 write-up done (split.md §0-§8). Status partial: Thm1 (split A_5 => abs #887 phase-uniform, #886 range), Thm2 (Pell
hypersurface W of 15/4 curves; parity lemma excludes in twisted). (b) undecided. (c) no mechanism.
