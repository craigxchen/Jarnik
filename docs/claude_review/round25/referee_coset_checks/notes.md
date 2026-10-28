# referee notes (coset.md)
## Proof checks (done by hand)
- Lemma 1.1 OK. Lemma 1.2 OK (exact Gaussian factorisation argument).
- Lemma 2.1 (i)-(iv) re-derived: constant c_2 = pi(1+sqrt2) correct; |mu|_1/|mu|_inf <= m; 2|mu|_1 eps/pi < b*|v|/2 (eps<pi/4). OK.
- Thm 2.2 OK (q_j = preceding split prime; Prop B for 5).
- Cor 2.3 OK (Rankin z^2=1/tau; f + z f' = f(1+2rz/(1-z^2)); safety from y_v >= P^{-|v|/2}). const 2.543 at tau=2 OK.
- Rem 2.3(2): refutes FULL-TORUS Poisson extrapolation (fakes §5.1) asymptotically. "open point 3 false for coset T^r":
  T^r is not a realisation coset of any profile with M>=2 (L != 0). So open point 3 as intended (realisation cosets) untouched. wording issue.
- Thm 3.1 OK (trivial given Lemma 2.1). Note: converse "positivity => slack existence" NOT proved; sandwich is honest.
- Thm 4.1 OK. Cor 4.2 formula re-derived OK; factor 1/2 by v->-v symmetry OK.
- Thm 5.1 OK; scope: bounds of the form sum_{mu in F_k} U(mu) (union bound over mu within a step). Does not cover
  exact computation of mu_1(A_k), non-uniform (not inf over theta_<k) integration, or non-Fubini arguments. Correctly labelled "union-bound certificate".
- Cor 5.2: f(M)=cM^{d/2}-log(1+2M^2/k) decreasing-then-increasing; "last zero" OK. Numbers to check.
- Prop 6.1 OK, needs N(beta_j) <= p_j (as stated). Summary/§6 "non-arithmetic" must mean "not a transplant with N<=p" (arg of Gaussian integers of large norm are dense).
- Prop 7.1 OK. Lemma 7.2 OK (incl. odd pi/4 shift swaps families; bound 1/(sqrt2|gamma|) >= tau^{|v|} y_v uses 2^{|v|/2}/sqrt2>=1).
- Cor 7.3 OK as stated (rho explicit).
- ISSUE: Reading of 7.3, §0 item 6 and (Q_seed): heavy mass "e^{-W/4+O(r)}" and heavy chars "w(v)>=W/2-O(log W)"
  drop log(C kappa_inv) (right-inverse norm: denominator cost). Hadamard allows log kappa_inv up to O(r log r).
- ISSUE: summary "A_0 window question reduces to (Q_seed)": only sufficient direction, for seeding mechanism.
- No reduction to route 5 or one-point caps claimed in coset.md.
## Computations: see below
## Computations done
- constants.py: M_0 reproduced (4.48e3, 4.55e5, 6.99e8, 4.93e17, 4.33e20). Cor 2.3 reproduced with author's tau,
  BUT tau for 13..73 is sqrt(41/37)=1.0527 not sqrt(73/61)=1.094 -> n0=1859, log mu >= -3.2e4 (minor error §8.1/§8.3).
- graph_check.py: Cor 4.2 formula vs brute force rel err <= 1.4e-7 (truncation); Thm 5.1: sum U_k (truncated, lower bd) = 7.23 >= 2.53.
- cover_check.py: 60 random lines certified (min overlap 6.0e-4), eta=1e-3 line certified (4.39e-6) at T=20 and 24. Same digits as author.
- char_check.py: eta=1e-3 line characters hold to T=40 (min margin 1.1e-6), eta=0 to T=40 (2.0e-6). Truncated: evidence.
- density_ref.py: 0.001,0.0915,0.4845,0.9195,1.0 at eps=1e-3..1e-7 (normalised units); uniform 0. Matches §8.1.
- seed_check.py: Lemma 7.2 exact-Gaussian test, r=4,m=2, all 256 g0 in (pi/4 Z)^4, 588 non-chars w<=22: min rel margin 0.026>0.
## Extra findings
- Borel-Cantelli: nu(G(p) ∩ C)>0 <=> exists phi in C with characters (weights p) and non-characters slack (p/tau^2). Equivalence;
  coset.md's "by Thm 3.1 such a point is necessary" is not what Thm 3.1 proves, but true by B-C.
- open point 4 parenthetical (genuine background fibre) NOT settled by Thm 2.2.
- "supersedes A'" only qualitatively.
