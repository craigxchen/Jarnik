# notes for coset.md (round 25, coset positivity)

## Key insight 1 (density points of transplants)
transplant phi_beta (beta_j Gaussian primes of distinct split-prime norms q_j < p_j) has slack factor tau^{|v|_1},
tau^2 = min p_j/q_j > 1. Unsafe monomials in a cube of radius eps have fraction <= c sqrt2 r (1+|v|_1)/(tau^{|v|_1}-1),
summable over v; each v safe for eps small => density point => mu(G(p)) > 0 iff 5 not in p (Prop B for converse).
Supersedes Prop A' on the full torus; answers fakes open point 4.

## Key insight 2 (sandwich on cosets)
G(p/tau^2) cap coset nonempty => nu(G(p) cap coset) > 0 => G(p) cap coset nonempty.
Coset version of density lemma uses frequencies mu = B^T v; characters constant on coset.

## Key insight 3 (barrier)
Sum over steps of fibre-union masses >= kappa (Pi_J - 1)/2 for ANY basis; method needs each < 1;
window with m >= cM => impossible for M >= M_0(c,delta).

## Key insight 4 (circularity of arithmetic seeds on realisation cosets)
arithmetic point on Real(delta) with N(beta_j) <= p_j => actual cluster on norm N_q <= N_p, arc const C (N_q/N_p)^{1/4}.

## Heuristic (cell first moment)
log E[# surviving cells at scale e^{-X/2}] ~ m X/2 - S(X). Full torus threshold Lambda ~ (1/4) log P + O(loglog P).
In window: S ~ e^{2 Lambda_J} >> m Lambda log P: no Haar-generic survivors on cosets. Full-torus positivity
comes from arithmetic density points, which on realisation cosets are actual clusters.

## Results of computations so far
- density_check: transplant (q=pred split prime) for 13..73: survival in cubes eps=1e-3:0.00025, 1e-4:0.07, 1e-5:0.47, 1e-6:0.915, 1e-7:1 (T=24); uniform: 0.
- line_cosets (b=(1,-1,2,0,1,-2,1,1)): lines through transplant keep length ~1.7e-5 at T=28 (gaps all within theta in (-1.6e-4,5e-5) of transplant);
  eta=1e-3 line: characters OK (margin 5.6e-5 at T=20,24; 3.1e-6 at T=28) but non-character windows cover (cover margin 4.4e-6) => G cap line empty.
  S_nc ~ 14.3-15.3; n_intervals 5e4-7e5.
- line_batch (T=22): 60 random (b,g): 9 pass characters, 0 have survivors; 30 lines through transplant: all survive (min len 1.6e-5).

## Quantitative full-torus bound (Cor 2.3)
mu(G) >= exp(-C_tau r^2 log P): refutes Poisson extrapolation exp(-c(Pi-1)) once Pi >> r^2 log P (window!).
So open point 3 (doubly exponential upper bound) is FALSE on the full torus at large Lambda.

## seeded base cosets (section 7)
- realisation cosets are Delta-translates of a component g_t + H_0 of H (base coset). Light non-character v: shift <v,eta> << slack
  => protected; heavy v (w >~ W/2): union bound, total mass e^{-W/4+O(r)}.
- Prop 7.1 (rigorous): nu(G cap C_eta) >= (4eps/pi)^m (1 - F_1(eps)) - sum_{heavy} omega. Positivity needs ~ C m^2 log P < W/4 - O(r): small m only.
- seeded_base.py: b=(1,..,1), H_0 seed arg(2+i); non-character survival 0.30 (=mu(G_1(13)), fibre alignment), 0.29 at eta=1e-3, 0.18 at 1e-2, 0 at 0.1.
  CAVEAT: toy has light characters (b^perp) which would fail; real profiles must have characters heavy or torsion.
- tension: Gaussian seeds need base-coset shifts in (pi/4)Z; light characters at such torsion fail; torsion denominators >=3 break Gaussian Liouville.
  grouped (2-blocks) profiles: K = span{e_i - e_j}; characters = representative profile (pinned, open point 5); seeds need Phi Gaussian-compatible.
