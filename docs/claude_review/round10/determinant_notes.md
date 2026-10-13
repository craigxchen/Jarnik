# Notes for determinant.md (working file)

## Key realisation
construct.md's barrier (C3) is the *uniform* cube measure on {+-1}^M, i.e. 2^M columns.  For a FIXED F a
random r-column profile approximates the cube; but F may be chosen AFTER the profile (per cluster).  The
restricted LP  max_A ord / sum_p w_p h_A(beta_p)  is then >= c sqrt(n/log 4r): slab-box support.

## Slab-box parameters
k' in [-T,T]^{n-1}, k_n = -sum k'.  <k,beta> = sum_{j<n} k_j (beta_j - beta_n), |coef| <= 2.
Hoeffding: P(|X| >= s) <= 2 exp(-s^2/(8 T^2 (n-1))).  s = w := T sqrt(8 (n-1) ln(4r)) gives 1/(2r).
Union => |A| >= (2T+1)^{n-1}/2, h_A(beta_p) <= 2w.
E(m) = binom(m+n-2, n-1) = dim Poly_{<m} in n-1 vars.  Need U = |A| >= 4E for the counting Siegel.
Moments bounded by B = T^{m-1}.

Counting Siegel: (2H+1)^d >= (H+1)^U / (2UHB+1)^E, d = dim span of kernel vectors with |c|_inf <= H.

## Liouville/divisibility
F(z) != 0  =>  (m/2 - w) W <= 2 ln|c|_1 + 2 m ln(D C/2) - 2 ln m!,  D = 2(n-1)T.

## Hadamard cap
eps = (1/n) sum_p <eps,h_p> h_p ;  E|S_n| <= sqrt n  => Omega(A) <= sqrt(n) * (1/n) sum_p h_A(h_p).

## Rank lemma (squarefree)
cut Gram G_xy = 1 - 2 d_xy/W' <= eta = 4 log^+ C'/W', G_xx = 1. rank >= n^2/||G||_F^2 >= n/(2+2n eta).

## Obstructions to nonvanishing
* box/zonotope subsets => product of binomials => character inequalities (critical).
* level-sum lemma: if h_{supp c}(a_p) < m then all level sums of c in direction a_p vanish
  => no unique-minimiser (Newton) certificate at p.  With residues: face sum sum c_k rho^k mod p.
* grid lemma tuples: E_t h_A(beta_p o t) = E_{x ~ mu_p^{(x)M}} h_A(x) = Omega(A) for balanced sign columns.
* low rank: T_S subset Z_A if reg(B_S A) < m; B_S A in box of sides e_p w.
* sparse delta (|supp| < m/(2T)) in Z_A.
