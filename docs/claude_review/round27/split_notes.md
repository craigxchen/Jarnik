# split lane notes (round 27)

## N1. Literature (2026-10-29)
- formal-conjectures 887.lean (raw fetch): parts.i, parts.ii research open; rosenfeld_infinite solved; rosenfeld_4 solved (negated)
  with formal proof link mnn31/erdos-lean Erdos887.lean#L538. 886.lean: erdos_886 open; rosenfeld_infinite, rosenfeld_bound solved.
- mnn31/erdos-lean ErdosLean/Erdos887.lean (downloaded to split_checks/dl/, data): theorem five_divisors_infinite with C=7/2:
  n = u(u+2)(u+3)(u+5), (u+2)(u+3) = 3(w+1)(w+2) [(2u+5)^2 - 3(2w+3)^2 = -2], start (1843,1064), divisors
  (u+2)(u+3), 3(w+1)(w+3), (u+2)(u+5), (u+3)(u+5), 3(w+2)(w+3). Also card_window_le_four for C=2 (at most 4 divisors when C=2).
- czech72.py: along 40 Pell points, log(dist^{7/2} H) periodic in {4.064, 4.267, 4.96}; gamma -> 7/2; lcm/n -> 7/8; C_act -> 2 sqrt 3.
  Torus profile: 8 rational places of the conic, P1=(0,0) constant; others give B = {25},{34},{1}^c,{123},{135},{124},{145}:
  7 cuts once each, B.C = 7, two real-quadratic branches (u/w -> +-sqrt3) of order mu_0=1. MR Thm 2.16: alpha = 7/(2*1) = 7/2.
  => UNBALANCED A-curve at alpha = 7/2 < 15/4 (but > tau = 10/3). Must lie in Z of split A_5(delta) for every delta <= 1/14.

## N2. Structure facts (derived here)
- Near-tight balanced split: x_B = log n_B/log X, x_B + x_{B^c} = 1, P_ab := sum_{B ni a, not b} x_B; Liouville tight iff all P_ab = 4.
  Symmetric solution x = 1/2 (Gaussian-like); W has x_B in {0,1} (one side trivial on each cut). Solution space of
  antisymmetric y with all P_ab(y)=0: >= 5-dim.
- T = G_m^5/G_m is an open subset of M_{0,7} (points 0, inf, u_1..u_5). K_{Mbar_{0,7}} = -(1/3) sum_{|T|=2} D_T.
  Vojta with D = D_{[5]} (all light points collide) at bounded residues: h_K + lambda = (2/3)lambda_inf - (1/3) sum h_{D_T}
  = [(5/2)(1+2 delta) - 5/3] log X > 0 for balanced bounded-residue points: predicts Z even at delta > -1/6.

## N3. C=7/2 family = quadric cone (fam72.py, exact)
- Family: conic u(u+a+b) = k w(w+c+d), k = ab/(cd); D1=(u+a)(u+b), D2=(u+a)(u+a+b), D3=(u+b)(u+a+b), D4=k(w+c)(w+c+d), D5=k(w+d)(w+c+d).
  Second factorisation (u+a)(u+b) = k(w+c)(w+d) holds identically on the conic. Jacobian rank 3 at 4 random rational points (exact).
- With y_j = D_j/D_1: Y2=a/(u+b), Y3=b/(u+a), Y4=c/(w+d), Y5=d/(w+c) (Y=y-1) => Y2*Y3 = Y4*Y5 (from (u+a)(u+b)=k(w+c)(w+d), k=ab/cd).
  So W_{7/2} = quadric cone {(u2-u1)(u3-u1) = (u4-u1)(u5-u1)} (vertex e). S_5-images: 15 such cones (base index x pairing).
- Twisted Galois (u -> 1/conj u) does not preserve the cone: not defined over Q in T_5.

## N4. Section exponent depends on split asymmetry (key new idea)
- For a profile x (30 weights, near-e: X_a = sum_{B ni a} x_B common), a Laurent poly f with ord_e = mu and toric pole orders pi_B
  gives: A-points at gamma < 15 mu / sum_B x_B pi_B (log H = sum x = 15) lie on {f=0}.
- Binomial chi^m - 1: cost c_x(m) = sum_B x_B (m(B))^+ ; pairs give P_ba = 4 on Pi_nt (constant).
- 1/gamma_sec(x) is concave (inf of linear) and invariant under x_B <-> x_{B^c}; hence symmetric (twisted-type) profiles are the
  hardest for sections. At one-sided realizable vertex x_W, binomials cannot beat 15/4 (else W in a subtorus; contradiction).
- Integrality at the root arrangement (bounded residues) with poles along Delta-bar: f = g/V^m; V (Vandermonde) has ratio exactly 15/4
  at x=1/2 => capped (mediant argument).

## N5. Section cost on Pi_nt (profiles.py, search_beat.py, qf_reg.py, cone_section.py)
- dim Pi_nt = 11 (not 5): cut balance gives P_ab + P_ba = 8, so P_ab = 4 <=> X_a = X_b; only 4 extra conditions.
- IDENTITY: for x in Pi_nt and support S, cost_x(S) = sum_cuts pi_{B^c}(S) - 2 sum_cuts x_B mid_B(S), mid_B = midpoint of <S,1_B>.
  If S is centrally symmetric (mid_B(S) = c(B) linear) then cost_x(S) is CONSTANT on Pi_nt (uses sum_B x_B c(B) = X_c sum c_a = 0).
  => binomials, zonotopes, Vandermonde products: same cost at every near-tight profile => capped at 15/4 by Theorem B (x=1/2).
  Binomial cost c(m) = (1/4) sum_{30 B} |m(B)| >= 4 (pair B with B xor {a}).
- Random search: 8000+ point sets S, LP min over Pi_nt of cost vs 4 reg(S): max ratio exactly 1.0000 (no beat).
- One-sided free polytopes Q_F(n): reg/n = 3, 3.5, 3.33, 3.5 (W class), <= 3.33 (1^5 3^10) < 15/4. Volumes 4!vol in [89, 118] < 197.75.
- Quadric-cone section G: min cost over Pi_nt = 8.333 > 8 (no beat).
- => Conjecture WTB: 4 reg(S) <= sum_B x_B pi_B(S) for all S, x in Pi_nt. Proved for centrally symmetric S and at x=1/2.

## N6. Phase-0 reformulation of #887 (derived)
- divisors d_j in [sqrt n, sqrt n + C n^{1/4}] <-> s_j = d_j + n/d_j, z_j = d_j - n/d_j, s_j^2 - z_j^2 = 4n, k_j = s_j - s_1 in [0, C^2+1].
  => z_j^2 - z_1^2 = 2 k_j s_1 + k_j^2 (j=2..M). For fixed k-tuple: CI of (M-1) quadrics in P^{M+1}:
  M=3: dP4; M=4: K3; M=5: K = O(1) general type; M=6: K=O(2).
  m-czech: k = (0,224,225,324,900) (offset^2), mnn31: k=(0,3,4,9,12). Their families are Pell conics on S_k.
  #887 parts.ii with K <=> (for each fixed C finitely many k) each S_k with K+1 entries has finitely many integral points.

## N7. W-vertex decision (Prop 3.3)
- W code profile F_R = {B : |B cap {2,3,4}| even} (labels F_0..F_4). Exactly 278 one-sided near-tight profiles (onesided_enum.py).
- WTB at x_R holds for all f not divisible by g_W (W-curves, mult 4, E_B.C = [B in F_R]).
- deg W > 14 (w_equation.py, F_1000003). Tangent cone contains cone over irreducible octic Sigma_W; branches are cusps
  Y = xi s^2 + e3 (1,1,1,1) s^3 + ... (F_1..F_4 even in s, F_0 odd part -e3 x) => two sheets, ord_e g_W >= 16 (derived).
- iota(W) != W: germ (xi, eta) -> (xi, -i eta); numerically membership residual 0.14 / 0.05 vs 1e-15 controls (w_member2.py).
  => g_W satisfies the cap at the complementary vertex iota x_R. Inversion trick would have closed Q-split-1 if iota(W) = W.
- Numerical counts (w_count.py, random Newton, lower bounds): line -> deg W, linee -> deg W - mult_e W, curve -> cost_{x_R}(g_W).
- RESOLVED (numerically): splitting the cusps of a W-curve gives the expected 4-dim family of profile-x_R curves through
  e with 4 smooth branches (solution-set dim 6 = 4 + 2 affine), dominant (8/8 random targets; W-family 0/8).
  Prop 3.7 => cap at x_R for all f. Code vertices (20) span 10-dim hull containing 1/2 => Cor 3.8 (WTB on Q_code).
  Random-start solver fails everywhere (also at W-class): need start curves at other vertices (e.g. even-set x_[5]).
