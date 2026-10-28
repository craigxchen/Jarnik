# intra notes (parking derivations)

## Key observation
All y_S = pr_S(Z), Z = [z] in T_M(Q). pr_S : T_M -> T_5 coordinate projection (homomorphism over Q;
commutes with conj and inversion). Phi = (pr_S1..pr_Sn): T_M -> T_5^n. G = Phi^*F.
- G(Z) = F(y).  Q_G subset Phi^*Q_F, Phi^*(s_1..s_n) = sum iota_{S_i}(s_i).
- <Phi^* s, v_w(Z)> = sum <s_i, v_w(y_Si)>  (v_w(y_S) = v_w(Z)|_S mod 1).  => widths equal.
- ord_e G >= min|J| (Taylor of F composed with linear map in log coords).
- => (*) for F implied by (TL_M) for G at Z.  = route 5's class (TL on T_n at cluster point, poly chosen after profile)
- route 5 slab monomials x^k, k' in {-1,0,1}^d, k_n = -sum: anchored star blocks S_j = {n} u 4 indices,
  s_j = sum k_x (e_x - e_n) in 4R_5 (sum|s| <= 8). So F_c = Phi^*F with block degree 4, ceil(d/4) blocks.

## Balanced cap (T_M MW)
Theorem B (threshold/construct.md): reg(S) <= 2^{1-n} sum_{P subset [n]} w_P(S) for S in Z^n.
On sum-zero lattice: w_P = w_{P^c}, so reg <= 2^{2-M} sum_cuts width(1_A).
Simple eta-balanced profile on all 2^{M-1}-1 cuts: RHS_F = 1/2 sum_p width_{Phi^*Q_F}(1_{A_p}) log p
 >= 1/2 (1-eta) w sum_cuts width_{supp G}(1_A) >= 2^{M-3}(1-eta) w ord_e G.
LHS(*) <= max_i lambda_i * min|J| <= max_i lambda_i * ord_e G.
=> (*) holds whenever max_i lambda(y_Si) <= 2^{M-3}(1-eta) w.
Subtuple scales: w_S = 2^{M-5} w (1 +- eta) per S-cut; so condition = mu(y_Si) <= 4(1-eta)/(1+eta) roughly.
Cluster: lambda(Z) ~ W/4 = (2^{M-1}-1)(1 +- eta) w /4. Need (2^{M-1}-1)(1+eta) <= 2^{M-1}(1-eta) - ...: eta < 1/(2^M-1).
mu(y_S) >= lambda(Z)/w_S ~ 4 - 2^{3-M} - O(eta).
This is Cor C3 (cube profile) transported. Cube profile needs r >= 2^{M-1} primes: outside route-5 regime d >= 316 ln 4r.

## TL normalization check
h_B(y) = 1/2 sum_cuts width_{R_5}(1_A) log n_A = sum log n_A since width_{R5}(1_A)=2. OK.

## Route 5 identification (round10/determinant.md §1)
G_k(z) = N^{(D-|k|_1)/2} prod z^{k+} conj z^{k-} = N^{D/2} Z^k ; F_c(z) = N^{D/2} G(Z), G = sum c_k x^k on T_n.
(I1) Norm F_c(z) >= N^D prod p^{-h_supp(a_p)}  <=>  log|G(Z)| >= -1/2 sum_p width log p  = TL_M exactly.
Star partition: blocks S_j = {n} u B_j, B_j partition [n-1] in 4-sets => Phi^*: M_5^q -> M_n iso (T_n ~= T_5^q).

## Lacunary dichotomy
width_{Q_F}(v) >= width_{pr_1 Q_F}(v_1) - sum_{i>=2} d_i width_{R5}(v_i) => RHS_F >= RHS_{pr1}(y_1) - sum_{i>=2} d_i h_B(y_i).
LHS <= R'_1 lambda_1 + sum_{i>=2} floor(10 d_i/3) lambda_i  (Prop 4.2, R'_i = reg(pr_i supp F) <= t_max(d_i)).
Roth lemma needs d_{i+1} <= omega d_i and d_i h(y_i) >~ d_1 h(y_1): heights of subtuples in [W/2, W] (Lemma 1.2) => fails.

## Gcd/integrality of reduced tuples
monomial m bihomog (a,b): m(z) = g^a gbar^b m(zeta), |m(zeta)|_pi <= 1 => toric bound at pi <= integrality bound. Toric TL dominates.

## Template-dependence of base locus
B(V_H) depends on profile, H, m' => Lemma 3.3 counting c(Z) not uniform. Owner constraint.

## Plan of computations
A. exact reduction identities on actual clusters (random F, exact Q(i) values, widths).
B. route-5 slab monomials as pullbacks (explicit).
C. balanced-cap numeric check; sparse-profile supercritical example and exact relation check G(Z)=0.
D. fakes: same profile real angles; G(Z_fake) != 0.
