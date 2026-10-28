# notes for search.md (round 24)

N1 (reduction). Curve C ~ P^1_Q over Q, twisted map u_j = c_j P_j/Pbar_j (Hilbert 90 over Q(i)(t)/Q(t)).
 Balanced B.C=30: E_A.C = 1 for all 30 rays (parity => E_A.C = E_{A^c}.C, cut weight 2).
 If 30 distinct places with ray vectors: Gal(Qbar/Q(i)) preserves vectors => each place Q(i)-rational:
 r_T in P^1(Q(i)), T subset [4] nonempty (vector 1_T), conj place rbar_T (vector 1_{T^c}).
 P_j = prod_{T contains j} (t - r_T), deg 8. Degenerate (chain vectors) = coincidences r_A = r_B.
N2 (branch condition). Branch pair t0, t0' = +-sqrt(m), m in Q nonsquare (PGL2(Q) normalisation).
 ord_{t0}(u_j - 1) >= 4 for all j  <=>  gamma_j P_j - conj(gamma_j P_j) = 0 mod (t^2-m)^4 (c_j = gamma_j/gammabar_j)
 <=> Im(gamma_j P_j) = kappa_j (t^2-m)^4, kappa_j = Im gamma_j (coefficientwise Im). kappa_j != 0 else u_j == 1.
 Cayley form: u_j = (f_j + i)/(f_j - i) with f_j = Re(G_j)/(kappa_j W) real rational, poles only at +-sqrt m order 4.
N3 (count, real). unknowns 30 (r_T) + 4 (phases) = 34; gauge: stabiliser of pair in PGL2(R): 1-dim.
 equations 4 x 8 = 32. expected dim of curve family = 34 - 1 - 32 = 1.  (single rational branch, toric 7: 34-2-32 = 0.)
 split Thm 2 check of method: 15 + 4 - 1 - 16 = 2 = actual dim. k=4 twisted 2-branch: 14+3-1-12 = 4.
N4 alpha: dist ~ |t-t0|^4, H ~ H(t)^30, |t - t0| >> H(t)^-2 => gamma = 30/8 = 15/4.
N5 argument principle: wlog kappa_j>0, arg G_j in [0,pi] on R => each R_j has 4 roots in each half plane.
## plan (after reading expert note s3, r23 README, split s1,s4-6, toric s1,s3,s6)
Verdict expected: reduces_to_split_cap for all 4 candidates. Key facts:
- X_{Q(i)} = LM_5; Q-structure sigma_X = iota o conj (iota: x->x^{-1}). y in X(Q) iff conj(y) = iota(y).
- point parity: v_{bar pi}(y) = -v_pi(y); no contact at inert q or 2 (T_5(Q_q) compact).
- Binary-form model (Hilbert 90): T^Q = (P^1 \ {Q=0})^5 / SO(Q); split Q=ab, twisted Q=a^2+b^2.
  x_a/x_b - 1 = CR - 1 ~ det(P_a,P_b)/...; |Q(P_j)| = prod p^{|eps_j(p)|} (both forms, up to 2).
  H = prod_p p^{range_j eps_j(p)}; dist ~ |det|/sqrt(Q Q).
- C1 isotypic: iota-even/odd; orders: ord f = min(ord f+, ord f-); t_par = t_max, t_other = t_max-1 by explicit
  products V, V_abc (odd, ord3, H0(B)), D_ab=x_a/x_b+x_b/x_a-2 (even, ord2, H0(B)). sigma-semilinear: V^sigma, i V^sigma same orders.
  Reality: Q-rational f => f(y) in Q; forced divisibility already symmetric (width). no gain.
- C2 conjugate: (y, ybar)=(y, iota y) graph pullback -> single section of (m1+m2)B. Thue on X^n at rational target:
  weighted index additive (filtered tensor product), forced div. of product polytope <= sum of projections.
- C3 parity of points: built into Gaussian profiles => fake-satisfied (toric s6.3: balanced 5-row fakes at 7/15).
  Pair: split min(X,Y) vs twisted (X+Y)/2 -> parity weakens Liouville.
- C4 Hilbert 90: det divisibility = pair gcd; Pluecker = Ptolemy (abc vacuous); bihomogeneous => K1.
- general: Pic rank 15 (twisted) vs 26 split: only log powers in counting.
## C4 derivation (Hilbert 90)
x_j = P_j/Pbar_j, P_j=a_j+ib_j primitive. [P,Q]:=Im(Pbar Q)=a_P b_Q - b_P a_Q. x_a/x_b - 1 = -2i[P_b,P_a]... |x_a/x_b-1| = 2|[P_a,P_b]|/(|P_a||P_b|).
eps_j(p)=v_pi(P_j)-v_pibar(P_j); N(P_j)=2^{e_j} prod p^{|eps_j|}; H_B = prod p^{range eps}.
v_p([P_a,P_b]) >= c_ab = max(0,min(eps_a,eps_b),min(-eps_a,-eps_b)); c - (|ea|+|eb|)/2 = -|ea-eb|/2 => dist_ab >= c e^{-d_ab/2}: pair bound exactly.
KEY: G(P,Pbar) bihomog deg_j in (P_j,Pbar_j) => G = prod Pbar_j^{deg_j} * g(x), g(x)=sum coef x^c, c_j in[0,deg_j].
 forced div of G at pi,pibar == toric forced div of g. => Hilbert-90 polys == Laurent polys (TL).
 reality: Gbar=±G identity <=> Newton set symmetric c<->d <=> f_pi=f_pibar => max(f_pi,f_pibar) no gain.
 NOTE: MW fails for half-integer sets (S={0,1/2}: reg 1 > MW 1/2) but the P-dictionary lands in integer exponents c, so no issue.
## C2: tensor filtration -> index additive; forced div of product polytope <= sum of projections.
## isotypic: exact t_+/t_- m<=4 match prediction (t_par(tmax)=tmax, other=tmax-1). sections.py ok.
## C2 cap subtlety: width_{Q_F}(v_1..v_n) <= sum width_{Q_i}(v_i) (wrong direction for cap). Fix: cap on data with
 pairwise disjoint prime supports (then equality) + simple (squarefree) eta-balanced profiles. F=g(x/x') = (TL) at y1/y2.
 (y, iota y): Delta_iota^* F single section of (d1+d2)B. Note width sublinear on chains -> use simple profiles for caps.
## C3: twisted RHS_f(v) = mean(split RHS_f(v), split RHS_f(-v)) -> twisted feasible ⊇ split feas(v) ∩ split feas(-v).
## t_ab = [P_a,P_b] / (prod p^{c_ab} 2^{(e_a+e_b)/2}) up to sign -> check in h90.py
