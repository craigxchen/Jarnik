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
(notes.md is shared with the theory agent; my notes live here from now on.)
N6 (elimination of a_j; reduced system). g_j = a_j + i W (a_j real deg<=8) has roots R_j. Let S_j = {r_T : T contains j, |T|>=2} (7 pts),
 Pi_j = prod_{S_j}(t - r). Then g_j = (lam t + mu) Pi_j, and Im((lam t+mu)Pi_j) = W.
 <=> W in span_R{Re Pi_j, Im Pi_j, Re tPi_j, Im tPi_j}  (codim 5 in R^9). 11 complex unknowns (22 real), 20 equations, gauge 1 => 1.
 Gauge fix: hyperbolic flow fixing +-1 moves r_1234 along flow lines; each crosses the imaginary axis once => Re r_1234 = 0.
 Singleton root r_j = -mu/lam.
 First LM run on full 34-var system (num_lm.py, 7 converged): all degenerate (Im r_T -> 0 or kappa -> 0 or coincident u_j).
N7 (pair structure, proved): g_j - g_k = a_j - a_k real, vanishes at r_T (T >= {j,k}) and conj => a_j - a_k = c_jk N_jk,
 N_jk = prod_{T>={j,k}} |t-r_T|^2 > 0 on R. => f_j = a_j/W never cross on R: the 5 circle points keep cyclic order for t in R.
 Cayley/cot form: f_ab = cot(theta_a - theta_b); condition (A_b - A_a) | (A_a A_b + W^2) for all pairs (tangent addition).
 S5-symmetric: S_{a|b} = prod_{A contains a, not b}(t - r_A) (deg 8), Im(e^{i phi_ab} S_{a|b}) in R W for all 10 pairs.
N8 (x-structured ansatz = twisted Thm-2 analogue): projection forgetting coord 4 pulled back from x-line (conic double cover):
 then it's a k=4 Thm 4.1 curve (rational branch x=inf, mu=4, B.C_x = 14), u_4 = H/Hbar, H = alpha + beta z,
 contact => beta real, Im alpha const. Count: per Thm4.1 curve codim 3 (Im-degree cond codim 1; square factor codim 2),
 Thm 4.1 moduli 3 => expected dim 0 (isolated). [gauges: H -> sH, (beta,q)->(beta/lam, lam^2 q)]
 Reduced LM (num_red.py, 2x1500 starts): converged ones all degenerate (coincident r_T, u_j = u_k).
N9 (clean s-form, proved equivalence): real Mobius s = (1+t)/(1-t) (UHP->UHP) sends branches +-1 -> inf, 0.
 Condition: Q_j(s) = prod_{T contains j}(s - rho_T) has Im Q_j = kappa_j s^4, i.e. e_k(R_j) real for k in {1,2,3,5,6,7,8}.
 <=> Im p_k(R_j)=0 (k=+-1,+-2,+-3) and prod R_j real. 30 real unknowns, 28 eqs, gauge s->lam s => dim 1.
 Single branch (toric 7) = e_1..e_7 real (gauge 2) => dim 0. k=4 two-branch: e_1,e_3,e_4 real: 14-9-1 = 4.
 Calibration (redk.py k=4, 200 starts): 194 converged, 57 nondegenerate. k=5 (num_red, num_defl): 0 nondegenerate so far.
N10 (unit-circle ansatz = all r_T in iQ = symmetry t->-t with inversion): 15 unknowns (y_T), 16 eqs e_{2l}(y on R_j) = C(4,l) m^l:
 overdetermined by 1. Real-arg constraint: sum_{R_j} Arg rho in {0, +-pi} exactly (arg Q_j stays in [0,pi] on R).
N11 numerics: LM in s-form converges (36/40) but ALL to degenerate (real roots / coincidences / kappa=0). Deflation kills convergence.
 Sign-pattern parametrisation (15 S4-classes of half-plane patterns, sum_{T∋j} eps_T = 0): pattern 0: 2/20 converge, degenerate.
 Next: symmetric ansaetze sigma = (g in {-s,1/s,-1/s}, pi in S5 involution, eps=+-1) -> smaller systems; compute expected dims.
N12 (ordering refinement): on (0,inf) u_j never passes 1, leaves on side sign A_j(0), returns on side sign lc(A_j) => prod R_j > 0,
 i.e. sum_{R_j} Arg rho = 0 exactly; ordering of 1/A_j(0) and of kappa_j coincide.
N13 (pair-structured reduction): (A): c_1 N_14 = c_2 N_24 = c_3 N_34 at rho_123 (2 real conds; c's real), rho_jk = extra root of
 c_j N_j4 - c_k N_k4. (B): Q_4 (with rho_4) real except s^4 (5 net), g_j = A_4 + sigma c_j N_j4 + i s^4 vanish at rho_123,
 rho_12, rho_13, rho_23 (8). 27 unknowns, 25 eqs, gauge 1 => 1.
 Symmetric ansaetze (symm.py): single involutions (g, pi, eps) with g in {-s,1/s,-1/s}; LM (symm_lm.py) on each.
N14 (S5-symmetric column form, proved): column polys C_a = prod_{A contains a}(s - rho_A) (deg 15, 30 places A).
 Condition <=> C_a - C_b = kappa_ab s^4 G_ab (G_ab = prod_{A>= {a,b}}, deg 7) <=> all C_a agree in coeffs s^14,s^13,s^12 and s^3..s^0
 <=> column moment vectors M_a = sum_{A contains a} (rho^k)_{k=+-1,+-2,+-3} equal and column products equal (5-column two-sided PTE).
 Cube/fiber form: Rbar_j is a full fibre of Z_j(s) = s^-4 prod_{R_j}(s - rho).
N15 structured sub-ansaetze: (X1) x = s+1/s (unramified at branches): Thm4.1 + conic coordinate, codim 3 -> 0.
 (X2) x = s^2 (ramified): 3 coords pulled back from k=4 two-branch mu=2 curve (4-dim) + u_4: 4 + 2 - 7 = -1.
 Complex Newton (cnewton.py, 21/30 converged): all degenerate (rho'_T = rho_S coincidences, kappa = 0).
 Symmetric ansaetze LM (symm_lm, 14 cases x 40 starts): many converged, all degenerate.
N16 contact-pair scan (mu_inf, mu_0): K = {1..mu_inf-1} U {9-mu_0..8}. kset/testK/fastK.
 (3,3): nondeg found directly. (4,0): found (14/30 fast). (4,1),(4,2): reached by staged continuation; (4,2) direct 2/10.
 (4,3): staged collapses to massive coincidences; (5,0): 0/40 direct+staged (slow solver). fastK_all.sh running all pairs, N=2000.
N17 Ptolemy identities (proved, from C_a - C_b = kappa s^4 G_ab): for each 4-subset abcd (e = 5th):
 k_ab k_cd N(rho_ab)N(rho_cd) - k_ac k_bd N(rho_ac)N(rho_bd) + k_ad k_bc N(rho_ad)N(rho_bc) = 0, N(rho)=|s-rho|^2,
 k_xy = v_x - v_y (Pluecker automatically). Only the 10 places of 2|3 cuts enter. Affine invariant (branch-free).
N18 CORRECTION: half-plane split 4/4 is forced only when Im(gamma Q_j) has constant sign: (4,4) [kappa s^4] and (8,0) [const].
 In general |n+ - n-| <= #sign changes of Im(gamma Q_j) (Cauchy index). Unrestricted patterns: (5,0) 29/400 nondeg (fastK).
N19 doubled-k=4 (dbl4.py): V-null (4,4) curves of class 28 on 4 points exist, non-pullback, local dim 7 = 6 + gauge (exact count).
 k=5 = doubled k=4 + choice of one place per pair (2^7) + rho_4, 5 net conditions (lift.py). Lifting+polish (96 starts):
 all polished solutions degenerate; dominant degeneration = perfect matching of 30 places into 15 coincident pairs with
 disjoint labels A,B (merged vector 1_{A u B}), e.g. (1,245),(2,145),(12,45),(13,2345)...
N20 (half-angle / star form, proved): z_a = (s-rho_a) prod_{x!=a}(s-rho_ax) (deg 5; 15 places = vertices+edges of K_5,
 one per cut, small-side label). C_a = z_a * Omegabar / zbar_a, Omega = prod of all 15. So u_a/u_b = (z_a zbar_b)/conj.
 V-null (4,4) <=> Im(g_ab z_a zbar_b) = k_ab s^4 N_ab (all pairs) <=> star sums Im p_k(Z_a), k=+-1,+-2,+-3, independent
 of a, and Arg-sums equal mod pi. All stars have the same # of UHP points.
N21 (Grassmannian / Moebius form): det(z_a,z_b) = k s^4 N_ab: Pluecker coords (w_a-w_b)N_ab(s) give a conic in Gr(2,5);
 phi: P^1 -> Mbar_{0,5} meets each boundary line twice (class -2K of dP5). Z(s) = M(s) Zhat(s), Zhat linear,
 det M = c s^4. Equivalently: 5 real Moebius maps zhat_a (all zhat_b^{-1} zhat_a elliptic, fixed pts rho_ab),
 a degree-4 map gamma=[g1:g2] with gamma(rho_ab) = zhat_a(rho_ab) (10 interpolation conds, one sign choice per pair),
 and Im(g1 g2bar) = c (s-p)^4 (s-q)^4. Count: 9 - 2 - 6 = 1. Vertex places = 5th coincidence of gamma and zhat_a.
 scan: (4,4) fastK 400 starts: 251 converged, 0 nondeg; (5,0) unrestricted: 29/400 nondeg.
N12 is WRONG (u_j may make a full turn on (0,inf)); only Sum Arg in {0,pi} (2.3) holds. Not used in search.md.
N22 (17:1x) xstruct.py (targeted x-structured ansatz, square 23x23) found 2 NONDEG real solutions (it 51, 76) among ~39 converged!
 Features: rho_123 = i (gauge), rho_4 = -q1 + i (?), Im n has degree 3 (x^4 coeff 0) so q2 = 1 in my normalisation.
 Square class of branch field: m ~ -1/sigma (sigma = c kappa^2 gauge-invariant). sigma ~ -0.535, -0.549.
 MUST VERIFY: construct conic, H, parametrise, check 30 places/vectors and contact 4 at both branches.
 Big fastK (4,4): 2000 starts, 1195 converged, 0 nondeg. split calib 2/400.
N23 !!! verify_x.py: xstruct solution 76 is GENUINE: 30 boundary places, all 30 rays once, all 15 cuts weight 2, contact 4 at
 both branches (w0 and inf), mpmath 60 digits (residual 2e-51). Solution 51 degenerate (16 places).
 (Earlier 44-place output was a bug in Hpoly; replaced by direct evaluation of H at the 16 points over the 8 x-values.)
 => real curves in the class EXIST (x-structured, X1). General searches missed them. Next: s-form conversion, Jacobian rank
 of general system there (local family dim), continuation of the family, algebraic identification (PSLQ), Q-points.
N24 s-form conversion (xs_to_s.py): 15 places satisfy (3.1) to 5e-51; general Jacobian rank 28 -> local dim 1 (mod gauge).
 identify_x.py (130 digits): coordinates NOT algebraic of deg<=12 (findpoly), and x-structured 23x23 Jacobian has one
 singular value 1.7e-81 => x-structured solutions form a 1-PARAMETER family (count 3.4 was off by one). So near this point
 the general family = the x-structured family. Symmetry: rho_{[3]-S} = -q1 - conj(rho_S), rho_4 = -q1 + i (rho_123 = i),
 q2 = 1 (my normalisation), Im n_4 = 0. Next: continue the family, fit implicit equations, look for Q-points.
N25 xfam.py: continuation gives 382 points (sigma in [-3.80,-0.51], q1 in [0.066,1.01]); Im n4 = 0 and q2 = 1 identically.
 xfam_sweep.py: sweep map (lambda, w) -> T_5(R) has rank 2 at 45/45 samples => the family sweeps a real SURFACE.
 Numerical algebraic relations between coordinate pairs of degree 5-8 (xfam_rel.py) -> model of the parameter curve
 not low degree. Symmetric structure (centered X): pairs p_r, -conj(p_r), r=0..3 (p_0=rho_123,p_j=rho_j), sum Im p_r = 0.
N26 EXACT STRUCTURE (xsym*.py): symmetric x-structured ansatz, centred X, p_r = a_r + i b_r (r=0..3), partners -conj(p_r).
 Thm 4.1 part <=> sum b_r = 0, a_r b_r = t (all r), t^2 = b0 b1 b2 b3. (II) automatic: Im n = c3 X^3 + c1 X.
 sigma = -prod |p_r|^2 (from L(0)=0), beta = X(X^2+e), q = X^2 + D (even). (III) <=> l(Y)/Y has a double root:
 ONE condition on (b1,b2) [b0=1, b3=-1-b1-b2] -> plane curve C. Q-points need t in Q (b1 b2 b3 a square).
N27 EXACT: (III) <=> disc_Y(l(Y)/Y) = 0; family lies on the degree-24 factor f24(b1,b2) (b0=1), which is S4-invariant:
 F(e2,e3,e4) = 32e2^6e4^3 - 96e2^5e3^2e4^2 + 37e2^4e3^4e4 - 256e2^4e4^4 - 4e2^3e3^6 + 1024e2^3e3^2e4^3
   - 992e2^2e3^4e4^2 + 288e2e3^6e4 - 27e3^8 - 256e3^4e4^3   (e_k elementary symm. of b0..b3, e1 = 0).
 u = e4/e2^2, v = e3^2/e2^3: quintic Phi(u,v) with triple point at origin.
N28 RATIONALITY: lines v = k u through the triple point: disc_u = (k-4)^5 (9k-4)^3 => curve Phi=0 birational to the conic
 w^2 = (k-4)(9k-4), which has the rational point (4,0): k = 4(m^2-1)/(m^2-9), sqrt = 32m/(m^2-9). So the INVARIANT curve is
 RATIONAL over Q. Q-curves need: quartic W^4 + v W^2 - v^2 W + u v^2 splits over Q, and u in Q^2 (t^2 = e4).
N29 u(m) = -(m+3)(m^3-11m^2+39m-13)/(16(m-3)^2(m-1)^2), v(m) = -(m+1)(m^3-11m^2+39m-13)/(4(m-3)^3(m-1)) (other branch m->-m).
 b-quartic (W = vZ): R = A Z^4 + B (Z^2 - Z) + D, A = m^5-9m^4+18m^3+54m^2+13m-13, B = -4m^5+36m^4-104m^3+72m^2+108m-108,
 D = -m^5+9m^4-18m^3-54m^2+243m-243; irreducible over Q(m); disc = -2^15 (m-3)^12 (m+1)^6 (m+3)(cubic)(m^4+70m^2-48m+9);
 resolvent cubic irreducible => Gal S_4. 8 simple branch points => S_4-closure genus >= 25 => finitely many m with all
 b's rational (Faltings) => FINITELY MANY Q-curves in this family. (Need: geometric irreducibility check; construct-from-m test.)
N30 construct_m.py: curves built from rational parameter m (m=-2,-1.8,-2.5): 30 places, every ray once, contact 4 at both
 branches (mpmath 50 digits). The real family is m in an interval of the branch with u>0, v<0, 4 real b's.
 Geometric monodromy: normal subgroup of S_4 containing a transposition (simple zero of disc, A != 0 there) => S_4.
 Riemann-Hurwitz for S_4-closure: >= 8 transposition branch points => 2g-2 >= -48+96 => g >= 25 => finite (Faltings).
N31 FINITENESS (rigorous): Y = b^2/(-e2) (scale & sign invariant) satisfies G = (Y^2 - Y + u)^2 + v Y = 0, monic-ish (LC = 256(m-3)^9(m-1)^7),
 disc_Y = c (m-3)^42 (m-1)^34 (m+1)^2 (m+3) cub^3 qua: 8 points (m=-3, 3 roots of cub, 4 roots of qua) with ODD order => odd local
 monodromy (e>=2). Galois closure (group N != 1): 2g-2 = |N|(-2 + sum(1-1/e)) >= |N|(-2+4) => g >= |N|+1 >= 3. Faltings =>
 finitely many m in Q with G split => finitely many Q-curves in the component. Exact search |p|,|q|<=60: no m with even one rational b-root.
 NB A = cub*(m+1)^2 (Res(A,cub)=0) so the earlier Z-form analysis was replaced by the Y-form.
N32 FINAL: search.md rewritten (v1 backup search_v1_backup.md). Status anisotropy_gain: Theorem A (real family M, m in (-3,-1),
 sweeps a surface), Theorem B (finitely many Q-members; Faltings), blind-search 'emptiness evidence' withdrawn.
