# kill_dp2 working notes
## Plan (killer, own code only)
1 geometry of S: smooth (Q,F3), T, F3, F4', Res, cubic field, disc, no bitangent through p
2 Frobenius traces at 3 n=1..8 -> char poly on NS
3 net |-3K-4P|: square tangent cone members
4 Gamma_3 global section unique; D-points; Gamma_3 inequality (T5)
5 S_0: smooth, P0 off lines, C_0 points
6 usearch reproduction
7 toric orders: L'=2(-K)-B14 on L5 (k<=3, maybe 4), -K calibration, k=4 control pair on L4
8 literature (WebSearch)
## Desk checks done (in head, recorded)
- Theorem A cone step: NE(St)∩F^perp = R+F via Zariski decomposition argument (no vertical neg curves) OK.
- tau=k+1/2 case OK. Projection = Galois average OK (avg of isometries is orthogonal proj onto invariants).
- beta=17/18 recomputed by hand OK. eps=4/3 decomposition OK.
- disc(19t^3-10t^2+20t+8) = -1707008 = -1667*1024 OK, 1667 prime.
- Coxeter element of E7 (exponents 1,5,7,9,11,13,17) gives traces t1..t8 = 0,2,3,2,0,-1,0,2: matches author's char poly.
- Theorem C: K_Y+Dt+E = L-E, (L-E)^2=1 nef big OK. Y\(Dt+E) = S\D. Bounded u -> Sigma_U-integral OK. Z_K depends on K (non-uniform).
- Y_U and double cover W_U of S branched on D are same K3 (zeta = w'); boundary 5 curves match D4 + ramification.
- (UE) for m=1 is a restatement (tautological) of the window question, not a reduction.
## Results so far (own code)
- k_geom.out: smooth over Q,F3,F5; T,F3,F4' as claimed; Res=-80353/1024; vertical dir F3=19/8; cubic field disc -2^10*1667, 1 real root; tangent-line disc t^2*(irred deg10). ALL CONFIRMED.
- k_frob.out: traces 0,2,3,2,0,-1,0,2 ; char poly (X-1)(X+1)Phi18 (Coxeter class of W(E7)). Pic(S)=ZK CONFIRMED.
- k_deg6.out: net dim 3, contains xD,yD; square-tangent-cone members = intersection of two conics, eliminant irreducible quartic. CONFIRMED (no Q-member).
- k_gamma3.out: dims (2,3,0)=1,(6,8,3)=1,(6,9,3)=0,(6,8,4)=0,(6,3,8)=1 => Gamma_2,Gamma_3,Gamma_-2 unique. D-points x+nu+ -> 3/4 from below (0.664..0.733). CONFIRMED.
- by hand: x_sec(0,nu-)=2/3 iff nu-<=2/9 (min_k (k-1)/(3k) at k=3; k=-2 needs nu-<=1/2; fibre nu-<=1/3). CONFIRMED.
- k_s0.out: S_0 smooth, Res=1, C_0 points u=1, X^1.5 H -> 1 CONFIRMED. BUT f3 = (x+y)(3x^2-xy+2y^2): D_0 has a RATIONAL branch at P_0 (alpha_{D_0}=4), unlike S. S_0 special in several ways (Q-lines, hyperflex line z=0, cusp = P_0 on parabolic curve of cubic surface, rational tangent of D_0). "agrees with S in every fibre-local statistic" overclaimed (minor).
- Prop A3 = MR-L Thm 4.5 (split cubic surface, tangent-plane section cuspidal => alpha=3/2) transported to dP2 via contracting sigma(E_0). Mechanism known; instance new.
## Transfer: one-line proof for ALL k and all multiples (closes author's open point 4)
Kapranov: L_k = iterated blowup of P^{k-1} along coordinate strata; singleton rays = hyperplanes zeta_i=0.
General line l through e avoids all codim>=2 coordinate strata => l~.E_A = 0 for |A|>=2, =1 for |A|=1.
Z_ab = strict transform of hyperplane {zeta_a=zeta_b} (contains e) => Z_ab.l~ = 1, (Z_ab)~.l~ = 0 on Bl_e.
(K_{Bl} + sum Zt_ab + E).l~ = (-k + (k-2)) + 0 + 1 = -1 < 0; lines through e cover => not pseudo-effective, every k>=3.
At k=5: L'.l~ = sum_{|A|=1} d_A = 5 => tau_e(L') <= 5; polynomial proof: F in H^0(tL') restricted to zeta=e+tv times prod(1+tv_j)^t is a poly of degree <= 5t, so ord_e F <= 5t. With order-5 section at t=1: tau_e(L') = 5 EXACTLY.
k=4 control (T2): same -1 => no Lang-Vojta prediction at k=4 => consistent with gauss Thm 4.1 (dense near-tight 4-tuples). PASS.
Also: lines through e minus one point are A^1's in the complement => complement log-uniruled.
- k_toric3.out: -K 3,6,10; L' 5,10,15; K+sum Z (orientation rep) 5,10,15; k=4: -K 3,6,9; K+sum Z 2,4,6. Line bound tight at k=4,5.
- k_usearch.out (after fixing my own scaling bug): exactly the author's 4 points. CONFIRMED.
- INCIDENT: k_usearch.py was first written to /home/user/Jarnik (backgrounded cd); moved out immediately; git status clean.
## Minor issues found
- 6.3 "genus <= 9 ... Faltings": needs genus >= 2 (rough RH: 10 I_1 + one II gives R=12, g=3; unverified).
- 6.1 Ru-Vojta threshold: for L=pi^*(-K) near-tight h_E=(2/3)h_L so need beta>3/2, not >1. Conclusion unchanged.
- T1 analogue overclaims sameness of S_0 (D_0 rational branch, f1=f2=0, Q-lines).
- (UE) for m=1 is a verbatim restatement of the window question.
- 8.1 family claim on X is a sketch stated as fact.
- Theorem A restricted set non-vacuity not shown; heuristic count C^{1-2x-nu} < 1 for x>1/2: no richer than generic.
