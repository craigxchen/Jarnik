# notes for lift.md (round 18)
## plan
A. counts: curve level (genus 1 / any genus, beta_7, 2beta_7, beta_8), integer level (lattice index).
B. Theorem L (rigorous): no class-beta_7 curve (any genus) in Mbar07 projects onto E_t (t generic).
   Gamma in |L (x) O(2)| on E_t x P^1 through 25 points (B_A, v_A), v_A in {x1,x2,x3,x4}.
   Lagrange in v at x1,x2,x3: sigma_j in H^0(L(-sum B_j)) (1-dim, tau_j); need rank [tau_j(B_{4,k})] (4x3) <= 2.
   Model: E_t = {u^2 = kb (a^2 - s h^2), w^2 = kc (a^2 - t' h^2)} in P^3; b = (-A1+a u)/(2A2), c = (-A1c + w)/(2A2c).
   L_q = 3H - r1 - r2 - r3 - q.  Spurious q: r_i, B_{j,k} (j<=3) -> handle directly with double-vanishing.
   Specialisation: rank-drop locus Z in Pic^8 closed, proper over Z_(p) -> F_p-computation suffices.
C. Corollary: beta_8 lift of E_t => both 7-projections birational of class 2 beta_7, degree 4 over E_t.
D. integer lattice test on Phi_N: (u,v) with u beta_j = alpha_j v mod E_j, need |u|,|v| <= e^{8w}; index e^{25w}.
## counts derived
- one point (k=6->7), genus 1 beta_7 over fixed E_t: x6 deg 8 on C' (deg 2 unram. cover): 16 params - 25 = -9; family 3-9=-6.
- any genus beta_7: Gamma: 24 - 25 = -1 per curve (linear!).
- 2beta_7 genus1 over E_t: 32 - 50 = -18 (Gamma: 80-75 nodes =5, minus 23 extra nodes = -18).
- beta_8 genus1 over E_t: 2*(-18) - 6 (D_{i67}) = -42.
- integer: one point per 6-config: box e^{16w} vs index e^{25w}: -9 (L7 sec8). two points: -42 (with 3 inflation pieces).
- E_t's own lattice: mover a: index e^{15 w6} vs need e^{8 w6}: -7 w6 -> integer heuristics fail on curves (structure beats count).
## KEY FINDING (beta7lift.py p=1009,10007, t0)
- only L with a Gamma: q = B_bc in old param (r_j = S=empty pts over y_j) -> nullity 2 = the PENCIL:
  Gamma_lam = {P_lam(v) - T Q_lam(v) = 0}, L = T^*O(1) = 2H. every member passes through all (B, x_j).
  => beta_7 lifts of E_t EXIST: C'_lam = {f_a(a)=f_b(b)=f_c(c)=f_lam(d)}, genus 9 (fibre product, 6 branch pts, (Z/2)^4),
     genus 5 if lam in {lam_s', lam_t'} (shares branch value s or t'). = item 112 mechanism. 56 contacts all over y_j.
  So Theorem: beta_7 lifts over E_t = pencil lifts only (genus >= 5).  [need openness-in-t argument]
- every rational d over Phi lies on pencil member lam(d) = -(d^2-T)/(P1(d)-T Q1(d)).
- beta_8 pencil lift: 5 movers, (Z/2)^5 cover, 8 branch pts: genus 33.
- genus<=1 beta_8 lift => both 7-projections birational, class 2beta_7 (else degree 2 onto genus>=5 pencil curve: RH impossible).
## openness argument plan
 R_i = div(m_i) - forced (r_i:3, B_{j<=3}:1, B_4 row: 1 for minors containing row, q*: 2). gcd of pushforwards = 1 at t0 => open nbhd.
## DONE
- beta7lift2.py over Q: 9 t's PASS (gcd rho = 1, specials 0 except q*=B_bc nullity 2). coplanar.py symbolic. branch_st.py 9 t's.
- intlift.py N=3..8: unsplittable 6..16 contacts; max skew 0.05w at N=8; mover skew ~7w.
- lift.md written (sections 0-7).
