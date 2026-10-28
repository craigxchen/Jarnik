# notes round 26 (ze)
## record facts
- [D]=round10/determinant.md. Thm 2.1 slab box (T=n), Thm 4.1 rigidity: F_c(z)=0 all c in V_H, i.e. sum_{k in A} c_k w^k = 0, w=e^{i delta} in T^{n-1}; forms P_c(w) have support A (Laurent), |A|>=(2T+1)^{n-1}/2, dim V_H >= |A|/2, ord m.
- ZE(n): no nondegenerate delta, |delta|<=1/(4D), in Z_P. Thm 5.7 ZE => n=O(loglog R) squarefree.
- R11: ZE-braid/SL/SL* false via cheap rational subtori (Sylvester-Hadamard, log2 n coordinate columns cost k*2floor(w/2)<m). ZE-col* => M=O((loglog R)^2) all norms. Prop 6.2: dim V=1 has nondeg zeros.
  T=1: m ~0.64n max, or m'=ceil(4w). Step 0; actual |delta|<=C e^{-W/4}.
- R12: T=1 reduction to MIN^pair; separation lemma dist(u,L)>= 1/2 exp(-N2 abar/(2P(L))); per-pair gaps |u_x-u_y|>=e^{-eps_xy/2}. Prop 2.7 nearstrata: no bound depending only on dim/degree/height/qualitative zero set reaches exponent d/2.
  J6: TC(A,m) ⊇ Bs of fat-point system; cheap spans = linear base cycles.
- R25: Z in T_M(Q), T_M=(R^1 G_m)^M/diag; ZE at actual point the one missing lemma.
## R16: (R_H) min_{G_1} Psi_{A,H'} > d^{2m'} eps_0^2; "arithmetic spent": box-level proof works only if NO admissible regime box (fakes incl.) has near-zero in G.

## IDEAS (round 26)
A. Extrapolation lemma (new?): rigidity at powers Z^j. (I1) holds for z^j (Gaussian ints, norm N^j, a->j a, e->j e, beta same).
   Liouville: ln Norm F_c(z^j) >= jDW - j W Lam(c)/2, Lam(c)=sum_p om_p h_{supp c}(beta_p).
   Taylor: ln Norm <= jDW + 2ln|c|_1 + 2m ln(DjC/2) - mW/2 - 2 ln m!.
   => P_c(Z^j)=0 if (m - j Lam)W/2 > err. Derivatives: c*k^alpha in V(A,m-|alpha|) => ord_{Z^j} P_c >= m - j Lam - O(err/W).
   So vanishing data: {Z^j : |j| < m/Lam}, order ~ m - |j| Lam. Extrapolation range J* < ord/Lam = LP ratio.
B. Spectral LP cap (generalise [D] Thm 2.3): G = sum_p om_p bt_p bt_p^T on 1-perp (bt = P beta). If lam_min(G)>=lam>0 then
   Omega(S) <= (sqrt(n)/lam) Lam(S) => ord(c) <= (sqrt n/lam) Lam(supp c) for all c. So J* <= sqrt(n)/lam.
   proof: P eps = sum om_p bt_p <bt_p, G^+ eps>; h subadditive, even; E|<G^+bt_p,eps>| <= |G^+ bt_p|_2 <= sqrt(n)/lam.
   careful: |bt_p|_2 <= |beta_p|_2 <= sqrt n.
C. Philippon needs Sigma(n) (n = dim G = d): with Sigma in <Z> (rank-1 group; only points where Liouville works), need progression of length >= d with orders >= d*T0 -> m/Lam >= d/2 + d/Lam. Contradiction with B for n >~ 4/lam^2.
D. Single-polynomial barrier: dimension count: for ANY point Y, exists nonzero P in V(A,m) vanishing at Y^j to orders T_j (sum binom(T_j+d-1,d) < dim V). So no single-polynomial zero estimate can separate Z from fakes. Need the whole system (overdetermination) - consistent with R11.
E. Dual picture: Z in Bs(V(A,m)) <=> k->Z^k on A in Poly_{<m}|_A. Exponential polynomial on G_a^d x G_m, phi(k)=(k,Z^k). Obstruction subgroups: V x G_m (projection-cheap: |pi_V(Sigma)| small) and V x {1} (only if Z^v=1 on V cap Z^d: excluded by nondegeneracy). Needs Sigma(d+1) subset A: at T=1 Sigma={0}; at T=n Sigma = T=1 slab box! But degree vacuity: |Sigma|<=3^d vs (d+1) m^d.
F. Obstruction subtori vs cheap subtori: K_F = exp(i L_F) in Bs; Z not on any proper subgroup (nondegenerate) => card(Sigma mod H)=|Sigma| for all H; obstruction = H={e} counting deficit, not a subtorus.
## record precedent
- round11/theory.md l.465 (summary remark, no proof): multiplicity estimates need many points; Thm 4.1 forces w^j into Z_P only for j <~ m/(2w) ~ sqrt(n/ln r), "far too few". Round 26 must make this rigorous & go beyond.
- web: arxiv/numdam blocked (network). Philippon statement = memory-quote.
## plan of theorems
T1 Extrapolation lemma (proved): ord_{Z^j} P_c >= t if (m-t+1) - |j|Lam(c) > E (E explicit).
T2 Spectral LP cap (proved): ord(c) <= (sqrt n/lam) Lam(supp c). Cor: J*(c) < sqrt(n)/lam.
T3 Torus-picture vacuity (Philippon, memory-quote form, + generous variant): H={e} satisfies inequality.
   general full-rank profile: Lam(c) >= 2 log5/W => J* <= T d W/log5; vacuity unless W >= (c d)^d.
T4 Single-polynomial barrier (Hilbert level, any zero estimate): for EVERY Y, nonzero c in V_H with route-5 vanishing data at Y^j (needs decay from T2).
T5 Dual picture G_a^d x G_m: subgroups V x {1}, V x G_m; Z-blindness for nondegenerate Z; degree version vacuous (|A| <= (2T+1)^d <= (d+1)(m-1)^d since m-1>=2T+1); Hilbert-sharp version false at cheap-subtorus points; any sharp Philippon-shape statement => ZE(n) for all nondegenerate points (delta-free) = barrier question.
T6 Z on obstruction subtori: none (nondegenerate); distance via character inequality |<v,delta>| >= 2 e^{-X_v/2}.
Conclusion: no ZE-col*/MIN^pair. Status no_new_result or obstruction_found (barrier for zero-estimate method) -> "obstruction_found".
