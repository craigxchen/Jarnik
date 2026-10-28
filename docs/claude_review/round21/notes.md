# Round 21 notes (working file; not the write-up)

## Imported facts (verified by reading)
- stability.md Thm 3.3(2): e_p log p <= 2 beta' < 4W'/M + 8 log^+ C;  beta' = max_{ceil(M/2)<=s<=M} b'_C(s),
  b'_C(s) = W'/s + (4(s-1)/s) log C - 8F*(s)/s^2 ; (3.0): beta' < 2W'/M - 3 log(M/2) + K'.
- Prop 3.5: W'_{<=Y} < pi(Y)(4W'/M + 8 log^+ C).  With Y = P_+(N'): W' < pi(Y)(...)  => M < 4 pi(Y) (C<=1).
  Since N' has only split primes: omega(N') <= pi(Y;4,1). => M < 4 pi(P_+(N');4,1) for C<=1.
- Cor 3.6: near-extremal: omega <= r <= (2/eps+o(1)) M; P_+ <= M^{4/eps-3+o(1)}.
- analytic.md Lemma 3.2(1): graft LS change <= sqrt(8 s_X (2K+1)) ||a||; Thm B needs s << X^{1/2}: fails in window.
- Known fakes: primes in [P, P e^{1/(4r)}] -> need P >~ r^2 log r.

## Key observations (to be proved / written)
1. Lower boundary of window: M < 4 pi(P_+(N');4,1) [C<=1]; general C: factor (1 + 2 M log^+C / W').
   => growth for P_+(N) = o(log R). Window in exponent of log R: [1, 2-delta).
   Near extremal (A = W'/(M log M)): omega >= (A/(4A-6))M(1-o(1)) -> M <= ((4A-6)/A) pi(Y;4,1)(1+o(1)).
2. Counting cannot go past exponent 1: near-extremal clusters have Theta(M) primes; pi(M^a) >> M for a>1.
3. Large sieve content in window = spreading of single prime angles at scale 1/K (lattice count already gives
   #{cluster angles in arc eta} <= eta X + O(sqrt X)). Short-arc condition constrains only linear forms
   <2(a_x-a_y), theta> at precision Delta = C e^{-W/4}; never single angles.
4. Spread-graft lemma: if grafted angles (and removed genuine ones) have <= m_K points per 1/K arc, change in
   LS sqrt-norm <= O(sqrt(K m_K))||a|| ... gives 1+o(1) even for s >> X^{1/2}. So Thm B's window failure is an
   artefact of worst-case Cauchy-Schwarz; LS is not what keeps the window fake-free.
5. The actual obstruction to window fakes: monomial-gap union bound sum_v e^{-w(v)/2} = prod_j(1+2 sum_m p_j^{-m/2}) - 1
   ~ exp(2 sum_j p_j^{-1/2}) >= exp(2 s P^{-1/2}) -> infinity exactly when s >> P^{1/2} (= LS non-triviality).
   Lattice-derived angles (arg beta, N beta <= p) satisfy gaps but give actual clusters (Thm D).
6. Real is invariant under subtorus H = ker of row-difference map; random point of H is pairwise uniform -> spread.

## Plan fixed (after scans launched)
Thm 1 counting boundary; Lemma 2 spreading (lattice sector count) with constants; Prop 3 Real contains spread
points (H-orbit, pairwise uniform coords); Prop 4 spread-graft extension of Thm B (Gallagher split into 2m
delta-separated classes + Montgomery-Vaughan (N-1+1/delta)); Cor 4.1 conditional cap; Prop 5 gap first-moment
threshold sum_j p_j^{-1/2}; heuristic Poisson => window fakes likely exist; computations: smoothscan.
LS split lemma: points with <= m in any half-open arc of length delta -> <= 2m delta-separated classes
(wrap-around handled: remainder block singletons). Constant 2m(N-1+1/delta), N=2K+1, delta=1/(2K): 8mK.

## Sec 4 design (Prop 4.1)
H = {h in T^s : B h in (pi/2)Z^{M-1}}, T=R/(pi/2)Z; Real(Delta)+H = Real(Delta) (shifts absorbed in k_x / theta_0).
Dual of H = Z^s / B^T Z^{M-1} =: Cok. D_j = order of e_j in Cok (inf = free). Pair (i,j) independent iff <e_i> cap <e_j> = 0 in Cok
-> (h_i,h_j) Haar on product of images. P(h_j+phi_j in I) <= 2 eta/pi + 1/D_j.
Grid arcs level l (length (pi/2)2^-l), 2^L <= D_0, t_l = sqrt(2 L n_J (1+2^l/D_0)) -> union prob <= 1/2.
Result: exists h: every arc length eta >= (pi/2)2^-L has <= (16/pi) s eta + 4 sqrt(L s) + (s - n_J).
LS (Lemma 3.2, K <= 2^{L-1}): Lambda <= 32 s + 32 K sqrt(L s) + 8K(s - n_J).
Full-rank remark: if rank B = s then H finite and Real(Delta) subset H + O(Delta ||B^+||): theta_j near (pi/2)Z/D_j;
monomial gap at v = D_j e_j forces D_j log p_j >~ W/2 - O(log): D_j >~ (A/(2a)) M -- stronger than LS's D >= s/(2.6 sqrt Y).

## Sec 5 derivation (Prop 5.1 spread graft)
Graft points: {phi_p, pi/2-phi_p} (2 per prime). Sp(L): every theta-arc of length eta >= (pi/2)2^-L has
<= (16/pi) s eta + 4 sqrt(Ls) of the s angles, and the 2s points are pairwise >= (pi/2)2^-L apart (both families).
K <= 2^{L-1}: m (for 2s points, x-arc 1/(2K) = theta-arc pi/(4K)) <= 8s/K + 8 sqrt(Ls); Lemma3.2 -> 8Km = 64 s + 64 K sqrt(Ls).
||D|| <= 2*8 sqrt(s + K sqrt(Ls)) ||a|| = 16 sqrt(...). t = 256 (s + K sqrt(Ls))/Lambda <= (256/c)(s/X + sqrt(Ls/X)).
K >= 2^{L-1}: MV with delta = 2^-L (x-units): 2K + 2^L <= 4K; ||D|| <= 2 sqrt(4K)||a|| = 4 sqrt K ||a||; t <= 16/(c sqrt X).
X >= P_0 always (else no change). Window: P_0 >= M^{1+eta}, s <= (2/eps+o(1))M, L = O(log M) => t -> 0.
R2/R3: Thm B proof verbatim, needs s log X_1 / X_1^gamma -> 0: gamma(1+eta) > 1.
## Sec 6: events B_v = B_{-v}: distinct events = (1/2) sum_{v != 0}. mu(B_v) = (4/pi)(arcsin(e^{-w/2}) + arcsin(e^{-w/2}/sqrt2)) >= 2.17 e^{-w/2}.
sum_{m in Z} p^{-|m|/2} = 1 + 2/(sqrt p - 1); log(1+x) >= x/(1+x) gives prod >= exp(sum 2/(sqrt p_j + 1)).
Monte Carlo plan: s ~ 6 primes, s/sqrt(P) ~ 0.5-1, test Poisson heuristic P(good) ~ exp(-E).
