# L4K working notes
N1. Budget: per-unit height m_(v,s) <= nu_s^{-1} e^{w(v)/2}, nu_s = 1 (s even), 2^{-1/2} (s odd).
    Collision any unit iff <v,ell> = 0 mod n_q, n_q = m_q/4. Unit s = <v,ell>/n_q mod 4.
N2. Pigeonhole forced region: S = {||v||_1 <= floor(K/2)} values mod n_q; collision forced iff |S| > n_q
    => q <~ 4|S_{K/2}| ~ (2r)^{K/2}. README §6 sentence "pigeonhole forces #V_K/q collisions at q<=Q_K" is the
    random expectation, not forced.
N3. Free region (collision-free parity-respecting vector exists):
    (a) greedy, sharpened Prop H: m_q > 2 F(r,K), F = 4 sum_{v: v_r>0, ||v||<=K} v_r  (~ (2r)^{K-1}).
    (b) explicit: a_j = 2K'(Lambda + b_j) + nu_j, b B_h set (h=floor(K/2)), with leveling => m > C_K r^{floor(K/2)}.
        K=2: a_j = 4j + nu_j works for m/4 > 8r+2.
N4. Literal L4(K) on [1,Q_K] refuted: genuine below free threshold, random greedy above; union bound on
    coherence over zeta (norm <= H), u_q, w, subsets.
N5. Forced region: the real remaining question. Coherence non-vacuous there (H vs e^{theta(T)}).
N6. Thm 2 (refutation): genuine at n_q <= A_K, (2.1) with random sigma_q above. P(coherent at q | zeta) <= n_q/r!.
    Union: (13H)^r 2^{pi_7(Y)} (Y/(4 r!))^{k0}, k0 >= delta theta_7(Y)/log Y. Needs r -> infinity (r = Theta(M) in
    near-extremal regime, round8 stability). K >= 2 (K=1 trivial, A_1 > Q_1).
N7. Finite instance: genuine below q*, random-DFS free vectors above; exact goodness check all levels <= Q_K;
    pair-coherence max mass over zeta pairs norm <= H via unit-orbit table T_q[a,b].
N8. forced_sig r=4 first K=2 (tightest budgets): 5.07e15 good level-1 systems on F={3,7,11,19,23,31,41};
    1.37e11 (=2^37) have genuine signature at every q. Loose sets: DFS over signature tuples explodes -> sample.
N9. Levels: level-2 lift l+mt; colliding v (level 1, <v,l> = n d_v) collides at level 2 iff q | d_v + 4<v,t>.
    Same unit as level 1. One forbidden residue of <v,t> mod q per colliding v.
N10. Coherence check on forced moduli: DFS over q choosing u_q, buckets residue->zeta set per column.
N12. K=3 tight one-modulus: rigid q<=41 (=q*(K=2)), flexible from 43. K=4 tight: rigid q<=103 (=q*(K=3)) except 83;
     flexible 107..199. core_sig: core [3..47] K=4 with rest genuine -> ONLY genuine class (1 signature tuple).
     Without rest: 87.5% non-genuine at q>=11. K=3 core [3..41] w/o rest: 0.35% genuine everywhere.
     prefix: genuine q<=103, top 107..199 searched -> 3/5 systems non-genuine at ALL 18 top moduli, not coherent
     (H=2000) on top Q_7 [127,151,167,191,199] (~half of forced Q_7 mass).
