# L4 working notes (not the write-up)

## N1. Lattice reformulation
At q (top level, ignore higher levels first): n_q = m_q/4. Collision of v (any unit) iff <v,ell(q)> in (m/4)Z/m
iff <v,c(q)> = 0 mod n_q, c = ell mod n_q. Kernel Lambda_q = {v: <v,c(q)> = 0 mod n_q}.
Q_1 = {q = +-1 mod 8}: n_q even, parity of c = parity of ell = nu(q); Lambda_q subset H_q = {v: <v,nu(q)> even}
(when c(q) surjective: Lambda_q + 2Z^r = H_q). Q_3: n_q odd, no constraint (Thm 1.1(1)).

## N2. Certificates from any number field are parity-locked (GRH)
Parity of x/x^sigma in T_q = (N x/q). Need p_j N(x) pseudosquare on Q_1 -> p_j N x = square or 2 square.
=> any global certificate has height >= sqrt(p_j) on e_j. Zero slack.

## N3. Random range
Random parity-respecting labels: union bound OK only for k >= M^{c/beta} (Poisson hits h log h ~ k log M).
Hard range: 1 <= k <= M^{O(1/beta)}.

## N4. Piecewise certification
Twisted = piecewise certified with T+1 pieces sharing one twist family. Cost quadratic (cap).
Parity-free: P pieces each certified => demand <= P * cap: linear for P = O(1).
Parity: pieces must realise chi_{p_j}: need P ~ log M pieces per column (pigeonhole), cap then counts
realized witness combinations.

## N5. Parity height bound only on odd coordinates
If ell = B lambda, lambda certified by xi, GRH: N(Xi_u) prod p_j^{v_j} = (1 or 2) square => |Xi_u|^2 >= prod_{v_j odd} p_j.

## N6. Plan of results (decided)
- Regime: stability Cor 3.6 => near-extremal (W <= c0 M log M) actual clusters have p_j <= M^{h0}, h0 = 2c0-3+o(1).
  For growth it suffices to prove (T') : for each c0 there is A=A(c0) ... So X = 3M^A >> p_j^C allowed.
- Thm A (uncond, Linnik): exact coherence at q = 7 mod 8 (q not | N automatically since q=3 mod 4) forces pN(zeta) in {t^2, 2t^2}
  once X >= C (8 p N zeta)^L.  Key: (n/q) = sigma_m (q/m) for q = 3 mod 4, m odd sqf part.
- Thm B (uncond, Heath-Brown quadratic large sieve): all but << X^eps delta^-2 log D columns: every sub-genuine zeta
  agrees in parity on <= (1+delta)/2 of theta(Q_7(X)).
- Thm C (uncond Linnik / GRH): P sub-genuine pieces cannot cover column j if X >= C(8 prod n_c)^L; GRH: need
  P >= (1/2 - o(1)) log_2 X. Duality: relations prod_S psi_c = 1 have |S| even and sign product +1.
- Prop D (k<=2 free): labels: Q_3 unit shift; S0={17,31,47} parity-matched zeta; b_q=3 CRT (parity bit is direct factor Z/8);
  b_q>=4: odd injective orbit-encoding of comp_1 mod N/4 (N=(q-chi)/4>=16), q-part homomorphic.
  For ||v||_1<=2: collision at q^a => q^a | Im(Z_v^16).  => m_v <= |Z_v|^16.
- Prop E (known, p3Dconstruct L3.2): classes b<=B cost 2^{B+1}; obstruction lives in q = +-1 mod 2^B, B->infty.
- Additive labels die even at bounded norm via highly composite moduli (q-1 | L): remark.
