# notes (scratch)
N1. Root-of-unity parity repair: at q with k=v_2(m_q), the odd-label elements of T_q are exactly those whose
    2-part generates the 2-Sylow (order 2^k). f_q = primitive 2^k root (i if k=2 i.e. q=+-3 mod 8; zeta_8 if k=3).
    Design ell_j = lam(zeta_j) + eps_j(q) h_q. Collision of v (any unit) => (Z_v/conj Z_v)^(2^k) = 1 mod q
    => q | Z_v^(2^k) - conj(Z_v)^(2^k). For k=2,3 cost O(log|Z_v|). For large k cost 2^k log|Z_v| (bad).
    Non-root fixers xi (norm nonresidue): collision depends on E_v(q)=sum v_j eps_j(q) in [-kk,kk] -> quadratic.
N2. Coordinates: ell_s(x)=lam'_x - t_x, ell_f(x)=t_x. kappa = P + eps*t, P = Atilde ltilde.
    Monomial label = <u, ltilde> + <c', kappa>, c'_x = eps_x(v_f(x)-v_s(x)), u = vtilde - Atilde^T c'.
    v in L <=> u=0 and sum c'=0. Twin single e_f(x): u=-eps a~_x (weight ~r/2), c'=eps e_x.
N3. Triangle: pair label (kappa_x-kappa_y) = twin-diff label + P-diff label. At each q, not exactly two of the three
    collide (units compose). So structured pairs => twin-diff collisions = P-diff collisions off pair-collision primes.
N4. Robustness of Poisson tail: independent uniform labels on >= eps M columns at primes of reciprocal mass
    lambda >= (log M)^(-C) in [M^(1-eta),M] => height fails below c log^2 M/loglog M. Applies to twin columns for
    ANY fixed point design kappa if non-twin labels free (P uniform on image, independent across x).
N5. Counting heuristic: coherent (small zeta) + uniform margins: pair demands = smooth parts of huge integers,
    max over M^2 pairs ~ 2 log^2 M/loglog M; entropy b d M log M << M^2 rho(mu) needed => coherent fails (heuristic).
N6. Coherent + nonuniform margins (cut LP, p3lower Thm A): demand matrix random-like, lambda_max ~ M log M =>
    real-weight O(M log M). Open: discreteness (Q), parity at v_2>=4 primes, rigorous demand statistics.
N3 CORRECTED: labels pair = twin-diff + P-diff; ultrametric only (min attained twice). No forced inheritance.
N7. T4 (Paley b copies, height-certified coherent design G): rows t with psi_t != 0 (class sums) satisfy
    ||phi_t||_1 + ||F_t - med||_1 >= M/4 (Hadamard: ||y||_1 >= M||z||_1/||z||_0, z=(-m,psi/2)).
    => #T* <= 4(C+2C') ; pigeonhole on s-sparse 0/1 beta gives degenerate v=sum beta_a s_a, ||v||<=4b T0
    => collides at all q=+-3 mod 8 => need ||v|| >= X/(4 tau_max) => C + 2C' >= M/(64 b c log M).
N8. Parity balance heuristic: root fixer at class k certified cost 2^k log H; pseudo-random above K costs
    beta ~ A log M/(K log2); total 2^K d + A log M/K -> log M/loglog M. Same rate.
N9. Union-bound principle: any pseudo-random component with per-v mean collision count >= (log M)^(-C) * kk
    gives rate log^2/loglog; certification needed up to kk ~ M^eps.
PLAN: T1 (robust Poisson, rigorous), T2 (root repair lemma), T3(Paley height barrier = N7), heuristics H1/H2,
      numerics: (a) T2 check, (b) single-column moduli under coherent+root repair vs random, (c) T4 numbers.
PROGRESS: doc sections 1 (setting, (1.1)-(1.4)), 2 (Thm 1 + Cor 1), 3 (root repair Lemma 3.2, Prop 3.3,
Cor 3.4, H1, other repairs) written. check_root.py PASSED (89 prime powers, 17038 collisions).
NEXT: section 4 = Theorem 3 (Paley height barrier, from N7) + H2 + reduction remark (coherent + cut LP);
then check_tail.py + check_T3.py numerics; then section 0 verdict + section 6 open points + 5 checks.
PROGRESS 2: section 4 written (Def 4.1, Thm 3, Cor 3.6, scope, 4.2 open + H2).
CHECK Cor 3.6 constants: (4.1) holds -> K M log M >= (M-1)(M/4)/(8s) -> s >= M/(32 K log M); C+2Hbar >= s/8.
  => C + 2Hbar >= M/(256 K log M) (I wrote 330: fine since weaker? 1/330 < 1/256 ok). (4.1) fails ->
  4 s^2 (2sbC+1) >= M^2 -> with s <= 8(C+2Hbar)+1: 8 s^3 b C ~ M^2 -> (C+2Hbar)^3 C >= M^2/(8*512*b)=M^2/(4096 b)
  I wrote 2^11 b: need 2^12. FIX.
NEXT: numerics check_T3.py (uncertainty ineq on Paley 43/47; pigeonhole demo), check_tail.py; then sec 0,5,6.
