# notes (p3Dlower)
N1. Pigeonhole over k-sparse 0/1 vectors among r columns vs prod_Q |G_q| gives v (norm<=2k) colliding on Q with
    log prod_Q q ~ log C(r,k) ~ k log(r/k); budget w(v)/2 >= k log p_min >= k log(r log r). => automatically satisfied
    (distinct primes). Entropy mechanism gives NOTHING (not even const*M log M beyond PNT). Parity removes 1 bit per (j,q):
    negligible for counting.
N2. Collision at q of v: <v,ell(q)> in mu_4 (4-torsion, = 2-part only). At q = +-3 mod 8 (m_q = 4 mod 8), mu_4 = full 2-Sylow,
    so parity (2-part bit) is INVISIBLE to collisions: odd-part labels free => coherent design (logs of small zeta_j) works there.
N3. At q = +-1 mod 8: twist design ell_j = log(zeta_j/zbar_j) + delta_j(q) log(rho/rhobar), rho non-residue norm, delta = mismatch.
    collision <=> q | D_s(Z_v rho^d), d=<v,delta(q)> in [-k,k]  => cap quadratic (2k+1)*k*log. Random-model: same Poisson tail
    truncated at quadratic cap => again tau ~ log^2M/loglogM (k ~ logM/loglogM is critical).
N4. 2-part-only twists lose 2^e (deep 2-adic levels) => also log^2/loglog by tradeoff.
N5. Uncertainty principle: kappa = A ell, A dense. Sparse v touching twin columns: <v,ell> = <c',kappa> + <dense, ell_free>.
    Column-coherent => pairs genuine-random (Poisson tails; maybe absorbed by cut LP of p3lower).
    Point-coherent => twin-column monomials random.
N6. Exact coherence with right parity: (Norm zeta_j / q) = (p_j/q) all q<=X => p_j Norm(zeta_j) square or pseudosquare.
N7. Any collision set lattice per q has index m_q/4 exactly; parity does not force collisions; lower bound must be an
    inverse theorem (linear caps => coherent) -- not available.
N8. TWIN DECOMPOSITION (key): basis labels (ell_J, U_y=ell_f+ell_s, Lam_y=eps_y ell_f(y)). kappa_y = F_y + Lam_y with
    F_y = sum_J a_yj ell_j + sum_z sigma_yz U_z. <v,ell> = <v_J,ell_J> + <v_s,U> + <c',Lam>, c'_y = eps_y(v_f-v_s).
    => sparse v = sparse comb of basis labels. Pairs see only kappa. Free labels decoupled from pairs.
    Mechanism must go through Lam = kappa - F.
N9. Hall/SDR: at q >> r, kappa class values can avoid all k<=2 hits. For k>=3, |C_{k-1}| >> q => hits forced per modulus;
    need cross-modulus balancing (discrepancy).  Random modifications at t<<q moduli add negligible agreement.
N10. Additive (integer-coherent) labels: Siegel relation of norm O(h) => collides at all q => tau >= ~X/k. dead.
N11. GRH: coherent at all q=+-1 mod 8 with log(N zeta p) << sqrt X => N zeta in p*sq or 2p*sq => |zeta|>=sqrt p:
    certificate saturated (no slack for completion's (4er)^k).
N12. Height-certified designs (torus exponent matrix N, A N small) + RIP => contradiction (uncertainty principle); parity-free.
