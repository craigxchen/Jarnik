# L5 working notes (not the write-up)

Plan of results:
1. Def of L5 with quantifiers; unit-1 only for pairs; log weights.
2. Translate form (kappa_y avoids <phi_y - eps w, f> + T); class-mate transfer identity; additivity over moduli.
3. Lemma: per-modulus avoidance (greedy) => counting cannot prove L5.
4. Thm: large moduli (m_q/2 >= c M) are free under IFL + Hall spread: switching among t shifts
   Lambda_y = lambda_y - 2<w_i^(y), f>; hits <= sum over i. Cor: L5 <=> L5 restricted to q <= 2cM+1.
   Lemma: random candidates satisfy Hall w.p. >= 1 - 1.5/n when c^(t-1) >= 2e^2.
5. Thm (volume): Gamma-coherent, G injective, col l1 <= B, pair l1 from x0 <= K_x:
   prod_x K_x * B^(r-M+1) * binom(r,M-1)^(1/2) >= covol(L) = sqrt(det Gram(a_x - a_x0)).
   Proof: singular values; prod of all sigma >= 1; prod of top k <= ||Lambda^k G||_F <= binom(r,k)^(1/2) B^k.
   covol(G L) >= covol(L) * prod_{M-1 smallest} sigma.
   Check: covol(G(L)) >= sigma_min-product: for a k-dim subspace V, vol_k(G|V) >= prod of k smallest sigma. yes.
6. Non-twin: decoupling iff gcd(Delta_M, m/2) = 1 (per q); section s; interface c'(v) = s^T v;
   sparse sections; no unimodular minor: image H contains d_M (Z/n)^M; injective kappa' exists if n/gcd(d_M,n) >= M.
   Non-split example A=(2), n=4.

Key quantities: n_q := m_q/2 = size of a parity class in Z/m_q.
Parity class of kappa_y fixed: kappa_y = F_y + nu_f(y) mod 2.
