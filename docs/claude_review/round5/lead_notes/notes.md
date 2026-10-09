# Lead notes for round 5/6 (2026-10-09)

## 1. Which angle-model class decides the growth question

Notation as in `round4/outside.md`: squarefree-normalised profile, `W = log N`, characters
`v = sum lambda_x a_x`, margin `m(v) = w(v) - W/2`, arc width `Delta = C e^(-W/4)`.

**Worst-case residue class (A.3, property 3) is too strong for a lower bound.** Property 3 asks every
character, *pairs included*, to satisfy `dist >= arcsin(Q_M e^(-w/2))`, `log Q_M ~ 2.08 M`. For a pair,
`dist <= |delta_x - delta_y|/2 <= Delta/2`, so `w(a_x - a_y) >= W/2 + 2 log Q_M - 2 log(C/2)`.
Summing over pairs and using the cut bound `sum_{x<y} w(a_x - a_y) <= W M^2/4`:

```text
W >= 2(M-1)(2 log Q_M - 2 log(C/2)) ~ 8.3 M^2.
```

So in the A.3 class `8.3 M^2 <~ W_min <= (30+o(1)) M^2 log M`. This says nothing about actual circles,
because actual circles do **not** satisfy worst-case strengthening: only the moduli at which the pair
actually collides enter (the 2/3 proof uses exactly this). A.3 satisfies the worst case, which makes
its *barrier* stronger; a *lower bound* must use only inequalities valid for actual circles.

**Consistent-residue class (the right one for the growth question).** Angle model plus, for every
prime modulus `q` (inert: group of order `q+1`; split `q` not dividing `N`: `q-1`) and prime powers,
an assignment of prime residues `ell_j^(q)`, inducing the class map
`kappa_x^(q) = <a_x, ell^(q)> + (unit term) mod m_q`. A character `v = sum c_x a_x` is strengthened by
`prod {q : sum_x c_x kappa_x^(q) = 0 mod m_q}` (pairs: actual collisions). Actual clusters lie in this
class (genuine residues), every input of `two_thirds/proof.md` is valid in it, and
`3 M log M - O(M) <= W_min^cons <= W_min^(A.3)`.

## 2. Twins make every class map realisable

If row `x` has a twin pair (columns `j(x)`, `f(x)` differing only in row `x`), then the column span of
the exponent matrix contains `+-e_x` for every `x`, so it is all of `Z^M`. Hence for every modulus the
map `ell -> (<a_x, ell>)_x` is onto `(Z/m)^M`: any residue design (class maps) is realisable. Consistent
residues therefore impose **no** structure on which pairs collide, beyond the pigeonhole counts.

## 3. The obstacle for an O(M log M) fake: pair demands vs uniform margins

Pair `(x,y)` needs `m_xy >= D_xy - O(log C)`, where `D_xy = 2 sum_{q: collision} log q` (prime powers
alike). The pigeonhole forces average `D ~ 2 log M`. But:

* random designs give `max D ~ 4 log^2 M / loglog M` (Poisson tail of collisions at `q ~ M`);
* dyadic Reed-Solomon / affine-geometry designs give `<= 1` collision per dyadic range of moduli, but
  there are `~log M` ranges, so again `max D ~ log^2 M`;
* Hadamard-type profiles (forced by Plotkin tightness when `W = O(M log M)`) have *uniform* pair margins
  `m_xy = -(1/2)<s_x,s_y>_w`, `~ (t/2) log P` for `t` blocks.

So the multi-block construction gives `W = O(M log^2 M / loglog M)` at best unless (a) a residue
design with `max D = O(log M)` exists, or (b) the profile's Gram matrix can be adapted to the design
(`G_xy <= -2 D_xy`, realised by +-1 columns with near-equal weights). Open.

## 4. No growth improvement from the max demand

For actual circles the PSD bound with a nonnegative test vector gives only
`W >= 2 rho(D + P)` (Perron root), and `rho ~ M * average` unless the residue partitions for different
moduli are correlated. Consistent residues force no correlation (Section 2). So the max-demand
phenomenon obstructs the *construction*, not actual circles.

## 5. Multi-block random-permutation profile (data, `blocks.py`)

Per Hadamard block, excess `E(c) = ||H'^T Pi^T c||_1 - (M-1)`; random row subsets:

| base | M | s=2 | s=3 | s=4 | s=5,6,8 |
| --- | --- | --- | --- | --- | --- |
| Sylvester | 32, 64 | 1 always | M+1 always | 1 w.p. ~1/M (aligned quadruple), else M/2+1 | >= M/2+1 |
| Paley | 32, 60 | 1 always | M+1 always | >= 13 / 25 | >= 17 / 33 |

For three rows of any Hadamard matrix the four sign classes each occur `M/4` times, so `s = 3`
characters `2e_1 - e_2 - e_3` have excess exactly `M+1` (deterministic).

## 6. A residue design with max pair demand O(log M) (fixes Section 3 for small moduli)

Index the points by `x in {0, ..., M-1}`. Map every small modulus `q` (group order `m_q in {q-1, q+1}`)
injectively to a prime `p_q <= m_q` with `log p_q = log q + o(1)` (inert `q -> q`; split primes
`q -> previous split prime`; `5 -> 2` or `3`). Use the class map `x -> x mod p_q` (only `p_q <= m_q`
classes used), and for `q^a`: `x -> x mod p_q q^(a-1)`. A pair collides mod `q` iff `p_q | (x - y)`,
so its weighted demand is `<= sum_{p | (x-y)} 2 log(next prime after p) + (prime-power levels)
<= O(log M)`. With twins every class map is realisable.

**Caveat (class definition).** Strengthening must be restricted to moduli an argument can access by
pigeonhole (`(q-1) q^(a-1) < M`, or any fixed `<= M^O(1)` range). If *all* moduli are allowed with
free per-modulus assignments, a generic assignment makes some characters collide modulo infinitely
many `q` (Borel-Cantelli, `sum 1/q = infinity`), and any fixed integer labelling `lambda` has a short
`c` with `<c, lambda> = 0` (Siegel), which collides everywhere. Making the residues those of genuine
Gaussian integers built from the profile brings back genuine-residue statistics (max pair demand
`~ log^2 M / loglog M` for random-like primes). So:

* **class III-small** (pigeonhole-accessible residues): the sharp question; contains every input of
  the 2/3 proof; candidate fake below gives `W = O(M log M)`.
* **genuine-residue class**: max-demand effects, but no argument is known that can exploit them
  (Section 4).

## 7. Candidate theorem: `W_min^(III-small) = Theta(M log M)`

Construction: `t = O(1)` independent uniformly random row permutations of Sylvester `H'` (non-constant
columns), plus one twin column `f(x)` per row (copy of a block-`t` column flipped at row `x`, label map
`x -> a(x)` onto the labels). Primes `= 1 mod 4` in `[P, P(1 + M^-2)]`, `P ~ M^5` (Huxley-type short
intervals), `e = 1`, `k = 0`, residues by Section 6.

Margins (equal weights `w = log P`): `m(c) >= (w/2)[sum_i E_i(c) + E_t(c) - 1 - 2||c||_1]`,
`E_i(c) = ||H'^T Pi_i^T c||_1 - (M-1) >= 1`. Pairs: `m >= (w/2)(t - 4)`.

Needed per-block lemma (probability over `Pi`): for every non-pair `c`, `sum c = 0`, `n = ||c||_1 < 2M`,
`P(E(Pi c) < K n + 3) <= M^(O(L log M)) / multinomial(c)` outside explicit deterministic cases.
Tools:
* `F >= ||c||_1` (since `c = H T / M`), so `n >= 2M` is deterministic;
* small support `s <= sqrt(M/(2K))`: Hoelder chain gives `>= M/(2s)` exactly aligned columns, and for
  Sylvester these form a coset of `span(S-S)^perp`, so `S` lies in an affine subspace of dimension
  `<= log2(2s)`; aligned quadruples have probability `~1/M` per block (data agrees);
* large support: `||f||_A = F/M <= K + 1` for `f = Pi c` on `F_2^m`, so by the Green-Sanders quantitative
  idempotent theorem (via level sets `(f + f^2)/2` etc., algebra-norm submultiplicativity) `f` is a
  `+-` sum of `L = L(K)` coset indicators; there are only `M^(O(L log M))` such functions.
Union bound: `sum_c p(c)^t <= sum_mu multinomial(mu)^(1-t) M^(O(t L log M))`, small for `t >= 3` once
the small-support and spiky cases are removed. Then the delta union bound, Theorem A.1 completion
(`P >> r^2`), and pair demands `O(log M)` (Section 6) give `W = (t+1) M log P = O(M log M)`.
Constants are astronomical (Green-Sanders), but `t` is absolute.

## 8. Checks (`design.py`)

* Design of Section 6, moduli with group order `< 2M`: max weighted pair demand
  `3.56, 4.03, 3.79, 4.18` times `log M` for `M = 64, 256, 1024, 4096` (stable `O(log M)`).
* Multi-block Sylvester, all `+-1` characters of support 4, minimum of `sum_i E_i`:
  `t = 1`: 1; `t = 3`: 3 (flat in every block, about `M^4 M^-3` of them) or 19; `t = 5`, `M = 64`: 69
  (none flat in all blocks). Matches the per-block probability `~1/M` of an aligned quadruple.

## 9. Per-block lemma: proof route (Sylvester, `G = F_2^m`)

For `f = Pi c` on `G`, `||f||_A := sum_a |f^(a)| = F/M` (`f^(a) = E_x f(x)(-1)^<a,x>`). Then:
1. `|c_x| <= ||f||_A` (so near-flat characters have bounded entries) and `F >= ||c||_1`.
2. Small support `s` (`s E < M`): Hoelder chain `sum_b |T_b|(n - |T_b|) <= n E` gives
   `#aligned >= M/n - E/2`; aligned columns form a coset of `span(S-S)^perp`, so `|span(S-S)| < 2s` is
   the least power of two `>= s`: `S` lies in a coset of size `2^ceil(log2 s)`. Non-+-1 small-support
   characters have `E >= 2M/n - 1`. Probability per block `<= M^(d+1-s) s^O(s)`, `d = ceil(log2 s)`;
   support 4 needs `t >= 5`, `s >= 6` needs `t >= 4` (odd `s` cannot be +-1).
3. Large support: level sets `1_{f=v} = prod_{u != v}(f-u)/(v-u)` have bounded algebra norm
   (Wiener algebra is a Banach algebra), so by Green-Sanders (Ann. Math. 2008; Sanders 2019 for
   `F_2^n`) each is a `+-` sum of `L = L(K)` coset indicators: `<= M^(O(L log M))` structured `f`.
4. Union bound over value multisets `mu`: `sum_mu multinomial(mu)^(1-t) M^(O(t L log M))`, small once
   `s >= C L log M` (multinomial `>= binom(M, min(#neg, #pos))`, and bounded entries with `sum c = 0`
   force both signs to have a positive proportion of the support).
Result: for `M >= M_0` (astronomical), `t = 5` blocks suffice; rows can be deleted (sub-Hadamard), so
every large `M` works.
