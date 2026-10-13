# Round 10: the auxiliary-polynomial (determinant) method in the growth regime

Task: in the regime of near-extremal clusters (`Theta(M)` primes, all `<= M^{O(1)}`, near-Plotkin-tight), does an
auxiliary polynomial, combined with forced Gaussian divisibility, small-prime residues and the rank structure,
become supercritical (`W >> M log M`)? Analyse the LP/tropical criterion of `threshold/construct.md` restricted to
`r = O(M)` columns, find the extremal profile, and test against the fakes. Throughout, `n` is the number of
cluster points (the task's `M`); `M` is kept for the number of variables of a form, as in `threshold/`.

## 0. Summary and verdict

**No growth improvement is proved.** What is proved changes the picture recorded in the growth table
(`round8/README.md` §2.4), which lists forced divisibility as "exactly critical".

1. **The LP becomes strongly supercritical once columns are restricted (Theorem 2.1, proved).** The cube barrier
   of `threshold/` needs `2^n` prime columns. For *every* profile with `r` columns, a support chosen from the
   profile (a box cut by one slab per column) has `ord / sum_p omega_p h(beta_p) >= 0.12 sqrt((n-1)/ln(4r))`. No
   hypothesis on balance, Plotkin tightness, prime sizes or residues is needed. For the Plotkin-tight Hadamard
   profile the best ratio is at most `sqrt n` (Theorem 2.3), so `sqrt n` is the right order. An LP-critical
   profile needs `r >= e^{0.0144(n-1)}/4` (Corollary 2.4). **In the growth regime no profile is LP-critical**, and
   there is no extremal profile at the LP level (§3).
2. **Rigidity theorem for actual clusters (Theorem 4.1, proved, unconditional).** If `n - 1 >= 300 ln(4r)` and
   `W >= 16 ln n + 8 log^+ C + 14`, which holds for every cluster with `n >= c log W`, then the cluster
   satisfies an exact algebraic system. Its relative-angle vector `delta` makes
   `sum_{k in A} c_k e^{i<k,delta>} = 0` for every `c` in an explicit space `V_H` of integer measures of dimension
   `>= |A|/2`. Equivalently, `F_c(z) = 0` for all auxiliary forms in `V_H`.
3. **This input is true for actual circles and false for the fakes (Corollary 4.3).** The solution set `Z_P` is
   Lebesgue-null (Proposition 4.2). The round-5 fakes and the skew double core are drawn from positive-measure
   sets of positions, so almost every one of them violates the rigidity relation at the same `(n, r, W)`.
4. **The gap is a zero estimate (§5).** A contradiction needs one `c in V_H` with `F_c(z) != 0`. Every natural
   certificate fails or is capped:
   * parallelepiped (product) supports reduce exactly to single characters (Proposition 5.3);
   * valuation (Newton-polygon) certificates are impossible at primes with `e_p <= 3`, by the level-sum lemma.
     With residues they become a mod-`p` copy of the same rigidity, and residue collisions work *against*
     nonvanishing (Lemma 5.4, Proposition 5.5);
   * grid-lemma/tuple and interpolation-determinant arguments see the cube on balanced profiles and are
     critical (Proposition 5.6);
   * low-rank clusters (`rank <~ sqrt(n/log r)`) lie in `Z_P` automatically (Proposition 5.2).
5. **Squarefree clusters have linear rank (Lemma 3.1, proved).** Every normalised cluster on a squarefree norm has
   exponent rank `>= n/4 - 2` and contains a nondegenerate full-rank subcluster (Lemma 3.2).
6. **Conditional growth (Theorem 5.7).** Under a zero-estimate hypothesis ZE (no *nondegenerate* small `delta`
   lies on the base locus `Z_P`), clusters on squarefree norms satisfy `n = O(log W) = O(loglog R)`. ZE is an
   unlikely-intersection statement (a toric variety against a linear space of codimension `>= |A|/2`). It is not
   a standard conjecture, and I cannot prove it; small-case evidence is in §6.

| question in the task | answer |
| --- | --- |
| supercritical LP with `r = O(n)` columns? | **yes**, ratio `>= 0.12 sqrt(n/ln 4r)` (Thm 2.1); `<= sqrt n` at Hadamard (Thm 2.3) |
| extremal profile | none at the LP level; the true adversaries are low-rank profiles (§3, Prop. 5.2), which are not squarefree (Lemma 3.1) |
| does forced divisibility + auxiliary polynomial see something the fakes violate? | **yes**: the rigidity relation (4.2), true for all actual clusters in the regime, false for a.e. fake (Thm 4.1, Cor. 4.3) |
| small-prime residue collisions | enter only through face residue sums (Prop. 5.5); collisions favour vanishing, so they do not certify |
| rank structure (`weakest.md` Thm 1.3) | rank `<~ sqrt n` profiles are invisible (Prop. 5.2); squarefree clusters have rank `>= n/4 - 2` (Lemma 3.1) |
| `W >> M log M`? | **not proved**; it would follow (as `n = O(loglog R)`, squarefree norms) from Hypothesis ZE (Thm 5.7) |

Fake-test, as required: the mechanism's decisive input is the exact vanishing (4.2). It is **true** for actual
circles (Theorem 4.1) and **false** for almost every fake of size `O(n log n)` (Corollary 4.3). The missing step
(nonvanishing) is not satisfied by the fakes either, and is not required of them.

## 1. Setting, notation, and what is imported

Conventions. *Theorem / Proposition / Lemma*: proved here with all quantifiers, or quoted with an exact
source. *Exact check*: integer or exact-rational computation. *Evidence*: floating-point or finite
experiments that no proof uses. Scripts are in `round10/determinant_checks/` (outputs `*_out.txt`).

**Clusters.** `N >= 2`, `W = log N`, `R = sqrt N`. A *cluster* is a set of `n` distinct Gaussian integers
`z_1, ..., z_n` with `|z_x|^2 = N` lying on an arc of length `<= C sqrt R` of the circle `|z| = R`. Write
`z_x = R e^{i theta_x}`. The angular width of the arc is `Theta <= C sqrt R / R = C N^{-1/4}`. Put
`theta_0` = midpoint of the arc and `delta_x = theta_x - theta_0`, so `|delta_x| <= Theta/2`.

For a split prime `p = pi conj(pi)` with `p^{e_p} || N`, put `a_{p,x} = v_pi(z_x)` (so `0 <= a_{p,x} <= e_p`),
`a_p = (a_{p,x})_x in Z^n`, and the *centred column* `beta_p = (2 a_p - e_p 1)/e_p in [-1,1]^n`. The
*weights* are `omega_p = e_p log p / W`; `sum_{p split} omega_p <= 1`. Let `r` be the number of split
primes dividing `N`. The *profile* of the cluster is `(omega_p, beta_p)_{p <= r}`.

**Monomials on the zero-sum slice.** Let `Z^n_0 = {k in Z^n : sum_x k_x = 0}`. For `k in Z^n_0` and an even
integer `D >= |k|_1` put

```text
G_k(z) = N^{(D - |k|_1)/2} prod_x z_x^{k_x^+} conj(z_x)^{k_x^-}  in Z[i].
```

(`|k|_1 = 2 sum_x k_x^+` is even.) On a cluster, `G_k(z) = N^{D/2} e^{i <k, theta>} = N^{D/2} e^{i <k,delta>}`,
because `sum k = 0` removes `theta_0`. For a finite `A` in `Z^n_0` and integers `(c_k)_{k in A}` put
`F_c(z) = sum_{k in A} c_k G_k(z)`; this is the restriction to the cluster of a form of degree `D` on the
equal-norm variety `X_n` of `threshold/`.

**Order.** `ord(c) >= m` means `sum_k c_k P(k) = 0` for every polynomial `P` of total degree `< m` (on `Z^n_0`,
i.e. in the coordinates `k' = (k_1, ..., k_{n-1})` with `k_n = -sum k'`). `V(A, m)` is the space of such `c`.
`Poly_{<m}|_A` is the space of restrictions to `A` of polynomials of degree `< m`; `V(A,m)` is its
annihilator. `h_S(b) = max_{k in S} <k,b> - min_{k in S} <k,b>` is the width of `S` in direction `b`.

**Imported.**

* (I1) `threshold/construct.md` Theorem C1 (refereed): for every split `p`, `v_pi(G_k(z)) = e_p D/2 + <k, a_p>`
  and `v_{conj pi}(G_k(z)) = e_p D/2 - <k, a_p>`; inert primes and `1+i` divide every `G_k(z)` to the full
  power. Consequently, if `F_c(z) != 0` then `Norm F_c(z) >= N^D prod_{p split} p^{-h_{supp c}(a_p)}`.
  (The proof is four lines and is repeated in §4.1 for the reader's convenience.)
* (I2) `threshold/construct.md` Theorem B and Corollary B1 (refereed): for every finite `A` in a slice
  `{sum k = s}` of `Z^M` and every `c != 0` on `A`,
  `ord(c) <= reg(A) <= Omega(A) := 2^{-M} sum_{eps in {+-1}^M} h_A(eps)`.
* (I3) `round4/outside.md` Lemma 0.1 (gaps) and the angle-model definition; `round5/construct.md` §5 and
  `round5/exact.md` Theorem 5.1 (the fakes are obtained by drawing the positions `delta` from a set of
  positive Lebesgue measure).
* (I4) The pair inequality `d_xy := sum_p |a_{p,x} - a_{p,y}| log p >= W'/2 - 2 log C'` for normalised
  clusters (`round8/stability.md` §1, from `two_thirds/proof.md`).

**The question, made precise.** `threshold/construct.md` §2.2 shows that, with forced divisibility, one
auxiliary form `F_c` proves `F_c(z) = 0` exactly when the profile violates

```text
sum_p omega_p h_{supp c}(beta_p)  >=  ord(c) - O((ord(c) log C + log K)/W),                         (1.1)
```

and that the *uniform cube measure* on `{+-1}^n` satisfies (1.1) for every `c` (Theorem B). The cube
measure needs `2^n` distinct columns, hence `W >= 2^n log 5`, which is far outside the growth regime
`W = O(n log n)`. The task asks what happens when the adversary is restricted to `r = O(n)` columns
(more generally `log r = o(n)`), whether the method then becomes supercritical, and whether it sees
anything the fakes violate.

## 2. The LP/tropical criterion restricted to r = O(M) columns, primes <= M^O(1)

The cube barrier of `threshold/construct.md` (Corollary C3) is a statement about a *fixed* form against the
*uniform* measure on `{+-1}^n`. A cluster, however, is given first, and the form may be chosen afterwards,
from the cluster's own profile. Against a profile with `r` columns the right support is a box cut down by
`r` slabs, one orthogonal to each column. Its regularity grows like `n T` (dimension times side) while its
widths in the `r` profile directions grow only like `T sqrt(n log r)` (Hoeffding). So the method is
supercritical by a factor `sqrt(n / log r)`.

**Theorem 2.1 (slab-box supports).** Let `n >= 64` and `r >= 1`, and let `beta_1, ..., beta_r in [-1,1]^n`
and `omega_1, ..., omega_r >= 0` with `sum_p omega_p <= 1` be arbitrary. Put `T = n`,

```text
B_T = {k in Z^n_0 : |k_x| <= T for x = 1..n-1},     w = T sqrt(8 (n-1) ln(4r)),
A   = {k in B_T : |<k, beta_p>| <= w for p = 1..r},
m   = the largest integer with binom(m+n-2, n-1) <= (2T+1)^{n-1}/8,       D = 2(n-1)T.
```

Then:

* (a) `|A| >= (2T+1)^{n-1}/2 >= 4 binom(m+n-2, n-1) = 4 dim Poly_{<m}`;
* (b) `m >= 0.7 n(n-1)`, every `k in A` has `|k|_1 <= D`, and `h_A(beta_p) <= 2w` for every `p`;
* (c) with `U = |A|`, `H = (U T^{m-1})^2`, the integer vectors `c in V(A,m)` with `|c|_inf <= H` span a
  subspace `V_H` of `V(A,m)` of dimension `>= U/2`;
* (d) if moreover `n - 1 >= 300 ln(4r)`, then `m >= 4w`, so every `c in V_H \ {0}` satisfies

  ```text
  sum_p omega_p h_{supp c}(beta_p)  <=  2w  <=  m/2  <=  ord(c)/2 .                                   (2.1)
  ```

  In general (no condition on `r`), `m/(2w) >= 0.12 sqrt((n-1)/ln(4r))`.

*Proof.* (a) Draw `k` uniformly from `B_T`: the coordinates `k_1, ..., k_{n-1}` are independent and uniform on
`{-T, ..., T}`, and `k_n = -(k_1 + ... + k_{n-1})`. For a fixed `p`,
`X_p = <k, beta_p> = sum_{x<n} k_x (beta_{p,x} - beta_{p,n})` is a sum of independent mean-zero terms with ranges
of length `2T |beta_{p,x} - beta_{p,n}| <= 4T`. Hoeffding's inequality gives
`P(|X_p| > w) <= 2 exp(-2 w^2 / ((n-1)(4T)^2)) = 2 exp(-ln(4r)) = 1/(2r)`. A union bound over the `r` columns
gives `P(k notin A) <= 1/2`, so `|A| >= |B_T|/2 = (2T+1)^{n-1}/2`. The second inequality is the definition of `m`.

(b) `binom(a+b, b) <= (e(a+b)/b)^b` with `b = n-1` shows that
`m_0 = floor((n-1)(2T+1) 8^{-1/(n-1)}/e) - (n-2)` satisfies the defining inequality, so `m >= m_0`. For `n >= 64`,
`e 8^{1/(n-1)} <= e^{1 + ln 8/63} < 2.8097`, hence
`m >= (n-1)((2n+1)/2.8097 - 1) >= (n-1)(0.7118 n - 0.644) >= 0.7 n(n-1)` (the last step needs `n >= 55`).
For `k in B_T`, `|k|_1 <= sum_{x<n} |k_x| + |k_n| <= 2(n-1)T`. The width bound is the definition of `A`.

(c) Let `M` be the `E x U` moment matrix with rows indexed by `alpha in N^{n-1}`, `|alpha| < m`, entries
`k'^alpha` (`E = binom(m+n-2, n-1)`); its entries are integers of absolute value `<= T^{m-1} =: B`, and
`V(A,m) = ker M`. The `(H+1)^U` vectors `c in {0..H}^U` have images `Mc` in `[-UHB, UHB]^E`, so some fibre has
`>= (H+1)^U/(2UHB+1)^E` elements; their differences with one fixed element are distinct integer vectors of
`ker M` with `|.|_inf <= H`, hence lie in `V_H`. A subspace of dimension `d` meets `[-H,H]^U cap Z^U` in at most
`(2H+1)^d` points (it projects injectively onto some `d` coordinates). So
`d ln(2H+1) >= U ln(H+1) - E ln(2UHB+1)`. With `L = ln(UB) >= ln U > 300`, `ln(H+1) >= 2L`,
`ln(2H+1) <= 2L + ln 3`, `ln(2UHB+1) <= 3L + ln 3` and `E <= U/4` (by (a)):
`d >= U(1.25 L - 0.275)/(2L + 1.1) >= U/2`.

(d) Since `ln(4r) <= (n-1)/300`, `4w = 4n sqrt(8(n-1) ln(4r)) <= 4 sqrt(8/300) n(n-1) < 0.66 n(n-1) < m`.
Every `c in V_H` is supported in `A`, so `h_{supp c}(beta_p) <= h_A(beta_p) <= 2w`, and `sum omega_p <= 1` gives
(2.1); `ord(c) >= m` because `c in V(A,m)`. For the general ratio,
`m/(2w) >= 0.7 n(n-1)/(2 n sqrt(8(n-1) ln 4r)) = (0.7/(2 sqrt 8)) sqrt((n-1)/ln(4r)) >= 0.12 sqrt((n-1)/ln(4r))`.
QED

*Exact check* (`determinant_checks/params.py`, output `params_out.txt`): the exact `m` (computed from
`lgamma`, binary search) satisfies the closed-form lower bound on a grid of `(n, T)`, and the exact threshold
for `m >= 4w` at `T = n` is `n >= 2137` (`r = n`), `n >= 2679` (`r = n ln n`), `n >= 4278` (`r = n^2`), all
below the sufficient condition `n - 1 >= 300 ln(4r)`. The constants are poor and irrelevant.

**Remark 2.2 (what Theorem 2.1 does and does not say).** It is a statement about the linear programme of
`threshold/construct.md` §2.2 only. It does not need the pair inequality, balance, Plotkin tightness, the size of
the primes, or residues: *every* profile with `ln(4r) <= (n-1)/300` violates (1.1) for some admissible support,
with a factor-2 margin. The cube measure escapes only because it has `2^n` columns. What Theorem 2.1 does not
give is `F_c(z) != 0`; that is §4-§5.

**Theorem 2.3 (Plotkin-tight Hadamard profiles: the order sqrt n is sharp).** Let `H` be a Hadamard matrix of
order `n` whose first row is all ones, with rows `h_1, ..., h_n`. The *Hadamard profile* has columns
`beta_p = h_p` and weights `omega_p = 1/(n-1)` (`p = 2, ..., n`); every pair of points is separated by exactly `n/2`
of these columns, so it meets the pair inequality and the Plotkin bound with equality. For every finite `A` in a
slice `{sum k = s}` of `Z^n` and every `c != 0` on `A`,

```text
ord(c)  <=  sqrt(n) * Lambda_H(A),         Lambda_H(A) := sum_{p>=2} omega_p h_A(h_p).
```

*Proof.* `H^T H = n I`, so every `eps in {+-1}^n` equals `(1/n) sum_{p=1}^n <eps, h_p> h_p`. The width `h_A` is even
and subadditive, and `h_A(h_1) = 0` on a slice. Hence `h_A(eps) <= (1/n) sum_{p>=2} |<eps,h_p>| h_A(h_p)`. For `eps`
uniform, `<eps, h_p>` is a simple random walk of length `n`, so `E|<eps,h_p>| <= sqrt n`. Averaging,
`Omega(A) <= (sqrt n / n) sum_{p>=2} h_A(h_p) <= sqrt n Lambda_H(A)`, and (I2) gives `ord(c) <= Omega(A)`. QED

So for the Plotkin-tight Hadamard profile the best ratio `ord/Lambda` lies between `0.12 sqrt((n-1)/ln(4n))` and
`sqrt n`.

**Corollary 2.4 (where an LP barrier can live).** Call a profile *LP-critical* if every `c != 0` satisfies (1.1)
without error terms, `sum_p omega_p h_{supp c}(beta_p) >= ord(c)`. By Theorem 2.1 an LP-critical profile with
`r` columns on `n >= 64` points has `0.12 sqrt((n-1)/ln(4r)) <= 1`, i.e. `r >= e^{0.0144 (n-1)}/4`, and hence
`W >= r log 5` is exponential in `n`. The cube measure (`r = 2^n`) shows that this is the right order.
**In the growth regime `W = O(n log n)`, and in particular for every profile with `Theta(n)` primes, no
profile is LP-critical.** The single-place auxiliary-polynomial barrier of `threshold/` therefore does not
transfer to the growth question at the level of the linear programme; the entry "forced Gaussian divisibility:
divisibility is exactly critical" in the growth table of `round8/README.md` §2.4 is correct only for a form
fixed before the profile (Proposition 5.6 below makes this precise for the grid-lemma route).

## 3. Extremal profile in the restricted regime

**Answer at the level of the linear programme: there is none.** By Theorem 2.1 every profile with `r` columns
has `sup_c ord(c)/sum_p omega_p h_{supp c}(beta_p) >= 0.12 sqrt((n-1)/ln(4r))`. By Theorem 2.3 the Plotkin-tight
Hadamard profile has ratio `<= sqrt n`, so it is extremal up to the factor `sqrt(log n)`, but so is every
`O(n)`-column profile. Near-Plotkin tightness, balance and small primes play no role in the LP.

The real extremal objects appear only once nonvanishing is taken into account (§5.1): they are profiles of
**low rank**. What matters there is the dimension of the space spanned by the columns, not their number. The
next two statements show that for squarefree norms (in particular for every fake of rounds 5-6, all of which are
squarefree) the rank is linear in `n`, and that a full-rank subcluster of that size always exists.

**Lemma 3.1 (squarefree clusters have rank >= n/4 - 2).** Let `z_1, ..., z_n` be a normalised cluster
(gcd 1; `round8/stability.md` §1) on a squarefree norm `N'`, with arc constant `C'` and `W' = log N'`. Suppose
`W' >= 4 n log^+ C'`. Then the exponent lattice `span_Z{a_x - a_y}` has rank `k >= n/4 - 2`.

*Proof.* All `e_p = 1`. Let `S` be the `n x r` sign matrix `s_{x,p} = 2 a_{p,x} - 1` and `w_p = log p / W'`
(so `sum_p w_p = 1`). Let `G = S diag(w) S^T`, a positive semidefinite matrix with `G_xx = 1` and
`G_xy = 1 - 2 d_xy/W'`. By (I4), `G_xy <= eta := 4 log^+ C'/W' <= 1/n` for `x != y`. Also
`sum_{x != y} G_xy = 1^T G 1 - n >= -n`. Let `P` and `Q` be the sums of the positive and negative parts of the
off-diagonal entries. Then `P <= n(n-1) eta <= n` and `Q <= n + P <= 2n`. Since `|G_xy| <= 1`,
`||G||_F^2 <= n + P + Q <= 4n`. For a positive semidefinite matrix,
`rank G >= (tr G)^2/||G||_F^2 >= n^2/(4n) = n/4` (Cauchy-Schwarz on the eigenvalues).

Next, `rank G <= rank S <= rank(A) + 1`, where `A = (a_{p,x})` and `S = 2A - J`. Subtracting the first row from
the others changes the rank by at most one, so `k >= rank(A) - 1`. Hence `k >= n/4 - 2`. QED

(`round8/stability.md` Prop 3.4 proves `omega(N') > M/4` for the number of primes; Lemma 3.1 is the
corresponding statement for the rank, which is stronger. For non-squarefree norms the rank can be as small as
`log n/log 5.5` (`round8/weakest.md` Thm 1.3 is the matching upper bound on `n`), and I do not know a linear
lower bound.)

**Lemma 3.2 (full-rank subclusters are nondegenerate).** Let a normalised cluster have exponent rank `k`. Then
it has a subcluster of `k + 1` points whose exponent vectors are affinely independent. For any such subcluster
`J` and every nonzero `v in Z^J` with `sum v = 0`, `sum_{x in J} v_x delta_x != 0`.

*Proof.* A maximal affinely independent subset of `{a_x}` has `k+1` elements, because its affine span contains
every `a_x`. Write (1.1) of `round4/outside.md`:
`<2a_x - e, phi> = theta_0 + delta_x + (pi/2) kappa_x` with integer windings `kappa_x`. Then
`sum_x v_x delta_x = <2 u, phi> - (pi/2) <v, kappa>` with `u = sum_x v_x a_x`. Here `u in Z^r` is nonzero by affine
independence, so `<2u, phi> notin (pi/2) Z` by the gap lemma (`outside.md` Lemma 0.1). QED

**Conjectural extremal object.** In a rank-`k` profile, every quantity the method uses factors through the
`k`-dimensional lattice `Lambda = {(<k, a_p>)_p : k in Z^n_0}` (§5.1). There, Theorem B of `threshold/` (in
dimension `k`) makes a profile whose `r` functionals look like the `k`-dimensional cube critical, which needs
`r >= e^{ck}`. Such *low-rank cube-like* profiles need large exponents. A sign profile of rank `k` has at most
`O(k)` points by Lemma 3.1, so `e_p >= 2` is necessary. I did not construct one with pair separation at
`W = O(n log n)`, and I do not know whether one exists. This is the natural adversary for any
continuation of this route (§7).

## 4. Adding forced Gaussian divisibility: a rigidity theorem for actual clusters

### 4.1 The theorem

**Theorem 4.1 (rigidity).** Let `n >= 64`. Let `z_1, ..., z_n` be a cluster (§1) on `x^2 + y^2 = N` with arc
constant `C > 0`, whose norm has `r` split prime factors, with

```text
n - 1 >= 300 ln(4r)      and      W = log N >= 16 ln n + 8 log^+ C + 14.                          (4.1)
```

Let `(A, m, V_H)` be the slab-box of Theorem 2.1 built from the cluster's own profile `(omega_p, beta_p)`. Then

```text
F_c(z) = 0,   equivalently   sum_{k in A} c_k e^{i <k, delta>} = 0,      for every c in V_H.             (4.2)
```

Equivalently, the restriction to `A` of the character `f_delta(k) = e^{i <k,delta>}` lies in the subspace
`V_H^perp` of `C^A`, which has dimension `<= |A|/2` and contains `Poly_{<m}|_A`.

*Proof.* Let `c in V_H` be an integer vector and suppose `F_c(z) != 0`. Below, `D = 2(n-1)n` and
`G_k(z) = N^{D/2} e^{i<k,delta>}` as in §1.

*Upper bound (Liouville).* For real `x`, `|e^{ix} - sum_{t<m} (ix)^t/t!| <= |x|^m/m!`. For `t < m`,
`k -> <k,delta>^t` is a polynomial of degree `t` in `k'`, so `sum_k c_k <k,delta>^t = 0`. Hence

```text
|F_c(z)| = N^{D/2} |sum_k c_k (e^{i<k,delta>} - sum_{t<m} (i<k,delta>)^t/t!)| <= N^{D/2} |c|_1 (D Theta/2)^m / m!,
```

because `|<k,delta>| <= |k|_1 max|delta_x| <= D Theta/2`.

*Lower bound (forced divisibility, (I1)).* For a split `p` and `pi | p`,
`v_pi(G_k(z)) = e_p(D - |k|_1)/2 + sum_{k_x>0} k_x a_{p,x} + sum_{k_x<0} |k_x|(e_p - a_{p,x}) = e_p D/2 + <k, a_p>`
(using `sum k = 0`). Symmetrically, `v_{conj pi} = e_p D/2 - <k, a_p>`. Every `c_k` is a rational integer, so
`v_pi(F_c(z)) >= e_p D/2 + min_{supp c} <k,a_p>` and `v_{conj pi}(F_c(z)) >= e_p D/2 - max_{supp c} <k,a_p>`.
An inert `q` with `q^{2f} || N` divides every `z_x` exactly `f` times, hence every `G_k(z)` at least `fD` times;
likewise `1+i`. Taking norms,
`Norm F_c(z) >= N^D prod_{p split} p^{-h_{supp c}(a_p)}`. On `Z^n_0`, `<k, a_p> = (e_p/2) <k, beta_p>`, so
`h_{supp c}(a_p) <= (e_p/2) h_A(beta_p) <= e_p w` (Theorem 2.1(b)). Therefore
`ln Norm F_c(z) >= D W - w sum_p e_p log p >= (D - w) W`.

*Comparison.* `Norm F_c(z) = |F_c(z)|^2`, and `ln Theta <= ln C - W/4`. Therefore

```text
(m/2 - w) W  <=  2 ln|c|_1 + 2m ln(DC/2) - 2 ln m!.
```

By Theorem 2.1(d), `m/2 - w >= m/4`. Next, `|c|_1 <= U H = U^3 n^{2(m-1)}` with `U <= (2n+1)^{n-1}`, and
`ln m! >= m ln(m/e)`, and `e D/(2m) <= e/0.7 < 3.884` because `m >= 0.7 n(n-1)`. So

```text
W  <=  24 (n-1) ln(2n+1)/m + 16 ln n + 8 ln(3.884 C)  <=  2.61 + 16 ln n + 10.85 + 8 ln C.
```

The first term is at most `24 ln(2n+1)/(0.7 n) <= 2.61` for `n >= 64`. If `C < 1`, replace `C` by `1`; the arc
only gets shorter. This contradicts (4.1). Hence `F_c(z) = 0` for every integer `c in V_H`. By linearity the
same holds on the span. The dimension statement is Theorem 2.1(c). `Poly_{<m}|_A` is contained in
`V_H^perp` because `V_H` is contained in `V(A,m)`. QED

**Remarks.**

* Theorem 4.1 is unconditional. Its hypotheses hold for every cluster in the growth regime: `two_thirds` gives
  `W >= 3 n log n - O(n)`, and `r <= W/log 5`, so (4.1) holds as soon as `n` exceeds an absolute constant times
  `log W`. It holds a fortiori for every near-extremal cluster (`Theta(n)` primes, all `<= n^{O(1)}`).
* Theorem 4.1 is the honest output of "auxiliary polynomial + forced divisibility" in the restricted regime. The
  method is supercritical (Theorem 2.1). It does not produce a contradiction. It produces an exact
  algebraic relation that the cluster must satisfy.
* Only the column *directions* `beta_p` enter the construction of `A`. The weights enter only through
  `sum omega_p <= 1`, and the sizes of the primes not at all. Residues and rank are not used. §5 examines
  whether adding them closes the gap.

### 4.2 The fakes violate the rigidity relation

For a profile `P` (primes, exponents, rows `a_x`, windings), let `Z_P` be the set of `delta in R^n` with
`sum_k c_k e^{i<k,delta>} = 0` for all `c in V_H` (with `(A, V_H)` built from `P`).

**Proposition 4.2.** `Z_P` is a closed real-analytic subset of `R^n` of Lebesgue measure zero.

*Proof.* `V_H != 0` by Theorem 2.1(c). Fix `c in V_H \ {0}`. The function `g(delta) = sum_k c_k e^{i<k,delta>}`
is a trigonometric polynomial with distinct frequencies `k` and not all coefficients zero, hence not
identically zero. The zero set of a nonzero real-analytic function on `R^n` is Lebesgue-null, and
`Z_P` is contained in it. QED

**Corollary 4.3 (the input is of the mandatory type).** Let `G` be a set of positions `delta` of positive Lebesgue
measure such that `(P, delta)` is an angle model for every `delta in G`. The fake-circle constructions of
`round5/construct.md` §5 (step 4: the good `delta` have measure `> 1/2` of `[0, Delta]^M`) and `round5/exact.md`
Theorem 5.1 (skew double core; same mechanism) produce exactly such sets. Suppose `(P, n, r, W)` satisfies (4.1);
this holds for those fakes once `n` exceeds an absolute constant times `log W`, since they have `W = O(n log n)`.
Then almost every `delta in G` violates (4.2). The rigidity relation (4.2) is therefore

* **true for every actual cluster** in the regime (4.1) (Theorem 4.1), and
* **false for almost every fake** of rounds 5-6, including the skew double core, at the same `(n, r, W)`.

It escapes the reduction `round4/outside.md` Corollary A.2 for a simple reason. That corollary covers properties
of the angles that hold *off* a Haar-null set; (4.2) is a property that holds only *on* a null set.

*Proof.* `A` and `V_H` depend only on `P`, which is fixed as `delta` ranges over `G`. Apply Proposition 4.2. QED

So the determinant method, combined with forced divisibility and restricted to `log r = o(n)`, does **see
something the fakes violate**. What it does not do is turn this into a bound. For that one would need
`F_c(z) != 0` for at least one `c in V_H`, which is the content of §5.

## 5. Nonvanishing: what the residue, rank and determinant inputs can and cannot certify

To turn Theorem 4.1 into a bound one needs **one** `c in V_H` with `F_c(z) != 0`. This section examines each
natural source of such a certificate. Throughout, `(A, m, V_H)` is as in Theorem 2.1 and `Z_P` is as in §4.2.

### 5.1 Degenerate points of the base locus

**Proposition 5.1 (sparse positions).** If `delta_x = delta_n` for all `x` outside a set `S` of `{1..n-1}` with
`2T|S| < m`, then `delta in Z_P`. With `T = n` and `m >= 0.7 n(n-1)`, this covers every `delta` in which more
than `65%` of the points coincide.

*Proof.* With `sum k = 0`, `f_delta(k) = prod_{x in S} zeta_x^{k_x}`, `zeta_x = e^{i(delta_x - delta_n)}`, is a
function of `k_S in [-T,T]^S`. Every function on that grid is a polynomial of degree `<= 2T|S|` (Lagrange
interpolation in each variable). So `f_delta|_A in Poly_{<m}|_A`, which is contained in `V_H^perp`. QED

**Proposition 5.2 (low-rank subtori are invisible).** Let `S` be a set of split primes and let `kappa in Z^n`. If
`delta = sum_{p in S} t_p a_p + (pi/2) kappa + t_0 1` for some reals `t`, and
`sum_{p in S} h_A(a_p) + h_A(kappa) < m`, then `delta in Z_P` (mod `2 pi`).

*Proof.* `f_delta(k)` is a function of `y = ((<k,a_p>)_{p in S}, <k,kappa>)`, which ranges over an integer box with
sides `h_A(a_p)` and `h_A(kappa)`. On that box every function is a polynomial in `y` of degree at most the sum of
the sides (Lagrange), and `y` is linear in `k`. QED

For an actual cluster, `delta_x = 2 <a_x, phi> - (pi/2) kappa_x + const` (`outside.md` (1.1)). If the exponent
lattice has rank `k`, then `span{a_p} + R1` is spanned by `1` and `k` of the columns, so `delta` has the form of
Proposition 5.2 with `|S| = k`. In the slab-box, `h_A(a_p) <= e_p w`. **So Theorem 4.1 carries no information
for clusters of rank `k` with `sum_{p in S} e_p + h_A(kappa)/w < m/w`.** With bounded exponents this means
`k <~ 0.25 sqrt(n/ln 4r)`. Such clusters lie in `Z_P` automatically. By `weakest.md` Thm 1.3 they have
`n <= (5.5^k+1)/2`, which leaves room for them when `log n <~ k <~ sqrt n`. By Lemma 3.1 they are not squarefree.
This is the precise sense in which the extremal profiles of this route are the low-rank ones (§3).

### 5.2 Product (parallelepiped) certificates are exactly the character inequalities

The only supports on which nonvanishing is automatic are lattice parallelepipeds. On them the method reduces to
single characters, which are capped.

**Proposition 5.3.** Let `v_1, ..., v_d in Z^n_0` be linearly independent, `L_i >= 1`, and
`P = {k_0 + sum_i t_i v_i : 0 <= t_i <= L_i}`. Put `zeta_i = e^{i<v_i,delta>}` and `I = {i : zeta_i != 1}`.

* (i) `f_delta|_P in Poly_{<m}|_P` iff `sum_{i in I} L_i < m`. So some `c in V(P,m)` has `F_c(z) != 0` iff
  `sum_{i in I} L_i >= m`.
* (ii) `h_P(a_p) = sum_i L_i |<v_i, a_p>|`. Hence if `sum_{i in I} L_i >= m`, the forced-divisibility
  inequality fails at leading order,
  `sum_p log p h_P(a_p) < (m/2) W`, only if some `i in I` has `X_{v_i} := sum_p log p |<v_i, a_p>| < W/2`.
* (iii) For an actual cluster every `v` with `<v,delta> != 0` satisfies
  `X_v >= W/2 - log(C |v|_1/2) - 1`, the character (pair) inequality.

So every contradiction obtainable from a parallelepiped support is, up to the `O(log(C|v|_1))` error terms, a
contradiction already given by a single character. Characters are capped at `log R/loglog R` by the refereed
barriers (`round8/README.md` §2.4), and every fake satisfies all character inequalities.

*Proof.* (i) `t -> k_0 + sum t_i v_i` is an injective affine map, so `Poly_{<m}|_P` corresponds to polynomials of
degree `< m` in `t` on the box `Pi = prod [0, L_i]`, and `f_delta = e^{i<k_0,delta>} prod_i zeta_i^{t_i}`. On `Pi`,
`zeta^t = sum_{a <= L} binom(t,a) (zeta - 1)^a`, so

```text
f_delta|_Pi = e^{i<k_0,delta>} sum_{alpha <= L} prod_i (zeta_i - 1)^{alpha_i} binom(t, alpha),     binom(t,alpha) = prod_i binom(t_i, alpha_i).
```

The `binom(t,alpha)`, `alpha <= L`, form a basis of functions on `Pi`, and `binom(t,alpha)` has degree `|alpha|`.
Every monomial `t^beta` with `|beta| < m` restricts to a combination of `binom(t,alpha)` with `alpha <= beta`
(Stirling expansion; `binom(t_i, a) = 0` on `Pi` for `a > L_i`). Hence `Poly_{<m}|_Pi = span{binom(t,alpha) : |alpha| < m}`.
So `f_delta|_Pi` lies in it iff `prod_i (zeta_i-1)^{alpha_i} = 0` for every `alpha <= L` with `|alpha| >= m`. That
holds iff no such `alpha` is supported in `I`, i.e. iff `sum_{i in I} L_i < m`. The second sentence follows
because `V(P,m)` is the annihilator of `Poly_{<m}|_P`.

(ii) The width of a zonotope is the sum of the widths of its generators. Suppose every `i in I` has
`X_{v_i} >= W/2`. Then `sum_p log p h_P(a_p) = sum_i L_i X_{v_i} >= (W/2) sum_{i in I} L_i >= (m/2) W`.

(iii) Let `u = sum_x v_x a_x in Z^r`. Then `<v,delta> = 2<u,phi> - (pi/2)<v,kappa>`, and `w(u) = X_v`. By
`outside.md` Lemma 0.1 applied to `u`, `dist(<u,phi>, (pi/2)Z) >= arcsin(e^{-X_v/2})`. Hence
`|<v,delta>| >= 2 arcsin(e^{-X_v/2}) >= 2 e^{-X_v/2}` whenever `<v,delta> != 0` (mod `2 pi`; `delta` is tiny). Also
`|<v,delta>| <= |v|_1 Theta/2 <= (|v|_1 C/2) e^{-W/4}`. Taking logarithms gives the claim. QED

### 5.3 Valuation (Newton-polygon) certificates, and residues

A Gaussian integer is nonzero as soon as one of its valuations is finite. This is the natural way to use
"forced divisibility plus small-prime residues" as a *nonvanishing* tool.

**Lemma 5.4 (level sums vanish).** Let `c in V(A, m)` and let `lambda` be a linear form taking `s <= m` distinct
values on `supp c`. Then `sum_{k : lambda(k) = l} c_k = 0` for every `l`.

*Proof.* The push-forward of `c` under `lambda` is a measure on `s` reals whose moments of order `0, ..., s-1`
vanish (`lambda^t` has degree `t < m`). The Vandermonde matrix is invertible, so the measure is zero. QED

**Proposition 5.5 (valuation certificates).** Let `p` be split with `pi | p`. Write `z_x = pi^{a_x} conj(pi)^{e-a_x} y_x`,
`sigma = conj(pi) mod pi in F_p^*`, and `rho_x = y_x sigma^{-a_x} mod pi in F_p^*` (the point residues at `pi`).
For integer `c`, write `c_k = p^{v_k} c'_k` with `p` not dividing `c'_k`. Put `mu = min_k (v_k + <k, a_p>)` and let
`Phi` be the set of minimisers.

* (i) `F_c(z) / pi^{e D/2 + mu} == sigma^{eD/2} (N/p^e)^{D/2} sum_{k in Phi} c'_k sigma^{v_k} rho^k  (mod pi)`,
  with `rho^k = prod_x rho_x^{k_x}`. In particular, if this *face residue sum* is nonzero mod `pi`, then
  `F_c(z) != 0`.
* (ii) If `h_{supp c}(a_p) < m`, then `Phi` is never a single point. So no unique-minimiser certificate exists
  at such `p`. In the slab-box, `h_A(a_p) <= e_p w <= e_p m/4`, so this covers every `p` with `e_p <= 3`.
* (iii) Under (4.1), Theorem 4.1 forces every face residue sum to vanish mod `pi`, for every `c in V_H` and every
  split `p` (and the conjugate statement at `conj pi` with maxima). Residue collisions push towards this
  vanishing. If all `rho_x` are equal then `rho^k = 1` on `Z^n_0`, and when moreover every `v_k` is zero on
  `Phi`, the face sum reduces to a level sum, which vanishes by Lemma 5.4.

*Proof.* (i) From the proof of Theorem 4.1, `G_k(z) = pi^{eD/2 + <k,a>} conj(pi)^{eD/2 - <k,a>} (N/p^e)^{D/2} prod_x y_x^{k_x}`.
This uses `y_x conj(y_x) = N/p^e` and `sum_x k_x^- = |k|_1/2`. Reduce mod `pi`, writing
`p^{v_k} = pi^{v_k} conj(pi)^{v_k}`. Terms outside `Phi` have strictly larger `pi`-valuation.
(ii) Suppose `Phi = {k*}`. All `k != k*` on the level `l* = <k*, a_p>` have `v_k > v_{k*}`. By Lemma 5.4 (`lambda = <., a_p>`
takes at most `h + 1 <= m` values) `c_{k*} = -sum_{k != k*, level l*} c_k`, which is divisible by `p^{v_{k*}+1}`. That
contradicts the definition of `v_{k*}` (or `c_{k*} = 0`). (iii) If `F_c(z) = 0`, its `pi`-valuation is infinite, so
the reduction in (i) vanishes. QED

So the "small-prime residue" input does enter, through Proposition 5.5(i), as a polynomial condition mod `p` on
the point residues. Under rigidity it holds *with* vanishing. Turning it into a certificate requires showing
that an actual cluster's residues make some face sum nonzero. That is a mod-`p` version of the same rigidity
problem, and collisions (which are forced when `p < n`) work against it.

### 5.4 Grid-lemma and interpolation-determinant certificates

The classical way to get nonvanishing is the grid lemma, as in `threshold/verify.md`. Fix a form `F` in `M`
variables, nonzero on `X_M`; some tuple `t in [n]^M` of cluster points has `F(z_t) != 0`. The tuple, however, is not
under control. Laurent's interpolation determinants have the same feature: some minor is nonzero, but not a
chosen one.

**Proposition 5.6 (tuples see the cube).** Let the cluster's columns `beta_p` be balanced sign vectors (`n/2`
entries `+1`, `n/2` entries `-1`), with `sum omega_p = 1`. For every finite `A` in a slice of `Z^M` and every `c`
on `A`, `E_{t uniform in [n]^M} sum_p omega_p h_A(beta_p o t) = Omega(A) >= ord(c)`. In particular some tuple
`t` satisfies the leading-order forced-divisibility inequality for `c`. A grid-lemma argument, which must handle
the unknown nonvanishing tuple, is therefore at best critical on balanced profiles. This includes the Hadamard
profile and, up to the balance defects, every near-Plotkin-tight profile.

*Proof.* For `t` uniform, the coordinates `beta_p(t_1), ..., beta_p(t_M)` are independent fair signs, so
`beta_p o t` is uniform on `{+-1}^M` and `E h_A(beta_p o t) = Omega(A)`. Sum over `p` and apply (I2). QED

For `M = 1` (one circle, interpolation determinants with frequency set `v`), the support is the permutation
orbit of `v` and `h(eps) = 2(T_s(v) - B_s(v)) >= sum_{i<j} |eps_i - eps_j|` with `s = #{eps = +1}`. So against
*every* pair-separated profile the ratio is `<= 1`. This is the exact form of `growth/auxiliary.md`'s remark that
one-circle interpolation determinants are Vandermonde times a factor carrying no extra arithmetic.

### 5.5 The zero-estimate hypothesis and a conditional growth theorem

Everything above reduces the growth question for this method to one statement. Call `delta in R^n`
*nondegenerate* if `sum_x v_x delta_x != 0` for every nonzero `v in Z^n_0`. By Lemma 3.2, full-rank subclusters
have nondegenerate positions. The degenerate points of Propositions 5.1-5.2 are all degenerate in this sense.

**Hypothesis ZE(n).** For every profile on `n` points with `n - 1 >= 300 ln(4r)`, no nondegenerate `delta` with
`max_x |delta_x| <= 1/(4D)` lies in `Z_P`.

**Theorem 5.7 (conditional).** Let `n_1 >= 64` and assume ZE(`n'`) for all `n' >= n_1`. Let `z_1, ..., z_n` be a
normalised cluster on a squarefree norm `N'`, with arc constant `C'` and `W' = log N'`, satisfying
`W' >= 4n log^+ C'` and `W' >= 16 ln n + 8 log^+ C' + 14`. Then

```text
n  <  4 max(n_1, 300 ln(4r) + 1) + 4  <=  4 max(n_1, 300 ln(4 W'/log 5) + 1) + 4,
```

so `n = O(log W') = O(loglog R)`. Both conditions on `W'` hold for `n >= n_2(C)` by `two_thirds`.

*Proof.* By Lemma 3.1 the rank is `k >= n/4 - 2`. By Lemma 3.2 there is a subcluster `J` of `n' = k+1 >= n/4 - 1`
points with nondegenerate positions. `J` is a cluster on the same arc, and its profile has the same `r` columns
restricted to `J`. Suppose `n' >= max(n_1, 300 ln(4r) + 1)`. Then Theorem 4.1 applies to `J` ((4.1) holds since
`n' <= n`), so `delta_J in Z_{P_J}`. Moreover `max|delta_x| <= C' N'^{-1/4}/2 <= 1/(4D')`, because
`W' >= 4 ln(4C' n'(n'-1))`, which is implied by the hypothesis. This contradicts ZE(`n'`). Hence
`n/4 - 1 <= n' < max(n_1, 300 ln(4r) + 1)`, and `r <= W'/log 5`. QED

**Status of ZE.** It is an *unlikely-intersection* statement. `Z_P` is the preimage of a linear space of
codimension `dim V_H >= |A|/2` under the monomial map `T^{n-1} -> P^{|A|-1}`, `w -> (w^k)_{k in A}`. Its expected
dimension is negative, and all of its known positive-dimensional pieces near the identity lie in rational
subtori (Propositions 5.1-5.2). ZE asserts that there are no others near the identity. I cannot prove it. Neither
standard zero estimates (which concern groups, not one point) nor arithmetic Bezout (whose bound
`D^{n-2} h` is exponential in `n`, against `W = O(n log n)`) reaches it. Evidence from small cases is in §6.
Note also that ZE is a statement about *all* nondegenerate small `delta`, so it is stronger than what is needed.
For actual clusters `w = e^{i delta}` is an algebraic point with coordinates in `Q(i)` of height `~ W/2`, and only
such points must avoid `Z_P`.

### 5.6 The fakes, summarised

| input used by the auxiliary-polynomial argument | fakes (round 5, skew double core) |
| --- | --- |
| LP violation of (1.1) by the cluster's own profile (Thm 2.1) | **shared** by the fakes: it depends on the profile only, and the fakes' profiles have `log r = O(log n)` |
| rigidity (4.2) (Thm 4.1) | **false for a.e. fake** (Cor. 4.3) |
| parallelepiped certificates (Prop. 5.3) | satisfied: they reduce to characters, which the fakes satisfy |
| valuation certificates (Prop. 5.5) | not defined for fakes (no Gaussian integer `F_c`); for actual clusters forced to vanish |
| grid-lemma tuples (Prop. 5.6) | critical on balanced profiles: no contradiction |
| ZE (Hyp.) | false for a.e. fake *as a conclusion* (the fake's `delta` avoids `Z_P`, so no contradiction arises); ZE itself concerns `Z_P`, not the fakes |

The input that separates actual circles from fakes is the exact vanishing (4.2). What is missing is a way to
show that actual clusters cannot vanish there.

## 6. Checks and evidence

* `params.py` (exact binomials via `lgamma`; output `params_out.txt`). This checks the closed-form lower bound
  for `m` used in Theorem 2.1(b) on a grid. At `T = n` it finds the exact thresholds for `m >= 4w`: `n >= 2137`
  for `r = n`, `2679` for `r = n ln n`, `4278` for `r = n^2`. For `n <= 4096` the Theorem 4.1 threshold
  `W_1(n, C)` is at most `~150` (versus `16 ln n + 8 log^+ C + 14` proved). The ratio `m/(2w)` grows like
  `sqrt(n/ln r)`: `0.44, 0.59, 0.79, 1.07, 1.45, 1.96, 2.67` at `n = 64, ..., 4096` (`r = n`).
* `exact_checks.py` (exact Gaussian-integer arithmetic; output `exact_checks_out.txt`). It takes 40 windows of
  4 consecutive lattice points on `N = 5*13*17*29` and random supports `A` (40 points of the `T = 2` box), and
  for each it uses three kernel vectors of the `m = 3` moment matrix (exact nullspace). It checks:
  * the valuation formula `v_pi(G_k(z)) = e D/2 + <k,a_p>` of (I1): 19,200 checks;
  * the bound `v_pi(F_c(z)) >= eD/2 + mu` and the face-residue congruence of Proposition 5.5(i): 480 checks;
  * the vanishing of level sums (Lemma 5.4) whenever there are at most `m` levels: 13 checks.

  All pass. These are small instances of identities that are proved in general.
* **Evidence on ZE** (`baselocus2.py`, `baselocus3.py`; floating point; outputs `baselocus2_out.txt`,
  `baselocus3_out.txt`). Small slab-boxes with a random sign profile, keeping the half of the box with the
  smallest slab values, and the full space `V(A,m)` (not `V_H`). Method: grid search for zeros of
  `||(I - P) f_delta|_A||^2` on the whole torus `[-pi,pi)^{n-1}`, BFGS refinement of the best grid points, an
  integer-relation test, and a minimisation of the tangent-cone residual
  `tau(u) = ||(I-P) <k,u>^m||^2/||<k,u>^m||^2` over unit `u`.
  * `n = 4, T = 3, m = 5` (`|A| = 203`, `dim V = 168`): the only zero on the torus is `delta = 0`; `min tau = 7.9e-3`.
  * `n = 4, T = 4, m = 6` (`|A| = 369`, `dim V = 313`): the only zero is `delta = 0`; `min tau = 2.8e-3`.

  * `n = 5, T = 2, m = 5` (`|A| = 423`, `dim V = 353`; grid `20^4`, 150 refinements): 85 distinct zeros. Every one
    of them has exactly one nonzero coordinate (one point moved, all others coincident). That is precisely the
    sparse stratum of Proposition 5.1 (`2T|S| = 4 < 5`). No zero lacks an integer relation, and the tangent cone
    is attained at sparse directions (`min tau ~ 1e-30` at `u = e_3`).

  In all three cases the base locus consists of the identity and the degenerate strata predicted by
  Propositions 5.1-5.2, and nothing else was found, so ZE holds in these instances, even globally on the torus.
  These sizes are far below the regime of Theorem 2.1, the search is numerical, and it cannot exclude isolated
  zeros between grid points. It is evidence, not support for any proof.

## 7. Open points

1. **ZE (Hypothesis of §5.5), or any substitute.** The cheapest form that would suffice is this. *For the slab-box
   `(A, V_H)` of a squarefree profile, no algebraic point `w = e^{i delta} in (Q(i) cap S^1)^{n-1}` whose
   coordinates are multiplicatively independent and of height `<= W`, and which lies within `e^{-W/4}` of the
   identity, satisfies `sum_k c_k w^k = 0` for all `c in V_H`.* By Theorem 5.7 this gives
   `M = O(loglog R)` on squarefree norms, which is far beyond `o(log R/loglog R)`. It is an unlikely-intersection
   problem (toric variety against a linear space of huge codimension), with the arithmetic of `w` entering
   through its height. A first test is the tangent cone at the identity: the common real zeros of the degree-`m`
   forms `Phi_c(u) = sum_k c_k <k,u>^m`, `c in V_H`. Are all of them degenerate (in rational subspaces)?
2. **Freedom not used.** Theorem 4.1 holds for *every* admissible support. Examples are other `T`, other slab
   widths, unions over many slab-boxes, and supports adapted to sub-clusters. So an actual cluster in the regime
   lies in the intersection of all the corresponding base loci. A proof of nonvanishing may exploit this
   intersection rather than a single `Z_P`. The multi-scale structure near the sparse points of Proposition 5.1
   (a dense core plus outliers) suggests a recursion on sub-arcs.
3. **Non-squarefree norms.** Low-rank profiles are invisible (Proposition 5.2). Combining the method with the
   rank bound `n <= (5.5^k+1)/2` would need `k <= theta log W` with `theta < 1/log 5.5`. The method, granted ZE,
   gives `k <= c_0 log r` only with `c_0 ~ 300`. Better constants in Theorem 2.1 are possible (the counting loses
   `e`, and Hoeffding loses the variance), but probably not below `1/log 5.5`. A linear lower bound on the rank
   for near-extremal non-squarefree clusters, analogous to Lemma 3.1, would remove the issue.
4. **Low-rank cube-like profiles (§3).** Do pair-separated profiles of rank `k ~ C log n` with `r >= e^{ck}`
   columns exist at `W = O(n log n)`? If not, the LP is supercritical in the effective dimension for every
   near-extremal cluster.
5. **Status of the growth table.** `round8/README.md` §2.4 lists "single-place auxiliary polynomials ... forced
   Gaussian divisibility: divisibility is exactly critical" among the methods that give no growth. In the growth
   regime this is accurate only for a form fixed before the cluster (Proposition 5.6). With the form chosen from
   the cluster's profile the linear programme is supercritical (Theorem 2.1). The obstruction is nonvanishing
   (§5), not criticality.
