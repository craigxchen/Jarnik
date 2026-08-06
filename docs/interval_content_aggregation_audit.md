# Interval-content aggregation: exact cancellation and surviving constraints

## Status

This note audits whether the transition-content cocycle can be summed over all
ordered intervals to beat the endpoint cut bound.  The answer for scalar
content is negative and exact: after summing every interval, the content term
cancels prime by prime and one recovers precisely the old pairwise cut
identity plus the angular Vandermonde.

Several non-scalar consequences survive:

- a prime word with many switches creates cubic aggregate interval content;
- almost every three-vertex split has zero quadratic carry;
- more than four critical points force linear-in-height Gaussian overlap
  energy; and
- finite Gaussian character codes can realize all cocycle identities only by
  paying excessive lcm height.

These are rigorous paper-level statements in the odd squarefree split model,
but they do not prove the uniform endpoint bound and are not yet formalized in
Lean.

---

## 1. All-interval product identity

Let the ordered vertices be `0,...,n-1`.  Write

```text
E_r = B_(r,r+1),
q_r = Norm(E_r),
W_r = (r+1)(n-1-r).
```

For every interval `a<b`, primitively reduce the consecutive product:

```text
product_(r=a)^(b-1) E_r = C_(a:b) B_ab,             (1.1)
```

where `C_(a:b)` is its positive rational content.  The edge `r` occurs in
exactly `W_r` intervals, so multiplying (1.1) over all `a<b` gives

```text
product_r E_r^(W_r)
  = (product_(a<b) C_(a:b))(product_(a<b) B_ab).     (1.2)
```

Taking norms,

```text
sum_r W_r log q_r
  = 2 sum_(a<b) log C_(a:b)
      +sum_(a<b) log q_ab.                          (1.3)
```

Let

```text
B_ab = X_ab+i delta_ab,
Theta_ab = arg B_ab > 0.
```

Then `delta_ab=sqrt(q_ab)sin Theta_ab`, and (1.3) yields

```text
sum_(a<b) log delta_ab
  = (1/2)sum_r W_r log q_r
      -sum_(a<b)log C_(a:b)
      +sum_(a<b)log sin Theta_ab.                   (1.4)
```

At first sight, a lower bound on total content in (1.4) appears capable of
forcing a residual deficit.  The prime-word calculation below shows that its
entire contribution is already paid by the weighted edge norms.

---

## 2. Primewise cancellation

Fix one conductor prime and write its binary word as

```text
e_0,e_1,...,e_(n-1).
```

Put

```text
s_r = |e_(r+1)-e_r|,
T_(a:b) = sum_(r=a)^(b-1) s_r.
```

The transition-content cocycle gives

```text
v_p(C_(a:b))
  = floor(T_(a:b)/2)
  = (T_(a:b)-|e_b-e_a|)/2.                          (2.1)
```

If `n_0,n_1` are the two class sizes of the word, summing (2.1) over all
intervals gives

```text
2 sum_(a<b) v_p(C_(a:b))
  = sum_r W_r s_r-n_0 n_1.                          (2.2)
```

The first term on the right is exactly the prime contribution to the
weighted edge norm in (1.4); the second is exactly the prime contribution to
all pairwise transition norms.  Substituting (2.2) into (1.4), and then
summing over conductor primes, gives the exact identity

```text
sum_(a<b) log delta_ab
  = (1/2)sum_p n_(0,p)n_(1,p) log p
      +sum_(a<b)log sin(theta_b-theta_a).            (2.3)
```

Thus total interval content contains no scalar information beyond the old
pairwise cut metric.  Any improvement must use how contents and coordinates
are coupled across overlapping intervals, not their first moment.

---

## 3. Exact growth of aggregate content

Suppose the fixed prime word has `t` switches and consecutive run lengths

```text
ell_0,...,ell_t,
sum_i ell_i=n.
```

Its total exponent in all interval contents is

```text
K_p = sum_(i<j) ell_i ell_j floor((j-i)/2).          (3.1)
```

For fixed `n,t`, the exact minimum is

```text
K_min(n,t) = F(t)+(n-t-1)A(t),                      (3.2)
```

where

```text
F(2u)   = u(u+1)(4u-1)/6,
F(2u+1) = u(u+1)(4u+5)/6,

A(2u)   = 2 floor(u^2/4),
A(2u+1) = u(u+1)/2.
```

To see this, set `ell_i=1+x_i`.  The expansion of (3.1) is `F(t)` plus
linear terms in the `x_i` and nonnegative quadratic cross terms.  The minimum
is obtained by putting all extra length in a central run.  In particular,

```text
K_p >= F(t) >= (t^3-t)/12.                          (3.3)
```

So a word that switches linearly often creates cubic aggregate content.  This
large quantity is real, but (2.2) shows why it does not yield a residual gain:
the same repeated switches create exactly the weighted edge-norm term needed
to cancel it.

---

## 4. A quadratic re-entry inequality

Let there be `L` positive consecutive transition angles `g_r`, with

```text
sum_r g_r <= eta.
```

For a prime `p`, let `t_p^+` and `t_p^-` count its `0->1` and `1->0`
switches.  Alternation gives

```text
t_p^+ t_p^- = floor(t_p^2/4),
t_p=t_p^++t_p^-.
```

For two edges `r<s`, let `P^-_(r,s)` be the product of primes that switch in
opposite directions on those edges.  This product is the rational content of
`E_r E_s`.  The primitive imaginary coordinate after division is a positive
integer, so

```text
P^-_(r,s)
  <= sqrt(q_r q_s) sin(g_r+g_s)
  <= sqrt(q_r q_s)(g_r+g_s).                        (4.1)
```

Multiply (4.1) over all edge pairs.  Since

```text
sum_(r<s)(g_r+g_s)=(L-1)eta,
```

AM--GM gives

```text
sum_p floor(t_p^2/4)log p
  <= (L-1)/2 sum_p t_p log p
      +binom(L,2)log(2eta/L).                       (4.2)
```

This is an exact quadratic re-entry/gap constraint.  At endpoint scale its
negative term is `-binom(L,2)H/4`, but balanced random words satisfy it with
room to spare.  It controls the switch distribution without removing the
diagonal `LH` slack.

---

## 5. Quadratic carry and the near-complement matching

For `a<b<c`, factor the interval product as

```text
B_ab B_bc = gamma_abc B_ac,                         (5.1)
```

where

```text
gamma_abc = C_(a:c)/(C_(a:b)C_(b:c))
```

is the cross-content.  Taking imaginary coordinates and norms gives the
exact identity

```text
delta_ab delta_bc/gamma_abc
  = sqrt(q_ac) sin Theta_ab sin Theta_bc.            (5.2)
```

If the whole angular diameter is at most `C exp(-H/4)`, then

```text
delta_ab delta_bc/gamma_abc
  <= (C^2/4)sqrt(q_ac/N).                            (5.3)
```

Consequently a nonzero integer carry

```text
floor(delta_ab delta_bc/gamma_abc) >= 1
```

forces

```text
q_ac >= 16N/C^4.                                    (5.4)
```

Such outer pairs are near complements in the weighted Hamming metric.  Two
near complements of the same vertex are within `O_C(1)` of each other, which
contradicts the pair lower bound `H/2-O_C(1)` for large `H`.  Hence the graph
of outer pairs admitting a nonzero carry has maximum degree one.  Outside at
most a matching, every split `a<b<c` satisfies

```text
delta_ab delta_bc < gamma_abc.                      (5.5)
```

This is a genuine simultaneous restriction, but direct summation still
fails: one prime with word `010/101` can contribute its cross-content to
quadratically many triples.  At exact pairwise-critical character balance,
every triple has cross-content height `H/4`, so all carries may vanish while
the same blocks are reused cubically.

---

## 6. Lcm height and overlap energy

Let `w_1,...,w_n` be primitive Gaussian direction vectors, and set

```text
h_i = log Norm(w_i),
s_ij = log Norm(gcd(w_i,w_j)),
H = log Norm(lcm_Gauss(w_1,...,w_n)).
```

Prime-power-layer inclusion--exclusion gives the universal inequality

```text
H >= sum_i h_i-sum_(i<j)s_ij.                       (6.1)
```

The raw determinant of each distinct pair is a nonzero integer.  Endpoint
diameter therefore gives

```text
h_i+h_j >= H/2-2log C.                              (6.2)
```

Summing (6.2) and using (6.1) yields the necessary overlap bound

```text
sum_(i<j)s_ij
  >= (n/4-1)H-n log C.                              (6.3)
```

Thus a critical family with more than four points must have overlap energy
linear in the full conductor height.  In particular, comparable-size vectors
with subpower pairwise Gaussian gcds admit at most four critical points.

This cleanly excludes the weak-overlap or junk-prime regime.  It is not sharp
enough in the balanced character regime, where the pairwise overlap energy is
quadratic in `n` and (6.3) has substantial slack.

---

## 7. Exact stress tests

### 7.1 Content need not enter the residuals

Let

```text
pi=2+i, rho=3+2i, sigma=4+i,
w=(1,pi,rho,pi sigma).
```

These are compatible divisors of one squarefree oriented conductor of norm
`1105`.  Their primitive consecutive blocks are

```text
2+i, 8+i, 33+4i.
```

The nontrivial interval contents are `5,65,65`, whose product is `21125`,
whereas the six pair residual ordinates are

```text
1,2,6,1,1,4,
```

with product `48`.  Neither divisibility nor an inequality comparing total
content directly to total residuals can hold.

### 7.2 Arbitrary-length continuant flexibility

Put

```text
f(k)=k^2-k+1,
z_k=f(k)+i=conjugate(k+i)(k-1+i).
```

Then, for `ell<k`,

```text
product_(j=ell+1)^k z_j
  = [product_(j=ell+1)^(k-1)(j^2+1)]
      ((k ell+1)+i(k-ell)).                         (7.1)
```

After primitive reduction, every interval residual is

```text
(k-ell)/gcd(k ell+1,k-ell) <= k-ell,               (7.2)
```

despite the enormous interval contents.  Cocycle and continuant algebra alone
therefore permit arbitrary length with only linear residual growth.

### 7.3 Finite character codes and height inefficiency

Every fixed split-prime binary character code can be realized by genuine
common-oriented Gaussian divisors in an arbitrarily small absolute sector.
One may choose

```text
w_j=T+d_j+i
```

with fixed CRT offsets `d_j`, add dummy split primes to avoid conjugate
orientation conflicts, and then choose `T` in a suitable residue class.  The
reduced pair residuals stay bounded as `T` grows, and all content cocycles are
honest.

What fails is precisely critical height efficiency.  If there are `n`
vectors and their pairwise gcds remain bounded, (6.1) gives

```text
H >= 2n log T-O(1),                                 (7.3)
```

while their angular span is only `Theta(T^(-2))`.  The critical scale is at
most `T^(-n/2+o(1))`, so the construction cannot be critical once `n>4`.

---

## 8. Consequence for the proof strategy

The exact all-interval identities close off a tempting but false route:
large scalar content, even cubic aggregate content, cannot by itself force an
extra residual penalty.  The same switch process pays for it through the edge
norms.  Quadratic carries rule out only a matching of near-complement pairs,
and lcm height rules out only weak-overlap realizations.

The surviving target is a **critical height-efficiency theorem** for a
balanced overlap code.  It must show that linearly many simultaneous
endpoint-small signed Gaussian products cannot reuse the same one-use
conductor atoms with optimal lcm height.  The companion note
`high_return_pattern_descent.md` proves that any hypothetical counterexample
is forced into exactly this dense, diffuse regime.
