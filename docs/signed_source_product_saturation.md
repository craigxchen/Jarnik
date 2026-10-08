# Signed source products recover the full contact-jet lattice

The product-height barrier in
[higher_degree_gram_jet_barrier.md](higher_degree_gram_jet_barrier.md)
does not reflect a large missing arithmetic index. Under the median
source incidence, the integer span of all canonical source products
recovers the full jet lattice after multiplication by a power of the
three small anchor determinants. This is an exact saturation result.
It gives no coefficient-height bound for the signed representations
and does not exclude a first short frame.

## 1. Fixed arithmetic data

Let `v_1,...,v_m` be primitive integer binary vectors, with the first
three pairwise independent. For every contact `e`, fix a nonempty
cluster `C_e` containing at most one of the three anchors, a primitive
representative `lambda_e` in that cluster, and a positive integer
`d_e`. Assume the `d_e` are pairwise coprime and

```text
d_e divides det(v_i,lambda_e) for every i in C_e.      (1)
```

These are ordinary full-depth congruences. No hypothesis is made that
outside rows remain distinct modulo primes of `d_e`.

Put

```text
L_i(X,Y)=det(v_i,(X,Y)),
I_e=(d_e,L_(lambda_e)) subset Z[X,Y],
Gamma_(d,h)=Z[X,Y]_d intersect intersection_e I_e^h,
delta=lcm(|det(v_1,v_2)|,|det(v_1,v_3)|,
          |det(v_2,v_3)|).
```

For every nonnegative integer multiplicity vector `n` with `sum n_i=d`,
set

```text
n(C_e)=sum_(i in C_e)n_i,
P_n= [product_e d_e^max(h-n(C_e),0)] product_i L_i^n_i.
```

Let `Lambda_(d,h)` be the integer span of **all** these polynomials.
Both degree `d` and jet order `h` are fixed nonnegative integers. The
basic divisibility calculation gives `Lambda_(d,h) subset Gamma_(d,h)`.

## 2. Exact saturation theorem

For the data above,

```text
delta^d Gamma_(d,h) subset Lambda_(d,h)
                    subset Gamma_(d,h).              (2)
```

In particular

```text
[Gamma_(d,h):Lambda_(d,h)] <= delta^(d(d+1)).          (3)
```

All conclusions remain valid with private contacts omitted. They use
only the displayed cluster and anchor hypotheses.

**Proof.** Work first over `Z_p`. If `p` divides no `d_e`, all the
multipliers in `P_n` are units. Products of any two independent anchor
linear forms span a sublattice containing `delta^d Z_p[X,Y]_d`, by
clearing their inverse linear coordinate matrix in each of the `d`
factors.

Suppose now `p^a || d_e`. Every other contact norm is a unit in
`Z_p`. Choose a primitive inside row `v_i`, and choose two anchors
`v_j,v_k` outside `C_e`, which is possible because the cluster contains
at most one anchor. The determinant identity

```text
det(v_j,v_k) v_i
 =det(v_i,v_k) v_j-det(v_i,v_j) v_k
```

and primitivity of `v_i` imply

```text
min(v_p det(v_i,v_j),v_p det(v_i,v_k))
 <= v_p det(v_j,v_k) <= v_p delta.                    (4)
```

Choose an outside anchor attaining the first inequality. Since `v_i`
and `lambda_e` are primitive and congruent projectively modulo `p^a`,

```text
I_e Z_p[X,Y]=(p^a,L_i).
```

Use a unimodular `Z_p` coordinate change with `L_i=Y`. The selected
outside anchor form becomes

```text
L_j=uX+vY,       v_p(u)<=v_p(delta).
```

A basis of the local degree-`d` jet lattice is

```text
p^(a max(h-k,0)) X^(d-k)Y^k,      0<=k<=d.            (5)
```

The canonical products supported only on rows `i,j` give the same
expressions with `X` replaced by `L_j`, up to units from all the other
contact norms. Multiply (5) by `u^d` and expand:

```text
u^d p^(a max(h-k,0)) X^(d-k)Y^k
 =u^k p^(a max(h-k,0)) (L_j-vY)^(d-k)Y^k.
```

Every resulting term has `Y` exponent `l>=k`. Its integer multiplier
is divisible by the required factor `p^(a max(h-l,0))`, because
`max(h-k,0)>=max(h-l,0)`. It belongs to the local canonical span.
Since `delta^d/u^d` is integral over `Z_p`, this proves (2) locally.

The local inclusions at every prime imply the integral inclusion in
(2). The coefficient lattice has rank `d+1`, so its quotient is killed
by `delta^d` and (3) follows. The proof retains all exponents `a`; it
never discards a prime merely because an outside row collides there.

## 3. Full-profile consequences and the absence of a height bound

In the median frame, all three anchor determinants have absolute value
at most the primitive residue bound `T`. Thus

```text
delta<=T^3,
log[Gamma_(d,h):Lambda_(d,h)] <=3d(d+1)log T=o(w)     (6)
```

for fixed `d` in a full profile with `log T=o(w)`. Consequently the
source-product lattice has the same leading determinant exponent as
the full jet lattice. The large individual product heights cannot be
attributed to an exponentially large saturation defect.

For an admissible Gram frame and `d>=h>=1`, apply the exact evaluation
map `epsilon_d(P)=P(w_0,-z)` from the higher-degree note. It gives

```text
delta^d (Delta^h) subset epsilon_d(Lambda_(d,h))
                       subset (Delta^h).             (7)
```

The Gaussian image has index at most `delta^(2d)` in `(Delta^h)`.
Thus there are signed integer source-product representations of both
`delta^d Delta^h` and `i delta^d Delta^h` after evaluation. Neither (2)
nor (7) controls the integer coefficients of these representations or
the coefficient height of a lift with that evaluation.

The Gram kernel itself illustrates why signed cancellation cannot be
excluded merely from the product-height gap. Since `Q_U` belongs to
every `I_e`,

```text
delta^(2h) Q_U^h in Lambda_(2h,h).                   (8)
```

Its full-profile coefficient height is
`exp(2h A_m w+o(w))`. The canonical products can individually be much
larger; the representation (8) then necessarily uses cancellation.
It gives a kernel polynomial, not the nonkernel witness needed for a
first-frame contradiction. Saturation therefore identifies signed
combinations as the same underlying jet-lattice problem, up to a
subpower factor, rather than an independent arithmetic shortcut.

## 4. Even the quotient determinant scale remains above the target

There is a useful exact comparison after removing the Gram kernel.
Use the coefficient norm

```text
||sum p_k X^(d-k)Y^k||_B^2=sum p_k^2/binomial(d,k),
H=||U||_F.
```

Let `K_d=Q_U R[X,Y]_(d-2)` and let `pr` be orthogonal projection onto
its two-dimensional orthogonal complement. The real Jacobian of
`epsilon_d` on that complement is

```text
J_d=1/2 sqrt(H^(4d)-(H^4-4)^d)
   ~ sqrt(d) H^(2d-2)                                (9)
```

for fixed `d` as `H` grows. To verify it, in orthonormal coefficient
coordinates the complex evaluation vector has squared Hermitian norm
`H^(2d)` and complex bilinear square `(w_0^2+z^2)^d`. Also

```text
|w_0^2+z^2|^2=H^4-4
```

because `Im(bar(z)w_0)=1`. The Gram determinant of the two real evaluation
rows is therefore the square of (9).

Since the evaluation image of `Gamma_(d,h)` is exactly `(Delta^h)`,
whose real covolume is `D^h`,

```text
covol(pr Gamma_(d,h))=D^h/J_d.                       (10)
```

The corresponding covolume for the canonical span is between this
value and `delta^(2d)` times it, by (7).

The two-dimensional geometric-mean scale of (10) is consequently

```text
~ d^(-1/4) D^(h/2)/H^(d-1).
```

It is still a factor of order `H` above the nonkernel ideal threshold
`D^(h/2)/H^d`. In the full profile that difference has logarithm
`A_m w+o(w)`. This is only a quotient determinant-scale comparison.
It is not a lower bound on the shortest projected vector, and lifting
a projected vector back into the integral coefficient lattice also
requires controlling the Gram-kernel component.

## 5. Exact certificate

Run `python3 docs/check_signed_source_product_saturation.py`. The
checker constructs the integer spans for all fifteen median clusters
of five fixed primitive source rows and directly verifies the inclusion
`delta^d Gamma subset Lambda`. It includes a contact of depth two and
depth five at 5, while outside rows coincide at depth one; the observed
saturation defect is unchanged by the increase in prescribed depth.
It also checks 400 exact real Gram determinants of the Bombieri
evaluation map against (9). These examples verify the stated lattice
identities, not asymptotic full-profile heights or coefficient bounds
for signed representations.
