# Primitive content bound for the signed six-factor family

The later [ramified refinement](sixfactor_moving_content_c_half_exclusion.md)
improves the bound here to `|G|<=5 sqrt(2)|P|`. It excludes this entire
representative family at normalized arc constant `C<=1/2`, including
moving rational coefficients. The proof below remains the source of its
odd-prime bound; neither result treats arbitrary endpoint configurations.

The signed five-point construction in
[five_point_affine_shape_cubic_family.md](five_point_affine_shape_cubic_family.md)
has a uniform primitive-normalization bound even without the divisibility
assumptions used in Section 6 of that note.

Let

```text
(a_1,...,a_6)=(u,2u,v,u+v,2u+v,3u+v),
P=product_r a_r,
w_j=(-1)^|S_j| product_(r in S_j)(a_r+iT)
                    product_(r notin S_j)(a_r-iT),
```

where `u,v,T` are integers, `gcd(u,v,T)=1`, and `P!=0`.  Let `G` be
a Gaussian gcd of the five `w_j`.  Then

```text
G divides (1+i)^6 5 P,
N(G) <= 2^6 5^2 P^2.                                (1)
```

In particular `|G|<=40|P|`.  This estimate is deliberately not
optimized.  It is absolute and is enough to show that no choice of
signed coefficients, relative parameter rates, or rational denominator
clearing in this six-factor construction makes the normalized primitive
arc constant tend to zero.

## 1. Inert and ramified primes

Fix a rational prime `p`.  If `p=3 mod 4`, it remains prime in
`Z[i]`, and

```text
v_p(a+iT)=v_p(a-iT)=min(v_p(a),v_p(T)).               (2)
```

The exponent contributed to every row is therefore
`sum_r min(v_p(a_r),v_p(T))`, which is at most `v_p(P)`.  Its
contribution to `N(G)` is at most `p^(2v_p(P))`.

For the ramified prime, put `s=v_2(a)` and `t=v_2(T)`.  When `T!=0`,

```text
v_(1+i)(a plus-or-minus iT)
    =2 min(s,t) + indicator(s=t).                    (3)
```

Indeed, after removing the smaller rational power of two, the real and
imaginary parts have different parity unless `s=t`; in the latter case
both are odd and the remaining Gaussian valuation is exactly one.
The value in (3) does not depend on the sign.  Summing over the six
factors gives

```text
v_(1+i)(G) <= v_(1+i)(P)+6.                          (4)
```

If `T=0`, all five products are associates of `P`, so (1) is immediate.

## 2. Split primes dividing `T`

Now let `p>=5` split, and fix one Gaussian prime `pi` above `p`.  Put
`t=v_p(T)>0`.  Every oriented factor has the sign-independent base
valuation

```text
min(v_p(a_r),t).
```

The sum of these six base valuations is at most `v_p(P)`.  An additional
valuation is possible only for an index with `v_p(a_r)=t`; after the
base power of `p` is removed, exactly one of its two orientations can
vanish modulo `pi`.

Primitivity of `(u,v,T)` says that `u` and `v` are not both zero modulo
`p`.  If they are both units, the first three forms `u,2u,v` have zero
valuation, and at most one of `u+v,2u+v,3u+v` vanishes modulo `p`.
Thus at most one index can have an additional valuation, and one signed
support cannot cover all five rows.  If `p|u`, then `v` is a unit and
only indices `1,2` can be active.  Across the five rows their membership
pairs include all four choices `(0,0),(1,0),(0,1),(1,1)`, so for either
choice of their two active orientations some row avoids both.  If
`p|v`, only index `3` can be active.  Hence there is no common valuation
beyond the base in any case, and

```text
v_pi(G) <= v_p(P).                                   (5)
```

The same argument applies to the conjugate prime.

## 3. Split primes not dividing `T`

Suppose next that `p>=5` and `p` does not divide `T`.  Reduction modulo
`pi` identifies `iT` with a nonzero residue.  A factor can vanish only
when its coefficient form equals one of two opposite nonzero residues.
Tag such a form `+` or `-` according to the vanishing orientation.  Its
five-row support is `SUPPORTS[r]` for one tag and its complement for the
other.  A positive common `pi`-valuation requires the tagged supports to
cover all five rows.

The exact signed-cover calculation has only the following possibilities
for `p>=7`:

```text
condition       tagged indices and orientations
u=0             3+,4+,5+,6+   or   3-,4-,5-,6-
u+v=0           2-,6-         or   1+,3-,5+
2u+v=0          2-,3+         or   1+,4-,6+ .        (6)
```

This is a characteristic-free finite calculation.  Enumerate the
`3^6` signed tag vectors.  Every full-rank covering system has nonzero
three-by-three determinant content supported on `{2,3,5}`.  Closing the
rank-deficient covering systems under all six coefficient equations
gives exactly the six patterns in (6).  The accompanying checker performs
both integer calculations, rather than inferring (6) from samples over
finite fields.

The valuation in each boundary pattern is charged to `P`.  For `u=0`
modulo `p`, write the four active same-orientation factors as

```text
f_j=v+ju plus-or-minus iT,       j=0,1,2,3.
```

Put `e=v_p(u)`.  The difference of two such factors is `(j-k)u`, and
`j-k` is a unit for `p>=5`.  Thus at most one `f_j` has `pi`-valuation
larger than `e`.  Every row uses exactly two of these four factors, and
for each one of them some row uses two of the other three.  A row avoiding
the possible exceptional factor has total valuation at most `2e`.
Consequently

```text
v_pi(G)<=2e<=v_p(a_1a_2)<=v_p(P).                    (7)
```

For `u+v=0`, the two-factor pattern has factors `2-` and `6-`; rows
witness each factor alone, while their difference has coefficient
`a_6-a_2=u+v`.  The three-factor pattern has `1+`, `3-`, and `5+`;
again there are rows witnessing each alone, and the sum of the first two
oriented factors is `a_1+a_3=u+v`.  In either pattern the common
valuation is at most `v_p(u+v)`.  The two patterns for `2u+v=0` are
identical after using `a_2+a_3=2u+v` or `a_1+a_4=2u+v`.  This proves
(5) also when `p` does not divide `T`, apart from one exceptional
mod-five cover.

At `p=5`, exact enumeration of all residues adds only

```text
2+,4+,5-.                                            (8)
```

Rows `1`, `2`, and `3` respectively witness factors `5`, `4`, and `2`
alone.  If the common valuation is at least `h`, all three oriented
factors in (8) vanish modulo `pi^h`.  Subtracting factors `4+` and `2+`
gives `v=u mod pi^h`; adding factors `5-` and `2+` gives
`v=-4u mod pi^h`.  The cover (8) has `u` a unit modulo five, so
`5u=0 mod pi^h` forces `h<=1`.  Therefore for every split prime

```text
v_pi(G) <= v_p(P) + indicator(p=5).                  (9)
```

## 4. Global content and the arc constant

Multiplying the local bounds, remembering that an inert prime has norm
`p^2`, a split prime has two factors of norm `p`, and `N(1+i)=2`, gives

```text
N(G) <= 2^6 5^2 product_p p^(2v_p(P))
     = 1600 P^2,
```

which proves (1).

The exact chord identity from the earlier note implies, after primitive
division by `G`,

```text
max(|w_1-w_4|,|w_2-w_4|)/(|G| sqrt(|w_0|/|G|))
    >=2 sqrt(|P|/|G|).
```

Thus every arc containing the primitive five-tuple has

```text
arc length / sqrt(primitive radius) >=1/sqrt(10).    (10)
```

Finally, a common integer scaling of an integral triple multiplies both
`P` and `G` by the sixth power of that integer.  Given a rational triple,
clear a common denominator and divide the resulting integral triple by
`gcd(u,v,T)`.  Every integral clearing is a common integer multiple of
this primitive representative, so `|P|/|G|` is independent of the
clearing.  The primitive Gaussian tuple on that rational projective ray
is unique up to a unit. To verify the last assertion, two such tuples
differ by a Gaussian rational factor: any one nonzero coordinate gives
that factor. Bezout applied to the first primitive tuple makes the
factor Gaussian integral, and Bezout applied to the second makes its
inverse Gaussian integral. Thus it is a unit. This also treats common
Gaussian rational clearings of the five coordinates. Hence (10) covers
arbitrary common denominator clearings of this particular rational
six-factor family.

Run `python3 docs/check_five_point_signed_tangent.py` for the exact chord
identities, the signed-cover determinant calculation, the six rational
boundary closures, and the complete mod-five tag enumeration.
