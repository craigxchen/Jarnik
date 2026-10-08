# The common scalar discriminant fails for complex endpoint jets

The first-order divisor proof for a real common jet does not extend by
replacing the common derivative square with another polynomial of the
same degree. In the existing primitive five-row example with `d=3`,
**no polynomial `B` of degree below 20** has

```
y_j | 4 B P - P'^2       for all five rows.             (1)
```

The natural degree of a derivative square is `4d-2=10`. This obstruction
allows arbitrary `B`, not merely squares, and already concerns individual
divisibility; requiring a product with full root multiplicities cannot
repair it. It is an obstruction to this precise scalar ansatz, not to
every rotation-covariant differential construction or a uniform row bound.

## Why the radial constants stop being constant

Work in a constant real orthogonal coordinate system in which
`deg x_0=2d` and `deg y_0=q<2d`. Suppose `q>d` and all row differences
have degree at most `d`. Write `x_j=x_0+a_j`, `y_j=y_0+b_j`. Equal norms give

```
2x_0 a_j + 2y_0 b_j + a_j^2 + b_j^2 = 0.
```

If `deg a_j>q-d`, the first term has degree greater than `q+d`, while
the other terms have degrees at most `q+d, 2d, 2d`. Therefore

```
deg a_j <= q-d.                                      (2)
```

For `q<=d`, the same argument gives degree at most zero. This explains
the constant shifts in the real-jet proof. For general `q`, the shifts
are polynomials. With `H=x_0`, differentiating the actual norm identity
and reducing modulo the full polynomial `y_j` gives the exact congruence

```
4H'^2 P-P'^2 = -4P a_j'(2H'+a_j')  (mod y_j).        (3)
```

No squarefree replacement is used in (3). The new derivative term is
the obstruction even before shared-root multiplicity accounting.

## Exact primitive fixture and CRT certificate

Use the unscaled five rows from
[the cubic affine-shape family](five_point_affine_shape_cubic_family.md):

```
x_j=-T^6-45T^4-c_j T^2+e_j,
y_j=-T^5-d_j T^3-f_j T,
P=product_(a=1)^6 (T^2+a^2),

j    c_j    e_j    d_j    f_j
0    404     720     85    1164
1    428    -720     61     396
2    464    -720     25     324
3    484    -720      5    -276
4    524     720    -35   -1236.
```

These have common norm `P`, degree six, and difference degree at most
three. Their common Gaussian polynomial gcd is one: in the factorization
by `a+iT` and `a-iT`, every factor is omitted in at least one row. Here
`q=5` and the radial degree two reaches (2).

Set `U=T^2` and `m_j(U)=U^2+d_j U+f_j`. The five quadratics are pairwise
coprime, have nonzero constant terms, and satisfy `m_j(-a^2)!=0` for
`a=1,...,6`. Consequently `P` is invertible modulo every `y_j`.
The differentiated norm identity makes (1) equivalent to

```
B(T) = x_j'(T)^2  (mod y_j(T)).                     (4)
```

In particular `B(0)=0`. Modulo `m_j(T^2)`, (4) requires

```
B(T) = 4U(3U^2+90U+c_j)^2  (mod m_j(U)).            (5)
```

Let `M(U)=product_j m_j(U)`. The unique CRT remainder `b(U)` of degree
below ten for the five congruences (5) has

```
deg b=9,
[U^9]b=389/13440000,
b(0)=370384723467/4375,
M(0)=50947247932416.
```

The unique polynomial of degree below eleven satisfying those
congruences and the additional condition at zero is

```
C(U)=b(U)-b(0) M(U)/M(0),
deg C=10,          [U^10]C=-67/40320000.             (6)
```

Thus `C(T^2)` satisfies (4) and has degree twenty. The full least common
multiple of the `y_j` is, up to a nonzero scalar,

```
L(T)=T M(T^2),                 deg L=21.
```

Every solution to (4) differs from `C(T^2)` by a multiple of `L`.
The unique remainder of degree below 21 therefore has degree twenty;
no solution has smaller degree. This also rules out odd polynomial
terms as an escape: the CRT statement is in the full ring `Q[T]`, not
only its even subring.

The [exact standard-library checker](check_complex_jet_scalar_discriminant_obstruction.py)
verifies the row norms and degrees, CRT inverses (hence coprimality),
the displayed rational coefficients, and divisibility by every complete
`y_j`. Run it with `python3 docs/check_complex_jet_scalar_discriminant_obstruction.py`.

The minimum degree twenty is for the fixed leading-real coordinate
system above. It does not give a degree-independent upper bound for
the number of general complex-jet rows.

## Common rational rotation and exact primitive normalization

A moving rational rotation cannot straighten a primitive polynomial
tuple for free. In fact, let `gcd_j z_j=1` in `Q(i)[T]`, let
`U` be a rational function of norm one, and let a rational function `F`
clear its denominators so that every `F U z_j` is polynomial. Polynomial
Bezout gives `sum_j h_j z_j=1`, hence

```text
B=F U=sum_j h_j(F U z_j) is a polynomial.
```

The cleared rows are exactly `B z_j`, whose common polynomial gcd is
`B` up to a constant. Their raw norm degree is `4d+2 deg B`; removing
this exact content restores the original tuple and degree `4d`, up to
one common nonzero constant. In particular, if `U z_j` were already
polynomial for all rows, then `U` would be polynomial; `U bar U=1`
would force it to be constant. This applies to rational common
similarities as well, since the Bezout argument itself does not use
the norm-one hypothesis.

The displayed fixture's coefficients of `T^6` and `T^5` are `-1` and
`-i`. No common nonzero constant can make both real. Thus neither a
constant change of coordinates nor a rational rotation followed by
full polynomial content removal supplies the real-jet hypothesis for
that primitive tuple. A differential calculation in a rational frame
would have to retain its denominators explicitly.

This rigidity does not apply to general nonconformal maps of the
half-angle parameter. For example, the primitive pair
`z_0=T^2-i, z_1=T^2+i` has common norm `T^4+1`, and its relative
half-angle parameters are `0,1/T^2`. The map `t -> T^2 t` takes these
to `0,1`; their primitive realization is the constant pair `1-i,1+i`
of norm two. Thus a degree-monotonicity claim for arbitrary projective
maps would already fail on two rows. No construction retaining
arbitrarily many rows in the real-jet class has been proved.

Sol independently checked the scalar CRT by rational linear algebra
and expanded all five factor products. The permanent checker now
links those exact expansions to the displayed rows, so its primitive
factor-support verification does not rely on an untested identification.
