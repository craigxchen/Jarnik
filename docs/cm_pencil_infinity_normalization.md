# Common-root-at-infinity normalization for the CM pencil

This is a bounded algebraic normalization and search aid.  It does not
assert a global nonexistence result.

Start with

```text
f0(x)=x^2/(x-1),
v=x/(x-2).
```

Scaling `v` by a rational `a` and returning to the `x` coordinate gives

```text
phi_a(x)=a x/[1+(a-1)x/2].
```

The corresponding target map is

```text
M_a(y)=a^2 y/[1+(a^2-1)y/4],
M_a^(-1)(y)=4y/[4a^2-(a^2-1)y].
```

Direct substitution gives

```text
f0(phi_a(x))=M_a(f0(x)).
```

Therefore the transformed pencil is

```text
f_new(x)=M_a^(-1)(f(phi_a(x))),
```

and `f0` is fixed.  The map `M_a` fixes the branch values `0` and `4`,
although it generally moves the other branch values.  Since `phi_a` and
`M_a` are rational Möbius maps, rationality, distinctness, ramification,
and contact multiplicities are preserved on the projective lines.
Being a pole is a chart-dependent property and is not preserved by a
target Möbius map.

For a finite rational common argument `rho` with `rho != 0,2`, choose

```text
a=rho/(rho-2).
```

Then

```text
phi_a(infinity)=2a/(a-1)=rho,
```

so the old common argument `rho` is the new common argument at infinity.
The exclusions `rho=0` and `rho=2` are the two critical points of `f0`.
The normalized common argument at infinity has common target infinity:
it is a simple pole of `f0`, hence an unramified projective contact.
One must retain this contact rather than discard it as an affine pole.

After this normalization the common-root cross-product has homogeneous
degree four and a simple factor at infinity, so its affine degree is
three.  All denominator quadratics have zero leading coefficient.
Subtracting a multiple of the base pair removes the other numerator's
quadratic coefficient; scaling the remaining non-base pair then gives

```text
P1=bx+c,  Q1=x+e,
```

and hence

```text
x^2 Q1-P1(x-1)
 =x^3+(e-b)x^2+(b-c)x+c
 =product_(j=1)^3 (x-r_j).
```

Writing `S1=sum r_j` and `S2=sum_(i<j) r_i r_j` gives

```text
c=-r1r2r3,  b=S2+c,  e=b-S1.
```

For

```text
f_k=(x^2+k(bx+c))/(x-1+k(x+e)),
```

the exact branch discriminant is

```text
Delta_k(y)=(1+k)^2 y^2
 +[-2kb(1+k)-4(1-ke)]y+k^2b^2-4kc.
```

Consequently the parameters whose pencils have branch value `0` and `4`
are

```text
lambda_b=4c/b^2,
lambda_c=4(2b-8+c-4e)/(b-4)^2.
```

Put `H=2b-8+c-4e`.  The other branch values at those two parameters are

```text
s=4b[b^3+2cb(b-2e)+8c^2]/(b^2+4c)^2,

t=8H[b^2(b-4-2e)+4c(b-2)]
  /[b^2+4c-16(e+1)]^2.
```

The unordered branch set is therefore `{0,4,s,t}`.  Its harmonic tests are

```text
st-2s-2t=0,
st+4s-8t=0,
st-8s+4t=0.
```

The exact checker applies only necessary rational denominator filters before
testing these equations.  A surviving harmonic tuple would still require:

```text
* lambda_b and lambda_c finite, distinct, and away from denominator drops;
* the corresponding numerator/denominator pairs to be coprime;
* both maps to have degree exactly two;
* 0,4,s,t to be distinct finite branch values;
* the common contacts and relevant fibers to be unramified and to have
  the required multiplicities, using projective charts at poles.
```

Run

```text
python3 docs/check_cm_pencil_infinity_normalization.py
```

The reproduced necessary-filter search has no harmonic hits for 3,617
integer triples `r_j` in `[-15,15]` (excluding `0,2`) and no hits for
186,481 distinct rational triples with numerators bounded by 15 and
denominators at most 5 (again excluding `0,2`).  The partial degree and
branch-distinctness checks
are reported separately and also have zero hits in these searches.  These
bounded results do not claim that all rational triples have been tested.

There is also a residual symmetry preserving the normalized common
argument at infinity.  In the `v` coordinate it is `v -> 1/v`; in the
source and target coordinates it is respectively `x -> 2-x` and
`y -> 4-y`.  Its exact action is

```text
(b,c,e) -> (4-b, c+2b-4e-8, -2-e),
(r1,r2,r3) -> (2-r1,2-r2,2-r3),
(lambda_b,lambda_c) -> (lambda_c,lambda_b),
(s,t) -> (4-t,4-s).
```

This is an involution and preserves the unordered harmonic condition.
It can identify equivalent candidates, but supplies no independent CM
equation.  The general transformation `v -> a/v` uses the same one
continuous parameter as scaling; after fixing a common root at infinity,
only a finite residual symmetry remains.  There is no second free
normalization of a generic common-root set.

An independent symbolic audit checked the normalizer identity, the full
discriminant, both remaining-branch formulas, and this residual
involution.  The audit also corrected the affine-pole wording above.
