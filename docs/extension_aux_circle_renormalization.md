# Extension auxiliary-circle renormalization

## Status

This note develops a global algebraic renormalization hidden in the existing
four-point certificate machinery.

The key identity is now formalized in
`GaussianChain/ExtensionAuxCircle.lean` and passes the full Lean build.

It does **not** yet prove the uniform endpoint bound.  It gives a genuine
circle-to-circle transformation and identifies the exact arithmetic parameter
that determines whether the transformation is a descent or a critical
self-similarity.

---

## 1. The discriminant is an auxiliary circle

For a geometric certificate `C`, the extension discriminant is a quadratic in
`W`:

```text
Delta(W) = -A^2 W^2 + L W + C0,
A = 2 y z V.
```

The negative-square leading coefficient is the formalized identity

```text
discLeading = -4 y^2 z^2 V^2.
```

If `Delta(W)=Y^2`, completing the square gives

```text
(2 A Y)^2 + (2 A^2 W - L)^2 = L^2 + 4 A^2 C0.       (1.1)
```

Thus every possible fifth-point extension maps to an integral point on one
fixed auxiliary circle.  This is not an approximation.

---

## 2. Integral slope reduction

The geometric certificate contains integers `n3` and `J`.  Suppose first that

```text
n3 = J r
```

for an integer `r`.

The certificate identities imply two exact equations:

```text
N V^2 - 2 h y z q = -2 y z V r,                     (2.1)

z (r^2+1) = 2 h N.                                  (2.2)
```

For an actual extension with parameters `(t,W)`, define

```text
X = W - h z r,
Y = h(2t-z) - r W.
```

Substituting (2.1) into the extension quadratic and using (2.2) gives

```text
X^2 + Y^2 = 2 h^3 N z.                              (2.3)
```

This is formalized as

```text
GeomCertificate.extension_reducedCircle
```

in `GaussianChain/ExtensionAuxCircle.lean`.

So every later point after the initial certificate becomes an integral point on
a second circle under one affine-linear transformation.

---

## 3. General rational slope

The integrality assumption is not essential conceptually.  Put

```text
g = gcd(n3,J),
r0 = n3/g,
q0 = J/g.
```

Then `gcd(r0,q0)=1`.  Multiplying the reduced coordinates by `J` gives integral
coordinates

```text
XJ = J W - h z n3,
YJ = h J(2t-z) - n3 W,
```

satisfying

```text
XJ^2 + YJ^2 = 2 h^3 N z J^2.                        (3.1)
```

The base point `(t,W)=(z,0)` is

```text
h z g (-r0,q0).
```

After removing this common scale, the primitive norm of the auxiliary root is

```text
N_new = r0^2+q0^2.
```

The certificate equations give the exact relation

```text
z N_new = 2 h N q0^2.                               (3.2)
```

Since `gcd(N_new,q0^2)=1`, equation (3.2) implies

```text
q0^2 divides z.
```

Writing

```text
z = lambda q0^2
```

therefore yields

```text
N_new = 2 h N / lambda.                             (3.3)
```

The integer `lambda` is the renormalization factor.

---

## 4. Primitive chord interpretation

The same parameter appears without the four-point certificate.

Let a rooted circle point have coordinates `(k,s)` in the standard rooted
parameter plane:

```text
k^2+s^2 = 2 h N k.
```

Put

```text
g = gcd(k,s),
k = g q,
s = g r,
gcd(q,r)=1.
```

The circle equation implies `q` divides `g`.  Write

```text
g = lambda q.
```

Then

```text
k = lambda q^2,
s = lambda q r,
q^2+r^2 = 2 h N/lambda.                             (4.1)
```

Thus `lambda` is also the primitive chord-norm parameter.  The auxiliary
circle is the primitive circle of the chord direction.

---

## 5. Descent versus critical self-similarity

Equation (3.3) gives a genuine primitive-norm descent whenever

```text
lambda > 2 h,
```

because then

```text
N_new < N.
```

The non-descent regime is

```text
lambda <= 2 h.                                      (5.1)
```

At first this appears finite, but `h` is not uniformly bounded.

For a primitive root (`h=1`), only two non-descent values remain:

```text
lambda=1 or lambda=2.
```

They are explicitly critical:

```text
lambda=1:  q^2+r^2 = 2N,

lambda=2:  q^2+r^2 = N.
```

The endpoint restriction gives `q=O_C(N^(1/4))`.  Hence each fixed lambda
places `(q,r)` in an axis-aligned square-root arc on the new circle.  Its
vertical coordinate varies through an interval of bounded length, so each
fixed lambda contributes only `O_C(1)` points.

This is a real partial theorem:

> For one fixed rooted point and one fixed renormalization factor `lambda`, the
> endpoint arc contains only `O_C(1)` extensions.

---

## 6. Why this does not yet finish the proof

Different cluster points may have different values of `lambda`.  Equation
(4.1) only shows

```text
lambda divides 2 h N,
```

and the number of possible divisors is not uniformly bounded.

Consequently, the transformation replaces one endpoint cluster by a union of
bounded fibers indexed by the arithmetic parameter `lambda`.  To finish the
proof one would need a global bound on the number of renormalization factors
represented inside one short arc.

That remaining problem is not artificial.  It is equivalent to controlling
primitive chord directions of many different norms in one tangent sector, the
same global compatibility that defeated the earlier local approaches.

The renormalization therefore gives:

```text
one lambda  -> uniformly bounded,
large cluster -> many distinct lambda values.
```

The next target is a theorem coupling distinct lambda values.  Promising
quantities are:

1. gcd and divisibility relations between the corresponding primitive chord
   Gaussian integers `q+i r`;
2. the fact that all `lambda(q^2+r^2)` equal the same integer `2hN`;
3. cyclic order of the original points, which imposes an order on the rational
   slopes `r/q`;
4. the affine extension map, which couples every later point to the same
   initial certificate rather than treating lambda fibers independently.

---

## 7. Verdict

The auxiliary-circle construction is not merely another discriminant bound.
It reveals an exact critical renormalization of the endpoint problem.

It succeeds in two respects:

- it converts all extensions of one certificate into points on one explicit
  integral circle;
- it gives a primitive-norm descent except for explicitly parameterized
  critical fibers.

It fails to prove uniformity because the number of critical fibers `lambda` is
not yet controlled.  Any continuation of this route should focus on
interactions between distinct renormalization factors, not on counting points
within one factor.
