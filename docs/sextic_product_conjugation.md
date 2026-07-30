# Conjugation and rational normalization of the three sextic products

## 1. The three products on the norm-one torus

Work in the affine chart

```text
[1:x:y],        x conjugate(x)=y conjugate(y)=1.
```

For `b,c >= 0` with `b+c <= 6`, put

```text
m_{b,c}=x^b y^c.
```

The three explicit star-basis products are

```text
P_0(x,y)
  = product_{(b,c)!=(0,0)} (x^b y^c-1),

P_1(x,y)
  = product_{(b,c)!=(6,0)} (x^b y^c-x^6),

P_2(x,y)
  = product_{(b,c)!=(0,6)} (x^b y^c-y^6),
```

where every product ranges over the 28 pairs `b,c >= 0`, `b+c <= 6`, with
one root omitted.  Each product has 27 factors.

The exponent sums over all degree-six monomials are

```text
sum b = sum c = 56.
```

---

## 2. Exact conjugation formulas

For the first product,

```text
conjugate(x^b y^c-1)
  = x^(-b)y^(-c)-1
  = -(x^b y^c-1)/(x^b y^c).
```

There are 27 factors, so

```text
conjugate(P_0) = -P_0 x^(-56)y^(-56).                 (2.1)
```

For `P_1`,

```text
conjugate(x^b y^c-x^6)
  = x^(-b)y^(-c)-x^(-6)
  = -(x^b y^c-x^6)/(x^(b+6)y^c).
```

After omitting `(b,c)=(6,0)`,

```text
sum (b+6)=212,       sum c=56.
```

Hence

```text
conjugate(P_1) = -P_1 x^(-212)y^(-56).                (2.2)
```

Similarly,

```text
conjugate(P_2) = -P_2 x^(-56)y^(-212).                (2.3)
```

Equivalently,

```text
P_0/conjugate(P_0) = -x^56 y^56,
P_1/conjugate(P_1) = -x^212 y^56,
P_2/conjugate(P_2) = -x^56 y^212.                     (2.4)
```

---

## 3. Rational normalization

All exponents in (2.4) are even.  Define

```text
R_0 = P_0 x^(-28)y^(-28),
R_1 = P_1 x^(-106)y^(-28),
R_2 = P_2 x^(-28)y^(-106).
```

Then

```text
conjugate(R_i)=-R_i.                                  (3.1)
```

Since `R_i in Q(i)` and the anti-fixed subspace of conjugation on `Q(i)` is
`i Q`, there are rational numbers `T_i in Q` such that

```text
R_i=i T_i.                                             (3.2)
```

Thus every product-section value splits into

```text
P_0 = i T_0 x^28 y^28,
P_1 = i T_1 x^106 y^28,
P_2 = i T_2 x^28 y^106.                               (3.3)
```

The point-dependent nonarchimedean behavior is therefore the sum of:

1. an explicit linear function of the valuations of `x` and `y`; and
2. valuations of three rational numbers `T_0,T_1,T_2`.

This is much narrower than a general family of `S`-unit linear forms.

---

## 4. Conjugate-prime pairing

At a split rational prime `p=pi conjugate(pi)`, write

```text
a=v_pi(x),        b=v_pi(y).
```

If `A_i=v_pi(P_i)` and `B_i=v_conjugate(pi)(P_i)`, then (2.4) gives

```text
A_1-B_1-(A_0-B_0)=156a,
A_2-B_2-(A_0-B_0)=156b.                               (4.1)
```

Equivalently, after the normalization (3.3),

```text
v_pi(T_i)=v_conjugate(pi)(T_i)=v_p(T_i).               (4.2)
```

The two orientations of a split prime therefore share the same residual
rational valuations.  All orientation dependence is carried by the explicit
monomial exponents in (3.3).

If the tropical chamber at `pi` selects coordinate `i` and the conjugate
chamber selects coordinate `j`, the paired product-section valuation is

```text
v_pi(P_i)+v_conjugate(pi)(P_j)
```

and (3.3) rewrites it as

```text
v_p(T_i)+v_p(T_j) + an explicit linear conductor term. (4.3)
```

Thus the remaining max problem is a three-rational-number problem over the
rational primes, with the conductor term separated exactly.

---

## 5. Why this may help

The unrestricted uniform-in-`S` statement failed because the exceptional
behavior could be hidden independently at arbitrarily many Gaussian prime
places.  Formula (4.2) removes that freedom for the residual sextic products:
prime and conjugate-prime contributions are forced to be equal and come from
ordinary rational valuations of only three numbers.

A successful completion would prove a conductor-relative bound for the
selected sum in (4.3), using simultaneously:

- the rational product formula for `T_i`;
- the explicit archimedean products of sines represented by `T_i`;
- the fact that the selected pair `(i,j)` is the pair of extreme coordinate
  valuations, so the omitted index is the tropical middle coordinate.

The missing estimate is now a rational three-number inequality rather than a
moving family of local adapted bases.
