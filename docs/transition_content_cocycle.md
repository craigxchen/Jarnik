# The transition-content cocycle

## Status

This note records an exact identity coupling prime re-entry to the actual
Gaussian transition product. In the odd squarefree rooted model, the weighted
re-entry count on an interval is not merely a support statistic: it is exactly
the rational coordinate content removed from the product of its primitive
transition blocks.

The identity supplies the phase/conductor coupling missing from the abstract
width models. It does not by itself prove the uniform endpoint bound. Exact
five-transition data show that very large content can coexist with a final
primitive residual equal to one.

---

## 1. Rooted Gaussian divisors

Fix one oriented Gaussian prime `pi_p` above each prime in an odd squarefree
split conductor. Write

```text
w_i = unit_i product_p pi_p^(e_(i,p)),
e_(i,p) in {0,1},
d_i = Norm(w_i).
```

For two vertices put

```text
g_ij = Norm(gcd(w_i,w_j)),
B_ij = conjugate(w_i) w_j / g_ij
     = X_ij+i delta_ij,
q_ij = Norm(B_ij).
```

Then `B_ij` is coordinate-primitive away from the suppressed ramified factor,
and

```text
v_p(q_ij) = |e_(i,p)-e_(j,p)|.                     (1.1)
```

The sign of `e_(j,p)-e_(i,p)` selects `pi_p` or its conjugate in `B_ij`.

---

## 2. The three-vertex content cocycle

Direct multiplication gives

```text
B_ij B_jk = C_ijk B_ik,                            (2.1)
```

where

```text
C_ijk = d_j g_ik/(g_ij g_jk).                      (2.2)
```

This quotient is a positive integer. At a conductor prime,

```text
v_p(C_ijk)
  = (|e_i-e_j|+|e_j-e_k|-|e_i-e_k|)/2.             (2.3)
```

For binary exponents it is one precisely on the two re-entry words

```text
010, 101,
```

and zero on the other six words. Thus `C_ijk` is exactly the product of the
primes that switch twice on the triple.

Taking coordinates in (2.1) yields the exact continuant identities

```text
C_ijk delta_ik
  = X_ij delta_jk+delta_ij X_jk,                   (2.4)

C_ijk X_ik
  = X_ij X_jk-delta_ij delta_jk.                  (2.5)
```

Because `B_ik` is coordinate-primitive, `C_ijk` is also the ordinary rational
content of the product:

```text
C_ijk
  = gcd_Z(X_ij X_jk-delta_ij delta_jk,
          X_ij delta_jk+delta_ij X_jk).            (2.6)
```

This identity recovers the re-entry core directly from the Gaussian
coordinates themselves.

---

## 3. Arbitrary path intervals

For an ordered interval `a<b`, multiplication telescopes to

```text
product_(r=a)^(b-1) B_(r,r+1)
  = C_(a:b) B_ab,                                  (3.1)
```

with

```text
C_(a:b)
  = g_ab product_(a<r<b) d_r
      / product_(r=a)^(b-1) g_(r,r+1).             (3.2)
```

At one prime let

```text
T_p(a:b) = sum_(r=a)^(b-1) |e_(r+1,p)-e_(r,p)|.
```

Then

```text
v_p(C_(a:b))
  = (T_p(a:b)-|e_(b,p)-e_(a,p)|)/2.                (3.3)
```

The right side counts completed pairs of switches: every departure followed
by a return contributes one rational factor `p` to the content. Consequently

```text
C_(a:b)
  = content_Z(product_(r=a)^(b-1) B_(r,r+1)).      (3.4)
```

The norm identity accompanying (3.1) is

```text
product_(r=a)^(b-1) q_(r,r+1)
  = C_(a:b)^2 q_ab.                                (3.5)
```

Thus the conductor-weighted re-entry statistic, the rational content, and the
norm lost when a transition product is primitively reduced are the same
quantity.

---

## 4. Positivity and what it does not give

Inside an endpoint arc all consecutive transition arguments are positive and
their sum is small. Equation (3.1) preserves the full lifted phase:

```text
arg B_ab = sum_(r=a)^(b-1) arg B_(r,r+1).
```

Writing the product coordinates out gives a continuant whose common rational
divisor is exactly `C_(a:b)`. However, the primitive ordinate after division
is `delta_ab`; the content need not divide it or force it to be large. The
content is removed simultaneously from the real and imaginary coordinates.

The same phenomenon is visible in the known four-point cluster through the
exact identity

```text
(987+i)(843+i) = 610(1364+3i).
```

Thus a two-switch content `610` accompanies primitive residuals `1,1,3`.
More systematic exact product data are given next.

---

## 5. Exact five-transition calibration

For the six-point endpoint example, take the five primitive consecutive
blocks

```text
157+i, 993+i, 307+i, 302+i, 278+i.
```

The contents of the successive prefix products are

```text
1,
50,
18850,
11856650,
70487784250,
```

while their primitive ordinates are

```text
1, 23, 27, 17, 1.
```

In particular, the full product has exact rational content `70487784250` and
primitive ordinate `1`. The consecutive two-block contents are

```text
50, 650, 29, 145.
```

Therefore neither

```text
C_(a:b) <= product of primitive residuals
```

nor any divisibility of the content into those residuals can be true. A proof
must compare content and residual efficiency over many overlapping intervals,
not bound either one interval at a time.

---

## 6. Corrected next target

The exact remaining problem can now be phrased as follows.

> Prove that an unbounded positive transition path cannot have endpoint-small
> primitive ordinates on every relevant interval while the contents
> `C_(a:b)` simultaneously realize the re-entry formula (3.3) for one common
> squarefree Gaussian conductor.

The abstract diagonal countermodel does not satisfy (2.1)--(3.5) with its
assigned conductor characters, so this is genuinely stronger than support,
entropy, or projective-collision information. The finite example shows that a
valid theorem must be aggregate or asymptotic: four and five transitions can
remain maximally content-efficient.
