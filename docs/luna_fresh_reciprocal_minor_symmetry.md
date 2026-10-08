# Reciprocal minors and the exact reanchoring obstruction

This note records a backwards calculation for four points in a one-sided
endpoint arc. It gives an integer lower bound for the norm of the Gaussian
least common multiple, obtained from reciprocal chord
minors. The calculation is exact, including arbitrary gcds. It also explains
why iterating the usual anchor reflection cannot create a smaller-radius
tuple: on the full labelled tuple it is only conjugation followed by a common
rational rotation.

The result does not prove the uniform point bound. The missing step is a
lower bound for the lcm in terms of the smallest cotangent; the minor bound
below is a concrete necessary condition that any such argument must improve.

## 1. Reciprocal minor lemma

Let `d>=1` and let

```text
t_0<t_1<...<t_{k-1},       w_j=-d+i t_j,
A=lcm_G(w_0,...,w_{k-1}),       q_j=A/w_j.
```

Here `lcm_G` is a Gaussian least common multiple, chosen up to a unit, and
`q_j` are Gaussian integers. Put `n_j=d^2+t_j^2` and, for `i<j`, reduce

```text
d(t_j-t_i)/(n_i n_j) = a_ij/b_ij,       gcd(a_ij,b_ij)=1,
```

with `b_ij>0`. Then

```text
B_k := lcm_(i<j) b_ij  divides  N(A).                   (1)
```

Consequently

```text
|A| >= sqrt(B_k),
|A| >= max_(i<j) sqrt(n_i n_j/(d(t_j-t_i))).            (2)
```

The last inequality is the individual reciprocal-minor estimate; the lcm in
(1) retains all `binom(k,2)` constraints simultaneously.

### Proof

For two indices, take the ordinary integer determinant of the two reciprocal
vectors:

```text
delta_ij = det(q_i,q_j) = Im(conjugate(q_i) q_j).
```

Since `q_i=A/w_i`,

```text
w_i conjugate(w_j) = d^2+t_i t_j+i d(t_j-t_i),
```

and therefore

```text
delta_ij = N(A)d(t_j-t_i)/(n_i n_j).                   (3)
```

The determinant is a nonzero integer: the `q_i` and `q_j` cannot be real
collinear, since their quotient `w_j/w_i` is nonreal when `t_i!=t_j`.
Writing the fraction in lowest terms in (3), integrality gives
`b_ij | N(A)`. Taking the lcm proves (1). Since
`b_ij = n_i n_j/gcd(n_i n_j,d(t_j-t_i))`, (2) follows from
`b_ij >= n_i n_j/(d(t_j-t_i))` and `N(A)=|A|^2`.

No coprimality of the `w_j`, no restriction on prime powers, and no
primitivity of the complementary quotients `q_j` was used.

## 2. What the lcm bound does and does not control

For an endpoint cotangent chart, `t_j/d` is the cotangent parameter. Thus
the quadratic height target for four finite points would ask for a lower bound
of the shape

```text
|A| >= c d (t_0/d)^2 = c t_0^2/d.                     (4)
```

The reciprocal-minor lemma reduces this to the explicit arithmetic question

```text
B_4 >= c^2 t_0^4/d^2.                                 (5)
```

Equation (5) is not proved here. The gap is real: the proof of (1) only uses
that the six reciprocal minors are integers, and their denominators can share
substantial factors. Replacing the lcm by the product asserts information
that is not present in (3).

For orientation, exact calculations give the following values. In each line
`B_4` is the lcm in (1), while `N(A)` is the norm of the actual Gaussian lcm.

```text
d=2, t=(6,8,16,42):   B_4=8840,   N(A)=8840,
d=1, t=(3,5,13,31):  B_4=81770,  N(A)=81770,
d=2, t=(6,8,10,36):  B_4=44200,  N(A)=44200.
```

The first line is the bounded four-factor endpoint configuration used in the
vertical-cofactor checks. Equality here shows that the new minor lcm can be
sharp; it is not a disposable denominator estimate. For 30,000 random
ordered quadruples in total, with `1<=d<=30`, direct integer Gaussian lcm
computation gave `B_4 | N(A)` in every case. This check is only a sanity
check; (1) is proved by (3). An exact checker is saved as
`check_luna_fresh_reciprocal_minor_symmetry.py`.

## 3. Exact backwards reflection and its limitation

There is a related symmetry of the rational unit circle. Let
`q_j=z_j/z_0` be a tuple of distinct rational unit-circle points, represented
by a primitive Gaussian integral tuple `z_0,...,z_{k-1}` of radius `R`. If
`q_a` is chosen as a new anchor, the standard divisor-complement reflection
acts by

```text
q_j  ->  q_a/q_j.                                      (6)
```

After reanchoring at the image of the old anchor, the new ratios are

```text
(q_a/q_j)/(q_a/q_0) = 1/q_j = conjugate(q_j).          (7)
```

Thus the full labelled tuple is the complex conjugate of the old tuple,
followed by the common rational rotation `q_a`. Conjugation preserves
Gaussian primitivity and the radius. A common rational Gaussian multiplier
only changes the chosen projective representative; after primitive
normalization it gives the same radius. The normalized angular width is also
unchanged: if the old lifted angles are `theta_j`, (6) has angles
`theta_a-theta_j`, whose width is the old width.

This remains true after any finite sequence of reanchorings: every step
changes the projective tuple by conjugation and a common rotation, so no
backwards iteration can produce a radius descent. The exact transformation is
still useful because it transports every pair-minor condition and every
prime-power denominator condition to the reflected chart. To close the
endpoint problem one would need to combine these transported minors with a
new cross-chart inequality; (1) by itself leaves precisely the lcm-height
gap (5).

## 4. Remaining gap

The lemma supplies a different necessary condition from direct chord products:
all reciprocal pair minors are integral after one common normalization, and
their reduced denominators divide one integer `N(A)`. What remains is to
show that four (or another fixed number of) positive cotangents force `B_k`
to have fourth-power size in the smallest cotangent, or to find an explicit
family violating that inequality while still realizing the full Gaussian
endpoint tuple. The three-finite-coordinate Pell family only tests `k=3`; it
does not decide (5) for four finite coordinates.
