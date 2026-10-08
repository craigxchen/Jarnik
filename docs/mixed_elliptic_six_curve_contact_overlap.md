# The contact-overlap obstruction on the first mixed elliptic six-curve

The curve in
[mixed_elliptic_six_curve_example.md](mixed_elliptic_six_curve_example.md)
has all twenty-five boundary pullbacks of degree one, but it cannot
support a full independent arithmetic contact profile.  The obstruction
is local and pointwise: six pairs of different compatible boundary
divisors have identical pullback divisors on the source.  Their exact
integral contact contents must consequently be equal up to a bounded
factor, whereas an independent cut profile makes their large cores
coprime.

The cleanest proof uses the five-row contact dictionary already proved
in [five_row_del_pezzo_arithmetic.md](five_row_del_pezzo_arithmetic.md),
Sections 2--3.  It needs neither an elliptic logarithm estimate nor a
Mordell--Weil argument.

## 1. The repeated source divisors

Write the six labels as

```text
infinity, 0, 1, a, b, c.
```

For a nonempty proper subset `A` of `{a,b,c}`, let `P_A` be the point
over `T=0` at which precisely the coordinates in `A` equal zero.  The
two clusters at this point are

```text
{0} union A,             {1} union A^c.
```

They are disjoint and hence define two compatible, distinct boundary
divisors.  Both contacts are transverse, and every individual boundary
pullback has degree one.  Therefore, as effective Cartier divisors on
the source elliptic curve `E`,

```text
f^* D_({0} union A)=f^* D_({1} union A^c)=[P_A].       (1)
```

There are six choices of `A`, explaining exactly the difference between
twenty-five boundary incidences and nineteen supporting source points.

For the established arithmetic dictionary, forget `a`.  The resulting
degree-two quotient `E'` is the smooth anticanonical five-point image
with coordinates

```text
(infinity,0,1,b,c).
```

The two points of `E` with `b=0,c=1` and the two choices of `a` map to
one point `R_01` of `E'`.  This point lies on the two boundary lines

```text
D_e,  e={0,b},             D_f,  f={1,c}.
```

Each line has intersection degree one with `E'`, and both intersections
are transverse at `R_01`.  Thus

```text
D_e|_(E')=D_f|_(E')=[R_01].                            (2)
```

The point with `b=1,c=0` similarly pairs `D_{1,b}` with `D_{0,c}`.
One pair is already enough for the contradiction.

## 2. Equal pullback divisors give commensurable contact contents

The precise local statement used here is the following.

**Lemma.**  Let `C/Q` be a fixed smooth projective curve, let
`phi:C -> X` be a fixed morphism, and let `D_1,D_2` be fixed effective
Cartier divisors on `X`.  Suppose

```text
phi^*D_1=phi^*D_2=Z                                      (3)
```

as effective Cartier divisors.  For `Q in C(Q)-supp(Z)`, let
`q_i(Q)` be positive integral contact contents obtained from fixed
integral projective models, so that `v_p(q_i(Q)) log p` is the normalized
finite local Weil function for `D_i` at `phi(Q)`, up to a fixed
`M_Q`-constant.  Then there are constants `c_p>=0`, zero for all but
finitely many primes, such that

```text
|v_p(q_1(Q))-v_p(q_2(Q))| log p <= c_p                 (4)
```

for every such `Q`.  Consequently

```text
log gcd(q_1(Q),q_2(Q))=log q_1(Q)+O(1)
                      =log q_2(Q)+O(1),                (5)
```

where the constants are independent of `Q`.

Indeed, by functoriality the two pulled-back local Weil functions are
Weil functions for the same divisor `Z`.  Uniqueness of local Weil
functions makes their difference bounded by an `M_Q`-constant.  After
spreading out the fixed data, the two equal generic-fiber Cartier
divisors have the same model away from finitely many primes, so the
bound is zero there.  At each remaining prime, uniqueness on the proper
curve gives a fixed bound; equivalently, the two arithmetic divisors
can differ only by a fixed vertical divisor.  This also handles bad
reduction primes and primes introduced by the chosen integral bases.
Finally,

```text
log q_1-log gcd(q_1,q_2)
 =sum_p max(v_p(q_1)-v_p(q_2),0) log p
 <=sum_p c_p,
```

which proves (5).

The hypothesis that the integers are the exact contact contents is
essential.  Equality of divisor heights alone, or independently chosen
large divisors of the contact integers, does not imply (4).

## 3. Conflict with the independent five-row cores

Use the notation of the five-row arithmetic dictionary.  For a boundary
line indexed by `e`, put

```text
d_e=n_e n_(e^c),             q_e=g_e/g,
B=product_(ten edges ij) b_ij.
```

The exact graph-content divisibilities are

```text
D_0 d_e | g_e,       D_0 | g,       g | D_0 B^2,
g_e | D_0 d_e B^2.                                    (6)
```

Writing `m=g/D_0`, equations (6) give

```text
m | B^2,       d_e | q_e m,       q_e m | d_e B^2,
q_e | d_e B^2.                                        (7)
```

For `e={0,b}` and `f={1,c}`, the four cuts
`e,e^c,f,f^c` are distinct.  In an independent profile their blocks
have disjoint prime support, so

```text
gcd(d_e,d_f)=1.                                        (8)
```

The last divisibility in (7), applied prime by prime with (8), yields
the particularly useful exact bound

```text
gcd(q_e,q_f) | B^2.                                    (9)
```

If the block scale is `w_5`, the full profile assumptions give

```text
log d_e=log d_f=2w_5+o(w_5),       log B=o(w_5).
```

Equation (8) of the five-row arithmetic note, or directly (6)--(7),
then gives

```text
log q_e=2w_5+o(w_5).                                  (10)
```

On the other hand, (2) and the lemma give

```text
log gcd(q_e,q_f)=log q_e+O(1)=2w_5+o(w_5),             (11)
```

while (9) gives `log gcd(q_e,q_f)=o(w_5)`.  This is a contradiction.

A full six-row independent profile restricts under forgetting `a` to
the required five-row profile: the two six-row blocks whose labels
differ only by membership of `a` are grouped together.  If each original
block has logarithmic norm `w+o(w)`, every grouped block has

```text
w_5=2w+o(w).
```

The grouped blocks retain disjoint prime support.  The full-profile
hypothesis must also include the already used condition that inherited
primitive-residue and correction contents have logarithm `o(w)`; this
is exactly what gives `log B=o(w_5)` above.

Thus the overlapping-contact curve is geometrically balanced but cannot
realize the full independent arithmetic profile.  The conclusion is
specific to this curve and its repeated boundary support.  It does not
exclude a balanced curve whose twenty-five boundary contacts occur at
twenty-five distinct source points, and it supplies no global uniform
lattice-arc bound.
