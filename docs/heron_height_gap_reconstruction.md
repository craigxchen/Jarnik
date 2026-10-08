# Heron coherence and the height gap between its two branches

This note studies the simultaneous arithmetic of the spin equations in
[graph_trace_and_spin_defects.md](graph_trace_and_spin_defects.md).
There is a positive reconstruction result: in the oriented core height
range, each anchor Heron equation forces the exact primitive Gaussian
pair numerator, including its entire gcd. All remaining triangle and
four-index identities then follow. An exact four-point example shows
why the unsquared branch information cannot simply be omitted outside
that height range.

This does not prove a positive lower exponent for a residue. It reduces
one proposed route to the existence of the simultaneous anchor Heron
solutions themselves, rather than additional identities among their
correctly reconstructed pairs.

## 1. The normalized Heron equation has exactly two branches

Let

```text
h_i=x_i+i t_i,       h_j=x_j+i t_j,       b=y+i s
```

be conjugate-primitive Gaussian integers with all six coordinates
positive. Assume their arguments `alpha,beta,gamma` lie in `(0,pi/4)`.
In fact only `alpha+beta<pi/2` and `0<gamma<pi/2` are needed.
Write

```text
q_i=N(h_i),       q_j=N(h_j),       q_b=N(b),
A=t_i/sqrt(q_i),  B=t_j/sqrt(q_j), C=s/sqrt(q_b).
```

The symmetric normalized Heron equation is

```text
A^4+B^4+C^4-2A^2 B^2-2A^2 C^2-2B^2 C^2
                 +4A^2 B^2 C^2=0.                         (1)
```

It involves only rational numbers `A^2,B^2,C^2`. Equivalently,

```text
(A^2+B^2-C^2)^2=4A^2 B^2(1-C^2).
```

Regarding (1) as a quadratic in `C^2`, its two roots are exactly

```text
(A sqrt(1-B^2)+B sqrt(1-A^2))^2,
(A sqrt(1-B^2)-B sqrt(1-A^2))^2.
```

The angle restrictions therefore prove the exact alternative

```text
gamma=alpha+beta       or       gamma=|alpha-beta|.           (2)
```

No approximation or field-degree assumption enters this step. If the
two anchors are distinct and primitive with positive coordinates,
their arguments differ: two primitive Gaussian integers on the same
positive real ray are equal.

## 2. Primitivity identifies the entire Gaussian numerator

Put

```text
g_-=gcd_G(h_i,h_j),       g_+=gcd_G(h_i,bar(h_j)),
d_-=h_j bar(h_i)/N(g_-),  d_+=h_i h_j/N(g_+).                (3)
```

Both quotients are Gaussian integers and conjugate-primitive. This is
the exact primitive transition identity; for `d_+`, apply that identity
after conjugating one anchor. Their moduli are

```text
|d_-|=|h_i h_j|/N(g_-),
|d_+|=|h_i h_j|/N(g_+).                                    (4)
```

Two conjugate-primitive Gaussian integers on the same unoriented line
differ only by sign; on the same ray they are equal. Indeed their
quotient in the latter case is a positive rational number, and
primitive integer coordinates force its numerator and denominator to
be one. Therefore (2), together with the positive coordinates of `b`,
has the stronger exact arithmetic consequence

```text
b=d_+,
or b is d_- or bar(d_-), whichever has positive imaginary part.
                                                               (5)
```

In particular, the strict inequality

```text
|b| < |h_i h_j|/N(g_+)                                     (6)
```

excludes the sum branch and forces the primitive difference numerator.
In that case

```text
N(g_-)=sqrt(q_i q_j/q_b)                                    (7)
```

is an exact positive integer equality. Thus neither the integer square
root in (7) nor the Gaussian gcd equality needs to be postulated in
addition to (1), primitivity, and the height gap.

There is also an immediate arithmetic consequence before selecting
either branch. Clearing the denominators in (1) gives

```text
(q_j q_b t_i^2+q_i q_b t_j^2-q_i q_j s^2)^2
 =4 q_i q_j q_b (t_i t_j y)^2.                              (8)
```

All coordinates on the right are nonzero integers. Hence
`sqrt(q_i q_j q_b)` is rational, and therefore an integer. Once
(5) is known, the three quotients of this integer by the respective
pair norms are integers as well, by the Gaussian gcd valuation lemma
in the spin note. This derives the integral Heron coefficients from
the normalized equation and primitive inputs.

## 3. The oriented core makes the height gap automatic

Suppose the anchors have the actual form

```text
h_i=K_i A_i,
gcd_G(A_i,bar(A_j))=1       for every i,j,
kappa_i=log|K_i|.
```

The independent oriented core blocks give the displayed coprimality.
At each Gaussian prime, the minimum of the valuations in `h_i` and
`bar(h_j)` is at most the sum of the valuations in `K_i` and
`bar(K_j)`. Therefore the exact divisibility

```text
gcd_G(h_i,bar(h_j)) divides K_i bar(K_j)                     (9)
```

retains all prime powers, and implies

```text
log|d_+| >= log|h_i|+log|h_j|-2kappa_i-2kappa_j.             (10)
```

For a finite version, assume

```text
|log|h_i|-W/4|<=E,       |log|h_j|-W/4|<=E,
log|b|<=W/4+E,
W/4>3E+2kappa_i+2kappa_j.                                  (11)
```

Then (10) is strictly greater than `W/4+E`, so (6) follows.
Consequently an integral primitive solution of the Heron equation in
this height range is forced to be the exact primitive difference.

In central truncation, `E,kappa_i,kappa_j=O(kW/sqrt(M)+1)`.
Thus (11) holds simultaneously for all selected pairs when `M` is
large. The competing sum branch has logarithmic modulus
`W/2-o(W)`, whereas every prescribed pair modulus is
`W/4+o(W)`. This is a positive height gap, rather than a matched
trace or determinant estimate.

## 4. Simultaneous reconstruction and what remains independent

Choose the anchor at an endpoint, so that all anchor arguments are
positive. Let `h_1,...,h_m` be distinct primitive anchors as above.
For every pair choose a primitive candidate `b_ij` with positive
coordinates and the prescribed pair height. Suppose each anchor
triple satisfies (1), and suppose (11) holds for every pair.

Then (5)--(11) identify every candidate with the absolute-argument
primitive difference. Ordering the anchors by argument fixes its
orientation. In that ordering, write

```text
b_ij=x_ij+i t_ij,       t_ij>0       for i<j,
G_ij=N(gcd_G(h_i,h_j)).
```

One obtains the exact simultaneous identities

```text
h_j bar(h_i)=G_ij b_ij,
(b_ij/bar(b_ij))(b_jl/bar(b_jl))=b_il/bar(b_il),
G_ij t_ij=x_i t_j-x_j t_i.                                 (12)
```

Let `L=lcm_G(h_1,...,h_m)`. Then the Gaussian lattice points

```text
z_0=bar(L),       z_i=bar(L) h_i/bar(h_i)
```

give the exact integral realization, with least radius `|L|`.
All the prescribed `b_ij` are its actual primitive pair numerators.
Every remaining triple Heron equation, every four-index Pluecker
identity, and the common-factor-normalized versions of those identities
now follow from (12). Their coefficients are the actual Gaussian gcds;
they have not been replaced by freely chosen norm variables.

For example, the unsquared identity for an ordered triple is

```text
t_13^2 Z=t_12^2 U+t_23^2 V+2t_12 t_23 x_13,                 (13)
```

with the integral coefficients of the spin note. It is the selected
branch of the Heron equation, and follows by the elementary sine
addition identity. Four-index comparisons of (13) therefore yield
the already reconstructed Gaussian cocycle, rather than a new
independent condition.

This is a sufficiency theorem for the indicated system. It does not
show that the simultaneous integral Heron solutions exist with the
full large core profile, or that they cannot exist. Establishing a
positive residue exponent still requires an arithmetic obstruction to
those anchor solutions, or a new divisibility statement not already
implied by (12).

## 5. A four-point countermodel to unsigned local Heron tests

Take

```text
a=10+i,       b=12+i,       c=ab=119+22i,
N(a)=101,     N(b)=145,     N(c)=101*145.
```

All are conjugate-primitive, their norm supports are disjoint for
`a,b`, and their arguments satisfy
`arg(c)=arg(a)+arg(b)<pi/4`. On a four-cycle assign

```text
h_01=h_23=a,
h_12=h_03=b,
h_02=h_13=c.                                               (14)
```

Every triangle has the same three angular distances
`alpha,beta,alpha+beta`. Thus every unsigned Heron equation holds,
and its three integral coefficients are `101,145,1` in some order:
they are pairwise coprime. The prime valuation metrics are themselves
binary cuts. Primes in `101` see the cut `{0,3}|{1,2}`, and those
in `145` see `{0,1}|{2,3}`. Hence pairwise coprime norm coefficients
and even a coherent **unoriented** norm-cut profile do not repair the
missing angular branch information.

Nevertheless these six distances do not embed in a real interval.
Anchor zero has distances `alpha` to vertex one, `beta` to vertex
three, and `alpha+beta` to vertex two. The first two vertices must
lie between zero and vertex two, forcing their mutual distance to be
`|alpha-beta|`, not the prescribed `alpha+beta`.

There is an exact arithmetic witness. The actual difference numerator
between the anchor directions `a,b` is `a bar(b)=121+2i`; its
norm equals `N(c)`, but its residue is `2`, not `22`. The two
Gaussian branches have equal modulus here, because their relevant
gcds are both units. Thus the strict gap (6) correctly fails.

This countermodel is not an endpoint counterexample and does not have
the full uniform cut profile. It isolates a specific failure of local
unsigned Heron and norm-only tests. Sections 2--4 explain precisely
how the actual oriented core height range removes that failure.

## Verification

The exact cleared Heron equation, both primitive branches, their full
gcd norm formulas, and the displayed countermodel passed 1,000
Gaussian-integer checks with seed 82614. Two independent audits checked
the branch alternatives, denominator clearing, core-gcd bound, finite
gap, simultaneous realization, and the scope of the countermodel.
No new Lean formalization is claimed.
