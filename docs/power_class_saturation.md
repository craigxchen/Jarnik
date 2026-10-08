# Independent correction cuts saturate every norm power-class test

## Status

This note constructs correction and core data with three simultaneous
properties. The correction height is arbitrarily small relative to each
core block. Products of primitive pair norms have exactly the perfect-power
relations that are universal for lattice-circle configurations, for every
power order at once. Finally, the desired equations `Im h_i=1` can be made
valid modulo any prescribed finite integer modulus.

The construction does **not** prove `Im h_i=1`, or any sufficiently small
bound on its actual integer value. It is therefore not an endpoint cluster
construction. Rather, it isolates an obstruction to excluding the systems
in [independent_block_reconstruction.md](independent_block_reconstruction.md)
using only power classes and a prescribed finite set of local congruences.
No positive correction-height exponent or uniform endpoint bound is proved.

## 1. The universal power relations among pair norms

Let the vertex set be `V={0,...,m}`, and let `E` be the edges of its
complete graph. For a cut `S subset V`, write

```text
delta_S(ij)=1  if exactly one of i,j belongs to S,
delta_S(ij)=0  otherwise.
```

Complementary cuts are identical. Thus it suffices to use the nonempty
subsets of `[m]`, with vertex zero always outside the displayed cut.

For a fixed integer `q>=2` and an edge coefficient vector
`c=(c_ij) in (Z/qZ)^E`, the following description is exact:

```text
sum_(ij crossing S) c_ij=0 mod q for every cut S

iff

q odd:   c_ij=0 mod q for every edge;

q even:  c=(q/2)b mod q,
         where b is a binary edge set with even degree at every vertex.
                                                               (1)
```

To prove this, the integer vector identity

```text
delta_{ {i} }+delta_{ {j} }-delta_{ {i,j} }=2 e_ij
```

shows that vanishing of all cut sums implies `2c_ij=0 mod q`.
Singleton cuts also give zero sums of the coefficients incident to each
vertex. Conversely, for any cut `S`,

```text
sum_(i in S) sum_(j!=i) c_ij
  =sum_(ij crossing S)c_ij+2sum_(i<j, i,j in S)c_ij.
```

Thus those two conditions imply vanishing of every cut sum. If `q` is
odd, multiplication by two is invertible. If `q` is even, write
`c_ij=(q/2)b_ij`; the vertex conditions say exactly that the binary
degrees of `b` are even. This proves (1), including composite `q`.

These are universal perfect-power relations for primitive pair norms
in the common-unit split-prime allocation setting used here. Indeed,
at a split prime with allocation exponents `a_i`,
the exponent in the pair norm is `|a_i-a_j|`. Its threshold expansion is

```text
|a_i-a_j|=sum_t delta_{ {v:a_v>=t} }(ij).        (2)
```

Therefore a coefficient vector satisfying all cut sums modulo `q`
always makes the corresponding product of pair norms a rational
`q`th power. Negative coefficients are allowed, so this statement
includes quotients of pair norms. For example, products around binary
cycles are always rational squares. Formula (1) describes the analogous
universal consequences for every other power order.

## 2. Correction data realizing exactly the universal relations

Choose, for every nonempty cut `S subset [m]`, a Gaussian integer
`J_S` with the following properties:

* every `J_S` is coprime to its conjugate;
* the rational norms `N(J_S)` are pairwise coprime;
* each `N(J_S)` has a private rational prime divisor `p_S` whose
  exponent is exactly one.

Choose independent conjugate-coprime core blocks `H_S`, including an
empty block, whose prime supports avoid every `N(J_S)`. Put

```text
A_i=product_(S contains i) H_S,
K_i=product_(S contains i) J_S,
h_i=K_i A_i,
h_0=1.
```

All `h_i` are Gaussian integers coprime to their conjugates. Their
primitive change-of-anchor numerators are exactly

```text
h_ij=
 product_(j in S, i notin S) (J_S H_S)
 product_(i in S, j notin S) conj(J_S H_S).      (3)
```

This follows either directly by cancellation, or from the gcd identity
in [all_anchor_primitive_transition.md](all_anchor_primitive_transition.md).
Let `n_ij=N(h_ij)`. At the private correction prime `p_S`,

```text
v_(p_S)(n_ij)=delta_S(ij).                      (4)
```

Every other prime exponent is a positive integer multiple of one of
these same cut vectors. Consequently, simultaneously for **every**
integer `q>=2`,

```text
product_(ij in E) n_ij^(c_ij) is a rational qth power

iff  every cut sum sum_(ij crossing S)c_ij vanishes mod q.       (5)
```

Necessity follows from the private exponent-one primes. Sufficiency
follows because all remaining prime rows are multiples of cut rows.
Thus (1) lists every power relation in this construction. No
non-universal perfect-power relation remains for a separation argument
to exploit, even if it tests all pair norms and all orders at once.

In particular, for any anchor, its `m` primitive norm squareclasses
are linearly independent. A graph supported in an anchor star has no
nonempty binary even-degree edge set, so (1) proves independence;
the same argument gives independence modulo every `q`. Hence the
six-wise squareclass restriction is satisfied with full independence,
uniformly over all anchors.

## 3. Signed phase power classes are saturated as well

The singleton correction block `J_{ {i} }` has a private Gaussian
prime appearing in `h_i` with exponent one, and in no other `h_j`.
For an integer vector `(c_1,...,c_m)`, consider the primitive Gaussian
numerator of

```text
product_i (h_i/bar(h_i))^c_i.
```

If all its Gaussian prime exponents are divisible by `q`, then each
`c_i` is divisible by `q`, by examination of these private primes.
Conversely, that divisibility of the coefficients makes every prime
exponent divisible by `q`. The same conclusion holds after any change
of anchor, since changing from vertex-zero differences to another
anchor is an invertible integer change of basis.

There is also a useful quantitative check for phase tests using the
uniform cut load. Define

```text
mu(c)=average_(S subset [m]) |sum_(i in S)c_i|.
```

Conditioning on all bits except the `i`th and applying
`|y|+|y+c_i|>=|c_i|` gives

```text
mu(c)>=max_i |c_i|/2.                           (6)
```

Thus a nonzero coefficient vector divisible by `q` has load at least
`q/2`. This construction eliminates nontrivial low-load power-class
relations without requiring any increase in core block height.

## 4. Exact integral corrections of one common modulus

Set

```text
D=product_(S nonempty) N(J_S),
B_0=D,
B_i=D K_i^2/N(K_i).
```

Since `N(K_i)` divides `D`, every `B_i` is a Gaussian integer, and

```text
|B_i|=D,
B_i/B_0=K_i/bar(K_i).                           (7)
```

Hence these are actual equal-norm integral correction factors of the
type required in the centrally truncated setting. If a prime occurs
in `J_S` to exponent `f`, the exponents in `B_i` are `2f,0` in the
two Gaussian orientations for inside rows, and `f,f` for outside
rows. The anchor is an outside row. At this prime the central
truncation with `r=f` removes the full exponent `2f` from the core.
There is no fractional correction or uncharged denominator.

This describes the selected correction factors themselves; it does
not assert that the resulting full points have small angular span.
In particular, no small-angle conclusion follows from (7).

## 5. Every prescribed finite modulus can also be accommodated

Fix an even positive integer `Q`, and any desired positive integer
`L`. The preceding construction can additionally be arranged so that

```text
J_S = i mod Q        when S is a singleton,
J_S = 1 mod Q        otherwise,
H_S = 1 mod Q        for every S,
H_S is a Gaussian Lth power for every S.        (8)
```

It follows that

```text
h_i=i mod Q,
Re h_i=0 mod Q,
Im h_i=1 mod Q,
N(h_i)=1 mod Q.                                (9)
```

Thus the proposed small-value equations `Im h_i=1` pass all
congruence tests modulo this prescribed `Q`. The actual integer
imaginary parts have not been bounded.

Here is an elementary construction of the correction blocks, which
does not require primes in a prescribed Gaussian residue class.
Choose them successively. For a singleton use `J=a+i`, with
`a=0 mod Q`. For any other cut use `J=a+iQ`, with `a=1 mod Q`.
In both cases write the fixed imaginary part as `b`, and impose

```text
a=0 mod product_(earlier cuts T) N(J_T).
```

Each earlier norm is coprime to `Q`, hence to `b`. Therefore the
new norm is coprime to every earlier norm. Choose a fresh split
prime `p=1 mod 4` avoiding all these moduli, and prescribe `a`
modulo `p^2` so that

```text
a^2+b^2=0 mod p,       a^2+b^2!=0 mod p^2.
```

Such a residue exists: either nonzero square root of `-b^2` modulo
`p` has exactly one lift that is a root modulo `p^2`, and any other
lift works. The Chinese remainder theorem supplies `a` satisfying
all conditions. The coordinates of `J` are coprime, and have
opposite parity, so `J` is coprime to its conjugate. Its norm is
coprime to `Q` and has the required private exponent-one prime.
All rational prime divisors of its norm are odd split primes.

Now choose the core blocks on fresh split Gaussian prime supports.
Their classes modulo `Q` are units in the finite ring `Z[i]/(Q)`.
Take each exponent to be a multiple of `L` and of the order of
the corresponding residue class in that finite unit group. This
gives both core conditions in (8).

For different nonreference rows `i,j`, their common factors contain
no singleton correction blocks. Hence their common Gaussian divisor
is `1 mod Q`, and (3) also gives `h_ij=1 mod Q`. These pairwise
congruences are consistent with the exact change-of-anchor equations;
they do not assert a zero actual pair residue.

### Finite depths at every chosen core prime can also be imposed

The congruence construction need not avoid the eventual core support.
Choose first the distinct oriented core primes `gamma_S`, including
the empty block's prime, and arbitrary positive finite depths `a_S`.
Write `p_S=N(gamma_S)`; these labels in this paragraph refer to
core primes, not the private correction primes of Section 2.

Choose a common exponent period divisible by `L`, by the orders of
all the core primes modulo two, and by their orders in all the
finite unit groups

```text
(Z[i]/gamma_S^(a_S))^*,
(Z[i]/bar(gamma_S)^(a_S))^*
```

where the prime in question is a unit. Take every core exponent to
be a sufficiently large multiple of this period. Thus `H_T=1`
in both quotients belonging to `S` whenever `T!=S`, while

```text
H_S=0 mod gamma_S^(a_S),
H_S=1 mod bar(gamma_S)^(a_S).
```

For row `i`, prescribe the following unit residues for `K_i`:

```text
i in S:     K_i=1   mod gamma_S^(a_S),
            K_i=2i  mod bar(gamma_S)^(a_S);

i notin S:  K_i=i   in both quotients.
```

Also prescribe `K_i=i mod 2`. The Gaussian Chinese remainder
theorem combines these into one unit residue class modulo

```text
Q_core=2 product_S p_S^(a_S).
```

When `i in S`, the corresponding `h_i=K_i A_i` has residues zero
and `2i` in the two conjugate quotients. Therefore
`(h_i-bar(h_i))/(2i)=1` in both quotients. When `i notin S`,
`h_i=i` in both quotients, giving the same result. Hence

```text
Im h_i=1 mod p_S^(a_S)  for every row i and every core prime S.   (10)
```

These prescribed `K_i` residues can be implemented by giving the
singleton correction block `J_{ {i} }` the desired residue class
and every nonsingleton correction block residue one. The following
elementary strengthening of the previous construction ensures that
the correction blocks still have all properties of Section 2.

Given any unit Gaussian residue `u mod Q_core` and the previously
chosen correction norms, impose `J=u mod Q_core` and `J=1`
modulo every previous norm. These conditions give a unit residue
`a_0+i b_0` modulo their product `M`. At a fresh split prime `p`
choose a residue modulo `p^2` with imaginary part one and norm
valuation exactly one, as above. Combine it with the residue
modulo `M`, obtaining prescribed coordinates modulo `Mp^2`.
Choose a nonzero integer `b` in the prescribed imaginary class.
For each prime dividing `b` but not `Mp^2`, further impose
`a=1` modulo that prime. The Chinese remainder theorem now
chooses `a` in its prescribed class modulo `Mp^2` and in these
additional classes.

The resulting coordinates are coprime. A prime dividing both `b`
and `Mp^2` cannot divide the prescribed real coordinate: the
residue is a unit modulo `M`, and the imaginary coordinate is one
modulo the fresh prime. All other prime factors of `b` were
explicitly avoided. The coordinates have opposite parity since
the prescribed residue modulo two is a unit. Thus `J` is
conjugate-coprime. Its norm is coprime to all old norms and to
the core support, and its fresh private prime has exponent one.

The correction blocks and their heights can therefore be fixed
after the core supports and finite depths are selected, but before
the core exponents grow. All norm power-class conclusions remain
valid. Congruences (10) can hold to arbitrarily prescribed finite
depths at **every** core prime with correction heights negligible
relative to each sufficiently large core block.

This does not handle depths growing proportionally to the full
core exponents while retaining the same height budget. In
particular it does not replace the actual integer equation
`Im h_i=1`, or the required subpower upper bound on that integer.

## 6. Uniform core weights and negligible correction height

For each fixed `m,Q,L`, the correction blocks can be fixed before
choosing the core exponents. Write

```text
E=sum_i log|K_i|
 =(1/2)sum_(S nonempty) |S| log N(J_S).
```

This is a finite constant at that stage. Choose a target block log
norm `w` tending to infinity. For every `S`, select an allowed core
exponent whose norm logarithm is closest to `w`. The allowed exponent
spacing is a fixed positive integer for each core prime, so

```text
log N(H_S)=w+O_(m,Q,L)(1).
```

Consequently `E=o(w)`, all core weights are asymptotically equal,
and the row marginals satisfy

```text
log|A_i|=W_core/4+o(w),
W_core=sum_S log N(H_S).
```

The same statements can be arranged along `m` tending to infinity,
even if `Q` and `L` grow arbitrarily with `m`. At each stage choose
`w` sufficiently large that both `E/w` and the total accumulated
core-exponent rounding error divided by `w` tend to zero. There is
no upper bound on the block scale in this construction. Equal-norm
correction height `log D`, total coefficient height `E`, and all
fixed congruence parameters can therefore be negligible relative
to every block simultaneously.

The reconstruction theorem would require the **actual** maximum
`T=max_i |Im h_i|` to satisfy `E+log T=o(w)`. That requirement has
not been established here. Equations (8)--(9), maximal power-class
separation, and the equal-norm integral corrections all coexist
without supplying it. Any exclusion of these data must therefore
use more than universal perfect-power relations and congruences at
prescribed finite depths, even across the chosen core support;
in particular it must retain
information that controls the actual small integer values.

As a finite supplement to the proof of (1), exhaustive integer checks
verified 67,630 coefficient vectors on complete graphs with two through
four vertices and power orders two through six. The calculation agreed
with the displayed kernel in every case. No Lean formalization is claimed.
