# High-return pattern descent and the dense cancellation branch

## Status

This note proves a new structural reduction in the odd squarefree split
model.  A hypothetical unbounded endpoint cluster cannot put a positive
fraction of its conductor in a bounded number of recurrent prime patterns.
Every complete orientation-pattern block has logarithmic weight
`O(H/n)`, while the Farey/run inequality forces a positive fraction of the
conductor into patterns that return linearly many times.  Consequently the
hard branch contains linearly many pairwise-disjoint, one-use Gaussian
pattern blocks.

More precisely, if the cluster has `n` points and `H=log N`, then for every
fixed `0<a<1/4` the blocks whose words return at least `a(n-1)` times include
at least

```text
((1/4-a)/(1-2a)) n-o(n)                            (0.1)
```

distinct complete words.  At `a=1/8` this is `n/6-o(n)`.  Moreover, a
positive proportion of the consecutive transitions are both short and each
flip a linear number of these blocks.

The argument is rigorous at paper level but is not yet formalized in Lean.
It does not prove the uniform endpoint bound.  It reduces the remaining
case to a dense signed phase-cancellation problem that the abstract
character code can model at the support level.

---

## 1. Rooted squarefree model

Let

```text
z_0,z_1,...,z_(n-1)
```

be distinct ordered points in an endpoint arc.  Put

```text
m = n-1,
H = log N,
angular diameter <= C exp(-H/4).
```

After the standard removal of common, inert, and ramified factors, assume
that `N` is an odd squarefree product of split rational primes.  For every
conductor prime `p`, record its orientation along the ordered points by a
binary word

```text
e_p(0),e_p(1),...,e_p(m),
```

normalized by `e_p(0)=0`.  If `q_ij` is the norm of the forced-reduced
Gaussian transition between points `i` and `j`, then

```text
log q_ij = sum_(p:e_p(i) != e_p(j)) log p.          (1.1)
```

Its primitive imaginary coordinate is a nonzero integer.  The endpoint
diameter therefore gives the pair-separation estimate

```text
log q_ij >= H/2-kappa_C,                            (1.2)
```

where one may take `kappa_C=2 log(max(C,1))` after harmless normalization
of the angular parameter.  Only the fact that `kappa_C=O_C(1)` is used
below.

For a word `e_p`, let

```text
V_p = number of adjacent switches,
r_p = floor(V_p/2).
```

Thus `r_p` is the number of completed returns to the root orientation.

---

## 2. Forced conductor mass at high return count

Let `c_p` be the number of nonempty one-runs of `e_p` among the nonroot
indices `1,...,m`, and define

```text
Lambda = H^(-1) sum_p c_p log p.
```

The half-angle Farey product inequality gives

```text
(m-1) log(m-1)
  <= H(Lambda-(m+1)/4)+(m+1)log C.                 (2.1)
```

Because the word starts at zero,

```text
r_p >= c_p-1.                                      (2.2)
```

Summing (2.2) over the entire conductor and using (2.1) gives the exact
lower bound

```text
sum_p r_p log p
  >= (m-3)H/4+(m-1)log(m-1)-(m+1)log C.            (2.3)
```

Fix `0<a<1/4`, and put

```text
P_a = {p : r_p >= a m},
W_a = sum_(p in P_a) log p.
```

Since every binary word on `m+1` positions has `r_p<=m/2`,

```text
sum_p r_p log p
  <= a m(H-W_a)+(m/2)W_a.                          (2.4)
```

Combining (2.3) and (2.4) yields

```text
W_a >= [ (1/4-a)mH-3H/4
          +(m-1)log(m-1)-(m+1)log C ]
        / [ (1/2-a)m ].                            (2.5)
```

Along a hypothetical unbounded family, the already proved sublogarithmic
bound implies `m/H -> 0`.  Hence

```text
W_a/H >= (1/4-a)/(1/2-a)-o(1).                    (2.6)
```

For `a=1/8`, this becomes

```text
W_(1/8) >= H/3-o(H).                               (2.7)
```

Thus a fixed positive fraction of the conductor is carried by prime words
that oscillate linearly many times.

---

## 3. The pattern-block descent bound

Group conductor primes by their complete binary word on all `n` points.  A
pattern block is the product of all Gaussian primes in one word class.  Let
`w` be its logarithmic norm.

Choose the larger of the zero and one classes of this word, and call its
cardinality `M`.  Then

```text
M >= ceil(n/2).                                    (3.1)
```

Every point in that class contains the same oriented Gaussian block: the
block itself on one side and its conjugate on the other.  Dividing it from
all `M` points is an injective similarity to a circle whose remaining
conductor height is `H-w`.  Equivalently, all pair-difference supports inside
the class lie in the remaining conductor.

Sum (1.2) over the `binom(M,2)` pairs.  At each remaining conductor prime,
the corresponding binary cut separates at most `floor(M^2/4)` pairs.
Therefore

```text
binom(M,2)(H/2-kappa_C)
  <= floor(M^2/4)(H-w).                            (3.2)
```

The parity cases give particularly clean exact forms.  If `M` is even,

```text
M w <= H+2(M-1)kappa_C.                            (3.3)
```

If `M` is odd,

```text
(M+1)w <= H+2M kappa_C.                            (3.4)
```

In either case,

```text
w <= H/M+2kappa_C
  <= H/ceil(n/2)+2kappa_C.                         (3.5)
```

Since the sublogarithmic theorem also gives `H/n -> infinity` along an
unbounded cluster sequence,

```text
w <= (2+o(1))H/n.                                  (3.6)
```

This estimate is important for two reasons.

1. It applies separately to every complete pattern block, not merely to the
   heaviest one.
2. It uses the original height `H` in the pair lower bound (1.2), even after
   descending to height `H-w`.  The angular diameter is unchanged by the
   common-factor similarity.

Globally constant words should simply be divided from the whole cluster and
removed before applying the lemma.

---

## 4. Linearly many high-return patterns

Let `K_a` be the number of distinct complete words represented among the
primes in `P_a`.  The blocks are pairwise disjoint in rational-prime support,
their total weight is `W_a`, and each block obeys (3.6).  Consequently

```text
K_a >= W_a / ((2+o(1))H/n)
    >= ((1/4-a)/(1-2a))n-o(n).                     (4.1)
```

At `a=1/8`,

```text
K_(1/8) >= n/6-o(n).                               (4.2)
```

Letting a fixed `a` tend to zero makes the coefficient in (4.1) tend to
`1/4`, while every selected word still has `Omega_a(n)` returns.  Hence, for
any fixed positive margin, the hard branch contains almost `n/4` distinct
linearly oscillatory pattern blocks.

Each block has norm at most

```text
N^(2/n+o(1)),                                      (4.3)
```

so these are genuinely diffuse `N^o(1)` atoms, not a bounded family of heavy
compensators.

---

## 5. Dense short transitions

Use the threshold `a=1/8`, and abbreviate its total block weight by `W`.
For the consecutive edge `k`, let

```text
S_k = total weight of high-return pattern blocks
      whose word switches across k.                (5.1)
```

Every selected word has

```text
V_p >= 2r_p >= m/4.                                (5.2)
```

Therefore

```text
sum_(k=0)^(m-1) S_k
  = sum_blocks V_block w_block
  >= mW/4.                                         (5.3)
```

Since `0<=S_k<=W`, at least `m/7` edges satisfy

```text
S_k >= W/8.                                        (5.4)
```

Indeed, if `t` edges satisfy (5.4), then the left side of (5.3) is at most

```text
tW+(m-t)W/8 = mW/8+7tW/8.
```

Using (2.7) and the individual block bound (3.6), every edge in (5.4) flips
at least

```text
n/48-o(n)                                          (5.5)
```

distinct high-return blocks.

At least `7m/8` of the `m` positive gaps are at most eight times their mean.
Intersecting this set with the `m/7` dense edges leaves at least

```text
m/56-O(1)                                          (5.6)
```

edges that are simultaneously dense and short.  For each such edge, the
primitive-transition integrality estimate gives

```text
q_k >= m^2 sqrt(N)/(64C^2),                        (5.7)
```

while its high-return part is a signed product of at least `n/48-o(n)`
pairwise-disjoint one-use blocks.

Thus the hard branch is not a sparse repair process.  It contains linearly
many short positive transitions, each built from linearly many small
recurrent Gaussian atoms.

---

## 6. A finite-body phase obstruction

Distinct pattern blocks have disjoint rational-prime support.  If `G_j` and
`G_k` are two such nonunit blocks, their primitive Gaussian rays cannot be
equal.  Hence

```text
1 <= |Im(conjugate(G_j)G_k)|
   <= sqrt(Norm(G_j)Norm(G_k)) |sin(beta_j-beta_k)|.
```

Using (3.6),

```text
|sin(beta_j-beta_k)|
  >= exp(-2H/n-o(H/n)).                             (6.1)
```

This is enormously larger than the endpoint width `exp(-H/4)` when `n`
grows.  More generally, if two disjoint products of high-return atoms have
directions within the endpoint width, their total logarithmic norm is at
least `H/2-O_C(1)`.  Since every atom has weight at most
`(2+o(1))H/n`, the two products together must use at least

```text
n/4-o(n)                                           (6.2)
```

atoms.

Consequently no bounded-body extraction can explain the small transition
angles.  Cancellation must be genuinely collective.  This also explains why
the four-coprime-block sector lemma does not immediately apply: the individual
pattern-block rays are separated, while only large signed aggregate products
are endpoint-close.

There is nevertheless a sharpened short-gap sector consequence.  Let
`G_1,...,G_4` be four nonunit effective products made from pairwise-disjoint
sets of pattern blocks, with the one-use budget

```text
sum_(j=1)^4 log Norm(G_j) <= H.
```

Suppose their rays all lie in a sector of width

```text
s <= A_C exp(-H/4)/m.
```

For each of the six pairs, primitive determinant spacing gives

```text
log Norm(G_i)+log Norm(G_j)
  >= 2log(1/s)
  >= H/2+2log m-O_C(1).                             (6.3)
```

Summing (6.3) over all pairs counts every block height three times, so

```text
sum_(j=1)^4 log Norm(G_j)
  >= H+4log m-O_C(1),                               (6.4)
```

contradicting the one-use budget for large `m`.  Thus a short-gap sector
contains at most three pairwise-disjoint nonunit effective products.  At the
ordinary endpoint width without the factor `1/m`, the extra `4log m` is absent
and four bodies remain possible.  The short-gap surplus is therefore real,
although dense transitions evade it by overlapping rather than by producing
four disjoint bodies.

### A coherent-core iteration across short transitions

The dense short-edge set also supports a growing common Gaussian core.  Let
`E_0` be the set from (5.6), and let `W` be the total high-return weight.  Since
every edge in `E_0` flips high weight at least `W/8`, weighted averaging finds
one high pattern block that occurs on at least `|E_0|/8` of those edges.
Fixing one of its two switch orientations retains a set `E_1` with

```text
|E_1| >= |E_0|/16 >= m/896-O(1).                   (6.5)
```

Suppose inductively that a product `K_d` of `d` distinct complete pattern
blocks, in fixed Gaussian orientations, divides every primitive transition
block `B_i` for `i` in `E_d`.  Put

```text
c_d = log Norm(K_d),
B_i = K_d D_i.
```

The pattern-block cap (3.6) gives

```text
c_d <= (2+o(1))dH/n.                               (6.6)
```

Every original edge in this family is short, so

```text
log Norm(D_i)
  >= H/2+2log m-O_C(1)-c_d.                        (6.7)
```

The tail supports lie in the remaining one-use conductor of height `H-c_d`.
If `lambda_i` is any fractional matching on these supports, then

```text
[sum_i lambda_i]
  [H/2+2log m-O_C(1)-c_d]
    <= H-c_d.                                      (6.8)
```

For `d=o(n)`, the fractional matching number is therefore

```text
at most 2+O(d/n)+o(1).                              (6.9)
```

Linear-programming duality supplies a remaining prime occurring on at least
`|E_d|/(2+O(d/n)+o(1))` tails.  All primes with its complete word have the
same occurrence pattern; grouping them into their full pattern block and
fixing one switch orientation loses at most another factor two.  This extends
the coherent core and gives

```text
|E_(d+1)| >= |E_d|/(4+O(d/n)+o(1)).                (6.10)
```

Consequently, uniformly for `d=O(log m)`,

```text
|E_d| >= (m/896)(4+o(1))^(-(d-1)).                 (6.11)
```

Taking `d=o(log m)`, for example `d=floor(log log m)`, produces an unbounded
coherent core of distinct Gaussian pattern blocks that divides
`m^(1-o(1))` short transitions.  Taking `d=c log m` retains a polynomial-size
family `m^(1-c log 4+o(1))` whenever `c<1/log 4`.

This transition statement also gives a genuine point descent.  After the
switch orientation of every core block has been fixed, all left endpoints of
the edges in `E_d` have one common `d`-bit signature, and all right endpoints
have the complementary signature.  Dividing `K_d` from either endpoint class
therefore gives `|E_d|` distinct points on conductor height `H-c_d`, with
endpoint constant

```text
C_d = C exp(-c_d/4).                                (6.12)
```

The descent is quantitatively real but does not yet close the argument.  In
the diffuse character model, `c_d` is only of order `dH/n` while the retained
cardinality loses a factor about `4^d`; the gain in (6.12) is exactly too small
for the known endpoint inequalities to beat that loss.

After division, all tails `D_i` lie in the same translated sector of width
`O_C(exp(-H/4)/m)`.  The four-body lemma shows that this tail family contains
no four members with pairwise-disjoint one-use supports.  The character model
realizes the factor-four loss at every step, so support theory alone cannot
turn the coherent core into a larger descent.  Its value is to isolate a
large family of actual Gaussian coordinate products on which a future
residual or resultant theorem can act.

---

## 7. Sharp support-level obstruction

Take `n=2^r` rows indexed by `F_2^r`, and one equal-weight conductor block for
every nonzero linear character.  Pairwise weighted Hamming distance is

```text
nH/(2(n-1)) > H/2.
```

There are orderings of the rows in which every balanced character switches
linearly often.  For example, a random permutation has this property
simultaneously for all `n-1` characters with positive probability by standard
bounded-difference concentration and a union bound.  In this model:

- every conductor block is high-return;
- there are `K=n-1` distinct complete words;
- every bounded set of revealed independent characters leaves large fibers;
- the joint signature map is injective.

Thus all conclusions above are compatible with the diagonal character model.
Return counts, pattern entropy, support matching, and bounded-depth descent
cannot finish the proof on their own.  The Gaussian Walsh theorem rules out a
large character-compatible affine cube only after the required phase and
tail-compatibility hypotheses have been extracted.

There is an even sharper weighted obstruction when `H` is sufficiently large.
Use every nonzero rooted Boolean word as an abstract atom, and mix a small
amount of the nontrivial-character distribution with the uniform distribution
on these words.  For a nonempty subset `J` of transition edges, its boundary
vector on `F_2^r` is nonzero and has total sum zero.  Fourier inversion gives

```text
sum_(nontrivial characters chi) |hat boundary(chi)| >= n.
```

Consequently the character part improves the normalized positive-subset width
from `1/2` to `1/2+1/(2m)`, where `m=n-1`.  Giving it weight

```text
lambda = 4m log(4m)/H
```

supplies

```text
W(J) >= H/2+2log(4m)                                (7.1)
```

simultaneously for every nonempty `J`.  The uniform Boolean part preserves
the dense high-return behavior outside exponentially small weight.  Every
atom switches somewhere, so the formal one-use lcm height of all transition
products is still exactly `H`; distinct transition supports have enormous
overlap.  Because all rooted words occur, the nonroot row vectors are linearly
independent over `F_2`, so this model contains no affine parallelogram, let
alone an affine five-cube.

The transition matrix has full row rank because the word set contains the
standard-basis columns.  Hence, at the level of real torus angles, it can
realize arbitrary distinct positive gaps of the required total width after a
generic kernel perturbation, while avoiding exact multiplicative collisions.
What it does not provide is a realization by pairwise-coprime,
coordinate-primitive Gaussian blocks with total lcm height `H`.  This proves
that dense returns, all positive-subset width inequalities, one-use lcm
bookkeeping, and exact real phase solvability still do not suffice.  The next
input must use simultaneous Gaussian integrality or resultants, not another
support or scalar-height inequality.

The exact six-point Gaussian example provides a finite calibration.  In its
threshold-layer coordinates, the high-return layers carry all but one factor
`5` of the conductor and form seven distinct nonconstant words whose joint
signatures separate all six points.  Heavy recurrent mass therefore need not
produce even a two-point complete-signature fiber at finite scale.

---

## 8. Corrected next target

Any hypothetical unbounded endpoint family now has a sharply defined hard
subsystem:

```text
linearly many disjoint Gaussian atoms G_j,
log Norm(G_j) <= (2+o(1))H/n,
each atom switches linearly many times,
linearly many short edges use linearly many atoms,
every full signed edge product has positive endpoint-small phase.
```

A successful continuation must use the actual Gaussian phases to exclude
this dense cancellation matrix.  Sufficient forms would include:

1. extracting a character-compatible affine five-cube, so Walsh rigidity
   applies;
2. extracting disjoint **aggregate** products in one endpoint sector while
   preserving the one-use conductor budget; or
3. proving a height-efficiency theorem showing that this many simultaneous
   signed near-relations force extra Gaussian lcm height beyond `H`.

The last formulation is consistent with the oriented CRT stress test: every
finite character pattern can be realized by genuine Gaussian divisors with
bounded residuals and arbitrarily small absolute angles, but the lcm height of
that realization is too large for the critical `exp(-H/4)` scale.
