# Majority freezing: the 63-cut audit and an optimal 49-cut anchor design

This note audits the proposed eight-vertex freezing argument. The
63-cut Fano design has the claimed strict majority on every pair.
For determining the primitive configuration, however, the seven
anchor pairs suffice. An explicit 49-cut design gives strict
majorities on those pairs, and 49 is optimal for that combinatorial
criterion. The conclusion is uniqueness of the labelled primitive
configuration, not uniqueness of every block factorization and not
a uniform bound on the number of configurations.

## 1. Exact Gaussian uniqueness with a frozen divisor

Let `Q` be a nonzero Gaussian integer coprime to its conjugate.
Suppose

```text
Im(Q U)=t,       Im(Q U')=t',
U,U' in Z[i],       t,t' in Z minus {0}.
```

Then `Q(t'U-tU')` is real. Conjugating that equality shows that
`bar(Q)` divides `Q(t'U-tU')`; conjugate coprimality gives

```text
bar(Q) divides t'U-tU'.                         (1)
```

If

```text
|t'||U|+|t||U'|<|Q|,
```

the Gaussian integer in (1) vanishes. Consequently `U/t=U'/t'`.
In particular, among solutions satisfying `|U|<=U_max` and
`0<|t|<=T_max`, the condition

```text
2 U_max T_max<|Q|                               (2)
```

determines at most one ratio `U/t`, and therefore at most one
phase ratio `(Q U)/bar(Q U)`.

The nonzero-target assumption is necessary. No primitivity of `U`
is required. It is the frozen factor `Q` that must be coprime to
its conjugate, and independent oriented core blocks guarantee this.

## 2. The eight-vertex block model and the quantitative margin

Use vertices `0,1,...,7`, with each cut represented by a subset
`S` of `[7]` excluding the anchor zero. There are 128 such blocks.
Every fixed pair is separated by exactly 64 cuts. Suppose the
independent Gaussian blocks satisfy

```text
(1-eta)w <= log N(H_S) <= (1+eta)w.              (3)
```

For an ordered pair `i,j`, its core product uses `H_S` when
`j in S,i notin S`, and `bar(H_S)` when `i in S,j notin S`.
Suppose its actual primitive numerator is

```text
h_ij=K_ij product_(cuts separating i,j) oriented(H_S),
t_ij=Im h_ij !=0,
log|K_ij|<=kappa w,       log|t_ij|<=tau w.      (4)
```

For anchor pairs this is exactly the already proved central
factorization `h_i=K_i A_i`. For all pairs, (4) must come from the
audited all-anchor factorization or its trimmed version; it cannot
be inferred just by assigning arbitrary pair factors.

Freeze a labelled set of Gaussian blocks, including their actual
complex values and the fixed orientation relative to vertex zero.
If `a_ij` of its cuts separate `i,j`, put their oriented product
in `Q_ij` and the remaining factors, including `K_ij`, in `U_ij`.
Then

```text
log|Q_ij| >= a_ij(1-eta)w/2,
log|U_ij| <= (64-a_ij)(1+eta)w/2+kappa w.
```

Condition (2), uniformly for all candidate systems with the same
frozen data and bounds, follows from

```text
(a_ij-32-32eta-kappa-tau)w>log 2.               (5)
```

Thus every strict integer majority `a_ij>=33` suffices when
`eta,kappa,tau ->0` and `w->infinity`. Conjugating a factor changes
neither its modulus nor the required coprimality. Freezing just its
norm or prime support would not supply the fixed Gaussian `Q_ij`
needed by (1).

## 3. The proposed 63-cut Fano design is correct

Freeze all 35 subsets of size four, all 21 subsets of size five,
and the seven lines of a Fano plane on `[7]`. For example its lines
can be the cyclic translates of `{0,1,3}` modulo seven, with one
added to each vertex label. Each vertex lies on three lines and
each vertex pair lies on one line.

For an anchor pair `0,j`, the separating frozen cuts number

```text
binom(6,3)+binom(6,4)+3=20+15+3=38.
```

For two nonanchor vertices, the count is

```text
2 binom(5,3)+2 binom(5,4)+(3+3-2)=20+10+4=34.
```

These are both strict majorities of 64. The nominal logarithmic
`Q/U` margins are `6w` for anchor pairs and `2w` for the other
pairs. The uniform finite sufficient condition is

```text
(2-32eta-kappa-tau)w>log 2.                    (6)
```

There are 63 frozen blocks, fewer than half of all 128 blocks.
No empty-block factor occurs in any pair product.

## 4. Only anchor ratios are needed: a 49-cut construction

Fixing the phase ratios from zero to each other vertex fixes every
pair phase ratio. Thus it is enough to have a strict majority on
the seven anchor pairs. The following design freezes only 49 cuts.

Let `F_1` consist of the translates of `{0,1,3}` and `F_2` the
translates of `{0,2,3}` modulo seven, again relabelled by `[7]`.
These are disjoint Fano planes. In explicit labels they are

```text
F_1: 124, 137, 156, 235, 267, 346, 457,
F_2: 126, 134, 157, 237, 245, 356, 467.
```

Their 14 complementary four-sets are distinct. Each vertex lies
in four complements per Fano plane, hence in eight of these 14
four-sets. Add the further omitted four-set `{1,2,3,4}`, which
is not among those complements.

Now freeze all 29 subsets of sizes five, six, and seven, and all
the 20 four-sets that were not omitted. The total is `29+20=49`.
The large subsets contribute `15+6+1=22` anchor incidences per
vertex. The retained four-sets contribute eleven incidences for
vertices one through four and twelve for vertices five through
seven. Thus the seven separating counts are exactly

```text
33,33,33,33,34,34,34.                           (7)
```

All anchor ratios are therefore fixed under

```text
(1-32eta-kappa-tau)w>log 2.                    (8)
```

This design needs only the anchored equations, so no separate
all-pair correction factorization is required for its conclusion.

### Optimality for strict whole-cut anchor majorities

Every frozen cut `S` contributes exactly `|S|` to the sum of the
seven anchor-incidence counts. Among 48 cuts, this total is at
most the sum of the 48 largest subset sizes:

```text
7+7*6+21*5+19*4=230<7*33=231.
```

Thus some anchor edge has at most 32 separating frozen cuts.
No family of at most 48 whole cuts meets the strict-majority
criterion on all seven anchor edges. This proves optimality only
within this whole-cut counting criterion at asymptotically equal
block weights. It is not a lower bound on the amount of information
needed by every arithmetic method, by partial-factor freezing, or
by some different injectivity argument.

## 5. Exactly what is determined

The uniqueness lemma fixes `U_i/t_i`, hence `h_i/t_i` and
`h_i/bar(h_i)`, for every anchor pair. The targets may change sign
or scale, but that does not change these phase ratios.

Given all the phase ratios, the labelled **primitive** Gaussian
configuration is unique up to multiplication by one common
Gaussian unit. Indeed any two such configurations differ by a
common factor `c in Q(i)`. Primitivity of the first configuration
and Gaussian Bezout imply `c in Z[i]`; primitivity of the second
implies `c^(-1) in Z[i]`. Thus `c` is a unit. Equivalently, the
Gaussian-lcm construction gives the unique least radius for the
fixed ratios, as proved in
[independent_block_reconstruction.md](independent_block_reconstruction.md).

The complete block-data tuple need not be unique. A small divisor
of an unfrozen block can be transferred into every incident `K_i`
without changing any `h_i`. Such transfers can preserve all the
asymptotic height and block-weight bounds. Consequently a claim
that the frozen data determine every unfrozen factor requires an
additional canonical factorization convention. The argument here
determines the primitive point configuration itself.

Exact finite enumeration verified both Fano systems, all 28
separating counts in the 63-cut design, all seven counts in the
49-cut design, and the 48-cut incidence upper bound. This is a
combinatorial and arithmetic injectivity result; it does not by
itself bound how many different frozen data sets can occur as
the radius grows.
