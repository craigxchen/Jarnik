# Midpoint reflection has a quantitative endpoint scaling cost

The reflection integrality domain was already computed in
[the rational reflection route](rational_reflection_route.md).
This note adds a different, quantitative conclusion: in the growing
balanced independent-block regime, making an entire cluster and its
pair-midpoint reflection simultaneously integral forces the normalized
arc constant to grow at least as `R^(1/4-o(1))`. This holds uniformly
over the choice of pair. The very small original endpoint constant
does not pay for the required scaling.

This is an obstruction to that particular symmetrization operation,
not an exclusion of the original integer configurations. It does not
prove the uniform endpoint bound or a new balanced-cut bonus.

## 1. The exact cost of adjoining a reflected primitive cluster

Let `Z={z_0,...,z_m}` be Gaussian integers of common nonzero modulus
`R`, with common Gaussian gcd one. Let `rho in Q(i)`, `|rho|=1`,
and write its reduced Gaussian fraction as

```text
rho=a/d,       a,d in Z[i],       gcd_G(a,d)=1.
```

Necessarily `|a|=|d|`. The reflection is `z -> rho bar z`.
If a common complex multiplier `gamma` makes both configurations

```text
gamma Z,       gamma rho bar Z                       (1)
```

Gaussian integral, then

```text
gamma in d Z[i],       |gamma|>=|d|.                 (2)
```

Indeed Bezout coefficients for the primitive original cluster show
that `gamma` is a Gaussian integer. Bezout for its conjugate likewise
shows that `gamma rho` is a Gaussian integer. These two statements
are equivalent to the divisibility in (2). Conversely `gamma=d`
works. Its union `d Z union a bar Z` is primitive, since any common
Gaussian divisor would divide both `d` and `a`. Thus common-gcd
normalization of the enlarged configuration cannot remove this cost.

Write `D_G(rho)=|d|`. For norm-one Gaussian rationals one has

```text
D_G(alpha beta)<=D_G(alpha)D_G(beta),
D_G(alpha beta)>=D_G(alpha)/D_G(beta),
D_G(alpha^(-1))=D_G(alpha).                           (3)
```

The first inequality follows by multiplying reduced fractions, and
the others follow by inversion. If `c in Z[i]` is conjugate-primitive,
then `D_G(c/bar c)=|c|`; for arbitrary nonzero Gaussian `c`, the same
quantity is at most `|c|`.

In particular unscaled reflection of a whole primitive cluster
requires `D_G(rho)=1`, so `rho` must be a Gaussian unit. Arbitrary
angular recentering followed by conjugation does not preserve its
integrality.

## 2. A lower bound for every pair midpoint in the reconstruction

Use the actual independent-block data and least-radius reconstruction
from [independent_block_reconstruction.md](independent_block_reconstruction.md):

```text
B=product_(S nonempty) H_S,        P=|B|,
A_i=product_(S contains i)H_S,    h_i=K_i A_i,
L=lcm_G(h_1,...,h_m)=BJ,          J in Z[i],
z_0=bar L,                       z_i=bar L h_i/bar h_i,
R=|L|=P|J|.
```

The blocks and their conjugates have disjoint prime support, each
`h_i` is conjugate-primitive, and the phase ratios are distinct.
Set `A_0=K_0=h_0=1`. For `0<=i<j<=m`, put

```text
B_sep(i,j)=product_(S separating i,j) H_S,
B_agree(i,j)=product_(S not separating i,j, S nonempty) H_S,
s_ij=|B_sep(i,j)|,
k_ij=|K_i K_j|,
P=s_ij |B_agree(i,j)|.                               (4)
```

The anchor `0` is outside every subset. Reflection about the
midpoint axis of `z_i,z_j` has coefficient

```text
rho_ij=z_i z_j/R^2
      =(h_i h_j/L)/conjugate(h_i h_j/L).              (5)
```

The core part of `h_i h_j/L` is
`A_i A_j/B`, the product of the both-inside blocks divided by
the product of the both-outside blocks. Its conjugate quotient
has reduced Gaussian denominator modulus `|B_agree(i,j)|`.
This uses the full independent, conjugate-coprime block hypothesis.
The remaining factor is `K_i K_j/J`. By (3),

```text
D_G(rho_ij) >= |B_agree(i,j)|/(k_ij |J|)
             =P/(s_ij k_ij |J|).                    (6)
```

We next retain an actual nonzero integer residue. Let
`g=gcd_G(h_i,h_j)` and

```text
v_ij=h_i bar(h_j)/N(g) in Z[i].
```

Distinctness of the phases says `Im v_ij!=0`, so
`|Im v_ij|>=1`. The product of all blocks containing both rows
divides `g`. Therefore

```text
|v_ij|=|h_i h_j|/|g|^2 <= k_ij s_ij,
|z_i-z_j|=2R |Im v_ij|/|v_ij|
           >=2R/(k_ij s_ij).                        (7)
```

Now multiply the original cluster and its reflection by any common
`gamma` making the union integral. Its radius is `R'=|gamma|R`.
Every containing arc has length at least the chord of the scaled
original pair. Its normalized arc constant `C_ref` consequently
satisfies

```text
C_ref >= sqrt(|gamma|) |z_i-z_j|/sqrt R
      >= sqrt(D_G(rho_ij)) 2sqrt R/(k_ij s_ij)
      >= 2P/(k_ij s_ij)^(3/2).                       (8)
```

The factor `|J|` cancels exactly between (6) and the actual radius
`R=P|J|`. Thus there is no uncharged correction or common-denominator
cost in (8). Any subsequent primitive normalization leaves the same
lower bound, by the minimality assertion in Section 1.

For the anchor-row reflection this specializes to

```text
C_ref >= 2P/(|K_i A_i|)^(3/2).                       (9)
```

This conclusion uses integer reconstruction and `Im v_ij` being
a nonzero integer. It does not treat a formal real cofactor merely
as if it were integral.

## 3. The growing fair core makes the reflected constant diverge

Suppose the number `m` of nonanchor rows tends to infinity, and write

```text
w_S=log N(H_S),          W_c=sum_(S subset [m]) w_S,
w_empty=w_(empty),      E=sum_i log|K_i|.
```

Assume the near-uniform block estimate and correction budget

```text
max_S |w_S-W_c/2^m|=o(W_c/2^m),
E=o(W_c),               W_c -> infinity.             (10)
```

These are weaker than the corresponding output already supplied by
central truncation. Each pair is separated by exactly half the
subsets, none of them empty. Hence, uniformly over every pair,

```text
log s_ij=W_c/4+o(W_c),
log k_ij=o(W_c),
log P=(W_c-w_empty)/2=W_c/2+o(W_c).
```

Also `1<=|J|<=exp E`, so `log R=W_c/2+o(W_c)`. Taking logarithms
in (8) proves

```text
log C_ref >= W_c/8+o(W_c)
           =(1/4-o(1)) log R,
C_ref >= R^(1/4-o(1)).                              (11)
```

Thus no pair-midpoint reflection can be adjoined to these reconstructed
clusters and made primitive integral while staying on an endpoint arc
with any fixed constant. This includes a hypothetical family whose
original normalized constants tend to zero as in the reconstruction
theorem: the mandatory scaling and the nonzero chord residue together
already imply (11).

The conclusion concerns reflection of the whole point set and adjoining
it on one integer circle. It does not obstruct the separate algebraic
reflection of direction rows used to preserve invariant relations in
[complement_reflection_invariant_relations.md](complement_reflection_invariant_relations.md).
That argument does not assert preservation of a small-arc circle after
the present adjoining-and-scaling operation.

## 4. Remaining scope

An arbitrary rational reflection axis is governed by (2), but need
not have the special coefficient (5); no analogue of (11) is asserted
for every such axis. Nor does (11) force a third original point into
one pair's reflection integrality domain. The later
[domain-capacity theorem](pair_reflection_domain_capacity.md) rules out
forcing a positive fraction into one such domain: it contains only O_C(1)
source points. The literal third-point question is restricted there to
denominator norm `D<=C^4/4`; its general existence is not proved.

The new obstruction is therefore operational and quantitative: pair
midpoints cannot be used to double the actual balanced reconstructed
cluster at bounded endpoint scale by clearing their denominators.
It supplies no contradiction to the original cluster itself.

## Exact verification

`python3 docs/check_midpoint_reflection_scaling_obstruction.py` checks
actual Gaussian reconstructions with prime powers and conjugate
correction factors. It verifies the reduced reflection denominator,
primitive integral reflected union, bounds (6)--(8), and the exact
pair-residue identities using integer arithmetic. The asymptotic
conclusion (11) follows from the displayed inequalities, not from a
finite search.

Independent reviews checked the primitive union, arbitrary correction
supports, integer pair residue, and conversion to the actual reconstructed
radius. No new general point-count estimate is claimed.
