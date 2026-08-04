# Short-gap transition matching and prime-reveal renormalization

## Status

This note strengthens the ordered half-angle/Farey process.  Among a positive
fraction of the shortest consecutive gaps, the corresponding Gaussian
transition blocks are so large that their prime supports have fractional
matching number at most two.  Consequently one rational split prime switches
on a positive fraction of those transitions, and both of its orientation
classes contain a positive fraction of the cluster.

Dividing the two orientation classes by the corresponding Gaussian prime gives
a genuine circle-to-circle renormalization.  A hypothetical unbounded family at
one endpoint constant therefore produces unbounded families at arbitrarily
small endpoint constants and, for every fixed depth, an adaptive binary tree of
renormalized endpoint clusters.

This is rigorous in the odd squarefree model.  It does not by itself prove the
uniform endpoint bound: the adaptive prime labels may change from branch to
branch, so the resulting tree need not contain a fixed-coordinate Boolean
cube.

---

## 1. Short consecutive gaps

Use the rooted half-angle notation

```text
0=t_0<t_1<...<t_(M-1)<=T,
T<=C N^(-1/4).
```

For the transition from `t_i` to `t_(i+1)`, remove the maximal common
Gaussian divisor of the two half-angle divisors.  The remaining Gaussian
transition block is

```text
B_i = X_i+iY_i,
q_i = Norm(B_i),
```

and

```text
arg(B_i)=arg(A_(i+1))-arg(A_i).
```

Since the transition is nontrivial, `Y_i` is a nonzero integer.  Hence

```text
1 <= |Y_i|
   = sqrt(q_i) |sin(arg B_i)|
   <= sqrt(q_i) (t_(i+1)-t_i).                       (1.1)
```

At least

```text
L = floor((M-1)/2)
```

of the consecutive gaps satisfy

```text
t_(i+1)-t_i <= 2T/(M-1).                             (1.2)
```

For every such short transition,

```text
q_i >= (M-1)^2/(4T^2)
    >= (M-1)^2 N^(1/2)/(4C^2).                       (1.3)
```

Thus the shortest half of the gaps have residual Gaussian norm larger than the
ordinary half-conductor threshold by a factor quadratic in the cluster size.

---

## 2. Fractional matching bound

Let `D_i` be the rational-prime support of `B_i`.  Since the conductor is
squarefree,

```text
q_i = product_(p in D_i) p.
```

Let `lambda_i>=0` be a fractional matching on the short-transition support
hypergraph:

```text
sum_(i:p in D_i) lambda_i <= 1
```

for every conductor prime `p`.  Then

```text
product_i q_i^(lambda_i) <= N.                       (2.1)
```

Combining (1.3) and (2.1), with

```text
Lambda_f = sum_i lambda_i,
```

gives

```text
Lambda_f
  <= log N /
      ( (1/2)log N + 2log(M-1)-log(4C^2) ).          (2.2)
```

In particular,

```text
Lambda_f <= 2+o(1)                                   (2.3)
```

along any unbounded endpoint family.

By finite hypergraph linear-programming duality, the fractional vertex-cover
number has the same bound.  Hence there are nonnegative prime weights `x_p`
with

```text
sum_p x_p <= 2+o(1),

sum_(p in D_i) x_p >= 1
```

for every short transition `i`.

Summing the latter inequalities over the `L` short transitions shows that one
conductor prime occurs in at least

```text
L/(2+o(1)) >= (M-1)/4-o(M)                           (2.4)
```

of the short transition supports.

---

## 3. A prime with two large orientation classes

For a fixed prime `p`, let `V_p` be the number of adjacent transitions on
which its orientation switches, and let `n_0(p),n_1(p)` be the numbers of
cluster points in the two orientation classes.  Every switch alternates the
binary state, so

```text
V_p <= 2 min(n_0(p),n_1(p)).                         (3.1)
```

The prime supplied by (2.4) therefore satisfies

```text
min(n_0(p),n_1(p)) >= (M-1)/8-o(M).                  (3.2)
```

Thus both orientations contain a positive fraction of the original endpoint
cluster.

A slightly weaker version follows without the short-gap selection: summing the
ordinary pair separation over consecutive gaps gives

```text
sum_p V_p log p
  >= (M-1)((1/2)log N-2log C),                       (3.3)
```

and hence some prime switches at least `(M-1)/2-o(M)` times.  The short-gap
argument is stronger because it ties the frequent switching to residual blocks
of norm `M^2 N^(1/2+o(1))`.

---

## 4. Exact geometric renormalization

Fix one of the two orientation classes at `p`.  Every point in that class is
divisible by the same Gaussian prime `pi` above `p` (or by its conjugate).
Dividing every point by this common factor is injective and sends the class to
lattice points on the circle

```text
|z|^2 = N/p.
```

It is a similarity of ratio `p^(-1/2)`.  Hence an original arc of physical
length

```text
C N^(1/4)
```

becomes an arc of length

```text
C N^(1/4) p^(-1/2)
  = (C p^(-1/4)) (N/p)^(1/4).                        (4.1)
```

Therefore the renormalized endpoint constant is

```text
C' = C p^(-1/4).                                     (4.2)
```

By (3.2), each of the two children has cardinality at least

```text
M/8-o(M).                                             (4.3)
```

---

## 5. Self-improvement of a hypothetical bad family

Assume endpoint clusters of constant `C` have unbounded cardinality.  Apply the
preceding construction and choose either large orientation class.  For every
fixed integer `d`, iterating `d` times gives a descendant family with

```text
M_d >= 8^(-d) M-o(M),                                (5.1)
```

and

```text
C_d <= C 5^(-d/4),                                   (5.2)
```

because every odd split prime is at least `5`.

Since `d` is fixed while `M->infinity`, the descendant cardinalities still tend
to infinity.  Thus:

> Unboundedness at one positive endpoint constant implies unboundedness at
> every smaller positive endpoint constant.

More precisely, it produces for every fixed depth an adaptive full binary tree
of renormalized endpoint clusters, each internal node labelled by a split prime
and each child retaining a fixed positive fraction of the parent cluster.

---

## 6. Why the adaptive tree is the remaining obstruction

If the same prime were used at every node of a given level, the depth-`d`
renormalization tree would give a fixed-coordinate Boolean cube.  The existing
Walsh-character argument rules out depth five for sufficiently large
conductor.

The present construction is adaptive: the frequent-switch prime may depend on
the branch.  Abstract binary codes can support arbitrarily deep adaptive trees
without containing a fixed five-dimensional subcube.  The prime labels in
different branches may therefore drift, and the Walsh cancellation no longer
isolates five fixed Gaussian blocks.

This is the precise obstruction to closing the renormalization argument.  A
successful continuation must prove one of the following arithmetic statements:

1. a bounded-depth adaptive tree contains a strong subtree whose prime labels
   synchronize by level; or
2. branchwise drift creates five pairwise-coprime effective Gaussian blocks in
   one endpoint sector, contradicting the four-block sector lemma.

The second alternative is especially natural because labels in different
branches are supported on disjointly renormalized factor histories.  What is
missing is a canonical choice of branch compensator whose prime support stays
inside its own subtree.
