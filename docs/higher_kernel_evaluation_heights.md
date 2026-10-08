# Evaluation gcds for the higher cut kernels

This note extends the evaluation-content calculation from the balanced-cut
kernel to every active depth of the cut hierarchy. The lower gcd bound comes
from the monomial support forced by the vanishing cuts. It does not require a
graph-monomial spanning theorem, so cancellations between invariant
polynomials cannot reduce the bound. Graph monomials are used only to attain
the cut order and to prove the matching upper gcd bound.

## 1. Cut kernels and the forced monomial support

Let `m=2q` and let `V_m` be the integral lattice of simultaneous `SL_2`
invariants that are homogeneous of degree two in each binary row. For
`0<=k<=floor(m/6)`, let `K_k` be the subspace whose restrictions vanish when
any at least `q-k+1` rows are collapsed to one direction. Thus `K_0=V_m`;
for `k>=1` this is the hierarchy in
[all_cut_invariant_relation_hierarchy.md](all_cut_invariant_relation_hierarchy.md).
For `k=0`, cut restrictions of size at least `q+1` vanish automatically by
weight balance. The Specht dimensions of these kernels are established in
[kernel_specht_filtration.md](kernel_specht_filtration.md), but are not needed
for the evaluation-content proof here.

Use coordinates `(P_i,Y_i)` on row `i`, with the collapsed direction given
by `P_i=0`. The chart used to define a cut restriction sets rows in the cut
to `(0,1)` and the other rows to `(1,z_j)`. By rowwise homogeneity, vanishing
on this chart extends to all choices of the outside rows on the dense set
where their `P_j` are nonzero, and hence to the polynomial identity

```text
Q(P_i=0 for every i in T) = 0,       |T|>=q-k+1.          (1)
```

The diagonal `SL_2` torus forces every monomial of `Q` to have total
`P`-degree `m`. Write a monomial as

```text
P^a Y^(2-a) = product_i P_i^a_i Y_i^(2-a_i),
0<=a_i<=2,          sum_i a_i=m.                           (2)
```

If `a_i=0` for every `i` in some `T` of size `q-k+1`, that monomial
survives the specialization (1). Its outside exponent vector is unique, so
its coefficient cannot cancel against a different monomial. Therefore every
nonzero monomial coefficient in `Q` has `P`-support meeting every such `T`.
Equivalently, its support has at least `q+k` labels.

For any subset `S` of size `s`, the support condition and the total-degree
condition give

```text
sum_(i in S) a_i >= s-q+k,
sum_(i in S) a_i >= m-2(m-s)=2s-m,
sum_(i in S) a_i >= 0.
```

Thus every monomial of every `Q in K_k` has inside `P`-degree at least

```text
g_k(s)=max(0, 2s-m, k+s-q).                               (3)
```

For `k=0`, the support condition is automatic from `sum a_i=m` and
`a_i<=2`, and (3) reduces to `g_0(s)=max(0,2s-m)`. The proof is
coefficientwise, so arbitrary linear combinations can only keep or increase
the order.

## 2. Sharp cut order from triangle and pair graphs

Assume `k>=1` and `m>=6k`. Form a bracket graph from `2k` disjoint
triangles and `q-3k` disjoint pairs, using one bracket on each triangle
edge and a squared bracket on each pair. Its polynomial has degree two in
every row. Its largest set containing no internal graph edge has size

```text
2k+(q-3k)=q-k.
```

Consequently the graph polynomial belongs to `K_k`. At `k=0`, use a perfect
matching with every edge doubled; it belongs to `V_m`.

For either graph, let `I` be the number of internal edge factors in a cut
`S`, counted with multiplicity, and let `X` be its number of crossing edge
factors. The graph is 2-regular, so `2s=2I+X`. With
`L=q-k` components, each component contributes at most two crossing factors,
and the smaller side of the cut has at most `min(s,m-s)` vertices. Hence

```text
X<=2 min(L,s,m-s),
I>=s-min(L,s,m-s)=g_k(s).                               (4)
```

Equality can be attained for every `s`. If `s<=L`, put one vertex from
each of `s` components inside. If `L<=s<=m-L`, split every component and
add the extra inside vertices among the triangles; there are `2k` triangles,
so this covers the full interval of width `2k`. For `s>=m-L`, apply the
construction to the complement. When `k=0`, splitting pairs and then
complementing gives `g_0(s)`.

This establishes sharpness of (3) among actual graph witnesses. It does not
assume that those witnesses span `K_k`.

## 3. Exact evaluation gcd bounds

Use the core factorization from
[endpoint_central_truncation.md](endpoint_central_truncation.md):

```text
P_i=K_i product_(S containing i) H_S,
n_S=Norm(H_S),
|Delta_ij|=b_ij product_(S containing i,j) n_S,
B=product_(i<j) b_ij.                                  (5)
```

The Gaussian blocks `H_S` and their conjugates have disjoint prime support;
corrections and bracket residues are included in the positive integers
`b_ij`. Empty-set factors have exponent zero throughout. Let `G_k` be the
positive gcd of the integer evaluations of the saturated integral lattice
`K_k` at the given distinct rows, and set

```text
D_k=product_S n_S^g_k(|S|).                              (6)
```

Then

```text
D_k | G_k | D_k B^2.                                    (7)
```

For the lower divisibility, write `P_i=X_i+iY_i`. The change
`(X_i,Y_i)->(P_i,Y_i)` is a determinant-one shear, so invariance gives
`Q(P,Y)=Q(X,Y)`, an ordinary integer. By (3), each monomial term in this
evaluation contains at least `g_k(|S|)` copies of `P_i` across the rows of
each `S`. The factorization (5) therefore makes
`product_S H_S^g_k(|S|)` divide every term and every evaluation. Its
conjugate also divides the real integer evaluation; coprimality of the
oriented blocks gives the norm divisor `D_k`.

For the upper divisibility, consider a rational prime `p`. If `p` divides
one core norm `n_S`, choose the graph witness in Section 2 attaining
`g_k(|S|)` for that cut. The graph evaluation has exactly that core
valuation at `p`. Its bracket-residue contribution is at most
`2 sum_(i<j) v_p(b_ij)`, since each edge has multiplicity at most two.
If `p` divides no core norm, any graph witness has only this residue
contribution. The core supports are disjoint, so there is at most one such
`S` for a given rational prime. Since `G_k` divides every graph evaluation,
the primewise upper bound in (7) follows.

## 4. Primitive evaluation height and relation-lattice covolume

Let

```text
A_0=m*2^(m-2)-q*binomial(m,q),
Delta_k=sum_(s=0)^m binomial(m,s) max(0,k-|s-q|),
A_k=A_0-Delta_k.                                        (8)
```

The identity

```text
g_k(s)-g_0(s)=max(0,k-|s-q|)
```
shows that the leading logarithmic evaluation gcd exponent increases by
`Delta_k` from `V_m` to `K_k`. In a fixed full core profile with
`log n_S=w+O(epsilon w)` for every `S`, and
`log B=o(w)`, (7) gives

```text
log G_k = (sum_s binomial(m,s) g_k(s)) w
          +O_m(epsilon w)+O(log B).                      (9)
```

For fixed `m,k`, every bracket has logarithmic size
`2^(m-2)w+O_m(epsilon w)+O(log B)`. A graph witness in `K_k` has exactly
`m` bracket factors counted with multiplicity, and a fixed integral basis
changes coordinate norms by only a constant depending on `m,k`. Hence the
raw evaluation vector on `K_k` has logarithmic height
`m*2^(m-2)w+O_m(epsilon w)+O(log B)`. After division by its gcd, the
primitive evaluation vector has logarithmic height

```text
log ||f_k|| = A_k w+O_m(epsilon w)+O(log B)+O_m(1).      (10)
```

The evaluation functional on `K_k` is nonzero because the graph witness is
nonzero at distinct directions. If `d_k=rank(K_k)` and
`Lambda_k={x in Z^d_k : f_k dot x=0}` in a fixed integral basis, then
`Lambda_k` has rank `d_k-1` and covolume exactly `||f_k||_2`. Thus its
logarithmic covolume is also
`A_k w+O_m(epsilon w)+O(log B)+O_m(1)`.

For reference, the first values of the coefficient `A_k` are

| `m` | `k` | `A_k` |
|---:|---:|---:|
| 6 | 1 | 16 |
| 8 | 1 | 162 |
| 10 | 1 | 1,048 |
| 12 | 1 | 5,820 |
| 12 | 2 | 3,312 |
| 18 | 3 | 357,528 |

The exact checker
[`check_higher_kernel_evaluation_heights.py`](check_higher_kernel_evaluation_heights.py)
checks the monomial-support minimum, the graph cut minimum, the binomial
identities, and two exact split-prime valuation examples. No representation
formula for `K_k` is needed in the proof of the gcd exponents. These are
heights of the restricted evaluation functionals on `K_k`; the result does
not put short numerical relations into successive kernels and does not prove
a uniform endpoint bound.
