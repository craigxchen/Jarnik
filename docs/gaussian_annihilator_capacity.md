# Gaussian annihilator capacity

## Purpose

This note introduces a new global arithmetic tool for endpoint-scale lattice-point clusters on circles.

The recurring failure of earlier approaches was that they projected the actual Gaussian divisor configuration onto pairwise distances, cut metrics, Fourier magnitudes, additive energy, or one exceptional set. Those projections admit abstract extremizers that need not be realizable by one Gaussian conductor.

The present tool keeps the full realization.

It works in the integral group algebra of the exponent lattice and measures whether a signed combination of actual divisor monomials can have simultaneously:

1. high-order cancellation of the lifted angular phases;
2. small weighted support width in the Gaussian-prime coordinates;
3. controlled coefficient mass.

A positive capacity certificate forces a contradiction by comparing an analytic upper bound with a Gaussian-integer lower bound.

The framework does not yet prove the uniform endpoint theorem. It creates a precise object whose existence would prove it, and whose nonexistence has a dual structural meaning.

---

## 1. Realized divisor monomials

Write the split part of the common norm as

```text
N = product_j p_j^(e_j),
```

and choose one Gaussian prime `pi_j` above each `p_j`. Put

```text
rho_j = pi_j / conjugate(pi_j).
```

The elements `rho_j` are norm-one elements of `Q(i)^*`.

After choosing a root point, every lattice-point ratio has the form

```text
u_a = product_j rho_j^(a_j),
```

for an exponent vector `a` in an integer box. Translating all exponents by one fixed vector only multiplies every `u_a` by a common norm-one factor, so support widths are translation-invariant.

Let `A` be the exponent set of a cluster. A signed group-algebra element supported on `A` is

```text
P = sum_(a in A) c_a [a],
```

with coefficients `c_a in Z`. Its arithmetic evaluation is

```text
P(rho) = sum_(a in A) c_a u_a.
```

This is a linear combination of the actual Gaussian divisor ratios, not an abstract model.

---

## 2. Gaussian support width

For the support of `P`, define

```text
L_j(P) = min {a_j : c_a != 0},
U_j(P) = max {a_j : c_a != 0}.
```

The conductor-weighted support width is

```text
W(P) = sum_j (U_j(P)-L_j(P)) log p_j.
```

Multiply `P(rho)` by

```text
D(P) = product_j pi_j^(-L_j(P)) conjugate(pi_j)^(U_j(P)).
```

Every term becomes a Gaussian integer:

```text
D(P) P(rho)
 = sum_a c_a product_j
     pi_j^(a_j-L_j(P))
     conjugate(pi_j)^(U_j(P)-a_j).
```

Therefore, if `P(rho) != 0`, then

```text
|D(P) P(rho)| >= 1.
```

Since

```text
|D(P)| = exp(W(P)/2),
```

we obtain the exact arithmetic separation lemma

```text
|P(rho)| >= exp(-W(P)/2).                 (2.1)
```

This lower bound is global. It depends on the coordinate width of the whole support, not on a sum of independently estimated pairwise heights.

There is a second branch:

```text
P(rho) = 0.                                (2.2)
```

That is an exact additive relation among realized Gaussian divisor monomials. Such zero relations must be treated as arithmetic degeneracies and reduced by minimal support, rather than discarded.

---

## 3. Lifted angular coordinates

Suppose the cluster lies in an angular interval of width `delta`.

Choose a center angle `alpha_0` and real lifts `t_a` satisfying

```text
u_a = exp(i(alpha_0+t_a)),
|t_a| <= delta.
```

The exponent formula is only valid modulo `2 pi`, so attach a winding integer `n_a` and define the lifted exponent vector

```text
b_a = (a,n_a) in Z^(r+1)
```

such that

```text
t_a = 2 sum_j a_j theta_j - 2 pi n_a - alpha_0.   (3.1)
```

This extra winding coordinate is essential. It prevents false moment cancellations caused by silently choosing incompatible logarithm branches.

---

## 4. Tensor annihilators

For `q >= 0`, call `P` a lifted tensor annihilator of order `q` when

```text
sum_a c_a b_a^(tensor m) = 0
```

for every `0 <= m < q`.

The `m=0` condition is

```text
sum_a c_a = 0.
```

The `m=1` condition says that both exponent coordinates and winding numbers cancel:

```text
sum_a c_a a = 0,
sum_a c_a n_a = 0.
```

Tensor cancellation implies phase-moment cancellation for every real linear functional on the lifted lattice. Applying the functional in (3.1) gives

```text
sum_a c_a t_a^m = 0
```

for every `0 <= m < q`.

Taylor expansion of `exp(i t)` therefore yields

```text
|P(rho)|
 <= ||c||_1 exp(delta) delta^q / q!,       (4.1)
```

where

```text
||c||_1 = sum_a |c_a|.
```

This is the analytic side of the tool.

Affine cubes are one special source of tensor annihilators, but the definition is much broader. Arbitrary signed circuits, finite-difference schemes, and higher-order cubature identities are all included.

---

## 5. Capacity

Define the annihilator gain

```text
Gain(P)
 = q log(1/delta)
   - W(P)/2
   - log ||c||_1
   + log(q!)
   - delta.
```

At the endpoint scale

```text
delta = C exp(-H/4),
H = log N,
```

this is

```text
Gain(P)
 = q H/4 - W(P)/2 - log ||c||_1 + log(q!) - O_C(q).
```

If `P` is a tensor annihilator and

```text
Gain(P) > 0,
```

then (4.1) is strictly smaller than the nonzero lower bound (2.1). Hence

```text
P(rho) = 0.                                (5.1)
```

If one can also rule out exact zero relations of the relevant support type, a positive-gain annihilator gives an immediate contradiction.

Define the Gaussian annihilator capacity of the realized cluster by

```text
Cap(A;rho)
 = sup_P Gain(P),
```

where `P` ranges over nonzero integer lifted tensor annihilators supported on `A`.

The endpoint theorem would follow from a dichotomy of the following form:

> For every sufficiently large cluster, either `Cap(A;rho)>0`, or the absence of positive-capacity annihilators forces a bounded-complexity arithmetic structure that admits descent.

This is the intended primal-dual formulation.

---

## 6. Why this tool is different

### 6.1 It remembers the actual Gaussian realization

The arithmetic lower bound is taken after evaluating at

```text
rho_j = pi_j / conjugate(pi_j).
```

Abstract balanced codes or random sign systems do not automatically satisfy the Gaussian-integer separation lemma.

### 6.2 It is genuinely global

The support width is measured once for the whole signed relation. It is not a nonnegative sum of pairwise cut distances, so the max-cut barrier does not apply automatically.

### 6.3 It handles winding numbers

The lifted coordinate `n_a` records the branch of the angle. This was missing from several informal additive reductions.

### 6.4 It separates analytic and arithmetic complexity

The three competing quantities are explicit:

```text
annihilation order q,
weighted support width W(P),
coefficient mass ||c||_1.
```

The research problem becomes an optimization problem rather than a vague request for a stronger inequality.

---

## 7. Exact-zero relations

The branch `P(rho)=0` cannot be ignored. Multiplicative independence of the `rho_j` does not imply linear independence of their Laurent monomials.

The correct response is to choose a zero relation of minimal support and study its Newton polytope.

Possible reduction principles are:

1. if the support lies in a proper affine sublattice, reduce exponent rank;
2. if all terms share a Gaussian monomial factor, divide it out;
3. if the Newton polytope decomposes as a Minkowski sum, attempt factorization in the group algebra;
4. use nondegenerate `S`-unit equation bounds only after the support size has been fixed by the annihilator construction;
5. exploit the fact that all support monomials themselves lie in one tiny angular interval, which is far stronger than an arbitrary vanishing `S`-unit sum.

A complete theory should therefore classify minimal zero relations with small angular diameter and small weighted Newton width.

---

## 8. Dual obstruction

Suppose no low-coefficient positive-gain annihilator exists.

The tensor-annihilation conditions form an integer linear system whose columns are the lifted moment vectors

```text
(1,b_a,b_a^(tensor 2),...,b_a^(tensor(q-1))).
```

The absence of a short integer kernel vector should have a geometry-of-numbers dual: there exists a low-complexity polynomial functional that separates the lifted points.

Such a dual certificate would say that the cluster is unisolvent for low-degree polynomials with respect to a weighted lattice norm. That is a rigid global condition and may force one of:

```text
large affine rank,
large weighted simplex volume,
heavy conductor coordinates,
large winding dispersion.
```

Each of these has a potential descent or arithmetic consequence.

This dual side is what allows the tool to address sparse, cube-free configurations rather than merely rediscover affine cubes.

---

## 9. First concrete theorem to attack

The first target should be:

> **Positive-capacity or weighted-simplex theorem.** Fix `M`. For every `M` lifted exponent vectors in an endpoint cluster, either:
>
> 1. there is an integer tensor annihilator `P` of order at least `3` with
>    ```text
>    W(P) <= H,
>    log ||c||_1 = o(H),
>    ```
>    and hence positive capacity; or
> 2. some bounded subset of the lifted vectors has weighted simplex determinant at least `exp(cH)`, forcing a bounded set of prime or winding coordinates to carry a positive proportion of the conductor.

Order `3` is the first meaningful threshold: at endpoint scale it supplies analytic decay `exp(-3H/4)`, while the worst full-box arithmetic lower bound is `exp(-H/2)`.

The margin is `H/4`, enough to absorb subexponential coefficient growth.

The statement must be proved with the actual weighted lattice norm, not ordinary Euclidean rank.

---

## 10. Falsification tests

Before treating this as the main route, test the capacity on:

1. the known four-point endpoint example;
2. random squarefree Gaussian conductors;
3. balanced Hadamard-style sign systems realized with artificial angles;
4. sparse Sidon-like exponent sets;
5. configurations with one dominant prime power;
6. minimal exact zero relations among Gaussian divisor monomials.

The tool should be abandoned if large genuine endpoint-like configurations can have simultaneously:

```text
Cap(A;rho) <= 0,
no heavy coordinates,
no low-rank support,
no reducible exact zero relation.
```

---

## 11. Computational implementation

A search routine can be built around the following data:

```text
input:
  exponent vectors a,
  winding integers n,
  prime weights log p,
  maximum order q,
  coefficient bound B.
```

For each `q`:

1. build the integer moment matrix through degree `q-1`;
2. compute a short integer kernel basis using Smith normal form and LLL;
3. enumerate short combinations;
4. compute `W(P)`, `||c||_1`, and `Gain(P)`;
5. evaluate `P(rho)` exactly in `Q(i)` when feasible;
6. classify positive-capacity certificates and exact zero relations.

This would be a genuinely useful discovery tool: it searches for global arithmetic cancellations that all previous pairwise and magnitude-based methods were incapable of seeing.
