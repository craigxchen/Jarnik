# Minimal-counterexample balance and polarization

## Status

This note proves a quantitative global consequence of the endpoint hypothesis.
It does not prove the uniform bound.

For a hypothetical sequence of endpoint clusters with cardinality tending to
infinity, almost all conductor weight must lie on prime-power threshold layers
that split the cluster almost exactly in half. Moreover, on every prime block
with small average defect, almost all points are polarized near exponent `0` or
near exponent `e`.

This justifies reducing the remaining problem to an almost-binary balanced
orientation code before applying affine-cube or arithmetic rigidity arguments.

---

## 1. Setup

Let

```text
A = {a_1,...,a_M} subset product_j {0,...,e_j}
```

be exponent vectors for distinct Gaussian lattice points of common norm

```text
N = product_j p_j^(e_j),
H = log N.
```

Assume their angles lie in an interval of width

```text
delta = C exp(-H/4).
```

For two points define their conductor distance

```text
D(i,k) = sum_j |a_(i,j)-a_(k,j)| log p_j.
```

The corresponding norm-one Gaussian quotient has height `D(i,k)`. The elementary
nonzero-Gaussian-numerator estimate and endpoint closeness give

```text
D(i,k) >= H/2 - O_C(1)                                 (1.1)
```

for every distinct pair.

---

## 2. Threshold layers

For a fixed prime block `p^e`, define

```text
k_t = #{i : a_i >= t},       1 <= t <= e.
```

The scalar identity

```text
|a_i-a_k|
  = sum_(t=1)^e 1_{exactly one of a_i,a_k is >= t}
```

implies

```text
sum_(i<k) |a_i-a_k|
  = sum_(t=1)^e k_t (M-k_t).                            (2.1)
```

After summing over prime blocks,

```text
sum_(i<k) D(i,k)
  = sum_(j,t) k_(j,t) (M-k_(j,t)) log p_j.              (2.2)
```

For every integer `0 <= k <= M`,

```text
k(M-k) <= floor(M^2/4).                                 (2.3)
```

Define the global cut defect

```text
Def(A)
  = sum_(j,t)
      [floor(M^2/4)-k_(j,t)(M-k_(j,t))] log p_j.
```

Every summand is nonnegative.

---

## 3. Global near-equality

The upper bound (2.3) gives

```text
sum_(i<k) D(i,k)
  <= floor(M^2/4) H.
```

The pairwise endpoint lower bound (1.1) gives

```text
sum_(i<k) D(i,k)
  >= binom(M,2) [H/2-O_C(1)].
```

Therefore

```text
Def(A)
  <= floor(M^2/4) H - binom(M,2) H/2 + O_C(M^2).        (3.1)
```

For even `M`, the main difference on the right is exactly

```text
M H / 4.
```

For odd `M` it has the same order. Uniformly,

```text
Def(A) = O(M H +_C M^2).                                (3.2)
```

Relative to the total possible scale `M^2 H`,

```text
Def(A)/(M^2 H) = O(1/M +_C 1/H).                        (3.3)
```

Thus a hypothetical family with `M -> infinity` and `H -> infinity` is forced
arbitrarily close to the Plotkin equality configuration globally.

---

## 4. Almost every weighted layer is balanced

For simplicity take `M` even and write `m=M/2`. Then

```text
m^2-k(M-k) = (k-m)^2.
```

Fix `epsilon>0`. Call a threshold layer `epsilon`-unbalanced if

```text
|k-m| >= epsilon M.
```

Every such layer contributes at least

```text
epsilon^2 M^2
```

to the unweighted defect. If `W_bad(epsilon)` is the total conductor weight of
those layers, then

```text
epsilon^2 M^2 W_bad(epsilon) <= Def(A).
```

Using (3.2),

```text
W_bad(epsilon)/H
  = O(1/(epsilon^2 M) +_C 1/(epsilon^2 H)).              (4.1)
```

Hence for every fixed `epsilon>0`, the proportion of conductor height carried by
`epsilon`-unbalanced layers tends to zero.

The odd-cardinality version uses the two balanced sizes `(M-1)/2` and `(M+1)/2`
and gives the same conclusion.

---

## 5. Polarization inside one prime block

The threshold counts

```text
k_1 >= k_2 >= ... >= k_e
```

form a monotone chain.

Suppose, for parameters `epsilon,eta>0`, all but at most `eta e` values of `t`
satisfy

```text
M/2-epsilon M <= k_t <= M/2+epsilon M.                  (5.1)
```

Then:

- at least `M/2-epsilon M` points satisfy

  ```text
  a_i >= (1-eta)e;
  ```

- at least `M/2-epsilon M` points satisfy

  ```text
  a_i <= eta e;
  ```

- consequently at most `2 epsilon M` points can remain in the middle interval

  ```text
  eta e < a_i < (1-eta)e.
  ```

Proof: if more than `eta e` early thresholds had count above
`M/2+epsilon M`, they would already be exceptional in (5.1). Thus by threshold
`eta e` the count has fallen below the upper boundary. The analogous argument
from the far end gives the lower boundary at `(1-eta)e`. Monotonicity then gives
the displayed endpoint populations.

So a low-defect prime block acts, on almost all points, like one balanced binary
orientation variable:

```text
roughly half near pi^e,
roughly half near conjugate(pi)^e,
very few in the middle.
```

---

## 6. Global binary reduction

Combining Sections 4 and 5 gives the following asymptotic description of every
hypothetical unbounded endpoint family.

After discarding conductor weight `o(H)` and a point fraction `o(1)` on the
remaining good blocks, the exponent code is a balanced binary orientation code.
Every good prime block contributes a sign column

```text
epsilon_j(i) in {+1,-1}
```

with approximately half the points in each sign, and the point angles have the
form

```text
alpha_i
  = constant + sum_j epsilon_j(i) beta_j + small error,
```

where

```text
beta_j = e_j arg(pi_j)
```

modulo the harmless choice of orientation.

This is exactly the regime analyzed by the affine-cube/Walsh note

```text
docs/affine_cube_walsh_rigidity.md.
```

---

## 7. What this proves and what it does not

The balance theorem rules out attempts to obtain a contradiction merely from
one visibly unbalanced prime or from a positive average determinant surplus.
Any counterexample is already forced into the balanced equality regime.

Walsh rigidity then eliminates the high-additive-energy part of that regime.
The only surviving case is an almost-binary, balanced, additive-Sidon-like code
whose Gaussian angle sum is exponentially concentrated.

Thus the global contradiction program has now reduced the endpoint conjecture
to a sparse arithmetic statement rather than a general weighted box theorem.
