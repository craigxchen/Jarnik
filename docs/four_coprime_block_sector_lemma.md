# Four-coprime-block sector lemma

## Status

This note proves a global endpoint obstruction that does not reduce to a
pairwise cut bound.

At most four pairwise-coprime Gaussian conductor blocks can have their
arguments in one endpoint-scale sector while their norm product divides the
common conductor.  The proof multiplies all pairwise determinant lower bounds;
for five blocks the conductor exponents are incompatible.

The lemma explains structurally why four-point configurations survive many of
the preceding arguments.  It also sharpens the weighted-sunflower obstruction:
five disjoint petals after one common core are already impossible.

---

## 1. Statement

Let

```text
Gamma_1,...,Gamma_k in Z[i]
```

be nonzero Gaussian integers with pairwise disjoint nonunit Gaussian prime
supports.  Put

```text
n_j = Norm(Gamma_j).
```

Assume

```text
product_j n_j <= N.                                  (1.1)
```

Suppose their arguments admit real lifts contained in one interval of width

```text
delta <= C N^(-1/4).                                 (1.2)
```

Assume also that no two `Gamma_i,Gamma_j` lie on exactly the same real ray.
For primitive pairwise-coprime split-prime blocks this last condition is
automatic, apart from units and rational factors that can be removed first.

Then, for sufficiently large `N` depending only on `C`,

```text
k <= 4.                                               (1.3)
```

---

## 2. Pairwise integer determinant

For `i != j`, put

```text
Delta_ij = |Im(conjugate(Gamma_i) Gamma_j)|.
```

This is a positive integer, hence

```text
Delta_ij >= 1.                                        (2.1)
```

On the other hand,

```text
Delta_ij
  = |Gamma_i| |Gamma_j|
      |sin(arg Gamma_j-arg Gamma_i)|
  <= sqrt(n_i n_j) delta.                             (2.2)
```

Therefore

```text
n_i n_j >= delta^(-2)
          >= C^(-2) N^(1/2).                          (2.3)
```

The important point is that the same disjoint block norms occur in every pair.

---

## 3. Multiply all pair inequalities

Multiplying (2.3) over all unordered pairs gives

```text
(product_j n_j)^(k-1)
  >= C^(-k(k-1)) N^(k(k-1)/4).                       (3.1)
```

Using (1.1),

```text
N^(k-1)
  >= C^(-k(k-1)) N^(k(k-1)/4).                       (3.2)
```

Equivalently,

```text
N^((k-1)(1-k/4)) >= C^(-k(k-1)).                     (3.3)
```

For `k>=5`, the exponent on the left is negative.  Hence (3.3) fails for all
sufficiently large `N` depending only on `C`.

Thus `k<=4`.

---

## 4. Why this is genuinely global

One pair only gives

```text
n_i n_j >= N^(1/2)/C^2,
```

which is exactly critical and allows arbitrarily large individual block norms.
The gain appears only after multiplying the inequalities for every pair:
`k` disjoint blocks contribute to `k(k-1)/2` pair constraints, while the total
conductor pays for each block only once.

This is the type of nonlocal accounting missing from the earlier determinant,
Markov, and entropy estimates.

The exponent comparison is

```text
pair requirements:  k(k-1)/4,
conductor budget:    k-1.
```

They agree at `k=4` and become inconsistent at `k=5`.

---

## 5. Sunflower corollary

Let a rooted exponent-support family contain five sets

```text
S_j = K union P_j,
```

where the petals `P_1,...,P_5` are pairwise disjoint.  Factor the common
Gaussian block associated with `K`.  The remaining petal products

```text
Gamma_1,...,Gamma_5
```

have pairwise disjoint prime support.

All five original circle points lie in one endpoint arc, and removing the same
common block rotates and rescales all of them equally.  Hence the five petal
block arguments remain in one interval of width `C N^(-1/4)` in the original
conductor normalization.

Their norm product divides `N`.  The four-block lemma therefore gives a
contradiction for sufficiently large `N`.

Thus:

> An endpoint cluster contains no five-petal exact weighted sunflower after a
> common Gaussian core is removed.

This improves the previous elementary weight-only obstruction, which required
six petals.

---

## 6. Near-equality model

The lemma also explains the persistent four-point examples.  A model family

```text
Gamma_j = (X+j)+i,
```

has pairwise determinant of order one and norms of order `X^2`.  The product of
`k` norms is of order `X^(2k)`, so the endpoint angular scale is of order
`X^(-k/2)`.  The pairwise angle spacing is of order `X^(-2)`.

The inequality

```text
X^(-2) <= X^(-k/2)
```

is compatible exactly for `k<=4`; `k=4` is the critical case.

This toy calculation is not itself a circle construction, but it isolates the
same exponent balance as the rigorous proof.

---

## 7. Remaining extraction problem

To use the lemma for the full theorem, one must extract five pairwise-disjoint
Gaussian blocks whose arguments are all endpoint-close.

An exact five-petal sunflower supplies them immediately.  General balanced
orientation codes need not contain such a sunflower, so ordinary set-system
extraction is insufficient.

The prime-reveal/renormalization program can now be targeted more sharply:

```text
large reveal entropy
  -> five disjoint arithmetic completion blocks
  -> four-block sector contradiction.
```

Unlike the earlier support-only matching goal, the blocks must carry their
actual Gaussian products and must be pairwise angle-close.  The present lemma
is the arithmetic endpoint once such a five-block extraction is obtained.
