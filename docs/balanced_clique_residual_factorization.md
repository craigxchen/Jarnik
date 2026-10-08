# Balanced cliques and residual Gaussian factorization

## Status

This note gives a global extraction theorem for a hypothetical unbounded
endpoint family in the squarefree regime.

It proves that, after passing to any fixed number of points, every pairwise
chord can be divided by its forced common Gaussian conductor factor so that the
remaining Gaussian integer has norm `N^o(1)`.

For three extracted points the conductor splits into four pairwise-coprime
blocks, and the three residual chord equations become an exact small-coefficient
Gaussian relation.  Equivalently, the three pairwise products of the
nonconstant blocks are all `N^o(1)`-close to perfect squares.

This is a genuine global reduction.  It does not yet prove the uniform endpoint
bound: the remaining three-block system is an integer analogue of a balanced
`abc`/simultaneous-almost-square configuration, and the available Pell bounds do
not control its large coefficients.

---

## 1. Squarefree orientation distance

Let

```text
N = product_(p in S) p
```

be squarefree, with every odd prime in `S` split in `Z[i]`.  Choose one Gaussian
prime above each `p`.  A lattice point of norm `N` chooses one of the two
conjugate factors over every `p`.

For two points `z_i,z_j`, let

```text
P_ij = product of the rational primes at which their orientations differ,
G_ij = N/P_ij.
```

There is a Gaussian common divisor `g_ij` of `z_i,z_j` with

```text
Norm(g_ij)=G_ij.
```

Write

```text
z_i-z_j = g_ij w_ij,

m_ij = Norm(w_ij).
```

Then, exactly,

```text
|z_i-z_j|^2 = G_ij m_ij.                            (1.1)
```

If the points lie in an arc of length `C N^(1/4)`, then

```text
|z_i-z_j|^2 <= C^2 N^(1/2).
```

Since `m_ij` is a positive integer for distinct points, (1.1) gives

```text
P_ij >= N^(1/2)/C^2.                                (1.2)
```

Put

```text
H = log N,
D_ij = log P_ij,
K_C = 2 log C.
```

Then

```text
D_ij >= H/2-K_C.                                    (1.3)
```

The residual norm satisfies the complementary upper bound

```text
m_ij <= C^2 exp(D_ij-H/2).                          (1.4)
```

Thus the excess of the weighted orientation distance over `H/2` measures
exactly how much arithmetic remains after the forced common divisor is removed.

---

## 2. Global average excess

Let the cluster contain `M` points.  At one rational prime, suppose `k` points
choose one Gaussian orientation and `M-k` choose the conjugate orientation.
That prime contributes to exactly

```text
k(M-k) <= floor(M^2/4)
```

unordered pair distances.  Summing over all conductor primes gives

```text
sum_(i<j) D_ij <= floor(M^2/4) H.                   (2.1)
```

For each pair define

```text
E_ij = max(0,D_ij-H/2),
F_ij = max(0,H/2-D_ij).
```

Equation (1.3) gives `F_ij <= K_C`, while

```text
D_ij-H/2 = E_ij-F_ij.
```

Since

```text
floor(M^2/4)-binom(M,2)/2 <= M/4,
```

(2.1) yields

```text
sum_(i<j) E_ij
  <= M H/4 + binom(M,2) K_C.                        (2.2)
```

This is stronger than merely saying that the average cut is balanced: it
controls the positive excess that appears in the residual norm bound (1.4).

---

## 3. Fixed-size balanced clique extraction

Assume a hypothetical bad sequence has

```text
M -> infinity,
H -> infinity,
M/H -> 0.
```

The last condition is supplied by the already-proved sublogarithmic arc bound.

Fix an integer `t>=3`.  Call a pair bad if

```text
E_ij > ((t-1)/M) H.
```

By (2.2), the number of bad pairs is at most

```text
M^2/(4(t-1)) + o(M^2).                              (3.1)
```

If the graph of good pairs had no clique `K_t`, Turan's theorem would force at
least

```text
M^2/(2(t-1))-O(M)
```

bad pairs.  This contradicts (3.1) for sufficiently large `M`.

Therefore every hypothetical unbounded family contains, for each fixed `t`,
`t` points satisfying

```text
D_ij <= H/2 + ((t-1)/M)H                            (3.2)
```

for every pair in the selected set.  Combining (1.4) and (3.2),

```text
m_ij <= C^2 exp((t-1)H/M) = N^o(1).                 (3.3)
```

This is the balanced-clique theorem:

> After removal of the pairwise forced Gaussian conductor divisor, all chords
> in any fixed extracted clique have subpower Gaussian norm.

The exponent tends to zero with the cardinality of the original cluster, not
with the fixed extracted clique size.

---

## 4. Shared-prime compatibility of half-angle vectors

Let `(q_i,r_i)` be primitive half-angle vectors, with

```text
D_i=q_i^2+r_i^2,
Delta_ij=q_i r_j-q_j r_i.
```

The chord identity has the form

```text
d_ij D_i D_j = 4 N Delta_ij^2.                      (4.1)
```

If an odd prime `p` divides both `D_i,D_j`, occurs only once in squarefree `N`,
and divides neither the factor `4` nor `N/p`, then (4.1) forces

```text
p | Delta_ij.                                       (4.2)
```

This is now formalized in

```text
GaussianChain/HalfAngleSharedPrime.lean.
```

The theorem means that every common conductor prime places all half-angle
vectors containing that prime on one common isotropic line modulo `p`.  It is
the local compatibility needed for the global block factorization below.

---

## 5. Exact three-point block factorization

Take three points from a balanced clique.  Every squarefree conductor prime has
one of four orientation types:

1. all three points use the same orientation;
2. point 1 is the unique odd orientation;
3. point 2 is the unique odd orientation;
4. point 3 is the unique odd orientation.

Group the chosen Gaussian prime factors into pairwise-coprime Gaussian integers

```text
Gamma_0, Gamma_A, Gamma_B, Gamma_C
```

for these four types.  Up to a common unit, the points have the exact form

```text
z_1 = Gamma_0 Gamma_A conjugate(Gamma_B) conjugate(Gamma_C),
z_2 = Gamma_0 conjugate(Gamma_A) Gamma_B conjugate(Gamma_C),
z_3 = Gamma_0 conjugate(Gamma_A) conjugate(Gamma_B) Gamma_C.  (5.1)
```

Let

```text
G = Norm(Gamma_0),
A = Norm(Gamma_A),
B = Norm(Gamma_B),
C = Norm(Gamma_C).
```

Then

```text
G A B C = N.                                        (5.2)
```

Define the pair residuals

```text
W_AB = Gamma_A conjugate(Gamma_B)
       -conjugate(Gamma_A) Gamma_B,
```

and cyclically.  Subtracting the points in (5.1) gives

```text
z_1-z_2 = Gamma_0 conjugate(Gamma_C) W_AB,
z_2-z_3 = Gamma_0 conjugate(Gamma_A) W_BC,
z_3-z_1 = Gamma_0 conjugate(Gamma_B) W_CA.           (5.3)
```

The three chord vectors sum to zero, so after cancelling `Gamma_0`:

```text
conjugate(Gamma_C) W_AB
 + conjugate(Gamma_A) W_BC
 + conjugate(Gamma_B) W_CA = 0.                     (5.4)
```

The coefficients `Gamma_A,Gamma_B,Gamma_C` are pairwise coprime, while by the
balanced-clique theorem

```text
Norm(W_AB), Norm(W_BC), Norm(W_CA) <= N^o(1).        (5.5)
```

Thus every hypothetical unbounded family produces an exact three-term Gaussian
relation with pairwise-coprime large blocks and subpower coefficients.

---

## 6. Simultaneous almost squares

Write

```text
Gamma_A conjugate(Gamma_B) = X_AB+i Delta_AB.
```

Then

```text
A B = X_AB^2+Delta_AB^2,
Norm(W_AB)=4 Delta_AB^2.                             (6.1)
```

Equations (5.5) and (6.1) give

```text
A B = X_AB^2+N^o(1),
B C = X_BC^2+N^o(1),
C A = X_CA^2+N^o(1).                                (6.2)
```

In addition, endpoint separation for the three point pairs gives

```text
A B, B C, C A >= N^(1/2)/C^2.                       (6.3)
```

Multiplying (6.3) shows

```text
A B C >= N^(3/4)/C^3,
```

so the common block satisfies

```text
G <= C^3 N^(1/4).                                   (6.4)
```

The canonical remaining local object is therefore:

> Three pairwise-coprime Gaussian blocks whose three pairwise norm products are
> all subpower-close to squares, and which satisfy the exact vector relation
> (5.4).

---

## 7. Relation to the Erdos--Rosenfeld/Chan method

Chan's perfect-square divisor argument succeeds because the Pythagorean
parameter `mu_i` and the gap `y_i-x_i` are bounded solely in terms of the fixed
interval constant.  Three divisors then yield simultaneous Pell equations with
bounded coefficients and bounded inhomogeneous terms.

In (6.2), the errors are subpower rather than bounded, which is still
potentially usable.  The decisive difference is that the Pell coefficients
`A,B,C` can have size `N^Theta(1)`.  Rewriting, for example,

```text
C X_AB^2-A X_BC^2 = A * error_BC-C * error_AB
```

leaves both coefficients and the right-hand side at polynomial scale.  The
available simultaneous-Pell estimates are then ineffective.

So the balanced-clique extraction reaches an almost-square system, but not the
bounded-coefficient system in Chan's theorem.

---

## 8. Concrete failure of the elementary continuation

The three-term equation (5.4) alone cannot yield a contradiction from pairwise
coprimality and small coefficients.  Integer and Gaussian equations

```text
A x+B y+C z=0
```

with pairwise-coprime, comparable large coefficients and small nonzero
`x,y,z` can occur at arbitrarily large scales.  The inequalities (6.2)--(6.4)
add strong compatibility, but no elementary gcd, triangle-inequality, or
Ptolemy argument forces one of `A,B,C` to be bounded.

The route therefore stops at a precise global arithmetic problem rather than at
the old pairwise `1/2` inequality:

```text
classify balanced pairwise-coprime Gaussian block triples
whose every pair product is subpower-close to a square.
```

A successful theorem here would combine naturally with the fixed-clique
extraction and could prove the endpoint conjecture.  At present, neither the
available simultaneous-Pell machinery nor standard `S`-unit estimates control
this large-coefficient regime.
