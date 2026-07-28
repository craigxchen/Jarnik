# A uniform generalized-GCD route to the square-root arc problem

## Status

This note separates three logically distinct statements.

1. **Circle-side reduction.**  A uniform bounded-degree exceptional-curve theorem for the two-variable generalized-GCD inequality implies a uniform bound for lattice points on arcs of length `C * sqrt R`.  This implication is proved below.
2. **Finite-level Ru--Vojta calculation.**  On the blowup of `P^2` at `[1:1:1]`, the finite section counts can be evaluated exactly.  In particular, the asymptotic beta threshold used in Yasufuku's argument is attained at an explicit finite level.
3. **Remaining gap.**  The cited quantitative Subspace Theorem for the max-over-local-bases system has an explicit dependence on the number of places `|S|`.  Therefore the published results do **not yet** give a total exceptional-curve degree bounded uniformly as `S` varies.

Accordingly, this document is a rigorous **conditional proof and gap audit**, not a claim that the endpoint conjecture has already been proved.

---

## 1. The uniform exceptional-curve hypothesis

Let `K = Q(i)`.  Fix a real number `eta < 1/2`.

### Hypothesis UEC(eta)

There exist constants `D, H0` depending only on `K` and `eta` such that for every finite set of places

```text
S superset M_K^infty
```

there is a nonzero polynomial

```text
F_S(X,Y) in K[X,Y],       deg F_S <= D,
```

with the following property.  For all `S`-units `x,y in K^*` with

```text
h([1:x:y]) > H0
```

and `F_S(x,y) != 0`,

```text
log GCD^+(x-1,y-1) < eta * h([1:x:y]).                 (UEC)
```

The polynomial may depend on `S`; its total degree may not.

Yasufuku proves a two-variable generalized-GCD inequality with coefficient at most `0.45` outside a proper hypersurface.  What is not supplied in the published statement is a bound on the **total degree** of the exceptional hypersurface independent of `S`.

---

## 2. Conditional uniform bound for square-root arcs

### Theorem 2.1

Assume `UEC(eta)` for some `eta < 1/2`.  Then for every `C > 0` there is a constant `B(C)` such that every arc of length at most

```text
C * sqrt R
```

on a centered circle of radius `R`, with `R^2 in N`, contains at most `B(C)` points of `Z[i]`.

### Proof

Let

```text
z_0, z_1, ..., z_{M-1} in Z[i]
```

be distinct, have modulus `R`, and lie on one such arc.  Put

```text
W   = log R,
u_j = z_j / z_0        (1 <= j < M).
```

#### Step 1: one common `S`

Let `S` consist of the archimedean place of `Q(i)` and all finite places lying over rational primes dividing `R^2`.  Since every `z_j` has norm `R^2`, every quotient `u_j` is an `S`-unit.  The same `S` works for the entire cluster.

#### Step 2: archimedean proximity

Chord length is at most arc length, hence

```text
|u_j - 1| = |z_j-z_0| / R <= C * R^(-1/2).            (2.1)
```

For distinct `i,j`, the complex-place contribution to the exceptional-divisor local height at `[1:1:1]` is therefore at least

```text
W/2 - log C.                                           (2.2)
```

At every nonarchimedean place,

```text
|x-1|_v <= max(1,|x|_v),
```

so the corresponding local contribution to `GCD^+` is nonnegative.  Consequently

```text
log GCD^+(u_i-1,u_j-1) >= W/2 - log C.                 (2.3)
```

#### Step 3: joint-height upper bound

We claim

```text
h([1:u_i:u_j]) <= W.                                   (2.4)
```

At a split rational prime `p = pi * conjugate(pi)`, write

```text
E = v_p(R^2),
a_k = v_pi(z_k),
v_conj(pi)(z_k) = E-a_k.
```

The combined normalized contribution of the two conjugate places to
`h([1:u_i:u_j])` is

```text
(1/2) * (max{a_0,a_i,a_j} - min{a_0,a_i,a_j}) * log p,
```

which is at most `(E/2) log p`.  Inert and common ramified valuations cancel in the ratios, and the complex place contributes zero because `|u_i|=|u_j|=1`.  Summing gives

```text
h([1:u_i:u_j]) <= (1/2) sum_p E_p log p
                 = (1/2) log(R^2)
                 = W.
```

#### Step 4: every off-diagonal pair is exceptional

Choose `R` sufficiently large that

```text
log C < (1/2-eta) W                                   (2.5)
```

and the joint height is above `H0` whenever needed.  If
`F_S(u_i,u_j) != 0`, then `UEC(eta)`, (2.3), and (2.4) imply

```text
W/2 - log C < eta * W,
```

contradicting (2.5).  Hence

```text
F_S(u_i,u_j) = 0       for every i != j.               (2.6)
```

The same polynomial works for every pair because `S` is common to the circle.

#### Step 5: grid-zero bound

Let

```text
A = {u_1,...,u_{M-1}}.
```

The elements are distinct.  A nonzero bivariate polynomial of total degree at most `D` has at most `D * |A|` zeros on `A x A`: view it as a polynomial in the second variable and account for the at most `D` first-coordinate values at which the specialization vanishes identically.

Equation (2.6) supplies every ordered off-diagonal pair, so

```text
|A| (|A|-1) <= D |A|.
```

Thus `|A| <= D+1` and

```text
M <= D+2.
```

The finitely many radii not covered by the large-height argument can be absorbed into `B(C)`.  This proves the theorem.  `square`

---

## 3. Exact finite-level beta calculation on `Bl_P P^2`

Let

```text
V = Bl_P P^2,       P = [1:1:1],
```

and write `H` for the pullback of a line and `E` for the exceptional divisor.  Fix integers `N >= 1` and `k >= 0`; the plane degree is

```text
d = N+k.
```

The space of degree-`N+j` plane forms vanishing to order at least `N` at `P` has dimension

```text
h_j = binom(N+j+2,2) - binom(N+1,2).                   (3.1)
```

The finite numerator in the beta approximation is

```text
sum_{j=0}^{k-1} h_j.
```

A direct summation gives

```text
sum_{j=0}^{k-1} h_j
  = k(k+1)(3N+k+2)/6.                                  (3.2)
```

The finite beta ratio is therefore

```text
beta_{N,k}
  = [k(k+1)(3N+k+2)/6]
    / [N * (binom(N+k+2,2)-binom(N+1,2))].             (3.3)
```

### Explicit level

Take

```text
N = 100,
k = 247,
d = 347,
c = d/N = 347/100.
```

Then

```text
h^0(V,347H-100E) = 55,676,
finite numerator = 5,604,924,
beta_{100,247}   = 45,201 / 44,900 > 1.                (3.4)
```

These identities are formalized in `GaussianChain/FiniteBeta.lean`.

This calculation shows that the finite approximation level in the Ru--Vojta/Autissier filtration argument can be chosen explicitly; the ambient projective dimension need not be left implicit behind asymptotic Riemann--Roch.

---

## 4. What fixed finite level does and does not prove

Fixing `(N,k)` has two important consequences.

1. The vector space of sections is finite-dimensional with an explicit dimension.
2. The direct Ru--Vojta proof uses only finitely many adapted bases once the discretization parameter is fixed.  Hence only a fixed finite collection of algebraic linear forms occurs.

This resolves Yasufuku's stated problem of not knowing the finite approximation dimension.

It does **not**, by itself, prove `UEC(eta)` uniformly in `S`.

The relevant published quantitative theorem for the generalized max-over-local-bases inequality has a bound on the number of exceptional subspaces containing a factor of the form

```text
(constant)^(dimension * |S|).
```

Thus its total exceptional degree is not uniform as the number of places varies.  The fixed-form parametric theorem has stronger uniformity, but the direct Ru--Vojta inequality takes a local maximum over adapted bases, and passing to that generalized system is precisely where the published quantitative count reacquires `|S|`.

Therefore the logical state is:

```text
explicit finite beta level                         proved;
finite adapted-basis family                        proved;
conditional circle endgame                         proved/formalized;
uniform exceptional degree independent of |S|      open in this argument.
```

---

## 5. The precise remaining theorem

The remaining research target is not another circle estimate.  It is the following special quantitative statement.

### Fixed-family, variable-place exceptional theorem

For the fixed finite family of linear forms produced by the explicit Ru--Vojta filtration on `Bl_P P^2`, prove that the solutions to the corresponding max-over-bases Subspace-Theorem inequality lie in a number of proper subspaces bounded independently of the finite set of active places.

Any proof must exploit additional structure absent from the general Corollary 3.2 estimate, for example:

- the forms arise from one fixed filtration geometry;
- the local weight vectors lie in finitely many rational polyhedral chambers;
- the points under consideration are `S`-units and the coordinate divisors contribute no outside-`S` term;
- only the distinguished Harder--Narasimhan exceptional space may be needed, rather than the full quantitative covering by arbitrary subspaces.

A uniform theorem of this form would establish `UEC(eta)` and, by Theorem 2.1, the square-root arc conjecture.

---

## 6. Lean formalization boundary

The repository contains two new modules.

### `GaussianChain/FiniteBeta.lean`

Formalizes:

- the Hilbert section-count polynomial;
- the closed finite numerator formula;
- the finite beta ratio;
- the explicit values at `(N,k)=(100,247)`;
- the strict inequality `1 < beta_{100,247}`.

### `GaussianChain/UniformExceptionalReduction.lean`

Formalizes, as a theorem with explicit hypotheses:

- the numerical `1/2-eta` contradiction;
- the conclusion that every off-diagonal pair is exceptional;
- an abstract grid/fiber bound;
- the resulting cluster-cardinality estimate `M <= D+2`.

The deep algebraic-geometry and Diophantine approximation input is represented as a hypothesis.  It is not available in mathlib, and it has not been introduced as an axiom.

---

## References

- K. Yasufuku, *GCD inequalities inspired by Vojta's conjecture*, Monatshefte fuer Mathematik (2025).
- M. Ru and P. Vojta, *A birational Nevanlinna constant and its consequences*.
- J.-H. Evertse and R. Ferretti, *A further improvement of the quantitative Subspace Theorem* / quantitative absolute and parametric Subspace Theorem formulations.
- C. Chen, *A Sub-Logarithmic Bound for Lattice Points on Small Circular Arcs*.
