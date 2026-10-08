# Markov graph route: flat holonomy and entropy versus conductor switching

## Status

This note tests the proposal to replace angular ordering by a weighted graph on
all lattice points in an endpoint arc and to study a Markov chain on that graph.

Two conclusions are rigorous.

1. The two most natural cycle labels are flat coboundaries, so ordinary return
   probabilities do not create new arithmetic equations.
2. A stationary walk satisfies an exact information-theoretic inequality
   relating its entropy rate to the prime-by-prime switching probabilities.

The resulting inequality is a genuine stochastic formulation of the endpoint
problem, but by itself it does not prove a uniform bound because Shannon
entropy counts prime coordinates while conductor distance weights them by
`log p`.

---

## 1. Weighted graph

Let `V` be a finite set of Gaussian integers of common norm `N`, all lying in
one endpoint arc.  For `i != j`, let

```text
P_ij = product of squarefree conductor primes at which z_i,z_j
       choose opposite Gaussian orientations,
E_ij = log P_ij - (1/2) log N.
```

After removing the maximal common Gaussian conductor factor, write

```text
z_j-z_i = g_ij w_ij,
Norm(g_ij)=N/P_ij,
r_ij=Norm(w_ij) >= 1.
```

Then

```text
|z_j-z_i|^2 = N^(1/2) exp(-E_ij) r_ij.
```

A natural conductance is

```text
c_ij = N^(1/4)/|z_j-z_i|
     = exp(E_ij/2) r_ij^(-1/2).
```

Any symmetric nonnegative function of `(E_ij,r_ij)` gives a reversible chain
when rows are normalized.

---

## 2. Flat multiplicative holonomy

The norm-one quotient edge label is

```text
q_ij = z_j/z_i.
```

It is an exact multiplicative coboundary.  For every closed walk

```text
i_0 -> i_1 -> ... -> i_k=i_0,
```

one has identically

```text
product_a q_(i_a,i_(a+1)) = 1.
```

Prime-by-prime, every orientation switch occurs an even number of times around
the cycle.  This is forced by the vertex labels and contains no additional
arithmetic information.

Choosing Gaussian representatives `alpha_ij` with

```text
q_ij = alpha_ij/conjugate(alpha_ij)
```

only restates the same fact as

```text
product_a alpha_(i_a,i_(a+1)) is real.
```

Thus a high return probability does not by itself create a rare multiplicative
constraint.

---

## 3. Flat additive holonomy

The chord label

```text
c_ij = z_j-z_i
```

is an exact additive coboundary.  Every closed walk satisfies identically

```text
sum_a c_(i_a,i_(a+1)) = 0.
```

Writing `c_ij=g_ij w_ij` produces nontrivial-looking weighted relations among
residuals, but these are simply the same additive coboundary after pairwise
factor extraction.  One triangle gives

```text
g_ij w_ij + g_jk w_jk + g_ki w_ki = 0,
```

and overlapping cycles give Plucker/Ptolemy compatibility already analyzed in
other notes.

Therefore ordinary graph traces `tr(K^k)` count many identities that are
algebraically automatic.

---

## 4. Stationary prime-switch process

Now let `(X_n)` be any stationary Markov chain on `V`.  For each squarefree
conductor prime `p`, let

```text
epsilon_p(i) in {-1,+1}
```

be the Gaussian orientation of vertex `i` at `p`, and put

```text
q_p = Pr(epsilon_p(X_(n+1)) != epsilon_p(X_n)).
```

The expected conductor distance of one transition is exactly

```text
E[D(X_n,X_(n+1))] = sum_(p|N) q_p log p.
```

Hence its expected excess over the critical half-conductor level is

```text
E[E_(X_n,X_(n+1))]
  = sum_(p|N) (q_p-1/2) log p.                       (4.1)
```

This is the stochastic version of ordered transition excess.

---

## 5. Entropy-rate inequality

Let `h(u)=-u log u-(1-u)log(1-u)` be binary entropy in natural units.
Since the next vertex is determined by its complete prime-orientation vector,
subadditivity of conditional entropy gives

```text
H(X_(n+1) | X_n)
  <= sum_(p|N) H(epsilon_p(X_(n+1)) | X_n).
```

Conditioning only on the current orientation can increase entropy, so

```text
H(epsilon_p(X_(n+1)) | X_n)
  <= H(epsilon_p(X_(n+1)) | epsilon_p(X_n)).
```

For a stationary two-state process with average switch probability `q_p`,
concavity of binary entropy gives

```text
H(epsilon_p(X_(n+1)) | epsilon_p(X_n)) <= h(q_p).
```

Therefore every stationary walk satisfies

```text
H_rate := H(X_(n+1)|X_n)
  <= sum_(p|N) h(q_p).                                (5.1)
```

This is the exact Markov tradeoff:

- a walk with many possible next vertices requires large entropy rate;
- switching one prime much more or much less than half the time creates an
  entropy deficit at that coordinate.

Using

```text
h(1/2+x) <= log 2 - 2x^2,
```

one obtains

```text
sum_p (q_p-1/2)^2
  <= (r log 2-H_rate)/2,                              (5.2)
```

where `r` is the number of active squarefree prime coordinates.

By Cauchy--Schwarz,

```text
E[E_edge]
  <= sqrt(sum_p (log p)^2)
       sqrt((r log 2-H_rate)/2).                      (5.3)
```

Equation (5.3) is a rigorous entropy-versus-excess inequality.

---

## 6. Why this does not yet finish the endpoint theorem

The entropy deficit in (5.2) counts coordinates.  The desired excess in (4.1)
weights coordinates by `log p`.  A single large prime can contribute large
conductor excess while costing at most `log 2` entropy.

Thus the same information mismatch seen in the Shannon-entropy audit remains:

```text
probabilistic information per prime <= log 2,
arithmetic information per prime = log p.
```

In a strongly diffuse regime, where `sum (log p)^2` is small relative to
`(log N)^2`, (5.3) is useful.  In the unrestricted problem it cannot yield an
absolute bound without a separate arithmetic mechanism controlling large
prime weights.

---

## 7. Directed residual-edge lift

One can enlarge the state space from vertices to directed edges `(i,j)` and
let the next state be `(j,k)`.  The transition then sees the residual Gaussian
integers `w_ij,w_jk` and the exact triangle relation

```text
g_ij w_ij + g_jk w_jk + g_ki w_ki = 0.
```

This removes the trivial one-step coboundary description.  However, closed
walks in the line graph still assemble from overlapping triangle and Plucker
relations.  No uniform bound on the number of admissible continuations follows
from the relation alone; variable-support Gaussian three-term equations remain
flexible.

The edge lift is therefore a valid bookkeeping device, not yet a source of
spectral decay.

---

## 8. Honest conclusion

A Markov chain is useful for converting a dense short-edge graph into an
entropy-rate question.  The rigorous global invariant is

```text
H_rate <= sum_p h(q_p),
```

coupled with

```text
mean conductor excess = sum_p (q_p-1/2) log p.
```

But graph recurrence itself adds no arithmetic constraint because the natural
edge labels are flat coboundaries.  To finish the theorem via a stochastic
method, one needs a conductor-weighted entropy or a new arithmetic bound saying
that high-weight prime switches cannot be concentrated in a low-entropy part
of the transition process.
