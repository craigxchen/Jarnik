# Independent audit of arbitrary rational CM contacts

The arithmetic theorem in
[cm_arbitrary_rational_contact_splitting.md](cm_arbitrary_rational_contact_splitting.md)
passes this independent proof audit. It extends the Boolean denominator
calculation in [rational_cm_inert_mass.md](rational_cm_inert_mass.md) to
arbitrary finitely many rational contacts. No estimate for p-adic
approximation to a noncyclic translate is required. The following
details make the compactness and normalization arguments explicit.

The existing exact checker `check_rational_cm_inert_mass.py` also
passed all nine cases (`m=5,13,65`). Those finite checks support the
cell identities; the infinite subsequence is established by the
proof below.

## CRT and simultaneous local bounds

The integers `c_f-1_(f=e)` are divisible by `m_e` for every `f`.
Consequently the displayed `S_e` are actual rational points, and
`Q_0-B_e=[m_e]S_e` holds in the abelian group `E(Q)`, including its
torsion subgroup. There is no division in `E(Q)` hidden in this step.
The identity `Q(n)-B_e=[m_e]R_e(n)` is exact.

At a fixed stage choose the finite rational prime set `Sigma` to
contain all model/CM-action exceptions, all primes dividing `2M`,
and all contact-section intersections. Omitting a rational prime
means omitting every prime above it in `Q(i)`; this keeps the
denominator decomposition invariant under conjugation. This finite
enlargement is harmless even when the initially prescribed
exceptional set is smaller.

Each map `z -> [alpha_S]([(M/m_e)]z+S_e)` has a nonzero isogeny as
its linear part. It can send `[n]P` to the identity for at most one
integer `n`: two such integers would give a nonzero multiple of a
nontorsion point in the finite kernel. Thus an integer `n_0` avoiding
all the finitely many zero conditions exists. A nonidentity algebraic
point remains nonidentity in every completion.

Here is the precise recurrence needed. In the compact group

```text
G = E(R) times product_(p in Sigma) E(Q_p),
```

there are arbitrarily large positive integers `q` for which the
diagonal point `[q]P` is arbitrarily near zero. For example, take a
convergent subsequence of `[n]P` with increasingly large gaps between
its indices, and subtract successive terms. Hence infinitely many
positive `n=n_0+q` return to any prescribed neighborhood of `[n_0]P`.
No density assertion in the entire product is needed.

The finitely many translated CM maps are continuous into the real,
complex, and appropriate finite extensions of these local fields.
Small neighborhoods of their nonidentity values have bounded
x-coordinates. Taking their simultaneous inverse images therefore
bounds every omitted local denominator term and both complex
contributions along the returning subsequence. The bounds are
stage-dependent constants, which is exactly what the proof needs.

## Cells, heights, and prime-power bounds

Away from `Sigma`, the reduction and formal-group argument applies
to any rational `R_e(n)`, not only to a multiple of one fixed point.
At a prime dividing the full denominator its reduction annihilator
contains `(m_e)`, so is a squarefree ideal divisor of `(m_e)`. The
denominator incidence is one upper Boolean interval. On that
interval, multiplication by the complementary CM factor has unit
differential and preserves the entire formal-parameter valuation.
Boolean inversion thus produces integral, pairwise coprime cells
with the full prime-power exponents. Rationality gives
`bar(C_T)=C_(bar(T))`; in particular inert primes occur only in
conjugation-stable cells.

Write `h_e=hat h(R_e(n))`. With the degree-one canonical-height
normalization used in both notes, bounded local terms give

```text
log Norm D_S = 2 Norm(alpha_S) h_e + O_stage(1),
log Norm C_T = 2 h_e product_(pi in T)(Norm(pi)-1) + O_stage(1).
```

The finitely many Boolean errors can be absorbed into the stage
constant. The full weight is `m_e^2`, the stable-cell weight is
`product_(p|m_e)(1+(p-1)^2)`, and the largest cell weight is
`product_(p|m_e)(p-1)^2`. Since
`hat h(Q(n)-B_e)=m_e^2 h_e`, these give exactly the factors `rho(m_e)`
and `2L(m_e)` in (3). An inert rational factor `p^a` has Gaussian
ideal norm `p^(2a)`, accounting for cancellation of the factor two
in the inert-mass estimate. For a split prime, either Gaussian
prime with its entire depth fits in one cell; the factor two in
the stated maximum-prime-power bound is a valid upper bound.
The omitted primes contribute only `O_stage(1)`.

The same height normalization gives
`log b(Q(n)-B_e)=hat h(Q(n)-B_e)+O_stage(1)`. Quadraticity gives
`hat h(Q(n))=M^2 n^2 hat h(P)+O_stage(n)` and changes this by only
`O_stage(n)` under any fixed contact translation.

## Diagonalization, allocation, and scope

At each stage the return times are unbounded. Choose one so large
that all stage constants and linear height errors are negligible
relative to `H=hat h(Q(n))`, and also require `H` to exceed the
previous stage's height. Taking `rho(m_e)` to zero proves all
three limits in (1), including the initially exceptional primes.
There is no missing uniformity requirement on the growing `m_e`.

Removing stable cells and `Sigma` costs `o(H)` rational logarithmic
mass. Choose one cell from each remaining conjugate pair. Their
norms multiply to the retained rational denominator, and their
largest logarithmic norm is `o(H)`. Cumulative allocation therefore
realizes any fixed positive proportions summing to one, with
`o(H)` errors. Coprimality within a contact follows from the cell
decomposition; coprimality across contacts follows because two
simultaneous reductions to their distinct sections are impossible
outside the fixed section-intersection primes.

The geometric consequence requires “degree-one boundary contact”
in the divisor sense: the pullback of each boundary divisor must
be one reduced rational point. Merely having a rational point of
higher intersection multiplicity would multiply its local and
global contact mass and would not give the stated normalization.
With reduced degree-one pullbacks and fixed integral models, the
local parameter comparison is valid outside finitely many primes,
as asserted in the source note.

This proves a denominator-and-norm allocation theorem. It neither
identifies the verified six-point curve as Gaussian CM nor provides
a common oriented Gaussian lift, the required private factors,
small primitive residues, or a short determinant-one frame. The
uniform centered-circle `C sqrt(R)` arc bound remains unproved.
