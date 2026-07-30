# Uniform endpoint problem: research reset

## Status

The uniform bound for lattice points on arcs of length `C * R^(1/2)` is not proved.

The Ru--Vojta/generalized-GCD track and the determinant equality-case track are no longer active proof strategies. They produced useful diagnostics, but every attempted completion either required a uniform exceptional-set theorem essentially as strong as the original problem or returned to the same endpoint equality with no strict gain.

The old files are retained only as an audit trail. They are no longer imported by `GaussianChain.lean` and should not be treated as the current roadmap.

## What happened to the counterexample framing

The intended contradiction setup was sound, but it was absorbed into a local determinant-surplus analysis.

Assume there is a sequence of endpoint clusters with cardinality tending to infinity. After removing common Gaussian factors, write the split part of the common norm as

```text
N = product_j p_j^(e_j),
```

choose one Gaussian prime `pi_j` above each `p_j`, and represent every lattice point by an exponent vector

```text
a = (a_j),  0 <= a_j <= e_j.
```

Up to a fixed unit, its angular coordinate is

```text
Phi(a) = sum_j (2 a_j - e_j) theta_j  mod 2 pi,

theta_j = arg(pi_j).
```

An endpoint cluster is therefore a set `A` of exponent vectors for which all values `Phi(a)` lie in an interval of length

```text
O(N^(-1/4)).
```

The previous work replaced this global concentration statement by equality conditions in a Ramana determinant. That was too lossy: abstract balanced cut systems can imitate the determinant equality case without respecting the arithmetic of Gaussian prime angles.

The counterexample framing should instead remain global from the start.

## New active program: inverse concentration plus arithmetic rigidity

The proposed route has two stages.

### Stage 1: inverse concentration

Prove that if a large finite set of box points

```text
A subset product_j [0,e_j]
```

is mapped by `Phi` into one interval of width `N^(-1/4)`, then the active angle multiset has low-complexity additive structure.

The desired conclusion is not merely that many threshold cuts are balanced. It should say that, after discarding a bounded number of coordinates or a bounded fraction of conductor height, the angles `theta_j` lie near a generalized arithmetic progression of bounded rank, or equivalently that the difference set `A-A` is controlled by a bounded-rank approximate kernel of `Phi`.

This is an inverse Littlewood--Offord/Freiman-type statement adapted to:

- a non-random subset of an integer box;
- a circular target interval of exponentially small width;
- weighted coordinates with multiplicities `e_j`;
- a conclusion measured by conductor height `sum e_j log p_j`, not by the raw number of coordinates.

### Stage 2: Gaussian arithmetic rigidity

Exploit that

```text
exp(2 i theta_j) = pi_j / conjugate(pi_j)
```

is a norm-one element of `Q(i)` attached to a Gaussian prime.

A bounded-rank additive model for many `theta_j` should force one of the following:

1. a large block of conductor height is carried by a bounded set of Gaussian primes, allowing direct descent;
2. many Gaussian primes lie in a bounded collection of exceptionally narrow angular progressions or sectors;
3. there is an exact multiplicative relation among distinct `pi_j / conjugate(pi_j)`.

The third alternative is impossible by unique factorization unless the relation is trivial. The first is already compatible with existing heavy-block descent. The second becomes an analytic number-theory problem about adversarial Gaussian primes in narrow sectors, rather than a moving exceptional-set problem.

## Minimal-counterexample normalization

A future proof should begin with a counterexample minimizing, in order:

1. cluster cardinality above a fixed threshold;
2. total conductor height;
3. number of active split-prime coordinates;
4. total exponent mass.

This normalization gives immediate reduction rules:

- no coordinate is constant across the cluster;
- no nonempty prime block divides every point in the same orientation;
- no proper coordinate projection preserves the full cluster;
- any exact multiplicative relation among active norm-one prime ratios reduces support and contradicts minimality.

These are the replacement for determinant equality conditions.

## First concrete theorem to attack

The first target is the following weighted inverse-concentration statement.

> **Gaussian box inverse theorem.** For every sufficiently large fixed `M`, there exist constants `r(M)` and `eta(M)>0` such that whenever `M` exponent vectors from a Gaussian conductor box have angular diameter at most `N^(-1/4)`, either:
>
> - some set of at most `r(M)` prime coordinates carries at least `eta(M)` of the conductor height; or
> - the corresponding Gaussian prime angles admit a bounded-complexity approximate additive relation strong enough, after algebraic separation, to become an exact multiplicative relation.

The theorem is deliberately stated as a dichotomy with a heavy-coordinate outcome. Pure inverse Littlewood--Offord statements based only on the number of coordinates are not adequate in the diffuse weighted setting.

## Methods worth testing

The live tools are:

- weighted inverse Littlewood--Offord theory;
- Balog--Szemeredi--Gowers and Freiman compression on selected exponent differences;
- concentration inequalities for product measures on conductor boxes;
- lower bounds for linear forms in logarithms after rank reduction;
- Gaussian-prime counting in narrow sectors for the remaining structured case;
- numerical searches for minimal high-concentration exponent boxes, used to identify the correct inverse statement.

The following are not active proof routes:

- uniform-in-`S` Ru--Vojta exceptional sets;
- further optimization of sextic adapted bases;
- universal or averaged determinant bases;
- balanced-cut surplus refinements;
- Segre re-embeddings intended only to alter the height coefficient.

## Repository policy after the reset

- `GaussianChain.lean` imports only the proved core development.
- Experimental uniform-bound modules and notes remain as historical audits but are not part of the umbrella build.
- New work should enter only after it proves a statement used by the inverse-concentration program.
- No new file should merely rename the missing theorem or transfer it to another exceptional-set or determinant lemma.
