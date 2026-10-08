# Uniform endpoint problem: research reset

## Status

The uniform bound for lattice points on arcs of length `C * R^(1/2)` is not proved.

The Ru--Vojta/generalized-GCD track and the determinant equality-case track are no longer active proof strategies. They produced useful diagnostics, but every attempted completion either required a uniform exceptional-set theorem essentially as strong as the original problem or returned to the same endpoint equality with no strict gain.

The old files are retained only as an audit trail. They are no longer imported by `GaussianChain.lean` and should not be treated as the current roadmap.

## Global counterexample framing

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

The counterexample must be studied globally. Local determinant surplus and pairwise cut metrics discard too much of the Gaussian-angle structure.

## Binary threshold-layer expansion

Each exponent coordinate may be expanded into binary threshold layers

```text
sigma_(j,t)(a) = +1 if t <= a_j, and -1 otherwise,
1 <= t <= e_j.
```

This embeds the cluster into a binary orientation cube, with monotone chains inside each prime block. Every layer carries conductor weight `log p_j` and angle `arg(pi_j)`.

This representation is now the active state space for affine-cube and additive-energy arguments.

## Proved structured-regime theorem

The active detailed note is:

- `docs/affine_cube_walsh_rigidity.md`

It proves, modulo the standard imported form of Roth's theorem, that a sufficiently large endpoint cluster cannot contain a five-dimensional Boolean affine subcube on which every Gaussian layer restricts to an affine character.

The mechanism is exact Walsh isolation. On a `32`-point character cube, multiplying the points with the signs of one Walsh character isolates one Gaussian character block. Angular concentration forces that block into an `O(N^(-1/4))` neighborhood of one of the fixed `32` algebraic rays. Roth's theorem then forces every active character block to carry more than one fifth of the total conductor height. At most four blocks can be active, but four binary characters distinguish at most `16` points, contradicting the `32` distinct cube vertices.

Consequences:

1. no character-compatible Boolean affine cube of dimension at least five occurs for sufficiently large conductor;
2. after importing Balog--Szemeredi--Gowers, Freiman theory in `F_2^n`, and the fixed-density affine-subspace theorem, every hypothetical bad family must have

   ```text
   additiveEnergy(A) / |A|^3 -> 0;
   ```

3. the high-energy/small-doubling regime is therefore eliminated.

## Exact remaining regime

An unbounded counterexample must now be simultaneously:

- diffuse in conductor support;
- nearly balanced at almost every prime layer;
- free of five-dimensional affine subcubes;
- additive-Sidon-like, with almost every pair difference distinct;
- mapped by the Gaussian-prime angle form into an interval of width `O(N^(-1/4))`.

This sparse regime is not forced to contain affine cubes. Arbitrarily large abstract cube-free subsets of binary cubes exist, so no purely combinatorial cube-extraction statement can finish the proof.

The next theorem must be a sparse arithmetic inverse theorem:

> Bound an additive-Sidon-like family of Gaussian divisor exponent vectors whose pairwise ratios all lie within `O(N^(-1/4))` of `1`.

Possible inputs are:

- simultaneous Roth/Subspace-Theorem estimates for the many distinct exponent differences;
- a higher common-divisor theorem for their Gaussian numerators;
- sector packing for a sparse family of divisors of one Gaussian conductor, stronger than pairwise slope separation.

## Minimal-counterexample normalization

A proof by contradiction should choose a bad family minimizing, in order:

1. cluster cardinality above a fixed threshold;
2. total conductor height;
3. number of active split-prime coordinates;
4. total exponent mass.

This gives the following exact reduction rules:

- no coordinate is constant across the cluster;
- no nonempty oriented factor divides every point;
- no proper coordinate projection preserves an unbounded subcluster at endpoint scale;
- any exact multiplicative relation that reduces the active support contradicts minimality;
- any high-energy subset is excluded by Walsh rigidity.

These are the current replacements for determinant equality conditions.

## Live methods

The active tools are now:

- affine-cube/Walsh analysis in the high-energy regime;
- additive-energy decomposition into structured and Sidon-like cases;
- rational approximation to fixed algebraic rays after Walsh isolation;
- sparse simultaneous logarithmic-form estimates;
- Gaussian-divisor sector packing;
- numerical searches for minimal sparse high-concentration exponent codes.

The following are not active proof routes:

- uniform-in-`S` Ru--Vojta exceptional sets;
- further optimization of sextic adapted bases;
- universal or averaged determinant bases;
- balanced-cut surplus refinements;
- Segre re-embeddings intended only to alter the height coefficient.

## Repository policy

- `GaussianChain.lean` imports only the proved core development.
- Experimental uniform-bound modules and notes remain historical audits and are not part of the umbrella build.
- New active work must prove a statement used by the global inverse-concentration program.
- No new file should merely rename the missing theorem or transfer it to another exceptional-set or determinant lemma.
