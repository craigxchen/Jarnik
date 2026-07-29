# Uniform closed exceptional set versus the finite remainder

## Main distinction

Several formulations of the Ru--Vojta arithmetic general theorem state that,
for fixed geometric data, there is a proper Zariski-closed set `Z` independent
of the number field and of the finite place set `S`, such that the target
inequality holds for all but finitely many rational points outside `Z`.

This is stronger than an arbitrary `S`-dependent closed exceptional set, but it
is not yet enough for the common-conductor endpoint problem.

The quantifiers are

```text
exists Z, independent of S,
for every S,
there exists a finite set E_S,
such that the inequality holds outside Z union E_S.
```

The uniform endpoint argument needs either

```text
E_S is empty above a height threshold H0 independent of S,
```

or a uniform bound on the size or algebraic complexity of `E_S`.

Neither conclusion follows merely from the independence of `Z`.

## Why the finite remainder matters

For one Gaussian conductor `G`, the ordered divisor-pair grid is finite.  An
`S(G)`-dependent finite exceptional set could contain the entire grid.  Thus
"all but finitely many points" gives no information at one fixed conductor,
even if every pair has arbitrarily large conductor height as `G` varies.

To use the qualitative uniform-`Z` theorem, one would need a sequence argument
showing that infinitely many bad conductor grids force a single infinite
sequence outside `Z` to which one fixed `S` applies.  That is unavailable:
`S(G)` varies with `G`, and passing to the union of all supports usually gives
an infinite set of places, outside the theorem's hypotheses.

## Quantitative higher-degree theorem

The older quantitative higher-degree Subspace Theorem gives an explicit height
threshold and explicit bounds for exceptional subvarieties.  However its stated
bounds include the cardinality `s = |S|` in the exponent of the component-count
bound.  Therefore it does not immediately combine with the uniform-`Z`
qualitative refinement.

The desired missing theorem would have the mixed strength:

1. a closed exceptional set or bounded-degree exceptional family independent of
   `S`; and
2. a height threshold independent of `S`, eliminating the finite remainder for
   all sufficiently large common-conductor points.

## Current conclusion

The statement

```text
Z can be chosen independently of S
```

is genuine and useful, but it does not prove the uniform endpoint bound.
The exact remaining issue is no longer the closed exceptional set itself; it is
uniform elimination or uniform control of the finite `S`-dependent remainder.
