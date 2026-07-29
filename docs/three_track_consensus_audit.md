# Three-track consensus audit

## Constraint

The current environment does not provide independent subagent processes. To preserve the user's requested decision discipline without misrepresenting the process, each proposed step is evaluated through three separately developed proof tracks:

1. direct finite sextic-basis construction;
2. product-section / higher-degree Subspace Theorem encoding;
3. direct quantitative Ru--Vojta exceptional-set analysis.

A step is retained only when at least two tracks support it after an explicit objection pass.

## Consensus 1: do not reconstruct the qualitative Ru--Vojta theorem

All three tracks reject rebuilding the entire adapted-basis argument merely to recover the known beta inequality. The Ru--Vojta theorem already packages the filtration and basis choices into the intrinsic beta constants.

The relevant fixed geometric data are:

- the blowup of `P^2` at `[1:1:1]`;
- the fixed boundary divisors;
- the fixed line bundle, at the minimal finite level represented by `6H-E`;
- a fixed approximation margin supplied by `6 - 3 beta = 4/9 < 1/2`.

Therefore the load-bearing question is not whether one can choose finitely many adapted bases. It is whether a quantitative version of the Ru--Vojta conclusion gives exceptional subvarieties whose number and degree are uniform when the finite place set varies but the geometric divisors remain fixed.

## Consensus 2: product sections are not currently load-bearing

The product-section identity

```text
max_B sum_{s in B} log |s(P)|_v = log max_B |prod_{s in B} s(P)|_v
```

is formally correct. It does not by itself produce a fixed finite family because the set of adapted bases can be infinite. No track currently has a proved finite domination lemma with uniform constants. The product-section route is therefore suspended unless the direct quantitative route fails.

## Consensus 3: the decisive parameter is the number of distinct forms/divisors, not the raw number of places

Quantitative Subspace Theorem bounds are often stated in terms of the number of distinct linear or higher-degree forms occurring across all places. In the common-conductor application, the same fixed boundary divisors are used at every place. If the Ru--Vojta reduction can be matched to a quantitative theorem without introducing place-dependent new forms, the relevant distinct-form count is fixed even though `|S|` varies.

This is the strongest remaining simplification and must be checked theorem-by-theorem, with exact quantifiers.

## Exact theorem to verify or prove

Let `X` be the fixed blowup and let `D_1,D_2,D_3` be the fixed properly intersecting divisors used in the generalized-GCD inequality. Fix a big line bundle `L` and `epsilon > 0` such that the beta coefficient gives a strict endpoint margin.

Prove that there exist constants `T`, `Delta`, and `H_0`, depending only on `(X,L,D_1,D_2,D_3,epsilon)`, such that for every finite set of places `S`, all points of height at least `H_0` violating the Ru--Vojta inequality lie in a union of at most `T` proper subvarieties of `X`, each of degree at most `Delta`.

The constants must not depend on `S`.

## Decision rule for the next step

The next step is accepted only if two tracks independently verify all of the following:

1. the forms/divisors entering the quantitative theorem are fixed as `S` varies;
2. the exceptional number and degree bounds depend on the number of distinct forms, not on `|S|` separately;
3. the height threshold is uniform in `S`;
4. pulling the exceptional subvarieties down to `P^2` preserves a uniform degree bound;
5. the numerical coefficient remains strictly below `1/2` after all error terms.

Until these five items are verified, the uniform bound is not proved.
