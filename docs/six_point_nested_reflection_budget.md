# A sharp nested profile for one-coordinate reflection

The full fair threshold-cut profile does not force a positive conductor
cost for the critical replacement

```text
e_c -> e_* = e_a+e_b-e_c.
```

The general positive-support estimate gives only
`Delta log N >= W_{ab}-W_c=0` when every unordered cut has weight one. The
following nested profile attains equality for every number `m>=4` of source
coordinates. Thus zero is the exact optimum of the corresponding
piecewise-linear allocation problem.

## The sharp profile

Fix distinct labels `a,b,c` and let `Z` be the nonempty set of the other
`m-3` labels. Give every column below unit weight.

For `T=empty`, use the binary valuation vector

```text
e_a=e_b=1,       e_c=e_Z=0.
```

Its source width is one. The replacement has `e_c'=2`, so its width is two
and its cost is `+1`.

For every nonempty subset `T subset Z`, use the three-level vector

```text
e_c=e_(Z\T)=0,       e_a=e_b=1,       e_T=2.          (1)
```

Its two unit threshold gaps supply exactly the unordered cuts

```text
T,       {a,b} union T.                               (2)
```

If `T` is proper, a zero remains after deleting `c`, while the replacement
has value two. Both widths are two, so the cost is zero. If `T=Z`, deleting
`c` removes the sole zero level; the width falls from two to one, giving
cost `-1`.

The cuts in (2), together with the initial cut `{a,b}`, are all cuts whose
chosen side away from `c` contains either both or neither of `a,b`. Supply
every remaining cut by its binary valuation vector, oriented away from
`c`. Such a cut contains exactly one of `a,b`; the replacement gives `c`
the same binary value one, so its width cost is zero. Every unordered cut
has now received weight exactly one, and the total cost is

```text
+1-1=0.                                               (3)
```

Together with the general lower bound, (3) proves sharpness. Independent
binary primes alone have a positive cost; the three-level nested columns
are what reach the boundary.

## The replacement remains full fair

The same replacement permutes the leading cut profile. For nonempty proper
`T`, the two source cuts in (2) become

```text
{a,b,c} union T,       {c} union T,
```

whose unordered complements are `Z\T` and
`{a,b} union (Z\T)`. This is the profile pair indexed by `Z\T`. The
`T=Z` column becomes the single cut `{a,b}`, while the initial binary
column becomes the two cuts `{c}` and `Z`. On every remaining binary cut,
which has exactly one of `a,b`, the operation permutes that family by
toggling the side containing `c`. Hence all replacement cut weights are
also exactly one.

There is an asymptotic version using fixed split rational primes and
varying powers. Assign a different fixed prime `p=1 mod 4` to each column,
choose a positive integer exponent `k_p`, and multiply that column's levels
by `k_p`. Use these as valuations at one Gaussian prime above `p`, with complementary
valuations at its conjugate. Every coordinate then has the same norm, and
the minimum valuation at both conjugate primes is zero, so the tuple is
primitive.

Explicitly, if `h_p=max_k e_(p,k)` and `p=pi_p conjugate(pi_p)`, take

```text
z_k = product_p pi_p^(k_p e_(p,k))
                 conjugate(pi_p)^(k_p(h_p-e_(p,k))).
```

Then every `Norm(z_k)` is `product_p p^(k_p h_p)`. Primitive normalization changes
its exponent at `p` to the valuation width, so (3) is exactly the exponent
change in the common norm under `z_c -> z_a z_b/z_c`.

Choose the initial binary exponent so that
`k_empty log p_empty=w+O(1)`, the `T=Z` exponent so that

```text
k_Z log p_Z = w+sqrt(w)+O(1),
```

and every other exponent so that `k_p log p=w+O(1)`. Integer rounding
provides all these choices because the finitely many primes are fixed.
Both source and replacement profiles remain full fair to leading order,
while the exact squared-radius ratio is

```text
N_new/N_old = p_empty^k_empty/p_Z^k_Z
            = exp(-sqrt(w)+O(1)) -> 0.                    (4)
```

This is a subpower reduction and an exact nested valuation countermodel to
strict positivity. It is not a constructed endpoint configuration, has no
short-arc control, and does not make an iterable descent: applying the same
reflection again returns to the source.

The checker
[check_six_point_nested_reflection_budget.py](check_six_point_nested_reflection_budget.py)
enumerates all unordered cuts for `4<=m<=8`, verifies unit source and
replacement cut weights, checks the cut permutation, and confirms total
width change zero.
