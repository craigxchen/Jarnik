# The inert radical forced by an integer cotangent clique

This note extracts an aggregate radius-free consequence
from the local valuation bounds. It is a necessary condition on `L`, but it
does not improve the established `(1+epsilon) log R/log log R` growth rate.

Let `A={X_1,...,X_m}` be distinct integers and `L>0`, with

```text
(X_i X_j+L^2)/(X_i-X_j) in Z                 (1)
```

for every pair. For an odd inert prime `p=3 mod 4`, write `e=v_p(L)`. The
local capacity theorem gives

```text
m <= (p+1)p^e-1.                              (2)
```

Consequently every inert prime below `m` is forced into `L`, and its exact
required exponent is bounded below:

```text
v_p(L) >= max(0, ceil(log_p((m+1)/(p+1)))).     (3)
```

Equivalently, with `I_m={p<m:p=3 mod 4}`,

```text
L_m := product_(p in I_m) p^ceil(log_p((m+1)/(p+1))) divides L.  (4)
```

The divisibility in (4) is an ordinary integer divisibility statement; no
pairwise coprimality or squarefree assumption on `L` is being made.

## Proof and the occupancy content

For completeness, partition `A` by `s=v_p(X)`. If `s<e`, write
`X=p^s u` with `u` a unit. Since `p` is inert, valuation of `X^2+L^2` is
`2s`; (1) therefore implies

```text
v_p(X-Y)=s+v_p(u-v)<=2s,
```

so the normalized residues `u` are distinct modulo `p^(s+1)`. This stratum
has at most `(p-1)p^s` members. The stratum `v_p(X)>=e`, after writing
`X=p^e u`, has at most `p^(e+1)` members. Summing gives (2), exactly as in
[the local capacity note](integer_cotangent_local_bounds.md).

If `p<m` and `e=0`, (2) says `m<=p`, a contradiction. Thus every such `p`
divides `L`. Keeping the exponent in (2) gives (3), and multiplying these
independent prime-power requirements proves (4).

The dyadic capacity also gives

```text
v_2(L) >= max(0, ceil(log_2((m+1)/4))).           (5)
```

This is weaker than the inert radical for large `m`, but records the full
local information at the ramified prime.

## Size and why it does not close the height target

The radical part alone yields

```text
log L >= sum_(p<m, p=3 mod 4) log p.             (6)
```

The classical prime number theorem in arithmetic progressions gives
the asymptotic
`log L >= (1/2+o(1))m`. The extra exponents in (3) contribute only lower-order
weight relative to this primorial term.

This is not a radius lower bound. In the cotangent normalization, clearing
all anchor and pair imaginary residues only gives

```text
log L <= binom(m+1,2) log T,                       (7)
```

when the relevant residues have absolute value at most `T`. Combining (6)
with (7) gives at most a relation between `m` and the residue height; it does
not force the all-edge lcm `N` to dominate `(min X_i/L)^4`. In particular it
does not improve the current radius-dependent point bound.

The local conditions are also not a hidden factorial bound. At a fixed prime,
the normalized residue classes may occupy all allowed classes in each
valuation stratum, and local assignments at different primes can be combined
by the Chinese remainder theorem. This only constructs compatible local
residue data; it does not assert that those data form an actual global
cotangent clique. The actual pair equation can still rule out many such
assignments, but the local inert restrictions alone
provide no additional cross-prime loss beyond (4). The ordered Gaussian
family in [the joint inert sieve note](ordered_inert_sieve_route.md) further
shows that simultaneous prime-power occupancy constraints can be realized
with balanced collisions to leading order; angular order does not change the
leading `m^2 log m` cost.

## Bounded exact check

The companion script
[check_integer_cotangent_inert_radical.py](check_integer_cotangent_inert_radical.py)
enumerates small `L` and divisor-generated cliques, verifies (1), and checks
(3)--(5) for every clique found. It also prints the largest clique seen for
each `L`; this is a finite sanity check and is not used as evidence for an
unbounded theorem.
