# Fixed-weight reformulation audit

## Purpose

The previous one-scale audit treated the twisted weights `c` as pair-dependent.
This note isolates the more basic question: is that dependence intrinsic, or is
it only the local maximizing basis that depends on the point?

## What the absolute parametric theorem actually fixes

Evertse--Ferretti fix:

- a tuple of local systems of linearly independent forms `L_i^(v)`;
- one normalized weight tuple `c_{iv}`;
- the parameter `Q`.

The point `x` enters only through the local maximum

```text
max_i ||L_i^(v)(x)||_v Q^(-c_{iv}).
```

Thus the maximizing index may depend on `x` without changing the twisted-height
system. Point-dependent maxima are built into the definition.

## Consequence for the Ru--Vojta matching

The load-bearing question is therefore not automatically

```text
for every point P choose a new weight tuple c(P).
```

It is instead:

```text
Can the finite max-over-adapted-bases inequality be written using one fixed
finite collection of local forms and one conductor-dependent but
point-independent normalized weight tuple c_G?
```

If yes, then all ordered divisor pairs from one conductor are solutions of one
fixed twisted-height system at one fixed scale `Q_G`. The quantitative
parametric theorem then gives a uniformly bounded number of exceptional
subspaces for the whole pair set; no union over pair-dependent weight systems is
needed.

## Why the previous fixed-Q argument did not settle this

Proposition 4.2 of Evertse--Ferretti constructs, for fixed `L,c` and a fixed
parameter interval, a proper space by spanning all small twisted-height points.
That space may vary with `c`. This is harmless if the common-conductor
construction supplies one `c_G` for the entire pair grid, but fatal if `c` is
chosen separately for each pair.

## Exact next lemma

For the sextic space `H^0(P^2,I_[1:1:1](6))`, construct finite local form lists
`L_i^(v)` and normalized conductor weights `c_{iv}(G)` such that for every
ordered divisor pair `P=(u_A,u_B)` from the same Gaussian conductor,

```text
Ru--Vojta local max at P
  <= max_i ||L_i^(v)(P)||_v Q_G^(-c_{iv}(G))
```

at every place, with a total loss smaller than the available coefficient margin

```text
1/2 - 4/9 = 1/18.
```

The tuple `c(G)` may depend on the conductor, but it must not depend on `A,B`.

This is now the shortest plausible completion route. It supersedes attempts to
control a union over pair-dependent exceptional subspaces before first proving
that such pair-dependence is genuinely necessary.
