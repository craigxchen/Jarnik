# Local reduction of central moments at an unbalanced cut

This note records the exact Vandermonde statement at one good split prime.
It prevents an incorrect inference from a rectangular moment matrix: maximal
rank leaves a one-dimensional cofactor kernel and does not force any one
coefficient times a maximal minor to vanish.

Let `n=2h`, and let the ordinary primitive central vector `u` satisfy

```text
sum_i u_i r_i^j=0,          -h+1<=j<=h-1.            (1)
```

Fix a split rational prime `p` not dividing `2L`, and one Gaussian prime
`pi` above it. Suppose the entire allocation profile at this prime has
exactly two levels, separated by `a>0`, with minority size `s<h`.
Choose the orientation and an anchor
in the majority so that

```text
v_pi(r_i)=a  on the minority,
v_pi(r_i)=0  on the majority.                        (2)
```

The exact central-coefficient valuation formula gives

```text
v_p(u_i)=(h-s)a  on the minority,
v_p(u_i)=0       on the majority.                    (3)
```

## The surviving majority system

For a minority index, the valuation of its term in (1) is

```text
v_pi(u_i r_i^j)=(h-s+j)a.
```

It is therefore positive for

```text
s-h+1<=j<=h-1.                                      (4)
```

There are exactly `n-s-1` exponents in (4).  Reduction modulo `pi`
deletes the minority terms and gives that many consecutive moment equations
on the `n-s` majority nodes.

The majority nodes are units and have pairwise distinct reductions.  Indeed,
the good-prime difference estimate gives

```text
v_pi(r_i-r_j)=min(v_pi(r_i),v_pi(r_j))=0
```

for two majority indices.  Multiplying each column by a suitable unit turns
the reduced matrix into the standard `(n-s-1)` by `(n-s)` Vandermonde matrix.
It has rank `n-s-1`, hence a one-dimensional kernel.  Its signed maximal
minors are all nonzero; explicitly, a kernel vector is proportional to

```text
1/product_(j!=i)(r_i-r_j)
```

times a harmless unit power of `r_i`.  Thus every kernel entry is a residue
unit.  This is fully compatible with all majority coefficients in (3) being
units.  It creates no additional divisibility by `p` and no lower bound for
the cut height `a`.

## Why the remaining moments cannot be added

At `j=s-h`, the minority and majority terms both have valuation zero, so the
equation couples the two blocks and neither block can be deleted.  For

```text
-h+1<=j<=s-h-1,
```

write `r_i=pi^a xi_i` and `u_i=pi^((h-s)a) eta_i` on the minority.  After
scaling the equation by its lowest valuation and reducing modulo `pi`, the
majority terms disappear and one obtains `s-1` consecutive equations on the
`s` minority unit nodes `xi_i`.  The same good-prime difference estimate
shows that their reductions are distinct, so this is another rank `s-1`
Vandermonde system with a one-dimensional unit cofactor kernel. The boundary
equation `j=s-h` then couples the two one-dimensional kernels nontrivially.
For the majority let `m=n-s` and `j_0=s-h+1`. Its kernel has entries
`c r_i^(-j_0)/product_(ell!=i)(r_i-r_ell)`. At exponent `j_0-1`, its
moment is `c (-1)^(m-1)/product_i r_i`, a nonzero residue. On the
minority, the boundary is the first exponent after its `s-1` consecutive
vanishing moments, so its cofactor moment is likewise a nonzero multiple
of its kernel parameter. The boundary therefore imposes exactly one
relation between the two parameters, with unit coefficients.

Consequently all `n-1` central moment equations have exactly the expected
local rank behavior.  The high equations alone do not yield a vanishing
maximal minor or coefficient growth.  Any growth conclusion must use the
already known coefficient valuations in (3), interaction among several
prime layers, or additional global arithmetic.

For example, over `F_7` the nodes `1,2,3,4` and weights
`(-1,3,-3,1)` satisfy the three moments of degrees `0,1,2`.  The `3` by `4`
Vandermonde matrix has rank three and nonzero maximal minors; its nonzero
cofactor vector is precisely what rank three predicts.
