# A balanced nonlinear Walsh assignment with no low-height coset character

This finite audit records a limitation of the affine-coset phase method for
arbitrary balanced flip assignments. It neither constructs endpoint arcs nor
gives an asymptotic obstruction to them. The assignment is on `M=32` rows,
uses every nonzero Walsh label once and label `3` twice, and admits distinct
physical flipped columns for any `b>=2` copies per label:

```text
[28,21,11,30,24,23,13,3,12,7,25,4,31,1,22,3,
 19,18,10,9,2,16,17,27,26,6,15,5,14,8,20,29].
```

The exact coset calculation extends the fully aligned calculation in
[the linear-support character note](linear_support_character_height_obstruction.md).
Let `P=x0+V`, `h=|V|>=2`, and choose a nonzero restriction `l in V*`.
Put `c_(x0+v)=(-1)^l(v)` and `c=0` off `P`; its sum is zero. Let `k`
be the number of rows `x in P` for which the assigned label `a(x)`
restricts to `l`. The unflipped coefficient is

```text
u_a=(1/2) sum_x c_x (-1)^(a dot x)
   =(h/2)(-1)^(a dot x0) if a|V=l, else 0.
```

Exactly `M/h` nonzero Walsh classes meet `a|V=l`, so their `b`
physical copies have total absolute coefficient sum `bM/2`. Each
matched row flips a distinct copy, lowers its absolute coefficient by
one, and each unmatched row raises an initially zero coefficient by
one. Thus the **exact physical-column formula**, without assuming full
alignment, is

```text
L(P,l)=sum_(physical j) |v_j|=bM/2+h-2k.               (1)
```

At equal prime weights the height lies below `W0/2` only if
`2k>h+b/2`. For `b=5`, the exhaustive counts for this nonlinear map
are:

| `dim V` | `h` | directions | maximum `k` | minimum `L(P,l)` | `L-5(M-1)/2` |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 31 | 2 | 78 | 0.5 |
| 2 | 4 | 155 | 3 | 78 | 0.5 |
| 3 | 8 | 155 | 4 | 80 | 2.5 |
| 4 | 16 | 31 | 3 | 90 | 12.5 |
| 5 | 32 | 1 | 2 | 108 | 30.5 |

There are 155 two-dimensional directions, eight cosets per direction,
and three nonzero restrictions per coset. Among their 3720 triples,
exactly 128 have `k=3` and none have `k=4`: this assignment has no
aligned affine plane. The checker enumerates subspaces, cosets, and
restrictions rather than relying on a search heuristic. It also builds
the physical coefficients of one `k=3` plane directly.

The obstruction survives the near-equal **actual prime logs** used in
the prior notes for primitive source tuples. If `r=b(M-1)` and
`ell<=log p_j<ell+1/r`, then
`V(Beta)>=L(P,l) ell` and `W0<r ell+1`. Hence

```text
V(Beta)-W0/2 > (b/2+h-2k) ell-1/2.                  (2)
```

For this map with `b=5` and `ell>=log 5`, the right side is at least
`(log 5-1)/2>0` for every affine-coset character. Increasing `b`
only increases the positive gap. A large common Gaussian factor adds
to the source radius but cancels from the character; equation (2) is
therefore a statement about the primitive prime weight `W0`, and it
does not exclude a certificate that uses such common content.

The four-row limitation is sharper than a failed plane search. Any
zero-sum vector supported on four distinct rows with coefficients
`+1,+1,-1,-1` either has affine-plane support, in which case it is the
nonzero character of that plane and obeys (1), or has nonplanar
support. In the latter case translate one row to zero. The other three
row differences are independent, so their Walsh signs run uniformly
through eight patterns. For each of the six balanced choices of row
signs, the sum of the absolute half-coefficients over these eight
patterns is six. The constant Walsh coefficient is zero because the
row signs sum to zero. This gives old absolute coefficient sum
`3bM/4`; four distinct flips can lower it by
at most four. At `M=32,b=5` this is at least `116`, far above
`W0/(2ell)=77.5` at equal weights. Thus this balanced nonlinear map
blocks every four-row ±1 route to a character below half the primitive
weight, not merely the four-row aligned-coset construction.

The elementary rational-angle phase gap in the prior notes gives a
*lower* bound on any endpoint character height. Heights above `W0/2`
do not contradict that bound. This finite example therefore identifies
the precise combinatorial hurdle for extending the affine certificate:
one must find a different short-row character, exploit common content
or actual phase relations, or prove some large partial alignment for
all balanced nonlinear maps. It does not settle the uniform endpoint
count.

The [checker](check_walsh_nonlinear_coset_height_obstruction.py) verifies
the balance, all exact counts, all minimum coefficient sums, and a
physical-column fixture.
