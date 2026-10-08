# A universal large isotropic subspace for binary bilinear forms

Let `B` be any bilinear form on an `n`-dimensional vector space over
`F_2`; `B` need not be symmetric or nonsingular. For `n>=2`, there is a
subspace `U` of dimension at least

```text
d=max(0,floor((n-2)/2))                                  (1)
```

such that `B(u,v)=0` for every `u,v in U`. The one-dimensional case is
trivial with `d=0`. For the Walsh affine assignment
`a(x)=Lx+b`, the relevant form is `B(u,v)=u dot Lv` (or its transpose),
so this improves the isotropic dimension `floor(t/3)` in
[the affine-assignment note](walsh_affine_assignment_plane_certificate.md)
to `floor((t-2)/2)` for sufficiently large `t`. The improvement uses the
binary field; the symmetric Witt-index theorem alone does not prove it
for this nonsymmetric `B`.

## Nondegenerate polar form

First assume the alternating polar form
`A(x,y)=B(x,y)+B(y,x)` is nondegenerate. Its dimension `n` is even.
Write `B(x,y)=A(Tx,y)` for a unique linear operator `T`. With `T*`
denoting the adjoint for `A`, the identity `B+B^T=A` gives

```text
T+T*=I.                                                  (2)
```

We prove by induction on `n` that, for every integer `k` with
`0<=k<n/2`, there is a `k`-dimensional `B`-isotropic subspace.
For `k=0` there is nothing to prove. For `k>=1`, define quadratic
moments `m_j(v)=A(v,T^j v)` and impose the `k` equations

```text
m_1(v)=m_3(v)=...=m_(2k-1)(v)=0.                          (3)
```

They have a nonzero common solution. Indeed, the number of common
zeros modulo two is the sum over all `v in F_2^n` of
`prod_i(1+m_(2i-1)(v))`. This polynomial has degree at most `2k<n`.
Every monomial misses some variable, so its sum over all `v` is even.
The zero vector is a common zero; an even total number of zeros
therefore forces another one. This is the elementary binary parity
form of the Chevalley–Warning argument, proved here without invoking an
external theorem.

All even moments through `m_(2k-2)` then vanish too. For `r>=1`,
alternation and (2) give

```text
0=A(T^r v,T^r v)
 =A(v,(I+T)^r T^r v)
 =sum_(s=0)^r binom(r,s) m_(r+s)(v).                     (4)
```

The last term is `m_(2r)`; every other index is smaller than `2r`.
Induction on `r` proves the assertion from the odd equations (3).

Let `j` be the cyclic length of this nonzero `v`, so
`v,Tv,...,T^(j-1)v` are independent and span a `T`-invariant space if
`j<k`. For any `i,s<k`,

```text
B(T^i v,T^s v)
 =A(v,(I+T)^(i+1) T^s v),                              (5)
```

which is a sum of moments with indices at most `i+s+1<=2k-1`.
Thus the first `k` Krylov vectors are mutually `B`-orthogonal whenever
they are independent. They are also mutually `A`-orthogonal, by the
analogous formula with exponent `i` instead of `i+1`.

If `j<k`, set `W=span(v,Tv,...,T^(j-1)v)`. It is `T`-invariant and
both `A`- and `B`-isotropic. Its `A`-orthogonal complement `W^perp`
is `T`-invariant: `(I+T)W` is contained in `W`, and
`A(Ty,w)=A(y,(I+T)w)=0` for `y in W^perp,w in W`.
Moreover `B` vanishes between `W` and `W^perp` in both orders.
Consequently `B` descends to `W^perp/W`; the induced polar form is
nondegenerate alternating, and the quotient has dimension `n-2j`.
Since `n-2j>2(k-j)`, induction supplies a `(k-j)`-dimensional
isotropic subspace there. Its full preimage in `W^perp` is a
`k`-dimensional `B`-isotropic subspace. Taking `k=n/2-1` proves (1)
when `A` is nondegenerate.

## Polar radical

Now let `R=rad A` be nonzero, and use induction on `n`, with `n<=3`
as the trivial base case because `d(n)=0`. On `R`, the
self-pairing `Q(r)=B(r,r)` is linear because its polar form `A`
vanishes there. If some nonzero `r in R` has `Q(r)=0`, then
`B(r,x)=B(x,r)` for every `x`. Put `f(x)=B(r,x)`.

* If `f=0`, `B` descends to `V/<r>`. Lift an isotropic subspace of
  dimension `d(n-1)` and include `r`, obtaining dimension
  `d(n-1)+1>=d(n)`.
* If `f!=0`, the hyperplane `H=ker f` contains `r`, and `r` is
  `B`-orthogonal to `H` in both orders. The restriction of `B` descends
  to `H/<r>`, whose dimension is `n-2`. Lift an isotropic subspace
  there and include `r`, obtaining `d(n-2)+1=d(n)`.

If no nonzero `r in R` has `Q(r)=0`, linearity of `Q|R` forces
`dim R=1`. Write `R=<r>` with `B(r,r)=1`. The rank of alternating `A`
is even, so `n` is odd. The hyperplane `H=ker B(r,.)` excludes `r`,
and `V=<r> direct_sum H` is `B`-orthogonal in both orders. The
restriction `A|H` is nondegenerate: any element of its radical is in
`R intersect H=0`. The nondegenerate-polar result on the even space
`H` gives dimension `(n-1)/2-1=d(n)`. This completes the proof.

The bound cannot in general be raised to `floor(n/2)` without further
hypotheses. At `n=4`, take `B(x,y)=x dot Ly` with binary columns of
`L` equal to `(1,2,4,9)`, so its matrix is `I+E_(1,4)` and `L` is
invertible. The nonzero self-isotropic vectors are exactly
`{3,5,6,10,11,12,13}` (binary row indexing), and no distinct pair
among them is orthogonal in both orders. Its largest `B`-isotropic
subspace has dimension one.

## Affine Walsh consequence

Suppose `a(x)=Lx+b` on `F_2^t`, `rank L>=t-1`, and a controlled number
of rows may subsequently have their labels changed. For `t>=6`, take
the direction `V` from (1), of dimension
`d=floor((t-2)/2)>=2`. Every restriction `a(x+v)|V` is constant on
each coset `x+V`, since `r dot Lv=0` for all `r,v in V`. The induced
affine map from the `2^(t-d)` cosets to `V*` has nonzero linear part:
otherwise `V subset ker L^T`, contradicting
`dim ker L^T<=1<d`. At least half the cosets therefore have a nonzero
common restriction. Thus there are at least

```text
2^(t-d-1)                                                  (6)
```

disjoint aligned cosets, and changing fewer than that many row labels
leaves one intact, provided the changed assignment still uses nonzero
labels and distinct physical flip columns with enough copies.

For the repeated-Walsh nearby-prime endpoint profile and any fixed arc
constant `C`, the character calculation already proved in
[the affine-assignment note](walsh_affine_assignment_plane_certificate.md)
forces `b_copy>h=2^d` for sufficiently large `t`. Here
`h>=sqrt(M)/(2 sqrt(2))`, `M=2^t`. Distinct physical prime norms then
give `r=b_copy(M-1)>=c M^(3/2)` and
`log N>=log((r+1)!)>=c' M^(3/2) log M` for absolute positive constants.
Consequently, within this affine or sufficiently lightly modified
nearby-prime profile,

```text
M=O_C((log N/loglog N)^(2/3)).                            (7)
```

The prime-log spread, copy capacity, assignment rank, and surviving
aligned coset are essential scope conditions. This result does not
prove a uniform endpoint count or apply to arbitrary nonlinear maps.

The [checker](check_walsh_arbitrary_binary_bilinear_isotropic_bound.py)
tests exhaustive four-dimensional forms, random higher-dimensional
forms, and the explicit half-dimensional counterexample.
