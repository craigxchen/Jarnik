# Weighted symplectic phase averaging gives unconditional two-thirds radius growth

For the symplectic assignment of one flipped entry per Walsh row, every
actual fixed-constant endpoint realization satisfies

```text
M = O_C(1+(log R/loglog R)^(2/3)),
```

with the logarithmic expression interpreted for sufficiently large `R`.
The physical split primes may have arbitrary distinct sizes and arbitrary
Gaussian-prime orientations. No flipped-prime mass assumption is needed:
phase averaging forces the original logarithmic radius to be at least
order `sqrt(M)` times the sum of the `M-1` distinct flipped-prime logs.
The factorial lower bound for that sum gives the conclusion directly.

This is a theorem for the specified symplectic repeated-Walsh family,
not the uniform endpoint theorem. It extends the
[aligned-coset phase argument](linear_support_character_height_obstruction.md)
by averaging exact heights instead of imposing nearly equal weights.
A relative mass assumption is retained only for the optional stronger
conclusion that the number of physical copies grows with `sqrt(M)`.

## 1. Source model and an optional weight hypothesis

Let `t=2d>=2`, `M=2^t`, `h=2^d=sqrt(M)`, and `b>=5`. Index rows by
`x in F_2^t`, and let `J` be a nonsingular alternating binary matrix.
Thus `B(x,y)=(Jx) dot y` is a nondegenerate alternating bilinear form.
There are `b` physical copies of each nonzero Walsh label `a`, giving
`r=b(M-1)` columns with unflipped entries `chi_a(x)=(-1)^(a dot x)`.

For each nonzero row `x`, flip one physical column of label `a_x=Jx`.
These labels are all distinct. Row zero flips a distinct physical column
of any nonzero label; a second copy makes that possible. Each physical
column is flipped at most once. Attach arbitrary distinct split primes
`p_j` to these physical columns, and choose Gaussian primes `pi_j` above
them. With the resulting signs `s_xj`, form the actual Gaussian rows

```text
z_x = g epsilon_x product_j pi_j^((1+s_xj)/2)
                            bar(pi_j)^((1-s_xj)/2),
epsilon_x in {1,-1,i,-i},        g != 0 in Z[i].
```

All rows have the same norm. Put

```text
W_0=sum_j log p_j,       D=log Norm(g),
log R^2=W_0+D,
F=sum_(x!=0) log p_(flipped column at row x).
```

In particular `F` excludes the flipped column at row zero. The optional
weight assumption used only for the copy-count corollary is

```text
F >= eta W_0/b,                    eta>0 fixed.                 (1)
```

It is not automatic for arbitrary prime assignments. The radius growth
theorem below does not require (1).

## 2. Symplectic cosets and their exact characters

A Lagrangian is a `d`-dimensional subspace `V` on which `B` vanishes.
It has `h` elements and equals its `B`-orthogonal complement. Choose a
coset `P=x_0+V` with `x_0 notin V`, and define

```text
c_(x_0+v)=(-1)^B(x_0,v) for v in V,       c_x=0 outside P.
```

The character `B(x_0,.)` on `V` is nonzero. Hence `sum c_x=0` and
`sum |c_x|=h`. For an unflipped Walsh label,

```text
u_a=(1/2)sum_x c_x chi_a(x)
   =(h/2)chi_a(x_0) if a in J(P), and 0 otherwise.             (2)
```

Indeed the nonzero restriction condition is `a dot v=B(x_0,v)` for
all `v in V`, whose solution set is `Jx_0+JV=J(P)`. It has `h` elements
and excludes zero.

For each selected row `x in P`, its assigned label `Jx` lies in `J(P)`.
Moreover `chi_(Jx)(x)=1` and `chi_(Jx)(x_0)=c_x`. Flipping its entry
therefore changes its physical character coefficient from `(h/2)c_x`
to `(h/2-1)c_x`, reducing its absolute value by one. This happens for
both signs of `c_x`. Every source row has one assigned flipped entry;
the rows with `c_x=-1` are not a separate set of assigned flips. All
unselected rows have coefficient zero in this correction, and row zero
is never in `P`.

Write the resulting integral physical character as `v_j=(sum c_x s_xj)/2`.
Its Gaussian composite `Beta` is conjugate-primitive, because the physical
primes are distinct and each appears in just one orientation in that
composite. It is a nonunit: its absolute physical coefficient sum is
`bM/2-h>0`. Its exact logarithmic norm is therefore

```text
V_Beta=sum_j |v_j| log p_j.                                  (3)
```

## 3. Averaging uses combinatorial symmetry, not equal prime sizes

Average these certificates over the images of one pair `(V,P)` under
the finite symplectic group. That group acts transitively on nonzero
vectors: a nonzero vector can be paired with a vector of bilinear product
one, the pair can be extended to a symplectic basis, and a map between
two such bases is an isometry. It also acts transitively on the nonzero
labels, since `JT=T^(-transpose)J` for every symplectic transformation `T`.

Every support `P` and every selected-label set `J(P)` has `h` elements
and excludes zero. Thus each nonzero row and each nonzero label is
selected with probability `h/(M-1)`. Equation (2) and the unit reductions
at the selected flips give

```text
E V_Beta = M W_0/[2(M-1)] - h F/(M-1),
E(V_Beta-W_0/2) = (W_0/2-hF)/(M-1).                         (4)
```

All prime weights remain fixed during this averaging. Only valid
character certificates for the same actual source are varied; no
symmetry of the prime values or physical Gaussian arguments is assumed.

If the source lies on an arc of angular width
`Delta<=C exp(-(W_0+D)/4)`, the exact signed product is

```text
product_x z_x^c_x = unit * Beta/bar(Beta).
```

Lift arguments in the containing arc. The zero row sum gives
`|sum c_x theta_x|<=h Delta/2`, so the argument of `Beta` is within
`h Delta/4` of the grid `(pi/4)Z`. A conjugate-primitive nonunit Gaussian
integer lies on none of the axis or diagonal lines of that grid. Its
distance to such a line is at least `1/sqrt(2)`, giving the necessary
height inequality for every certificate

```text
V_Beta-W_0/2 >= D/2+log 8-2log(hC).                          (5)
```

This step uses integer coordinates and no transcendence estimate.
Combining (4) and (5), and using `h^2=M`, proves the exact necessary
weight condition

```text
hF <= W_0/2+(M-1)(log(M C^2/8)-D/2).                         (6)
```

The common Gaussian content strengthens (6). Arbitrary row units have
already been accounted for by the half-unit phase grid.

## 4. Unconditional growth with the original radius

Rearrange (6) using `W=log R^2=W_0+D`, retaining the original radius:

```text
W >= 2 sqrt(M) F - 2(M-1)log(M C^2/8) + M D.                 (7)
```

The `M-1` flipped nonzero rows use pairwise distinct physical split
primes. On ordering these primes, the j-th is at least `j+1`. Hence

```text
F >= log(M!).                                               (8)
```

Combining (7)--(8) gives the exact necessary radius bound

```text
log R^2 >= 2 sqrt(M) log(M!)
             - 2(M-1)log(M C^2/8) + M log Norm(g).           (9)
```

This uses actual source phase through (5), not only the norm profile.
It is uniform in the number of copies, in the prime sizes, and in their
orientations. The common Gaussian content strengthens the inequality.

For an explicit large-`M` threshold put
`K_C=1+2log^+(C)/log 4`. For `M>=4`,
`log(M C^2/8)<=K_C log M`, while the last `M/2` factors of `M!` give

```text
log(M!) >= (M/2)log(M/2) >= (M/4)log M.
```

Consequently (9) implies

```text
W >= M log M (sqrt(M)/2-2K_C) + M D.
```

When `sqrt(M)>=8K_C`, this yields

```text
W >= (1/4) M^(3/2) log M.                                   (10)
```

Inverting gives `M=O((W/log W)^(2/3))`. Explicitly, for sufficiently
large `W`, if `M>W^(1/3)` then `log M>(log W)/3`, so (10) gives
`M<=(12W/log W)^(2/3)`. If `M<=W^(1/3)`, the same order bound holds.
The bounded range `sqrt(M)<8K_C` is absorbed by the constant depending
on `C`. Thus, without the optional hypothesis (1),

```text
M = O_C(1+(log R/loglog R)^(2/3)).                            (11)
```

Small flipped-prime mass relative to `W_0` is harmless for this radius
conclusion: the phase inequality still forces a large `W_0`, because
there are `M-1` distinct primes in `F`. It need not force a large copy
count. The structured symplectic flip assignment remains an essential
hypothesis and has not been obtained for arbitrary endpoint tuples.

## 5. Conditional copy-count corollary

If the additional mass hypothesis (1) does hold, the original averaging
also forces the copy count itself to grow. Distinct prime norms give

```text
W_0 >= log((r+1)!) >= (r/2)log(r/2)
    >= (bM/4)log M,                                         (12)
```

using `r=b(M-1)>=bM/2` and `b>=5`. Apply (1) to (6), discard `D>=0`,
and use (12) to obtain

```text
eta sqrt(M) <= b/2+4K_C,
b >= 2eta sqrt(M)-8K_C.                                     (13)
```

In particular `sqrt(M)>=8K_C/eta` implies `b>=eta sqrt(M)`.
This conditional conclusion is stronger than the unconditional radius
bound in a different direction; it is not needed to prove (11).

The [checker](check_walsh_symplectic_weighted_phase_growth.py) enumerates
the small Lagrangian cosets, verifies row and label incidence counts,
checks every selected physical coefficient reduction, and tests the
exact linear weight identity. Its literal Gaussian fixtures include
arbitrary row units and nontrivial common content. Root and Astra
independently audited the averaging and the radius-growth deduction.
