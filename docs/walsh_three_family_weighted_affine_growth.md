# Three certificate families give arbitrary-weight affine radius growth

For repeated Walsh profiles, arbitrary physical prime weights can be
handled for **all binary affine flip assignments of bounded corank**, with
a quantitative allowance for exceptions. Symmetry of the affine linear
part and nearly equal prime logarithms are unnecessary.

Let `M=2^t>=2`, let `b>=5` physical copies of each nonzero Walsh label
be present, and flip one physical column at each row. All assigned
physical columns must be distinct. Suppose the actual assigned labels
agree with `a(x)=Lx+a_*` except at at most `K` rows, where
`rank L>=t-k`. Every actual assigned label is nonzero; rows where the
affine formula gives zero are necessarily exceptions. Attach arbitrary
pairwise distinct split primes and choose arbitrary Gaussian-prime
orientations. Independent row units and common Gaussian content are
allowed. If the actual Gaussian rows lie on an arc of length at most
`C sqrt(R)`, then

```text
log R^2 >= c_C M^2 log M / (2^k sqrt(M)+K).                  (1)
```

The constant depends only on `C`, not on `b,k,K`, the prime sizes, or
orientations. This is a theorem for the stated repeated-Walsh profiles;
it does not establish the uniform endpoint lattice-point count.

There is also a rank-free consequence for every affine assignment:

```text
log R^2 >= c'_C M^2 log M / (M^(3/4)+K).                    (1a)
```

Thus exact affine assignments, or `K=O(M^(3/4))` exceptions, satisfy
an unconditional within-profile exponent `4/5` with arbitrary prime
weights, without any rank hypothesis. A fixed corank improves this to
`2/3` when the edit allowance is `O(sqrt(M))`.

For fixed `k`, (1) gives the usual two-thirds radius growth when
`K=O(sqrt(M))`. If `K=O(M^delta)` with fixed `1/2<delta<1`, it gives

```text
M = O_(C,k,delta)(1+(log R/loglog R)^(1/(2-delta))),          (2)
```

with dependence also on the constant in `K=O(M^delta)`. Thus every fixed
polynomially sublinear edit allowance gives an exponent below one.

## Actual source phase controls individual label masses

Write `s_xj` for the physical signs after the flips and form

```text
z_x = g epsilon_x product_j pi_j^((1+s_xj)/2)
                            bar(pi_j)^((1-s_xj)/2),
W_0=sum_j log p_j,  D=log Norm(g)>=0,  W=log R^2=W_0+D.
```

For each nonzero label `a`, let `W_a` be its total physical-prime log
weight and `F_a` its total assigned flipped-prime log weight. Put
`F_tot=sum_a F_a`, so `F_a<=W_a` and `F_tot<=W_0`.
For a row, `f_x` denotes its assigned flipped-prime log.

Every row-pair certificate `c=e_x-e_y` has a conjugate-primitive nonunit
Gaussian composite `Beta`. Its signed source product is exactly
`unit * Beta/bar(Beta)`. Lift row arguments in the containing arc. The
integer-coordinate gap from `(pi/4)Z` gives

```text
log Norm(Beta) >= (W_0+D)/2+log 2-2log C.
```

The weighted physical row Gram `G_xy=sum_j s_xj s_yj log p_j` satisfies
`log Norm(Beta)=(W_0-G_xy)/2`. Hence the actual phase premise implies

```text
G_xy <= kappa=4log C-log 4-D   for every x!=y.               (3)
```

For `d!=0`, average `G_(x,x+d)` over all rows. Each physical column
is flipped at most once, so its two affected pair incidences give

```text
E_x G_(x,x+d) = sum_(a!=0) q_a (-1)^(a dot d),
q_a=W_a-4F_a/M > 0       (M>4).
```

Write this transform as `qhat(d)`, and put
`u_d=kappa-qhat(d)>=0` for `d!=0` and `Q=sum_a q_a`.
Walsh inversion, with the absent zero-label coefficient set to zero,
gives the exact identities

```text
sum_(d!=0) u_d = Q+(M-1)kappa,
q_a+kappa = (2/M)sum_(d:a dot d=1) u_d.
```

Using `Q<=W_0`, `F_a<=W_a`, and
`kappa_+=max(0,4log C-log 4)` yields the uniform class bound

```text
W_a <= B_max := 2W_0/(M-4)+(M-2)kappa_+/(M-4).              (4)
```

This is a consequence of actual pair phases. It is not a free
regularity assumption on the physical prime weights.

## Isotropic directions and the affine zero fibre

By the [binary bilinear isotropic theorem](walsh_arbitrary_binary_bilinear_isotropic_bound.md),
the form `B(v,u)=v dot Lu` has a totally isotropic subspace `V_0` of
dimension at least `max(0,floor((t-2)/2))`. Choose a subspace
`V⊂V_0` complementary to `V_0 intersect ker L^T`, and, if necessary,
shrink it to dimension

```text
d=max(0,floor((t-2)/2)-k),  h=2^d.                          (5)
```

Then `V` is still totally isotropic and `L^T|V` is injective. Thus the
map `x -> (Lx+a_*)|V` is onto `V*`. Choose a row `x_0` in its zero
fibre and translate row coordinates by `x_0`.

This does **not** require solving `Lx_0=a_*`, which need not have a
solution when `L` is singular. After translation the residual affine
term `a_*'=Lx_0+a_*` only has to annihilate `V`. The baseline physical Walsh columns
acquire fixed orientation signs `(-1)^(a dot x_0)`, absorbed by
conjugating the corresponding Gaussian primes; all weights stay fixed.
The good-coset calculation below is unaffected by these orientations.

In translated coordinates put

```text
H={x:(Lx)|V=0},   |H|=M/h=nh,   n=M/h^2.
```

Surjectivity gives the stated size; isotropy gives `V⊂H`. Assume for
now `M>=16` and `h>=2`. The choice (5) ensures `n>1` (indeed `n>=4`).
No condition on actual assigned labels inside `H` will be needed.

For any signed row certificate used below, let
`v_j=(sum_x c_x s_xj)/2` and `V_Beta=sum_j |v_j|log p_j`.
Every certificate has coefficients in `{0,1,-1}`, zero row sum, and
integral `v_j`. With `m=sum|c_x|`, the actual source phase gives

```text
V_Beta >= W_0/2+D/2+A_m,
A_m=log 8-2log(m C).                                        (6)
```

Indeed `product_x z_x^c_x=unit*Beta/bar(Beta)`, and the argument of
`Beta` is within `m C exp(-W/4)/4` of `(pi/4)Z`. Its distance from
that grid is at least `exp(-V_Beta/2)/sqrt(2)` by integer coordinates.
All physical split primes are distinct, so the composite is
conjugate-primitive. Nonunit status is checked below.

## Three families and their exact baseline coefficients

Partition nonzero labels into

```text
I = {a:a notin V^perp},
II = V^perp minus H^perp,
III = H^perp minus {0}.
```

Let `W_I,W_II,W_III` be their total actual physical-prime weights.
Put `F_out=sum_(x outside H) f_x`, `F_H=sum_(x in H) f_x`, and let
`F_bad` be the flipped-prime weight of exceptional rows outside `H`.

**Family g: good affine cosets of `V`.** Each coset `P=x+V` outside
`H` has a common nonzero predicted restriction
`ell=(Lx+a_*')|V`. Put
`c_(x+v)=(-1)^ell(v)` on `P`, and zero elsewhere. There are
`n(h-1)` such cosets, and every nonzero restriction occurs on exactly
`n` of them. Each label in `I` therefore has expected baseline absolute
coefficient `h/[2(h-1)]`; other labels have zero coefficient.

At a nonexceptional selected row the assigned flip reduces its absolute
physical coefficient by one. Any exceptional flip changes absolute
coefficient by at most one, so it loses at most two relative to this
reduction. The exact baseline and correction bound are

```text
E_g V_Beta <= h W_I/[2(h-1)]
              - F_out/[n(h-1)] + 2F_bad/[n(h-1)].           (7)
```

**Family p: differences of two `V`-cosets inside `H`.** Choose an
unordered pair of distinct cosets among the `n` such cosets; use `+1`
on one and `-1` on the other. The orientation choice does not matter
for absolute coefficients. A nonzero character of `H/V` separates
exactly `n^2/4` unordered pairs, a proportion `n/[2(n-1)]`.
Hence a label in `II` has expected baseline coefficient
`hn/[2(n-1)]`; every other label has zero coefficient.

A row of `H` is selected with probability `2/n`. Bound each selected
flip's change in absolute coefficient by `+1`, without any assumption
on its actual assigned label, including the translated origin:

```text
E_p V_Beta <= hn W_II/[2(n-1)] + 2F_H/n.                   (8)
```

**Family f: differences of two entire `H`-cosets.** Choose an unordered
pair among the `h` cosets of `H` in the full row space, with coefficients
`+1` and `-1`. Only labels in `III` can occur. Their expected baseline
absolute coefficient is `M/[2(h-1)]`. Each row is selected with
probability `2/h`, so arbitrary actual flips obey

```text
E_f V_Beta <= M W_III/[2(h-1)] + 2F_tot/h.                 (9)
```

The respective row coefficient sums are `h,2h,2nh`. Each unflipped
physical absolute-coefficient sum is `bM/2`. At most one unit can be
lost per selected flip. Thus every resulting composite is nonzero and
nonunit: even for the largest support,

```text
bM/2-2nh = nh(bh/2-2)>0.
```

This also verifies (6) for every member of all three families. The
origin and any other exceptional row can occur in families p and f;
its actual flip is explicitly included in their error bounds. No
exceptional-row prime is discarded or charged to fictitious common
content.

## The mixture closes the missing label block

Use the positive, not quite probability, weights

```text
lambda_g=(h-1)/h,
lambda_p=(n-1)/(hn),
lambda_f=(h-1)/M.
```

They satisfy

```text
lambda_g+lambda_p+lambda_f=1-1/M.                           (10)
```

Each of the three baseline coefficients in (7)--(9), multiplied by
its lambda, is exactly `1/2`. Thus the mixed baseline is `W_0/2`.
The mixed correction is at most

```text
-hF_out/M + 2hF_bad/M
 + 2(n-1)F_H/(hn^2) + 2(h-1)F_tot/(Mh).                    (11)
```

For the mixed lower bound, apply (6) to each family and multiply by
its lambda. Since every row coefficient sum is at most `2M/h<=M`,

```text
-sum lambda A_m <= 2log M+2log^+ C.                         (12)
```

The common-content contribution is nonnegative and may be dropped.
Combining (10)--(12), multiplying by `M`, and retaining all weights gives

```text
hF_out <= W_0/2 + 2hF_bad
 + [2M(n-1)/(hn^2)]F_H + [2(h-1)/h]F_tot
 + M(2log M+2log^+ C).                                    (13)
```

The row sets contributing to `F_H` and `F_bad` contain at most `M/h`
and `K` rows, respectively. Their flipped physical primes are distinct.
Even if many rows use the same label, the total flipped weight at that
label is at most its full class weight. Therefore (4) implies

```text
F_H <= (M/h)B_max,   F_bad<=K B_max,   F_tot<=W_0.          (14)
```

No injectivity of the actual label assignment is being assumed here.
Only physical flipped-column distinctness is used.

## Quantitative exceptions and growth

Set

```text
q=1+hK/M,
L_C=2+[(7/3)kappa_+ + 2log^+ C]/log 16.
```

For `M>=16`, use `M/(M-4)<=4/3` and
`(M-2)/(M-4)<=7/6` in (13)--(14). The `F_H` term costs at most
`(16/3)W_0+(7/3)M kappa_+`; the bad-row term costs at most
`(16/3)(q-1)W_0+(7/3)(q-1)M kappa_+`. Together with `W_0/2` and
`2F_tot`, the `W_0` coefficient is at most `8q`. Thus

```text
hF_out <= 8q W_0 + L_C q Mlog M.                           (15)
```

There are `M-M/h>=M/2` distinct flipped physical primes outside `H`.
Ordering them, and using a factorial lower bound, gives

```text
F_out >= log((M-M/h+1)!) >= (M/8)log M.                     (16)
```

If `h/q>=16L_C`, equations (15)--(16) yield

```text
W>=W_0 >= (h/q) Mlog M/128.
```

For the complementary range, distinctness of all
`r=b(M-1)>=5(M-1)` physical primes already gives

```text
W_0>=log((r+1)!) >= Mlog M.
```

Therefore `W_0 >= (h/q)Mlog M/(16L_C)` in that range. With
`c_C=min(1/128,1/(16L_C))`, both cases prove

```text
W >= c_C M^2 log M/(M/h+K).                                (17)
```

When `h=1` or `M<16`, the same inequality follows directly from the
baseline factorial bound, with the definition (5); for `M<16`, (5)
always has `h=1`. Finally (5) implies

```text
M/h <= 2sqrt(2) 2^k sqrt(M).
```

Substitution in (17), absorbing `2sqrt(2)` into the constant, proves
(1). For fixed corank and `K=O(M^delta)`, ordinary monotone inversion
of `M^(2-max(1/2,delta)) log M` gives (2) and the two-thirds case.
The full original radius has been retained throughout.

## Removing the rank hypothesis by physical copy capacity

In (1), take `k` to be the actual corank of `L` and put `x=2^k`.
The affine map has exactly `M/x` image labels. At least `M-K` actual
rows agree with it; their distinct assigned physical columns must lie
among those image classes, with at most `b` physical columns per class.
Therefore

```text
b >= x(1-K/M).                                             (18)
```

If zero is among the image labels, it has no physical columns; counting
it only weakens the displayed necessary bound. For `K<=M/2`, (18)
and the full-prime factorial bound give

```text
W>=W_0 >= (bM/4)log M >= (xM/8)log M.
```

If `x>=M^(1/4)`, this is at least `M^(5/4)log M/8`, which is at
least a constant times the right-hand side of (1a). If `x<M^(1/4)`,
use (1) and `x sqrt(M)<M^(3/4)` instead. If `K>M/2`, the baseline
`W_0>=Mlog M` proves (1a) directly because
`M^2/(M^(3/4)+K)<=2M`. This proves (1a) in all ranks.

Consequently, with no rank assumption, `K=O(M^delta)` gives

```text
M = O_C(1+(log R/loglog R)^(1/(2-max(3/4,delta))))
```

for fixed `delta<1`, with dependence also on `delta` and the edit-count
constant. At `delta<=3/4` the exponent is `4/5`; at `3/4<delta<1`
it is `1/(2-delta)`. These remain results for the repeated-Walsh affine
or edited-affine family, not arbitrary endpoint profiles.

The [checker](check_walsh_three_family_weighted_affine_growth.py) checks
all three baseline families, signed source products and flip corrections,
Walsh Gram inversion, the mixture coefficients, deficient-rank affine
zero fibres, arbitrary exception rows, and the quantitative inequalities.
Its literal Gaussian fixtures are not asserted to be endpoint clusters;
it verifies the exact identities on which the necessary phase argument
rests.
