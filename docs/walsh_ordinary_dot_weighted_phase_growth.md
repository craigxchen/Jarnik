# Unconditional weighted radius growth for ordinary-dot Walsh flips

The repeated-Walsh profile with `a_x=x` for every nonzero row also
satisfies unconditional two-thirds radius growth in every binary dimension.
The canonical assignment also follows by restriction to even-index rows.
Arbitrary distinct physical split primes, their Gaussian-prime
orientations, row units, and common Gaussian content are allowed. The
ordinary binary dot product is nonalternating, so the uniform incidence
identity from [the symplectic note](walsh_symplectic_weighted_phase_growth.md)
does not apply. Pair certificates supply exactly the missing parity
constraint.

This is a source-phase theorem for this specified profile, not a bound
for general endpoint tuples.

## Even binary dimension: source and parity weights

Let `t>=4` be even, `M=2^t`, `h=sqrt(M)`, and let `b>=5` physical
copies of each nonzero Walsh label have signs `chi_a(x)=(-1)^(a dot x)`.
Each nonzero row `x` flips one physical column of label `x`, with all
physical flipped columns distinct. Row zero flips any further unused
physical column. Attach pairwise distinct split primes and construct the
actual Gaussian rows as in the symplectic note, including arbitrary row
units and a common factor `g`. Write

```text
W_0=sum_j log p_j,  D=log Norm(g),  W=log R^2=W_0+D,
w=(1,...,1),  E=w^perp={even-parity binary vectors}.
```

Partition physical-prime weight by label:

```text
O = sum_(a odd) sum_(physical copies of a) log p_j,
A = sum_(a in E\{0,w}) sum_(physical copies of a) log p_j,
Q = sum_(physical copies of w) log p_j,
W_0=O+A+Q.
```

Let `F_o` be the sum of flipped-prime logs at odd rows, and `F_e` the
sum at rows in `E\{0,w}`. These definitions use the fixed actual prime
assignment. No prime or angle is permuted during averaging.

Assume the actual rows lie on an arc of length at most `C sqrt(R)`.
Every nonzero integral signed row certificate `c` below has `sum c=0`
and integral physical coefficients `v_j=(sum_x c_x s_xj)/2`. Its
composite is conjugate-primitive and nonunit. For `m=sum|c_x|`, the
exact product and lifted source arguments give

```text
product_x z_x^(c_x)=unit * Beta/bar(Beta),
log Norm(Beta) >= W_0/2+D/2+A_m,
A_m=log 8-2log(m C).                                        (1)
```

The last inequality is the integer-coordinate distance from `(pi/4)Z`,
not a norm-only assumption.

## Even-coset certificates

The dot product restricted to `E` has radical `<w>`, and the quotient
`E/<w>` is a nondegenerate alternating space of dimension `t-2`.
Lift any Lagrangian of this quotient to `V⊂E`. Then
`dim V=t/2`, `w in V`, and `V=V^perp` for the full dot product.
For an even nonzero coset `P=x0+V⊂E`, with `x0 notin V`, put

```text
c_(x0+v)=(-1)^(x0 dot v),  c=0 outside P.
```

The character on `V` is nonzero; hence the coefficients sum to zero
and their absolute sum is `h`. The unflipped Walsh coefficient is
`(h/2)chi_a(x0)` exactly on labels `a in P` and zero elsewhere.
For each row `x in P`, its flip at label `x` reduces its coefficient's
absolute value by one. All selected rows and labels are even and avoid
both `0` and `w`.

Choose a Lagrangian of `E/<w>` uniformly and then an even nonzero coset
of its lift uniformly. The quotient has `h/2` elements per Lagrangian.
Each nonzero quotient vector belongs to a random Lagrangian with
probability `1/(h/2+1)=2/(h+2)`. There are `h/2-1` eligible cosets.
Thus every even row or label other than `0,w` is selected with
probability

```text
q = (1-2/(h+2))/(h/2-1) = 2h/(M-4).
```

Averaging exact physical heights yields

```text
E log Norm(Beta) = M A/(M-4)-2h F_e/(M-4).
```

Apply (1) with `m=h`:

```text
M A-2h F_e >= (M-4)(W_0/2+D/2+A_h).                        (2)
```

## Odd-pair certificates cancel the exceptional parity mass

For every odd row `x`, its mate `x+w` is also odd. Use
`c=e_x-e_(x+w)`. The unflipped character has coefficients of absolute
value one on odd labels and zero on even labels. Both assigned flipped
labels `x,x+w` are odd, and each flip cancels its physical coefficient.
Therefore the exact composite height is

```text
log Norm(Beta)=O-f_x-f_(x+w),
```

where `f_x` is the flipped-prime log at row `x`. This composite is a
nonunit: there are `bM/2-2>0` remaining physical exponents. Equation
(1) with `m=2` and the average over the `M/4` disjoint odd pairs give

```text
O >= W_0/2+D/2+A_2+4F_o/M,
A <= W_0/2-Q-D/2-A_2-4F_o/M.                               (3)
```

Insert the upper bound (3) into (2). Rearranging with the original
radius `W=W_0+D` gives the exact necessary inequality

```text
W >= h F_e + 2F_o + (M/2)Q + (M/2)D
       + [M A_2+(M-4)A_h]/2.                               (4)
```

All terms involving actual prime weights on the right are nonnegative.
The odd-pair constraint is what permits this weighted conclusion; even
cosets alone do not eliminate arbitrary concentration on odd labels.

## Radius growth

There are `n=M/2-2` distinct flipped primes in `F_e`, so
`F_e>=log((n+1)!)`. For `M>=16`, use
`n+1=M/2-1>=M/4`, and retain the last half of these factorial factors:

```text
F_e >= (M/8)log(M/8) >= (M/32)log M.
```

The second inequality uses `M>=16`. Set
`L_C=1+4log^+(C)/log 16`. Since

```text
[M A_2+(M-4)A_h]/2
 = (M/2)(log 2-2log C)
   + ((M-4)/2)(log 8-log M-2log C)
 >= -(M/2)log M-2M log^+C
 >= -L_C M log M,
```

(4) gives

```text
W >= M log M (sqrt(M)/32-L_C).
```

In particular `sqrt(M)>=64L_C` implies
`W >= M^(3/2)log M/64`. Inversion, with the bounded remaining `M`
absorbed in the constant, proves

```text
M = O_C(1+(log R/loglog R)^(2/3)).                            (5)
```

The case `t=2` has bounded `M=4` and is harmless for (5). No relative
mass assumption or nearly equal prime-size hypothesis is required.

## Odd binary dimension

Now let `t>=3` be odd. Set `m=M/2`, `h=sqrt(m)`, and retain
`w=(1,...,1)` and `E=w^perp`. This time `w dot w=1`, so `w` lies
outside `E` and the restricted dot product on `E` is a nondegenerate
alternating form. The rows of `E` are a symplectic space of size `m`.

Two global labels have the same restriction to `E` precisely when they
differ by `w`. The only nonzero global label constant on `E` is `w`.
Thus every nonzero restricted Walsh label has `2b` physical columns.
Let `Q` denote the total physical-prime log weight of global class `w`,
and let `F_E` denote the flipped-prime log sum over `E\{0}`.
For every nonzero row `x in E`, its assigned global label `x` restricts
to its symplectic dual on `E`. These assigned physical columns are
pairwise distinct. Flips outside `E` have no effect on these rows.

Use the symplectic Lagrangian coset certificates on `E`. Their supports
avoid zero. In particular any exceptional flip at row zero contributes
nothing, even when its label is `w`. Their exact average height is

```text
E log Norm(Beta)
 = m(W_0-Q)/[2(m-1)]-h F_E/(m-1).
```

The original source phase still gives (1) with `m` there replaced by
the row coefficient sum `h`; the original logarithmic radius remains
`W=W_0+D`. Consequently

```text
W >= 2 sqrt(m) F_E + m(Q+D) + 2(m-1)A_h.                    (6)
```

There are `m-1` distinct flipped primes in `F_E`, so
`F_E>=log(m!)`. Since `A_h=log 8-log m-2log C`, (6) is exactly the
symplectic radius inequality at size `m`, with extra nonnegative term
`m Q`. It gives `W=Omega_C(m^(3/2)log m)` for sufficiently large `m`.
As `m=M/2`, (5) follows for odd `t` as well.

## Constant columns and an unused zero row

Every certificate in both parity cases has `c_0=0`. The proofs therefore
remain valid if the row at zero is omitted entirely, or if its signs are
arbitrary. They require the asserted flip pattern only at nonzero rows
used by the certificates. Physical columns constant on those rows have
zero character coefficient and may be included with any total weight
`Q_const>=0`: replace `D` in the displayed inequalities by
`D_eff=D+Q_const`, and define `W_0` there using the nonconstant columns.
Then the original radius remains `W=W_0+D_eff` and every deduction
continues unchanged. This use of `D_eff` is an exact height/phase
accounting identity; it need not be a common Gaussian factor of an
unused exceptional zero row.

The same observation holds for the symplectic proof, whose certificates
also avoid zero. Reorienting any whole physical column merely conjugates
its Gaussian prime, leaving all stated heights and phase arguments valid.

## Canonical assignment with arbitrary prime weights

Consider the canonical assignment
`a_x=1+(x mod(M-1))` on `M=2^t` rows, with distinct assigned physical
columns as in the earlier profile construction. Restrict the source to
the even-index rows `x=2y`, where `y` ranges over `F_2^(t-1)`.
For every nonzero such row, `a_(2y)=2y+1`. Its restriction to these
rows is the ordinary Walsh label `y`.

The two global labels `2a,2a+1` restrict to the same nonzero label `a`,
so each nonzero restricted label has `2b` physical copies. Global class
`1` is constant on all these rows except for the assigned flip at row
zero. The certificates avoid row zero, so class `1` contributes zero
to every certificate. Its entire logarithmic weight can be charged to
`D_eff` as just described. Flips at odd-index rows do not affect the
restricted rows.

The ordinary-dot theorem in dimension `t-1` therefore applies with
`M/2` selected rows and the unchanged original radius. Both parity cases
are available, so the canonical assignment also satisfies (5), uniformly
in arbitrary distinct physical prime sizes and orientations. The finitely
many small dimensions are absorbed by the constant depending on `C`.
This removes the nearly equal prime-log assumption from the canonical
and relabelled families in the earlier
[linear-support note](linear_support_character_height_obstruction.md).

The [checker](check_walsh_ordinary_dot_weighted_phase_growth.py) checks
the exact even-coset incidences, pair cancellations, arbitrary-weight
identities, and literal Gaussian signed products. It does not assert
that its source fixtures are endpoint clusters.
