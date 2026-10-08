# Three mixed cuts retain content without the old anchor

The old core cuts `124`, `135`, and `145` each force at least one half
of their logarithmic norm into the common Gaussian matching factor of
the six transformed points alone, up to negligible correction loss.
This assertion does not retain the original anchor and does not use the
experimental valuation-grid minima.

## 1. A standalone sorted-valuation inequality

Let the six signed phase valuations, measured in any common reference,
be ordered as `d_1<=...<=d_6`. For the least circle of these six phases
alone the exponent is `n=d_6-d_1`. Its common matching-content contribution,
divided by `log p`, is at least

```text
Gamma=(3n+d_1+d_2+d_3-d_4-d_5-d_6)/2
     =((d_6-d_5)+(d_2-d_1)+(d_6-d_4)+(d_3-d_1))/2.       (1)
```

All four terms in the last expression are nonnegative. Consequently
four phases at least `e`, together with one phase at most zero, force
`Gamma>=e/2`. The same holds for four phases at most zero and one phase
at least `e`. With four phases at least `e-E` and one at most `E`, or
four at most `E` and one at least `e-E`, the lower bound is `e/2-E`.

## 2. The exact Gaussian expressions

Use the two forms of the raw vectors from
[joint_circle_conductor_lift.md](joint_circle_conductor_lift.md):

```text
Q_4=H P_4-2 Delta_34 Delta_56 Delta_41 P_2
   =V P_4-2 Delta_34 Delta_56 Delta_42 P_1,
Q_5=F P_5-2 Delta_35 Delta_46 Delta_51 P_2
   =V P_5-2 Delta_35 Delta_46 Delta_52 P_1,
Q_6=H P_6+2 Delta_36 Delta_45 Delta_61 P_2
   =F P_6+2 Delta_36 Delta_45 Delta_62 P_1.                (2)
```

At an old cut `S` of exponent `e`, an inside `P_i` has Gaussian orders
`(e,0)`, an outside one `(0,0)`, before small corrections are added.
An old bracket has rational order `e` exactly when both its labels lie
inside `S`, and zero otherwise.

## 3. The cut `124`

The relevant matching orders give `v_p(F),v_p(H),v_p(V)>=e`.
For `Q_5` and `Q_6`, all three-bracket coefficients of the `P_1` and
`P_2` expressions in (2) have order zero. The reference vectors `P_1,P_2`
are inside the cut, and `P_5,P_6` are outside. Thus each first Gaussian
orientation has order at least `e`; the conjugate orientation has order
zero, since the scalar-times-reference term strictly dominates there.
Therefore the moving phases `d_5,d_6` are at least `e`.

Together with the fixed `P_1,P_2`, there are four phases at least `e`.
The fixed `P_3` has phase zero. Equation (1) gives `Gamma>=e/2`.
No restriction on the phase of `Q_4` is required.

## 4. The cut `135`

The two identities

```text
F=(A+C)-(A+B),       V=(A+C)+(A+B),
A+C=Delta_13 Delta_24 Delta_56,
A+B=Delta_12 Delta_35 Delta_46
```

show that both `F` and `V` have valuation at least `e`, because both
matching summands have order `e` at this cut.

Use the `P_1` expressions for `Q_4,Q_6` in (2). Their scalar-times-old-row
terms have valuation at least `e`, their triple-bracket coefficients have
order zero, and `P_1` is inside while `P_4,P_6` are outside. As in the
previous case, the moving phases `d_4,d_6` are at least `e`.
Together with fixed `P_1,P_3` this gives four high phases; fixed `P_2`
has phase zero. Hence `Gamma>=e/2`, independently of `Q_5`.

## 5. The cut `145`

Here `v_p(F),v_p(H)>=e`. The original `P_2` expressions for `Q_4,Q_5`
have three-bracket coefficients of order `e`; their old moving rows
`P_4,P_5` are inside the cut, whereas `P_2` is outside.
In the first Gaussian orientation the scalar-times-reference term has
order `e` and strictly dominates the first term, whose order is at least
`2e`. In the conjugate orientation both terms have order at least `e`.
Thus the moving phases `d_4,d_5` are at most zero.

Fixed `P_2,P_3` also have phase zero, producing four low phases. Fixed
`P_1` has phase `e`. Equation (1) again gives `Gamma>=e/2`.
The phase of `Q_6` is unrestricted.

## 6. Uniform transfer of all correction valuations

As in the anchored five-cut note, put

```text
E_p=sum_i v_p(N(K_i))+sum_(all15 edges) v_p(b_ij)+v_p(2).
```

Each old row orientation error and each bracket-product error is at
most `E_p`, and `sum_p E_p log p<=12s+15beta+log2=o(w)`. More precisely,
the sum of the errors in any product of three distinct brackets and
one row is at most `E_p`, because its bracket errors and row error are
charged to separate summands of this total. This local argument permits
bounded nonnegative errors at both Gaussian orientations of a row;
it does not require conjugate primitivity.
If `e>2E_p`, the strict comparisons above persist.

For `124` and `135`, the two designated moving rows have first
orientation at least `e` and conjugate orientation at most `E_p`.
The two fixed inside phases are at least `e-E_p`; the indicated outside
fixed phase is at most `E_p`. There are therefore four phases at least
`e-E_p` and one at most `E_p`.

For `145`, the dominant reference terms have first orientation at most
`e+E_p`, while the conjugate orders remain at least `e`. The designated
moving phases are at most `E_p`; the two fixed outside phases are at
most `E_p`, and the inside fixed phase is at least `e-E_p`.

Equation (1) gives `Gamma>=e/2-E_p` in all three cases. If `e<=2E_p`,
that same bound follows from `Gamma>=0`, so it is universal. Summing over
the three blocks proves a retained standalone matching-content lower
bound of

```text
(1/2)(log n_124+log n_135+log n_145)-sum_p E_p log p
 >=(3/2)(1-eta)w-(12s+15beta+log2).                      (3)
```

Higher prime powers are retained exactly. The
[exact complement covariance](joint_involution_standalone_local_symmetry.md)
sends every output signed phase `d_i` to `e-d_i` and preserves `E_p`.
Because Section 6 permits both orientation errors, it applies after this
local reflection. The same bounds therefore hold for the complementary
cuts `356`, `246`, and `236`. These six mixed-cut blocks together
contribute at least `3(1-eta)w-(12s+15beta+log2)` to the standalone
common matching content.

Together with the fixed-pair and two-pair cases, these six blocks enter
the [standalone inflation theorem](joint_involution_standalone_inflation.md).
Its sixteen distinct cuts give total common matching content
`8w-o(w)`, which rules out bounded normalized endpoint scale for every
labeling of this fixed-three joint pullback.

The root agent and fresh-algebraic agent independently checked the three
direct cases and the error accounting; the complement covariance was
independently audited in its companion note. No numerical grid is used
to infer any actual valuation in this argument.
