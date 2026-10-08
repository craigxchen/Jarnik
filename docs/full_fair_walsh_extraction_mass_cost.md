# Mass cost of extracting a Walsh core from full fair cuts

Selecting rows from a full fair cut system does not cheaply leave a small
Walsh core. Restriction makes every target cut class have the same number of
source preimages. A complete `m`-row Walsh core uses only `m-1` of the
`2^(m-1)-1` nonconstant target classes, so the discarded nonconstant mass is
an exact fixed fraction of the source mass. Grouping equal restrictions and
retaining constant restrictions as common content is already optimal.

This is a combinatorial mass calculation. It does not claim that the formal
Gaussian phase assignments below occur on a short endpoint arc.

## 1. Restriction fibers

Let `X` be the `M` source rows and let every nonempty proper subset of `X`,
modulo complementation, be a source cut class. There are

```text
C_M = 2^(M-1)-1
```

classes, each with equal squared-norm logarithmic weight `w`. Select
`I subset X` with `|I|=m`, and restrict every source cut to `I`. The target
has `2^(m-1)-1` nonconstant unoriented classes.

For each target nonconstant class, the restriction fiber has exactly

```text
2^(M-m)
```

source classes. Indeed, each of its two oriented representatives has
`2^(M-m)` lifts, and complementing a lift identifies the two lists. The
constant target class has `2^(M-m)-1` source classes: the omitted endpoint
lifts are the empty and full source cuts.

This count is independent of the selected rows. It is also invariant under
any permutation or nonlinear relabeling of the selected rows, since such a
relabeling only permutes target sign patterns.

## 2. Exact optimal mass loss

Choose `m-1` target classes that form a complete Walsh core. Every other
nonconstant target class is an extra. Equal restrictions can be merged into
one core block, and constant restrictions can be retained as one common
Gaussian factor. The retained and removed weights are therefore

```text
W_keep_common = [m 2^(M-m)-1] w,
W_removed_extra = [2^(m-1)-m] 2^(M-m) w,
W_source = [2^(M-1)-1] w.                              (1)
```

The exact removed fraction is

```text
W_removed_extra/W_source
 = (2^(m-1)-m) 2^(M-m) / (2^(M-1)-1).                  (2)
```

If one also divides out the common factor to obtain a primitive core, the
additional common-content loss is `(2^(M-m)-1)w`, and

```text
W_core = (m-1) 2^(M-m) w,
W_source-W_core = [2^(M-1)-1-(m-1)2^(M-m)]w.           (3)
```

For `M` large with fixed `m`, retaining common content gives retained
fraction `m/2^(m-1)`, while the primitive core gives `(m-1)/2^(m-1)`. The
exact small-core values are

| selected rows `m` | core classes | extra classes | removed fraction, common retained | removed fraction, primitive |
|---:|---:|---:|---:|---:|
| 4 | 3 | 4 | `1/2+o(1)` | `5/8+o(1)` |
| 8 | 7 | 120 | `15/16+o(1)` | `121/128+o(1)` |

Allowing at most `k` extra target classes changes the retained coefficient
in (1) from `m` to `m+k`, provided `m-1+k<=2^(m-1)-1`. Thus a fixed
number of extras does not change the near-total loss when m grows.
For near-fair weights between `(1-eta)w` and `(1+eta)w`, every counted
mass lies between those multiples of its displayed count. In particular,
the limiting removed fractions persist as `eta` tends to zero. Additional
source correction mass `o(W_source)` cannot change those limiting
fractions either; no phase assertion follows from this observation.

For a positive-density selection `m>=delta M`, the common-retained
fraction is at most `m/2^(m-1)`, which tends to zero exponentially in `M`.
Thus the removed mass is `W_source*(1-o(1))`, not `o(W)`.
No choice of the `m-1` target Walsh classes improves this: every target
nonconstant fiber has the same size. Pairing complementary classes is already
accounted for by the unoriented quotient. Multiplying blocks from distinct
restricted classes changes their row-sign pattern and cannot make one Walsh
core block unless an additional phase or valuation relation is supplied.

## 3. Radius and phase bookkeeping

Write `W=log R^2`. If the extra blocks are removed while the constant
restrictions are retained, then

```text
W_keep = W - W_removed_extra,
R_keep = R exp(-W_removed_extra/2).                    (4)
```

For a removed block `E`, let `phi_E` be its actual Gaussian argument and
`s_E(i)` its restricted row sign. Removing it changes a lifted row argument
by

```text
eta_i = -sum_(E removed) s_E(i) phi_E.                 (5)
```

The pairwise error is `eta_i-eta_j`. A constant restriction has the same
sign on every selected row, so dividing out common content creates no
pairwise phase error. A nonconstant extra has two rows with opposite signs;
its phase contribution to such a pair is `2 phi_E`. Equal norm weights do not
control these arguments. With unrestricted block phases, the aggregate extra
error can be order one even when the removed weight is a small fraction of
`W`.

Even one literal conjugate-primitive Gaussian block can have growing
height and an order-one phase contribution: take `E=n+i(n+1)`, `n>=1`.
Its coordinates are coprime and of opposite parity, so it is coprime to
its conjugate. Its argument lies strictly between `pi/4` and `pi/2`,
while `log Norm(E)` tends to infinity. This example concerns one summand
in (5), not a whole endpoint realization or a lower bound on an aggregate
that could have cancellations.

The original endpoint precision is of order `e^(-W/4)=R^(-1/2)`. After
removal, a core theorem at the new radius would allow only order
`e^(-W_keep/4)`, which still tends to zero whenever `W_keep` tends to
infinity. Therefore a valid extraction needs an independent phase-cancellation
estimate at that exponentially small scale. Mass counting alone supplies
none. In particular, the `m=4` and `m=8` costs are already fixed fractions,
and positive-density selections lose almost all source mass before phase
control is considered.

The obstruction is about literal block removal. Restriction and exact
grouping within one fiber are harmless algebraically and preserve the source
rows. They do not reduce the number of target nonconstant classes: the other
fibers remain extras and must be removed, grouped into some other profile, or
controlled by a new phase relation.

## Verification scope

The companion [checker](check_full_fair_walsh_extraction_mass_cost.py)
exhausts selected-row subsets for `m=4` and small ambient orders and all
choices of three target classes. For `m=8`, it uses a valid Walsh core and
samples arbitrary seven-class sets. Those arbitrary sets need not be
Walsh cores; the mass formula holds for all of them. It also checks the
fiber formulas with exact integers, an order-one formal phase diagnostic,
and the primitive growing Gaussian blocks above. No endpoint realization
is asserted.
