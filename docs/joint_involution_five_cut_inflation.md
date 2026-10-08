# Five retained cuts force radius inflation for the anchored joint map

Keep the original anchor and apply the arc-preserving ordered joint
involution to six other points from a full sixty-four-cut extraction.
Then its least integral output circle satisfies

```text
log R_new >= (193/6)w-o(w),
R_new >= R_old^(193/192-o(1)).                            (1)
```

Its normalized arc length is at least `exp(w/12-o(w))`. Thus this
anchored seven-point operation cannot give an endpoint-scale radius
descent. Retaining the anchor is essential: this is not a statement
about the least circle of the six transformed points alone. The later
[standalone theorem](joint_involution_standalone_inflation.md) proves a
stronger obstruction for those six points by using sixteen other cut
bounds and exact complementary-cut symmetry.

## 1. The raw Gaussian identities and common matching factor

Write `P_i=K_i product_(S containing i) H_S` for the six primitive
half-angle numerators, with core block norms `n_S`, and retain the
original anchor of relative phase one. The core weights satisfy
`(1-eta)w<=log n_S<=(1+eta)w`. The selected corrections and residues
are negligible compared with `w`.

Let `Q_4,Q_5,Q_6` be the exact degree-seven vectors from
[joint_circle_conductor_lift.md](joint_circle_conductor_lift.md).
Both useful forms of the identities are

```text
Q_4=H P_4-2 Delta_34 Delta_56 Delta_41 P_2
   =V P_4-2 Delta_34 Delta_56 Delta_42 P_1,
Q_5=F P_5-2 Delta_35 Delta_46 Delta_51 P_2
   =V P_5-2 Delta_35 Delta_46 Delta_52 P_1,
Q_6=H P_6+2 Delta_36 Delta_45 Delta_61 P_2
   =F P_6+2 Delta_36 Delta_45 Delta_62 P_1.                (2)
```

The ordered real-arc construction retains the first three rows, moves the
other three inside the old arc, and gives `|Q_i|=exp(64w+o(w))`.
The matching determinants of `W=(P_1,P_2,P_3,Q_4,Q_5,Q_6)` equal `T(M)`
literally. If `gamma` is the common Gaussian factor of the fifteen
three-chord matching products on the least output circle, then

```text
log|gamma|=3log R_new-240w+log gcd(T(M))+o(w),
log gcd(T(M))<=150w+o(w).                                (3)
```

This is an exact scalar identity before the displayed asymptotic estimates;
it does not discard rational contents of the Gaussian vectors.

## 2. A local content formula

At an odd split prime in a core block, let its core exponent be `e`.
For each of the six output points let `d_i` be the signed Gaussian phase
valuation relative to the retained anchor. Put the six `d_i` in increasing
order and set

```text
n=max(0,d_1,...,d_6)-min(0,d_1,...,d_6).
```

The old prime contributes at least

```text
(1/2)(3n+sum_(i<=3)d_(i)-sum_(i>=4)d_(i)) log p          (4)
```

to `log|gamma|`. This follows by pairing the three smallest and the three
largest valuations in a matching, separately at the two Gaussian primes.
The complete unconditional proof is in
[joint_circle_matching_content.md](joint_circle_matching_content.md).

Two elementary consequences will be used. If all six `d_i>=e`, then
`n=max d_i` and (4) is at least `(3/2)e log p`. If the first three
phases are `(e,e,0)` and the other three are nonpositive, (4) is at least
`(1/2)e log p`.

## 3. Five directly forced old cuts

First ignore the small correction valuations. The old phases at a cut
`S` are `e` inside and zero outside. The following deductions come from
strictly unequal valuations in the two summands in (2); no genericity
assumption or finite-grid argument is involved.

| Cut `S` | Form valuations needed | Forced moving phases |
|---|---|---|
| `123` | `F,H,V>=e` | `d_4,d_5,d_6>=e` |
| `1234` | `H=V=e`, `F>=e` | `d_4=e`, `d_5,d_6>=e` |
| `1235` | `V=e`, `F,H>=e` | `d_5=e`, `d_4,d_6>=e` |
| `1236` | `F=H=e`, `V>=e` | `d_6=e`, `d_4,d_5>=e` |
| `12` | `F=H=V=0` | `d_4,d_5,d_6<=0` |

Here a form valuation such as `V=e` means `v_p(V)=e`, not equality of
the form and the integer exponent. The indicated forms are ordinary
integers, so their two Gaussian valuations agree.

For `123`, all the relevant forms have valuation at least `e`, while
each triple-bracket coefficient in (2) has valuation zero. The reference
vectors `P_1,P_2` have Gaussian orders `(e,0)`. Thus each `Q_i` has
first orientation at least `e` and conjugate orientation zero.

For `1234`, the matching orders are
`(v_p(A),v_p(B),v_p(C))=(2e,e,2e)`. Consequently
`H=B+C` and `V=2A+B+C` have valuation `e`. In the expression for `Q_4`,
its first summand has Gaussian orders `(2e,e)`, while the other has
orders `(3e,2e)`. The first dominates at both orientations. For rows
five and six the previous `123` reasoning still applies.

For `1235`, use

```text
V=(A+C)+(A+B),
A+C=Delta_13 Delta_24 Delta_56,
A+B=Delta_12 Delta_35 Delta_46.
```

These two matching products have orders `e,2e`, respectively. The
`P_1` expression for `Q_5` in (2) therefore has the same strict
dominance as `Q_4` in the preceding case. For `1236`, the orders of
`C,B` are `e,2e`, so `F=C-B` and `H=C+B` have order `e`; use the
expression for `Q_6` directly. The other moving rows again have
first orientation at least `e` and conjugate orientation zero.

For `12`, the orders of `A,B,C` are `e,e,0`, so `F,H,V` are units.
Each moving old `P_i` is also a unit in the first Gaussian orientation,
whereas the `P_1` or `P_2` summand has positive valuation. Hence every
`Q_i` has first orientation zero and conjugate orientation nonnegative.

The four cuts containing `123` therefore each contribute at least
`(3/2)log n_S` to `log|gamma|`. The cut `12` contributes at least
`(1/2)log n_12`. Their total is `6.5w` in the equal-weight model.

## 4. Correction errors and arbitrary prime powers

Let the actual bracket factorization be
`|Delta_ij|=c_ij b_ij`, with `log b_ij<=beta`. At a rational prime put

```text
E_p=sum_(i=1..6) v_p(N(K_i))
       +sum_(all15 edges) v_p(b_ij)+v_p(2).              (5)
```

Each Gaussian orientation correction of an old row, and the extra
valuation of any matching product, is at most `E_p`. Suppose first
that `e>2E_p`. Every strict comparison in Section 3 persists.

At the four cuts containing `123`, each unaffected moving row has first
orientation at least `e` and conjugate orientation at most `E_p`.
The designated inside moving row, when present, is still dominated by
the same scalar multiple of its old `P_i` at both orientations, so its
signed phase is at least `e`. Thus all six output phases are at least
`e-E_p`. Formula (4) is at least `(3/2)(e-E_p)log p`.

At `12`, the units of the idealized calculation have valuations at most
`E_p`. The first term of each `Q_i` has first-orientation valuation at
most `2E_p`, while its other term has valuation at least `e`. Thus each
moving signed phase is at most `2E_p`. The first two old phases are at
least `e` and the third is at most `E_p`. In the sorted formula (4),
the top two are the first two old phases. If `m` is the least phase or
zero, the bottom three sum is at least `3m`, so (4) is at least
`(e/2-E_p)log p`.

If `e<=2E_p`, omit the prime's nonnegative contribution and charge
at most `3E_p log p` for the idealized contribution. This proves,
uniformly at all primes in the five selected blocks,

```text
log|gamma| >= (13/2)(1-eta)w-3sum_p E_p log p.            (6)
```

If `log|K_i|<=s`, then
`sum_p E_p log p<=12s+15beta+log2=o(w)`.
All higher prime powers have been retained; primes with a large relative
correction were charged to the same negligible total, rather than
silently discarded. In particular,

```text
log|gamma| >= (13/2)w-o(w).                              (7)
```

## 5. Radius inflation and exact scope

Substitute (7) and the upper common-content bound in (3):

```text
3log R_new >= (240-150+13/2)w-o(w)=(193/2)w-o(w).
```

This proves (1). The inherited input scale is `log R_old=32w+o(w)`.
The retained first three directions have separation
`exp(-16w+o(w))`, so the normalized output arc length satisfies

```text
Theta_out sqrt(R_new) >= exp(w/12-o(w)).                 (8)
```

Both statements apply to the actual least anchored output circle; they
do not assume beforehand that this output is at endpoint scale.

The anchor enters through `n` in (4). Removing it can decrease the local
range and remove common allocations of all six transformed points.
Therefore the five-cut argument alone does not settle the standalone
six-point map. That case is settled by the later
[sixteen-cut theorem](joint_involution_standalone_inflation.md), which
proves `R_new>=R_old^(49/48-o(1))` for the arc-preserving ordering and
normalized arc length at least `exp(w/3-o(w))` for every labeling of the
fixed-three pullback. Other joint rational transformations, subsequent
arbitrary fractional-linear changes, and endpoint uniformity remain
outside these results.

## Independent checks

The root agent and uniformity-audit agent independently checked the five
idealized cut cases; the latter explicitly corrected the earlier failed
cross-cut heuristic before this argument. The fresh-algebraic agent
checked the raw Gaussian identities, their sizes under arc containment,
and the Gaussian matching scalar. The proof uses only the five direct
valuation comparisons above; no experimental all-cut minimum is used.
