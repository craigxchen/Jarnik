# Full-support phase inequalities bound the nonuniform class mass by flip mass

The simultaneous phase constraints are stronger than pair-Gram slack
alone. This note records an exact positive residual measure extracted
from the full-support Walsh certificates, valid for arbitrary nonlinear
assignments and arbitrary actual distinct prime weights. It does not
close the general arbitrary-weight growth problem.

Write `W0=sum_a W_a`, `F=sum_a F_a`, `D>=0`, and

```
q_a = W_a-4F_a/M,
A_m = log 8-2log(m C).
```

For a full-support nonzero Walsh character at label `a`, the exact height
is

```
V_a = (M/2)W_a+F-2F_a = (M/2)q_a+F.
```

The actual phase lower bound therefore gives, for every nonzero label,

```
q_a >= l := (W0+D-2F+2A_M)/M.                            (1)
```

Thus the measure `r_a=q_a-l` is nonnegative. Its total is exactly

```
R = sum_(a!=0) r_a
  = W0/M+(2-6/M)F-(1-1/M)(D+2A_M).                       (2)
```

In particular, at fixed `C`,

```
R <= 2F+W0/M+O_C(log M).                                 (3)
```

This is a decomposition of class mass into a uniform signed background
`l` and a **positive** residual whose total is controlled by the flipped
prime mass. The background itself need not be positive. For `M>=4`,
`F_a<=W_a` implies `q_a>=0`; if the background is negative, this
supplies the additional lower bound `r_a>=-l`.

## Exact residual formulation for every coset character

For any row subspace `V` of size `h>=2`, nonzero `ell in V*`, and row
coset `P`, let `Rcell` and `Fcell` be the residual and flip-label masses
in the label fiber `a|V=ell`. Let `Fp` be the flip mass on `P` and
`Fmatch` the matching part of it. The fiber contains exactly `M/h`
nonzero labels. Substitution of `W_a=l+r_a+4F_a/M` into the exact height
and the phase lower bound gives

```
(h/2)Rcell+(2h/M)Fcell+Fp-2Fmatch
  >= F+2log(M/h).                                       (4)
```

All `C,D,W0` terms cancel from this formula after they have entered
`l` and `r`. In particular an aligned flat must satisfy

```
h Rcell >= 2F+2Fp-(4h/M)Fcell+4log(M/h).                 (5)
```

Averaging (4) over **all** `M/h` row cosets, without any alignment or
assignment assumption, uses

```
E Fp=hF/M,        E Fmatch=hFcell/M.
```

The two `Fcell` terms cancel and leave the further exact condition

```
Rcell >= 2F(1/h-1/M)+(4/h)log(M/h).                      (6)
```

Equations (2), (4), and (6) apply simultaneously to every nonzero
restriction fiber and every row subspace. They offer a positive-measure
formulation of part of the global certificate problem; a single fixed
aligned family does not enforce them.

For example, a proposed two-level label-weight construction can satisfy
all pair-Gram inequalities while violating (1) at a light label. In the
nominal balanced case with nine physical copies, half of the label
weights near `9 tau`, half near `18 tau`, and flip total near
`(3/2)M tau`, a light full-support character has

```
V_a-W0/2 = -(3/4)M tau+O(tau).
```

It is therefore excluded by actual full-support phase inequalities for
large `M`. This observation does not say that every assignment or every
variant of the fixed-family obstruction violates the same certificate.

The companion `check_walsh_full_support_residual_mass.py` checks 919
coset certificates at `M=8,16`, including direct physical-column height
calculations, the full-support identity, residual total, and all-coset
averaging. Its arbitrary weighted fixtures need not satisfy the phase
inequalities; the checked statements are the exact algebraic identities.

Applying these inequalities to all affine label flats gives the stronger
[average class-weight regularity theorem](walsh_phase_class_mass_regularity.md).
The later [four-row Ramsey theorem](walsh_arbitrary_weight_ramsey_growth.md)
uses the same simultaneous constraints to obtain arbitrary-weight growth
without an additional affinity or slack-background hypothesis.
