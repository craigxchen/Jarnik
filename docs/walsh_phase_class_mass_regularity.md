# All coset phase constraints force average class-weight regularity

In the repeated Walsh model, actual endpoint phases force the physical
prime-log masses of the nonzero labels to become nearly equal in total
variation. No affinity, balance, or nearly equal prime-size hypothesis is
needed. This is a necessary regularity theorem, not yet a growth theorem
or a solution of the uniform endpoint problem.

Use the model of the [full-support residual lemma](walsh_full_support_residual_mass.md):
`M=2^t>=16`, at least five physical copies of each nonzero Walsh label,
one distinct physical-column flip per row, and arbitrary actual distinct
split primes. Let `W_a` be the total log-prime weight of label a,
`W0=sum_a W_a`, `F_a` its flipped weight, `F=sum_a F_a<=W0`, and
`D>=0` the log norm of the full common Gaussian content. Assume the rows
lie on an arc of length at most `C sqrt(R)`.

Then

```text
sum_(a!=0) |W_a-W0/(M-1)|
 <= 40F/log M + 2W0/M + 8log M + 8log^+ C.                (1)
```

In particular the left side is `O_C(W0/log M)`. All logarithms below
are natural except where a binary logarithm is explicitly written.

## Density forces a contained affine flat

Let `A` be any subset of `F_2^t`, of density alpha, and let `K=2^d<=M`.
A sufficient condition for A to contain an affine flat of size K is

```text
alpha >= 2(K/M)^(1/K).                                   (2)
```

The assertion is vacuous when the right side exceeds one. For K=1
it follows just from nonemptiness. For larger K, start with the selected
singletons A. Suppose there are n selected cosets of one subspace V of
size h, all contained in A, with union density `rho=nh/M`. Choose
their representatives in a linear complement of V. If `n>=2`, there
are at least `n^2/2` ordered distinct pairs. Among fewer than `M/h`
directions, one direction occurs at least `n^2 h/(2M)` times. Each
pair in that direction has a unique partner, so after dividing by two
for reversal these pairs give disjoint cosets of a subspace of size 2h.
Their new union density is at least `rho^2/2`.

After i successful steps, with `h=2^i`, the retained density is at least
`2(alpha/2)^h`. Condition (2) implies

```text
(M/h)(alpha/2)^h >= (M/K)(alpha/2)^K >= 1    for h<=K.
```

Thus the number of selected cosets is at least two at every required
extension step. The iteration reaches size K. This elementary density
argument has no label-map or prime-weight hypotheses.

## No affine flat can have a large total negative residual

Write `A_M=log 8-2log(MC)`. The full-support phase inequalities give

```text
q_a=W_a-4F_a/M,
l=(W0+D-2F+2A_M)/M,
r_a=q_a-l>=0.
```

Put `s0=2F/M>0` and `mu_a=r_a-s0`, so `mu_a>=-s0`.
Every affine label flat `P` that avoids zero is a nonzero restriction
fiber: if its direction is U, use the row subspace `V=U^perp`.
Its size is `K=M/|V|`. The averaged residual inequality in the cited
lemma therefore gives

```text
sum_(a in P) mu_a >= -s0+(4K/M)log K >= -s0.              (3)
```

This applies to every such affine flat, including singletons. It is
stronger than merely knowing that the residual is nonnegative.

For `0<eta<=1`, let

```text
A_eta={a!=0:mu_a < -eta s0},
alpha_eta=|A_eta|/M,
K=2^ceil(log_2(1/eta)).
```

If a K-flat were contained in `A_eta`, its residual sum would be
strictly less than `-K eta s0<=-s0`, contradicting (3). Such a flat
automatically avoids zero, since `A_eta` does. Applying (2), for
`eta>=2/sqrt(M)` we have `K<=2/eta<=sqrt(M)` and obtain

```text
alpha_eta < 2(K/M)^(1/K)
          <= 2 exp[-eta log M/4].                        (4)
```

The case eta=1 has empty `A_eta` and satisfies the same bound. For
smaller eta use just `alpha_eta<=1`.

Integrating these level-set counts bounds the full negative mass:

```text
N_minus := sum_(a!=0) max(0,-mu_a)
 = M s0 integral_0^1 alpha_eta d eta
 <= M s0 [2/sqrt(M)+8/log M]
 = 4F/sqrt(M)+16F/log M
 <= 20F/log M.                                         (5)
```

The last step uses `log M<=sqrt(M)` for `M>=16`. In particular this
controls all negative residual levels at once; choosing one fixed
threshold eta would not give (5).

## Returning to actual physical class weights

Let `wbar=W0/(M-1)` and

```text
qstar=l+s0=(W0+D+2A_M)/M.
```

Then `W_a=qstar+mu_a+4F_a/M`. The last term is nonnegative, so

```text
sum_(a!=0) max(0,wbar-W_a)
 <= N_minus+(M-1)max(0,wbar-qstar)
 <= N_minus+W0/M+4log M+4log^+ C.                        (6)
```

Here `D>=0` and `-2A_M=4log(MC)-log 64` give the second inequality.
Since the deviations `W_a-wbar` sum to zero, their absolute sum is
twice their negative mass. Equations (5)--(6) prove (1).

Finally the distinct physical primes imply
`W0>=log((b(M-1)+1)!)`, so the last three terms in (1), as well as
the term involving F, are `O_C(W0/log M)`. No bound on common content
was assumed, and it kept its favorable sign throughout.

The proof also gives an unconditional restriction on that common content.
Summing mu and using (5) gives

```text
sum mu=(W0-4F)/M-(1-1/M)(D+2A_M) >= -20F/log M.
```

Consequently

```text
D <= (W0-4F)/(M-1)+20MF/((M-1)log M)+4log(MC)-log 64
   = O_C(W0/log M).                                     (7)
```

Thus `D/W0` tends to zero as M tends to infinity in actual configurations
of this model. This is a consequence of all the phase inequalities;
it is not a primitivity hypothesis.

This common-content corollary is redundant: the previously established
pair-slack identity `U=Q+(M-1)kappa>=0`, with `Q=W0-4F/M`, already gives
the stronger bound
`D<=4log C-log 4+Q/(M-1)<=4log C-log 4+W0/(M-1)`.
The new conclusion of this note is the class-weight regularity in (1).

## Remaining limitation

Total-variation regularity of the class weights does not give the
pointwise signed-slack control needed by the current nonlinear aligned
flat. In particular it does not assert that
`Psi=sum_(d in V minus {0})(-1)^(ell(d)+1)u_d`
is `O_C(W0/M+1)` on a selected flat. Sparse deviations can still affect
a certificate much more than their average suggests. The new theorem
is a necessary consequence of simultaneous actual phases, and the
uniform endpoint count remains unproved.

The [density checker](check_walsh_phase_class_mass_regularity.py) verifies
the exact coset-merging recurrence, the sufficient-density implication
on finite fixtures, and the rational residual-to-class bookkeeping.
The level-set integration and the asymptotic bound are proved above.

The subsequent [four-row and Ramsey argument](walsh_arbitrary_weight_ramsey_growth.md)
improves (1) to `O_C(W0/sqrt M)` by Parseval and proves an unconditional
growth improvement within the repeated-Walsh model. The density proof
above remains valid but is no longer the strongest regularity estimate.
