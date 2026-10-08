# Pair-Gram slack does not control one prescribed aligned family

The Fourier slack identity below is exact. It also has a concrete
limitation: after an aligned family with common `V,ell` has been chosen,
actual distinct split-prime weights of total size `O(M log M)` can make
**every** certificate in that family exceed primitive half-height while
all individual pair-Gram phase inequalities hold. This applies to the
dense family retained by the nonlinear extraction theorem.

The construction does not settle weight-adaptive extraction or the
simultaneous constraints from all coset characters. Its Gaussian rows
are not asserted to lie on an endpoint arc.

## Exact identity and necessary phase inequality

Use the notation of the [three-family weighted theorem](walsh_three_family_weighted_affine_growth.md):
`W_a` is total physical weight at label `a`, `F_a` its total flipped
weight, `W0=sum W_a`, `Ftot=sum F_a`, and

```text
kappa = 4 log C - log 4 - D,
q_a = W_a - 4 F_a/M,
u_d = kappa - qhat(d) >= 0       (d != 0).
```

Fix any subspace `V` of size `h>=2`, nonzero `ell in V*`, and coset
`P=x0+V`. Let `Fcell` sum `F_a` over labels with `a|V=ell`;
let `Fp` be total flipped weight at rows in `P`; and let `Fmatch`
be the part of `Fp` whose assigned labels restrict to `ell`. Put

```text
Psi = sum_(d in V minus {0}) (-1)^(ell(d)+1) u_d.
```

Writing `Q=sum q_a`, character orthogonality gives

```text
h sum_(a|V=ell) q_a = Q-kappa+Psi,
```

since the nonzero character has sum `-1` on `V minus {0}`.
The exact signed coset-character height is therefore

```text
V_Beta = W0/2-kappa/2+Psi/2
         +2h Fcell/M-2 Ftot/M+Fp-2 Fmatch.                 (1)
```

The usual actual-phase lower bound is
`V_Beta>= (W0+D)/2+log 8-2log(hC)`. Substitution into (1) gives

```text
Psi >= 2(2 Fmatch-Fp)+(4/M)(Ftot-h Fcell)
       -4log(h/2).                                       (2)
```

Both `D` and `C` cancel outside the definition of `u_d`. For an aligned
coset `Fmatch=Fp`. Nonnegativity of the individual `u_d` does not
provide an upper bound for this signed sum.

There is a useful sufficient target for a future weighted extraction.
Since `Fcell<=Wcell` and `Qcell=Wcell-4Fcell/M`, when `M>4` we have

```text
h Fcell <= M/(M-4) [W0-4 Ftot/M-kappa+Psi].
```

For an aligned coset, (2) consequently implies

```text
M Psi+4W0-4kappa-4Ftot
  >= (M-4)[2Fp-4log(h/2)].                               (2a)
```

The weaker inequality obtained by dropping `-4Ftot` on the left is
also valid. Thus extraction of a growing aligned flat with
`Psi=O(W0/M+|kappa|)` and large `Fp` would control the radius. The
problem is obtaining that upper bound on `Psi` for an appropriate
family or flat.

There is no hidden common-content loss in this sufficient target.
The identity `sum_(d!=0)u_d=Q+(M-1)kappa>=0` gives
`kappa>=-W0/(M-1)`, while its definition gives
`kappa<=4log C-log 4`. Hence
`|kappa|<=W0/(M-1)+|4log C-log 4|`. The target above is therefore
`Psi=O_C(W0/M+1)` at fixed implied constant. If a selected aligned
flat also has `Fp>=c h log M`, with `h` growing, (2a) would imply
`W0>=c' M h log M`. Obtaining both properties for one flat is unproved.

## Restricting to any chosen row subspace costs few patched rows

Suppose every assigned global label is nonzero and each label occurs
at most `b` times, as follows from distinct assigned physical columns
with `b` copies per label. Let `H` be **any** row subspace of size
`N=2^s`; it may be chosen using the weights. Exactly `M/N-1` nonzero
global labels annihilate `H`. Consequently

```text
#{x: a(x)|H=0} <= b(M/N-1).
```

Averaging over the `M/N` cosets of `H` gives a coset `T=x0+H` with
at most `b(1-N/M)<=b` such rows. This count concerns actual labels
at actual rows in `T`; no restrictions are imposed on other cosets.

Identify `H` with `F_2^s` and restrict each `a(x0+y)` to `H`.
Replace the zero restrictions on those at most `b` rows by arbitrary
nonzero functionals on `H`. For `s>=4096`, apply the
[nonlinear extraction theorem](walsh_arbitrary_nonlinear_aligned_flat_extraction.md)
to this patched vector field. It produces a common `V<=H`, nonzero
`ell in V*`, and disjoint aligned `V`-cosets with

```text
h=|V| >= s/[16(log_2(s+2))^2],
total retained rows >= N^(3/4).
```

Discard every retained coset meeting a patched row. The retained cosets
are disjoint, so at most `b` of them are discarded, losing at most
`bh` rows. The unpatched family thus contains at least

```text
N^(3/4)-bh                                               (2b)
```

rows. If this is positive, at least one aligned flat survives. More
generally, at any dimension `s` where extraction provides a family
of `L` retained rows, the same argument retains at least `L-bh` rows.

Every surviving flat becomes the **actual global row coset**
`x0+y0+V`. Since none of its rows was patched, its actual global
assigned labels restrict to `ell` on `V`. Its coefficient fiber,
`Fcell`, slack `Psi`, full source height, and phase inequalities
(1)--(2a) are all computed in the original `M`-row model. One must
not replace the global `M` or the full physical prime weight by `N`.

This permits weight-adaptive selection of `H` without assuming the
restricted labels are everywhere nonzero. It does not itself supply
the required slack control. In particular, the obstruction below
concerns a fixed family; it does not preclude finding a useful family
inside an appropriately chosen different subspace.

Any proposed background regularity must allow an error for actual primes.
For `M>4`, the values `u_d` at distinct nonzero directions are all
distinct: equality would give a rational linear relation among the
logarithms of the distinct physical rational primes. A prime's
coefficient in `qhat(d)-qhat(e)` is its positive factor `1` or
`1-4/M`, multiplied by the difference of its two Walsh signs.
At least one label distinguishes `d` and `e`. Unique factorization
therefore forbids that relation. In particular, an exactly constant
slack background outside a small set of directions is not a valid
actual-prime hypothesis here.

## A halfspace of heavier labels shields a prescribed family

Fix the common nonzero `V,ell` first and choose `w in V` with
`ell(w)=1`. Take `b=9` copies of every nonzero Walsh label. Initially
give every physical copy at label `a` the nominal weight

```text
tau                    if a dot w=0,
sigma                  if a dot w=1,
1.9 tau <= sigma <= 2.1 tau.                              (3)
```

Any one-flip-per-row assignment using distinct physical columns is
allowed; the following bounds do not assume affinity or alignment.
The nominal total physical weight is

```text
W0 = b[(M/2-1)tau+(M/2)sigma].
```

Every label in the coefficient fiber `a|V=ell` lies in the heavier
halfspace. Its baseline height is `bM sigma/2`. Each selected flip
can lower height by at most `sigma`. Consequently, for **every coset**
of this `V` and its fixed character `ell`,

```text
V_Beta-W0/2 >= bM(sigma-tau)/4+b tau/2-h sigma
             >= M[(b/4-1)sigma-b tau/4]+b tau/2
             >= M tau/8+b tau/2.                         (4)
```

The bound holds for every `h<=M`, and hence for an arbitrarily dense
retained family. Alignment only fixes the exact flip correction.

At the same time, for every nonzero direction `d`, the unflipped
weighted Gram is exactly

```text
-b tau - (bM/2)(sigma-tau) 1_(d=w).
```

Two distinct rows affect only their two assigned columns. The Gram
can increase by at most `4 sigma`, so every off-diagonal entry obeys

```text
G_xy <= -b tau+4 sigma <= -0.6 tau.                       (5)
```

For large `tau` this meets the actual pair-phase necessary condition
with `C=1,D=0`, namely `G_xy<=-log 4`. Thus it also gives nonnegative
Fourier slack `u_d`.

For the nominal weights, its exact slack is

```text
u_d = kappa+b tau+(bM/2)(sigma-tau)1_(d=w)
      +(4/M) Fhat(d).
```

The single positive spike at `w` contributes positively to `Psi` for
every member of the prescribed family. The identity (1) does not
remove this large contribution.

## Actual distinct split primes at the same radius scale

This is not limited to formal real weights. Set `X=M^4`, and partition
each of `[X,2X]` and `[X^2,2X^2]` into `r=9(M-1)` equal intervals.
The fixed-modulus prime number theorem, followed by the same
pigeonhole argument as in the [nearby-prime construction](strict_obtuse_prime_box_countermodels.md),
gives at least `r` primes congruent to one modulo four in one interval
of each partition, for all sufficiently large `M`.

If the left endpoints are `A,B`, take `tau=log A` and `sigma=log B`.
The actual light logs lie in `[tau,tau+1/r)`, and the heavy logs in
`[sigma,sigma+1/r)`. The two intervals are disjoint, so all required
physical prime norms can be chosen distinct. Also
`tau=4log M+O(1)`, `sigma=8log M+O(1)`, giving (3) eventually.

The perturbation of each Gram entry has absolute value less than one.
For a coset character, every physical half-coefficient has absolute
value at most `h/2`, so the perturbation of `V_Beta-W0/2` has absolute
value less than `(h+1)/2`. Thus the strict positive gap in (4) and the
pair inequality (5) persist as `M` grows. The full physical weight is
`W0=O(Mlog M)`. Arbitrary Gaussian-prime orientations give literal
primitive Gaussian rows realizing these norms and heights.

This rules out closing the arbitrary-weight gap using only pair-Gram
constraints and a **single family selected without reference to the
weights**. A proof may still use other restrictions, other subspaces,
a weight-adaptive selection, or the simultaneous actual phase equations.

The [checker](check_walsh_pair_slack_fixed_family_obstruction.py) verifies
(1) exactly on rational weights and on arbitrary coset characters,
checks the prescribed-family obstruction and every pair Gram, and
repeats the height and Gram checks with actual distinct split primes.
It also checks restriction-zero counts, coset averaging, and the
patched-coset discard bound in finite examples; the large-dimensional
extraction assertion follows from the cited theorem.
