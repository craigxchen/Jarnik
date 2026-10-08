# Reflection minimality does not automatically survive extraction

The exact reflection denominator acquires many new private factors after
rows are discarded. The one-step reflection inequalities supplied by
full-set minimality therefore do not pass to a fixed row window through the
existing uniform-profile extraction alone. There are two precise ways to
recover a useful statement: control the change of the omitted-row range
under each reflection, or require the selected valuation columns to be
binary up to scale.

## 1. Exact retained-gcd formula

Let a full primitive common-norm tuple have split-prime valuation vector
`e_p=(e_(p,1),...,e_(p,M))`. Fix a subset `I` of `k` rows containing the
deleted row `ell`, and put

```text
A=I\{ell},       O={1,...,M}\I.
```

Write `G_ell^I=gcd_G(z_r:r in A)` and
`G_ell^full=gcd_G(z_r:r!=ell)`. At one chosen Gaussian prime above `p`,

```text
v_pi(G_ell^I/G_ell^full)
  = min_(r in A)e_(p,r)-min_(r!=ell)e_(p,r),
v_barpi(G_ell^I/G_ell^full)
  = max_(r!=ell)e_(p,r)-max_(r in A)e_(p,r).            (1)
```

Consequently `log Norm(G_ell^I/G_ell^full)` is exactly the sum of the
threshold-cut weights for cuts on which all rows of `A` lie on one side and
at least one omitted row of `O` lies on the other. This is a statement about
actual nested prime layers; no replacement by independent binary primes is
made.

If every unordered full-set cut has weight exactly `w`, then

```text
log Norm(G_ell^full) = w,
log Norm(G_ell^I) = (2^(M-k+1)-1)w.                    (2)
```

There are `2^(M-k+1)-2` new cuts in (1). The common gcd `C_I` of all rows
of `I` has

```text
log Norm(C_I)=(2^(M-k)-1)w.
```

After intrinsic normalization of the selected tuple, its retained-row gcd
is therefore

```text
log Norm(G_ell^I/C_I)=2^(M-k)w.                        (3)
```

Thus the private factor relevant to the selected reflection is exponentially
larger than the original full-set private factor, even in the perfectly
uniform cut profile.

## 2. The exact commutation criterion

For any row set `J`, let

```text
R(J)=sum_p log(p) [max_(r in J)e_(p,r)-min_(r in J)e_(p,r)]
```

be the logarithm of its primitive common norm. For a reflection `T` using
only rows of `I`, define the omitted-range loss

```text
kappa_I(E)=R(full,E)-R(I,E).
```

Restriction and this row operation commute, and direct cancellation gives

```text
Delta_I(T)=Delta_full(T)
           +kappa_I(E)-kappa_I(TE).                    (4)
```

Hence `Delta_I(T)>=Delta_full(T)` precisely when
`kappa_I(TE)<=kappa_I(E)`. This is sufficient to pass a known lower bound
unchanged; a lower bound can also survive an increase in kappa if the
full-set inequality has enough additional slack. More generally, an
increase by at most `B` loses at most `B` in the inherited logarithmic
bound. The current profile
extraction controls the source cut weights, but it does not impose these
inequalities simultaneously for the transformed tuples.

This failure is visible in an actual Gaussian valuation matrix. At primes
`5,13,17`, take the five row valuations

```text
e_5  =(0,2,2,0,1),
e_13 =(2,2,0,0,2),
e_17 =(1,1,0,2,0).                                    (5)
```

Use complementary conjugate valuations `2-e_p`, so all five Gaussian
integers have the same norm and the tuple is primitive. Every
one-coordinate reflection of the full five-row tuple, including repeated
anchors, preserves or increases the primitive common norm. But on the
four-row window `I={0,1,3,4}`, replacing row 3 by `z_0 z_1/z_3` has local
width changes

```text
(Delta_5,Delta_13,Delta_17)=(0,0,-1),
N_new/N_old=1/17.                                     (6)
```

Thus even the complete collection of one-step full-row reflection
inequalities does not imply the corresponding inequalities after
restriction. This example is locally reflection-minimal; it is not claimed
to minimize the entire affine-unimodular orbit or to satisfy a short-arc
condition.

## 3. A basis-invariant sufficient criterion

Let `E` now be the valuation matrix of a selected row window. For each
prime put

```text
g_p=gcd_r(e_(p,r)-e_(p,0)).
```

If `A in GL_k(Z)` and `A 1=1`, then `A` acts invertibly on
`Z^k/Z 1`. The content `g_p` of the class of `e_p` is therefore invariant.
All coordinate differences of `Ae_p` are multiples of `g_p`, so

```text
width(Ae_p)>=g_p,
B_I=product_p p^g_p divides N_A.                       (7)
```

Define the orbit defect

```text
D_I=log(N_I/B_I)
   =sum_p (width(e_p)-g_p)log p.                        (8)
```

Every affine-unimodular monomial row transform then satisfies

```text
N_A >= N_I exp(-D_I).                                  (9)
```

In particular, `D_I=0` exactly when every nonconstant selected prime
column has two levels separated by `g_p`. Such a binary window is a genuine
least-radius representative of its entire row-basis orbit, and every
further row restriction retains that property. A bound `D_I<log 16` would
likewise give the strict constant-factor lower bound `N_A/N_I>1/16`
relevant to the short-arc minimality argument. At equality the conclusion
is only `N_A/N_I>=1/16`.

The uniform cut profile does not control (8), because it records threshold
gaps but forgets which nested gaps belong to the same prime. Independent
binary primes can realize a uniform profile with `D_I=0`; the three-level
construction in
[six_point_nested_reflection_budget.md](six_point_nested_reflection_budget.md)
realizes the same cut weights with positive defect. The existing Bessel
extraction therefore does not guarantee a window satisfying (8) with
bounded defect. Selecting such a window would be an additional arithmetic
theorem, not a consequence of full fairness.

The checker
[check_reflection_minimality_extraction.py](check_reflection_minimality_extraction.py)
verifies (2)--(3) for `3<=k<M<=8`, exhausts every full one-coordinate
reflection in (5), and checks the exact factor `1/17` in (6).
