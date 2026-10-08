# The joint involution loses endpoint scale even without the old anchor

The arithmetic audit of this particular joint symmetry is now complete.
For its arc-preserving ordering, every integral realization of the six
transformed directions, allowing a common rotation and omitting the old
anchor, requires

```text
log R_out >= (98/3)w-o(w),
R_out >= R_in^(49/48-o(1)),       log R_in=32w+o(w).       (1)
```

Its normalized arc length is at least `exp(w/3-o(w))`. In fact the
normalized-arc obstruction holds for every labeling of this fixed-three
pullback, whether or not that labeling preserves the old arc. Thus
dropping the old anchor does not make this construction an endpoint
descent. This result does not prove the desired uniform lattice-point
bound or exclude other transformations.

## 1. Actual inputs and the existing height estimate

Take six primitive Gaussian half-angle numerators `P_i` from a full
sixty-four-cut extraction. Their independent core block norms satisfy

```text
(1-eta)w <= log n_S <= (1+eta)w,
log|K_i| <= s,       log b_ij <= beta,
eta -> 0,           s+beta=o(w).
```

Here `P_i=K_i product_(S containing i)H_S`, `n_S=N(H_S)`, and
`|Delta_ij|=b_ij product_(S containing i,j)n_S`. Use the actual
central correction bound `beta=2s+log T`; the argument also works
with the weaker `4s+log T` when only that bound is available.

The rational joint involution fixes the first three directions and
replaces the other three by the degree-seven Gaussian vectors

```text
Q_4=H P_4-2 Delta_34 Delta_56 Delta_41 P_2,
Q_5=F P_5-2 Delta_35 Delta_46 Delta_51 P_2,
Q_6=H P_6+2 Delta_36 Delta_45 Delta_61 P_2.
```

The previously audited matching map and its six exceptional contents
give the actual projective height lower bound

```text
h_out >= 42w-342eta w-240beta-16log2.                    (2)
```

This is [the exact height theorem](joint_segre_involution.md); it does
not assume the unproved generic equality `h_out=48w+o(w)`. It holds
for every labeling of the full profile. No hypothesis on the output
residue sizes is used.

## 2. Sixteen cuts retain actual common matching content

Let `gamma` be the Gaussian gcd of the fifteen products of three
matching chords on an integral circle realizing these six directions.
For sorted signed phase valuations `d_(1)<=...<=d_(6)`, the primewise
contribution to `log|gamma|` is at least `Gamma(d) log p`, where

```text
Gamma(d)=(3(d_(6)-d_(1))+d_(1)+d_(2)+d_(3)
                             -d_(4)-d_(5)-d_(6))/2.
```

The range contains only the six output rows. No zero-phase anchor is
inserted. Equivalently, for the five successive gaps,

```text
Gamma=g_1+g_2/2+g_4/2+g_5.                              (3)
```

The following eight cuts each satisfy `Gamma>=e/2-2E_p`, with
arbitrary prime powers and cancellation depths:

| Cuts | Proof |
|---|---|
| `12,13,23` | [Fixed-pair dominance](joint_involution_standalone_local_symmetry.md) |
| `124,135,145` | [Four high or four low phases](joint_involution_mixed_cut_content.md) |
| `16,25` | [Continuous two-pair argument](joint_involution_two_pair_cuts.md) |

The local error budget is

```text
E_p=sum_i v_p(N(K_i))+sum_(all15 edges)v_p(b_ij)+v_p(2).
```

These proofs allow nonnegative errors in both Gaussian orientations,
including nonprimitive intermediate local rows. The exact reflection

```text
P_i -> pi^e bar(P_i)/p^(e 1_(i in S))
```

complements the old cut, preserves `E_p`, and sends every output phase
to `e-d_i`. By (3), it preserves `Gamma` exactly. Thus the same bounds
hold for the eight complementary cuts

```text
3456,2456,1456,356,246,236,2345,1346.
```

All sixteen cuts are distinct. A core rational prime belongs to only
one block, so the error sums once across these blocks. Contributions
at all other primes are nonnegative and may be discarded. Therefore

```text
log|gamma| >= 8(1-eta)w-2(12s+15beta+log2)
           = 8w-o(w).                                  (4)
```

This is an actual common-content bound. Its proof uses the displayed
dominance and reflection identities, not minimum orders from a generic
polynomial substitution or an experimental integer valuation grid.

## 3. Circle radius and angular scale

On any integral circle, the matching chord products have the same
complex phase up to sign. Their primitive coordinate vector is
integral, and its integer Bezout combination makes `gamma` a nonzero
Gaussian integer. If the angular diameter is `Theta_out`, every chord
has modulus at most `R_out Theta_out`. Consequently

```text
h_out+log|gamma| <= 3log(R_out Theta_out).                (5)
```

For the proved arc-preserving ordering,
`Theta_out<=exp(-16w+o(w))`. Combining (2), (4), and (5) yields

```text
3log R_out >= (42+8+48)w-o(w)=98w-o(w),
```

which proves (1). With canonical `w=W_original/64` and input angular
bound `Theta_in<=C exp(-16w)`, the same argument gives the finite form

```text
log R_out >= (98/3)w-(350/3)eta w-8s-90beta
                           -6log2-log C.                (6)
```

Common rotation does not change any signed phase difference or chord
length. The old anchor need not be on the new circle. If one writes
`z_i=lambda W_i/bar(W_i)`, the common factor `lambda` can lie in
`Q(i)` instead of `Z[i]`; (5) still follows from the actual integral
chords. Thus no integrality of a removed anchor is being assumed.

For any labeling, the retained first three directions give
`Theta_out>=exp(-16w+o(w))`. Equations (2), (4), and (5) also give
`log(R_out Theta_out)>=(50/3)w-o(w)`. Hence, even without the
arc-preserving upper bound,

```text
2log(Theta_out sqrt(R_out))
 =log(R_out Theta_out)+log Theta_out
 >=(2/3)w-o(w).
```

This proves the normalized-arc lower bound `exp(w/3-o(w))` for every
labeling of the fixed-three pullback and every common rotation of it.

## Scope

The theorem concerns the actual six directions obtained by the specified
joint involution and its inverse normalization fixing three input rows.
It also applies when the old anchor is retained, since retaining another
point cannot lower the least necessary radius. It does not address a
subsequent arbitrary fractional-linear change of the output directions,
other joint maps, or a different construction acting on a whole cluster.
The general uniform-bound goal remains open in this project.

The raw polynomial identities, projective height estimate, arc ordering,
local inequalities, error transfer, and complementary-cut symmetry were
independently checked across the root and three research agents. The
supporting exact certificates are linked from the component notes.
