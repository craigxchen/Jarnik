# Additive parallelograms among the 35 balanced characters

This note enumerates every nontrivial four-term additive relation among the
35 balanced sign characters with row zero fixed positive. It tests whether
factoring the corresponding torus sum produces a leading height beyond the
ordinary pair and determinant constraints on the full 127-cut equal-weight
profile.

Let `L` be the 35 vectors `lambda in {+1,-1}^8` with four plus signs and
`lambda_0=+1`. A nontrivial relation is

```
lambda + mu = nu + rho,                                  (1)
```

where the two unordered pairs are different. There are 357 pair-sum classes:
245 classes contain one pair, 105 contain three pairs, and seven contain ten
pairs. Choosing two pairs in one class gives exactly 630 relations. Every
relation uses four distinct characters; there are no repeated-label
exceptions.

## Exact classification

Set

```
p=(lambda-mu)/2,  q=(nu-rho)/2,
r=(p+q)/2,       s=(p-q)/2.
```

Then `p,q,r,s` are integral zero-sum vectors and

```
e^(i theta.lambda)+e^(i theta.mu)
 -e^(i theta.nu)-e^(i theta.rho)
 = -4 e^(i theta.(lambda+mu)/2)
     sin(theta.r) sin(theta.s).                          (2)
```

The 630 relations split as follows.

| pair Hamming supports | `p dot q` | number | factor supports `(|r|_0,|s|_0)` |
|---|---:|---:|
| `(4,4)` | `0` | 315 | `(2,2)` |
| `(6,6)` | `+2` | 210 | `(4,2)` |
| `(6,6)` | `-2` | 105 | `(2,4)` |

Thus every additive factor is either a pair vector with two nonzero entries
`(+1,-1)`, or a balanced four-support vector with two `+1` and two `-1`.
There are no factors with a new support type.

## Per-cut heights of the two factors

For a nonconstant sign column `c`, define the primitive exponent magnitude of a
factor vector `f` by

```
m_f(c) = |f dot c|/2.
```

This is integral for both factor types. Across all 254 nonconstant sign
columns, the exact histograms are

```
factor support       m=0   m=1   m=2   total height sum_c m_f(c)
pair (2)             126   128     0               128
balanced (4)          94   128    32               192
```

Refining by minority size `t=min(#plus,#minus)` gives

```
support 2:  t=1  (12,4,0), t=2  (32,24,0),
            t=3  (52,60,0), t=4  (30,40,0)
support 4:  t=1  (8,8,0),   t=2  (20,32,4),
            t=3  (40,56,16), t=4  (26,32,12),
```

Each triple records `(m=0,m=1,m=2)`. With one independent block of weight
`h` per unoriented cut, the full profile has `W=127h`, while the two sine
factors in (2) have total leading primitive heights

```
315 relations: 128h = 64h + 64h,
315 relations: 160h = 96h + 64h.                       (3)
```

The pair factor's 64h is the ordinary row-pair height; the balanced factor's
96h is the corresponding two-versus-two scalar height. The additive identity
therefore factors entirely into already available pair/balanced directions.

## Scope of the factorwise scalar test

On an arc of angular width `Delta`, (2) gives a product of two sine factors
bounded by a constant times `Delta^2`. Whenever both associated primitive
Gaussian factors are nonunits, the usual axis gaps ask for the sum lower bound

```
log Norm(A_r) + log Norm(A_s) >= W - O_C(1).
```

The formal leading values in (3) are `128h` or `160h`, both larger than
`W=127h`, so every one of the 630 additive patterns passes this test with
positive margin. The first class is exactly two ordinary pair directions;
the other two classes add one existing balanced direction. No pattern has
a new factor support type or a leading factor height that contradicts the
full 127-cut profile. If a factor is a unit or its sine vanishes, that
nonzero axis gap is unavailable; the displayed inequality is not being
asserted in that case.

This exhausts this factorwise scalar-height test for these four-term
identities. It does not compute the complete arithmetic content of a
denominator-cleared additive numerator, or exclude additional divisibility
of its evaluated value. It also does not address higher additive identities
or arithmetic information that does not factor through the two sine
directions. The separate [generic-content theorem](additive_balanced_character_generic_content.md)
states precisely which common-denominator test is covered for higher
homogeneous identities.

[The checker](check_additive_balanced_parallelogram_audit.py) enumerates all
630 relations, verifies (2) as a vector identity, checks every one of the 254
columns, and reproduces the leading heights in (3).
