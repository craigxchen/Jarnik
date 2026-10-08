# Exact triangle compatibility and binary residue branching

The setting is the primitive fixed-unit Gaussian configuration in
[kurlberg_wigman_shrinking_arc_audit.md](kurlberg_wigman_shrinking_arc_audit.md).
Its common norm N is odd. All conclusions below keep the actual prime-power
allocations. They do not prove a new general growth rate.

## 1. The exact rational cancellation factor

For three rows, put `x_p=e_(a,p)-e_(b,p)` and
`y_p=e_(b,p)-e_(c,p)`. With the pair factors defined by oriented exponent
differences, set

```text
q_abc = product_(p:x_p*y_p<0) p^min(|x_p|,|y_p|).
```

Comparing exponents at each split Gaussian prime gives

```text
A_ab A_bc = q_abc A_ac,
d_ab d_bc = q_abc^2 d_ac,
q_abc = Norm(gcd_G(A_ab,conjugate(A_bc))).              (1)
```

The gcd is taken up to a Gaussian unit. Oppositely oriented factors cancel
as `pi_p^t conjugate(pi_p)^t=p^t`; factors with the same orientation simply
add their exponents. Thus q is a positive ordinary divisor of N and is odd.
It can differ from one. One-prime examples illustrate that algebraic fact,
but examples with uncontrolled angles prove no obstruction for endpoint
clusters.

## 2. Ordered angles give positivity

Order three source arguments `alpha_a<alpha_b<alpha_c`, with total span
less than pi/2. Conjugate the factors if necessary and multiply by signs
so their arguments are respectively the positive half-gaps. Write

```text
A_ab=r+is,    A_bc=u+iv,    A_ac=w+it,
r,s,u,v,w,t > 0.
```

These choices preserve (1) with positive q: both sides have the same
argument, and the original identity differed only by signs. Therefore

```text
ru-sv=qw,             rv+su=qt.                       (2)
```

In particular the imaginary summands have the same sign. Cancellation
between them is unavailable for an ordered triple.

## 3. Equal imaginary parts are impossible

Every pair norm is odd, so its real and imaginary coordinates have opposite
parity. If `s=v=t=B`, equation (2) gives `q=r+u`. But r and u have the
same parity, making q even, a contradiction. Thus no triangle has all
three absolute imaginary coordinates equal, without any norm-band
assumption.

An earlier intermediate calculation excluded equal B in a factor-four
norm band under an additional real-coordinate size hypothesis. It was
valid as a general integer calculation, but is unnecessary for the actual
odd-N setting: the parity contradiction is stronger.

## 4. Exact binary branching at the prime two

Put `k_ab=v_2(|Im A_ab|)`. For an ordered triple, (2) and odd q imply

```text
k_ab != k_bc  =>  k_ac=min(k_ab,k_bc),
k_ab == k_bc  =>  k_ac>k_ab.                          (3)
```

If the two valuations are positive, their real coordinates are odd.
When the valuations differ, the two summands in the imaginary equation
have different valuations, so the lower one survives. When they agree,
their sum has higher valuation. If exactly one valuation is zero, its
imaginary coordinate and the other factor's real coordinate are odd, so
the sum is odd. If both valuations are zero, both real coordinates are
even, so the sum is even. This proves all cases.

Equivalently, in any triangle the minimum of its three k-values is attained
exactly twice. For each integer t>=0, define two rows to be equivalent
when they agree or `k_ab>=t`. The ordinary ultrametric implication in (3)
makes this an equivalence relation. A class at depth t has at most two
children at depth t+1: three children would give three representatives
with every pair value exactly t, contradicting (3).

Consequently, if EVERY pair has `k_ab<=K`, there are at most `2^(K+1)`
rows. In particular,

```text
if every pair satisfies 1<=|Im A_ab|<=H,
then M<=2^(floor(log_2 H)+1)<=2H.                     (4)
```

No assertion in (4) follows merely from a lower bound on the number of
pairs with small imaginary part. The near-critical divisor lemma controls
at least `ceil(M(M-2)/4)` pairs, not necessarily all pairs. It therefore
does not justify inserting its H into (4) for the full cluster. A subset
whose every pair is controlled would satisfy (4), but extracting a large
such subset is an additional requirement.

The existing [prime-incidence analysis](residue_prime_incidence.md) and
[content hierarchy](residue_content_descent.md) already prove more general
ultrametric restrictions for primitive residues. The refinement here is the
explicit binary branching law at two in the present pair-factor coordinates.
Its collision cost is assessed in
[binary_pair_residue_capacity.md](binary_pair_residue_capacity.md).
This is a local arithmetic restriction, not the missing radius-independent
bound.

## Verification

The exact near-real divisor checker tests the triangle identities, odd
cancellation factors, and binary valuation law on a nested-prime tuple and
three actual four-point Pell clusters. General positivity is proved by the
ordered-argument normalization above. Finite checks support the identities;
they do not replace the proofs or certify a general point-count bound.
