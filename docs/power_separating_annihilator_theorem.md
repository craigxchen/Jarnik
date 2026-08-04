# Power-separating annihilators: theorem and audit

## Status

This note corrects the earlier proposed low-width annihilator theorem.

The proposed statement

> every support-minimal order-q integer annihilator has Newton polytope of width bounded only by q

is false, even in one dimension.

A valid replacement is proved below. It eliminates the exact-zero branch by evaluating the same annihilator at several powers and using a Vandermonde argument.

The resulting theorem is rigorous and directly applicable to endpoint clusters. It does not by itself prove that every sufficiently large cluster contains a suitable annihilator.

---

## 1. Counterexample to the low-width theorem

Fix integers q >= 1 and L >= 1. On the one-dimensional lattice take support

```text
0, L, 2L, ..., qL
```

and coefficients

```text
c_k = (-1)^(q-k) binom(q,k).
```

For every polynomial f of degree less than q,

```text
sum_(k=0)^q c_k f(kL) = 0.
```

Equivalently,

```text
sum_(k=0)^q c_k (kL)^m = 0,       0 <= m < q.
```

Thus this is an order-q moment annihilator. Its lattice width is qL, which is arbitrarily large.

It is support-minimal: on any q distinct scalar points, the Vandermonde matrix for moments 0,...,q-1 is invertible, so no nonzero order-q annihilator can use only q support points.

Therefore support minimality alone cannot bound Newton width.

---

## 2. Why one arithmetic evaluation is insufficient

All realized divisor ratios lie in Q(i). As a vector space over Q, Q(i) has dimension two.

Hence any three elements of Q(i) satisfy a nontrivial rational, and therefore integer, linear relation. Exact relations

```text
sum_a c_a u_a = 0
```

are therefore not exceptional once the support has at least three points.

The earlier dichotomy

```text
small nonzero value  versus exact zero relation
```

cannot close the argument: the zero branch is pervasive for purely linear evaluations.

The repair is to use the whole power sequence

```text
S_m(P) = sum_a c_a u_a^m.
```

Distinct values u_a have an invertible Vandermonde matrix, so sufficiently many consecutive zero power sums force every coefficient to vanish.

---

## 3. Setup

Let

```text
u_a = product_j rho_j^(a_j),

rho_j = pi_j / conjugate(pi_j),
```

be distinct norm-one Gaussian rationals indexed by a finite support A.

Assume their arguments lie in one interval. Choose real lifts

```text
u_a = exp(i(alpha_0 + t_a)),
|t_a| <= delta.
```

Let

```text
P = sum_(a in A) c_a [a],
```

with nonzero integer coefficients, and let

```text
s = |supp(P)|.
```

Assume P annihilates phase moments through order q-1:

```text
sum_a c_a t_a^ell = 0,       0 <= ell < q.
```

A lifted tensor annihilator from the earlier note implies these scalar moment identities, so the theorem applies to that stronger notion.

Define the weighted support width

```text
W(P) = sum_j (U_j-L_j) log p_j,
```

where U_j and L_j are the maximal and minimal j-th exponents on the support.

---

## 4. Analytic estimate for every power

For m >= 1, put

```text
S_m(P) = sum_a c_a nu_a^m.
```

Since

```text
nu_a^m = exp(i m alpha_0) exp(i m t_a),
```

Taylor expansion and moment cancellation give

```text
|S_m(P)|
  <= ||c||_1 exp(m delta) (m delta)^q / q!.       (4.1)
```

Proof: subtract the Taylor polynomial of exp(ix) of degree q-1. Every Taylor term vanishes after summation, and the remainder is bounded by exp(|x|)|x|^q/q!.

---

## 5. Arithmetic estimate for every power

The exponent width of the monomials nu_a^m is m times the exponent width of the nu_a. Clearing Gaussian denominators exactly as in the original capacity note gives:

```text
S_m(P) = 0
```

or

```text
|S_m(P)| >= exp(-m W(P)/2).                       (5.1)
```

No pairwise decomposition is used.

---

## 6. Power-separating annihilator theorem

### Theorem

Let P be as above, with s distinct support values nu_a and phase-moment order q. Suppose that for every integer

```text
1 <= m <= s-1
```

one has

```text
||c||_1 exp(m delta) (m delta)^q / q!
  < exp(-m W(P)/2).                                (6.1)
```

Then no such nonzero P exists.

### Proof

For each m=1,...,s-1, equations (4.1), (5.1), and the strict inequality (6.1) force

```text
S_m(P)=0.
```

The zeroth moment condition gives

```text
S_0(P)=sum_a c_a=0.
```

Thus

```text
sum_a c_a nu_a^m = 0,       0 <= m <= s-1.
```

The s by s Vandermonde matrix

```text
(nu_a^m)_(0 <= m <= s-1, a in supp(P))
```

is invertible because the nu_a are distinct. Hence every c_a=0, contradiction.

---

## 7. Endpoint-scale corollary

At endpoint scale

```text
delta = C exp(-H/4),
H = log N,
```

condition (6.1) is implied, for all sufficiently large H, by the strict width inequality

```text
q H / 4 > (s-1) W(P) / 2
```

with enough margin to absorb

```text
log ||c||_1 + q log(s-1) - log(q!) + O_C(s delta).
```

Equivalently, asymptotically,

```text
W(P)/H < q / (2(s-1)).                              (7.1)
```

Therefore an endpoint cluster cannot support a nonzero order-q integer annihilator of support size s whose normalized Gaussian width lies strictly below q/(2(s-1)) and whose coefficient mass is subexponential in H.

This is a true low-width statement, but the bound necessarily depends on the support-size/order ratio. The false theorem omitted that dependence.

For a q-th finite difference on q+1 points, s=q+1, and the critical threshold becomes approximately

```text
W(P) < H/2.
```

So finite-difference annihilators occupying less than half of the conductor width are ruled out uniformly at the endpoint.

---

## 8. What remains

The theorem removes the pervasive exact-zero problem and gives a usable contradiction certificate.

The unresolved existence problem is now precise:

> Must every sufficiently large endpoint cluster contain an integer phase-moment annihilator P with
>
> ```text
> W(P)/H < q/(2(s-1))
> ```
>
> and subexponential coefficient mass?

This is not automatic. General support-minimal annihilators can have arbitrarily large width, as Section 1 shows, and high-dimensional sparse sets may have no suitable finite-difference pattern.

A successful next step must use endpoint concentration to force such a low-normalized-width annihilator, or else derive a dual structural certificate from its nonexistence.

The valid theorem should therefore be used as the arithmetic contradiction engine, while the combinatorial/geometric existence of a certificate is treated separately and honestly.