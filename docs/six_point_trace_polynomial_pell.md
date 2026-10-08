# The rational trace gives an exact polynomial Pell identity

This note retains the constants and squareclasses in the rational Belyi
trace from
[`six_point_belyi_period_map.md`](six_point_belyi_period_map.md).
The resulting identity is exact, but its specialization imposes no new
condition on a rational point already lying on the circle cover.

Fix a primitive integral nonresonant central vector `u`, a controlled
rational hyperbolic frame, and one ruling. Write `U=max |u_i|`. Use the
integer raw forms of that frame, so `Delta` has degree twelve and `K` has
degree twenty. Both are squarefree and coprime. Choose the ten integer
minor quadratics `M_I`, one for each unordered balanced partition, with
their actual contents and signs. Put

```text
m_I=|sum_(i in I) u_i|>0,
B=product_I M_I^(m_I),
D=deg B=2 sum_I m_I=deg F.
```

No monic or primitive normalization of `B` is understood. Its divisor is
exactly the pole divisor of the base trace `G=(F+F^(-1))/2`. Consequently
there is a rational homogeneous form `P` of degree `D` such that

```text
G=P/B,                    gcd(P,B)=1 in Q[s,t].     (1)
```

Here `P` need not be integral: fixing the actual integer product `B` can
force rational coefficients in its numerator.

## 1. Factorization with the actual squareclasses

The trace passport gives homogeneous forms `A,C` over `Q` and nonzero
rational constants `a,b` such that

```text
P+B=a A^2,
P-B=b Delta^5 C^2,                                 (2)
deg A=D/2,          deg C=(D-60)/2.
```

One can initially take `A,C` primitive integer forms. Unique
factorization over `Q` proves (2): every geometric root of `P+B` has
multiplicity two, while the roots of `P-B` have multiplicity five on
`Delta` and multiplicity two elsewhere. Homogeneous factorization includes
the parameter infinity. The passport also shows that `A,C,Delta,B` have
pairwise disjoint zero sets, and `A,C` are squarefree unless constant.

Let `H=(F-F^(-1))/(2 i d)` as in the trace note. Since
`G^2+Delta H^2=1`, multiplication of (2) gives

```text
H^2=-ab Delta^4 A^2 C^2/B^2.
```

It follows that `-ab=c^2` for some nonzero rational `c`, and the sign of
`c` may be chosen so that

```text
H=c Delta^2 A C/B,
a A^2-b Delta^5 C^2=2B.                            (3)
```

Indeed a rational function whose square is a constant is itself constant;
thus this squareclass conclusion holds over `Q`, not just over an
algebraic closure. Equivalently, replace `C` by `(c/a)C`, retaining this
rational rescaling, to obtain the particularly simple normalization

```text
P+B=a A^2,       P-B=-a Delta^5 C^2,
2B=a(A^2+Delta^5 C^2),
H=a Delta^2 A C/B.                                 (4)
```

The final `C` in (4) need not be primitive or integral. In particular,
discarding `a` or pretending both square factors retain their original
integer normalization would lose arithmetic information.

There is also an exact derivative identity. In a chart `r` with
`d log F=i kappa d^3 dr/K`, differentiating the trace gives

```text
G'=-kappa Delta^2 H/K,
P'B-PB'=-kappa c Delta^4 A C (B/K).                 (5)
```

The quotient `B/K` is a polynomial up to the fixed nonzero rational factor
relating `K` to `product M_I`. Both sides have homogeneous degree `2D-2`
when interpreted as the numerator of the corresponding projective
differential. Formula (5) expresses precisely the already known critical
divisor; it is not an additional arithmetic hypothesis.

## 2. A safe coefficient-height cost

The controlled frame has coefficient height `T<=C_0 U^C`. Hence the
coefficients of its fixed-degree forms `Delta,K,M_I,Z_i` are at most
`U^O(1)`, after the specified common clearing. Also `D=O(U)` and
`W=sum_(u_i>0) u_i=O(U)`.

For clarity one can bound the trace without any root approximation.
Using `Z_i Z_i_tilde=4AK` (the `A` in this sentence is the conic
coefficient, not the square factor in (2)), write

```text
F = product_(u_i>0) Z_i^(u_i)
    product_(u_i<0) Z_i_tilde^(-u_i) / (4AK)^W.
```

Reducing powers with `d^2=Delta`, and taking the real part, gives an
integer polynomial numerator and denominator for `G` of degree `O(U)`
and coefficient magnitude at most `exp(O(U log(2U)))`. This also directly
bounds the actual product `B` by the same quantity.

Cancellation and extraction of the square factors preserve an explicit,
though weaker, safe bound

```text
log arithmetic_height(P,A,C,a,b,c)
    <= C_1 U^2 log(2U).                            (6)
```

Here arithmetic height includes the absolute numerators and common
denominators of rational coefficients; (6) also holds after the rescaling
in (4), with an enlarged absolute constant. A factor estimate with no
accumulating degree loss suffices. For an integer polynomial `J`, define
its Mahler product as `|lead(J)| product_alpha max(1,|alpha|)`, including
root multiplicity. This product is multiplicative and is at least one
for every nonzero integer polynomial. Expansion into elementary
symmetric functions bounds each coefficient of a degree-`n` factor by
`2^n` times its Mahler product. Jensen's circle-average identity and
Parseval give `Mahler(J)<=sqrt(sum coefficients(J)^2)`. Gauss's lemma
therefore bounds every primitive integer factor of `J` by
`2^n sqrt(n+1) max|coefficients(J)|`. This estimate applies successively
without multiplying the logarithmic height by the degree. Scalar
contents and rational proportionality constants are ratios of nonzero
coefficients subject to the same bounds. These estimates even give a
sharper bound than (6), but (6) already makes
the essential cost explicit: this argument does not produce coefficient
height polynomial in `U` for the high-degree Pell identity.

## 3. What specialization does and does not add

At a rational cover point `d^2=Delta(s,t)` with `B(s,t)!=0`, (4) becomes

```text
2 B(s,t)/a = A(s,t)^2 + (d^5 C(s,t))^2.             (7)
```

Thus it is a rational Gaussian norm. For example, any prime `p=3 mod 4`
has even valuation in the nonzero rational number `2B(s,t)/a`. All
denominators and the squareclass of `a` must remain in that statement.

But (4) holds identically for every base parameter. Once `d^2=Delta` is
imposed, (7) follows by substitution, without any further divisibility,
integrality, or size requirement on `(s,t,d)`. Conversely, (7) as a norm
condition on `B` alone generally forgets the requirement that its second
coordinate be exactly `d^5 C(s,t)`. It therefore cannot be promoted to an
equivalent new height obstruction by treating the two norm coordinates
as independent integers.

In particular, the fixed twelve index-five points supply the factor
`Delta^5` and the prescribed derivative (5), while the degree `D`, pole
orders, constants, and coefficient heights still move with `u`. The
identities give a structured reformulation of the existing cover
equation. No estimate comparing the actual parameter height or angular
diameter with `U`, and no improvement in radius growth, follows here.
