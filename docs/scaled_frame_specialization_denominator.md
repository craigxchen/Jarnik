# The degree-seven frame has an unavoidable specialization denominator

The rational frame in [the doubled contact construction](scaled_split_contact_profile_q2.md)
meets its four aggregate norm-contact conditions, but it cannot produce
integer frames at rational parameters of unbounded height. Its denominator
has an elementary lower bound at **every** rational parameter. No generic
specialization, recurrence, or approximation theorem is needed for this
bound.

Write

```text
f(t)=-2+(t^4-t^3)/(t^7+1),
U_f=((f,f-1),(1,1)).
```

For a primitive integer pair `(a,b)`, put `H=max(|a|,|b|)`. Exclude the
pole `a+b=0`. The value at `b=0` is interpreted projectively as `f(infinity)=-2`.
If `q(a,b)` is the positive denominator of `f(a/b)` in lowest terms, then

```text
q(a,b) >= H^6/4.                                      (1)
```

The exponent six is sharp. In particular the only integral value of this
frame is the constant matrix `U_(-2)`, attained at parameters `0,1,infinity`.
This is a limitation of this particular rational frame, not a proof about
all possible frames or the original lattice-arc bound.

## 1. Exact cancellation is at most two

Homogenize to common degree seven:

```text
Q(a,b)=a^7+b^7,
S(a,b)=a^3 b^3(a-b),
P(a,b)=S(a,b)-2Q(a,b),
f(a/b)=P(a,b)/Q(a,b).
```

Primitivity gives `gcd(Q,a)=gcd(Q,b)=1`. Hence

```text
gcd(P,Q)=gcd(S,Q)=gcd(a-b,a^7+b^7) divides 2.          (2)
```

For the last assertion, reduce modulo `a-b` to get `2b^7`, and use
`gcd(a-b,b)=1`. Thus

```text
q(a,b)>=|a^7+b^7|/2.                                 (3)
```

## 2. The nonreal pole factor costs six powers of height

The exact factorization is

```text
a^7+b^7=(a+b)B(a,b),
B=a^6-a^5b+a^4b^2-a^3b^3+a^2b^4-ab^5+b^6.
```

For all real `a,b`,

```text
B(a,b)>=H^6/2.                                       (4)
```

If `ab<=0`, substitution of `|a|,|b|` makes every term of `B` nonnegative,
and one term is `H^6`. If `ab>=0`, put `x=|a|,y=|b|`. Then
`B=(x^7+y^7)/(x+y)>=H^7/(2H)`. The origin is irrelevant here.

Since `a+b` is a nonzero integer, (3)--(4) prove (1). The six complex
poles encoded by `B` cannot be approached along the real rational
parameter line; their full degree remains in the denominator.

If `f(a/b)` is integral, then `q=1`, so `H<=4^(1/6)<2`. The primitive
projective points with `H=1` are `0,1,-1,infinity`. Excluding the pole
`-1` leaves `0,1,infinity`, and each gives `f=-2`.

## 3. Sharpness and the frame consequence

Take `a=1-N,b=N`. These integers are coprime and have opposite parity,
so (2) gives no cancellation. As `N` tends to infinity,

```text
H=N,
q=N^7-(N-1)^7=7N^6+O(N^5).                           (5)
```

Thus six cannot be replaced by a larger uniform exponent in (1).
The degree-seven map also has `h(f(a/b))=7log H+O(1)`, so its reduced
denominator accounts for at least six sevenths of the frame-value
height, up to a bounded error. The sequence (5) attains that ratio.

Every positive integer clearing the entries of `U_f(a/b)` is divisible
by `q(a,b)`: its only nonconstant entries are `f` and `f-1`, with the same
lowest denominator. Consequently there is no height-diverging sequence
of rational parameters for which such clearing integers have logarithm
`o(log H)`. A fixed rational source Möbius reparametrization does not
repair this, since it changes parameter height only by `O(1)` and (1)
already covers all rational parameters.

## 4. No fiber has support of size at most two

There is a further fixed-map constraint useful with the
[rational-parameter pole classification](rational_parameter_quasi_integral_poles.md).
The derivative numerator factors as

```text
(t^4-t^3)'(t^7+1)-(t^4-t^3)(t^7+1)'
 =t^2(-3t^8+4t^7+4t-3).                              (6)
```

The degree-eight factor is squarefree and coprime to the denominator.
Squarefreeness is certified by reduction modulo five in the attached
checker. The ramification indices at `0` and `infinity` are three;
both points lie over `-2`. All other ramified points have index two.
The fiber over `-2` is exactly `3[0]+[1]+3[infinity]`, so it has three
support points. Every other degree-seven fiber has support size at
least four, since its multiplicities are at most two.

In particular no rational target Möbius change can turn this map into
one with either one rational pole or just two real quadratic-conjugate
poles. Once the separately stated fixed-map pole classification is
applied, even such a fixed target change cannot create rational
specializations with subheight denominators. This last consequence uses
Roth's theorem through that classification; the explicit bound (1) and
the integral-value classification above are elementary.

## Verification

`python3 docs/check_scaled_frame_specialization_denominator.py` checks
the exact gcd and bound on primitive pairs, the integral-value list,
the sharp sequence, and the derivative's squarefreeness certificate.
The general inequalities are proved above, not inferred from the finite
checks.
