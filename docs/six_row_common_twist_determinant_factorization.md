# The six-row common-twist determinant factors into scalar characters

For six actual common-unit Gaussian points, two balanced half-characters
have a common rational norm twist. Comparing their real coordinates seems
to give a simultaneous Pell approximation. Its **exact** factorization,
however, is the product of the imaginary coordinates of two ordinary
integer row characters. Thus this two-character determinant supplies no
height inequality beyond the scalar character bounds. This does not
exclude relations involving three or more characters or arithmetic
cancellation at the evaluated values.

The companion [odd trigonometric determinant theorem](common_twist_triple_determinant_factorization.md)
extends the exact factorization to a specific determinant family with
three, five, or any odd number of balanced characters. Other additive
identities and extra arithmetic information remain outside these results.

## Exact identities, including nested prime allocations

Remove the complete common Gaussian factor. At each varying split prime
`p`, choose `pi_p` above `p` and let `a_i(p)` be its allocation in row
`i`. Choose two distinct balanced vectors `lambda,mu in {±1}^6` modulo
global negation, and put

```text
c_lambda(p)=sum_i lambda_i a_i(p),
c_mu(p)=sum_i mu_i a_i(p),
L_+=(lambda+mu)/2,       L_-=(lambda-mu)/2.
```

Both `L_±` are nonzero integer vectors with coordinate sum zero. One has
`||L_+||_1, ||L_-||_1 = 2,4` in some order. Define the conjugate-primitive
Gaussian product `A_v` by taking `pi_p^c` when `c=v·a(p)>0`, and
`bar(pi_p)^(-c)` when `c<0`. Let

```text
A_lambda=x_lambda+i y_lambda,   A_mu=x_mu+i y_mu,
B_+=A_(L_+)=a_++i b_+,          B_-=A_(L_-)=a_-+i b_-.
```

All balanced signed exponents have the same parity at each prime. With
`ε_p=sum_i a_i(p) (mod 2)` and `t=product_(ε_p=1) p`, there are
positive integers `h_lambda,h_mu` satisfying

```text
Norm(A_lambda)=t h_lambda^2,
Norm(A_mu)=t h_mu^2.
```

Write `g=gcd(h_lambda,h_mu)`. The literal Gaussian identities are

```text
B_+ B_-          = (h_mu/g) A_lambda,
B_+ bar(B_-)     = (h_lambda/g) A_mu.                (1)
```

To check content prime by prime, use
`|q_+|+|q_-|=max(|c_lambda|,|c_mu|)` for
`q_±=(c_lambda±c_mu)/2`. The rational content in the first product is
`p^[(max(|c_lambda|,|c_mu|)-|c_lambda|)/2]`, exactly the local factor
of `h_mu/g`; the second identity is symmetric. This argument uses the
total allocation `a_i(p)` and permits arbitrary nested threshold layers
at the same prime.

Taking real parts of (1) gives the exact factorization

```text
h_mu x_lambda - h_lambda x_mu = -2g b_+ b_-,
h_mu x_lambda + h_lambda x_mu =  2g a_+ a_-.       (2)
```

Also, taking absolute values in (1),

```text
|B_+| |B_-| = sqrt(t) lcm(h_lambda,h_mu),
log(sqrt(t) lcm(h_lambda,h_mu))
  = [H(L_+)+H(L_-)]/2
  = (1/2) sum_p max(|c_lambda(p)|,|c_mu(p)|) log p,
H(L)=log Norm(A_L).                                (3)
```

The exact product in (2), rather than only divisibility by `g`, is the
decisive content calculation. Any extra size of this determinant comes
from the two scalar coordinates `b_+` and `b_-` themselves.
Here `x_lambda,x_mu` are the **signed** real parts of the canonical
Gaussian products. If one instead uses `|x_lambda|,|x_mu|` to compare
positive cosines on a short arc, its difference is the absolute value
of either the difference or the sum in (2), according to the real-part
signs. The latter factors through `a_+a_-`; it must not silently be
identified with the small determinant `-2g b_+b_-` without controlling
that sign branch. The scalar derivation below needs no such branch choice.

## Radius inequality and its scalar origin

Suppose the six points are distinct, have common modulus `R>1` and one
common Gaussian unit, and lie on an arc of length at most `C sqrt(R)`
with `C<=2`. Let `Delta<=C/sqrt(R)` be their angular width and
`W=log R^2`. Affine allocation rigidity makes both `B_±` nonunits.
They are conjugate-primitive with odd norm, so their imaginary
coordinates `b_±` are nonzero integers. In particular the real-part
determinant in (2) is nonzero, without a separate Pell argument.

For an integer row vector `L` with coordinate sum zero, the source
angles give

```text
dist(arg A_L, pi Z) <= (||L||_1/4) Delta.
```

The two `l1` norms here are two and four. Hence
`|B_short|>=2/Delta` and `|B_long|>=1/Delta`, from
`|b_±|>=1`. Equations (3) give

```text
sqrt(t) lcm(h_lambda,h_mu) >= 2/Delta^2 >= 2R/C^2,
log(sqrt(t) lcm(h_lambda,h_mu))
    >= W/2 - 2 log C + log 2.                      (4)
```

Equation (4) is exactly the average of the two ordinary scalar
row-character height bounds for `L_+` and `L_-`. A cruder estimate
based only on the nonzero integer determinant and the common conic
`x^2+y^2=t h^2` gives the same leading exponent with a weaker constant.
The common twist and the full rational content cannot turn this
two-character argument into a new radius exponent.

For comparison, on the formal equal-weight 31-cut six-row profile,
each `H(L_short)=16w` and `H(L_long)=24w`, so the left side of (4) is
`20w`; `W/2=31w/2`, leaving a leading margin `9w/2`. This is a formal
profile test, not an endpoint-arc realization. It shows why (4) cannot
recover the conditional aggregate-parity exclusion of that profile.

[The exact checker](check_six_row_common_twist_determinant_factorization.py)
verifies (1)–(3) over literal Gaussian products with arbitrary nested
allocations, and all 45 character pairs on the 31-cut profile.
