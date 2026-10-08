# Same-squareclass divisor solutions

Let `H,B,A` be positive integers, and let the positive divisor `n` of `H` satisfy

\[
 Y_n^2=2An-n^2-B^2,
\]

with `Y_n>=0` integral. Fix a squareclass `d`, so that two divisors
in it have the form `n=da^2`, `n'=db^2`. Multiplying the two equations
by `b^2` and `a^2`, respectively, and subtracting gives the exact identity

\[
 (bY_n-aY_{n'})(bY_n+aY_{n'})
 =(a^2-b^2)(B^2-d^2a^2b^2)
 =(a^2-b^2)(B^2-nn').
\tag{1}
\]

Assume `n!=n'`, `nn'!=B^2`, and `A>=H^2+B^2`. Then
\[
Y_n^2=2An-n^2-B^2\ge An,
\]
and similarly for `n'`, because `n<=H`. Hence
\[
bY_n+aY_{n'}\ge 2ab\sqrt{Ad}\ge2\sqrt A.
\]
The first factor in (1) is a nonzero integer, so its absolute value is at
least one.  On the other hand,
\[
|(a^2-b^2)(B^2-nn')|
 \le (n+n')(B^2+nn')
 \le 2H(B^2+H^2).
\]
Therefore every such pair obeys the explicit polynomial bound

\[
 A\le H^2(B^2+H^2)^2.                           \tag{2}
\]

If `A<H^2+B^2`, it is already bounded by a polynomial in `H+B`.
Thus, for `A>(H+B)^6`, each squareclass contains at most two solutions;
if it contains two, their product must be `B^2`. The exponent six is a
safe uniform threshold, since `H^2(B^2+H^2)^2<=(H+B)^6`.

Zero values `Y_n=0` cause no exception: they imply
`2An=n^2+B^2`, hence `A<=(H^2+B^2)/2`, below the large-`A` regime.
The assertion is conditional divisor structure and does not give a uniform
endpoint bound by itself.

The [checker](check_luna_divisor_squareclass_elimination.py) verifies exact
pair identities and the large-height consequence. The
[even-product argument](reciprocal_squareclass_runge_growth.md) extends this
two-label observation to squareclass relations of arbitrary even length.
