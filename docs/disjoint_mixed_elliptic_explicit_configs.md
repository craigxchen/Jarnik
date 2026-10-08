# Two exact rational six-point configurations

This note records two small exact orbit computations on the fixed
disjoint-support curve.  It is intended as data for later frame
experiments; it makes no Gaussian-prime or pairwise-coprimality claim.

Use the quotient quartic

```text
F: Y^2=T(67144321T-212248804)(T-4)(46513T-208624).
```

Put `K=-4*212248804*208624` and use the exact cubic conversion

```text
x=K/T,                 y=K*Y/T^2,
y^2=(x+4*67144321*208624)(x+212248804*208624)
       (x+4*212248804*46513).
```

Take

```text
R0=(-1/2,45299439/4)
```

on the quartic, with the branch point `(0,0)` as origin.  If `P=nR0`,
the quotient map from the source is `q=2P-R0=(2n-1)R0`.  The exact
group law gives the following two target values, for `n=2,3`:

```text
T_3 =
5637458512269621578866836997714629776400 /
55746306906090343697611717689527150819,

T_5 =
69602127761781652956842789689474871015540248434267661885313166304077864877141852015392606565051804401779065625650575961 /
11684797357306360493184597792210432544867502183611853033481491744644578189986105242549395021845832501101415061343663068.
```

For a target `T`, the two original-coordinate roots `x_old` for parameter
`lambda` are

```text
x_old_lambda,+/- =
 (73lambda + T - 11lambda*T +/- sqrt(Delta_lambda(T)))
 / (2*(1-lambda*T)),
```

where the three exact discriminants are those in Section 4 of the
[positive-rank certificate](disjoint_mixed_elliptic_six_curve.md#4-exact-quotient-quartic-and-positive-rank).
Normalize each root by

```text
X=(x+1)/(16-4x).
```

Choosing the `+` root for each of `lambda_a=0`,
`lambda_b=-240/5329`, and `lambda_c=8/257`, and using primitive rows
`(p,q)` for `X=p/q`, gives the following explicit six-point rows.  The
first three rows are the anchors `infinity,0,1`.

For `T=T_3`:

```text
(1,0), (0,1), (1,1),
(-75832963907715415721, 288332793673229656864),
(-47258107518159388601, 155390021886869913184),
(-251457826378440801557, 1154031153698118652288).
```

For `T=T_5`:

```text
(1,0), (0,1), (1,1),
(-160065375142999604105092096433964946346082445523858934249369,
 77175515359484945240842427490971365862846934153385962385406),
(-733070910594489668954419384627817088070699473358841497530727581,
 355275812690162480828075007809824042341448071994293848547651268),
(-116015272229964288194998137746868294085800989161644933738609,
 55060999269003577737560366529906994830474588893147156413492).
```

The omitted minus-root choices give the other rational lifts over the
same target.  The exact checker
[check_disjoint_mixed_elliptic_explicit_configs.py](check_disjoint_mixed_elliptic_explicit_configs.py)
recomputes the cubic additions, verifies `q=(2n-1)R0`, checks the
quartic equation, checks all three discriminants are rational squares,
and verifies every displayed row is primitive and distinct.  It does
not factor the very large determinants or assert any Gaussian lifting
property.
