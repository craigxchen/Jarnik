# K5 fixed-ratio algebra and bounded-search journal

The K5 fixed-ratio subsection contains one exact characteristic-zero theorem;
the modular, ten-factor, and gap-free sections remain bounded diagnostics. None
of these results is a claim about the circle theorem.

## K5 fixed-ratio modular Groebner diagnostic

Source: [check_k5_fixed_ratio_modular_diagnostic.py](check_k5_fixed_ratio_modular_diagnostic.py).
It fixes adjacent K5 edge weights `x_01=1,x_02=2`, substitutes the four exact
linear vertex-sum constraints, forms the four cubic vertex-difference polynomials,
and runs Buchberger over

```text
p = 1000003
```

with the displayed monomial order (`max` on exponent tuples). The completed run had
8,128 pair reductions and a 128-element basis. Reducing the partial product

```text
P9 = product_(j=1..9) (x_0^2-x_j^2)
```

by that modular basis gave zero at the factor `(0,9)`, which is the ninth factor
in this product. Thus the completed run verifies the finite-field diagnostic
`P9 = 0 mod I` for this fixed adjacent ratio. Since `P9` is a factor of the
full pairwise product, it implies that full product also vanishes in this
finite-field quotient. Neither conclusion transfers to the rationals or reals;
the disjoint fixed-edge orbit was not completed. No global ratio claim is made.

Interreducing the 128-element modular basis produced 16 polynomials, with 54
standard monomials. These are reproducibility data for the finite-field run,
not a characteristic-zero certificate.

An exact characteristic-zero check is available for this same fixed ratio in
[check_k5_fixed_ratio_qq_reduction.py](check_k5_fixed_ratio_qq_reduction.py).
The displayed linear substitution is checked exactly; its pivot matrix,
generated directly from graph incidence, has determinant one. Thus it
parametrizes every solution of the linear equations at the fixed ratio.
The checker uses SymPy's completed 28-element Groebner basis over `QQ`, variables
`(a,b,c,d)`, grevlex order, and method `f5b`. Exact reduction gives:

```text
linear_differences_zero=True pivot_det=1
zero_dimensional=True basis_len=28 domain=QQ order=grevlex
P9 remainder=0 zero=True
```

This checker requires SymPy. The independently replayed run used version
`1.14.0`, available locally with:

```text
/Users/cxc/.cache/uv/environments-v2/high-level-tutorial-14eec61d3278d4de/bin/python docs/check_k5_fixed_ratio_qq_reduction.py
```

Therefore every characteristic-zero common zero of these four fixed-ratio
cubic equations has `P9=0`, and hence has a repeated absolute edge value.
Because `x_01=1` and the `j=1` factor is `1-2^2=-3`, one of the
eight unfixed edge values has absolute value one. This proves the fixed-ratio
algebraic exclusion of ten distinct absolute values, but says nothing about
arbitrary edge ratios. The optional square-sum reductions have nonzero normal
forms, so those polynomials are not in this ideal; without a radicality claim,
this does not assert that they fail at every common zero. For example, the
first normal form is

```text
2*(a*b-a*c-b**2+b+2*c+d**2-3*d)
```

The completion prefix was:

```text
step 1000 basis 126 queue 6875
step 2000 basis 128 queue 6128
...
step 8000 basis 128 queue 128
done steps 8128 basis 128
factor 0 1 remterms 1
factor 0 2 remterms 54
factor 0 3 remterms 50
factor 0 4 remterms 50
factor 0 5 remterms 52
factor 0 6 remterms 50
factor 0 7 remterms 50
factor 0 8 remterms 50
factor 0 9 remterms 0
P9 rem terms 0
```

## Exact K5 triple-product identity

There is a useful algebraic reformulation, but it does not by itself give a
new nonexistence argument. For vertices (i,j), let (r,s,t) be the three
remaining vertices. Equality of the vertex sums and of the vertex cube sums
implies

```text
(x_ir+x_is)(x_ir+x_it)(x_is+x_it)
  = (x_jr+x_js)(x_jr+x_jt)(x_js+x_jt).
```

Indeed, for any triple (a,b,c), with (q=ab+ac+bc), (r=abc), and
(u=a+b+c), one has

```text
a^3+b^3+c^3-u^3 = -3u q + 3r,
(a+b)(a+c)(b+c) = u q-r.
```

Thus equal (u) and equal cube sum force equal pair-sum products. Conversely,
under equal (u), the two identities show that this product equality is just
the cube-sum equality in another form. Expanding the ten such identities after
the adjacent fixed-ratio substitution produces cubic polynomials without a
new common factor; consequently this identity alone does not classify the
real solutions.

## Ten-factor integer-box search

Source: [check_tenfactor_integer_box_search.cpp](check_tenfactor_integer_box_search.cpp).
Compile and run:

```text
g++ -O3 -std=c++17 docs/check_tenfactor_integer_box_search.cpp -o /tmp/search_tenfactor
/tmp/search_tenfactor
```

It exhausts all `binom(M,10)` ten-tuples for `M=10,12,14,16,18`, groups all
1,024 subset masks by the exact pair (first sum, third sum), and checks fifth-sum
order only for fibers of size at least five. Maximum fiber sizes were:

```text
M=10: 2
M=12: 3
M=14: 3
M=16: 3
M=18: 4
```

The full `M=18` run has 43,758 coefficient tuples and 44,808,192 masks. No
five-subset fiber occurred, so there were no source outer-cut candidates and no
associated normalized `C^2` value. The largest fiber was

```text
a=(1,2,5,7,8,9,10,14,16,18)
common first/cube sums=(45,7695)
rows=[408,421,602,615]
```

with fifth sums (1,635,975,1,648,575,2,039,175,2,051,775). This is bounded
integer evidence only; it makes no claim outside coefficient bound 18.
