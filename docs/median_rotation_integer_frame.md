# Median rotation and a small integer frame

A common Gaussian rotation can remove all pairwise Gaussian gcds among
three chosen primitive rows. The resulting three directions admit a
small integral frame, while one large unimodular matrix contains their
angular compression. This is an exact normal form, not a radius descent
or a proof of the uniform endpoint bound.

## 1. The rotation

Let `P_i=X_i+iY_i` be odd conjugate-coprime Gaussian integers with
distinct projective directions. Write

```text
D_ij=X_i Y_j-Y_i X_j,
t_ij=|D_ij|/N(gcd_G(P_i,P_j)),       1<=t_ij<=T.
```

Choose three rows, labelled 1,2,3, and set

```text
g=gcd_G(D_23 P_1,D_13 P_2),
A=D_23 P_1/g,       B=-D_13 P_2/g.
```

Cramer's identity gives `A+B=-D_12 P_3/g`. Thus `A,B,A+B` are
pairwise Gaussian-coprime. Rotate every row by `1/g` and choose its
primitive integral representative by the explicit formula

```text
c_i=gcd(|Re(P_i bar(g))|,|Im(P_i bar(g))|),
Q_i=P_i bar(g)/c_i.                                      (1)
```

The first three representatives are rational real multiples of
`A,B,A+B`. Primitive normalization only divides integral coordinate
content, so `Q_1,Q_2,Q_3` are pairwise Gaussian-coprime.

There is no individual quarter-turn hidden in (1). Since each `P_i`
is odd, the valuation of `g` at `1+i` is the minimum of two even
integers. The coordinate content removes its entire even ramified
valuation. Consequently every `Q_i` is again odd and conjugate-coprime.
The ratios `(Q_i/bar(Q_i))/(Q_j/bar(Q_j))` equal their original
counterparts. Their primitive Gaussian numerators agree up to sign,
and every `t_ij` is preserved exactly.

Here is an intrinsic description of the same rotation. For an odd
split prime with chosen orientation `pi`, put

```text
n_i=v_pi(P_i)-v_bar(pi)(P_i),
s=v_pi(g)-v_bar(pi)(g).
```

The signed valuations after (1) are `n_i-s`. Pairwise coprimality of
the first three representatives means at most one is positive and
at most one negative. Therefore

```text
s=median(n_1,n_2,n_3).                                  (2)
```

This also proves that the displayed Gaussian gcd realizes the
prime-by-prime median without any factorization algorithm.

## 2. The integer frame and its size

Choose `W_0 in SL_2(Z)` taking the primitive vector `Q_1` to `(1,0)`.
Write `W_0 Q_2=(a_0,b)` and `W_0 Q_3=(c_0,d)`. Pairwise Gaussian
coprimality gives

```text
1<=|b|=t_12<=T,   1<=|d|=t_13<=T,
1<=|a_0 d-b c_0|=t_23<=T.
```

An integral shear fixes `(1,0)` and changes `a_0` by a multiple of
`b`. Choose it so that `|a|<=|b|/2`; let `W` be the resulting matrix.
Then `c=(ad-D_23(Q))/b` gives `|c|<=3T/2`. In particular, with

```text
U=W^(-1),       v_i=W Q_i,
```

we have

```text
U in SL_2(Z),   Q_i=U v_i,
v_1=(1,0),     max_(a=1,2,3) ||v_a||_infinity<=2T.       (3)
```

If `H=max_(a=1,2,3)||Q_a||_infinity`, then

```text
H/(4T)<=||U||_infinity<=2TH.                            (4)
```

The lower bound follows from `Q_a=Uv_a`. For the upper bound,
the first column of `U` is `Q_1` and the second is `(Q_2-aQ_1)/b`.
Thus its norm is at most `H/|b|+H/2<=3H/2`, which is even stronger
than the claimed bound for `T>=1`.

## 3. Heights forced by a full independent cut profile

Suppose, as in the endpoint reduction, that there are `m>=3` rows
and independent odd split Gaussian blocks with

```text
P_i=K_i product_(S containing i) H_S,
log|K_i|<=sigma,
(1-eta)w<=log N(H_S)<=(1+eta)w
```

for every nonempty subset `S` of the `m` labels. Blocks and their
conjugates have disjoint prime support. The empty subset is irrelevant.
Set

```text
A_m=2^(m-3),       E=3 A_m eta w+4sigma.
```

For the three frame rows and all other rows respectively,

```text
|log|Q_a|-A_m w|<=E,       a in {1,2,3},
|log|Q_j|-2 A_m w|<=E,     j outside {1,2,3}.             (5)
```

To prove this, initially omit the corrections `K_i`. Formula (2)
subtracts the majority bit among the first three labels. A fixed
one of those labels differs from their majority on `2^(m-2)` cuts;
an outside label differs on `2^(m-1)` cuts. Each selected block
contributes one half of its logarithmic norm to the row modulus.
The median is Lipschitz in the sum of its three input errors.
Restoring `K_i` changes each row's logarithmic modulus by at most
`log|K_i|+sum_(a=1,2,3)log|K_a|<=4sigma`. This proves (5), with
room in the stated error.

For every pair, its primitive Gaussian numerator `h_ij` satisfies

```text
|log|h_ij|-2 A_m w|<=2 A_m eta w+2sigma<=E.              (6)
```

Indeed the separating cuts number `2^(m-1)`. The signed-valuation
formula for `h_ij` takes the absolute difference of the two row
valuations, and the corrections cost at most `2sigma`.

Since `|h_ij|=|Q_i Q_j|/N(gcd_G(Q_i,Q_j))`, (5)--(6) give

```text
|log N(gcd_G(Q_a,Q_j))-A_m w|<=3E,
|log N(gcd_G(Q_j,Q_l))-2 A_m w|<=3E                    (7)
```

for frame/outside and two distinct outside rows. Using the exact
determinant/residue identity and (3), this implies

```text
A_m w-3E <= log||v_j||_infinity
          <= A_m w+3E+log(2T^2),
|log||U||_infinity-A_m w|<=E+log(4 sqrt(2) T),
2 A_m w-3E <= log|det(v_j,v_l)|
            <=2 A_m w+3E+log T.                        (8)
```

For the first line, `v_j=(x_j,y_j)` has
`y_j=det(Q_1,Q_j)` and `x_j=(a y_j-det(Q_2,Q_j))/b`.
The lower bound follows from the nonzero integer residue `t_1j`;
the upper bound follows from `|a|<=T/2`, `|b|>=1`, and (7).
The second line uses (4)--(5) and comparison of Euclidean and
coordinate norms. The last line is just (7) and invariance of
determinants under `SL_2(Z)`.

## 4. What this does and does not give

With `eta->0` and `log T,sigma=o(w)`, the first three frame vectors
have subexponential height. Both the matrix `U` and every outside
frame vector have logarithmic height `A_m w+o(w)`. Two outside
vectors have determinant of logarithmic size `2 A_m w+o(w)`.
This records an actual integer realization of the geometry, rather
than only prescribing formal prime contents.

Applying `W` changes the Euclidean circle metric. Its height is
large, and the least circle radius of the transformed directions
has not been bounded. Thus (8) is an arithmetic normal form and
does not produce a smaller endpoint configuration. It complements
[the residue-content identities](residue_content_descent.md), whose
global Smith normalization does not remove all three pairwise gcds.

The exact checks in
[check_median_rotation_integer_frame.py](check_median_rotation_integer_frame.py)
verify the rotation, residues, frame inequalities, signed medians,
and cut counts. They supplement the proofs, and do not establish
the remaining uniform bound.
