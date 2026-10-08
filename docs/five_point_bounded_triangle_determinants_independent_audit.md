# Independent audit of five-point bounded triangle determinants

This note audits `five_point_bounded_triangle_determinants.md`.  The
determinant comparison, conic reconstruction, metric scale, and annular
counterexample are correct.  The result is a bound in terms of the adjacent
determinant alphabet, namely `R=O(K^120)` with the displayed crude
constants; it is not a uniform point bound because no radius-independent
bound on `K` is proved.

## 1. Chord comparison and the alternating case

For adjacent chord lengths `a,b` with half-gaps `x,y`,

```text
determinant=ab sin(x+y),       a=2R sin x, b=2R sin y.
```

On a quarter-circle arc, `x+y<=pi/4`, so both cosines in
`sin(x+y)=sin x cos y+cos x sin y` lie between `1/sqrt(2)` and one.
This gives exactly the two constants in equation (2) of the source note.
Comparing successive determinants sharing their middle chord gives

```text
1/(sqrt(2)K)<=l_0/l_2,l_1/l_3<=sqrt(2)K.
```

If the largest and smallest chords have opposite index parity, both chords
of the smaller parity have length at most `sqrt(2)Km`.  Their midpoint
directions are distinct and differ by less than `pi`, so their integer
determinant is nonzero.  Using

```text
Delta<=4sqrt(2)L/R,
1/R<=2sqrt(2)K/(L^2m)
```

bounds its modulus by `32K^3m/L`.  Integrality therefore gives
`L/m<=32K^3`.  Substitution into

```text
|det(e_i,e_j)|<=4sqrt(2)L^3/R,
1/R<=sqrt(2)K/m^3
```

gives `|det(e_i,e_j)|<=2^18K^10`.  No direction-parallelism or endpoint
case is omitted.

## 2. Rank-four conic and the metric determinant

Let `A=[e_0 e_1]`, `d=det A>0`, and
`U=dA^(-1)(P-P_0)`.  The two coordinates of `dA^(-1)e_j` are signed
integer chord determinants, so all `U_i` are integral with the stated
coordinate bound.  The transported Euclidean metric is

```text
G=A^T A/d^2,       det G=(det A)^2/d^4=1/d^2.       (1)
```

The rank argument for the conic is valid over `Q`.  For each of five
circle nodes, the product of lines pairing the other four nodes vanishes
at those four and not at the chosen node, because no line contains three
circle points.  Thus evaluation of all quadratic monomials has rank five.
At `U_0=0`, its evaluation row has only the constant entry nonzero;
removing that row and the constant column leaves rank four for the four
equations in the five coefficients of a conic through zero.

The signed maximal minors therefore produce a nonzero primitive integer
kernel vector bounded by
`H=4!(16B^2)^4`.  Its quadratic matrix `G_0` is proportional to (1):

```text
G=lambda G_0,       lambda>0.
```

Since `4det G_0=4ac-b^2` is a positive integer,
`det G_0>=1/4`, and (1) gives `lambda<=2/d`.  Completing the square gives

```text
R^2=(lambda/4)(u,v)G_0^(-1)(u,v)^T.
```

Every inverse entry is at most `4H` in modulus.  Hence the quadratic form
is at most `16H^3`, proving `R^2<=8H^3/d<=8H^3`.  Since
`B=O(K^10)` and `H=O(B^8)`, this is `R=O(K^120)`.

## 3. Annular comparison

For the integer-parabola family, direct expansion gives the exact norm
defect

```text
|P_j|^2-R^2=j^3(2Q+j).
```

For `j<=M<=Q`, division by `|P_j|+R>=2R>=Q^3` gives the radial error
`3M^3/Q^2`.  Also `|P_M-P_0|<=3QM`.  The positive dot product with
`P_0`, together with the strictly signed derivative determinant, places
all arguments monotonically in one interval of width below `pi/2`.
Therefore

```text
Delta<=2 arcsin(|P_M-P_0|/(2R))<=2|P_M-P_0|/R,
Delta sqrt(R)<=6sqrt(2)M/sqrt(Q).
```

With `M=floor(Q^(1/4))`, both the normalized angular span and radial error
tend to zero, while the number of primitive consecutive chords grows and
their adjacent determinants remain exactly two.  The error is still much
larger than the radius quantization scale `1/R=Theta(Q^(-3))`, so it does
not contradict the exact-concyclicity theorem.

## 4. Audit result

No missing constant, rank hypothesis, metric determinant, or annular error
factor was found.  The existing checker again passes all 56 exact conic
reconstructions and all integer-parabola identities.  The proof establishes
a useful finite rigidity statement but leaves the growth rate unchanged.
