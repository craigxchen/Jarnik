# An integer eight-point obstruction for every order-eight Hadamard cut system

The intrinsic eight-point target is

\[
D+W_4\le 2W+B.
\tag{1}
\]

The following exact integer theorem excludes a sparse critical profile
which the ordinary determinant permits at leading weight. It uses actual
Gaussian integrality and the containing arc, with no polynomial-family
assumption and no restriction on prime sizes or exponents. It does not
prove (1) for arbitrary profiles or for this profile with independent
row-dependent twists.

**Theorem.** Consider the eight-cut model defined below, with independent
conjugate-coprime Gaussian integer blocks and one common Gaussian unit
class. Its eight points cannot be distinct and lie on an arc of length
at most \(\tfrac12\sqrt R\). A common Gaussian integer multiplier of all
the points is allowed. The same conclusion holds with **any** order-eight
Hadamard matrix whose anchor row is all ones; Section 7 proves this
generalization without classifying Hadamard matrices or enumerating
switchings.

## 1. The model and its critical cut weights

Index rows \(x\) and columns \(y\) by \(\mathbb F_2^3\). Write

\[
H_{xy}=(-1)^{x_1x_2+x\cdot y},\qquad
A_{xy}=\frac{1-H_{xy}}2.
\tag{2}
\]

The rows of \(H\) are orthogonal and have squared length eight; the
row indexed by zero is all ones, so \(A_{0y}=0\). The column sums are

\[
S_y=\sum_xH_{xy}=
\begin{cases}
4(-1)^{y_1y_2},&y_3=0,\\
0,&y_3=1.
\end{cases}
\tag{3}
\]

Thus the four columns in \(U=\{y:y_3=0\}\) are cuts of sizes two or
six, while the four columns in \(B=\{y:y_3=1\}\) are balanced cuts of
size four. There are three size-two columns and one size-six column.

For each column let \(Q_y\in\mathbb Z[i]\setminus\{0\}\). Require

\[
\gcd(Q_y,\overline{Q_y})=1,
\qquad
\gcd(Q_y,Q_z)=\gcd(Q_y,\overline{Q_z})=1\quad(y\ne z),
\tag{4}
\]

up to Gaussian units. Unit blocks are allowed. With one common nonzero
Gaussian integer \(d\), form

\[
z_x=d\prod_yQ_y^{A_{xy}}\overline{Q_y}^{\,1-A_{xy}},
\qquad R=|d|\prod_y|Q_y|.
\tag{5}
\]

The common unit may be absorbed into \(d\). Row-dependent Gaussian
units are not being freely inserted: the common-unit convention is
needed for the phase identities below. The model permits arbitrary
Gaussian blocks, not just linear polynomials or prime factors of
bounded height.

If all eight block weights \(w_y=2\log|Q_y|\) were equal to \(w\),
the primitive cut statistics would be

\[
W=8w,\quad V_2=4w,\quad V_4=4w,\quad D=16w=2W,
\quad D+W_4-2W=4w.
\tag{6}
\]

This is a leading-weight obstruction to (1) in the cut relaxation.
The theorem excludes every exact integer realization in (5) at arc
constant \(1/2\), even with unequal weights or unit blocks. It therefore
supplies an arithmetic exclusion for an actual critical mixed profile,
as distinct from exact balance alone.

## 2. Two exact small-power separation bounds

Let \(\gamma=a+bi\) be a conjugate-coprime nonunit Gaussian integer.
It lies on neither coordinate axis. Therefore \(a,b\ne0\), and

\[
|\sin(2\arg\gamma)|
=\frac{2|ab|}{|\gamma|^2}
\ge\frac{\sqrt2}{|\gamma|}.
\tag{7}
\]

Indeed one coordinate has absolute value at least
\(|\gamma|/\sqrt2\), while the other has absolute value at least one.

It also lies on neither diagonal: if \(|a|=|b|\), the Gaussian prime
\(1+i\) divides both \(\gamma\) and its conjugate. Hence

\[
|\sin(4\arg\gamma)|
=\frac{4|ab(a-b)(a+b)|}{|\gamma|^4}
\ge\frac1{|\gamma|}.
\tag{8}
\]

For an elementary proof, put \(u=\max(|a|,|b|)\),
\(v=\min(|a|,|b|)\), and \(e=u-v\). Then \(v,e\ge1\), so
\(2ve\ge v+e=u\). Consequently

\[
4uve\ge2u^2\ge u^2+v^2,
\qquad u+v\ge\sqrt{u^2+v^2},
\]

which proves (8). In particular \(\gamma^2\) and \(\gamma^4\) are
nonreal whenever these bounds are used. Primitivity is essential in the
fourth-power statement; it excludes the rational diagonal directions.

## 3. Balanced columns become exact fourth powers

Put \(h_x=\prod_yQ_y^{A_{xy}}\) for \(x\ne0\). From (5),

\[
\frac{z_x}{z_0}=\frac{h_x}{\overline{h_x}}.
\tag{9}
\]

If all points lie on an arc of length at most \(C\sqrt R\), choose
representatives \(\phi_x=\arg h_x\pmod\pi\) satisfying

\[
|\phi_x|\le\delta:=\frac{C}{2\sqrt R}.
\tag{10}
\]

This follows directly from their angular distances to the anchor; it
does not replace the actual phases by formal small real variables.

For a balanced column \(y\in B\), orthogonality of \(H\) gives the
integer vector identity

\[
4e_y=-\sum_{x\ne0}H_{xy}A_x.
\tag{11}
\]

Therefore

\[
Q_y^4=\prod_{x\ne0}h_x^{-H_{xy}}
\tag{12}
\]

in \(\mathbb Q(i)^\times\). If \(Q_y\) is a nonunit, (8), (10),
and (12) imply

\[
\frac1{|Q_y|}\le |\sin(4\arg Q_y)|
\le7\delta,
\qquad
\boxed{|Q_y|\ge\frac{2\sqrt R}{7C}.}
\tag{13}
\]

The equality in (12), and not merely a valuation divisibility, is why
the elementary fourth-power separation is applicable.

## 4. At most one unbalanced block can be nonunit

Take distinct \(y,z\in U\), and put
\(\varepsilon=S_y/S_z\in\{1,-1\}\). Another exact integer identity is

\[
2(e_y-\varepsilon e_z)
=-\sum_{x\ne0}\frac{H_{xy}-\varepsilon H_{xz}}2 A_x.
\tag{14}
\]

Its coefficients are \(0,1,-1\), with at most four nonzero terms.
Let

\[
\gamma_{yz}=\begin{cases}
Q_y\overline{Q_z},&\varepsilon=1,\\
Q_yQ_z,&\varepsilon=-1.
\end{cases}
\tag{15}
\]

Its squared phase is the phase of the monomial in (14). By (4),
\(\gamma_{yz}\) is conjugate-coprime and is a nonunit whenever at least
one of its two blocks is a nonunit. Equations (7), (10), and (14) give

\[
\frac{\sqrt2}{|Q_yQ_z|}
\le |\sin(2\arg\gamma_{yz})|\le4\delta,
\qquad
\boxed{|Q_yQ_z|\ge\frac{\sqrt R}{\sqrt2 C}.}
\tag{16}
\]

Suppose at least two of the four \(U\)-blocks are nonunits. Partition
the four columns into two disjoint pairs, each containing a nonunit.
When unit columns occur, give the nonunits unit partners as needed;
when all four are nonunits, any two disjoint pairs work. Multiplying
(16) for the two pairs yields

\[
\prod_{y\in U}|Q_y|\ge\frac{R}{2C^2}.
\tag{17}
\]

On the other hand the left side is at most
\(\prod_y|Q_y|\le R\), since \(|d|\ge1\). For \(C=1/2\), (17)
would give \(R\ge2R\), a contradiction. Thus at most one unbalanced
block is a nonunit. This step in fact holds for every \(C<1/\sqrt2\).

## 5. Distinctness forces three balanced nonunits, then a radius contradiction

Two rows agreeing on all nonunit block coordinates give points whose
ratio is a Gaussian unit. In this model the unit ratio is in
\(\{1,-1\}\). A containing arc of angular width at most
\(C/\sqrt R\le1/2\) cannot contain an antipodal pair, so distinctness
forces distinct patterns on the nonunit coordinates.

If at most two balanced blocks and at most one unbalanced block were
nonunits, there would be at most six distinct patterns. With no unbalanced
coordinate there are at most four. With one unbalanced coordinate, its
two values occur on two and six rows, while two additional binary
coordinates can distinguish at most four rows within each value; the
total is at most \(2+4=6\).

Eight distinct points therefore require at least three balanced nonunit
blocks. Multiplying (13) over three such blocks and using
\(\prod_y|Q_y|\le R\) gives

\[
R\ge\left(\frac{2\sqrt R}{7C}\right)^3,
\qquad
R\le\left(\frac{7C}{2}\right)^6.
\tag{18}
\]

At \(C=1/2\), this is

\[
R\le(7/4)^6=117649/4096<29.
\tag{19}
\]

Distinct equal-norm Gaussian integers have squared distances that are
positive even integers, hence distances at least \(\sqrt2\). The seven
successive chords of the eight points along the containing arc have
total length at least \(7\sqrt2\). Thus

\[
\tfrac12\sqrt R\ge7\sqrt2,
\qquad R\ge392,
\tag{20}
\]

contradicting (19). This proves the theorem, including the optional common
Gaussian integer multiplier in (5).

## 6. Why small row-dependent twists remain a real issue

If the row numerators are instead \(K_xh_x\), equations (12) and (14)
gain products of the \(K_x\). Their small phase bounds then concern a
twisted fourth or second power, not the pure powers in (7) and (8).
Even a fixed conjugate-coprime twist can change the relevant exponent.

For example, with \(\alpha=2+i\),

\[
\operatorname{Im}\bigl(\alpha(a+bi)^2\bigr)
=a^2+4ab-b^2=(a+2b)^2-5b^2.
\tag{21}
\]

The Pell solutions \(x+b\sqrt5=(9+4\sqrt5)^n\), with \(a=x-2b\),
make (21) identically one while \(|a+bi|\to\infty\). These Gaussian
integers have coprime coordinates and odd norms: the Pell equation
forces coordinate coprimality, and this family has \(x\) odd and \(b\)
even. Thus the twist does not merely introduce a nonprimitive or
ramified exception. Its normalized imaginary part is of order
\(|a+bi|^{-2}\), whereas the pure-square bound (7) has order
\(|a+bi|^{-1}\).

This is not an eight-point counterexample. It identifies the exact step
that cannot be transferred to arbitrary small corrections without a new
simultaneous argument. Nor is the full nearly uniform 128-cut profile a
small perturbation of this eight-cut support. The theorem gives a genuine
integer exclusion of one critical mixed profile, with the remaining
general finite inequality (1) still unproved.

## 7. Every normalized order-eight Hadamard matrix

Now let \(H\) be any \(8\times8\) matrix with entries \(\pm1\),
\(HH^{\mathsf T}=8I\), and anchor row all ones. Define \(A\), the
Gaussian blocks, and the points exactly as in (2), (4), and (5).
No Walsh representation or equivalence classification is assumed.

Let \(S_y\) be its column sums. Any two columns disagree in exactly
four positions. Their sums therefore differ by twice a sum of four
signs, hence by a multiple of four. All column sums have the same
residue modulo four. Also

\[
\sum_y S_y^2=\|H^{\mathsf T}\mathbf1\|^2=64.
\tag{22}
\]

If their common residue is zero, the possible absolute values are
\(0,4,8\), and (22) forces either one absolute value eight and seven
zeros, or four absolute values four and four zeros. If their common
residue is two, the possible absolute values are two and six; writing
\(k\) for the number of sixes, (22) becomes
\(36k+4(8-k)=64\), so \(k=1\). Thus the only absolute spectra are

\[
(8,0,0,0,0,0,0,0),\quad
(4,4,4,4,0,0,0,0),\quad
(6,2,2,2,2,2,2,2).
\tag{23}
\]

Both phase identities have the same general form. For \(S_y=0\),
equation (11) follows from column orthogonality. For nonzero
\(|S_y|=|S_z|\), equation (14) holds with
\(\varepsilon=S_y/S_z\). Its coefficient support is still at most
four, because the two columns agree or disagree in four positions.
Consequently (13) and (16) hold unchanged.

In particular, in any group of at least four columns with the same
nonzero absolute sum, **at most one block is a nonunit**. If there
were two, choose four columns containing them and pair the four so
each pair contains a nonunit. Applying (16) twice would give the same
contradiction \(R\ge2R\) as before.

Apply these facts to the three cases in (23):

* In the first case the absolute-eight column is constant. Because
  its anchor entry is one, the column is all ones, contributes no
  row incidence, and can be absorbed into the common factor. Eight
  distinct points require at least three nonunit balanced columns.
  Equations (18)–(20) give the radius contradiction.
* The second case is exactly the counting argument in Sections 4–5:
  the four nonzero-sum columns have at most one nonunit; their cut
  sizes are two or six. Distinctness again forces three balanced
  nonunits, and the same radius contradiction follows.
* In the third case the seven absolute-two columns have at most one
  nonunit. The single absolute-six column can supply at most one
  additional nonunit coordinate. Therefore at most four incidence
  patterns, and hence at most four points in the short arc, are
  possible. Eight distinct points are excluded directly.

This proves the theorem for every normalized order-eight Hadamard cut
system. The extension is entirely elementary and retains the actual
integer blocks and arc geometry. The common-unit and no-independent-twist
qualifications remain unchanged.

The companion [exact check](check_semibent_eight_point_phase.py) verifies
the projection identities and support argument for the displayed
semibent matrix, and also checks all 128 row-sign switchings of the
normalized Walsh matrix. The three spectrum counts are respectively
8, 56, and 64. That enumeration illustrates the general proof; equation
(22) and the modulo-four argument establish the result for arbitrary
Hadamard matrices.
