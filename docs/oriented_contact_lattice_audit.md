# Audit of the oriented contact-lattice theorem

This is an independent arithmetic audit of
[`oriented_contact_lattice_rigidity.md`](oriented_contact_lattice_rigidity.md).
The exact CRT and shared-five-vector checks are in
[`check_oriented_contact_lattice_rigidity.py`](check_oriented_contact_lattice_rigidity.py).

Under the stated hypotheses—integer primitive coefficient vectors
\(\lambda_e=(r_e,s_e)\), odd pairwise coprime Gaussian moduli, and
\((G_e,\overline{G_f})=1\) for every pair—the existence argument is valid.
At each split prime power, \(\lambda_e\) is primitive over
\(\mathbf Z/p^k\mathbf Z\), so an \(SL_2\) matrix can send it to the
kernel direction \((-\iota,1)\) of \((1,\iota)\). CRT combines the local
matrices, and surjectivity of
\(SL_2(\mathbf Z)\to SL_2(\mathbf Z/D\mathbf Z)\) supplies an integral
determinant-one lift.

The parity normalization is also correct.  An odd Gaussian integer has
residue (1) or (i) modulo (2), so a unit makes a chosen generator
\(\Delta\) satisfy \(\Delta\equiv1\pmod2\).  Since every Gaussian integer
is congruent to its conjugate modulo (2),

\[
y=\frac{x+\Delta\overline{x}}2
\]

is integral.  The oddness of every contact modulus puts (y) in the same
contact lattice as (x).  Direct determinant and Hermitian calculations
give

\[
\det(x,y)=-i\Delta,
\qquad
\bigl(h(x_i,x_j)\bigr)=
\begin{pmatrix}2&1\\1&(1-N\Delta)/2\end{pmatrix}.
\]

Because conjugate-primitivity forces all rational prime factors of
\(N\Delta\) to be \(1\pmod4\), we have \(N\Delta\equiv1\pmod4\).
Consequently the equivalent norm equation with coefficient
\((1-N\Delta)/4\) is integral on
Gaussian coefficients.

The duality formula is consistent with the one-contact model
\(L=G\mathbf Z[i]\oplus\mathbf Z[i]\):

\[
L^\vee=(\overline\Delta)^{-1}\overline L,
\qquad
L^\vee/L\cong
\mathbf Z[i]/(\Delta)\oplus\mathbf Z[i]/(\overline\Delta).
\]

The second displayed quotient is a noncanonical decomposition; the two
ideals are coprime, so it is also cyclic as an abstract Gaussian module.
No local discriminant obstruction is being hidden here.

Finally, if two (q=1) vectors in (L) are independent, their determinant
is divisible by \(\Delta\), giving
\(\lVert x\rVert\lVert x'\rVert\geq|\Delta|\).  If they are dependent,
primitivity of a (q=1) vector forces the scalar to be a Gaussian unit.
Thus the short-lift uniqueness and orientation-overlap bounds follow with
the logarithmic constants stated in the source note.

As a finite shortest-vector diagnostic, the companion checker
[`check_contact_lattice_min_height.py`](check_contact_lattice_min_height.py)
uses four pairwise conjugate-primitive generators of norms (5,13,17,29)
and \(x=(1,i)\), so each contact equation is satisfied by taking
\(\lambda_e=(\Re G_e,\Im G_e)\).  Exhausting Gaussian basis coefficients
of size at most (12) finds exactly the four unit multiples of (x), with
minimum squared Euclidean norm (2).  This finite observation only checks
that example; it supplies no universal shortest-vector or height bound.

The hypotheses should remain explicit.  If moduli share a Gaussian prime,
share a conjugate prime, or if a coefficient vector is not primitive, the
index product, local CRT independence, duality decomposition, and clean
orientation conclusions need not hold.
