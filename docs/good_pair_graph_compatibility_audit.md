# Good-pair graph: the precise scope of binary valuations

For primitive fixed-unit source points of odd squared radius \(N\), the
ordered Gaussian pair factors satisfy
\[
 A_{ab}A_{bc}=q_{abc}A_{ac},\qquad
 q_{abc}b_{ac}=a_{ab}b_{bc}+a_{bc}b_{ab},
\]
where \(q_{abc}\mid N\) is an odd positive integer. Opposite coordinate
parity gives, for \(k_{ab}=v_2(|b_{ab}|)\),
\[
 k_{ab}\ne k_{bc}\Rightarrow k_{ac}=\min(k_{ab},k_{bc}),\qquad
 k_{ab}=k_{bc}\Rightarrow k_{ac}>k_{ab}.
\]
Thus the full pair set has a binary hierarchy. If every pair has
\(|b_{ab}|\le H\), the proved capacity bound is
\(M\le2^{\lfloor\log_2H\rfloor+1}\le2H\).

A positive-density selected set of edges does not inherit that conclusion
from the valuation hierarchy alone. Assign distinct integer labels
\(L_i=2i\) and \(R_j=2j+1\), and let each pair's valuation be the
2-adic valuation of the difference of its labels. These labels satisfy the
triangle rule exactly. Every cross edge has valuation zero; within either
side, the valuation is \(1+v_2(|i-i'|)\). Taking balanced side sizes gives
\(\lfloor M^2/4\rfloor\) cross edges, enough for the good-pair quota
\(\lceil M(M-2)/4\rceil\), but this selected graph has no triangle.

This is only a valuation model. It neither assigns the actual integer
\(b_{ab}\) nor realizes the real-part sizes, positive triangle equations,
or near-critical norm interval of a short circle arc. In particular it does
not rule out stronger extraction using those additional arithmetic and
angular hypotheses. It rules out an inference from the binary hierarchy
and the edge quota alone.

See [the exact triangle proof](pair_norm_triangle_compatibility.md) and
[the full-pair capacity calculation](binary_pair_residue_capacity.md).
The separate [Möbius transfer audit](mobius_determinant_free_transfer_audit.md)
uses explicit isotropic coordinates and supersedes an incorrect local
matrix calculation removed from an earlier draft of this note.
