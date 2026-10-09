# Round 4 context (lead author, 2026-10-09)

Goal: an unconditional proof that arcs of length C sqrt(R) on x^2+y^2=R^2 hold <= M(C) lattice points,
or failing that a general growth improvement M = o(log R/loglog R). Neither is known.

Refereed barriers established this session (all in /home/user/Jarnik/docs/claude_review/, read the READMEs):
* Local inputs (pair chords, gcds, small-prime residue collisions, weighted cut/Plotkin bound) are capped
  at order log R/loglog R; best constant 2/3 (two_thirds/proof.md).
* Multi-row character certificates (items 311-335 of the notes; growth/walsh-extraction.md): uniform for
  pure Hadamard cores (M <= 24 at C <= 1), growth for Walsh with arbitrary flips; FAIL on fully flipped
  capacity-two Paley profiles, which pass every pair, integral/rational character, aggregate collision and
  capacity inequality at W = Theta(M log M). Missing: "fully flipped Hadamard phase lemma" (Sec. 10 there).
* Single-place auxiliary polynomials (threshold/upper.md): every form F on X_M = {x_j^2+y_j^2 equal} has
  ord along the diagonal <= (2 - 1/ceil(M/2)) deg F < 2 deg F; so Liouville and Ru-Vojta with S={inf}
  (beta_M < 2) cannot give exponent 1/2. Forced Gaussian divisibility only makes the Vandermonde exactly
  critical (threshold/construct.md); the uniform cube profile satisfies all such inequalities.
* Subspace/Ru-Vojta with S = primes of N: exceptional sets depend on S (research notes:
  claimed_complete_proof_salvage_audit.md).
* Ptolemy route (growth/ptolemy.md, threshold/L8.md): reduces to (L_8) on 2x8 integer matrices; (L_5),
  (L_6) false; no auxiliary polynomial on Mbar_{0,k} proves (L_k); (L_7) open (heuristically false).
* Conditional: Vojta (D=0) on Bl_Y X_7 implies the uniform bound (growth/conditional.md).
* Data: exhaustive to R^2 <= 1e12: max 3 points on arcs (1/2)sqrt R, 4 on sqrt R; heuristically
  M(1/2)=M(1)=5; a 10-point cluster exists at C = 8.24 (growth/families.md).
Research notes (1170 files): /tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/research/docs/
