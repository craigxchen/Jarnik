# Round 8: growth rate (drafts, refereeing in progress)

Both write-ups here are unrefereed drafts; their referees are running. Only the growth rate matters in
this round: the question is whether any input gives `M = o(log R/loglog R)`.

* `lw.md`: single-form Lang-Waldschmidt gives no growth improvement (draft claim). A linear form in
  logarithms sees a character only through its support weight. The Paley skew double core of
  `round5/exact.md` keeps every character's support weight at least `W/4 + (1/2) log P - 1/16`, so no
  single LW inequality excludes it. For each fixed `eps` there are LW-consistent fake circles at
  `W = O(M log M)`. Open: the full conjecture with `eps -> 0`. The draft confirms that
  `round4/flipped.md` Theorem L excludes the genuine versions of the round-5/6 fakes.
* `weakest.md` (draft claims):
  * Rank reduction, unconditional: `M <= (5.5^k + 1)/2` for a cluster whose exponent lattice has rank
    `k`. So the uniform theorem is equivalent to bounded rank.
  * New simultaneous Lang-Waldschmidt hypotheses: five forms imply the uniform theorem, and
    logarithmically many forms imply `M = O((log R)^(1-eta))`.
  * Single-form LW, abc, Bombieri-Lang and the quantitative Subspace Theorem with known
    `S`-dependence do not suffice.
