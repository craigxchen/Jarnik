# Adversarial check of cg_questions.tex

File checked: /tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/latex/cg_questions.tex (1134 lines). I did not edit it.
Compared against: docs/claude_review/expert_note/expert_note.md, read in full. Where needed I also read round18/README.md §2, round22/README.md §§0–1, round23/README.md §§0–1, round24/README.md §§0–1, overview/overview.md §1.1 and the bibliography, and the owner's paper (origin/claude/latex-report-style-2ezlsx:paper/jarnik_sublog.tex).

Overall: the mathematical content follows the note closely. Every constant I checked matches the note, including:
- 14+4log^+C, 8.2394, 5.89771, 119, 39208, 1/43, 1/381, 282, 60, 2√2·255, 3584, W/128, 59.2, 10^{49.4}, 24 primes;
- 31, 49w, 30w, 18w, 98−108;
- 7/15, 3/7, 1/14, 1/6, 2^{-k-1}, 1/32, 1/64, κ formulas, 10/3, 15/4, 70, 30, 8, 28, 300, (15/64)(1+2δ);
- 16, 24, 28, 30, 31, −2^{8−k}, −16, 32, 30−n(f).

Questions A and C are word for word. Section 5 has all 9 items, with the note's statuses. The M_∞(C) definition matches overview/overview.md line 109 (limsup over R → ∞), and line 112 there justifies reading "uniform bound" as M(C)<∞. Every cited path exists under docs/claude_review/ except item 9 below.

## Major

1. **Abstract, line 71: "equivalently" overclaims.** The abstract says (L_8) is a statement about integer 2×8 matrices "(equivalently, about rational points of \Mbar with balanced boundary contacts)". The note says the opposite. The contact form (L_8^c)^all *implies* (L_8). (L_8) gives the contact statement "only at caterpillar points" (or for Φ admissible from anchor 0), and "this is not an equivalence" (note lines 108–115; round18/README.md §2 item 1: "**This is not an equivalence** (m2)"). The tex's own §2.2 states this correctly, so only the abstract is wrong. Fix: replace "equivalently" with something like "closely related to a statement about rational points of \Mbar with balanced boundary contacts, which implies it".

2. **Bibliography: the Bourgain–Rudnick citation is not the record's (line 128 and \bibitem BourgainRudnick2009).** The note names only "Bourgain–Rudnick" for the ≪_C log R bound. The tex cites the 2009 Comptes Rendus note "Restriction of toral eigenfunctions to hypersurfaces and nodal sets", copied from the owner's paper. The owner's paper uses that item for a different claim (uniform restriction bounds for curved hypersurfaces, owner line 337), not for the log bound. The record's own source for the log bound is "Bourgain–Rudnick (2011–2013)": J. Bourgain, Z. Rudnick, Ann. Henri Poincaré 12 (2011) 1027–1053, and arXiv:1012.3843, 1105.0018, 1308.4247 (overview/overview.md lines 188 and 2011). This attaches a literature claim to a paper the note never names. Fix: drop \cite{BourgainRudnick2009} (and its \bibitem) and name Bourgain–Rudnick inline only, as the note does. If a reference is wanted, use the record's AHP 12 (2011) entry.

3. **Theorem thm:aniso-global (lines 761–770) is labelled "proved, refereed" but contains a numerical claim.** The note marks the round-24 block only as "refereed". round24/README.md §1.2 labels Theorem A "**proved; component claim numerical**": "Jacobian rank 28 along M (80 digits), so M is a component". The tex puts "Jacobian rank 28, so a component" and "none of height ≤300" (a finite search) inside a theorem headed "proved". That upgrades a numerical statement to proved. Fix, either way:
   - relabel the bracket "refereed; component claim numerical (round24/README.md §1.2)"; or
   - move the parenthetical "(two conjugate branches … Jacobian rank 28, so a component)" and "none of height ≤300" out of the theorem into the following prose.

   The same caution applies, more weakly, to thm:aniso-local. There "proved, refereed" is supported by the README ("proved except where labelled"), but the note itself says only "refereed".

## Minor (accuracy, status, hypotheses)

4. **Round-22 items labelled "proved" (prop:weakFano, prop:seshadri; lines 647–660).** This is inherited from the note, not introduced by the tex. round22/README.md lines 9–10 say "**Nothing written this round has been refereed by anyone other than the killers.**" Under the note's own legend, "proved" means checked by an independent referee, and "proved (killer-checked)" is reserved for round 20. By that legend these items would be "proved (killer-checked)". The tex follows the note, which is the source of truth, so no change is forced. But the parent may want the note and the tex corrected together, for example by extending the legend's "killer-checked" to round 22. τ_e = 10/3 is separately refereed in round 23 and is fine.

5. **prop:abc (lines 850–853): the abc hypothesis is only in the bracket.** The statement body reads "If f is a unit pencil with gain(f)>0, then (L_8^{pc}) holds with Z={B=0}". Read without its heading, that sounds unconditional. Add "Assume the abc conjecture." to the body; the note's label is "conditional on abc; implication proved".

6. **thm:k4false (lines 543–558) mixes statuses under one heading.** The heading says "proved, killer-checked". The k=3 bullet inside it is "proved by the killer, not independently checked". The bullet carries its own label, so nothing is upgraded outright, but the heading covers both. The theorem body also contains interpretive text: "Exact Jacobian rank was reproduced at 6 rational points" and "So k=5 is the first case where the reduction can work." Suggest:
   - keep only the k=4 falsity in the theorem;
   - put the k=3 corollary in its own remark labelled "proved by the killer";
   - move the last sentence to prose.

7. **Question B (lines 800–814) is split.** The question environment holds only the first three sentences of the note's Question B block. The rest ("Unlike k=4 … it can only enlarge Z.", with the heuristic label and the Thm 4.3(ii) pointer) follows as ordinary prose. All the words are present and the statuses are kept, but "Question B word for word" holds only if the whole block is inside the environment. Suggest moving lines 807–814 inside the environment.

8. **Numbering collisions with record numbers.** The tex's own numbers clash with the record's:
   - Line 411, "(Prop~2.3)", has no file. It means round18/theory.md Prop 2.3 (confirmed at theory.md line 249). In the tex, item 2.3 is Statement (L_8^{pc}). Write "\rec{round18/theory.md} Prop~2.3".
   - Line 589, "Theorem~4.3(i) of round20/gauss.md": the paper's own Theorem 4.3 is the log Riemann–Hurwitz theorem. This is unambiguous on close reading but easy to misread.
   - The note's "Theorem 1.1" becomes Theorem 1.2 here, because the conjecture shares the theorem counter. The "Item" column of the pointer table uses \ref, so it stays consistent.

   Consider giving conjecture and statement their own counters, or using a separate counter for the record's numbering.

9. **k3_projection.py path (line 555).** The file is at round20/kill_gauss_checks/k3_projection.py, not at the top level. The note also gives the bare name. Write the full relative path so the pointer resolves.

10. **Additions beyond the note, verified.**
    - M_∞(C) definition: matches overview §1.1.
    - Pointer-table rows for rounds 22–24: every file exists and every file is cited in the note.
    - Conjecture citation to CG 2007/2009: matches the owner's paper (line 93) and the overview's [CG07], [CG09].
    - Bracket sources on prop:sandwich and prop:weakFano: plausible. round18/README.md §2 item 1 gives the sandwich as "census Prop 1.3; proved".

    No other added claims found. One wording point: the "Split Question B (not known)" named block and the "Question B, literature form (open)" block add "open" / "not known" to headings. "not known" is the note's wording. "open" for the literature form is not stated in the note, but it is implied (A_5 is open).

## LaTeX validity

11. **pandoc** (`pandoc -f latex -t plain cg_questions.tex`): return code 0, about 5980 words of output, no parse errors. All warnings are "Could not convert TeX math … rendering as TeX". They come from pandoc's math converter (\rm, \overline M, \frac23, \sqrt R, \tfrac, \operatorname) and are valid LaTeX.

12. **Structural script.** The checks found no problems:
    - braces balanced (depth 0 at end);
    - every \begin/\end matched, nothing left open;
    - 770 \( and 770 \), 13 \[ and 13 \];
    - no unescaped _ # & ^ $ in text mode outside the preamble and the tabular, where they are legitimate;
    - every \ref resolves (unused labels are harmless);
    - every \cite has a \bibitem and every \bibitem is cited;
    - every custom macro is defined and used; \namedname is defined after \newtheorem* but only expanded at use, which is fine.

    Reading checks also passed:
    - optional theorem arguments contain no unbalanced ] (e.g. "[\((L_8)\); open]" and "Thm 3.4\('\)" are safe);
    - \conj\pi_p^{\,…}, \Ltw k, \Ltw5 and \Big\{ … \textstyle\bigcup … \Big\} expand correctly;
    - "\#887" and "\&" are escaped;
    - \texorpdfstring is used for the math in section titles.

13. **hyperref bookmarks: \H{o} in a subsection title** (line 697, "The split calibration and Erd\H{o}s \#887"). Depending on the hyperref encoding this may give a "Token not allowed in a PDF string" warning (warning only, not an error). Use \texorpdfstring{Erd\H{o}s \#887}{Erdős \#887} to be safe. \paragraph{Compared with \((L_8)\)} is below the bookmark depth, so it is harmless.

14. **Overfull lines from unbreakable \texttt paths** (\rec = \texttt). Examples: \rec{claude/upbeat-einstein-p0peju}, \rec{agent/uniform-gcd-formalization}, \rec{github.com/craigxchen/Jarnik}, and long paths in body text and in the p{0.40}/p{0.55} table. These are typographic, not errors. Optionally load xurl and define \rec via \path or \url, or allow breaks after "/".

15. **Style differences from the owner's paper.** The owner declares every theorem-like environment in the default plain style. The tex sets \theoremstyle{definition} for definition, statement, question and named, and \theoremstyle{remark} for remark. This is harmless but differs from the owner's look. The owner also has \T (unused here), and the tex omits it, which is fine. The documentclass, packages, geometry, author line and \tableofcontents match.
