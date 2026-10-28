# referee notes (width.md)
- Lemma 2.1-2.4, Thm 2.5 re-derived by hand: all correct (see report §2).
- Possible error: §3.2 / §2.4 remark "slicing S_m by x_1 with C_bal(4) bound (3/2)(2m-|c|) and peeling returns exactly 15m/4".
  Hand computation: peel = max_k (3/2)(2m-k) + 2k = 3m + k/2 -> 7m/2 at k=m. 15m/4 arises with Thm B's 4-dim bound (7/4)(a+b). To verify by code.
- Thm 3.3 last claim: each semigroup element is a Laurent poly so the bound applies elementwise; fine.
- Twisted form: t_max over Q-sections = reg by descent (toric.md §3 says so). fine.
