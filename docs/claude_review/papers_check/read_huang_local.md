# Huang, local distribution on toric surfaces I and II: full-text reading against Question B

Reader: subagent (papers_check). Both PDFs were read in full (I: pp. 1-64; II: pp. 1-48) with the Read tool.
Page numbers are the printed page numbers, which coincide with the PDF page indices in both files.
Labels: **[quote]** = text as printed (French originals, with my translation where useful); **[paraphrase]** = close
restatement of a printed statement; **[inference]** = my own reasoning, not in the papers. Nothing in the PDFs was
treated as an instruction. The third uploaded Huang paper (arXiv 1809.02001, ANT) was **not** part of this assignment
and was not read here.

---

## 1. Bibliographic data as printed in the PDFs

**Part I.** arXiv:1703.01772v2 [math.NT], margin stamp "17 Dec 2018".
Title: *Approximation diophantienne et distribution locale sur une surface torique*. Author: Zhizhong Huang,
"Institut Fourier, UMR 5582, Université Grenoble Alpes, 38610 Gières, France", zhizhong.huang@yahoo.com. "Current
address: Institut für Algebra, Zahlentheorie und Diskrete Mathematik, Leibniz Universität Hannover" (p. 64). 64 pages.
French résumé plus an English abstract (p. 1). **No journal reference, volume or page range is printed anywhere in this
PDF.**

**Part II.** arXiv:1805.03920v2 [math.NT], margin stamp "7 May 2019".
Title: *Approximation diophantienne et distribution locale sur une surface torique II*. Author: Zhizhong Huang,
"Institut für Algebra, Zahlentheorie und Diskrete Mathematik, Leibniz Universität Hannover, Welfengarten 1, 30167
Hannover" (p. 1). 48 pages. **No journal reference is printed.** Its bibliography cites Part I as
"[19] Zhizhong Huang. Approximation diophantienne et distribution locale sur une surface torique. *à paraître dans Acta
Arith.*, 2019." (p. 47). The volume and page range that the record gives ("Acta Arith. 189 (2019) 1-94") **do not
appear in either PDF**.

Related items cited inside the PDFs: [15] in I = [18] in II: Huang, *Distribution locale des points rationnels de
hauteur bornée sur une surface de del Pezzo de degré 6*, IJNT 13 (2017), no. 07, 1895-1930. Pagelot's essential
constant: I [22] = II [28]: S. Pagelot, *Étude de la distribution asymptotique fine des points rationnels de hauteur
bornée*, unpublished manuscript (II adds "Université de Grenoble I, 2008"). McKinnon-Roth: I [21] = II [27], Invent.
Math. 200 (2015) 513-583.

---

## 2. Definitions (exact statements, with pages)

### 2.1 Part I, §2.1 (p. 7)

Setting **[quote, p. 7]**: "Soit X une variété projective définie sur un corps de nombres k. On fixe un point rationnel
Q ∈ X(k̄), ν une place réelle et d = d_ν une distance projective […]. Soit L un fibré en droites sur X muni d'une hauteur
de Weil absolue H_L. Pour V une partie constructible de X, on considère les ensembles

A(Q,V) = {γ > 0 | ∃(y_i) ∈ (V(k) \ {Q})^N, ∃C > 0, d(Q,y_i) → 0 et d(Q,y_i)^γ H_L(y_i) < C, ∀i},
B(Q,V) = {γ > 0 | ∃C > 0, d(Q,y)^γ H_L(y) ≥ C, ∀y ∈ V(k) \ {Q}}."

It adds: "On suppose que k(Q) le corps résiduel de Q est contenu dans k_ν". The target point is therefore allowed to be
**algebraic** (Q ∈ X(k̄) with real residue field). It is called "rationnel" throughout, but Thm 1.3 applies it to
θ ∈ P^1(Q̄) ∩ R.

* **Def. 2.1 (Northcott, p. 7).** L has the Northcott property for an open U if #{x ∈ U(k) : H_L(x) ≤ B} < ∞ for all B.
* **Prop.-Déf. 2.2 (p. 7, citing [21] Prop. 2.11).** "Soient Q ∈ X(k̄) et V une partie constructible de X. Supposons
  que L vérifie la propriété de Northcott pour un ouvert de Zariski U contenant Q. Alors on a
  inf A(Q, U ∩ V) = sup B(Q, U ∩ V). Cette quantité est appelée *constante d'approximation en Q dans V*, notée
  α_L(Q,V), ou simplement α(Q,V)". So α is the McKinnon-Roth constant; **small α means well approximable**. This is
  the record's convention (dist^γ H bounded).
* **Déf. 2.3 (p. 7, credited to Pagelot [22]).** "On définit *la constante essentielle* en Q comme la quantité
  α_ess,L(Q) = α_ess(Q) = sup_V α(Q,V) (2.1) où V parcourt toutes les parties constructibles denses de X. S'il existe
  une sous-variété Z de X de sorte que pour tout ouvert dense W ⊂ Z, α(Q,W) < α_ess(Q), on dit que Z est *localement
  accumulatrice*." Translation: α_ess(Q) is the sup of α(Q,V) over dense constructible V; Z is *locally
  accumulating* if α(Q,W) < α_ess(Q) for every dense open W ⊂ Z.
  [inference] Because a dense constructible set contains a dense open set, α_ess(Q) > γ holds iff some dense
  Zariski-open U has α(Q,U) > γ, i.e. iff Zariski-dense sequences cannot reach γ. This is exactly the record's
  "α^ess_e(-K) > 15/4 off a proper closed set".
* p. 7: "on suppose que le fibré en droites L est gros" (big), with U ≅ its image under the rational map given by
  H^0(X, L^m), so Northcott holds on U. "ce qui est le cas pour le fibré anticanonique des variétés toriques complètes
  lisses […] En outre le fibré −K_X est sans lieu de base sur l'orbite ouverte et donc la propriété de Northcott est
  vérifiée en tout point de l'ouvert."
* **Thm 2.4 (p. 7; "Théorème de Roth, Principe de Dirichlet, [21] Lemma 2.15").** For x ∈ P^1(Q̄) ∩ P^1(R) and
  O(d): α(x,P^1) = d if x ∈ P^1(Q), and d/2 if x is irrational.
* **Remarque 2.5 (p. 8, quoting [21] Theorem 2.16).** "Prenons une courbe rationnelle C définie sur k et L un faisceau
  inversible ample sur C. On fixe Q ∈ C(k̄). Soit f : P^1 → C le morphisme de normalisation. Alors
  α(Q,C) = min_{P ∈ f^{-1}(Q)} d/(m_P r_P), où d = deg_C(L), m_P est la multiplicité de la branche de C passant par Q
  correspondant à P et r_P = 0 si k(P) ⊄ k_ν ; 1 si k(P) = k ; 2 sinon." For a section of O(2,2) on P^1 × P^1 that
  is nodal at a rational Q with branch slopes τ_i(Q): α(Q,C) = d if k(Q_i) = Q; d/2 if Q ⊊ k(Q_i) ⊂ R; ∞ otherwise.
* **Conjecture 2.6 (McKinnon, p. 8).** "Soient L un fibré ample sur X avec une hauteur de Weil choisie et Q ∈ X(k).
  Supposons qu'il y a une courbe rationnelle passant par Q (ceci implique en particulier que α(Q,X) < ∞). Alors il
  existe une courbe rationnelle C telle que α(Q,C) = α(Q,X)."
* p. 8, the method used to compute α_ess **[quote]**: "démontrer une borne inférieure uniforme valide dans un ouvert
  pour le produit d'une certaine puissance de la distance avec la hauteur et trouver une famille de courbes rationnelles
  dont la réunion dans la variété est dense pour la topologie de Zariski et chacune atteint la constante
  d'approximation souhaitée. […] cette procédure donne seulement une condition « suffisante » pour trouver α_ess. En
  général, il n'est pas clair que α_ess puisse toujours être atteinte sur une famille de courbes rationnelles."

### 2.2 Zoom (Part I §2.2, pp. 8-9)

* "l'opérateur de zoom Ψ_{r,B} avec le facteur r > 0 est défini par Ψ_{r,B}(P) = B^{1/r} ρ(P)" (ρ a local
  diffeomorphism X(R) → T_Q X sending Q to 0). "avec cette convention, *plus r est grand, plus le zoom est faible*."
  The measure is δ_{U,Q,B,r}(f) = Σ_{x ∈ U(k): H_L(x) ≤ B} f(B^{1/r} ρ(x)) for f continuous with compact support on
  (T_Q X)_R.
  [inference] With f = χ(ε) (the ball of radius ε), δ_{U,Q,B,r}(χ(ε)) counts the x ∈ U(k) with H(x) ≤ B and
  d(Q,x) ≤ εB^{-1/r}; all of them satisfy d^r H ≤ ε^r. **The zoom factor r is an approximation exponent.** The critical
  zoom r = α_ess is the threshold scale of the essential constant, and the "natural" (Batyrev-Manin) scale is r = dim X.
* **Heuristic (2.3) (p. 9)**, for L = ω_X^{-1}, assuming Batyrev-Manin on U and "une forme forte de
  l'équidistribution […] pour des voisinages réels dont la taille dépend de B": δ_{U,Q,B,r}(χ(ε)) ≫≪_ε
  B^{1-n/r}(log B)^{κ-1}, hence "β = 1 − n/r, γ = κ − 1". "En fait, dans les (rares) exemples connus, la valeur de β est
  en accord avec cette prédiction, mais celle de γ ne l'est pas."
* **Prop. 2.7 (p. 9).** "Soit U une partie constructible dense de X pour laquelle L vérifie la propriété de Northcott
  telle que α(Q,U) = α_ess(Q). Alors pour tout 0 < r < α_ess et f ∈ C^b_Q(X), on a 1. δ_{U,Q,B,r}(f) = f(Q) si
  Q ∈ U(k̄) ; 2. δ_{U,Q,B,r}(f) = 0 sinon, pour B ≫_{f,r} 1." (The proof is the Liouville-type definition of α(Q,U).)
* **Déf. 2.8 (p. 9).** The zoom is *critique* (resp. *sous-critique*) if r = α_ess (resp. r > α_ess).
* **Déf. 2.9 (p. 9).** "Supposons que α_ess(Q) < ∞. On dit qu'une sous-variété fermée stricte W de X est *localement
  faiblement accumulatrice* si α_ess(Q) = α(Q,W) et s'il existe U_1 ouvert dense de X tel que pour tout ε > 0
  suffisamment grand et pour tout ouvert U_2 dense de X satisfaisant à α(Q,U_2) = α_ess(Q) et U_1 ∩ U_2 ∩ W dense dans
  W, on ait, en notant U = U_1 ∩ U_2, δ_{U∖W,Q,B,α_ess(Q)}(χ(ε)) = o(δ_{U∩W,Q,B,α_ess(Q)}(χ(ε)))." Translation: W
  has the same approximation constant as α_ess but dominates the count at the critical zoom.

### 2.3 Part II, §2 (pp. 9-12)

* **Déf. 2.1 (McKinnon-Roth, p. 10).** X is projective irreducible over **Q**, "Q_0 ∈ X(Q̄) un point algébrique réel
  fixé et L un fibré en droites gros tel qu'il induise une application birationnelle de X sur son image et bien définie
  sur un ouvert U_0 contenant Q_0". "La constante d'approximation α(Q_0,V) est le supremum des δ > 0 tels que
  l'inégalité du type de Liouville suivante ∃C(δ) > 0, ∀y ∈ (V ∩ U_0)(Q) ∖ {Q_0}, d(Q_0,y)^δ H(y) > C(δ), soit
  valide." (Here α is defined directly as the Liouville-type sup, i.e. sup B(Q,V) of Part I.)
* **Déf. 2.2 (Pagelot, p. 10).** "La constante essentielle (par rapport à V) est α_ess(Q_0,V) = sup α(Q_0,Y), Y ⊆ V
  partie constructible dense pour la topologie de Zariski induite. Par convention, α_ess(Q_0) = α_ess(Q_0,X). V est
  dite localement accumulatrice si elle vérifie α_ess(Q_0,V) < α_ess(Q_0)." (This is a slightly different, essentially
  equivalent, form of Part I's Def. 2.3.)
* **Thm 2.3 (McKinnon-Roth [27] Thm 2.16, pp. 10-11).** "Soient C une courbe rationnelle définie sur Q et L un faisceau
  inversible ample sur C et Q_0 ∈ C(Q̄). Soit φ : P^1 → C le morphisme de normalisation. Alors
  α(Q_0,C) = α_ess(Q_0,C) = min_{P ∈ φ^{-1}(Q_0)} d/(m_P r_P), où d = deg_C(L), m_P est la multiplicité de la branche
  de C passant par Q_0 correspondant à P et r(P) = 0 si k(P) ⊄ R ; 1 si k(P) = Q ; 2 sinon".
* **Zoom (2.1), p. 11:** δ_{Y,Q_0,B,r}(f) = Σ_{x ∈ Y: H_L(x) ≤ B} f(B^{1/r} ρ(x)). "Si r < α(U,Q_0), alors on a
  δ_{U,Q_0,B,r} → δ_{Q_0} […]. Le cas le plus intéressant est quand U est une partie dense de X telle que
  α(Q_0,U) = α_ess(Q_0) < ∞ et que l'on prend le facteur de zoom r = α_ess(Q_0), appelé critique."
* **Thin sets, Déf. 6.1 (p. 21, Serre [33] Def. 3.1.1)** and p. 22: type I = contained in a proper Zariski-closed set;
  type II = image f(V(Q)) for f : V → X generically finite and dominant of degree ≥ 2 with V integral and no rational
  section. "celle du type II est parfois dense dans X pour la topologie de Zariski."
* **l(X) (p. 3):** "Nous désignons par A_1(X) le groupe de Chow des 1-cycles modulo l'équivalence algébrique, et par
  l(X) le rang du sous-groupe de A_1(X) engendré par les classes des courbes rationnelles C vérifiant
  α(Q,C) = α_ess(Q)."

### 2.4 Heights, distances and fields in the examples

* Everything is over **Q**, and the varieties are **split** toric surfaces (blow-ups of P^1 × P^1 at torus-fixed
  points). The target is Q = [1:1] × [1:1], a **rational** point of the open orbit. The place is the real place.
* **Y_4** (I, §5.1, p. 36): H_{ω^{-1}}(P) = max(|x²st|, |y²st|, |t²xy|, |s²xy|) / (pgcd(x,s) pgcd(x,t) pgcd(y,s)
  pgcd(y,t)) for x,y,s,t ∈ Z with pgcd(x,y) = pgcd(s,t) = 1. The local chart is (w,z) = (y/x, t/s) and
  ρ(P) = (w−1, z−1) (5.3). The distance is d(P) = max(|z−1|, |w−1|) (5.7).
* **Y_3** (II, §3.2, p. 14): H(P) = max_{f ∈ S}|f(P)| / pgcd(f(P), f ∈ S) with
  S = {x²st, y²st, t²xy, s²xy, xyst, y²t²} (3.2). The local chart is ρ = (w,z) = (x/y − 1, s/t − 1) (3.3) and the
  distance is d((w,z)) = max(|w|,|z|) (3.4).
* **P^1** (I, (4.1), p. 13): H([u:v]) = max(|u|,|v|), pgcd(u,v) = 1, so this is O(1). With this normalisation
  α(θ,P^1) = 1/2 for a real quadratic θ, which is the O(1)-normalised version of "Liouville exponent 2".

---

## 3. Theorems (exact statements, with pages)

### 3.1 Part I: the surface Y_4 (blow-up of P^1 × P^1 at its 4 torus-invariant points)

**Thm 1.2 (p. 3; proved as Thms 5.4, 5.5, 5.6) [quote].** "On a α(Q,Y_4) = α_ess(Q) = 2. Soit U = Y_4 ∖ ∪_{i=1}^4 Z_i.
Soit r ≥ 2. Alors
1. Si r = 2, on a que pour toute f ∈ C^b_Q(X), δ_{U,Q,B,r}(f) = O_f(1). Il n'existe pas de mesure limite dans ce cas.
2. Si 2 < r < 144/55, alors il existe une constante η ∈ ]0,1[ telle que pour toute positive f ∈ C^b_Q(X),
   δ_{U,Q,B,r}(f) ≫ (∫ f dδ_r) B^{(1+η)(1/2 − 1/r)}(log B)^3 (1.5) où δ_r est une mesure absolument continue par
   rapport à la mesure de Lebesgue planaire."

Here the Z_i are the 4 rational curves through Q of minimal anticanonical degree 2, with equations x = y, s = t,
xs = yt, xt = ys (5.1) (p. 35). In the chart they are z = 1, w = 1, z = w, zw = 1 (5.4) (p. 36). On p. 5 the paper says
the lower bound (1.5) is proved "uniquement pour tous les η < 1/35 (cf. Théorème 5.6)", whereas the naive
Batyrev-Manin heuristic would give the power B^{2−1/r}… ("on souhaiterait que l'on voudrait pouvoir prendre η = 1").

* **Prop. 5.1 (p. 37).** "α(Q,Y_4) = 2." The proof, for P = [x:y] × [s:t] ≠ Q with t ≠ s, is
  H_{ω^{-1}}(P) d(P)² ≥ t²xy/(pgcd…) · (t/s − 1)² = xy t²(t−s)²/(pgcd(x,s)pgcd(x,t)pgcd(y,s)pgcd(y,t) s²) ≥ 1 (5.8).
  The remark after it: "Cette majoration uniforme indique un phénomène de « trou » […] Ceci peut aussi s'interpréter
  comme une inégalité du type Liouville : d(P) ≥ 1/H_{ω^{-1}}(P)^{1/2}, ∀P ∈ Y_4(Q) ∖ {Q}." (The same remark prints
  "δ_{U,Q,B,1/2}(χ(ε)) = 0", apparently meaning zoom factor 2.)
* **Geometry (p. 36).** "Le diviseur anticanonique ω^{-1}_{Y_4} = O(2,2) − E_1 − E_2 − E_3 − E_4 […] est gros et
  engendré par ses sections globales mais il n'est pas ample (car, les fibrés L_1 = O(1,0) − E_3 − E_4,
  L_2 = O(0,1) − E_1 − E_2 représentent des courbes effectives et ⟨ω^{-1}, L_1⟩ = ⟨ω^{-1}, L_2⟩ = 0)". So Y_4 is a
  **weak** del Pezzo surface of degree 4: −K is nef and big but not ample. §6 (p. 49) identifies it as the minimal
  toric desingularisation f : Y_4 → V of the singular toric dP4 "de type 4A_1" defined by x_0x_1 = x_2x_3 = x_4² in
  P^4.
* **The nodal family (p. 37).** C_{a,b} : axy(t−s)² = bst(y−x)², (a,b) ∈ N*², pgcd(a,b) = 1, of class −K_{Y_4} with
  deg_{ω^{-1}} C_{a,b} = 4. **Lemme 5.2:** for a ≠ b, C_{a,b} is a "courbe cubique géométriquement intègre et nodale
  en Q", with tangent slopes ±√(b/a). **Lemme 5.3:** "Pour tout point rationnel P = (w_0,z_0) ∈ R où z_0, w_0 ∈ Q, il
  existe une unique courbe C_{a,b} passant par ρ^{-1}(P)." (Here R = {z > w > 1} is one of 8 symmetric regions.)
* **(5.10)-(5.11) (p. 38).** Case I (a, b both squares): α(Q,C_{a,b}) = deg C_{a,b} = 4. Case II ((a,b) ∉ □²):
  approximating Q along a branch "est équivalente à l'approximation du point irrationnel quadratique réel Q_i", so
  α(Q,C_{a,b}) = deg/2 = 2.
* **Thm 5.4 (p. 38).** "On a α_ess(Q) = α(Q,Y_4) = 2. Par conséquent, il n'y a pas de sous-variétés localement
  accumulatrices (Définition 2.3)." Proof: the Case-II curves are dense (even in the analytic topology), and Prop. 5.1
  gives the lower bound.
* **Thm 5.5 (p. 41).** "Pour tout ε > 0, on a que, pour tout B ≫_ε 1, δ_{V,Q,B,2}(χ(ε)) =
  Σ_{(a,b) ∉ □², a < b ≤ ε²} #E(a,b,ε,B,r) = O_ε(1). Par conséquent, la sous-variété ∪_{i=1}^4 Z_i est celle
  localement faiblement accumulatrice (définition 2.9). Il n'existe pas de mesure limite pour le zoom critique." The
  proof uses D_1D_2d_3 | au² − bv² and b ≤ ε², so only finitely many nodal curves meet a fixed critical zoom window;
  each contributes O(1) by Prop. 4.3. On the Z_i, by Pagelot's theorem, δ_{Z_i,Q,B,2}(χ(ε)) ≫≪_ε B^{1/2}.
* **Thm 5.6 (p. 42).** "Pour 2 < r < 144/55, 0 < η < 1/35, on a δ_{V,Q,B,r}(χ(W(ε_1,ε_2,τ_1,τ_2))) ≥
  B^{(1+η)(1/2−1/r)}(log B)^3 (C_2 ∫ χ(W) dw'dz'/z' + O_{τ_i,ε_i}(1/log B)), où C_2 = C_2(r) = (6/π²)(η(1−2/r))³
  ∏_p (1−1/p)³(1 + 3/p − 1/p² − 18/(p(p+2)))." (W = {ε_2 < z' ≤ ε_1, τ_2 ≤ z'/w' ≤ τ_1} (5.27).) This is a **lower
  bound only**. §5.3.3.3 (p. 49) explains why no asymptotic was obtained: the large-(a,b) range, where "le zoom sur la
  courbe nodale C_{a,b} compte au plus un point".
* **Prop. 5.7 (p. 42).** For 2 < r < 144/55 and (a,b) with pgcd 1, a < b, (a,b) ∉ □², τ_2² < ba^{-1} < τ_1² and
  b^{17/2} ≤ 16τ_2/(τ_1−1) B^{1/2−1/r}: #E(a,b,ε,B,r) = (3/(2π²)) Φ(b)Φ(a)Ψ(b−a)/(b a^{1/2}) · εB^{1/2−1/r} +
  O_{τ_i,ε,δ}(b^{23/8+δ} B^{(3/4)(1/2−1/r)} log B). This is "équidistribution locale des points rationnels sur chaque
  C_{a,b}".
* **§6 (pp. 49-50), thin-set interpretation.** h = g∘f : Y_4 → P^1 × P^1 "est un revêtement lisse de degré générique
  2", with "ω^{-1}_{Y_4} = f^*(ω^{-1}_V) = h^*(O_{P^1×P^1}(1,1))". The image of V(Q) is the thin set
  {[u:v] × [s:t] : uvst ∈ □} (6.2). "Ceci pourrait réexpliquer le fait que l'approximation rationnelle des nombres
  rationnels est équivalente à l'approximation rationnelle des nombres quadratiques."
* **Prop. 7.1 (p. 50).** For V_n = Y_4 × P^n and Q ∈ V_n(Q) in the open orbit: "α_ess,ω^{-1}_{V_n}(Q) =
  α_ess(π_1(Q)) + α_ess(π_2(Q)) = n + 3, qui peut être atteinte sur une famille de courbes nodales couvrant un ouvert
  dense de V_n." The lower bound is the product of (5.8) with the P^n Liouville bound d_2^{n+1}H ≥ 1. The upper bound
  comes from the nodal curves D_{a,b,l} = Ψ_{a,b,l}(P^1), of degree 2n + 6, which have α = (2n+6)/2 = n + 3 and whose
  union over (a,b) ∉ □² and lines l is Zariski dense (p. 51).

### 3.2 Part I: quadratic points on P^1 and Pell-Fermat equations

* **Thm 1.3 (pp. 3-4; Thms 4.1, 4.8) [quote].** "On munit le fibré O(1) sur P^1 d'une hauteur de Weil (cf. (4.1)).
  Soit θ un nombre algébrique irrationnel et réel. On identifie θ à un point rationnel [θ:1] de P^1. Alors
  α(θ,P^1) = 1/2.
  — Si r = 1/2 et θ est un nombre quadratique, alors pour toute fonction f ∈ C^b_θ(P^1), δ_{P^1,θ,B,1/2}(f) =
  O_{θ,f}(1) et il existe certaines fonctions g […] telles que liminf_B δ_{P^1,θ,B,1/2}(g) <
  limsup_B δ_{P^1,θ,B,1/2}(g). En d'autres termes, il n'existe pas de mesure limite pour le zoom critique.
  — Si r > 1/2. Alors pour toute f ∈ C^b_θ(P^1), δ_{P^1,θ,B,r}(f) = B^{2−1/r}(6/(π² sup(1,|θ|)²) ∫ f(x)dx +
  o_{θ,f}(1))."
  (The constant 6/π² here differs from the 3/π² in Thm 4.8 and (4.47); this is an internal inconsistency of a factor 2.)
* **Thm 4.1 (p. 13).** For θ = √(b/a): "Pour C_1, C_2 > 0, il existe ε > η > C_1, ε − η < C_2 tels que
  0 = liminf_B δ_{P^1,θ,B,1/2}(χ(ε,η)) < limsup_B δ_{P^1,θ,B,1/2}(χ(ε,η)) = O_{ε,η}(1)."
* **(4.5) (p. 14), the Pell link.** If |u/v − √(b/a)| ≤ εB^{-2} and max(|u|,|v|) ≤ B, then
  |au² − bv²| ≤ 2ε√(ab) + aε²/B² = 2ε√(ab) + o(1). "On en conclut que |au² − bv²| ne prend qu'un nombre fini de valeurs
  entières et ce nombre est indépendant de B."
* **Prop. 4.3 (p. 14).** For every ε > Ξ(θ) and B ≫_ε 1: δ_{P^1,θ,B,1/2}(χ(ε)) ≤
  6 Σ_{|m| ≤ 2ε√(ab)+1} τ(A'm)(⌊(log ε − log Ξ(θ))/(2 log ε*_{A'B'})⌋ + 1).
  **Prop. 4.4 (p. 14):** there are sequences B_n → ∞ with δ_{P^1,θ,B_n,1/2}(χ(ε,η)) = 0 (gaps). **Prop. 4.5
  (p. 16):** sequences with ≥ 1 point.
* **§3 (pp. 10-13).** Prop. 3.2 and Cor. 3.3 (#I_{|m|} ≤ τ(m_1) ≤ τ(|m|)). Cor. 3.4: U^+_D acts freely on S_{D,m} and
  #(S*_{D,m}/U*_D) ≤ 3#I^+_m ≤ 3#I_{|m|} ≤ 3τ(m). Thm 3.5 (Dirichlet): the solutions of x² − Dy² = 1 form a
  finite-index subgroup of O*_{Q(√d)}. Prop. 3.6: A*_{ab} acts freely on the solutions of ax² − by² = c and preserves
  pgcd(x,y). Cor. 3.7: once S*_{a,b,c} ≠ ∅, there are (x_n,y_n) ∈ S*_{a,b,c} with x_n, y_n → ∞.
* **Thm 4.8 (p. 18; sub-critical zoom, any real α with m(α) = 2).** "Soient α un nombre réel et r > 1/2. Alors pour
  tout ε_1 > ε_2 ≥ 0 et τ > 0, — si 1/2 < r ≤ 1, supposons que m(α) = 2, on a δ_{P^1,α,B,r}(χ(ε_1,ε_2)) =
  B^{2−1/r}(3/(π² sup(1,α²))) ∫ χ(ε_1,ε_2)dx + O_{α,ε_i,τ}(B^{1−1/(2r)+τ}) (4.16)"; with the effective Liouville
  hypothesis (4.17) the error becomes O_{ε_i}(Δ(α)Ξ(α)^{-1/2}B^{1−1/(2r)} log B log(…)) (4.19); "si r > 1, […]
  O_{|α|,ε_i}(B log B)". Thm 4.17 (p. 27) is a lattice version (1/2 < r < 7/10) with the main term
  Θ(Λ)εK²/2 · B^{2−1/r}.
* **Appendix A, Thm A.1 (Pagelot, p. 52).** For Q ∈ P^1(Q): α(Q,P^1) = α_ess(Q) = 1. At r = 1,
  δ_{P^1,Q,B,1}(f) = B ∫ f(x) σ(x)dx/x² + O_f(1). For r > 1, δ = B^{2−1/r}(3/π²)∫f dx + O_f(B(log B)^{2/3}(log log
  B)^{4/3}).

### 3.3 Part II: the surface Y_3 (blow-up of P^1 × P^1 at 3 of its 4 torus-invariant points)

* p. 3 **[quote]**: Y_3 "est une surface de del Pezzo généralisée (cf. [8]) (c'est-à-dire la désingularisation minimale
  d'une surface de del Pezzo singulière de degré 5)." p. 13: P_1 = [1:0] × [1:0], P_2 = [0:1] × [1:0],
  P_3 = [1:0] × [0:1], with Q = [1:1] × [1:1]. The minimal-degree curves through Q are l_i : y = x, t = s, yt = xs
  (3.1), and [ω^{-1}_{Y_3}] = O(2,0) + O(0,2) − [E_1] − [E_2] − [E_3].
* **Formule empirique 1.2 (version faible), p. 3 [quote].** "En dehors d'une partie mince M, nous avons que pour toute
  fonction f continue à support compact définie sur (T_Q X)_R, δ_{X∖M,Q,B,α_ess}(f) =
  O_f(B^{1−dim X/α_ess}(log B)^{l(X)−1}) (1.3)."
  **Répartition empirique 1.3 (version forte), p. 3.** "Si α_ess(Q) > dim X, alors en dehors d'une partie mince M,
  nous avons B^{−(1−dim X/α_ess)}(log B)^{−(l(X)−1)} δ_{X∖M,Q,B,α_ess} → δ_{α_ess}, au sens de convergence vague pour
  une certaine mesure δ_{α_ess} sur (T_Q X)_R." These are **proposed empirical formulas, not theorems**. p. 3 says
  [18] and [19] prove 1.3 for the dP6 surface X_3 and 1.2 for Y_4.
* **Thm 1.4 (p. 4; Prop. 4.1, Cor. 4.3) [quote].** "– α(Q,Y_3) = 2. Elle s'obtient sur les trois courbes rationnelles
  lisses l_i (1 ≤ i ≤ 3) de degré minimal passant par Q. Ces courbes sont localement accumulatrices ;
  – α_ess(Q) = 5/2. Elle peut être calculée sur des courbes nodales et des courbes cuspidales passant par Q."
* **Thm 1.5 (pp. 4-5; Thms 7.6, 7.15, 7.16) [quote].** "Soient Z = ∪_{i=1}^3 l_i, U = Y_3 ∖ Z.
  1. Pour 2 = α(Q,Y_3) ≤ r < 5/2, nous avons B^{−(1−2/r)} δ_{Y_3,Q,B,r} → δ_r, où δ_r est une mesure à support dans
     Z ;
  2. Pour r = α_ess(Q) = 5/2, il existe une partie mince, qui est la réunion de Z du type I et M du type II, telle
     qu'en notant V = Y_3(Q) ∖ (Z(Q) ∪ M), nous ayons B^{−1/5} δ_{V,Q,B,5/2} → δ_{5/2}, où δ_{5/2} est une mesure
     absolument continue par rapport à la mesure de Lebesgue sur R²."
* **Prop. 4.1 (p. 15):** α(Q,Y_3) = 2, with H d² ≥ |xy|(s−t)²/(pgcd…) ≥ 1 (p. 16).
* **Prop. 4.2 (p. 16) [quote].** "Pour P = [x:y] × [s:t] n'appartenant pas aux trois courbes l_i, 1 ≤ i ≤ 3, i.e.
  x ≠ y, s ≠ t, xs ≠ yt, (4.1) tel que d(P) ≤ C avec 0 < C < 1, il existe alors D = D(C) > 0 tel que
  H_{ω^{-1}_{Y_3}}(P) d(P)^{5/2} ≥ D." The proof (pp. 16-17) bounds d^{5/2} ≥ |x/y − 1||s/t − 1||xs/yt − 1|^{1/2}
  (1+C)^{-1/2} and then
  H d^{5/2} ≥ ( |y|/(pgcd(y,t)pgcd(y,s)) · |t|/(pgcd(y,t)pgcd(x,t)) · |yt−xs|/(pgcd(x,t)pgcd(y,s)) )^{1/2}
  · |x−y||t−s|(1−C)²(1+C)^{-1/2} ≥ (1−C)²(1+C)^{-1/2}. Every factor in the bracket is a nonzero integer.
* **Cor. 4.3 (p. 17).** "Soit U l'ouvert Y_3 ∖ (Z = ∪_{i=1}^3 l_i). Nous avons α(U,Q) = α_ess(Q) = 5/2. Les courbes
  (3.1) sont les variétés localement accumulatrices." Upper bound: "Les courbes nodales C_{a,b}, ab ≠ □ vérifient
  α(C_{a,b},Q) = 5/2 d'après le Théorème 2.3 et leur réunion couvre l'ouvert dense U."
* **Prop. 4.4 (p. 17).** "Pour ε > 0 fixé, on a que, pour tout B suffisamment grand et pour tout δ > 0,
  #{P ∈ U(Q) : H_{ω^{-1}}(P) ≤ B, d(P) ≤ εB^{-2/5}} ≪_{δ,ε} B^{1/5+δ}."
* **The families (§3.4, p. 15).** Sections c(t−s)²xy + d(y−x)²st − eyt(y−x)(t−s) = 0 of ω^{-1}, (c,d,e) ∈ Z³_prem.
  The tangent slopes at Q coincide iff e² = 4cd, i.e. c = r_1², d = r_2², e = 2r_1r_2 (3.12). This gives the
  **cuspidal** family R_{r_1,r_2} (3.13), with tangent slope r_2/r_1 at Q, and the **nodal** family
  C_{a,b} : axy(s−t)² = bst(x−y)² (3.15), with slopes ±√(b/a), "irrationnelles si et seulement si ab ≠ □".
* **The thin set (§6.2, p. 22).** M = {(w,z) ∈ Q² : −w − z − wz = □ ≠ 0} ⊂ R_2, "dense dans S_2 pour la topologie
  réelle". It is the image of V : x_3² = (x_2² − x_1)(1 − x_1) under a generically degree-2 map π, so it is thin of type
  II. **Prop. 6.2:** "Un point P = (w,z) ∈ Q² est dans M si et seulement s'il existe une courbe cuspidale R_{r_1,r_2}
  (3.14) passant par P."
* **Prop. 7.5 (p. 26).** Under the nodal parametrisation, (ρ∘Ψ)(A_{C_3,D}) ∩ M ≠ ∅ iff −C_3D = □, "En particulier dans
  ce cas le polynomial C_3X² + D est réductible dans Z[X]".
* **Thm 7.6 (p. 27), region S_1 [quote].** "On a δ_{ρ^{-1}(S_1)∩U,Q,B,5/2}(χ_{ε,θ}) = B^{1/5}( ∫ χ_{ε,θ}(w,z)
  E(wz√(w+z))/(wz√(w+z)) dw dz + O_{ε,θ}((log log B)^{5/6}/(log B)^{(2−√2)/6}) ), où E : ]0,∞[ → [0,∞[ est une
  fonction croissante en escalier." Remark: "La mesure obtenue fait apparaître les courbes l_i […] et la somme de leur
  puissances est égale exactement à la constante essentielle 5/2." Thms 7.15 and 7.16 (p. 40) are the analogues on
  S_2, S_3, S_4 (outside M).
* **Lemme 7.17 (p. 40).** "Pour toute région R à support compact dans S_3 ∪ S_4, nous avons
  δ_{M,Q,B,5/2}(1_R) ≪_R B^{1/5} log B." p. 41: a lower bound ≫_ε B^{1/5} log B is expected (Prop. A.7 and §7.3), so the
  thin set carries more points at the critical zoom ("non-équirépartition des points dans M").
* **Thm 2.4 ([19] Thm 1.3) and Thm 2.5 (Pagelot), p. 12-13.** On the nodal cubic C : y² = x³ + ax² (a not a square):
  α(Q_0,C) = α_ess(Q_0,C) = 3/2 and δ_{C,Q_0,3/2,B}(χ(ε)) = O_ε(1). On the cuspidal cubic C' : y² = x³:
  α = α_ess = 3/2 and B^{-1/3} δ_{C',Q_0,3/2,B}(χ(ε)) → δ_{3/2}. "La manière dont on calcule α peut entraîner de
  différents ordres de grandeur dans le zoom."
* **Appendix A (pp. 41-46).** Erdős's mean value for ρ_F (Prop. A.1). Hooley's equidistribution of roots of an
  irreducible F mod m (Thm A.2, Prop. A.3, with discrepancy O_F((log log N)^{(d²+1)/3}/(log N)^{2δ_d/3})). The split
  case: D_F(s) = ζ(s)²Φ(s), Σ_{n≤X} ρ_F(n) ~ C_F X log X (Prop. A.7), with no uniform distribution mod 1 (Dartyge-Martin).

---

## 4. Mechanism

**Lower bounds (approximation constants), verified for both parts.** Every lower bound on α or α_ess in these papers
is a single explicit Liouville inequality: the height times a power of the distance is at least a product of nonzero
integers. The instances are Part I Prop. 5.1 (d²H ≥ 1 on Y_4), Part II Prop. 4.1 (d²H ≥ 1 on Y_3), Part II Prop. 4.2
(d^{5/2}H ≥ D off the l_i), Part I Prop. 7.1 (a product of two such bounds), and Thm 2.4 / Thm 4.7 on P^1 (Liouville).
No Roth-, Schmidt- or Vojta-type input is used for a rational target. Roth enters only for irrational *algebraic*
targets on P^1, in Thm 4.8 for degree ≥ 3.

[inference, checked by hand] Prop. 4.2 is toric Liouville for one section. F = (x−y)²(t−s)²(xs−yt)·y·t has bidegree
(4,4). Its orders at the blown-up points are: at P_1 (y = t = 0), 2 from yt; at P_2 (x = t = 0), 1 from t and 1 from
xs − yt ≈ x − t; at P_3 (y = s = 0), 1 from y and 1 from xs − yt. So F ∈ H^0(Y_3, 2(−K)), and ord_Q F = 2 + 2 + 1 = 5.
Hence τ_Q(−K_{Y_3}) ≥ 5/2. Conversely a section s of m(−K) with ord_Q s = r gives, along any C_{a,b} not contained in
{s = 0}, α(Q,C_{a,b}) ≥ r/m. The C_{a,b} are Zariski dense, so r/m ≤ 5/2. Therefore **τ_Q(−K_{Y_3}) = α_ess(Q) = 5/2**.
The same argument gives **τ_Q(−K_{Y_4}) = α_ess = 2** (the section xy(t−s)² of −K, order 2, against C_{a,b} at 2) and
τ = α_ess = n + 3 on Y_4 × P^n. In all three examples, α_ess = τ_Q(−K) ≥ dim X.

**Upper bounds (dense families).** α_ess is bounded above by exhibiting a Zariski-dense family of rational curves
through Q, each attaining the constant through McKinnon-Roth Thm 2.16, α(Q,C) = min d/(m_P r_P). There are two kinds.
* **Pell type (nodal):** m_P = 1 and r_P = 2, i.e. two branches with **real-quadratic conjugate tangent slopes**
  ±√(b/a).
* **Rational-cusp type:** m_P = 2 and r_P = 1.

Both give d/2: 4/2 = 2 on Y_4 and 5/2 on Y_3. In Part I, Lemma 5.3 gives *every* rational point of an open region on
exactly one C_{a,b}. This is a complete local parametrisation of Y_4(Q) near Q by a pencil of Pell-type curves.

**Role of Pell-Fermat.** Along a nodal C_{a,b}, approximating Q is approximating the real quadratic √(b/a) on P^1. At
the critical zoom, (4.5) forces |au² − bv²| ≤ 2ε√(ab) + o(1). So the approximants are the solutions of finitely many
generalised Pell-Fermat equations au² − bv² = m. Their primitive solutions form ≤ 3τ(m) free orbits of the unit group
(Cor. 3.4), whose heights grow geometrically. Three consequences follow.
* Each curve contributes O(1) points per critical window (Prop. 4.3).
* There are B-sequences with no points at all (Prop. 4.4).
* There is no limit measure at the critical zoom (Thm 4.1, Thm 5.5). By contrast, the minimal curves Z_i (rational
  slopes) carry ≍ B^{1/2} points at the same zoom. They have the same α but dominate the count, which is the "localement
  faiblement accumulatrice" phenomenon.

In Part II the Pell structure is used for counting rather than for gaps. The pair ((a,b),(u,v)) satisfies
C_3(au² − bv²) = D(b − a) (7.1), and a "descente à la Legendre" eliminates (a,b). The count then becomes one of roots of
the quadratic congruence F(l) ≡ 0 [m] (7.45). Hooley's equidistribution handles irreducible F. Reducible F
(−C_3D = □) is exactly the thin set M of points on cuspidal curves (Prop. 7.5), which must be removed. §6 of Part I
gives a structural reason: h : Y_4 → P^1 × P^1 is a degree-2 cover with ω^{-1} = h^*O(1,1). The quadratic irrationality
lives in the **branches of curves through a rational point of a split variety**, not in any torus.

**Scales.** The zoom of factor r counts x with H(x) ≤ B and d(Q,x) ≤ εB^{-1/r}. The papers prove the following at each
range of r.
* r < α_ess: nothing off Q, by the definition of α (Prop. 2.7).
* α ≤ r < α_ess: concentration on the accumulating curves (II Thm 1.5(1)).
* r = α_ess: asymptotics with **relative** error (II Thm 7.6: relative error (log log B)^{5/6}/(log B)^{(2−√2)/6}),
  or O(1) bounds and oscillation (I Thm 5.5, Thm 4.1).
* r > α_ess: lower bounds (I Thm 5.6, only B^{(1+η)(1/2−1/r)} with η < 1/35, not the expected B^{1−2/r}) or
  asymptotics (P^1).

The critical zoom of Y_4 is at r = 2 = dim, the natural Batyrev-Manin scale. The critical zoom of Y_3 is at
r = 5/2 > dim = 2, a **coarser** scale than the natural one. No result is at a scale finer than both α_ess and the
natural scale. Neither paper treats any example with α_ess < dim X.

---

## 5. Relevance to Question B

1. **Applicability.** Both papers are over Q, for **split** toric surfaces (Y_4, Y_3; also Y_4 × P^n), at the rational
   point Q = [1:1] × [1:1] of the open orbit, with anticanonical Weil heights and the McKinnon-Roth normalisation
   d^γ H. The approximation constant and essential constant are exactly the record's (Part I Def. 2.3 = record lit §1.1).
   Nothing concerns anisotropic tori, dimension 4, or the twisted Losev-Manin space. The record's dictionary
   (A_5(δ < 1/30) ⇔ α^ess_e(−K) > 15/4 off a proper closed set) is consistent with Huang's definitions.

2. **Do they give approximation constants of Zariski-dense sequences?** Yes. α_ess(Q) = 2 (Y_4), 5/2 (Y_3) and n + 3
   (Y_4 × P^n) are exactly the sup over dense constructible sets. Each is attained by a Zariski-dense union of
   infinitely many curves, not by one curve. [inference] In every case α_ess equals the effective threshold
   τ_Q(−K) (§4) and is ≥ dim X: Y_4 sits at τ = α_ess = dim, and Y_3 and Y_4 × P^n at τ = α_ess = dim + 1/2 and dim + 1.
   **These examples confirm the record's statement that in every known example α_ess equals the section/Liouville
   value.** They live in the regime τ ≥ dim, where Liouville already gives the Vojta/McKinnon-Satriano prediction. The
   record's X is in the regime τ_e = 10/3 (balanced 15/4) < dim = 4, which neither paper touches.

3. **Beyond Liouville at a rational point?** No. No theorem in either paper proves α^ess > τ at any point. The lower
   bounds are single-section Liouville inequalities; the upper bounds are explicit curve families. The counting results
   are at r ≥ α_ess with relative errors (or are lower bounds), so they cannot exclude a sparse Zariski-dense set below
   α_ess. That is a statement at r < α_ess, which Huang obtains only from Liouville through Prop. 2.7.
   [inference] Could they be adapted? The engine that could transfer is the "pencil parametrisation + Pell descent"
   (Lemma 5.3, Lemma 7.3, Prop. 7.4). It needs a family of rational curves through e such that every rational point near
   e lies on exactly one member. On X, any such covering family has α(C) ≥ 4 by the record's covering inequality.
   Proving Question B through it would then need a bound d^γ H ≥ c with c uniform in the member. In Huang this
   uniformity comes back from the global Liouville inequality (Prop. 5.1 is used in Thm 5.5, and Lemma 7.1 uses the
   Prop. 4.4 bound). So the adaptation reduces to a uniform Liouville statement, which is capped at 15/4 on X (round 22
   Prop 1.3, round 23 width Thm 2.5, round 24 theory). **No lane opens.**

4. **Pell-Fermat and split Theorem 2 (round23/split.md).** This is the **same mechanism**. Huang's quoted McKinnon-Roth
   formula α(Q,C) = d/(m_P r_P) (I Rem. 2.5 p. 8; II Thm 2.3 pp. 10-11), with r_P = 2 for a branch over a real quadratic
   field, gives the following.
   * Split Theorem 2: d = B·C = 15, two conjugate real-quadratic branches with m_P = 2, so α = 15/(2·2) = 15/4.
   * Round 24 class: d = 30, m_P = 4, so α = 30/(2·4) = 15/4.
   * Huang's nodal curves: d = 4 or 5, m_P = 1, r_P = 2, so α = 2 or 5/2.

   The record's derivation (`round24/search.md` §1.2: Liouville/Dirichlet at the quadratic t_0) is this formula. Huang
   adds two quantitative facts [inference: transferring Prop. 4.3 / Thm 4.1 from P^1 to the normalisation of such a
   curve].
   * Along one Pell-type curve, the approximants at the curve's own critical exponent form finitely many unit orbits
     with geometrically growing heights. So there are O(1) points per critical window, with gaps and no limit measure.
     This matches the record's "dist^{15/4}H bounded along the convergents" and its observation that W-points reach
     γ < 4 only at log H ≈ 106.
   * Such curves therefore contribute negligibly to counts. They matter only for approximation constants.

   The calibration lesson [inference]: in Huang's surfaces the Pell-type curves form a **covering** (Zariski-dense,
   one-parameter) family, and they sit **exactly at the Liouville value**. That is how α_ess = τ is realised there. In
   split LM_5 the Theorem-2 Pell curves also sit exactly at the balanced Liouville value 15/4 but sweep only a
   hypersurface W (a 2-parameter family in a 4-fold). A Huang-type *covering* Pell family at 15/4 would refute split
   Question B. The record's covering bound (α ≥ 4 for covering families) excludes it, and Huang's examples respect that
   bound (Y_4 attains it with equality: α = 2 = dim). So the papers neither decide split Question B nor add a new case.

5. **Pell-Fermat and the twisted search (round24/search.md).** The searched class (two conjugate real-quadratic branches
   of contact 4, B·C = 30) is exactly Huang's "Cas II" (r_P = 2) with m_P = 4. The quadratic irrationality in Huang is
   a property of a curve's branches at a rational point of a **split** variety. It is not a twisted-torus phenomenon,
   which agrees with the round-24 finding that Theorem B (finiteness of Q-members of M_C) is form-independent. The
   record's own open question (round22 lit §9 item 3: "a quadratic-irrational phenomenon inside a split variety might
   be the closest analogue of our anisotropic situation") is answered negatively by the full text. The quadratic field
   in Huang comes from the branch slopes ±√(b/a) of nodal curves (equivalently from the degree-2 cover
   Y_4 → P^1 × P^1 of I §6), not from a non-split torus. The anisotropy of T_5 enters elsewhere (the parity lemma
   E_A·C = E_{A^c}·C).

6. **One structural caution [inference].** Huang's empirical formulas (II Formule 1.2, Répartition 1.3) and II Thm 1.5(2)
   remove a **thin** set (type II, Zariski dense: M is dense in a real region). Question B asks for a **Zariski-closed**
   exceptional set. In Huang's Y_3 the thin set does not lower α_ess, because cuspidal curves have α = 5/2 = α_ess too:
   it changes counts (an extra log B), not approximation constants. So nothing in the papers suggests that approximation
   constants need thin exceptional sets. But the papers also give no evidence either way for the record's 4-fold.

**Verdict: no_help.** The papers neither settle nor narrow Question B. They (a) confirm the record's normalisation and
definitions verbatim, (b) add three more examples where α_ess equals the single-section Liouville value, all with
τ ≥ dim, and (c) identify the record's split Theorem 2 and twisted Pell-type curves as instances of the McKinnon-Roth
r_P = 2 branch phenomenon, which Huang studies at length but only at or above the Liouville threshold.

---

## 6. Corrections to the record

1. `round22/lit.md` §5 H2 (and README tables): "Acta Arith. 189 (2019) 1-94" is **not printed** in either PDF. Part I
   (arXiv v2, 17 Dec 2018) carries no journal data. Part II's bibliography says only "à paraître dans Acta Arith.,
   2019" (p. 47). Keep volume/pages as [S] (unconfirmed).
2. `round22/lit.md` §5 H2 says the papers "study a toric surface Y_4, using Pell-Fermat equations". This is true of
   **Part I only**. **Part II studies Y_3** (P^1 × P^1 blown up at 3 of the 4 torus-fixed points, a generalised dP5),
   with α(Q,Y_3) = 2 and α_ess(Q) = 5/2 (II Thm 1.4, p. 4). In Part II, Pell-Fermat equations enter through
   C_3(au² − bv²) = D(b − a) and the reduction to quadratic congruences (Hooley), not through orbit counting.
3. `round22/lit.md` §5 H2 definitions. The essential constant and "locally accumulating" are confirmed verbatim (I Déf.
   2.3, p. 7; II Déf. 2.2, p. 10, credited to Pagelot). Add the companion notion **"localement faiblement
   accumulatrice"** (I Déf. 2.9, p. 9): α(Q,W) = α_ess(Q), but W dominates the count at the critical zoom. The Z_i on
   Y_4 are of this kind, and **Y_4 has no locally accumulating subvariety** (I Thm 5.4, p. 38).
4. Record lit §5, item 2, and README row "Huang torsor lower bounds … section-type … **assessment**, proofs unread". For
   the Acta Arithmetica papers (H2) this is now **verified**. The lower bounds are explicit single-section Liouville
   inequalities in the Cox coordinates of P^1 × P^1: I Prop. 5.1 (p. 37), II Props. 4.1-4.2 (pp. 15-17), I Prop. 7.1
   (p. 50). [inference, hand-checked] II Prop. 4.2 is Liouville for F = (x−y)²(t−s)²(xs−yt)yt ∈ H^0(2(−K_{Y_3})),
   ord_Q F = 5. The assessment for H3 (ANT, 1809.02001) is unaffected by this reading.
5. Record lit §8 / README J3 ("in every computed example α_ess equals the section/Liouville value"). Add Y_4, Y_3 and
   Y_4 × P^n: [inference] α_ess = τ_Q(−K) = 2, 5/2 and n + 3 respectively, all ≥ dim X. Add the caveat that every Huang
   example is in the regime τ ≥ dim, unlike X (τ_e = 10/3 < 4).
6. Record lit §9 item 3 (Pell-Fermat "might be the closest analogue of our anisotropic situation"): **answered
   negatively**. The quadratic field is that of the branch slopes ±√(b/a) of nodal curves through a rational point of a
   split surface. This is McKinnon-Roth Thm 2.16 with r_P = 2, the same mechanism as split Thm 2 and the round-24 twisted
   class, and it is not tied to anisotropy. Part I §6 adds the interpretation through the degree-2 cover
   Y_4 → P^1 × P^1 with ω^{-1}_{Y_4} = h^*O(1,1).
7. Record wording "counting / local distribution works at the natural scale with relative errors" (lit §0 item 5;
   README row 198). Refine it: Huang's zooms are at radius εB^{-1/r} with r ≥ α_ess. The critical zoom is at the natural
   scale on Y_4 (r = 2 = dim) and **coarser than natural** on Y_3 (r = 5/2 > 2). The results are relative-error
   asymptotics (II Thm 7.6), O(1) bounds and oscillation at r = α_ess (I Thm 4.1, Thm 5.5), or lower bounds only for
   r > α_ess (I Thm 5.6, with η < 1/35). For r < α_ess the only statement is the trivial Prop. 2.7, which *is* the
   Liouville bound. The conclusion (inapplicable to Question B) stands.
8. Minor internal discrepancies in Part I (no consequence for the record): the constant is 6/π² in Thm 1.3 (p. 4) but
   3/π² in Thm 4.8 and (4.47) (pp. 18, 25). The remark after Prop. 5.1 (p. 37) prints zoom factor 1/2 where the context
   requires 2. Thm 1.2(2) states "η ∈ ]0,1[" but proves it only for η < 1/35 (p. 5, Thm 5.6).
9. Huang's examples are **weak** del Pezzo in at least one case. Y_4 has −K nef and big but not ample (p. 36), with
   −K_{Y_4} = f^*(−K_V) for the crepant contraction to the singular toric dP4 of type 4A_1 (p. 49). This parallels the
   record's X → X', but here the contraction is divisorial (four (−2)-curves) rather than small. Lines in the record that
   say Huang's constants are "stated for ample L" (kill_lit) should note that I §2.1 (p. 7) assumes only L big, with
   Northcott on the open orbit.
