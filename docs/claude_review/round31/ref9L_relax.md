# Referee report, block relax: Astra Part IX.4 and the m = 43 computation in IX.6

Scope. This block covers §IX.4 of `astra_lemma_IX/answer.md` ("A squarefree obstruction to improving the larger sieve", lines 2039-2065), the sentence on the m = 43 relaxation in the COMPUTED paragraph of §IX.6 (line 2118), and summary items 8 and 9. I did not referee IX.0-IX.3, IX.5 or the rest of IX.6. My code is in `ref9L_relax/relax.py` and `ref9L_relax/lsbound.py`, in the working directory `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9L_relax/`. It does not import Astra's script. I read her `fake()` function to compare conventions.

## 0. Verdict in brief

* Every mathematical assertion of IX.4 is **correct**: the prime choice, the cut identities, the Sidon signatures, the edge sizes and gcds, the local valuation witnesses, and the non-realisability theorem. I recomputed all of them for m = 11, 19 and 43.
* The m = 43 numbers are **reproduced** with my own primes: 86 prime factors, 903 distinct unordered pair sums (= C(43,2)), and 12,341 failed triangle identities. The count 12,341 equals C(43,3). So **every** unordered triple fails. That count is generic and carries no information. The information is in the theorem that no edge table with these cuts can satisfy the identities.
* **The constant 1/2 is the larger-sieve constant (Proved here).** For q all of whose prime factors exceed b, the pair-valuation form of the larger sieve (Vandermonde refinement) gives b <= (1/2 + o(1)) log q/loglog q at L = C sqrt q. The relaxation attains m = (1/2 + o(1)) log q/loglog q. It is therefore extremal for that argument **to the leading constant**, but only asymptotically. At m = 43 it reaches m/(log q/loglog q) = 0.37, while the sieve allows about 118 marks at the same q.
* **Correction to the interpretation (the main point of this report).** Saying the relaxation "isolates the triangle identity" is accurate but misleading. Two facts matter. (i) Assume |d_ij| < q/3, which holds here. Then all triangle identities over Z are **equivalent** to a *pair-level local* condition: for each separated pair, the signed residue d_ij mod p^f equals the difference of the two prescribed residues. (ii) The relaxation fails this pair condition maximally. At every prime p | q, the squares d_ij^2 mod p over the 462 separated pairs take 426-461 distinct values; a genuine configuration has 1. So the relaxation is a barrier only for arguments that read the edge table through magnitudes, gcds and valuation (collision-partition) patterns. It is not a barrier for arguments that use signed residues of single differences, such as the sum and squared-difference congruences used in IX.3. "A barrier for pair valuations and larger-sieve arguments" is a fair summary. "A barrier for pair inputs" would be false.
* **Novelty.** The record already has Paley/skew-Hadamard profiles as fakes on circles: angle models at W <= (8+o(1)) M log M (`round8/lw.md` Thm 6.1, expert note §1), against the genuine lower bound 3 M log M. Three things are new here. The barrier is in the twin, where no angle model exists. It attains the sharp constant of the pair-valuation bound. And Paley cuts are provably unrealisable in the twin, by an elementary argument valid even at lengths up to q^{1-o(1)}. On circles, the corresponding question (the fully flipped Paley phase lemma) is open in the record.

## 1. Definitions

* Twin: q >= 1 squarefree in this block, with two residues r_{p,0}, r_{p,1} mod p for each p | q. Marks are integers y_1 < ... < y_m in an interval of length L = C sqrt q with y_i ≡ r_{p,b_ip} (mod p). The bit b_ip is the class (signature) of mark i at p.
* An *edge table* is a family of nonzero integers d_ij (i < j). It is *realised* if d_ij = y_j - y_i for marks as above.
* h_ij is the product of the primes p | q at which b_ip = b_jp, and t_ij = d_ij/h_ij.
* Paley data: m is prime with m ≡ 3 (mod 8), and χ is the Legendre symbol mod m. K_iv = χ(v - i), S = K - I, and A = (11^t + S)/2 is the bit matrix (rows are marks, columns v ∈ F_m are semiprimes Q_v). Both primes of Q_v carry column v.

## 2. Results, item by item

### (1) Prime choice. Proved (with one quantitative remark).

Let P = c_0 m (log m)^2 and J = ceil(m log m). The primes ≡ 1 (mod 4) in [P,2P] number ~ P/(2 log P). They form ~ P^2/(8 log^2 P) products of two distinct such primes, all lying in [P^2, 4P^2]. That range splits into ceil(J log 4) <= 2J bins of ratio exp(1/J), so by pigeonhole some bin holds E >> P^2/(J log^2 P) products.

For a fixed first factor p, the admissible partners lie in an interval of length T(e^{1/J} - 1)/p <= 8P/J. So the maximum degree satisfies Δ <= 8P/J + 1 = O(1 + P/J). Greedy matching gives at least E/(2Δ - 1) edges, which is >> P/log^2 P once P/J >= 1. Here P/J ~ c_0 log m → ∞, as she says. Also P/log^2 P = c_0 m (log m)^2/(log m + 2 loglog m + log c_0)^2 ~ c_0 m. So a fixed absolute c_0 gives m disjoint edges, and no prime-tuple hypothesis is used. Her "matching size" is the bound edges/(2 max degree), and it is right.

Then q = ∏ Q_v, s = 2m, and sqrt T ∈ [P, 2P]. Next, W = log q = m log T + O(m/J) = 2m log P + O(m) = 2m log m + 4m loglog m + O(m). This is (2 + o(1)) m log m, so m = (1/2 + O(loglog m/log m)) log q/loglog q.

Remark (Computed): convergence is slow. At m = 43:

| P | c_0 | log q | m/(log q/loglog q) |
|---|---|---|---|
| 6000 | 9.9 | 778.2 | 0.368 |
| 12000 | 19.7 | 837.8 | 0.345 |
| 24000 | 39.5 | 897.1 | 0.326 |

At P = 6000 the best bin had 449 products, maximum degree 6, and a greedy matching of 126 >= 43.

**The larger-sieve constant (Proved, by me).** Let b marks be given, and suppose every prime of q exceeds b. Then:

* prod_{i<j}(y_j - y_i) is divisible by prod_{k<b} k! (Vandermonde).
* At each p | q, at least floor((b-1)^2/4) pairs agree.
* The two divisibilities are coprime.

Hence C(b,2)(W/2 + log C) >= W floor((b-1)^2/4) + sum_{k<b} log k!. This rearranges to b(log b - 3/2 - log C + o(1)) <= W/2, so b <= (1/2 + o(1)) log q/loglog q.

Gallagher's form gives the same constant: use moduli the prime powers up to Q = W/log W outside q, with ν = 2 at p | q. The relaxation has m = (1/2 + o(1)) log q/loglog q, so it is extremal **to the leading constant**. The ratio m/b_max is 1 - O(loglog m/log m). The second-order loss is the factor k ~ c_0 log^2 m in the cofactors: the mean log|t| exceeds the Vandermonde minimum log m - 3/2 by about 2 loglog m.

Computed: the inequality allows b_max = 118, 124 and 132 at the three m = 43 moduli above. The m = 43 example is therefore not extremal numerically. It is only asymptotically extremal.

The task's phrase "Gallagher's larger sieve gives B << log q/loglog q" is therefore right, with constant 1/2 when the prime factors of q exceed B. For q with small primes I did not work out the constant.

Comparison with circles: there the analogous inequality is 3M log M <= W + O(M) (Theorem 1.1, constant 2/3). The record's circle fakes reach only W ≈ 8 M log M or (23/3) M log M. The twin relaxation is the first construction in the record I know of that meets the pair-valuation constant.

### (2) Cuts. Proved and Computed.

* K^t = -K because χ(-1) = -1, which holds for m ≡ 3 (mod 4).
* KK^t = K^tK = mI - J, since the sum over v ≠ i,j of χ(v-i)χ(v-j) is -1.
* Hence S^tS = SS^t = (m+1)I - J and S1 = S^t1 = -1.
* Two rows (and two columns) agree in exactly (m-1)/2 places.
* A has row weight (m-1)/2 and pair overlap (m+1-3)/4 = (m-3)/4. For m ≡ 3 (mod 8) these are odd and even, so AA^t = I over F_2.
* Therefore the rows are independent over F_2, the pair sums A_i + A_j are distinct, and D = 1 + C(m,2). This persists after duplicating each column for its two primes (signatures in F_2^{2m}).
* Each column has (m-1)/2 plus-entries and (m+1)/2 minus-entries. So the number of agreeing pairs per column is ((m-1)/2)^2, which is the minimum possible.

All of this was verified exactly for m = 11, 19, 43 (row agreements 5/9/21, overlaps 2/4/10, AA^t = I mod 2, and 55/171/903 distinct xors).

Side observation (Computed): m ≡ 3 (mod 8) is sufficient but not necessary for the Sidon property. At m = 23 (≡ 7 mod 8) AA^t ≠ I, yet all 253 xors are distinct. At m = 7 only 7 of 21 are distinct.

(IX.2) is correct: h_ij/sqrt(q/T) lies in [e^{-m/(2J)}, e^{(m-1)/(2J)}] = 1 + O(1/log m). Computed at m = 43, P = 6000, the range was [0.99800, 1.00141].

### (3) Edges and local data. Proved and Computed.

Let k = floor(P/(4m)) and w_i = ki + 1_{i=m-1}. Then 0 < |t_ij| <= k(m-1) + 1 < P <= every prime of q, so gcd(d_ij, q) = h_ij. The joint gcd of the t_ij is 1, because t_01 = k and t_{m-2,m-1} = k + 1.

max|d_ij|/sqrt q <= e^{1/(2 log m)}(P/4 + 1)/sqrt T < 1. Computed at m = 43: 0.168, 0.171 and 0.172 for the three values of P.

The mean of log|t_ij| is log k + log m + O(1) = log m + O(loglog m). Computed at P = 6000: 5.885, with log m = 3.761 and log k = 3.526.

The real points x_i = w_i sqrt(q/T) give x_j - x_i = d_ij · sqrt(q/T)/h_ij, which is d_ij(1 + O(1/log m)). This configuration is not exact: an exact real realisation would force the triangle identities.

Precisely what is and is not realised locally (Proved):

* **Realised.** For every prime l and every depth a, the relation {(i,j) : l^a | d_ij} is the collision partition of an actual local point configuration. At l ∤ q this is because d_ij = h_ij(w_j - w_i) with h_ij a unit at l, so the configuration is w. At p | q the configuration b_iv + p w_i has v_p of its differences equal to 1 or 0, matching v_p(d_ij) = v_p(h_ij) at all depths. By CRT, any finite set of places and depths is realised simultaneously. In record vocabulary, this is the valuation pattern and point-residue *partition* data (Prop. 5.1(a) of `round9/L1L2.md` says this is all the cocycle gives at a single prime not dividing the modulus).
* **Realised as well.** The reduced cofactors t_ij are themselves an exact cocycle (t_ij + t_jk = t_ik). The failure comes entirely from the h_ij.
* **Not realised.** The signed residues: d_ij is not ≡ z_j - z_i modulo any l for local points z. This fails at small primes too. Computed at m = 43, P = 6000: the triangle identity fails mod 3, mod 5 and mod 7 on 6387, 9324 and 10393 of 12341 triples, and on at least 9467 triples mod each p | q. It fails even at the pair level for p | q. A genuine configuration has d_ij ≡ ±(r_{p,1} - r_{p,0}) (mod p) for every pair separated at p, with the sign fixed by the classes. In the relaxation, d_ij^2 mod p takes between 426 and 450 distinct values over the 462 separated pairs.

Astra's sentence "exact point-level witnesses for valuations at all powers, not equality of signed edges with those differences in Z_p" is accurate.

### (4) Non-realisability. Proved.

Let a = S_u - S_v. Then a^t1 = 0, and a^tS_w = -1 - (-1) = 0 for w ≠ u, v. So a has zero sum on each class of column w. Marks constant on each class modulo both primes of Q_w therefore give Q_w | a^t y, and hence q/(Q_uQ_v) | a^t y.

a has m+1 nonzero entries, each ±2, and a^t1 = 0. So after translating the interval to [0, L], |a^t y| <= (m+1)L. (Her bound 2mC sqrt q is also valid.) If (m+1)L < q/(Q_uQ_v), then a^t y = 0 for every pair u, v. Then S^ty = c1. det(S^tS) = (m+1)^{m-1} ≠ 0 and S^t1 = -1, so y = -c1 is constant, which is a contradiction.

**Strengthening (Proved).** The argument needs only L < q/((m+1) 16P^4), not L = C sqrt q. So Paley cuts on all m rows are unrealisable in every interval of length up to q^{1-O(1/m)}, for every choice of the residues. Computed at m = 43: log(q/(Q_uQ_v))/log q = 0.953, and q/(Q_uQ_v) > (m+1) sqrt q holds for all u, v.

Triangle consequences (Proved):

* If all triangle identities held, then y_0 = 0, y_j = d_0j would realise the table, with collision partition equal to the two-class cut and length < sqrt q. That is a contradiction.
* **Equivalence.** For any edge table with |d_ij| < q/3 and any residues r: all triangle identities over Z hold if and only if every d_ij ≡ ρ_j - ρ_i (mod q), where ρ_i is the CRT residue of row i. One direction is reduction mod q. For the other, the defect d_ij + d_jk - d_ik is ≡ 0 (mod q) and has absolute value < q.

So for these tables "the triangle identity" and "the signed pair residues at the primes of q" are the same missing input.

### (5) Input classes. Proved, with the statement corrected.

Corrected statement: the relaxation satisfies (a)-(f) and violates (g)-(i).

* (a) Pair magnitudes: |d_ij| < sqrt q.
* (b) Exact gcds: gcd(d_ij, q) = h_ij with a two-class cut at every prime, and minimal agreement counts.
* (c) Integer reduced cofactors t_ij, of mean log size log m + O(loglog m), forming an exact cocycle. The twin analogue of the Gram identity, d = h t with integral t, holds.
* (d) At every prime and depth, the valuation (collision-partition) pattern of a genuine local point configuration, simultaneously at finitely many places.
* (e) Sidon signatures, D = 1 + C(m,2), so the bound B << sqrt D of IX.3 and the no-affine-3-space condition of IX.1 hold.
* (f) An ordered real configuration matching every edge to relative error O(1/log m), with absolute errors up to about sqrt q/log m.
* (g) Signed residues of single edges modulo every p | q. This is a pair-level local input, and it is the one the IX.3 proof uses.
* (h) Signed local consistency (triangle identities) modulo small primes l ∤ q.
* (i) Triangle identities over Z, and therefore also any exact real realisation.

Relation to angle-model fakes. An angle model has exact real points, so its real triangle identities hold. On circles it also satisfies signed point residues at small moduli (axiom (R)). This relaxation has no points at all, neither real nor local signed ones; it is only an edge table. It is therefore a fake for a strictly smaller class of arguments, namely those that read the configuration through magnitudes, gcds, valuation partitions and cut signatures. The larger sieve, in both its Gallagher and Vandermonde forms, is in that class and runs verbatim on the table. So "a barrier for pair valuations and larger-sieve arguments" is fair.

The record has the analogous barrier on circles in the form of skew-Hadamard/Paley angle-model fakes. What is new is listed in §0: the twin setting, where no angle model exists; the sharp constant 1/2; and the unrealisability theorem for Paley cuts. Her remark that the triangle identity is what is missing should be read through the equivalence in (4): what is missing is signed local data, which already appears at the level of single pairs.

### (6) The m = 43 computation. Computed and reproduced.

My own prime choice: all primes ≡ 1 (mod 4) in [P, 2P], products binned by floor(J log(ab/P^2)) with J = 162. I took the bin with the largest greedy matching and the 43 smallest products in it, for P = 6000, 12000 and 24000. This differs from Astra's choice: she used P = 40000 and a window of ratio 1 + 1/(8m) at T = (3P/2)^2.

In all three runs:

* there are 86 distinct primes;
* the duplicated-column signatures in F_2^86 have 903 = C(43,2) distinct unordered pair sums, so D = 904;
* with the convention "unordered triples i < j < l, test d_ij + d_jl = d_il", 12,341 of 12,341 = C(43,3) triples fail.

I reproduce 12,341 exactly, and it means "all triples". All asserted equalities, gcds and size bounds were checked with exact integers. For m = 11 (P = 3000) all 165 triples fail.

## 3. Labels

| Statement | Label |
|---|---|
| Prime choice, matching size, W = (2+o(1)) m log m, m ~ (1/2) log q/loglog q | Proved (asymptotic; slow convergence Computed) |
| The larger sieve gives constant 1/2 (prime factors of q > b); the relaxation is extremal to leading order | Proved (by referee) |
| Cut identities, Sidon signatures, (IX.2) | Proved; Computed for m = 11, 19, 43 |
| Edges: sizes, gcds, mean log cofactor, approximate real configuration | Proved; Computed |
| Local valuation witnesses at all places | Proved; Computed |
| Non-realisability for L = C sqrt q, and in fact for L < q/(16(m+1)P^4) | Proved |
| Triangle identities ⇔ signed pair residues mod q, when |d| < q/3 | Proved (by referee) |
| "Isolates the triangle identity" | Correct only with the equivalence above; the missing input is signed local data, already visible at the pair level |
| m = 43: 86 primes, 903 pair sums, 12,341 failures (= all triples) | Computed, reproduced |
| Barrier for pair-valuation and larger-sieve arguments, not a counterexample | Correct |
