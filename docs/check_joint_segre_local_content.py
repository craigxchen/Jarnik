"""Exact core-valuation audit for the joint Segre involution.

This audits the valuation constraints, not the full integer height theorem.
There are six genuine boundary exceptions to the generic output gcd.

Why the finite choices bound real valuations: ceiling preserves every
constraint. If two real valuations differ, their minimum must equal the
integer order of their matching sum/difference; ceiling leaves that minimum
unchanged. If they agree, their common ceiling is still at most that integer
order. Individual unequal-branch constraints are already integer equalities.
All pair orders are at most3, so at most one form can exceed3. Raising that
one value to20 preserves the constraints and increases the objective; an
output monomial omitting that form has valuation at most12. Thus values
0,1,2,3,20 suffice for an upper bound, with no discretization assumption on
the original normalized valuations.
"""

from itertools import product
import json
from pathlib import Path

from check_segre_gradient_arithmetic import MATCHINGS, all_matchings, graph
from check_near_balanced_invariant_congruences import add, scale


def plus(*polynomials):
    result = {}
    for polynomial in polynomials:
        result = add(result, polynomial)
    return result


A, B, C, D, E = [graph(edges) for edges in MATCHINGS]
L = plus(A, B, C, D, E)
forms = [plus(C, scale(B, -1)), plus(C, B), plus(B, C, scale(D, 2)),
         plus(scale(A, 2), B, C), plus(scale(A, 2), B, C, scale(D, 2), scale(E, 2))]
individual = [(C, scale(B, -1)), (C, B), (plus(B, D), plus(C, D)),
              (plus(A, B), plus(A, C)),
              (plus(L, scale(B, -1)), plus(L, scale(C, -1)))]
for form, (left, right) in zip(forms, individual):
    assert form == plus(left, right)

# form_j + sign*form_i = 2*matching; order F,H,U,V,G.
pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2),
         (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
signs = [-1, -1, -1, -1, -1, -1, 1, 1, -1, -1]
matchings = [B, plus(B, D), plus(A, B), plus(A, B, D, E), D,
             A, L, plus(L, scale(E, -1)), plus(A, E), plus(D, E)]
for (i, j), sign, matching in zip(pairs, signs, matchings):
    assert plus(forms[j], scale(forms[i], sign)) == scale(matching, 2)

all_matching_polynomials = [graph(edges) for edges in all_matchings(tuple(range(6)))]
assert len(all_matching_polynomials) == 15
for polynomial in matchings + [p for pair in individual for p in pair]:
    assert any(polynomial == candidate or polynomial == scale(candidate, -1)
               for candidate in all_matching_polynomials)


def order(polynomial, inside):
    return min(sum(monomial[i] for i in inside) for monomial in polynomial)


rows = []
for mask in range(64):
    inside = [i for i in range(6) if mask >> i & 1]
    old = [order(polynomial, inside) for polynomial in [A, B, C, D, E]]
    lower = [order(polynomial, inside) for polynomial in forms]
    pair_orders = [order(polynomial, inside) for polynomial in matchings]
    individual_orders = [(order(left, inside), order(right, inside))
                         for left, right in individual]
    options = []
    for low, (left, right) in zip(lower, individual_orders):
        options.append([min(left, right)] if left != right
                       else list(range(max(low, left), 4)) + [20])

    def common(z):
        f, h, u, v, g = z
        a, b, c, d, e = old
        return min(a + f + u + g, b + u + g + v, c + u + g + v,
                   d + f + g + v, e + 3 * f)

    generic = common(lower)
    best = -1
    witness = None
    feasible = 0
    for z in product(*options):
        if any(min(z[i], z[j]) > q or
               (z[i] != z[j] and min(z[i], z[j]) != q)
               for (i, j), q in zip(pairs, pair_orders)):
            continue
        feasible += 1
        value = common(z)
        if value > best:
            best, witness = value, z
    rows.append({"S": [i + 1 for i in inside], "generic": generic,
                 "maximum": best, "witness_F_H_U_V_G": witness,
                 "feasible_words": feasible})

expected = {(3, 4), (2, 5), (1, 6), (2, 3, 4, 5),
            (1, 3, 4, 6), (1, 2, 5, 6)}
exceptions = {tuple(row["S"]) for row in rows if row["maximum"] > row["generic"]}
assert exceptions == expected
assert all(row["maximum"] == row["generic"] + (tuple(row["S"]) in expected)
           for row in rows)
assert sum(row["generic"] for row in rows) == 144
assert sum(row["maximum"] for row in rows) == 150

path = Path(__file__).with_name("joint_segre_local_content_audit.json")
path.write_text(json.dumps(rows, indent=2) + "\n")
print("All individual decompositions and ten pair identities verified.")
print("All64 cuts checked; generic sum144, valuation upper sum150.")
print("Six extra-content cuts:", sorted(exceptions))
print("Feasible valuation words checked:", sum(row["feasible_words"] for row in rows))
