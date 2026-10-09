"""Referee check R3: construct.md Remark 2.3 (Walsh cores of order >= 16 have no skew pairing).

Independent backtracking.  Condition for a map p on F_2^t minus {0}:  (x+y).(p(x)+p(y)) = 1 for all distinct x,y.
  (a) t = 2, 3: a bijection p : F_2^t\0 -> F_2^t\0 exists (skew pairing of the Walsh core).
  (b) t = 4: no INJECTIVE p : F_2^4\0 -> F_2^4 (value 0 allowed) exists (the relaxed statement used by
      the restriction argument), hence no skew pairing for t >= 4.
Also re-derives that the condition is the 2x2 minor condition for H(x,a) = (-1)^{x.a} (checked for t=3)."""
import itertools

def dot(x, y): return bin(x & y).count("1") & 1

def search(t, injective_only, allow_zero, find_all=False):
    pts = list(range(1, 2**t))
    vals = list(range(0 if allow_zero else 1, 2**t))
    assign = {}; used = set(); count = [0]; nodes = [0]
    def bt(i):
        nodes[0] += 1
        if i == len(pts):
            count[0] += 1
            return not find_all
        x = pts[i]
        for v in vals:
            if v in used: continue
            if all(dot(x ^ y, v ^ assign[y]) == 1 for y in pts[:i]):
                assign[x] = v; used.add(v)
                if bt(i + 1): return True
                del assign[x]; used.discard(v)
        return False
    found = bt(0)
    return found, count[0], nodes[0], dict(assign) if found else None

for t in (2, 3):
    found, cnt, nodes, a = search(t, True, False)
    print("t=%d: skew pairing (bijection onto nonzero labels) exists: %s  (example %s)" % (t, found, a))
    assert found
# minor-condition equivalence check for t = 3 on the found pairing
t = 3; _, _, _, a = search(3, True, False)
H = lambda x, y: -1 if dot(x, y) else 1
for x in range(1, 8):
    for y in range(1, 8):
        if x != y:
            assert (H(x, a[x]) * H(y, a[y]) != H(x, a[y]) * H(y, a[x])) == (dot(x ^ y, a[x] ^ a[y]) == 1)
print("t=3: minor condition <=> (x+y).(p(x)+p(y)) = 1 verified on the pairing")
found, cnt, nodes, _ = search(4, True, True, find_all=True)
print("t=4: injective p: F_2^4\\0 -> F_2^4 (0 allowed) with the condition: %d solutions (search nodes %d)" % (cnt, nodes))
assert cnt == 0
found, cnt, nodes, _ = search(3, True, True, find_all=True)
print("t=3 (for comparison): injective maps into F_2^3 with 0 allowed: %d solutions" % cnt)
print("R3 DONE")
