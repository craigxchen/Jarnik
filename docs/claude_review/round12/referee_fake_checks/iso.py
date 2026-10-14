import itertools
from rlib import *
n = 8
pal_full = pal(7); core = pal(7, core_only=True); syl = sylvester(3)
A1, _ = box(n, pal_full, 3, 0); A2, _ = box(n, core, 3, 0)
s1 = set(map(tuple, A1)); s2 = set(map(tuple, A2))
print('pal7 L=3 base 0: full-profile box == Paley-core box:', s1 == s2, len(s1))
for base in range(n):
    Ab, _ = box(n, pal_full, 3, base)
    print('  pal7 L=3 base', base, '|A| =', len(Ab))
for base in range(n):
    Ab, _ = box(n, pal_full, 2, base)
    print('  pal7 L=2 base', base, '|A| =', len(Ab))
# permutation mapping Paley-core partition set onto Sylvester partition set
def partset(pats):
    return set(frozenset([frozenset(S), frozenset(set(range(n)) - set(S))]) for S in pats)
PC, PS = partset(core), partset(syl)
found = []
for perm in itertools.permutations(range(n)):
    img = set(frozenset(frozenset(perm[x] for x in S) for S in part) for part in PC)
    if img == PS: found.append(perm)
print('#perms mapping Paley-core cuts onto Sylvester cuts:', len(found), ' with 0->0:', sum(1 for p in found if p[0] == 0),
      ' example:', found[0] if found else None)
# the author's permutation (0,1,2,4,3,6,7,5)
p = (0,1,2,4,3,6,7,5)
img = set(frozenset(frozenset(p[x] for x in S) for S in part) for part in PC)
print("author's permutation maps core onto Sylvester:", img == PS)
