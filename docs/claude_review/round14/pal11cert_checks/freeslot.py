"""For each feasible profile (profiles2.py) designate a 'free' slot (the type with the most representatives in the
universal pool) and list the types of the remaining 7 slots."""
from collections import Counter
exec(open('profiles2.py').read().split("print('feasible profiles:'")[0])
size = {(2, 6): 66, (3, 9): 1760, (4, 6): 990, (4, 8): 3960, (4, 9): 495, (4, 10): 40260, (6, 2): 110, (6, 4): 1320,
        (6, 6): 11245, (7, 5): 9240, (7, 6): 86867, (8, 4): 30470, (8, 5): 290065, (8, 6): 4038980, (12, 1): 11,
        (12, 2): 880}
need = Counter()
for p in full:
    free = max(set(p), key=lambda ty: size[ty])
    rest = list(p); rest.remove(free)
    print('free', free, 'rest', sorted(Counter(rest).items()))
    for ty in set(rest): need[ty] += 1
print('types needed for the 7 enumerated slots:', sorted(need))
