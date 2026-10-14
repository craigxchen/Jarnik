import sys
sys.path.insert(0, '.')
from fibration_check import run
tests = [
    (6, [{0, 1, 2, 3}, {2, 3, 4, 5}, {0, 2, 4}], [1, 1, 1]),
    (6, [{0, 1, 2, 3}, {0, 1, 4, 5}, {2, 3, 4, 5}], [1, 1, 1]),
    (6, [{0, 1, 2, 3}, {0, 1, 4, 5}, {2, 3, 4, 5}], [2, 2, 2]),
    (7, [{0, 1, 2, 3}, {0, 1, 4, 5}, {0, 2, 4, 6}], [1, 1, 1]),
]
for d, c, L in tests:
    run(d, c, L); sys.stdout.flush()
