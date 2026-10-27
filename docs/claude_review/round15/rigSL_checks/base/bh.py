# greedy B_h sets (sums of h elements with repetition distinct as multisets) and Theorem 1(b) constants.
import itertools
def is_Bh(S, h):
    seen = set()
    for c in itertools.combinations_with_replacement(S, h):
        s = sum(c)
        if s in seen: return False
        seen.add(s)
    return True
def greedy_Bh(r, h):
    S = []
    x = 1
    while len(S) < r:
        if is_Bh(S + [x], h): S.append(x)
        x += 1
    return S
def AK(K, B):
    return K * (2 * (K + 1) ** 2 * B + 2 * K + 1)
if __name__ == '__main__':
    for K in (4, 5, 6, 7):
        h = K // 2
        for r in range(4, 9):
            S = greedy_Bh(r, h) if h >= 2 else list(range(1, r + 1))
            A = AK(K, S[-1]); print(K, r, S, "A_K", A, "q bound 4A+5", 4 * A + 5)
