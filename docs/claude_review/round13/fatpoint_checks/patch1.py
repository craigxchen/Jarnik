s = open('prodsearch.py').read()
old = s[s.index('    def rec(start):'):s.index('    rec(0)')]
new = '''    supp = [frozenset(x for x in range(n) if pool[i][x]) for i in range(len(pool))]
    inc = [0] * n
    capacity = 2 * n

    def rec(start, used):
        k = len(chosen)
        if k == m:
            leaves[0] += 1
            Z = [list(pool[i][1:]) for i in chosen]
            if in_T(Z):
                found.append(tuple(chosen))
            return
        if used + 2 * (m - k) > capacity:
            return
        for i in range(start, len(pool)):
            S = supp[i]
            if used + len(S) + 2 * (m - k - 1) > capacity:
                continue
            if any(inc[x] >= 2 for x in S):
                continue
            row = NZ[i]
            cnt[row] += 1
            if (cnt <= tau).all():
                for x in S:
                    inc[x] += 1
                chosen.append(i)
                rec(i, used + len(S))
                chosen.pop()
                for x in S:
                    inc[x] -= 1
            cnt[row] -= 1

'''
s = s.replace(old, new).replace('    rec(0)\n', '    rec(0, 0)\n')
open('prodsearch.py','w').write(s)
