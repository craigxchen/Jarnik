"""Referee: character margins (v with <b,v> = 0) of the eta=1e-3 line, at weights beyond the T=28 used in coset.md.
Enumerates only characters: v in Z^8 with <b,v>=0 and w(v) <= T (recursion with pruning on the last free coordinate)."""
import math, json
import numpy as np
PS = [13, 17, 29, 37, 41, 53, 61, 73]; LOG = [math.log(p) for p in PS]
d = json.load(open('../coset_checks/line_cosets_out.json')); b = d['b']; u = np.array(d['u'])
split = [5, 13, 17, 29, 37, 41, 53, 61, 73]
def gp(p):
    for bb in range(1, int(math.isqrt(p)) + 1):
        a2 = p - bb * bb; a = math.isqrt(a2)
        if a * a == a2 and a > bb: return a, bb
qs = [split[split.index(p) - 1] for p in PS]
xb = np.array([math.atan2(gp(q)[1], gp(q)[0]) / (math.pi / 2) for q in qs])
# b = (1,-1,2,0,1,-2,1,1); coordinate 3 has b=0 -> free; solve for coordinate 0 (b_0 = 1): v0 = -sum_{j>0} b_j v_j
def chars(T):
    out = []
    idx = [1, 2, 3, 4, 5, 6, 7]
    def rec(t, cur, w):
        if t == len(idx):
            v0 = -sum(b[j] * cur[k] for k, j in enumerate(idx))
            ww = w + abs(v0) * LOG[0]
            if ww <= T and (v0 or any(cur)):
                out.append(([v0] + cur, ww))
            return
        j = idx[t]; k = 0
        while w + k * LOG[j] <= T:
            for s in ((1, -1) if k else (1,)):
                cur.append(s * k); rec(t + 1, cur, w + k * LOG[j]); cur.pop()
            k += 1
    rec(0, [], 0.0)
    return out
res = {}
for eta in (1e-3, 0.0):
    g = xb + eta * u
    for T in (28.0, 34.0, 40.0):
        C = chars(T)
        worst = (1e9, None)
        for v, w in C:
            c = float(np.dot(v, g)); y = math.exp(-w / 2)
            m0 = abs(c - round(c)) - (2 / math.pi) * math.asin(y)
            m1 = abs(c - 0.5 - round(c - 0.5)) - (2 / math.pi) * math.asin(y / math.sqrt(2))
            m = min(m0, m1)
            if m < worst[0]: worst = (m, v, w)
        res[f'eta={eta},T={T}'] = {'n_chars(both signs)': len(C), 'min_margin': worst[0], 'argmin_v': worst[1], 'w': worst[2]}
        print(f'eta={eta} T={T}', res[f'eta={eta},T={T}'], flush=True)
json.dump(res, open('char_check_out.json', 'w'), indent=1)
