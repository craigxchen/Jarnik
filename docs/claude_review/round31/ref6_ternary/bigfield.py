"""F_{3^r} with numpy tables, for direct products over F_{q^2} (q^2 = 729)."""
import numpy as np, itertools

class BF:
    def __init__(self, r):
        q = self.q = 3**r; self.r = r
        for coeffs in itertools.product(range(3), repeat=r):
            if coeffs[0] == 0: continue
            v = [1] + [0]*(r-1); seen = []
            for _ in range(q-1):
                seen.append(sum(c*3**j for j, c in enumerate(v)))
                d = [0] + v[:-1]; top = v[-1]
                v = [(d[j] - top*coeffs[j]) % 3 for j in range(r)]
            if len(set(seen)) == q-1: break
        self.exp = np.array(seen + seen, np.int64)
        self.log = np.zeros(q, np.int64); self.log[self.exp[:q-1]] = np.arange(q-1)
        dig = np.array([[(x//3**j) % 3 for j in range(r)] for x in range(q)])
        pw3 = 3**np.arange(r)
        self.add = ((dig[:, None, :] + dig[None, :, :]) % 3) @ pw3
        self.neg = ((-dig) % 3) @ pw3
    def mul(self, a, b):  # arrays or scalars
        a = np.asarray(a); b = np.asarray(b)
        r = self.exp[(self.log[a] + self.log[b]) % (self.q - 1)]
        return np.where((a == 0) | (b == 0), 0, r)
    def pow(self, a, e):
        if a == 0: return 0 if e else 1
        return int(self.exp[(int(self.log[a])*e) % (self.q-1)])
