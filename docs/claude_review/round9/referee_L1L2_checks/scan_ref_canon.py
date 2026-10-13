def canon(cl):
    best = None
    for c in (1, -1):
        for r in range(4):
            Q = []
            for (x, y) in cl:
                y = c * y
                for _ in range(r):
                    x, y = -y, x
                Q.append((x, y))
            k = tuple(sorted(Q))
            if best is None or k < best:
                best = k
    return best


