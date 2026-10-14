s = open('prodsearch.py').read()
old = """    for f in found[:20]:
        print('  ', [pool[i] for i in f])"""
new = """    for f in found[:3]:
        print('  ', [pool[i] for i in f])
    # zonotope test for every found certificate: does k0 + sum_i [0,1] z_i' fit in A for some k0?
    if found and not countonly:
        Aset = set(A)
        nz = 0
        for f in found:
            G = np.array([pool[i][1:] for i in f])
            verts = np.array([G[list(I)].sum(axis=0) if I else np.zeros(d, dtype=np.int64)
                              for r in range(len(G) + 1) for I in itertools.combinations(range(len(G)), r)])
            lo = verts.min(axis=0); hi = verts.max(axis=0)
            ok = any(all(tuple(v + np.array(k0)) in Aset for v in verts)
                     for k0 in itertools.product(*[range(-1 - int(l), 2 - int(h)) for l, h in zip(lo, hi)]))
            nz += ok
        print(f'certificates found: {len(found)}; of which genuine zonotopes (translate inside A): {nz}', flush=True)"""
assert old in s
s = s.replace(old, new)
open('prodsearch.py', 'w').write(s)
