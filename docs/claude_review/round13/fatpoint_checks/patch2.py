s = open('prodsearch.py').read()
old = """    vals = Aarr @ D.T  # |A| x #dirs
    tau = np.array([len(np.unique(vals[:, i])) - 1 for i in range(D.shape[0])])"""
new = """    tau = np.zeros(D.shape[0], dtype=np.int64)
    for s0 in range(0, D.shape[0], 2000):
        vals = Aarr @ D[s0:s0 + 2000].T  # |A| x chunk
        vs = np.sort(vals, axis=0)
        tau[s0:s0 + 2000] = (np.diff(vs, axis=0) != 0).sum(axis=0)"""
assert old in s
s = s.replace(old, new)
open('prodsearch.py','w').write(s)
