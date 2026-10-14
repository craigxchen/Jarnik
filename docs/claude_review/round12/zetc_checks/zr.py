import re,sys
for f in sys.argv[1:]:
    n=0; gap=0
    for line in open(f):
        m = re.search(r"'rho': (\d+).*'zeta_ub': ([\d.]+)", line)
        if m:
            n+=1
            if float(m.group(2)) < int(m.group(1)) - 1e-6: gap+=1; print(line.strip()[:300])
    print(f, 'printed', n, 'with zeta_ub<rho:', gap)
