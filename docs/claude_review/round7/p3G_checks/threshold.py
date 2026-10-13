"""Evidence: fraction of columns left unpinned by the certified closure, random actual primes of size M^K,
r = b(M-1) columns, X = 2M. Records the tight degree per column."""
import sys, time, json, math
from trend import run
out = open('threshold_out.jsonl', 'a')
Ms = [int(x) for x in sys.argv[1].split(',')]
Ks = [float(x) for x in sys.argv[2].split(',')]
b = int(sys.argv[3]); seeds = int(sys.argv[4]) if len(sys.argv) > 4 else 1
for M in Ms:
    for K in Ks:
        if M ** K > 2e18: continue
        for seed in range(seeds):
            t = time.time()
            d = run(M, b, K, seed=seed)
            d['seed'] = seed; d['time'] = round(time.time() - t, 1)
            d['free_frac'] = d['free'] / d['r']
            out.write(json.dumps(d) + '\n'); out.flush()
            print(f"M={M} b={b} K={K} seed={seed} r={d['r']} pairdeg={d['deg']:.2f} singles/col={d['S']/d['r']:.2f} "
                  f"free={d['free']} ({d['free_frac']:.2f}) half={d['half']} exact={d['exact']} t={d['time']}", flush=True)
