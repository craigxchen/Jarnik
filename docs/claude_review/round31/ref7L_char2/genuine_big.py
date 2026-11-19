# larger exact instances of VII.6 (her construction), n=(1,2,3,4), to reach h/D > 1/2 + 1/60 and toward 8/15
import sys
sys.argv = sys.argv[:1]
exec(open("genuine.py").read().split("for ns, L in (((1, 2, 3, 5), 6)")[0])
for ns, L in (((1, 2, 3, 4), 100), ((1, 2, 3, 4), 300), ((1, 2, 3, 4), 700)):
    her(ns, L)
