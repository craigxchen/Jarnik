import io,contextlib
from spec import *
from comb import analyse
for name,roots in (('C1',C1),('C2',C2)):
    with contextlib.redirect_stdout(io.StringIO()):
        dirs=analyse(name,roots)
    try:
        run(name+' ALL',roots,dirs,5)
    except AssertionError as e:
        print('assertion failed',e)
