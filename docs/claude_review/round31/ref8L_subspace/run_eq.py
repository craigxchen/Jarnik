import sys
sys.argv=[sys.argv[0]]
import check_viii4 as c
c.equal_norm_instance(1000, 3, 1, 2, 3, 10, 'pair5', nmax_E=16)
c.equal_norm_instance(1000, 3, 1, 2, 3, 10, 'pair6', nmax_E=16)
c.equal_norm_instance(1000, 3, 1, 3, 3, 10, 'lin4', nmax_E=16)
c.equal_norm_instance(1000, 3, 3, 2, 9, 10, 'pair5', nmax_E=6)
open('part2_out.txt','w').write("\n".join(c.out)+"\n")
