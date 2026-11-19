import sys,re
D=int(sys.argv[2]); KL=D//2
for line in open(sys.argv[1]):
    m=re.search(r"WIT g=(\d+) sf=(\d) N\(high->low\)=(\d+) head\(low->high\)=(\d+) lows\(base9 idx\)=([\d,]+)",line)
    if not m: continue
    head=[int(c) for c in m.group(4)]
    codes=[]
    for x in m.group(5).strip(',').split(','):
        low=int(x); code=low+sum(h*9**(KL+1+j) for j,h in enumerate(head))
        codes.append(str(code))
    print(f"D={D} M={m.group(1)} base9_marks={','.join(codes)}, sf={m.group(2)}")
