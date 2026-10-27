import sys, os, re
root = sys.argv[1]; start = sys.argv[2]
seen = set(); stack = [start]
def path(mod): return os.path.join(root, *mod.split('.')) + '.lean'
while stack:
    m = stack.pop()
    if m in seen: continue
    p = path(m)
    if not os.path.exists(p): continue
    seen.add(m)
    for line in open(p, encoding='utf-8', errors='replace'):
        mm = re.match(r'\s*(?:public\s+)?import\s+(\S+)', line)
        if mm and mm.group(1).startswith('OAI'):
            stack.append(mm.group(1))
tot = 0; bad = []
for m in sorted(seen):
    txt = open(path(m), encoding='utf-8', errors='replace').read()
    tot += txt.count('\n')
    # strip comments crudely
    t2 = re.sub(r'/-.*?-/', '', txt, flags=re.S); t2 = re.sub(r'--[^\n]*', '', t2)
    for pat in [r'\bsorry\b', r'^\s*(?:private\s+)?axiom\s', r'\badmit\b', r'native_decide', r'implemented_by', r'@\[extern']:
        for x in re.finditer(pat, t2, flags=re.M):
            bad.append((m, pat))
print('modules', len(seen), 'lines', tot)
print('flags', bad[:40], len(bad))
