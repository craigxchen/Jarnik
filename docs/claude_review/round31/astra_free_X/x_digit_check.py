"""Exact checks of the digit construction in X.2; run python3 x_digit_check.py."""
from itertools import product
from random import Random


def cell_allocation(base, h, bits, colors, parity, modulus, target):
    """Construct one block with a prescribed sum. All endpoints use integers."""
    g, length = len(colors), len(bits)
    cells = base**length
    assert modulus >= 100*g*cells
    step = 2 if modulus % 2 == 0 else 1
    assert step == 1 or (target - sum(parity)) % 2 == 0
    remaining = target*cells//modulus - g//2
    indices = [0]*g
    unit = 1
    for row in reversed(bits):
        floors = [base//(4*h)*row[c] for c in colors]
        increment = (remaining - sum(floors)) % base
        assert increment <= g
        digits = [v+int(i < increment) for i,v in enumerate(floors)]
        remaining = (remaining - sum(digits))//base
        indices = [v+unit*d for v,d in zip(indices,digits)]
        unit *= base

    logs, limits = [], []
    for z, desired in zip(indices,parity):
        lower = (modulus*z+cells-1)//cells
        upper = modulus*(z+1)//cells
        lower += (desired-lower) % step
        upper -= (upper-desired) % step
        logs.append(lower)
        limits.append(upper)
    residual = (target-sum(logs)) % modulus
    assert residual % step == 0
    for j in range(g):
        used = min(residual,limits[j]-logs[j])
        logs[j] += used
        residual -= used
    assert residual == 0 and (sum(logs)-target)%modulus == 0
    for j,v in enumerate(logs):
        assert indices[j]*modulus <= v*cells <= (indices[j]+1)*modulus
        if step == 2:
            assert v % 2 == parity[j]
        for t,row in enumerate(bits):
            # fractional coordinate tested without performing division
            residue = (base**t*v) % modulus
            lo = base//(4*h)*row[colors[j]]
            assert lo*modulus <= base*residue <= (lo+2)*modulus
    return logs


def main():
    rng=Random(310)
    h=3
    base=64*h*h
    g=2*base
    # The four one-hot rows detect every nonzero vector on four colors.
    bits=[[int(i==j) for j in range(4)] for i in range(4)]
    colors=[i%4 for i in range(g)]
    # Optionally sampled coordinate tuples are supplementary to checking
    # every digit on every assigned coordinate.
    count=words=0
    for prime in (2,3,5,7,17):
        mod=prime
        while mod < 100*g*base**len(bits):
            mod *= prime
        for _ in range(40):
            parity=[rng.randrange(2) for i in range(g)]
            target=rng.randrange(mod)
            if prime==2:
                target=(target+(sum(parity)-target)%2)%mod
            logs=cell_allocation(base,h,bits,colors,parity,mod,target)
            count += 1
            for size in range(1,h+1):
                for signs in product((-1,1),repeat=size):
                    for rep in range(12):
                        indices=[rng.randrange(g) for i in signs]
                        delta=[0]*4
                        for i,v in zip(indices,signs):
                            delta[colors[i]]+=v
                        if any(delta):
                            assert sum(logs[j]*v for j,v in zip(indices,signs))%mod
                            words += 1
    print("digit blocks",count,"slots per block",g,"signed word tests",words,
          "individual digit intervals",count*g*len(bits))


if __name__=="__main__":
    main()
