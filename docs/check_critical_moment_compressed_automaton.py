"""Exact bounded-height/permutation automaton for necessary moment patterns.

Reusable interface: Automaton(s_parity, q), initial, terminal,
phases(state), valid(state), transitions(state). A state is (heights, order),
where order lists normalized polynomial values in increasing order.
Transitions yield (plus_mask, next_state, 28-entry distance increment).
Exact distances and path length must be accumulated separately.
"""
from collections import deque
from itertools import combinations, product


PAIRS = tuple(combinations(range(8), 2))


class Automaton:
    def __init__(self, s_parity, q):
        assert s_parity in (0, 1) and 0 <= q <= 4
        self.s_parity, self.q = s_parity, q
        eps = {b: 1 if b < 4+q else -1 for b in range(4, 8)}
        self.delta = (0,)*4+tuple(-eps[b] for b in range(4, 8))
        plus = tuple(b for b in range(4, 8) if eps[b] == 1)
        minus = tuple(b for b in range(4, 8) if eps[b] == -1)
        low, high = (plus, minus) if s_parity else (minus, plus)
        self.kappa_order = low+tuple(range(4))+high
        ranks = {i: k for k, i in enumerate(self.kappa_order)}
        self.kappa_sign = tuple(1 if ranks[i] > ranks[j] else -1 for i, j in PAIRS)
        self.kind = tuple('C' if (i < 4) != (j < 4) else
                          'A' if i < 4 or eps[i] == eps[j] else 'B'
                          for i, j in PAIRS)
        # s may be replaced by its parity in every modulo-four formula.
        self.distance = tuple(2*s_parity+1 if kind == 'C' else 2*s_parity+2
                              for kind in self.kind)
        self.orientation = tuple(ks*(-1)**(s_parity if kind == 'A' else s_parity+1)
                                 for ks, kind in zip(self.kappa_sign, self.kind))
        comparisons = tuple(ks*(-1)**((self.delta[i]+self.delta[j]-d)//2)
                            for (i, j), ks, d in zip(PAIRS, self.kappa_sign, self.distance))
        initial_order = self.order_from_comparisons(comparisons)
        assert initial_order is not None
        self.initial = ((0,)*8, initial_order)
        self.terminal = (self.delta, self.kappa_order)
        assert self.valid(self.initial) and self.valid(self.terminal)

    @staticmethod
    def order_from_comparisons(comparisons):
        """Comparison + means value_i>value_j for i<j."""
        below = [0]*8
        for (i, j), comparison in zip(PAIRS, comparisons):
            below[i if comparison == 1 else j] += 1
        if sorted(below) != list(range(8)):
            return None
        order = tuple(sorted(range(8), key=below.__getitem__))
        ranks = {i: k for k, i in enumerate(order)}
        if any(comparison != (1 if ranks[i] > ranks[j] else -1)
               for (i, j), comparison in zip(PAIRS, comparisons)):
            return None
        return order

    def phases(self, state):
        r, order = state
        ranks = {i: k for k, i in enumerate(order)}
        result = []
        for (i, j), ks, d in zip(PAIRS, self.kappa_sign, self.distance):
            comparison = 1 if ranks[i] > ranks[j] else -1
            e = int(comparison != ks)
            result.append((d+r[i]+r[j]-self.delta[i]-self.delta[j]+2*e) % 4)
        return tuple(result)

    def valid(self, state):
        r, order = state
        if r[0] != 0 or sorted(order) != list(range(8)):
            return False
        if any(abs(r[a]) > 1 for a in range(4)) or max(r[:4])-min(r[:4]) > 1:
            return False
        if any(r[b] not in (0, self.delta[b], 2*self.delta[b]) for b in range(4, 8)):
            return False
        signs = tuple((-1)**(self.delta[i]-r[i]) for i in order)
        if signs != tuple(sorted(signs)):
            return False
        for (i, j), kind, orientation, k in zip(PAIRS, self.kind, self.orientation, self.phases(state)):
            h = orientation*(r[i]-r[j])
            if kind == 'C' and h != (0, 1, 2, 1)[k]:
                return False
            if kind == 'A' and h != (0, 1, 0, -1)[k]:
                return False
            if kind == 'B' and not ((h == 0 and k == 0) or h == (2, 1, 2, 3)[k]):
                return False
        return True

    def next_signs(self, state):
        """Unoriented pattern sign, including successor-B first occurrence."""
        r, _ = state
        result = []
        for (i, j), kind, orientation, k in zip(PAIRS, self.kind, self.orientation, self.phases(state)):
            if kind == 'C':
                value = (1, 1, -1, -1)[k]
            elif kind == 'A':
                value = (1, -1, -1, 1)[k]
            elif r[i] == r[j]:
                assert k == 0
                value = 1
            else:
                value = (-1, 1, 1, -1)[k]
            result.append(orientation*value)
        return tuple(result)

    def raw_transition(self, state, plus_mask):
        """Algebraic update for any mask; reject nontransitive value orders."""
        r, order = state
        bits = tuple(plus_mask >> i & 1 for i in range(8))
        rp = tuple(r[i]+bits[i]-bits[0] for i in range(8))
        ranks = {i: k for k, i in enumerate(order)}
        comparisons = tuple((1 if ranks[i] > ranks[j] else -1)*(-1)**(bits[0]+bits[i]*bits[j])
                            for i, j in PAIRS)
        new_order = self.order_from_comparisons(comparisons)
        if new_order is None:
            return None
        new_state = (rp, new_order)
        return new_state if self.valid(new_state) else None

    def transitions(self, state):
        """Include empty/full constant self-loops as distinct column labels."""
        r, order = state
        zero = sum((-1)**(self.delta[i]-r[i]) == -1 for i in order)
        masks = {0, 255}
        for left in range(zero+1):
            for right in range(zero, 9):
                if left < right:
                    masks.add(sum(1 << i for i in order[left:right]))
        for mask in sorted(masks):
            bits = tuple(mask >> i & 1 for i in range(8))
            if mask in (0, 255):
                new_state = state
            else:
                selected = [k for k, i in enumerate(order) if bits[i]]
                left, right = min(selected), max(selected)+1
                updated = order[:left]+tuple(reversed(order[left:right]))+order[right:]
                if bits[0]:
                    updated = tuple(reversed(updated))
                new_state = (tuple(r[i]+bits[i]-bits[0] for i in range(8)), updated)
                if not self.valid(new_state):
                    continue
            increments = tuple(int(bits[i] != bits[j]) for i, j in PAIRS)
            yield mask, new_state, increments


def verify_bounded_prefixes(depth=5):
    checked_states = checked_edges = 0
    for parity in (0, 1):
        for q in range(5):
            machine = Automaton(parity, q)
            frontier = {machine.initial}
            seen = set()
            for _ in range(depth+1):
                following = set()
                for state in frontier-seen:
                    seen.add(state)
                    assert machine.valid(state)
                    # Row zero is always positive in the normalized value ordering.
                    r, order = state
                    zero = sum((-1)**(machine.delta[i]-r[i]) == -1 for i in order)
                    assert order.index(0) >= zero
                    transitions = {mask: (out, inc) for mask, out, inc in machine.transitions(state)}
                    assert transitions[0][0] == transitions[255][0] == state
                    phases = machine.phases(state)
                    next_signs = machine.next_signs(state)
                    for mask in range(256):
                        direct = machine.raw_transition(state, mask)
                        assert (direct is not None) == (mask in transitions)
                        if direct is None:
                            continue
                        out, inc = transitions[mask]
                        assert direct == out
                        assert machine.phases(out) == tuple((k+x) % 4 for k, x in zip(phases, inc))
                        for (i, j), increment, expected in zip(PAIRS, inc, next_signs):
                            if increment:
                                assert ((mask >> i & 1)-(mask >> j & 1)) == expected
                        following.add(out)
                        checked_edges += 1
                    checked_states += 1
                frontier = following-seen
    print(f'PASS: {checked_states} compressed states through depth {depth}; {checked_edges} transitions.')
    print('PASS: all 256 masks/state, phase updates, direct pair rules, constant loops, and zero-spanning permutation updates.')


if __name__ == '__main__':
    verify_bounded_prefixes()
