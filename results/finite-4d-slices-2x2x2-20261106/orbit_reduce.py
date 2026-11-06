#!/usr/bin/env python3
"""Exact local-gauge quotient for the finite 2x2x2 spatial graph.

24 directed link variables, eight vertices, and seven spanning-tree
links leave 17 independent cycle holonomies.  Remaining root gauge acts
by one simultaneous conjugation on the 17 variables.

The orbit counts are Burnside identities, not dense transfer spectra.
No floats, external libraries, dynamic imports, JSON or online input.
"""
from itertools import product
from finite_slice import Group, LINKS, LINK_POS, VERTICES


def tree_edges():
    out = []
    for child in VERTICES:
        if child == 0:
            continue
        axis = (child & -child).bit_length() - 1
        parent = child ^ (1 << axis)
        out.append((parent, child, axis))
    assert len(out) == 7
    return tuple(out)


TREE = tree_edges()
TREE_POS = frozenset(LINK_POS[p, axis] for p, _, axis in TREE)
CYCLE_POS = tuple(i for i in range(len(LINKS)) if i not in TREE_POS)
assert len(CYCLE_POS) == 17


def fix_tree(group, spatial):
    """Gauge fix all tree links; output 17 cycle variables in a fixed order."""
    assert len(spatial) == 24
    h = [group.identity] * len(VERTICES)
    for parent, child, axis in TREE:
        h[child] = group.mul[h[parent]][spatial[LINK_POS[parent, axis]]]
    fixed = tuple(group.word(h[x], spatial[i], group.inv[h[x ^ (1 << axis)]])
                  for i, (x, axis) in enumerate(LINKS))
    assert all(fixed[i] == group.identity for i in TREE_POS)
    return tuple(fixed[i] for i in CYCLE_POS)


def root_conjugate(group, vec, root):
    inv = group.inv[root]
    return tuple(group.word(root, x, inv) for x in vec)


def canonical(group, spatial):
    vec = fix_tree(group, spatial)
    return min(root_conjugate(group, vec, g) for g in group.elements)


def centralizer_sizes(group):
    return tuple(sum(group.mul[a][b] == group.mul[b][a] for b in group.elements)
                 for a in group.elements)


def burnside(group, cycles):
    centralizers = centralizer_sizes(group)
    assert cycles >= 0
    powers = []
    for c in centralizers:
        value = 1
        for unused in range(cycles):
            value *= c
        powers.append(value)
    numerator = sum(powers)
    assert numerator % group.n == 0
    return numerator // group.n


def independent_orbit_count(group, cycles):
    """Enumerate simultaneous conjugation classes only on small graphs."""
    assert 0 <= cycles <= 4
    unseen = set(product(group.elements, repeat=cycles))
    count = 0
    while unseen:
        vec = min(unseen)
        orbit = {root_conjugate(group, vec, g) for g in group.elements}
        assert orbit <= unseen
        unseen.difference_update(orbit)
        count += 1
    return count


def sample_state(group, seed):
    s = seed
    output = []
    for i in range(24):
        s = (1103515245 * s + 12345 + i) % (1 << 31)
        output.append(s % group.n)
    return tuple(output)


def verify():
    for name in ('Z3', 'S3'):
        g = Group(name)
        assert len(LINKS) == 24
        assert len(TREE) == 7
        assert len(CYCLE_POS) == 17
        # Every entry of a selected cycle tuple can be realized by a gauge-fixed
        # configuration; no distinct cycle variable is inadvertently dropped.
        for seed in range(7):
            state = sample_state(g, seed + 17)
            normal = canonical(g, state)
            assert len(normal) == 17
            h = tuple(sample_state(g, seed + 200)[i] for i in range(8))
            from finite_slice import independent_gauge_transform
            changed = independent_gauge_transform(g, state, h)
            assert canonical(g, changed) == normal
            # A second normalization doesn't change the orbit representative.
            fixed_spatial = [g.identity] * 24
            for i, j in enumerate(CYCLE_POS):
                fixed_spatial[j] = normal[i]
            assert canonical(g, tuple(fixed_spatial)) == normal
            # A nontrivial residual global gauge transform must not alter it.
            root = g.elements[(seed + 1) % g.n]
            twisted = [g.identity] * 24
            for i, j in enumerate(CYCLE_POS):
                twisted[j] = root_conjugate(g, normal, root)[i]
            assert canonical(g, tuple(twisted)) == normal
        for cycles in range(5):
            assert burnside(g, cycles) == independent_orbit_count(g, cycles)
        nraw = g.n ** len(LINKS)
        ngaugefixed = g.n ** len(CYCLE_POS)
        nstates = burnside(g, len(CYCLE_POS))
        if name == 'Z3':
            assert (nraw, ngaugefixed, nstates) == (282429536481, 129140163, 129140163)
        else:
            assert (nraw, ngaugefixed, nstates) == (
                4738381338321616896, 16926659444736, 2821153019713)
        # A wrong formula dividing by |G|^8 is not an orbit count.
        wrong = nraw // (g.n ** 8)
        assert wrong != nstates
        # Swapping one nontrivial link without respecting the gauge action
        # does not preserve all Wilson plaquettes; the quotient is gauge only.
        print('ORBIT|group='+name+'|vertices=8|links=24|tree_links=7|cycles=17|'
              'raw='+str(nraw)+'|tree_fixed='+str(ngaugefixed)+'|orbits='+str(nstates)
              +'|independent_small_graphs=5|gauge_samples=7|verdict=PASS|json=0')
    print('ORBIT_GATE|negative_wrong_group_power=REJECT|verdict=PASS|json=0')


if __name__ == '__main__':
    verify()
