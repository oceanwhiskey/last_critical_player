"""Improved weighted residual evaluations for one-dimensional weighted games.

These are the paper's WeightedResidualShares readings for Shapley--Shubik
on the winning QOBDD and Deegan--Packel on the minimal-winning QOBDD.
The benchmark labels this method bolus_weighted (Bolus_improved in the figures).
"""
from power_indices import weighted_suffix_counts
from fractions import Fraction
from math import factorial, gcd
from sortedcontainers import SortedList
FALSE = 0
TRUE = 1

class _QOBDD:
    """A complete QOBDD with terminals 0 and 1.

    ``T`` and ``E`` map inner-node ids to their then- and else-child.
    ``levels[i]`` contains the inner nodes on level ``i``.  The two terminal
    nodes are considered to be on level ``n + 1``.
    """

    def __init__(self, n, root, T, E, levels, order):
        self.n = n
        self.root = root
        self.T = T
        self.E = E
        self.levels = levels
        self.order = order

    def inner_nodes(self):
        return set(self.T)

def build_winning_qobdd(quota, weights, stats=None):
    """Build the winning QOBDD using Bolus's interval computed tables.

    At each level, the table stores intervals ``(lb, ub]`` of residual
    quotas that induce the same residual game.  ``SortedList`` keeps the
    entries ordered by ``ub`` so the interval containing a quota can be
    located with a binary search, as in Bolus's AVL-tree construction.
    """
    weights = list(weights)
    n = len(weights)
    order = sorted(range(n), key=lambda player: -weights[player])
    sorted_weights = [weights[player] for player in order]
    T = {}
    E = {}
    levels = [[] for _ in range(n + 1)]
    unique = {}
    computed = [SortedList() for _ in range(n + 2)]
    next_id = 2
    suffix = [0] * (n + 2)
    for i in range(n, 0, -1):
        suffix[i] = suffix[i + 1] + sorted_weights[i - 1]
    constant_nodes = {}
    counters = {'interval_lookups': 0, 'interval_insertions': 0}

    def mk(level, else_child, then_child):
        nonlocal next_id
        key = (level, else_child, then_child)
        node = unique.get(key)
        if node is None:
            node = next_id
            next_id += 1
            unique[key] = node
            T[node] = then_child
            E[node] = else_child
            levels[level].append(node)
        return node

    def constant_chain(level, value):
        """Return Bolus's level-specific I_i or O_i dummy chain."""
        if level == n + 1:
            return value
        key = (level, value)
        node = constant_nodes.get(key)
        if node is None:
            child = constant_chain(level + 1, value)
            node = mk(level, child, child)
            constant_nodes[key] = node
        return node

    def build(level, remaining_quota, remaining_weight):
        if level == n + 1:
            if remaining_quota <= 0:
                return (TRUE, float('-inf'), 0)
            return (FALSE, 0, float('inf'))
        if remaining_quota <= 0:
            return (constant_chain(level, TRUE), float('-inf'), 0)
        if remaining_quota > remaining_weight:
            return (constant_chain(level, FALSE), remaining_weight, float('inf'))
        table = computed[level]
        counters['interval_lookups'] += 1
        pos = table.bisect_left((remaining_quota, float('-inf'), -1))
        if pos < len(table):
            (upper, lower, node) = table[pos]
            if lower < remaining_quota <= upper:
                return (node, lower, upper)
        weight = sorted_weights[level - 1]
        (then_child, then_lb, then_ub) = build(level + 1, remaining_quota - weight, remaining_weight - weight)
        (else_child, else_lb, else_ub) = build(level + 1, remaining_quota, remaining_weight - weight)
        lower = max(then_lb + weight, else_lb)
        upper = min(then_ub + weight, else_ub)
        node = mk(level, else_child, then_child)
        table.add((upper, lower, node))
        counters['interval_insertions'] += 1
        return (node, lower, upper)
    (root, _, _) = build(1, quota, suffix[1])
    if stats is not None:
        stats.update(counters, nodes=len(T), interval_widths=[len(table) for table in computed])
    return _QOBDD(n, root, T, E, levels, order)

class _MinWinBuilder:
    """Build the minimal-winning QOBDD from the winning QOBDD."""

    def __init__(self, source):
        self.source = source
        self.n = source.n
        self.T = {}
        self.E = {}
        self.levels = [[] for _ in range(self.n + 1)]
        self.unique = {}
        self.memo = {}
        self.false_memo = {}
        self.next_id = 2
        self.level_of = {FALSE: self.n + 1, TRUE: self.n + 1}
        for i in range(1, self.n + 1):
            for node in source.levels[i]:
                self.level_of[node] = i

    def mk(self, level, else_child, then_child):
        key = (level, else_child, then_child)
        node = self.unique.get(key)
        if node is None:
            node = self.next_id
            self.next_id += 1
            self.unique[key] = node
            self.T[node] = then_child
            self.E[node] = else_child
            self.levels[level].append(node)
        return node

    def false_chain(self, level):
        if level == self.n + 1:
            return FALSE
        node = self.false_memo.get(level)
        if node is None:
            child = self.false_chain(level + 1)
            node = self.mk(level, child, child)
            self.false_memo[level] = node
        return node

    def minwin(self, node):
        if node in (FALSE, TRUE):
            return node
        result = self.memo.get(node)
        if result is not None:
            return result
        level = self.level_of[node]
        else_child = self.minwin(self.source.E[node])
        then_child = self.minwin(self.source.T[node])
        if self.source.T[node] == self.source.E[node]:
            result = self.mk(level, else_child, self.false_chain(level + 1))
        else:
            result = self.mk(level, else_child, then_child)
        self.memo[node] = result
        return result

    def build(self):
        root = self.minwin(self.source.root)
        return _QOBDD(self.n, root, self.T, self.E, self.levels, self.source.order)

class WeightedResidualIndices:
    """The paper's improved weighted residual readings of Bolus's formulas."""

    def __init__(self, quota, weights):
        if quota <= 0:
            raise ValueError('quota must be positive')
        weights = list(weights)
        if not weights or any((not isinstance(w, int) or w < 0 for w in weights)):
            raise ValueError('expected nonnegative integer weights')
        self.quota = quota
        self.weights = weights
        self.winning = build_winning_qobdd(quota, weights)
        self._minimal = None

    @property
    def minimal(self):
        if self._minimal is None:
            self._minimal = _MinWinBuilder(self.winning).build()
        return self._minimal

    @staticmethod
    def _scalar_counts(diagram):
        counts = {FALSE: 0, TRUE: 1}

        def visit(node):
            if node not in counts:
                counts[node] = visit(diagram.T[node]) + visit(diagram.E[node])
            return counts[node]
        visit(diagram.root)
        return counts

    @staticmethod
    def _vector_paths(diagram):
        n = diagram.n
        paths = {node: [0] * (n + 1) for node in diagram.inner_nodes()}
        paths[FALSE] = [0] * (n + 1)
        paths[TRUE] = [0] * (n + 1)
        paths[diagram.root][0] = 1
        for i in range(1, n + 1):
            for node in diagram.levels[i]:
                value = paths[node]
                (e, t) = (paths[diagram.E[node]], paths[diagram.T[node]])
                for k in range(i):
                    if value[k]:
                        e[k] += value[k]
                        t[k + 1] += value[k]
        return paths

    def _remap(self, sorted_values):
        result = [Fraction(0)] * self.winning.n
        for (sorted_pos, player) in enumerate(self.winning.order):
            result[player] = sorted_values[sorted_pos]
        return result

    def shapley_shubik(self, exact=False):
        """Unpacked weighted residual reading with O(n M) coefficient work."""
        diagram = self.winning
        n = diagram.n
        if not n:
            return []
        denominator = factorial(n)
        weights = [denominator // n]
        for m in range(1, n):
            weights.append(weights[-1] * m // (n - m))
        counts = weighted_suffix_counts(weights, diagram.T, diagram.E, diagram.levels)
        paths = self._vector_paths(diagram)
        numerators = []
        for i in range(1, n + 1):
            numerator = 0
            for node in diagram.levels[i]:
                prefix = paths[node]
                then = counts[diagram.T[node]]
                other = counts[diagram.E[node]]
                for k in range(i):
                    if prefix[k]:
                        numerator += prefix[k] * (then[k] - other[k])
            numerators.append(numerator)
        result = self._remap([Fraction(x, denominator) for x in numerators])
        return result if exact else [float(v) for v in result]

    def deegan_packel(self, exact=False):
        """Same MinWin formula with O(n M_min) weighted suffix rows."""
        n = self.winning.n
        if n == 0:
            return []
        minimal = self.minimal
        paths = self._vector_paths(minimal)
        scalar_counts = self._scalar_counts(minimal)
        total_minimal = scalar_counts[minimal.root]
        values = [Fraction(0)] * n
        if total_minimal:
            denominator = 1
            for k in range(1, n + 1):
                denominator = denominator * k // gcd(denominator, k)
            counts = weighted_suffix_counts([denominator // k for k in range(1, n + 1)], minimal.T, minimal.E, minimal.levels)
            for pos in range(n):
                level = pos + 1
                numerator = 0
                for node in minimal.levels[level]:
                    weighted = counts[minimal.T[node]]
                    numerator += sum((paths[node][k] * weighted[k] for k in range(level)))
                values[pos] = Fraction(numerator, denominator * total_minimal)
        result = self._remap(values)
        return result if exact else [float(v) for v in result]
