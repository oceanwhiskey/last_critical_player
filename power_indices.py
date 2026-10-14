"""Power-index algorithms from The Last Critical Player (reviewer edition).

Exact arithmetic; QOBDD readings use decreasing weight order. Experimental
variants, packing, interpolation and compatibility dispatch are omitted.
"""
from fractions import Fraction
from math import gcd
from operator import index as integer_index
O, I = 0, 1

def _armijos_prefix_max_sizes(prefix_sums, size, upper):
    """Largest feasible cardinality at each weight for this sorted prefix.

    The k smallest prefix weights sum to w(N_size) - w(N_{size-k}).
    A monotone cursor finds the largest k with this sum <= c in
    O(upper + size) work, without retaining histories of coefficient rows.
    """
    cutoff = size
    total = prefix_sums[size]
    bounds = []
    for c in range(upper + 1):
        while cutoff and prefix_sums[cutoff - 1] >= total - c:
            cutoff -= 1
        bounds.append(size - cutoff)
    return bounds

def johnston_dp(q, weights):
    """Aggregated Armijos-Toro Johnston DP with combinatorial counting tables.

    Sort positive weights decreasingly and assign each quasi-minimal
    coalition to its last critical player j. Prefix table entries count
    the total weight and cardinality of its other critical players;
    suffix counts record its noncritical continuations. A shared backward
    reward recurrence obtains every player's contribution in one sweep.
    Extremal prefix sums bound weight and cardinality loops. Equal
    weights share their numerator, retaining every reverse transition.

    This adaptation of Procedure 1 takes O(n**2 * q) integer arithmetic
    operations and O(n * q) integer cells, rather than storing individual
    coalitions. Integer operations are not constant-cost bit operations.
    Return exact Fractions in original player order; zero weights and
    games without critical winning coalitions receive zero.
    """
    q = integer_index(q)
    weights = [integer_index(w) for w in weights]
    if any((w < 0 for w in weights)):
        raise ValueError('weights must be nonnegative integers')
    result = [Fraction(0)] * len(weights)
    if q <= 0 or sum(weights) < q:
        return result
    order = sorted((i for (i, w) in enumerate(weights) if w), key=lambda i: -weights[i])
    w = [min(q, weights[i]) for i in order]
    common = 0
    for wi in w:
        common = gcd(common, wi)
    q = (q + common - 1) // common
    w = [wi // common for wi in w]
    n = len(w)
    distinct = sorted(set(w))
    prefix_sums = [0]
    tail_limits = {}
    for wi in w:
        tail_limits[wi] = min(q - 1, prefix_sums[-1])
        prefix_sums.append(prefix_sums[-1] + wi)
    minimum_size = []
    k = 0
    for c in range(q):
        while prefix_sums[k] < c:
            k += 1
        minimum_size.append(k)
    tail = {wi: [0] * (tail_limits[wi] + 1) for wi in distinct}
    counts = [1] + [0] * (q - 1)
    multiplicities = {wi: 0 for wi in distinct}
    for wi in w:
        multiplicities[wi] += 1
    low = 0
    suffix_upper = 0
    for boundary in distinct:
        prefix = [0]
        for c in range(suffix_upper + 1):
            prefix.append(prefix[-1] + counts[c])
        cap = suffix_upper + 1
        for wi in distinct:
            if wi < boundary:
                continue
            row = tail[wi]
            first = max(0, low + q - wi - suffix_upper)
            end = min(tail_limits[wi] + 1, boundary + q - wi)
            start = low + q - wi - first
            stop = boundary + q - wi - first
            for c in range(first, end):
                row[c] += prefix[stop if stop < cap else cap] - prefix[start if start > 0 else 0]
                start -= 1
                stop -= 1
        for _ in range(multiplicities[boundary]):
            suffix_upper = min(q - 1, suffix_upper + boundary)
            for c in range(suffix_upper, boundary - 1, -1):
                counts[c] += counts[c - boundary]
        low = boundary
    prefix_counts = [[0] * (n + 1) for _ in range(q)]
    prefix_counts[0][0] = 1
    for (i, wi) in enumerate(w):
        upper = min(q - 1, prefix_sums[i])
        maximum_size = _armijos_prefix_max_sizes(prefix_sums, i, upper)
        for c in range(min(q - 1, prefix_sums[i + 1]), wi - 1, -1):
            (row, source) = (prefix_counts[c], prefix_counts[c - wi])
            for k in range(minimum_size[c - wi], maximum_size[c - wi] + 1):
                row[k + 1] += source[k]
    denominator = 1
    for k in range(1, n + 1):
        denominator = denominator * k // gcd(denominator, k)
    reciprocal = [denominator // (k + 1) for k in range(n)]
    backward_rewards = [[0] * (n + 1) for _ in range(q)]
    numerators = [0] * n
    for i in range(n - 1, -1, -1):
        wi = w[i]
        upper = min(q - 1, prefix_sums[i])
        maximum_size = _armijos_prefix_max_sizes(prefix_sums, i, upper)
        for c in range(wi, min(q - 1, prefix_sums[i + 1]) + 1):
            (row, source) = (prefix_counts[c], prefix_counts[c - wi])
            for k in range(minimum_size[c - wi], maximum_size[c - wi] + 1):
                row[k + 1] -= source[k]
        tail_row = tail[wi]
        numerator = 0
        reuse_numerator = i + 1 < n and wi == w[i + 1]
        for c in range(upper + 1):
            (row, reverse) = (prefix_counts[c], backward_rewards[c])
            downstream = backward_rewards[c + wi] if c + wi < q else None
            count = tail_row[c]
            sizes = range(minimum_size[c], maximum_size[c] + 1)
            if reuse_numerator:
                for k in sizes:
                    contribution = count * reciprocal[k]
                    if downstream is not None:
                        contribution += downstream[k + 1]
                    reverse[k] += contribution
            else:
                for k in sizes:
                    contribution = count * reciprocal[k]
                    if downstream is not None:
                        contribution += downstream[k + 1]
                    numerator += row[k] * contribution
                    reverse[k] += contribution
        numerators[i] = numerators[i + 1] if reuse_numerator else numerator
    total = sum(numerators)
    if total:
        for (i, player) in enumerate(order):
            result[player] = Fraction(numerators[i], total)
    return result

def build_winning_qobdd(q, sorted_weights):
    """Build a complete winning QOBDD for weights already sorted.

    A node on level ``i`` carries the clamped prefix weight
    ``s = min(q, sum of selected weights before level i)``.  Redundant
    levels are kept as dummy nodes with ``T(v) == E(v)`` so that the
    diagram remains complete.
    """
    n = len(sorted_weights)
    levels = [[] for _ in range(n + 1)]
    suffix = [0] * (n + 2)
    for i in range(n, 0, -1):
        suffix[i] = suffix[i + 1] + sorted_weights[i - 1]
    T = {}
    E = {}
    unique = {}
    next_id = 2
    memo = {}

    def mk(level, e, t):
        nonlocal next_id
        key = (level, e, t)
        node = unique.get(key)
        if node is None:
            node = next_id
            next_id += 1
            unique[key] = node
            T[node] = t
            E[node] = e
            levels[level].append(node)
        return node

    def build(i, s):
        if i == n + 1:
            return I if s >= q else O
        key = (i, s)
        if key in memo:
            return memo[key]
        if s >= q:
            child = build(i + 1, s)
            node = mk(i, child, child)
        elif s + suffix[i] < q:
            child = build(i + 1, s)
            node = mk(i, child, child)
        else:
            e = build(i + 1, s)
            t = build(i + 1, min(q, s + sorted_weights[i - 1]))
            node = mk(i, e, t)
        memo[key] = node
        return node
    root = build(1, 0)
    return (root, T, E, levels)

def weighted_suffix_counts(weights, T, E, levels):
    """Return C_v[k] = sum_l c'_v,l*w_(k+l), through the add/shift recurrence.

    A node on level i needs only k=0,...,i-2, because the current selected
    player is outside the suffix.  This avoids a convolution per node.
    """
    n = len(weights)
    counts = {0: [0] * n, 1: weights}
    for i in range(n, 1, -1):
        for v in levels[i]:
            (e, t) = (counts[E[v]], counts[T[v]])
            counts[v] = [e[k] + t[k + 1] for k in range(i - 1)]
    return counts

class QOBDDIndices:
    """Read off all four indices from a winning QOBDD."""

    def __init__(self, n, root, T, E, levels, order):
        (self.n, self.root) = (n, root)
        (self.T, self.E, self.levels, self.order) = (T, E, levels, order)
        self.cnt = {O: 0, I: 1}
        self.paths = self.P = self.alpha = None
        self.h = {}
        self.last_critical_suffix_counts = None
        self._global_denom = self._shapley_denom = self._shapley_weights = None
        self.shapley_rewards = None

    @classmethod
    def from_weighted_game(cls, q, weights):
        weights = list(weights)
        n = len(weights)
        order = sorted(range(n), key=lambda i: -weights[i])
        (root, T, E, levels) = build_winning_qobdd(q, [weights[i] for i in order])
        return cls(n, root, T, E, levels, order)

    def _all_nodes(self):
        return set(self.T) | set(self.E) | {O, I}

    def _reciprocal_denominator(self):
        """Return ``lcm(1, ..., n)`` for one exact accumulation per player."""
        if self._global_denom is None:
            denom = 1
            for k in range(1, self.n + 1):
                denom = denom * k // gcd(denom, k)
            self._global_denom = denom
        return self._global_denom

    def suffix_counts(self, v=None):
        if v is None:
            v = self.root
        if v in self.cnt:
            return self.cnt[v]
        value = self.suffix_counts(self.T[v]) + self.suffix_counts(self.E[v])
        self.cnt[v] = value
        return value

    def prefix_counts(self):
        paths = {v: 0 for v in self._all_nodes()}
        paths[self.root] = 1
        for i in range(1, self.n + 1):
            for v in self.levels[i]:
                paths[self.T[v]] += paths[v]
                paths[self.E[v]] += paths[v]
        self.paths = paths
        return paths

    def prefix_rows(self):
        n = self.n
        P = {v: [0] * (n + 1) for v in self._all_nodes()}
        P[self.root] = [1] + [0] * n
        for i in range(1, n + 1):
            for v in self.levels[i]:
                pv = P[v]
                (e, t) = (P[self.E[v]], P[self.T[v]])
                for k in range(i):
                    if pv[k]:
                        e[k] += pv[k]
                        t[k + 1] += pv[k]
        self.P = P
        return P

    def minimal_winner_rewards(self):
        zeta = {O: 0, I: 1}
        alpha = {}
        for i in range(self.n, 0, -1):
            for v in self.levels[i]:
                zeta[v] = zeta[self.E[v]]
                alpha[v] = 1 if zeta[self.T[v]] == 1 and zeta[self.E[v]] == 0 else 0
        self.zeta = zeta
        self.alpha = alpha
        return (zeta, alpha)

    def banzhaf(self):
        """Last-critical scalar preparation and reverse swing reading, O(M)."""
        if not self.n:
            return []
        self.prefix_counts()
        self.last_critical_rewards()
        reverse = self.last_critical_suffix_counts
        return [Fraction(sum((self.paths[v] * (self.h[v] + reverse[self.T[v]]) for v in self.levels[i])), 2 ** (self.n - 1)) for i in range(1, self.n + 1)]

    def last_critical_rewards(self):
        """Compute scalar last-critical counts in O(M) work and storage.

        For sorted weighted games, U(v) counts accepting suffixes with a
        critical suffix player. Every swing at v has either the current
        player as its last critical player, or a later critical player:
        d(v) = cnt(T(v)) - cnt(E(v)) = h(v) + U(T(v)). Consequently
        U(v) = d(v) + U(E(v)) and h(v) = d(v) - U(T(v)). This is an
        identity for the local seeds, not for arbitrary H(a, b) pairs.
        """
        if self.last_critical_suffix_counts is not None:
            return self.h
        self.suffix_counts(self.root)
        cnt = self.cnt
        reverse = {O: 0, I: 0}
        h = self.h
        for i in range(self.n, 0, -1):
            for v in self.levels[i]:
                (t, e) = (self.T[v], self.E[v])
                difference = cnt[t] - cnt[e]
                reverse[v] = difference + reverse[e]
                h[v] = difference - reverse[t]
        self.last_critical_suffix_counts = reverse
        return h

    def _shapley_shubik_weights(self):
        if self._shapley_weights is None:
            n = self.n
            weights = [0] * n
            if n:
                weights[0] = 1
                for k in range(2, n):
                    weights[0] *= k
                for m in range(1, n):
                    weights[m] = weights[m - 1] * m // (n - m)
            denominator = weights[0] * n if n else 1
            self._shapley_weights = weights
            self._shapley_denom = denominator
        return self._shapley_weights

    def prepare_shapley_rewards(self):
        """Prepare factorial-scaled last-critical rewards in O(nM) work/slots.

        Requires the scalar weighted game and decreasing player order
        certified by the weighted-game constructors. Let R_v(z) count
        accepting suffixes with a critical selected suffix player. Then

            d_v(z) = c'_T(v)(z) - c'_E(v)(z)
            R_v(z) = R_E(v)(z) + z*d_v(z)
            H(T(v), E(v); z) = d_v(z) - R_T(v)(z).

        For u_m = m!*(n-m-1)!, C_v[k] weights c'_v and V_v[k]
        weights R_v by u_(k+l). Propagating weights directly gives

            V_v[k] = V_E(v)[k] + C_T(v)[k+1] - C_E(v)[k+1]
            local_v[k] = C_T(v)[k] - C_E(v)[k] - V_T(v)[k].

        No H/D pairs, tagged polynomials, coefficient convolutions,
        packing, or floating-point arithmetic are used. C and V on
        level i need k=0,...,i-2; local rewards need k=0,...,i-1.
        The root's V row is empty and never queried. Terminal V rows
        are zero because the empty suffix has no critical player.
        """
        if self.shapley_rewards is not None:
            return self.shapley_rewards
        weights = self._shapley_shubik_weights()
        if not self.n:
            self.shapley_rewards = {}
            return self.shapley_rewards
        counts = weighted_suffix_counts(weights, self.T, self.E, self.levels)
        zero = [0] * self.n
        vulnerable = {O: zero, I: zero}
        local = {}
        for i in range(self.n, 0, -1):
            for v in self.levels[i]:
                (e, t) = (self.E[v], self.T[v])
                (ce, ct) = (counts[e], counts[t])
                (ve, vt) = (vulnerable[e], vulnerable[t])
                local[v] = [ct[k] - ce[k] - vt[k] for k in range(i)]
                vulnerable[v] = [ve[k] + ct[k + 1] - ce[k + 1] for k in range(i - 1)]
        self.shapley_rewards = local
        return local

    def shared_reward_sums(self, s):
        """Read all marked-path values with one backward polynomial sweep.

        For an unmarked path reaching v with k selected prefix players,
        R[v][k] sums all rewards at v and later levels.  A selected edge
        at level j receives the local reward at j plus the rewards below
        its then-child.  Summing these edge values against P gives every
        player's marked reading in one pass.
        """
        if self.P is None:
            self.prefix_rows()
        if s == 'alpha':
            if self.alpha is None:
                self.minimal_winner_rewards()
            statistic = self.alpha
        elif s == 'h':
            self.last_critical_rewards()
            statistic = self.h
        else:
            raise ValueError("expected 'alpha' or 'h'")
        n = self.n
        denominator = self._reciprocal_denominator()
        reciprocal_num = [denominator // k for k in range(1, n + 1)]
        local = {v: [statistic[v] * reciprocal_num[k] for k in range(i)] for i in range(1, n + 1) for v in self.levels[i]}
        return self._reverse_weighted_rewards(local, denominator)[1:]

    def _reverse_weighted_rewards(self, local, denominator):
        """Accumulate all selected-edge rewards from local node rows."""
        zero = [0] * (self.n + 1)
        reverse = {O: zero, I: zero}
        for i in range(self.n, 0, -1):
            for v in self.levels[i]:
                else_row = reverse[self.E[v]]
                then_row = reverse[self.T[v]]
                reward_row = local[v]
                reverse[v] = [reward_row[k] + else_row[k] + then_row[k + 1] for k in range(i)]
        values = [Fraction(0)]
        for i in range(1, self.n + 1):
            total = 0
            for v in self.levels[i]:
                reward_row = local[v]
                then_row = reverse[self.T[v]]
                prefix = self.P[v]
                for k in range(i):
                    if prefix[k]:
                        total += prefix[k] * (reward_row[k] + then_row[k + 1])
            values.append(Fraction(total, denominator))
        return values

    def deegan_packel(self):
        """Compute Deegan--Packel with O(nM) work and coefficient storage."""
        self.prefix_rows()
        self.minimal_winner_rewards()
        Delta = self.shared_reward_sums('alpha')
        Delta_sum = sum(Delta)
        return [Fraction(Delta[j - 1], Delta_sum) if Delta_sum else Fraction(0) for j in range(1, self.n + 1)]

    def johnston(self):
        """Reference Johnston reading with per-node coefficient rows."""
        self.prefix_rows()
        self.last_critical_rewards()
        J = self.shared_reward_sums('h')
        J_sum = sum(J)
        return [Fraction(J[j - 1], J_sum) if J_sum else Fraction(0) for j in range(1, self.n + 1)]

    def shapley_shubik(self):
        """Shapley--Shubik via weighted last-critical rewards, O(nM); sorted order with leading zero."""
        if not self.n:
            return [Fraction(0)]
        if self.P is None:
            self.prefix_rows()
        local = self.prepare_shapley_rewards()
        return self._reverse_weighted_rewards(local, self._shapley_denom)
