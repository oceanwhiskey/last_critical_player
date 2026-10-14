"""Python port of Staudacher's dynamic-programming power-index kernels.

The translation follows the EPIC_ C++ implementation
``include/index/RawBanzhaf.cpp``, ``RawShapleyShubik.cpp``,
``RawDeeganPackel.cpp`` and ``RawJohnston.cpp``.  The surrounding
normalization wrappers (``Banzhaf``, ``ShapleyShubik``, ``DeeganPackel``
and ``Johnston``) are included as well, so every public function returns
the same mathematical quantity as EPIC.

EPIC computes the indices for players sorted by decreasing weight and then
maps the result back to the original order.  The same convention is used
here.  Players of weight zero contribute nothing and receive index zero.
Arbitrary-precision Python integers and :class:`fractions.Fraction` replace
EPIC's ``lint::LargeNumber`` and ``bigFloat``; this makes the comparison
values exact and deterministic.

The public kernels preserve EPIC's recurrences, with Python arithmetic
optimizations: Deegan-Packel tables use a common integer denominator;
Shapley-Shubik and Johnston cardinality columns are packed into integers
with carry bounds; Johnston deficiency sums use prefix queries over
intervals with the same suffix row.  These operations do not treat
arbitrary-precision integers as constant-cost machine words.

.. _EPIC: https://github.com/jhstaudacher/EPIC
"""
from fractions import Fraction
from math import lcm

def _prepare(q, weights):
    """Validate and sort ``weights`` in non-increasing order.

    Returns ``(q, sorted_weights, order, nz, total)`` where ``nz`` is the
    number of strictly positive weights and ``order`` maps sorted positions
    back to their original player indices.
    """
    weights = [int(w) for w in weights]
    q = int(q)
    n = len(weights)
    order = sorted(range(n), key=lambda idx: (-weights[idx], idx))
    sorted_weights = [weights[idx] for idx in order]
    nz = sum((1 for w in sorted_weights if w > 0))
    total = sum(sorted_weights)
    return (q, sorted_weights, order, nz, total)

def _remap(values, order, nz, zero=0):
    """Map ``nz`` sorted values back to the original player order."""
    n = len(order)
    result = [zero] * n
    for pos in range(nz):
        result[order[pos]] = values[pos]
    return result

def _all_zero(weights):
    zero = Fraction(0)
    return [zero] * len(weights)

def _winning_coalitions_per_weight(q, sorted_weights, nz, total):
    """``numberOfWinningCoalitionsPerWeight`` from EPIC."""
    n_wc = [0] * (total + 1)
    n_wc[total] = 1
    for wi in sorted_weights[:nz]:
        for k in range(q + wi, total + 1):
            n_wc[k - wi] += n_wc[k]
    return n_wc

def _raw_banzhaf_sorted(q, sorted_weights, nz, total):
    n_wc = _winning_coalitions_per_weight(q, sorted_weights, nz, total)
    n_sp = [0] * nz
    by_weight = {}
    for i in range(nz):
        wi = sorted_weights[i]
        if wi in by_weight:
            n_sp[i] = by_weight[wi]
            continue
        helper = [0] * (total + 1)
        m = max(total - wi, q - 1)
        for k in range(total, m, -1):
            helper[k] = n_wc[k]
        for k in range(total - wi, q - 1, -1):
            helper[k] = n_wc[k] - helper[k + wi]
        end = min(q + wi, total + 1)
        for k in range(q, end):
            n_sp[i] += helper[k]
        by_weight[wi] = n_sp[i]
    return n_sp

def raw_banzhaf(q, weights):
    """Return raw Banzhaf swing counts in the original player order."""
    (q, sorted_weights, order, nz, total) = _prepare(q, weights)
    if nz == 0 or q <= 0 or q > total:
        return [0] * len(weights)
    swings = _raw_banzhaf_sorted(q, sorted_weights, nz, total)
    factor = 1 << len(weights) - nz
    return _remap([value * factor for value in swings], order, nz)

def _swings_per_player_and_cardinality(q, sorted_weights, nz, total):
    bits = nz + 1
    mask = (1 << bits) - 1
    n_wc = [0] * (total - q + 1)
    n_wc[total - q] = 1 << bits * nz
    for wi in sorted_weights[:nz]:
        for k in range(wi, len(n_wc)):
            n_wc[k - wi] += n_wc[k] >> bits
    ssi = [[0] * (nz + 1) for _ in range(nz)]
    by_weight = {}
    for i in range(nz):
        wi = sorted_weights[i]
        if wi not in by_weight:
            helper = n_wc.copy()
            for k in range(len(helper) - wi - 1, -1, -1):
                helper[k] -= helper[k + wi] >> bits
            swings = sum(helper[:wi])
            by_weight[wi] = [swings >> bits * p & mask for p in range(nz + 1)]
        ssi[i] = by_weight[wi]
    return ssi

def _raw_shapley_shubik_sorted(q, sorted_weights, nz, total):
    ssi = _swings_per_player_and_cardinality(q, sorted_weights, nz, total)
    factorial = [1] * (nz + 1)
    for i in range(2, nz + 1):
        factorial[i] = factorial[i - 1] * i
    raw = [0] * nz
    for i in range(nz):
        for k in range(1, nz + 1):
            raw[i] += ssi[i][k] * factorial[k - 1] * factorial[nz - k]
    return raw

def shapley_shubik(q, weights):
    """Return the normalized Shapley-Shubik index."""
    (q, sorted_weights, order, nz, total) = _prepare(q, weights)
    if nz == 0 or q <= 0 or q > total:
        return _all_zero(weights)
    raw = _raw_shapley_shubik_sorted(q, sorted_weights, nz, total)
    factorial_n = 1
    for i in range(2, nz + 1):
        factorial_n *= i
    values = [Fraction(value, factorial_n) for value in raw]
    return _remap(values, order, nz, zero=Fraction(0))

def _number_of_minimal_winning_coalitions(q, sorted_weights, nz):
    wc = [0] * (q + 1)
    wc[0] = 1
    w_first = sorted_weights[0]
    if w_first <= q:
        for x in range(q - w_first, -1, -1):
            wc[x + w_first] += wc[x]
    sum_minwc = [0] * nz
    for i in range(1, nz):
        wi = sorted_weights[i]
        if wi < q:
            for x in range(q - wi, q):
                sum_minwc[i] += wc[x]
        for x in range(q - wi + 1, 0, -1):
            wc[x + wi - 1] += wc[x - 1]
    for i in range(nz):
        if sorted_weights[i] >= q:
            sum_minwc[i] += 1
    return sum(sum_minwc)

def _raw_deegan_packel_sorted(q, sorted_weights, nz):
    denominator = lcm(*range(1, nz + 1))
    f = [[0] * (q + 1) for _ in range(nz)]
    f[0][0] = 1
    for i in range(nz - 1):
        wi = sorted_weights[i]
        for k in range(i, -1, -1):
            (source, target) = (f[k], f[k + 1])
            for x in range(q, wi, -1):
                target[x - 1] += source[x - wi - 1]
    b = [[0] * (q + 1) for _ in range(nz)]
    solution = [0] * nz
    w_last = sorted_weights[nz - 1]
    dp_unscaled = 0
    for k in range(nz):
        tmp = denominator // (k + 1)
        for x in range(1, min(w_last, q) + 1):
            b[k][x] = tmp
        for x in range(q + 1):
            dp_unscaled += b[k][q - x] * f[k][x]
    solution[nz - 1] = dp_unscaled
    for i in range(nz - 2, -1, -1):
        wi = sorted_weights[i]
        w_next = sorted_weights[i + 1]
        if wi <= q:
            for k in range(i + 1):
                (current, following) = (b[k], b[k + 1])
                for x in range(q + 1, wi, -1):
                    current[x - 1] = current[x - wi + w_next - 1] + following[x - wi - 1]
        else:
            for k in range(i + 1):
                b[k][q] = denominator
        for k in range(i + 1):
            (source, target) = (f[k], f[k + 1])
            for x in range(wi, q):
                target[x] -= source[x - wi]
        minimum = q if q <= wi else wi + 1
        for k in range(i + 1):
            tmp = denominator // (k + 1)
            for x in range(w_next + 1, minimum):
                b[k][x] = tmp
        dp_unscaled = solution[i + 1] if wi == w_next else 0
        if wi != w_next:
            for k in range(i + 1):
                (reward, counts) = (b[k], f[k])
                for x in range(q + 1):
                    if counts[x]:
                        dp_unscaled += reward[q - x] * counts[x]
        solution[i] = dp_unscaled
    return [Fraction(value, denominator) for value in solution]

def deegan_packel(q, weights):
    """Return the normalized Deegan-Packel index."""
    (q, sorted_weights, order, nz, total) = _prepare(q, weights)
    if nz == 0 or q <= 0 or q > total:
        return _all_zero(weights)
    values = _raw_deegan_packel_sorted(q, sorted_weights, nz)
    total_mwc = _number_of_minimal_winning_coalitions(q, sorted_weights, nz)
    if total_mwc == 0:
        values = [Fraction(0)] * nz
    else:
        values = [value / total_mwc for value in values]
    return _remap(values, order, nz, zero=Fraction(0))

def _find_index_in_rev(sorted_weights, nz, weight):
    idx = nz
    for i in range(nz):
        if sorted_weights[i] <= weight:
            idx = nz - i - 1
            break
    return idx

def _subset_sum_from_smallest(q, sorted_weights, nz):
    h = [[0] * q for _ in range(nz)]
    upper = [0] * nz
    upper[nz - 1] = sorted_weights[nz - 1]
    j = nz - 1
    while True:
        j -= 1
        if j < 0:
            break
        upper[j] = upper[j + 1] + sorted_weights[j]
        if upper[j] >= q:
            for k in range(j + 1):
                upper[k] = q
            break
    h[0][sorted_weights[nz - 1] - 1] = 1
    for i in range(nz - 2, -1, -1):
        idx = nz - i - 2
        wi = sorted_weights[i]
        for x in range(upper[i], wi, -1):
            h[idx + 1][x - 1] = h[idx][x - 1] + h[idx][x - wi - 1]
        h[idx + 1][wi - 1] = h[idx][wi - 1] + 1
        for x in range(wi - 2, -1, -1):
            h[idx + 1][x] = h[idx][x]
    return h

def _check_no_veto(q, sorted_weights, nz):
    for wi in sorted_weights[:nz]:
        if wi >= q:
            raise ValueError('The Johnston index in EPIC requires a game without veto players')

def _johnston_deficiency_sums(q, sorted_weights, nz, h):
    """Sum original deficiency intervals using prefix sums of suffix rows."""
    prefixes = []
    for row in h:
        prefix = [0]
        for value in row:
            prefix.append(prefix[-1] + value)
        prefixes.append(prefix)
    distinct = sorted(set(sorted_weights[:nz]))
    intervals = [(low, high - 1, _find_index_in_rev(sorted_weights, nz, low)) for (low, high) in zip(distinct, distinct[1:])]
    by_weight = {}
    for wi in distinct:
        values = [0] * q
        for (low, high, idx) in intervals:
            high = min(high, wi - 1)
            if high < low:
                continue
            prefix = prefixes[idx]
            for deficiency in range(q - low):
                values[deficiency] += prefix[min(q, deficiency + high + 1)] - prefix[deficiency + low]
        by_weight[wi] = values
    return [by_weight[wi] for wi in sorted_weights[:nz]]

def _forward_counting_packed(rows, wi, upper, bits):
    """One EPIC forward step, with cardinality-minus-one in packed digits."""
    upper += wi
    for x in range(min(len(rows) - 1, upper), wi - 1, -1):
        rows[x - 1] += 1 if x == wi else rows[x - wi - 1] << bits
    return upper

def _raw_johnston_sorted(q, sorted_weights, nz):
    """EPIC's Johnston recurrences with exact packed cardinality counts."""
    _check_no_veto(q, sorted_weights, nz)
    h = _subset_sum_from_smallest(q, sorted_weights, nz)
    surplus_sums = [0] * sorted_weights[0]
    h_sum = 1
    surplus_sums[0] = h_sum
    for l in range(1, len(surplus_sums)):
        idx = _find_index_in_rev(sorted_weights, nz, l)
        if idx < nz:
            h_sum += h[idx][l - 1]
        surplus_sums[l] = h_sum
    deficiency_sums = _johnston_deficiency_sums(q, sorted_weights, nz, h)
    bits = nz + 2 * nz.bit_length() + 1
    mask = (1 << bits) - 1
    forwards_array = [0] * q
    upper = 0
    raw_j = []
    for i in range(nz):
        wi = sorted_weights[i]
        if i > 0:
            upper = _forward_counting_packed(forwards_array, sorted_weights[i - 1], upper, bits)
            if wi == sorted_weights[i - 1]:
                raw_j.append(raw_j[-1])
                continue
        row = deficiency_sums[i][q - 1 - wi] if wi > 1 else 0
        if i > 0:
            contribution = 0
            for x in range(1, q - wi):
                contribution += forwards_array[x - 1] * deficiency_sums[i][q - x - wi - 1]
            for x in range(q - wi, q):
                contribution += forwards_array[x - 1] * surplus_sums[q - x - 1]
            row += contribution << bits
        if i < nz - 1:
            w = sorted_weights[:nz]
            del w[i]
            interm = [0] * q
            interm_upper = 0
            included = 0
            for j in range(i + 1, nz):
                wj = sorted_weights[j]
                while included < j - 1:
                    interm_upper = _forward_counting_packed(interm, w[included], interm_upper, bits)
                    included += 1
                if q - wi - wj <= 0:
                    row += surplus_sums[q - wi - 1] << bits
                else:
                    row += deficiency_sums[j][q - wi - wj - 1] << bits
                    contribution = 0
                    for x in range(1, q - wi - wj):
                        contribution += interm[x - 1] * deficiency_sums[j][q - x - wi - wj - 1]
                    for x in range(q - wi - wj - 1, q - wi - 1):
                        contribution += interm[x] * surplus_sums[q - x - wi - 2]
                    row += contribution << 2 * bits
        raw_j.append(sum((Fraction(row >> s * bits & mask, s + 1) for s in range(nz)), Fraction(0)))
    return raw_j

def johnston(q, weights):
    """Return the normalized Johnston index in the original order."""
    (q, sorted_weights, order, nz, total) = _prepare(q, weights)
    if nz == 0 or q <= 0 or q > total:
        return _all_zero(weights)
    values = _raw_johnston_sorted(q, sorted_weights, nz)
    total_raw = sum(values)
    if total_raw == 0:
        values = [Fraction(0)] * nz
    else:
        values = [value / total_raw for value in values]
    return _remap(values, order, nz, zero=Fraction(0))
