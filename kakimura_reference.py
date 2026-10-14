"""Exact reference implementation of Kakimura--Terai (2026).

The paper expresses Banzhaf and Shapley--Shubik through the truncated
generating functions

    prod_p (1 + x**w_p),       prod_p (1 + y*x**w_p).

Its arithmetic model obtains these products with fast formal-power-series
operations over the reals.  For reproducible benchmark checks we use the
same truncated series and read-off formulas, but multiply them exactly by
divide-and-conquer Kronecker substitution.  Python's arbitrary-precision
integer multiplication supplies the fast convolution primitive and avoids
rounding errors.  This implementation therefore measures a practical exact
realisation of the paper's formulas; it is not a claim that Python exposes
the paper's unit-cost real-RAM bound verbatim.

Only Banzhaf and Shapley--Shubik are covered by the paper.
"""
from fractions import Fraction
from math import factorial

def _encode_univariate(coefficients, block_bytes):
    data = bytearray(len(coefficients) * block_bytes)
    for (index, coefficient) in enumerate(coefficients):
        if coefficient:
            start = index * block_bytes
            data[start:start + block_bytes] = coefficient.to_bytes(block_bytes, 'little')
    return int.from_bytes(data, 'little')

def _multiply_univariate(left, right, left_players, right_players, quota):
    """Multiply nonnegative polynomials modulo x**quota exactly."""
    players = left_players + right_players
    block_bytes = max(1, (players + 1 + 7) // 8)
    product = _encode_univariate(left, block_bytes) * _encode_univariate(right, block_bytes)
    data = product.to_bytes(max(1, (product.bit_length() + 7) // 8), 'little')
    mask = (1 << 8 * block_bytes) - 1
    result = [0] * quota
    for exponent in range(quota):
        start = exponent * block_bytes
        if start >= len(data):
            break
        result[exponent] = int.from_bytes(data[start:start + block_bytes], 'little') & mask
    return (result, players)

def _univariate_product(weights, quota):
    factors = []
    for weight in weights:
        row = [0] * min(quota, weight + 1)
        row[0] = 1
        if weight < quota:
            row[weight] += 1
        factors.append((row, 1))
    if not factors:
        return [1] + [0] * (quota - 1)
    while len(factors) > 1:
        next_level = []
        for index in range(0, len(factors), 2):
            if index + 1 == len(factors):
                next_level.append(factors[index])
                continue
            (left, nl) = factors[index]
            (right, nr) = factors[index + 1]
            next_level.append(_multiply_univariate(left, right, nl, nr, quota))
        factors = next_level
    (row, _) = factors[0]
    return row + [0] * (quota - len(row))

def _encode_bivariate(rows, block_bytes, stride):
    if not rows:
        return 0
    last = (len(rows) - 1) * stride + max((len(row) for row in rows), default=1)
    data = bytearray(last * block_bytes)
    for (cardinality, row) in enumerate(rows):
        base = cardinality * stride
        for (weight, coefficient) in enumerate(row):
            if coefficient:
                start = (base + weight) * block_bytes
                data[start:start + block_bytes] = coefficient.to_bytes(block_bytes, 'little')
    return int.from_bytes(data, 'little')

def _multiply_bivariate(left, right, left_players, right_players, quota, rows_limit):
    """Multiply modulo x**quota and y**rows_limit by exact convolution."""
    players = left_players + right_players
    block_bytes = max(1, (players + 1 + 7) // 8)
    stride = 2 * quota - 1
    product = _encode_bivariate(left, block_bytes, stride) * _encode_bivariate(right, block_bytes, stride)
    data = product.to_bytes(max(1, (product.bit_length() + 7) // 8), 'little')
    row_count = min(rows_limit, len(left) + len(right) - 1)
    result = []
    for cardinality in range(row_count):
        row = [0] * quota
        base = cardinality * stride
        for weight in range(quota):
            start = (base + weight) * block_bytes
            if start >= len(data):
                break
            row[weight] = int.from_bytes(data[start:start + block_bytes], 'little')
        result.append(row)
    return (result, players)

def _bivariate_product(weights, quota, rows_limit):
    factors = []
    for weight in weights:
        rows = [[1]]
        if rows_limit > 1:
            selected = [0] * min(quota, weight + 1)
            if weight < quota:
                selected[weight] = 1
            rows.append(selected)
        factors.append((rows, 1))
    if not factors:
        return [[1] + [0] * (quota - 1)]
    while len(factors) > 1:
        next_level = []
        for index in range(0, len(factors), 2):
            if index + 1 == len(factors):
                next_level.append(factors[index])
                continue
            (left, nl) = factors[index]
            (right, nr) = factors[index + 1]
            next_level.append(_multiply_bivariate(left, right, nl, nr, quota, rows_limit))
        factors = next_level
    (rows, _) = factors[0]
    rows.extend([[0] * quota for _ in range(rows_limit - len(rows))])
    for (index, row) in enumerate(rows):
        if len(row) < quota:
            rows[index] = row + [0] * (quota - len(row))
    return rows

def _prefix_sums(row):
    total = 0
    result = []
    for value in row:
        total += value
        result.append(total)
    return result

def _quotient_coefficient(g, weight, target):
    """[x**target] g/(1+x**weight), as in Kakimura--Terai Lemma 3."""
    if target < 0:
        return 0
    total = 0
    sign = 1
    for exponent in range(target, -1, -weight):
        total += sign * g[exponent]
        sign = -sign
    return total

def _bivariate_quotient_coefficient(g, weight, cardinality, target):
    """[y**cardinality*x**target] g/(1+y*x**weight)."""
    if target < 0:
        return 0
    total = 0
    sign = 1
    max_steps = min(cardinality, target // weight)
    for step in range(max_steps + 1):
        total += sign * g[cardinality - step][target - step * weight]
        sign = -sign
    return total

class KakimuraReference:
    """Exact implementation of the two power-index formulas in the paper."""

    def __init__(self, quota, weights):
        if quota <= 0:
            self.quota = quota
        else:
            self.quota = int(quota)
        self.weights = [int(weight) for weight in weights]
        if any((weight <= 0 for weight in self.weights)):
            raise ValueError('Kakimura--Terai assumes positive integer weights')
        self._hat_banzhaf = None
        self._g_banzhaf = None
        self._hat_shapley = None
        self._g_shapley = None

    def _constant_game(self):
        return self.quota <= 0 or sum(self.weights) < self.quota

    def banzhaf(self):
        """Return the benchmark's normalized Banzhaf vector exactly."""
        n = len(self.weights)
        if not n:
            return []
        if self._constant_game():
            return [Fraction(0)] * n
        q = self.quota
        self._hat_banzhaf = _univariate_product(self.weights, q)
        self._g_banzhaf = _prefix_sums(self._hat_banzhaf)
        denominator = 1 << n - 1
        by_weight = {}
        for weight in set(self.weights):
            upper = _quotient_coefficient(self._g_banzhaf, weight, q - 1)
            lower = _quotient_coefficient(self._g_banzhaf, weight, q - weight - 1)
            by_weight[weight] = Fraction(upper - lower, denominator)
        return [by_weight[weight] for weight in self.weights]

    def shapley_shubik(self):
        """Return the exact Shapley--Shubik vector from equation (SSg2)."""
        n = len(self.weights)
        if not n:
            return []
        if self._constant_game():
            return [Fraction(0)] * n
        q = self.quota
        rows_limit = min(n, q) + 1
        self._hat_shapley = _bivariate_product(self.weights, q, rows_limit)
        self._g_shapley = [_prefix_sums(row) for row in self._hat_shapley]
        denominator = factorial(n)
        factorials = [factorial(k) for k in range(n + 1)]
        by_weight = {}
        for weight in set(self.weights):
            numerator = 0
            for cardinality in range(min(n - 1, rows_limit - 1) + 1):
                upper = _bivariate_quotient_coefficient(self._g_shapley, weight, cardinality, q - 1)
                lower = _bivariate_quotient_coefficient(self._g_shapley, weight, cardinality, q - weight - 1)
                numerator += factorials[cardinality] * factorials[n - cardinality - 1] * (upper - lower)
            by_weight[weight] = Fraction(numerator, denominator)
        return [by_weight[weight] for weight in self.weights]
