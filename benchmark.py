"""Remeasure exactly the points and timeout markers in the four paper figures."""
import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
TIMEOUT = 120


def calculation(method, index, q, weights):
    """Load the selected algorithm before starting the measurement clock."""
    if method == 'qobdd':
        from power_indices import QOBDDIndices

        def run():
            diagram = QOBDDIndices.from_weighted_game(q, weights)
            if index == 'shapley_shubik':
                values = diagram.shapley_shubik()[1:]
            else:
                values = getattr(diagram, index)()
            remapped = [Fraction(0)] * len(weights)
            for position, player in enumerate(diagram.order):
                remapped[player] = values[position]
            return remapped
        return run
    if method == 'armijos_toro_dp':
        from power_indices import johnston_dp
        return lambda: johnston_dp(q, weights)
    if method == 'staudacher_dp':
        import staudacher_dp as dp
        if index == 'banzhaf':
            return lambda: [Fraction(x, 2 ** (len(weights) - 1))
                            for x in dp.raw_banzhaf(q, weights)]
        function = getattr(dp, index)
        return lambda: function(q, weights)
    if method == 'kakimura_terai':
        from kakimura_reference import KakimuraReference
        return lambda: getattr(KakimuraReference(q, weights), index)()
    if method == 'bolus_weighted':
        from weighted_residual_indices import WeightedResidualIndices
        return lambda: getattr(WeightedResidualIndices(q, weights), index)(exact=True)
    raise ValueError(method)


def worker(job, case):
    sys.setrecursionlimit(10000)
    import gc
    run = calculation(job['method'], job['index'], case['q'], case['w'])
    samples, result = [], None
    for _ in range(job['repeats']):
        gc.collect()
        start = time.perf_counter()
        values = run()
        elapsed = time.perf_counter() - start
        vector = [str(v) for v in values]
        if result is not None and result != vector:
            raise ValueError('Result changed between repetitions')
        samples.append(elapsed)
        result = vector
    return {'status': 'ok', 'seconds': statistics.median(samples),
            'samples_seconds': samples, 'result': result}


def measure(job, case):
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker',
               json.dumps({'job': job, 'case': case})]
    # As in the paper benchmark, the limit applies to the complete worker,
    # including startup and all repetitions; reported times exclude imports.
    try:
        process = subprocess.run(command, cwd=ROOT, capture_output=True,
                                 text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {'status': 'timeout', 'seconds': '', 'result': None}
    if process.returncode:
        raise RuntimeError(process.stderr)
    result = json.loads(process.stdout)
    expected = job.get('expected_result')
    if expected is not None and result['result'] != expected:
        raise ValueError(f"Exact vector differs from the paper: {job['case']}, "
                         f"{job['index']}, {job['method']}")
    return result


def write_tables(output, figures):
    columns = ['case', 'index', 'method', 'n', 'q', 'repeats', 'status', 'seconds']
    for figure, rows in figures.items():
        with (output / (figure + '.csv')).open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        lines = ['| ' + ' | '.join(columns) + ' |',
                 '| ' + ' | '.join(['---'] * len(columns)) + ' |']
        lines += ['| ' + ' | '.join(str(row[c]) for c in columns) + ' |' for row in rows]
        (output / (figure + '.md')).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    data = json.loads((ROOT / 'data/cases.json').read_text(encoding='utf-8'))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure', choices=list(data['figures']), action='append',
                        help='restrict execution to a paper figure; default: all four')
    parser.add_argument('--case', choices=list(data['cases']), action='append',
                        help='restrict execution to a plotted case; default: all cases')
    args = parser.parse_args()
    output = ROOT / 'tables' / time.strftime('%Y%m%d_%H%M%S')
    output.mkdir(parents=True, exist_ok=False)
    figures, cache, vectors = {}, {}, {}
    for figure, jobs in data['figures'].items():
        if args.figure and figure not in args.figure:
            continue
        rows = []
        figures[figure] = rows
        for job in jobs:
            if args.case and job['case'] not in args.case:
                continue
            key = job['case'], job['index'], job['method']
            case = data['cases'][job['case']]
            if key not in cache:
                print(f"{figure}: {', '.join(key)}", flush=True)
                cache[key] = measure(job, case)
            result = cache[key]
            if result['status'] == 'ok':
                vector_key = job['case'], job['index']
                if vector_key in vectors and vectors[vector_key] != result['result']:
                    raise ValueError(f'Algorithms disagree: {vector_key}')
                vectors[vector_key] = result['result']
            row = {k: job[k] for k in ('case', 'index', 'method', 'repeats')}
            row.update(n=len(case['w']), q=case['q'], **result)
            rows.append(row)
            # Keep completed rows if a later computation fails or is interrupted.
            write_tables(output, figures)
            (output / 'results.json').write_text(json.dumps(figures, indent=2) + '\n', encoding='utf-8')
    print(f'Tables: {output}')


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        payload = json.loads(sys.argv[2])
        print(json.dumps(worker(payload['job'], payload['case'])))
    else:
        main()
