#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Run the whole test-suite without pytest (e.g. on a cluster node): executes
every tests/test_*.py through its own standalone runner in a subprocess and
prints a summary table.

    python tests/run_all.py               # everything
    python tests/run_all.py rpa trans     # only files whose name matches
    python tests/run_all.py -q            # per-file summary lines only

Exit code 0 iff every file passed.  With pytest installed, `pytest` from the
repository root remains the recommended runner (parallel-friendly, --runslow).
"""
import argparse
import os
import subprocess
import sys
import time

TESTS = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument('filters', nargs='*',
                    help='substring filters on the test file names')
    ap.add_argument('-q', '--quiet', action='store_true',
                    help='suppress per-test output, keep the summary')
    args = ap.parse_args()

    files = sorted(f for f in os.listdir(TESTS)
                   if f.startswith('test_') and f.endswith('.py'))
    if args.filters:
        files = [f for f in files if any(s in f for s in args.filters)]
    if not files:
        print('no test files matched')
        return 1

    results = []
    for f in files:
        print(f"=== {f} ===", flush=True)
        t0 = time.time()
        proc = subprocess.run(
            [sys.executable, os.path.join(TESTS, f)],
            stdout=(subprocess.DEVNULL if args.quiet else None),
            stderr=None,
        )
        results.append((f, proc.returncode, time.time() - t0))

    print('\n' + '-' * 60)
    nfail = 0
    for f, rc, dt in results:
        status = 'PASS' if rc == 0 else 'FAIL'
        nfail += rc != 0
        print(f"  {status}  {f:38s} ({dt:6.1f}s)")
    print(f"\n{len(results) - nfail}/{len(results)} files passed")
    return 0 if nfail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
