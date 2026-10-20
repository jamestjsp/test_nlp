# Haverly super pooling with CVXPY

Globally certified solutions to the supplied super-pooling problem, using
CVXPY and HiGHS.

| Stream domain | Global maximum profit |
|---|---:|
| Every stream at least 1 | 494 |
| Every stream zero or at least 1 | 500 |
| Every stream at least 25 (author's correction) | Infeasible |
| Every stream zero or at least 25 | 500 |

The mandatory-flow model is feasible for the diagram as supplied. Analytical
upper bounds match exact feasible witnesses for both domains.

The author corrected the intended minimum to 25. See [the correction and
exact infeasibility proof](Correction.md). Mandatory feasibility ends at
$L=200/9$; optional streams still attain 500. The stated local values
400, 300, 100, and 0 have verified support-specific witnesses.

- [Problem statement](Problem.md) and [source diagram](1792413976415.pdf)
- [Formulation, optimal flows, and global proofs](Solution.md)
- [CVXPY implementation](solve_pooling.py)
- [Numerical results](results.json) and [regression tests](test_pooling.py)
- [Independent SciPy cross-check at fixed pool sulfur](fixed_quality_scan.py)
- [Corrected minimum-25 results](results_lower_bound_25.json)

## Run

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```sh
uv venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python solve_pooling.py
.venv/bin/python fixed_quality_scan.py
.venv/bin/python solve_pooling.py --lower-bound 25 --output results_lower_bound_25.json
.venv/bin/python fixed_quality_scan.py --lower-bound 25
.venv/bin/python -m pytest -q
```

On Windows, use `.venv\Scripts\python.exe` in place of `.venv/bin/python`.

Validated with CVXPY 1.9.2, HiGHS 1.15.1 and SciPy 1.18.1; 43 tests pass.
The solver writes `results.json` and validates the original nonlinear balances. Relaxation
objective estimates are reported separately from analytical certificates.
