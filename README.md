# Haverly super pooling with CVXPY

Globally certified solutions to the supplied super-pooling problem, using
CVXPY and HiGHS.

| Stream domain | Global maximum profit |
|---|---:|
| Every stream at least 1 | 494 |
| Every stream zero or at least 1 | 500 |

The mandatory-flow model is feasible for the diagram as supplied. Analytical
upper bounds match exact feasible witnesses for both domains.

- [Problem statement](Problem.md) and [source diagram](1792413976415.pdf)
- [Formulation, optimal flows, and global proofs](Solution.md)
- [CVXPY implementation](solve_pooling.py)
- [Numerical results](results.json) and [regression tests](test_pooling.py)

## Run

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```sh
uv venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python solve_pooling.py
.venv/bin/python -m pytest -q
```

Validated with CVXPY 1.9.2 and HiGHS 1.15.1; 11 tests pass. The solver writes
`results.json` and validates the original nonlinear balances. Relaxation
objective estimates are reported separately from analytical certificates.
