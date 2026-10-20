"""Independent fixed-pool-sulfur scan of the super-pooling model using SciPy only.

For a fixed pool sulfur q the original bilinear model becomes an LP
(mandatory streams) or MILP (zero-or-at-least-one streams). The scan is a
cross-check on the CVXPY results; the analytical certificate is the proof.
"""

import argparse
from fractions import Fraction

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

SULFUR = {"A": 3, "B": 1, "C": 2}
COST = {"A": 6, "B": 16, "C": 10}
PRICE = {"P1": 9, "P2": 15}
CAP = {"P1": 100, "P2": 200}
SPEC = {"P1": Fraction(5, 2), "P2": Fraction(3, 2)}
FEEDS = ("A", "B", "C")
PRODUCTS = ("P1", "P2")

DIRECT = [(i, j) for i in FEEDS for j in PRODUCTS]
ARCS = [("d", i, j) for i, j in DIRECT] + [("a", i) for i in FEEDS] + [("y", j) for j in PRODUCTS]
N = len(ARCS)


PROFILE_Q = [1.0, 1.011, 1.02, 1.05, 1.5, 2.0, 2.5, 2.95, 2.975, 3.0]


def arc_upper(arc):
    if arc[0] == "d":
        return CAP[arc[2]]
    if arc[0] == "y":
        return CAP[arc[1]]
    return sum(CAP.values())


def solve_fixed_q(q, optional):
    idx = {arc: k for k, arc in enumerate(ARCS)}
    nvar = 2 * N if optional else N
    c = np.zeros(nvar)
    for arc in ARCS:
        k = idx[arc]
        if arc[0] == "d":
            c[k] = -(PRICE[arc[2]] - COST[arc[1]])
        elif arc[0] == "a":
            c[k] = COST[arc[1]]
        else:
            c[k] = -PRICE[arc[1]]

    rows, lo, hi = [], [], []

    def add(row, lower, upper):
        rows.append(row)
        lo.append(lower)
        hi.append(upper)

    row = np.zeros(nvar)
    for i in FEEDS:
        row[idx[("a", i)]] = 1
    for j in PRODUCTS:
        row[idx[("y", j)]] = -1
    add(row, 0, 0)

    row = np.zeros(nvar)
    for i in FEEDS:
        row[idx[("a", i)]] = SULFUR[i] - q
    add(row, 0, 0)

    for j in PRODUCTS:
        quality = np.zeros(nvar)
        volume = np.zeros(nvar)
        for i in FEEDS:
            quality[idx[("d", i, j)]] = SULFUR[i] - float(SPEC[j])
            volume[idx[("d", i, j)]] = 1
        quality[idx[("y", j)]] = q - float(SPEC[j])
        volume[idx[("y", j)]] = 1
        add(quality, -np.inf, 0)
        add(volume, -np.inf, CAP[j])

    if optional:
        for k, arc in enumerate(ARCS):
            row = np.zeros(nvar)
            row[k], row[N + k] = 1, -1
            add(row, 0, np.inf)
            row = np.zeros(nvar)
            row[k], row[N + k] = 1, -arc_upper(arc)
            add(row, -np.inf, 0)
        lb = np.zeros(nvar)
        ub = np.array([arc_upper(a) for a in ARCS] + [1] * N, dtype=float)
        integrality = np.array([0] * N + [1] * N)
    else:
        lb = np.ones(nvar)
        ub = np.array([arc_upper(a) for a in ARCS], dtype=float)
        integrality = np.zeros(nvar)

    res = milp(
        c,
        constraints=LinearConstraint(np.array(rows), lo, hi),
        bounds=Bounds(lb, ub),
        integrality=integrality,
    )
    if res.status != 0:
        return None
    return -res.fun


def exact_profit(flows):
    """Validate a witness against the original bilinear model in exact arithmetic."""
    f = {k: Fraction(v) for k, v in flows.items()}
    assert all(v == 0 or v >= 1 for v in f.values())
    pool_in = sum(f[("a", i)] for i in FEEDS)
    pool_out = sum(f[("y", j)] for j in PRODUCTS)
    assert pool_in == pool_out
    pool_sulfur_load = sum(SULFUR[i] * f[("a", i)] for i in FEEDS)
    q = pool_sulfur_load / pool_in if pool_in else Fraction(0)
    assert q * pool_out == pool_sulfur_load
    revenue = 0
    for j in PRODUCTS:
        t = sum(f[("d", i, j)] for i in FEEDS) + f[("y", j)]
        load = sum(SULFUR[i] * f[("d", i, j)] for i in FEEDS) + q * f[("y", j)]
        assert t <= CAP[j]
        assert load <= SPEC[j] * t
        revenue += PRICE[j] * t
    purchases = sum(
        COST[i] * (f[("a", i)] + sum(f[("d", i, j)] for j in PRODUCTS)) for i in FEEDS
    )
    return revenue - purchases


def witness(direct, pool_in, pool_out):
    flows = {("d", i, j): v for (i, j), v in direct.items()}
    flows.update({("a", i): v for i, v in pool_in.items()})
    flows.update({("y", j): v for j, v in pool_out.items()})
    return flows


MANDATORY_WITNESS = witness(
    {("A", "P1"): 51, ("A", "P2"): 1, ("B", "P1"): 1, ("B", "P2"): 101, ("C", "P1"): 47, ("C", "P2"): 96},
    {"A": 1, "B": 1, "C": 1},
    {"P1": 1, "P2": 2},
)

OPTIONAL_WITNESS = witness(
    {("A", "P1"): 50, ("A", "P2"): 0, ("B", "P1"): 0, ("B", "P2"): 100, ("C", "P1"): 50, ("C", "P2"): 100},
    {"A": 0, "B": 0, "C": 0},
    {"P1": 0, "P2": 0},
)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=float, nargs="+", default=PROFILE_Q)
    args = parser.parse_args()
    print("| Pool sulfur q | Mandatory profit | Optional profit |")
    print("|---:|---:|---:|")
    for q in args.q:
        cells = []
        for optional in (False, True):
            value = solve_fixed_q(q, optional)
            cells.append("infeasible" if value is None else f"{value:.2f}")
        print(f"| {q:.4f} | {cells[0]} | {cells[1]} |")


if __name__ == "__main__":
    main()
