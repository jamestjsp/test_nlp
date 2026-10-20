import numpy as np
import pytest

from fixed_quality_scan import (
    COST,
    FEEDS,
    MANDATORY_WITNESS,
    OPTIONAL_WITNESS,
    SULFUR,
    exact_profit,
    solve_fixed_q,
)


def test_mandatory_witness_is_feasible_with_profit_494():
    assert all(v >= 1 for v in MANDATORY_WITNESS.values())
    assert exact_profit(MANDATORY_WITNESS) == 494


def test_optional_witness_has_profit_500():
    assert exact_profit(OPTIONAL_WITNESS) == 500


def test_cost_identities_behind_certificate():
    for i in FEEDS:
        assert COST[i] == 18 - 4 * SULFUR[i] + 2 * (i == "B")
        assert COST[i] == 22 - 6 * SULFUR[i] + 2 * (i == "A")


Q_GRID = np.unique(np.concatenate([np.linspace(1, 3, 801), [2.0]]))


@pytest.mark.parametrize("optional,bound", [(False, 494), (True, 500)])
def test_fixed_q_scan_never_beats_certificate(optional, bound):
    values = [solve_fixed_q(q, optional) for q in Q_GRID]
    feasible = [v for v in values if v is not None]
    assert feasible
    assert max(feasible) <= bound + 1e-6
    assert solve_fixed_q(2.0, optional) == pytest.approx(bound, abs=1e-6)


def test_mandatory_pool_quality_range_is_limited_by_inlet_lower_bounds():
    assert solve_fixed_q(1.0, optional=False) is None
    assert solve_fixed_q(3.0, optional=False) is None
    q_min = 297 / 294
    assert solve_fixed_q(q_min - 1e-4, optional=False) is None
    assert solve_fixed_q(q_min + 1e-4, optional=False) is not None
    assert all(solve_fixed_q(q, optional=False) is not None for q in np.linspace(1.02, 2.97, 50))
