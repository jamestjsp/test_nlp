from fractions import Fraction

import cvxpy as cp
import numpy as np
import pytest

from fixed_quality_scan import ARCS, FEEDS, PRODUCTS, exact_profit, solve_fixed_q
from solve_pooling import (
    build_original,
    check_original,
    exact_witness,
    local_optimum_witnesses,
    mandatory_feasibility_certificate,
    solve_case,
)


def test_corrected_mandatory_case_has_an_exact_infeasibility_certificate():
    case = solve_case("mandatory", lower_bound=25)
    assert case["status"] == "infeasible"
    certificate = case["feasibility_certificate"]
    assert certificate["lhs_exact"] == "225"
    assert certificate["rhs_exact"] == "200"
    assert exact_witness("mandatory", 25) is None


@pytest.mark.parametrize("lower_bound", [25, 30, 49.999, 50])
def test_optional_corrected_case_keeps_global_profit_500(lower_bound):
    case = solve_case("optional", lower_bound=lower_bound)
    assert case["status"] == "globally_optimal"
    assert case["feasible_lower_bound"] == pytest.approx(500, abs=1e-7)
    assert case["absolute_gap"] <= 1e-7
    assert case["exact_witness"]["feasible"]
    assert case["exact_witness"]["profit"] == 500


def test_threshold_witness_is_feasible_in_exact_arithmetic():
    bound = Fraction(200, 9)
    flows = {arc: bound for arc in ARCS}
    flows[("d", "B", "P2")] = Fraction(11, 2) * bound
    flows[("y", "P1")] = flows[("y", "P2")] = Fraction(3, 2) * bound
    assert all(v >= bound for v in flows.values())
    assert exact_profit(flows, lower_bound=bound) == Fraction(500, 3)
    for j, capacity in [("P1", 100), ("P2", 200)]:
        assert sum(flows[("d", i, j)] for i in FEEDS) + flows[("y", j)] == capacity


def test_both_sides_of_threshold_with_independent_solver():
    threshold = 200 / 9
    assert mandatory_feasibility_certificate(threshold - 1e-5)["feasible"]
    assert not mandatory_feasibility_certificate(threshold + 1e-5)["feasible"]
    assert solve_fixed_q(2, optional=False, lower_bound=threshold - 1e-5) is not None
    assert solve_fixed_q(2, optional=False, lower_bound=threshold + 1e-5) is None


def test_cvxpys_fixed_quality_mandatory_model_is_also_infeasible():
    model = build_original("mandatory", fixed_quality=2, lower_bound=25)
    model.problem.solve(solver=cp.HIGHS)
    assert model.problem.status == cp.INFEASIBLE


@pytest.mark.parametrize("q", [1, 1.25, 1.5, 2, 2.5, 2.75, 3])
def test_independent_corrected_profile(q):
    assert solve_fixed_q(q, optional=False, lower_bound=25) is None
    assert solve_fixed_q(q, optional=True, lower_bound=25) == pytest.approx(500)


def test_corrected_checker_rejects_the_old_mandatory_witness():
    old = exact_witness("mandatory")
    result = check_original(old["direct"], old["pool_inlet"], old["pool_outlet"],
                            old["pool_sulfur"], "mandatory", lower_bound=25)
    assert not result["feasible"]
    assert result["residuals"]["flow_domain"] == 24


def test_optional_active_pool_can_attain_500_with_corrected_minimum():
    point = check_original([[50, 0], [0, 100], [0, 0]], [0, 0, 150], [50, 100],
                           2, "optional", lower_bound=25)
    assert point["feasible"]
    assert point["profit"] == 500


def test_claimed_local_values_have_exact_feasible_support_witnesses():
    cases = local_optimum_witnesses(25)
    assert [case["profit"] for case in cases] == [400, 300, 100, 0]
    for case in cases:
        point = case["witness"]
        flows = {("d", i, j): point["direct"][a][b]
                 for a, i in enumerate(FEEDS) for b, j in enumerate(PRODUCTS)}
        flows.update({("a", i): 0 for i in FEEDS})
        flows.update({("y", j): 0 for j in PRODUCTS})
        assert exact_profit(flows, lower_bound=25) == case["support_upper_bound"]


def test_local_support_bounds_follow_exact_sulfur_inequalities():
    sulfur, cost = [3, 1, 2], [6, 16, 10]
    for j, feeds, price, spec, margin, penalty in [
        (1, [1, 2], 15, Fraction(3, 2), 2, 6),
        (1, [0, 1], 15, Fraction(3, 2), Fraction(3, 2), 5),
        (0, [0, 2], 9, Fraction(5, 2), 1, 4),
    ]:
        for i in feeds:
            assert price - cost[i] == margin - penalty * (spec - sulfur[i])
        assert margin * [100, 200][j] in [400, 300, 100]


def test_unclosed_profit_bound_is_not_reported_as_global():
    case = solve_case("mandatory", lower_bound=20)
    assert case["status"] == "feasible"
    assert case["absolute_gap"] > 1e-7
    assert "gap remains open" in case["certification"]


@pytest.mark.parametrize("lower_bound", [0, -1, np.inf, np.nan])
def test_invalid_lower_bounds_are_rejected(lower_bound):
    with pytest.raises(ValueError, match="finite and positive"):
        build_original(lower_bound=lower_bound)
