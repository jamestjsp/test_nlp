from fractions import Fraction

import numpy as np
import pytest

from solve_pooling import (
    build_original,
    build_relaxation,
    check_original,
    exact_witness,
    solve_case,
    solve_checked,
)


@pytest.mark.parametrize("mode, optimum", [("mandatory", 494), ("optional", 500)])
def test_global_bound_meets_original_feasible_point(mode, optimum):
    result = solve_case(mode)
    assert result["feasible_lower_bound"] == pytest.approx(optimum, abs=1e-7)
    assert result["strengthened_relaxation_objective_estimate"] == pytest.approx(optimum, abs=1e-7)
    assert result["absolute_gap"] <= 1e-7
    assert result["recovered"]["max_residual"] <= 1e-7
    assert not result["original_is_dcp"]
    assert result["fixed_quality_is_dcp"]
    assert result["fixed_quality_is_mixed_integer"] == (mode == "optional")


@pytest.mark.parametrize("mode, optimum", [("mandatory", 494), ("optional", 500)])
def test_witness_with_exact_rational_arithmetic(mode, optimum):
    point = exact_witness(mode)
    direct = [[Fraction(v) for v in row] for row in point["direct"]]
    inlet = [Fraction(v) for v in point["pool_inlet"]]
    outlet = [Fraction(v) for v in point["pool_outlet"]]
    quality = Fraction(point["pool_sulfur"])
    sulfur = [3, 1, 2]
    costs = [6, 16, 10]
    assert sum(inlet) == sum(outlet)
    assert sum(s * f for s, f in zip(sulfur, inlet)) == quality * sum(outlet)
    production = [sum(row[j] for row in direct) + outlet[j] for j in range(2)]
    assert production == [100, 200]
    for j, spec in enumerate([Fraction(5, 2), Fraction(3, 2)]):
        assert sum(sulfur[i] * direct[i][j] for i in range(3)) + quality * outlet[j] <= spec * production[j]
    flows = [f for row in direct for f in row] + inlet + outlet
    assert all(f >= 1 or (mode == "optional" and f == 0) for f in flows)
    expenditure = sum(costs[i] * (sum(direct[i]) + inlet[i]) for i in range(3))
    assert 9 * production[0] + 15 * production[1] - expenditure == optimum


def test_pool_cut_closes_relaxation_gap():
    raw = build_relaxation("mandatory", strengthen=False)
    strengthened = build_relaxation("mandatory", strengthen=True)
    assert solve_checked(raw.problem) == pytest.approx(496)
    assert solve_checked(strengthened.problem) == pytest.approx(494)
    candidate = check_original(raw.direct.value, raw.inlet.value, raw.outlet.value,
                               float(np.array([3, 1, 2]) @ raw.inlet.value / sum(raw.inlet.value)),
                               "mandatory")
    assert not candidate["feasible"]
    assert candidate["residuals"]["product_sulfur"] > 1e-7


def test_optional_relaxation_bound_matches_certificate():
    raw = build_relaxation("optional", strengthen=False)
    strengthened = build_relaxation("optional", strengthen=True)
    assert solve_checked(raw.problem) == pytest.approx(500)
    assert solve_checked(strengthened.problem) == pytest.approx(500)


def test_pool_contamination_cut_holds_for_well_mixed_pools():
    rng = np.random.default_rng(0)
    for _ in range(10_000):
        inlet = np.array([1, 1, 0]) + rng.exponential(20, size=3)
        split = rng.uniform(0, 1)
        outlet = inlet.sum() * np.array([split, 1 - split])
        routed = np.outer(inlet, outlet) / outlet.sum()
        assert routed[1, 0] + routed[0, 1] >= 1 - 1e-12


def test_original_checker_rejects_broken_pool_sulfur_balance():
    point = exact_witness("mandatory")
    candidate = check_original(point["direct"], point["pool_inlet"], point["pool_outlet"],
                               1.9, "mandatory")
    assert not candidate["feasible"]
    assert candidate["residuals"]["pool_sulfur_balance"] == pytest.approx(0.3)


def test_original_checker_rejects_semicontinuous_flow_below_one():
    candidate = check_original([[0.5, 0], [0, 0], [0, 0]], [0, 0, 0], [0, 0], 2, "optional")
    assert not candidate["feasible"]
    assert candidate["residuals"]["flow_domain"] == 0.5


def test_optional_witness_can_use_an_active_pool():
    candidate = check_original([[50, 0], [0, 100], [0, 0]], [0, 0, 150], [50, 100], 2, "optional")
    assert candidate["feasible"]
    assert candidate["profit"] == 500


def test_mandatory_streams_are_continuous_not_mixed_integer():
    assert not build_original("mandatory").problem.is_mixed_integer()
    assert build_original("optional").problem.is_mixed_integer()


def test_nondcp_original_never_sent_to_highs():
    with pytest.raises(ValueError, match="Only DCP"):
        solve_checked(build_original().problem)


def test_invalid_domain_is_not_silently_reinterpreted():
    with pytest.raises(ValueError, match="Unknown flow domain"):
        build_original("invalid")
