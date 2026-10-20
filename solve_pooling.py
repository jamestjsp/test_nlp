"""Solve the supplied super-pooling diagram and check the original equations."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import cvxpy as cp
import highspy
import numpy as np


SULFUR = np.array([3.0, 1.0, 2.0])
COST = np.array([6.0, 16.0, 10.0])
PRICE = np.array([9.0, 15.0])
CAPACITY = np.array([100.0, 200.0])
SPEC = np.array([2.5, 1.5])
MODES = ("mandatory", "optional")
TOL = 1e-7


@dataclass
class Model:
    problem: cp.Problem
    direct: cp.Variable
    inlet: cp.Expression
    outlet: cp.Expression
    quality: cp.Expression | float | None = None
    routed: cp.Variable | None = None
    activation: cp.Variable | None = None


def validate_lower_bound(lower_bound):
    if not np.isfinite(lower_bound) or lower_bound <= 0:
        raise ValueError("Stream lower bound must be finite and positive")


def mandatory_feasibility_certificate(lower_bound):
    validate_lower_bound(lower_bound)
    lhs = 9 * Fraction(str(lower_bound))
    return {
        "feasible": lhs <= 200,
        "necessary_and_sufficient_condition": "9 * lower_bound <= 200",
        "threshold_exact": "200/9",
        "threshold": 200 / 9,
        "lhs_exact": str(lhs),
        "rhs_exact": "200",
    }


def flow_domains(direct, inlet, outlet, mode, lower_bound=1.0):
    validate_lower_bound(lower_bound)
    if mode not in MODES:
        raise ValueError(f"Unknown flow domain: {mode}")
    if mode == "mandatory":
        return [direct >= lower_bound, inlet >= lower_bound, outlet >= lower_bound], None
    activation = cp.Variable(11, boolean=True, name="stream_on")
    flows = cp.hstack([cp.vec(direct, order="C"), inlet, outlet])
    upper = np.concatenate([np.tile(CAPACITY, 3), np.full(3, 300.0), CAPACITY])
    return [flows >= lower_bound * activation, flows <= cp.multiply(upper, activation)], activation


def build_original(mode="mandatory", fixed_quality=None, lower_bound=1.0):
    """The full nonconvex model, or an exact restriction at fixed pool quality."""
    direct = cp.Variable((3, 2), nonneg=True, name="direct")
    inlet = cp.Variable(3, nonneg=True, name="pool_inlet")
    outlet = cp.Variable(2, nonneg=True, name="pool_outlet")
    quality = cp.Variable(name="pool_sulfur") if fixed_quality is None else float(fixed_quality)
    constraints, activation = flow_domains(direct, inlet, outlet, mode, lower_bound)
    production = cp.sum(direct, axis=0) + outlet
    constraints += [
        quality >= 1,
        quality <= 3,
        cp.sum(inlet) == cp.sum(outlet),
        SULFUR @ inlet == quality * cp.sum(outlet),
        SULFUR @ direct + quality * outlet <= cp.multiply(SPEC, production),
        production <= CAPACITY,
    ]
    profit = PRICE @ production - COST @ (cp.sum(direct, axis=1) + inlet)
    return Model(cp.Problem(cp.Maximize(profit), constraints), direct, inlet, outlet,
                 quality, activation=activation)


def build_relaxation(mode="mandatory", strengthen=True, lower_bound=1.0):
    """Drop common pool composition; optionally add its valid contamination cut."""
    direct = cp.Variable((3, 2), nonneg=True, name="direct")
    routed = cp.Variable((3, 2), nonneg=True, name="pool_source_to_product")
    inlet = cp.sum(routed, axis=1)
    outlet = cp.sum(routed, axis=0)
    constraints, activation = flow_domains(direct, inlet, outlet, mode, lower_bound)
    delivered = direct + routed
    production = cp.sum(delivered, axis=0)
    constraints += [
        production <= CAPACITY,
        SULFUR @ delivered <= cp.multiply(SPEC, production),
    ]
    if strengthen:
        if mode == "mandatory":
            constraints += [routed[1, 0] + routed[0, 1] >= lower_bound]
        else:
            constraints += [routed[1, 0] + routed[0, 1]
                            >= lower_bound * (activation[6] + activation[7] - 1)]
    profit = PRICE @ production - cp.sum(cp.multiply(COST[:, None], delivered))
    return Model(cp.Problem(cp.Maximize(profit), constraints), direct, inlet, outlet,
                 routed=routed, activation=activation)


def solve_checked(problem):
    if not problem.is_dcp():
        raise ValueError("Only DCP subproblems are sent to HiGHS")
    problem.solve(solver=cp.HIGHS, highs_options={
        "mip_rel_gap": 0.0,
        "mip_abs_gap": 1e-9,
        "primal_feasibility_tolerance": 1e-9,
        "dual_feasibility_tolerance": 1e-9,
        "mip_feasibility_tolerance": 1e-9,
    })
    if problem.status != cp.OPTIMAL:
        raise RuntimeError(f"HiGHS returned {problem.status}")
    return float(problem.value)


def check_original(direct, inlet, outlet, quality, mode, lower_bound=1.0):
    validate_lower_bound(lower_bound)
    direct, inlet, outlet = (np.asarray(x, dtype=float) for x in (direct, inlet, outlet))
    if direct.shape != (3, 2) or inlet.shape != (3,) or outlet.shape != (2,):
        raise ValueError("Incorrect flow dimensions")
    if mode not in MODES:
        raise ValueError(f"Unknown flow domain: {mode}")
    if not all(np.isfinite(x).all() for x in (direct, inlet, outlet, quality)):
        raise ValueError("Nonfinite candidate")
    flows = np.concatenate([direct.ravel(), inlet, outlet])
    production = direct.sum(axis=0) + outlet
    sulfur_load = SULFUR @ direct + quality * outlet
    if mode == "mandatory":
        domain_violation = float(np.maximum(lower_bound - flows, 0).max())
    else:
        domain_violation = float(np.minimum(np.abs(flows), np.maximum(lower_bound - flows, 0)).max())
    residuals = {
        "nonnegative": float(np.maximum(-flows, 0).max()),
        "flow_domain": domain_violation,
        "pool_flow_balance": float(abs(inlet.sum() - outlet.sum())),
        "pool_sulfur_balance": float(abs(SULFUR @ inlet - quality * outlet.sum())),
        "product_sulfur": float(np.maximum(sulfur_load - SPEC * production, 0).max()),
        "product_capacity": float(np.maximum(production - CAPACITY, 0).max()),
        "pool_quality_bounds": float(max(1 - quality, quality - 3, 0)),
    }
    total_pool = outlet.sum()
    routed = np.outer(inlet, outlet) / total_pool if total_pool > 0 else np.zeros((3, 2))
    delivered = direct + routed
    effective_sulfur = SULFUR @ delivered
    quality_slack = SPEC * production - effective_sulfur
    profit = float(PRICE @ production - COST @ (direct.sum(axis=1) + inlet))
    loss = float((CAPACITY[0] - production[0]) + 2 * (CAPACITY[1] - production[1])
                 + 4 * quality_slack[0] + 6 * quality_slack[1]
                 + 2 * (delivered[1, 0] + delivered[0, 1]))
    residuals["profit_identity"] = abs(500 - profit - loss)
    return {
        "feasible": max(residuals.values()) <= TOL,
        "profit": profit,
        "direct": direct.tolist(),
        "pool_inlet": inlet.tolist(),
        "pool_outlet": outlet.tolist(),
        "pool_sulfur": float(quality),
        "production": production.tolist(),
        "product_sulfur": np.divide(sulfur_load, production,
                                    out=np.zeros(2), where=production > 0).tolist(),
        "effective_deliveries": delivered.tolist(),
        "loss_from_500": loss,
        "residuals": residuals,
        "max_residual": max(residuals.values()),
    }


def exact_witness(mode, lower_bound=1.0):
    validate_lower_bound(lower_bound)
    if mode == "mandatory":
        if lower_bound == 1:
            return check_original([[51, 1], [1, 101], [47, 96]], [1, 1, 1], [1, 2], 2, mode)
        if not mandatory_feasibility_certificate(lower_bound)["feasible"]:
            return None
        bound = lower_bound
        return check_original([[bound, bound], [bound, 5.5 * bound], [bound, bound]],
                              [bound] * 3, [1.5 * bound] * 2, 2, mode, bound)
    if mode == "optional":
        direct = [[50, 0], [0, 100], [50, 100]] if lower_bound <= 50 else np.zeros((3, 2))
        return check_original(direct, [0, 0, 0], [0, 0], 2, mode, lower_bound)
    raise ValueError(f"Unknown flow domain: {mode}")


def local_optimum_witnesses(lower_bound=25.0):
    supports = [
        (400, [[0, 0], [0, 100], [0, 100]], "B,C to P2: profit <= 2*t2 <= 400"),
        (300, [[0, 50], [0, 150], [0, 0]], "A,B to P2: profit <= 1.5*t2 <= 300"),
        (100, [[50, 0], [0, 0], [50, 0]], "A,C to P1: profit <= t1 <= 100"),
        (0, np.zeros((3, 2)), "All streams off: isolated flow vector"),
    ]
    results = []
    for value, direct, proof in supports:
        point = check_original(direct, [0, 0, 0], [0, 0], 2, "optional", lower_bound)
        if not point["feasible"] or abs(point["profit"] - value) > TOL:
            raise ValueError("Local witness does not satisfy this lower bound")
        results.append({"profit": value, "support_upper_bound": value,
                        "support_certificate": proof, "neighborhood_radius": lower_bound / 2,
                        "witness": point})
    return results


def solve_case(mode, lower_bound=1.0):
    original = build_original(mode, lower_bound=lower_bound)
    feasibility = mandatory_feasibility_certificate(lower_bound)
    if mode == "mandatory" and not feasibility["feasible"]:
        return {
            "flow_domain": mode, "lower_bound": lower_bound,
            "status": "infeasible", "original_is_dcp": original.problem.is_dcp(),
            "certification": "Analytical infeasibility: 9 * lower_bound > 200",
            "feasibility_certificate": feasibility,
        }
    relaxation = build_relaxation(mode, lower_bound=lower_bound)
    raw_relaxation = build_relaxation(mode, strengthen=False, lower_bound=lower_bound)
    raw_value = solve_checked(raw_relaxation.problem)
    relaxation_value = solve_checked(relaxation.problem)
    restriction = build_original(mode, fixed_quality=2, lower_bound=lower_bound)
    fixed_quality_profit = solve_checked(restriction.problem)
    secondary = cp.Problem(cp.Minimize(cp.sum(restriction.inlet)),
                           restriction.problem.constraints
                           + [restriction.problem.objective.expr >= fixed_quality_profit])
    solve_checked(secondary)
    recovered = check_original(restriction.direct.value, restriction.inlet.value,
                               restriction.outlet.value, 2, mode, lower_bound)
    witness = exact_witness(mode, lower_bound)
    analytic_upper = 500.0 - 6 * lower_bound if mode == "mandatory" else 500.0
    if not recovered["feasible"] or not witness["feasible"]:
        raise AssertionError("Candidate failed original-equation validation")
    if recovered["profit"] > analytic_upper + TOL or relaxation_value > analytic_upper + TOL:
        raise AssertionError("Profit exceeds the analytic upper bound")
    gap = max(0.0, analytic_upper - recovered["profit"])
    if abs(witness["profit"] - analytic_upper) <= TOL:
        certification = "Analytic upper bound matched by exact rational feasible witness"
    elif gap <= TOL:
        certification = "Analytic upper bound matched numerically within tolerance"
    else:
        certification = "Feasible fixed-quality restriction; global bound gap remains open"
    binaries = restriction.activation
    integrality = 0.0 if binaries is None else float(np.abs(binaries.value - np.rint(binaries.value)).max())
    if integrality > TOL:
        raise AssertionError("Activation variables are nonintegral")
    return {
        "flow_domain": mode,
        "lower_bound": lower_bound,
        "status": "globally_optimal" if gap <= TOL else "feasible",
        "original_is_dcp": original.problem.is_dcp(),
        "original_is_mixed_integer": original.problem.is_mixed_integer(),
        "relaxation_is_dcp": relaxation.problem.is_dcp(),
        "fixed_quality_is_dcp": restriction.problem.is_dcp(),
        "fixed_quality_is_mixed_integer": restriction.problem.is_mixed_integer(),
        "raw_relaxation_objective_estimate": raw_value,
        "strengthened_relaxation_objective_estimate": relaxation_value,
        "analytic_global_upper_bound": analytic_upper,
        "feasible_lower_bound": recovered["profit"],
        "absolute_gap": gap,
        "certification": certification,
        "solver_status": restriction.problem.status,
        "activation_integrality_residual": integrality,
        "recovered": recovered,
        "exact_witness": witness,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results.json"))
    parser.add_argument("--lower-bound", type=float, default=1.0)
    args = parser.parse_args()
    results = {
        "cvxpy_version": cp.__version__,
        "highs_version": highspy.Highs().version(),
        "source": "Problem.md and 1792413976415.pdf; no unshown feed or pool capacities assumed",
        "lower_bound": args.lower_bound,
        "cases": [solve_case(mode, args.lower_bound) for mode in MODES],
    }
    if args.lower_bound <= 50:
        results["local_optimum_witnesses"] = local_optimum_witnesses(args.lower_bound)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    for case in results["cases"]:
        if case["status"] == "infeasible":
            print(f"{case['flow_domain']}: infeasible (9 * lower_bound > 200)")
            continue
        print(f"{case['flow_domain']}: profit={case['feasible_lower_bound']:.9f}, "
              f"upper={case['analytic_global_upper_bound']:.9f}, "
              f"gap={case['absolute_gap']:.2e}, "
              f"max_residual={case['recovered']['max_residual']:.2e}")


if __name__ == "__main__":
    main()
