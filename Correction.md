# Author's correction: minimum stream flow 25

The author clarified that the intended common stream lower bound is **25**, not
1. This changes the result: **mandatory streams are infeasible; optional
streams retain globally optimal profit 500**. The original result at minimum 1
remains valid.

The parenthetical “≥25 and <50” is interpreted as a range for the common lower
bound $L$, not as an upper bound of 50 on actual stream flows. The latter would
be a different model. The conclusions below hold throughout $25\le L<50$;
the optional witnesses also remain valid at $L=50$.

| Common minimum $L$ | Mandatory streams | Zero-or-at-least-$L$ streams |
|---|---|---|
| 1 | Global profit 494 | Global profit 500 |
| 25 | Infeasible | Global profit 500 |

## Exact mandatory-flow feasibility threshold

Use the notation from [the original solution](Solution.md). Set
$T=a_A+a_B+a_C=y_1+y_2$ and $r=y_1/T$, the fraction of pool output sent to P1.
For mandatory flows, $T\ge3L>0$ and $0\le r\le1$.

P1 has three direct streams and one pool stream. Its capacity implies

$$
100\ge t_1\ge3L+y_1=3L+rT\ge3L(1+r).
$$

Let $\sigma_2$ be the sulfur load delivered to P2. The pool's sulfur excess
above feed B's concentration of 1 is $2a_A+a_C\ge3L$. Consequently

$$
\sigma_2-t_2
=2d_{A2}+d_{C2}+(2a_A+a_C)(1-r)
\ge3L(2-r).
$$

P2's quality and capacity constraints imply
$\sigma_2-t_2\le0.5t_2\le100$. Adding the two bounds eliminates $r$:

$$
\boxed{9L\le200.}
$$

At **$L=25$**, this requires **$225\le200$**, a contradiction. Equivalently,
P1 requires $r\le1/3$ while P2 requires $r\ge2/3$. This proves global
infeasibility without relying on a solver status or a quality grid.

The threshold is tight. For any $0<L\le200/9$, set

- Every pool inlet to $L$ and both pool outlets to $3L/2$, giving pool sulfur 2.
- All three direct P1 streams to $L$.
- Direct A → P2 and C → P2 to $L$, and B → P2 to $11L/2$.

All streams satisfy their lower bounds. Product quantities are $9L/2$ and
$9L$, with sulfur concentrations 2 and 1.5. Thus the mandatory model is
feasible **if and only if $L\le200/9\approx22.2222$**. At the threshold,
both product capacities are reached. The threshold witness is checked with
exact rational arithmetic in the tests.

## Optional-flow global optimum: 500

Use the original bypass witness:

| Stream | Flow |
|---|---:|
| A → P1 | 50 |
| C → P1 | 50 |
| B → P2 | 100 |
| C → P2 | 100 |
| All other streams | 0 |

Every active flow is at least 50. Product quantities are 100 and 200, and
sulfur concentrations are 2.5 and 1.5. Revenue is 3900, feed cost 3400,
and profit **500**. The original universal profit bound $\Pi\le500$ is
unchanged, so this is globally optimal for every $0<L\le50$.

An active-pool witness is also available: send 150 of C into the pool, split
its output as 50 to P1 and 100 to P2, retain direct A → P1 = 50 and
B → P2 = 100, and turn off the remaining streams.

## The reported local values exist

The following witnesses all satisfy the corrected optional-flow model. Pool
flows and every unlisted direct stream are zero.

| Profit | Active direct streams | Global upper bound on this fixed support |
|---:|---|---|
| 400 | B → P2 = 100; C → P2 = 100 | $\Pi\le2t_2\le400$ |
| 300 | A → P2 = 50; B → P2 = 150 | $\Pi\le1.5t_2\le300$ |
| 100 | A → P1 = 50; C → P1 = 50 | $\Pi\le t_1\le100$ |
| 0 | All streams off | The all-zero flow vector is isolated |

For the 400 and 100 rows, the stated bounds follow from the cost identities
in the original proof. For the 300 row, writing $x=d_{A2}$ and
$t=t_2$ gives $\Pi=-t+10x$; the sulfur limit forces $x\le t/4$, hence
$\Pi\le1.5t\le300$.

Every active flow in the table is at least 50. In a flow-space neighborhood
of radius $L/2$, an unused stream cannot activate because it must jump from
zero to at least $L$. Therefore each support optimum is a local maximum of
the disjunctive problem. These are **nonstrict** local maxima when the
physically undefined idle-pool concentration is retained as a variable.
The zero-profit point is isolated in flow space, so its locality is trivial.

This verifies the existence of all four reported local values; it does not
assert that they exhaust every possible local optimum or reproduce the
author's particular solver solutions. These support arguments also apply at
$L=1$.

## Reproduce and validate

```sh
.venv/bin/python solve_pooling.py --lower-bound 25 --output results_lower_bound_25.json
.venv/bin/python fixed_quality_scan.py --lower-bound 25
.venv/bin/python -m pytest -q
```

The CVXPY model, original-space feasibility checker, activation bounds, and
contamination cuts now accept `lower_bound`; the independent SciPy model
accepts the same parameter. The command-line option is `--lower-bound`.
Default 1 preserves the original study. For other feasible lower bounds, a
fixed-quality solution with an unclosed analytical profit gap is explicitly
reported as feasible rather than globally optimal.

`results_lower_bound_25.json` records the exact mandatory infeasibility
certificate, optional global bound and recovered solution, and four local
witnesses. CVXPY/HiGHS also reports the mandatory $q=2$ restriction infeasible;
the independent SciPy model agrees at all tested pool qualities. Those
numerical checks support the implementation; the analytical contradiction is
the global proof. **43 tests pass**, including the original 19 regressions.
