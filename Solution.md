# Super pooling: globally solved with CVXPY

**Optional streams: global maximum profit 500. Mandatory streams ≥ 1: global
maximum profit 494, so the diagram's mandatory-flow problem is feasible.**

**Author's correction:** the intended minimum is 25. At that minimum,
mandatory streams are infeasible and optional streams still attain 500.
See [the correction, exact feasibility threshold, and local witnesses](Correction.md).
This document retains the original minimum-1 study below.

These results use only `Problem.md` and the supplied diagram in
`1792413976415.pdf`. No feed availability limits, pool capacity, fixed charges,
or other unshown constraints are assumed. Additional constraints in the post's
implementation could explain its different feasibility result.

| Stream domain | Problem class | Feasible profit | Proven upper bound | Exact gap |
|---|---|---:|---:|---:|
| Every stream ≥ 1 | Nonconvex continuous NLP | 494 | 494 | 0 |
| Every stream ∈ {0} ∪ [1, U] | Nonconvex MINLP | 500 | 500 | 0 |

A positive lower bound alone does **not** make a model mixed integer. The
zero-or-at-least-one choice introduces the binary decisions.

## Data and original model

The diagram supplies:

| Feed | Sulfur | Unit cost |
|---|---:|---:|
| A | 3 | 6 |
| B | 1 | 16 |
| C | 2 | 10 |

| Product | Maximum quantity | Maximum sulfur | Unit revenue |
|---|---:|---:|---:|
| P1 | 100 | 2.5 | 9 |
| P2 | 200 | 1.5 | 15 |

Sulfur numbers use the diagram's consistent concentration scale; no percentage
conversion is needed. No sulfur-removal operation or pool processing cost is shown.

Let $d_{ij}$ be direct feed-to-product flows, $a_i$ feed-to-pool flows,
$y_j$ pool-to-product flows, and $q$ the common pool sulfur concentration.
There are 11 independent stream flows: six direct, three pool inlets, and two
pool outlets. Aggregate feed, pool, and product streams follow by conservation.

With $s=(3,1,2)$, $c=(6,16,10)$, $p=(9,15)$,
$D=(100,200)$, and $\bar s=(2.5,1.5)$:

$$
\max\; \Pi=\sum_jp_jt_j-\sum_ic_i\left(a_i+\sum_jd_{ij}\right),
\qquad t_j=\sum_i d_{ij}+y_j,
$$

subject to

$$
\sum_i a_i=\sum_j y_j,\qquad
\sum_i s_i a_i=q\sum_jy_j,\qquad 1\le q\le3,
$$

$$
\sum_i s_i d_{ij}+q y_j\le\bar s_jt_j,
\qquad t_j\le D_j.
$$

Mandatory mode imposes $d_{ij},a_i,y_j\ge1$. Optional mode uses a binary
$b_f$ for each stream and the exact constraints

$$
b_f\le f\le U_fb_f,\qquad b_f\in\{0,1\}.
$$

Valid upper bounds are $U_{d_{ij}}=D_j$, $U_{y_j}=D_j$, and
$U_{a_i}=300$, from total product capacity and conservation. The aggregate
streams also obey the corresponding domain automatically.

The original bilinear equalities and product-quality constraints are non-DCP.
CVXPY's [DCP rules](https://www.cvxpy.org/tutorial/dcp/) do not accept products
of two unrestricted continuous decision variables as affine expressions.
The script checks this and never sends the non-DCP model to HiGHS.

## Exact optimal flows

Each entry below is a flow, not a concentration. The optima are not unique.
The table shows the hand-built exact witnesses (`exact_witness` in
`solve_pooling.py`). The numerical recovery in `results.json` is a different
mandatory optimum: C → P1 = 46, C → P2 = 97, Pool → P1 = 2, Pool → P2 = 1,
with the other flows as below. It also earns 494 with total pool flow 3.

| Stream | Mandatory optimum | Optional optimum |
|---|---:|---:|
| A → P1 | 51 | 50 |
| A → P2 | 1 | 0 |
| B → P1 | 1 | 0 |
| B → P2 | 101 | 100 |
| C → P1 | 47 | 50 |
| C → P2 | 96 | 100 |
| A → pool | 1 | 0 |
| B → pool | 1 | 0 |
| C → pool | 1 | 0 |
| Pool → P1 | 1 | 0 |
| Pool → P2 | 2 | 0 |

The mandatory pool has quantity 3 and sulfur $q=2$:

$$
3(1)+1(1)+2(1)=2(1+2)=6.
$$

Product balances are exact:

$$
t_1=51+1+47+1=100,\qquad
3(51)+1+2(47)+2(1)=250,
$$

$$
t_2=1+101+96+2=200,\qquad
3+101+2(96)+2(2)=300.
$$

Thus product sulfur is exactly 2.5 and 1.5. Feed quantities are
$(53,103,144)$, giving

$$
\Pi=9(100)+15(200)-6(53)-16(103)-10(144)=494.
$$

The optional solution bypasses the pool, purchases $(50,100,150)$, and earns
$3900-3400=500$. A functioning pool can also attain 500: send 150 of C into
the pool, split its outlets as $(50,100)$, set both direct C streams to zero,
and retain the direct A → P1 and B → P2 streams. All active streams exceed 1.
An idle pool has undefined physical concentration; assigning $q=2$ is harmless
because every pool flow and sulfur load is zero.

## Global upper bounds: analytical certificate

For a nonempty pool, define the effective amount of feed $i$ reaching product
$j$ as

$$
m_{ij}=d_{ij}+\frac{a_i y_j}{T},\qquad T=\sum_i a_i=\sum_j y_j.
$$

For an empty pool take $m_{ij}=d_{ij}$. This definition preserves feed costs,
product quantities, and product sulfur because the original pool is well mixed.
Let $\sigma_j=\sum_i s_i m_{ij}$ and the nonnegative sulfur slack be
$\delta_j=\bar s_jt_j-\sigma_j$.

The feed costs satisfy two exact identities:

$$
c_i=18-4s_i+2\mathbf1_{i=B}
    =22-6s_i+2\mathbf1_{i=A}.
$$

Using the first identity for P1 and the second for P2 gives

$$
\Pi=t_1+2t_2-4\delta_1-6\delta_2
       -2m_{B1}-2m_{A2}.
$$

Every subtracted term is nonnegative, so **all original-feasible solutions
have profit at most $100+2(200)=500$**. The optional witness attains this bound.

When every stream is mandatory, $d_{B1},d_{A2},a_A,a_B\ge1$, hence

$$
m_{B1}+m_{A2}
=d_{B1}+d_{A2}+\frac{a_B y_1+a_A y_2}{T}
\ge2+\frac{y_1+y_2}{T}=3.
$$

Therefore **$\Pi\le500-2(3)=494$**. The mandatory witness attains this bound.
These are global certificates for the nonconvex originals; no local solver or
initialization assumption is needed. Exact rational checks of both witnesses
are included in the tests.

## CVXPY convexification and recovery

The implementation follows the `convexify-with-cvxpy` skill:

1. Introduce nonnegative $w_{ij}$ for source-to-product flow routed through
   the pool. Impose row sums $a_i$ and column sums $y_j$, but drop the common
   composition condition $w_{ij}=a_iy_j/T$. This is an **outer relaxation**,
   not an exact reformulation. Product sulfur and profit become affine in
   $d+w$.
2. In mandatory mode add the valid inequality
   $w_{B1}+w_{A2}\ge1$, proved above. In optional mode the corresponding
   valid cut is $w_{B1}+w_{A2}\ge b_{a_A}+b_{a_B}-1$.
3. Solve the DCP LP/MILP relaxation through CVXPY/HiGHS. The mandatory raw
   relaxation estimates 496; the strengthened relaxation estimates 494.
   The optional relaxation estimates 500. Integer activation constraints
   remain exact in the optional relaxation.
4. Fix $q=2$ in the original equations and solve the resulting LP or MILP.
   This is an **exact restriction at a fixed quality**, not a general global
   reformulation. It happens to recover global optima here because its
   feasible profit matches the independent universal upper bound.
   A second solve minimizes total pool flow while preserving optimal profit
   to select a simple witness.
5. Recompute profit, flow balance, bilinear sulfur balance, product quality,
   capacity, flow domains, and activation integrality from the recovered
   original-space solution. Nothing is certified by relaxation status alone.

The numerical relaxation objectives in `results.json` are labeled estimates.
The upper-bound certificate comes from the analytical inequalities, with exact
integer witnesses. The optional numerical recovery has a gap below $10^{-12}$;
the exact optional witness has zero gap. No spatial branching is required.
The supplied local-optimum claims and solver timing comparisons are not used
as evidence for these global bounds. The local values are verified separately
in `Correction.md`; the author's solver runs and timings are not reproduced.

## Independent cross-check at fixed pool sulfur

`fixed_quality_scan.py` rebuilds the original model with SciPy only, without
CVXPY. For a fixed pool sulfur $q$ every constraint is linear, so each value of
$q$ is one LP (mandatory) or MILP (optional) solved by HiGHS through
`scipy.optimize.milp`. The tests scan 802 values of $q$ in $[1,3]$. No value
exceeds 494 or 500, and both values are reached at $q=2$. This scan is a bug
check. The analytical certificate above is the proof.

| Pool sulfur $q$ | Mandatory profit | Optional profit |
|---:|---:|---:|
| 1.000 | infeasible | 500 |
| 1.011 | −639.27 | 500 |
| 1.020 | 400 | 500 |
| 1.050 | 492.10 | 500 |
| 1.500 | 493 | 500 |
| 2.000 | 494 | 500 |
| 2.500 | 493 | 500 |
| 2.950 | 480 | 500 |
| 2.975 | 360 | 500 |
| 3.000 | infeasible | 500 |

Each row is the best profit with the pool sulfur held at one value. It is not
a local optimum of the full problem. The post's local optima (400, 100, 0) come
from the zero-or-at-least-1 version, where every row is 500. The 400 at
$q=1.02$ is unrelated to them.

In mandatory mode the best profit is $494-2|q-2|$ over most of the range, and
it falls fast near the ends. Every pool inlet is at least 1, and every product
also receives three direct streams of at least 1. Thus the pool outflow is at
most $97+197=294$, and the pool sulfur must lie in about
$[297/294,\;2.9757]$. Outside this interval the mandatory model has no
solution. The optional model can turn the pool off, so its profit stays 500
for every $q$.

## Reproduce

From this folder:

```sh
uv venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python solve_pooling.py
.venv/bin/python fixed_quality_scan.py
.venv/bin/python -m pytest -q
```

On Windows, use `.venv\Scripts\python.exe` in place of `.venv/bin/python`.

`solve_pooling.py` writes `results.json`. Tested with CVXPY 1.9.2, HiGHS
1.15.1 and SciPy 1.18.1: **19 tests passed**. Both original models pass
feasibility checks with tolerance $10^{-7}$; the mandatory numerical witness has zero residual and the
optional numerical witness has maximum residual below $10^{-12}$.

The regression suite also checks that a raw relaxation point can violate the
original product quality, that broken bilinear balances are rejected, that
optional flows between zero and one are rejected, and that an active-pool
optional solution attains 500. It checks that the optional relaxation bound is
500, and that the pool contamination cut holds for sampled well-mixed pools.

This resolves the diagram as supplied. A claim that its mandatory-flow model
is infeasible would require additional constraints or different data.
The author's subsequent change of the minimum to 25 supplies that different
data and does make the mandatory model infeasible, as proved in `Correction.md`.
