Can we apply our quick detection, identification and elimination infeasibility handling to the Haverly *super pooling problem*?

In our last post we described our IMPL(c) SLPQPE solver's infeasibility handling capabilities applied to an industrial MINLP crude-oil blend scheduling optimization (CBSO) problem with discrete and continuous tank-to-tank transfer degrees-of-freedom to find globally feasible quantity-logic-quality solutions.

Here we apply our IMPL-SLPQPE solver to what we call the *super pooling problem* extending the well-known and non-convex classic Haverly pooling problem (1978) by adding all possible streams of raw materials A, B and C to the products P1 and P2 transshipped or pooled, blended or mixed via the pool unit-operation shown attached where all stream flows have a lower bound 1.0 extending the problem to be a MINLP.

Solving the problem as-is using the NLP only solvers IPOPT 3.14.20, UNO 2.9.1 and IMPL-SLPQPE-HIGHS 1.15.1 without our zero-continuous variables (ZCV's), we get the following results:

a. IPOPT takes 4.20-sec, 3000 iterations and does not detect an infeasibility (Maximum Number of Iterations Exceeded),

b. UNO takes 0.50-sec, 507 iterations and does not detect an infeasibility (Algorithmic Error) and

c. SLPQPE-HIGHS takes 0.04-sec, 4 iterations and is either unconverged or detects an infeasibility (Converged or Infeasible) identifying the pool's bilinear balance for flow times sulfur (S) as the worst constraint closure.

As IPOPT and UNO do not (as yet) have the capability to convert all stream flows to be zero-continuous variables (ZCV's), running SLPQPE-HIGHS from several initial-value restarts as the super pooling problem is non-convex, we find one global optimum of 500.0 and three (3) local optimums: 400.0, 100.0 and 0.0.  The global optimum in the classic pooling problem is 400.0 with two (2) local optimums: 100.0 and 0.0.

Although we have not been able to conclusively detect and identify the infeasibilities caused by the non-zero lower bounds on all streams, we are able to detect non-convergence quicker than both IPOPT and UNO with SLPQPE-HIGHS and most importantly we are able to eliminate all lower bound infeasibilities via our zero-continuous variables (ZCV's).

In addition, this new super pooling problem can be used as a small but representative (non-convex) MINLP problem to explore other interesting solution algorithms :-)
