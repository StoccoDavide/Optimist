# Golden-section search

The `GoldenSection` class minimizes a unimodal function of one variable by shrinking the bracket \f$[a, c]\f$ with two interior points placed in golden ratio, so that one of them is reused at every iteration and the interval shrinks by \f$1/\varphi \approx 0.618\f$ per function evaluation. The bracket is obtained as for Brent's method. Convergence is linear but unconditional.
