# Brent's method for minimization

The `Brent` class minimizes a function of one variable by combining golden-section steps with successive parabolic interpolation through the three best points (Brent, 1973). Parabolic steps are accepted only when they fall inside the current bracket and are shorter than half the step before the last, which guarantees the robustness of golden-section search with superlinear convergence on smooth functions.

The initial bracket \f$a < b < c\f$ with \f$f(b) < f(a)\f$ and \f$f(b) \leq f(c)\f$ is the bounded interval if finite bounds are set, or is found by golden expansion downhill from the initial guess; if no bracket is found (e.g., the function is unbounded below) the solver reports failure.
