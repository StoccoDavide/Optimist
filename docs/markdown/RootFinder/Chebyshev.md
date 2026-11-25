# Chebyshev's method

Chebyshev's method is a third-order iteration for scalar functions based on the logarithmic convexity degree \f$L_k = f(x_k) f''(x_k) / f'(x_k)^2\f$:
\f[
  x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)} \left(1 + \frac{L_k}{2}\right) \text{.}
\f]
The correction is only meaningful for \f$|L_k| < 1\f$; outside this range (far from the root, where the second-order model is unreliable and the step may even vanish) the Newton step is taken instead, and the third-order convergence is recovered near the root, where \f$L_k \to 0\f$.
