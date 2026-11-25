# Halley's method

Halley's method uses the first two derivatives of a scalar function to obtain a third-order iteration,
\f[
  x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)} \left(1 - \frac{f(x_k) f''(x_k)}{2 f'(x_k)^2}\right)^{-1} \text{,}
\f]
which can be derived by applying Newton's method to \f$f/\sqrt{|f'|}\f$. Like Newton's method, it converges only locally; the optional damping relaxes the step until the residual decreases.
