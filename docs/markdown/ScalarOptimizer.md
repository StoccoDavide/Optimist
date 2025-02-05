# Scalar optimizer

This section describes the scalar optimizers implemented in `Optimist`. The available optimizers are derivative and non-derivative methods. Derivative methods employ the function's derivative to find the minimum with high accuracy, while non-derivative methods approximate the derivative for improved efficiency in certain scenarios.

Here, the solvers are implemented for solving problems of the form
\f[
  \min_{x} f(x) = 0 \quad \text{with} \quad f: \mathbb{R} \rightarrow \mathbb{R} \text{,}
\f]
which consist in finding the minimum of the function \f$f\f$ by iteratively updating the current iterate \f$x_k\f$ until convergence is achieved. The solvers require the function \f$f\f$ and its first derivative \f$f^{\prime}(x)\f$ to be provided by the user. Alternatively, the derivative can be approximated numerically using finite differences, depending on the problem's complexity and the user's preference.
