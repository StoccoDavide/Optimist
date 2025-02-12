# Optimizer

This section describes the scalar optimizers implemented in `Optimist`. The available optimizers are derivative and non-derivative methods. Derivative methods employ the function's derivative to find the minimum with high accuracy, while non-derivative methods approximate the derivative for improved efficiency in certain scenarios.

Here, the solvers are implemented for solving problems of the form
\f[
  \min_{\mathbf{x}} \mathbf{f}(\mathbf{x}) = 0 \quad \text{with} \quad \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R} \text{,}
\f]
which consist in finding the minimum of the cost function \f$\mathbf{f}\f$ by iteratively updating the current iterate \f$\mathbf{x}_k\f$ until convergence is achieved. The solvers require the cost function \f$\mathbf{f}\f$ and its first derivative \f$\mathbf{f}^{\prime}(\mathbf{x})\f$ to be provided by the user. Alternatively, the derivative can be approximated numerically using finite differences, depending on the problem's complexity and the user's preference.
