# Scalar root finder

This section describes the scalar root-finders implemented in `Optimist`. The available root-finders are derivative and non-derivative methods. Derivative methods employ the function's derivative to find the root with high accuracy, while non-derivative methods approximate the derivative for improved efficiency in certain scenarios.

Here, the solvers are implemented for solving problems of the form
\f[
  f(x) = 0 \text{,} \quad \text{with} \quad f: \mathbb{R} \rightarrow \mathbb{R} \text{,}
\f]
which consist in finding the root of the function \f$f\f$ by iteratively updating the current iterate \f$x_k\f$ until convergence is achieved. The solvers require the function \f$f\f$ and its first derivative \f$f^{\prime}(x)\f$ to be provided by the user. Alternatively, the derivative can be approximated numerically using finite differences, depending on the problem's complexity and the user's preference.

## Affine-invariant step

By default, the derivative scalar root-finders in `Optimist` use optionally an affine-invariant advancing step to ensure convergence. The generic advancing step \f$h_k\f$ is computed as
\f[
  x_{k+1} = x_k + \alpha_k h_k \text{,}
\f]
where \f$\alpha_k\f$ is a damping coefficient that ensures convergence by satisfying
\f[
  \left\|f^{\prime}_{x}(x_k)^{-1} f(x_{k+1})\right\|
  \leq \left(1 - \displaystyle\frac{\alpha_k}{2}\right) \left\|f^{\prime}_{x}(x_k)^{-1}
  f(x_k)\right\| = \left(1 - \displaystyle\frac{\alpha_k}{2} \right)
  \left\|h_k\right\| \text{.}
\f]
For detailed information on the affine-invariant Newton's method, refer to [this link](https://www.zib.de/deuflhard/research/algorithm/ainewton.en.html), or the works by P. Deuflhard.
