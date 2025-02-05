# Root finder

This section describes the multi-dimensional root-finders implemented in `Optimist`. The available root-finding solvers are *Newton* (derivative) and *quasi-Newton* (non-derivative) methods. Derivative methods employ the function's derivative to find the root with high accuracy, while non-derivative methods approximate the derivative for improved efficiency in certain scenarios.

Here, the solvers are implemented for solving problems of the form
\f[
  \mathbf{f}(\mathbf{x}) = \mathbf{0} \text{,} \quad \text{with} \quad \mathbf{f}: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n} \text{,}
\f]
which consist in finding the root of the function \f$\mathbf{f}\f$ by iteratively updating the current iterate \f$\mathbf{x}_k\f$ until convergence is achieved. The solvers require the function \f$\mathbf{f}\f$ and its Jacobian matrix \f$\mathbf{Jf}_{\mathbf{x}}\f$ to be provided by the user. The Jacobian matrix can be computed analytically or numerically, depending on the problem's complexity and the user's preference. Alternatively, the Jacobian can be approximated numerically using finite differences, depending on the problem's complexity and the user's preference.

## Affine-invariant step

By default, the multi-dimensional root-finders in `Optimist` use optionally an affine-invariant advancing step to ensure convergence. The generic advancing step \f$\mathbf{h}_k\f$ is computed as
\f[
  \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{h}_k \text{,}
\f]
where \f$\alpha_k\f$ is a damping coefficient that ensures convergence by satisfying
\f[
  \left\|\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}_k)^{-1} \mathbf{f}(\mathbf{x}_{k+1})\right\|
  \leq \left(1 - \displaystyle\frac{\alpha_k}{2}\right) \left\|\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}_k)^{-1}
  \mathbf{f}(\mathbf{x}_k)\right\| = \left(1 - \displaystyle\frac{\alpha_k}{2} \right)
  \left\|\mathbf{h}_k\right\| \text{.}
\f]
For detailed information on the affine-invariant Newton's method, refer to [this link](https://www.zib.de/deuflhard/research/algorithm/ainewton.en.html), or the works by P. Deuflhard.
