# Trust-region Newton method

The `TrustRegion` class minimizes the quadratic model
\f[
  m_k(\mathbf{p}) = f(\mathbf{x}_k) + \nabla f(\mathbf{x}_k)^\top\mathbf{p} + \frac{1}{2}\mathbf{p}^\top\nabla^2 f(\mathbf{x}_k)\mathbf{p}
  \quad \text{subject to} \quad \|\mathbf{p}\| \leq \Delta_k \text{,}
\f]
and accepts the step if the ratio between the actual and the predicted reduction, \f$\rho_k\f$, is positive enough. The radius is reduced when \f$\rho_k < 1/4\f$ and enlarged when \f$\rho_k > 3/4\f$ and the step reached the boundary. Unlike line-search Newton methods, the iteration is well defined when the Hessian is indefinite.

Two subproblem solvers are available:

- **Steihaug-Toint** truncated conjugate gradient (default), which stops at the boundary of the region or along a direction of negative curvature, and only needs Hessian-vector products;
- **Powell's dogleg**, which moves along the piecewise-linear path from the Cauchy point to the Newton point, and falls back to Steihaug-Toint when the Hessian is not positive definite.

For vector-valued functions the solver minimizes \f$\frac{1}{2}\|\mathbf{f}\|^2\f$ with the full Hessian \f$\mathbf{J}^\top\mathbf{J} + \sum_i f_i\nabla^2 f_i\f$.
