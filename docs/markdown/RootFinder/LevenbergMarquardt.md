# Levenberg-Marquardt method

The `LevenbergMarquardt` class solves \f$\mathbf{f}(\mathbf{x}) = \mathbf{0}\f$ by minimizing \f$\frac{1}{2}\|\mathbf{f}\|^2\f$ with the damped Gauss-Newton step
\f[
  (\mathbf{J}^\top\mathbf{J} + \lambda\mathbf{D})\,\mathbf{h} = -\mathbf{J}^\top\mathbf{f} \text{,}
\f]
where \f$\mathbf{D} = \mathrm{diag}(\mathbf{J}^\top\mathbf{J})\f$ (Marquardt's scaling, default) or \f$\mathbf{D} = \mathbf{I}\f$. The damping interpolates between Gauss-Newton (\f$\lambda \to 0\f$) and scaled steepest descent (large \f$\lambda\f$), and is updated with Nielsen's rule from the ratio between actual and predicted reduction: \f$\lambda \leftarrow \lambda\max(1/3, 1-(2\rho-1)^3)\f$ on success, \f$\lambda \leftarrow \nu\lambda\f$, \f$\nu \leftarrow 2\nu\f$ on failure.

The damped normal equations remain solvable when the Jacobian is singular, which makes the method more robust than Newton's method far from the solution and at singular roots (where the convergence becomes linear). The solver also accepts overdetermined systems of \f$m > n\f$ equations, `LevenbergMarquardt<Input, Output>`: the iteration then stops either at a root or at a least-squares solution, where \f$\|\mathbf{J}^\top\mathbf{f}\| < \epsilon\max(1, \|\mathbf{f}\|)\f$.
