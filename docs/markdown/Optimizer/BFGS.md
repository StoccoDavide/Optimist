# BFGS and L-BFGS methods

The `BFGS` class implements the Broyden-Fletcher-Goldfarb-Shanno quasi-Newton method, where the search direction \f$\mathbf{d}_k = -\mathbf{H}_k \nabla f(\mathbf{x}_k)\f$ uses an approximation \f$\mathbf{H}_k\f$ of the inverse Hessian, updated with the step \f$\mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k\f$ and the gradient change \f$\mathbf{y}_k = \nabla f(\mathbf{x}_{k+1}) - \nabla f(\mathbf{x}_k)\f$ as
\f[
  \mathbf{H}_{k+1} = (\mathbf{I} - \rho_k \mathbf{s}_k \mathbf{y}_k^\top) \mathbf{H}_k (\mathbf{I} - \rho_k \mathbf{y}_k \mathbf{s}_k^\top) + \rho_k \mathbf{s}_k \mathbf{s}_k^\top \text{,}
  \qquad \rho_k = \frac{1}{\mathbf{y}_k^\top \mathbf{s}_k} \text{.}
\f]
The initial approximation is rescaled after the first step by \f$\mathbf{s}^\top\mathbf{y}/\mathbf{y}^\top\mathbf{y}\f$, and the update is skipped whenever the curvature condition \f$\mathbf{s}_k^\top\mathbf{y}_k > 0\f$ fails.

The limited-memory variant (`Method::LBFGS`) stores only the last \f$m\f$ pairs \f$(\mathbf{s}_i, \mathbf{y}_i)\f$ and computes \f$\mathbf{H}_k\nabla f\f$ with the two-loop recursion, so that memory and cost per iteration are \f$O(mn)\f$.

The step is selected by one of the line searches of `Optimist::LineSearch` (strong Wolfe by default), whose sufficient decrease condition tolerates changes of the cost below its rounding error, so that the iteration can proceed until \f$\|\nabla f\| < \epsilon\f$ even when \f$|f|\f$ is large at the minimizer.
