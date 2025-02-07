# Greenstadt's method

Similarly to the Broyden's methods, Greenstadt's methods are quasi-Newton methods that approximates the Jacobian matrix \f$\mathbf{Jf}_{\mathbf{x}}\f$ an update rule. The generic Greenstadt's method is defined as
\f[
  \mathbf{H}_k(\mathbf{x}_k) \mathbf{h}_k = -\mathbf{f}(\mathbf{x}_k) \text{,}
\f]
where \f$\mathbf{H}_k\f$ is the (inverse) Jacobian approximation at the \f$k\f$-th iteration.

Based on the update rule, Greenstadt's method is classified into two methods: trivially *method 1* and *method 2*. The update of the Jacobian approximation is performed as
\f[
  \mathbf{H}_{k+1}^{-1} = \mathbf{H}_k^{-1} - \displaystyle\frac{\mathbf{H}_k^{-1} \Delta\mathbf{f}_k - \Delta\mathbf{x}_k}{\mathbf{g} \Delta\mathbf{f}_k} \mathbf{g} \text{,}
\f]
where \f$\Delta\mathbf{x}_k = \mathbf{x}_k - \mathbf{x}_{k-1}\f$, and \f$\Delta\mathbf{f}_k = \mathbf{f}(\mathbf{x}_k) - \mathbf{f}(\mathbf{x}_{k-1})\f$. The quantity \f$\mathbf{g}\f$ is \f$\mathbf{g} = \mathbf{f}_k^\top\f$ for *method 1* and \f$\mathbf{g} = \mathbf{H}_{k\top}^{-1} \mathbf{H}_k^{-1} \Delta\mathbf{f}_k^{\top}\f$ for *method 2*. The choice of the method is based on the convergence history, switching between the methods to adapt to the problem's behavior.

> For more details on the Greenstadt's methods refer to the reference: Spedicato E., Greenstadt J. *On some classes of variationally derived Quasi-Newton methods for systems of nonlinear algebraic equations*, Numerische Mathematik, 29 (1978), pp. 363-380.