# Broyden's method

The Broyden family includes the good, bad, and combined methods. These methods are quasi-Newton approaches that approximate the Jacobian matrix or its inverse. The combined method integrates the updates of the good and bad methods, offering an alternative approach to the classical Newton's method. The Broyden's methods are based on the linearization of the function \f$\mathbf{f}(\mathbf{x})\f$ around the current iterate \f$\mathbf{x}_k\f$, which leads to the linear system
\f[
  \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}_k) \mathbf{h}_k = -\mathbf{f}(\mathbf{x}_k) \text{.}
\f]
Here, the Jacobian matrix \f$\mathbf{Jf}_{\mathbf{x}}\f$ is supposed not to be explicitly available or computationally expensive to compute. For this reason, Broyden's methods approximate the Jacobian matrix by updating the inverse of the Jacobian matrix iteratively. The Broyden's method is defined as
\f[
  \mathbf{H}_k(\mathbf{x}_k) \mathbf{h}_k = -\mathbf{f}(\mathbf{x}_k) \text{,}
\f]
where \f$\mathbf{H}_k\f$ is the (inverse) Jacobian approximation at the \f$k\f$-th iteration. The generic advancing step is then defined as
\f[
  \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{h}_k \text{.}
\f]
where \f$\alpha_k\f$ is a damping coefficient that ensures affine-invariant criteria is satisfied.

What distinguishes Broyden's *combined* method from the generic Broyden's method is the update of the Jacobian approximation. The Broyden's *combined* method uses the following update rule
\f[
  \begin{cases}
    \mathbf{H}_{k+1}^{-1} = \mathbf{H}_k^{-1} - \displaystyle\frac{\mathbf{H}_k^{-1} \Delta\mathbf{f}_k - \Delta\mathbf{x}_k}{\mathbf{C}_g \Delta\mathbf{f}_k} \mathbf{C}_g  &
    \left\| \displaystyle\frac{\Delta\mathbf{x}_k^\top \Delta\mathbf{x}_{k-1}}{\Delta\mathbf{x}_k^\top \mathbf{H}_k^{-1} \Delta\mathbf{f}_k} \right\| < \left\| \displaystyle\frac{\Delta\mathbf{f}_k^\top  \Delta\mathbf{f}_{k-1}}{\Delta\mathbf{f}_k^\top\Delta\mathbf{f}_k} \right\| \\
    \mathbf{H}_{k+1}^{-1} = \mathbf{H}_k^{-1} - \displaystyle\frac{\mathbf{H}_k^{-1} \Delta\mathbf{f}_k-\Delta\mathbf{x}_k}{\mathbf{C}_b \Delta\mathbf{f}_k} \mathbf{C}_b &
    \text{otherwise}
  \end{cases} \text{,}
\f]
with \f$\mathbf{C}_g = \Delta\mathbf{x}_k^\top \mathbf{H}_k^{-1}\f$, \f$\mathbf{C}_b = \Delta\mathbf{f}_k^\top\f$, \f$\Delta\mathbf{x}_k = \mathbf{x}_{k+1} - \mathbf{x}_k\f$, and \f$\Delta\mathbf{f}_k = \mathbf{f}(\mathbf{x}_{k+1}) - \mathbf{f}(\mathbf{x}_k)\f$. Such a rule allows Broyden's *combined* method to adapt to the problem's behavior, switching between the good and bad methods based on the convergence history.

> For more details on the Broyden's methods refer to the references: Broyden C.G., *A class of methods for solving nonlinear simultaneous equations*, Mathematics of Computation, 19 (1965), pp. 577-593, 10.1090/s0025-5718-1965-0198670-6, and Martinez J., Ochi L., *Sobre dois métodos de broyden*, Matemática Aplicada e Computacional, IV Congresso Nacional de Matemática Aplicada e Computacional, Rio de Janeiro, Brasil, setembro de 1981.

> **NOTE:** In `Optimist`, the implemented Broyden's class for the *combined* method can be used as a *good* or *bad* solver by setting the appropriate parameters.
