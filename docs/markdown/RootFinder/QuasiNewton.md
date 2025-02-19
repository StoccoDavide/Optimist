# Quasi-Newton method

The quasi-Newton is a family of methods that approximate the Jacobian matrix or its inverse. Such methods are based on the linearization of the function \f$\mathbf{f}(\mathbf{x})\f$ around the current iterate \f$\mathbf{x}_k\f$, which leads to the linear system
\f[
  \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}_k) \mathbf{h}_k = -\mathbf{f}(\mathbf{x}_k) \text{.}
\f]
Here, the Jacobian matrix \f$\mathbf{Jf}_{\mathbf{x}}\f$ is supposed not to be explicitly available or computationally expensive to compute. For this reason, quasi-Newton methods approximate it iteratively.
