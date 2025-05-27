# Newton's method

Newton's method, including its damped variant with an affine-invariant step, is based on the linearization of the function \f$f(x)\f$ around the current iterate \f$x_k\f$, which leads to the linear system
\f[
  f^{\prime}_{x}(x_k) h_k = -f(x_k) \text{.}
\f]
The advancing step \f$h_k\f$ is then computed as
\f[
  x_{k+1} = x_k + \alpha_k h_k \text{,}
\f]
where \f$\alpha_k\f$ is a damping coefficient that ensures affine-invariant criteria is satisfied.
