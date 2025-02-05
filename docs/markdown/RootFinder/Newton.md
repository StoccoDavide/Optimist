# Newton's method

Newton's method, including its damped variant with an affine-invariant step, is based on the linearization of the function \f$\mathbf{f}(\mathbf{x})\f$ around the current iterate \f$\mathbf{x}_k\f$, which leads to the linear system
\f[
  \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}_k) \mathbf{h}_k = -\mathbf{f}(\mathbf{x}_k) \text{.}
\f]
The advancing step \f$\mathbf{h}_k\f$ is then computed as
\f[
  \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{h}_k \text{,}
\f]
where \f$\alpha_k\f$ is a damping coefficient that ensures affine-invariant criteria is satisfied.
