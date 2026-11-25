# Adam method

The `Adam` class implements adaptive moment estimation (Kingma and Ba, 2015):
\f[
  \mathbf{m}_k = \beta_1\mathbf{m}_{k-1} + (1-\beta_1)\mathbf{g}_k \text{,} \quad
  \mathbf{v}_k = \beta_2\mathbf{v}_{k-1} + (1-\beta_2)\mathbf{g}_k^2 \text{,} \quad
  \mathbf{x}_{k+1} = \mathbf{x}_k - \eta_k \frac{\hat{\mathbf{m}}_k}{\sqrt{\hat{\mathbf{v}}_k} + \varepsilon} \text{,}
\f]
with the bias-corrected moments \f$\hat{\mathbf{m}}_k = \mathbf{m}_k/(1-\beta_1^k)\f$ and \f$\hat{\mathbf{v}}_k = \mathbf{v}_k/(1-\beta_2^k)\f$. The AMSGrad variant uses the running maximum of \f$\hat{\mathbf{v}}_k\f$, and an optional decay \f$\eta_k = \eta/(1 + \gamma k)\f$ reduces the step over time.

Adam is designed for stochastic gradients. On deterministic problems it is robust to scaling but converges only linearly and to moderate accuracy, so it is mostly useful as a baseline or when the gradient is noisy.
