# Conjugate gradient methods

The `ConjugateGradient` class implements the nonlinear conjugate gradient iteration
\f[
  \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{d}_k,
  \qquad
  \mathbf{d}_{k+1} = -\mathbf{g}_{k+1} + \beta_k \mathbf{d}_k
  \text{,}
\f]
where \f$\mathbf{g}_k = \nabla f(\mathbf{x}_k)^\top\f$.

The update parameter \f$\beta_k\f$ can be selected among the standard formulas summarized in Table 1.1 of Hager and Zhang's survey:

- Hestenes-Stiefel
- Fletcher-Reeves
- Polak-Ribiere
- Polak-Ribiere+
- Conjugate Descent
- Liu-Storey
- Dai-Yuan
- Hager-Zhang
- Hager-Zhang+

The implementation uses a strong-Wolfe line search and supports a Powell-style restart safeguard whenever the conjugate direction loses reliability. The `Hager-Zhang+` variant applies the standard survey-style truncation bound controlled by the solver truncation parameter.

The initial trial step passed to the line search can also be selected. By default the solver reuses the last accepted step. An optional `Barzilai-Borwein` alpha update computes
\f[
  \alpha_{k+1}^{BB} = \frac{\boldsymbol{s}_k^T \boldsymbol{s}_k}{\boldsymbol{s}_k^T \boldsymbol{y}_k}
  	ext{,}
\f]
with \f$\boldsymbol{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k\f$ and
\f$\boldsymbol{y}_k = \mathbf{g}_{k+1} - \mathbf{g}_k\f$, and falls back to the last accepted step when the curvature estimate is not positive or becomes numerically unreliable.

> Reference: Hager W. W., Zhang H., *A Survey of Nonlinear Conjugate Gradient Methods*, Pacific Journal of Optimization, 2 (2006), pp. 35-58.