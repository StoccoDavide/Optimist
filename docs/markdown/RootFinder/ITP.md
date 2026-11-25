# ITP method

The `ITP` (Interpolate, Truncate, Project) class implements the bracketing method of Oliveira and Takahashi (2020). At each iteration on the bracket \f$[a, b]\f$:

1. **Interpolate**: compute the regula falsi point \f$x_f\f$;
2. **Truncate**: move it towards the midpoint \f$x_{1/2}\f$ by \f$\delta = \kappa_1(b-a)^{\kappa_2}\f$;
3. **Project**: project it onto the interval of radius \f$r_j = \epsilon\,2^{n_{\max}-j} - (b-a)/2\f$ around the midpoint.

The projection guarantees that the method never needs more than \f$n_{\max} = \lceil\log_2((b-a)/2\epsilon)\rceil + n_0\f$ iterations, i.e., it preserves the worst case of bisection, while the interpolation gives superlinear convergence (order up to \f$\kappa_2\f$) on smooth functions. The default parameters are \f$\kappa_1 = 0.2/(b-a)\f$, \f$\kappa_2 = 2\f$ and \f$n_0 = 1\f$.
