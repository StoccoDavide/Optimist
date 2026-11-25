# Line searches

Descent optimizers update the iterate as \f$\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{d}_k\f$, where \f$\mathbf{d}_k\f$ is a descent direction, \f$\nabla f(\mathbf{x}_k)^\top \mathbf{d}_k < 0\f$. The step \f$\alpha_k\f$ is chosen by a line search on \f$\phi(\alpha) = f(\mathbf{x}_k + \alpha\mathbf{d}_k)\f$. Four strategies are available in the `Optimist::LineSearch` namespace:

- **Armijo** backtracking enforces the sufficient decrease condition
  \f[
    \phi(\alpha) \leq \phi(0) + c_1 \alpha \phi'(0) \text{,}
  \f]
  shrinking the trial step by safeguarded quadratic interpolation. It needs function values only.
- **Goldstein** enforces \f$\phi(0) + (1-c)\alpha\phi'(0) \leq \phi(\alpha) \leq \phi(0) + c\alpha\phi'(0)\f$ with \f$0 < c < 1/2\f$, by expansion and bisection, which also rules out steps that are too short.
- **Strong Wolfe** enforces sufficient decrease together with the curvature condition \f$|\phi'(\alpha)| \leq c_2|\phi'(0)|\f$ (Nocedal and Wright, Algorithms 3.5 and 3.6), with a bracketing phase followed by a zoom phase based on safeguarded cubic interpolation. It guarantees \f$\mathbf{s}_k^\top\mathbf{y}_k > 0\f$, which keeps quasi-Newton updates positive definite, and it is the default choice of the BFGS solver.

- **Moré-Thuente** (1994) enforces the same strong Wolfe conditions with the safeguarded cubic and quadratic interpolation steps of MINPACK-2 (`dcsrch`/`dcstep`), working on the auxiliary function \f$\psi(\alpha) = \phi(\alpha) - \phi(0) - c_1\alpha\phi'(0)\f$ until a suitable bracket is found.

Close to a minimizer the decrease of \f$\phi\f$ falls below its rounding error, so that a literal sufficient decrease test would reject every step. All searches therefore accept \f$\phi(\alpha) \leq \phi(0) + c_1\alpha\phi'(0) + \epsilon_f|\phi(0)|\f$, and compare function values inside the bracket only above that level (below it, the derivative decides), in the spirit of the approximate Wolfe conditions of Hager and Zhang. The default \f$\epsilon_f = 100\,\epsilon\f$ can be raised for cost functions evaluated with heavy cancellation.

The parameters (\f$c_1 = 10^{-4}\f$, \f$c_2 = 0.9\f$, \f$c = 0.25\f$, \f$\epsilon_f\f$, maximum step and number of trials) are collected in `LineSearch::Settings`.
