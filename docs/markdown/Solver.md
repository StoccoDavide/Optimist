# Solver

The solvers implemented in `Optimist` are derived from the `Solver` class, which provides a common interface for solving nonlinear systems of equations. The `Solver` class is an abstract base class that defines the methods and attributes required by the derived classes. The derived classes implement the specific algorithms for solving the nonlinear systems of equations.

This class provides a base class for developing solvers that handle a generic function \f$\mathbf{f}(\mathbf{x}): \mathbb{R}^{n} \rightarrow \mathbb{R}^{m}\f$. In general, the derived classes are designed to solve root-finding problems of the form
\f[
  \mathbf{f}(\mathbf{x}) = \mathbf{0} \text{,} \quad \text{with} \quad \mathbf{f}: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n} \text{,}
\f]
as well as optimization problems of the form
\f[
  \min_{\mathbf{x}} f(\mathbf{x}) \text{,} \quad f: \mathbb{R}^{n} \rightarrow \mathbb{R} \text{.}
\f]
If \f$ n > 1 \f$ the solver is a multi-dimensional solver, otherwise, it is a scalar solver. The linear algebra operations are performed using the `Eigen` library, which provides a high-level interface for matrix operations. Otherwise, for \f$ n = 1 \f$, the base type changes to `Real`, a scalar type defined in the `Optimist` library.

