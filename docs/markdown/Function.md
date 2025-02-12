# Function

The `Function` class in `Optimist` provides a common interface for representing and evaluating mathematical functions. This class is designed to handle a generic function \f$\mathbf{f}(\mathbf{x}): \mathbb{R}^{n} \rightarrow \mathbb{R}^{m}\f$, where \f$\mathbf{x}\f$ is a vector of input variables and \f$\mathbf{f}(\mathbf{x})\f$ is a vector of output values.

The `Function` class serves as a base class for developing specific function representations that can be used in various numerical algorithms. In general, the `Function` class and its derived classes are designed to support:
- **Root-finding** problems of the form
\f[
  \mathbf{f}(\mathbf{x}) = \mathbf{0} \text{,} \quad \text{with} \quad \mathbf{f}: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n} \text{,}
\f]
- **Optimization** problems of the form
\f[
  \min_{\mathbf{x}} f(\mathbf{x}) \text{,} \quad f: \mathbb{R}^{n} \rightarrow \mathbb{R} \text{.}
\f]

If \f$n > 1\f$, the function is multi-dimensional, otherwise, it is scalar. The linear algebra operations required for multi-dimensional functions are performed using the `Eigen` library, which provides a high-level interface for matrix operations. For scalar functions (\f$n = 1\f$), the base type changes to `Real`, a scalar type defined in the `Optimist` library.

