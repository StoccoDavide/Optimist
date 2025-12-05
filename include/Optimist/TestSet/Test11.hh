/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef OPTIMIST_TESTSET_TEST11_HH
#define OPTIMIST_TESTSET_TEST11_HH

#include "Optimist/Function.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
     * \brief Class container for the Test11 function.
     *
     * Class container for the Test11 function, which is defined as:
     * \f[
     * \mathbf{f}(\mathbf{x}) = \begin{bmatrix} x^2 \end{bmatrix} \text{.}
     * \f]
     * The function has one solution at \f$\mathbf{x} = [0]\f$, with \f$f(\mathbf{x}) = [0]\f$.
     * The initial guesses are generated on the square \f$x_i \in [-10, 10]\f$.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    requires TypeTrait<Vector>::IsEigen && (TypeTrait<Vector>::IsDynamicSize ||
      (TypeTrait<Vector>::IsFixedSize && Vector::RowsAtCompileTime == 1))
    class Test11 : public Function<Vector, Vector, Test11<Vector>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar = typename Vector::Scalar;
      using typename Function<Vector, Vector, Test11<Vector>>::Input;
      using typename Function<Vector, Vector, Test11<Vector>>::Output;
      using typename Function<Vector, Vector, Test11<Vector>>::Matrix;
      using typename Function<Vector, Vector, Test11<Vector>>::Tensor;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Test11 function.
       */
      Test11()
      {
        this->m_solutions.emplace_back(-10.0);
        this->m_solutions.emplace_back(-5.0);
        this->m_solutions.emplace_back(0.0);
        this->m_solutions.emplace_back(5.0);
        this->m_solutions.emplace_back(10.0);
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "Test11";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Input const & x, Output & out) const
      {
        out << x(0) * x(0);
        return out.allFinite();
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Input const & x, Matrix & out) const
      {
        out << 2.0 * x(0);
        return out.allFinite();
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Input const & /*x*/, Tensor & out) const
      {
        out << 2.0;
        return out.allFinite();
      }

    }; // class Test11

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_TEST11_HH
