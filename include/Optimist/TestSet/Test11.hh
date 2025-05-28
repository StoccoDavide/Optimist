/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef OPTIMIST_TESTSET_TEST11_HH
#define OPTIMIST_TESTSET_TEST11_HH

#include "Optimist/TestSet.hh"

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
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Test11 : public Function<Real, 1, 1, Test11<Real>, true>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      using Vector = typename Function<Real, 1, 1, Test11<Real>, true>::InputVector;
      using Matrix = typename Function<Real, 1, 1, Test11<Real>, true>::Matrix;
      using Tensor = typename Function<Real, 1, 1, Test11<Real>, true>::Tensor;

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
      std::string name_impl() const {return "Test11";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(const Vector & x, Vector & out) const
      {
        out << x(0) * x(0);
      }
      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(const Vector & x, Matrix & out) const
      {
        out << 2.0 * x(0);
      }

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(const Vector & /*x*/, Tensor & out) const
      {
        out << 2.0;
      }

    }; // class Test11

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_TEST11_HH
