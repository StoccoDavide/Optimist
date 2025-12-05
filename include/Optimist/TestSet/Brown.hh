/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_BROWN_HH
#define OPTIMIST_TESTSET_BROWN_HH

#include "Optimist/Function.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
     * \brief Class container for the Brown badly scaled function.
     *
     * Class container for the Brown badly scaled function, which is defined as:
     * \f[
     * \mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 - a \\ x_2 - 2a \\ x_1x_2 - 2 \end{bmatrix} \text{,}
     * \f]
     * where \f$a = 10^{-6}\f$. The function has one solution at \f$\mathbf{x} = [a, 2a]^\top\f$,
     * with \f$f(\mathbf{x}) = 0\f$. The initial guess is generated at \f$\mathbf{x} = [1, 1]^\top\f$.
     * \tparam Input Input vector type.
     * \tparam Output Output vector type.
     */
    template <typename Input, typename Output>
    requires TypeTrait<Input>::IsEigen && TypeTrait<Output>::IsEigen &&
      (!TypeTrait<Input>::IsFixed || TypeTrait<Input>::Dimension == 2) &&
      (!TypeTrait<Output>::IsFixed || TypeTrait<Output>::Dimension == 3)
    class Brown : public Function<Input, Output, Brown<Input, Output>>
    {
    private:

      using VectorTraitInput = TypeTrait<Input>;
      using VectorTraitOutput = TypeTrait<Output>;
      using Scalar = typename Input::Scalar;
      using typename Function<Input, Output, Brown<Input, Output>>::Input;
      using typename Function<Input, Output, Brown<Input, Output>>::Output;
      using typename Function<Input, Output, Brown<Input, Output>>::Matrix;
      using typename Function<Input, Output, Brown<Input, Output>>::Tensor;

      Scalar m_a{1.0e-6}; /**< Scaling value (keep it low to guarantee bad scaling). */

    public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Brown function.
       */
      Brown()
      {
        this->m_solutions.emplace_back(this->m_a, 2.0*this->m_a);
        this->m_guesses.emplace_back(1.0, 1.0);
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "Brown";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(const Input & x, Output & out) const
      {
        out << x(0) - this->m_a,
               x(1) - 2.0*this->m_a,
               x(0)*x(1) - 2.0;
        return out.allFinite();
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(const Input & x, Matrix & out) const
      {
        out << 1.0, 0.0, x(1),
               0.0, 1.0, x(0);
        return out.allFinite();
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(const Input & /*x*/, Tensor & out) const
      {
        out.resize(this->output_dimension());
        out[0].setZero();
        out[1].setZero();
        out[2] << 0.0, 1.0,
                  1.0, 0.0;
        return out[0].allFinite() && out[1].allFinite() && out[2].allFinite();
      }

    }; // class Brown

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_BROWN_HH
