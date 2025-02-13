/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_SCHAFFER2_HXX
#define OPTIMIST_SCHAFFER2_HXX

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the Schaffer2 function.
    *
    * Class container for the Schaffer2 function, which is defined as:
    * \f[
    * f(\mathbf{x}) = 0.5 + \displaystyle\frac{\sin^{2}(x_1^2 - x_2^2) - 0.5}{(1 + 0.001(x_1^2 + x_2^2))^2} \text{.}
    * \f]
    * The function has global minima at \f$\mathbf{x} = (0, 0)\f$, with \f$f(\mathbf{x}) = 0.0\f$.
    * The initial guesses are generated on the square \f$x_i \in \left[-100, 100\right]\f$.
    */
    class Schaffer2 : public CostFunction<2, Schaffer2>
    {
    public:
      using Vector    = typename CostFunction<2, Schaffer2>::Vector;
      using RowVector = typename CostFunction<2, Schaffer2>::RowVector;
      using Matrix    = typename CostFunction<2, Schaffer2>::Matrix;

      /**
      * Class constructor for the Schaffer2 function.
      */
      Schaffer2()
      {
        this->m_solutions.emplace_back(0.0, 0.0);
        for (Real x{-100}; x < 100 + EPSILON; x += 100/25.0) {
          for (Real y{-100}; y < 100 + EPSILON; y += 100/25.0) {
            this->m_guesses.emplace_back(x, y);
          }
        }
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Schaffer2";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(const Vector & x, Real & out) const
      {
        Real xx_0{x(0)*x(0)};
        Real xx_1{x(1)*x(1)};
        Real xx_0_m_xx_1{xx_0 - xx_1};
        Real xx_0_p_xx_1{xx_0 + xx_1};
        out = 0.5 + (std::sin(xx_0_m_xx_1)*std::sin(xx_0_m_xx_1) - 0.5) /
          ((1.0 + 0.001*(xx_0_p_xx_1))*(1.0 + 0.001*(xx_0_p_xx_1)));
      }

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(const Vector & x, RowVector & out) const
      {
        Real xx_0{x(0)*x(0)};
        Real xx_1{x(1)*x(1)};
        Real xx_0_m_xx_1{xx_0 - xx_1};
        Real xx_0_p_xx_1{xx_0 + xx_1};
        Real tmp{1.0 + 0.001*(xx_0_p_xx_1)};
        Real tmp2{tmp*tmp}, tmp3{tmp2*tmp};
        out(0) = 2.0*x(0)*std::sin(xx_0_m_xx_1) / tmp2 -
          2.0*0.001*x(0)*std::cos(xx_0_m_xx_1) / tmp3;
        out(1) = -2.0*x(1)*std::sin(xx_0_m_xx_1) / tmp2 +
          2.0*0.001*x(1)*std::cos(xx_0_m_xx_1) / tmp3;
      }

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(const Vector & x, Matrix & out) const
      {
        Real xx_0{x(0)*x(0)};
        Real xx_1{x(1)*x(1)};
        Real xx_0_m_xx_1{xx_0 - xx_1};
        Real xx_0_p_xx_1{xx_0 + xx_1};
        Real tmp{1.0 + 0.001*(xx_0_p_xx_1)};
        Real tmp2{tmp*tmp}, tmp3{tmp2*tmp}, tmp4{tmp2*tmp2};
        out(0, 0) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
          4.0*x(0)*x(0)*std::sin(xx_0_m_xx_1) / tmp3 +
          6.0*0.001*x(0)*x(0)*std::cos(xx_0_m_xx_1) / tmp4;
        out(0, 1) = -2.0*x(0)*x(1)*std::sin(xx_0_m_xx_1) / tmp3 +
          6.0*0.001*x(0)*x(1)*std::cos(xx_0_m_xx_1) / tmp4;
        out(1, 0) = out(0, 1);
        out(1, 1) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
          4.0*x(1)*x(1)*std::sin(xx_0_m_xx_1) / tmp3 +
          6.0*0.001*x(1)*x(1)*std::cos(xx_0_m_xx_1) / tmp4;
      }

    }; // class Schaffer2

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_COS_HXX
