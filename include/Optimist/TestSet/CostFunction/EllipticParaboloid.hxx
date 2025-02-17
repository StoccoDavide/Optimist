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

#ifndef OPTIMIST_COST_FUNCTION_ELLIPTIC_PARABOLOID_HXX
#define OPTIMIST_COST_FUNCTION_ELLIPTIC_PARABOLOID_HXX

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the paraboloid function.
    *
    * Class container for the paraboloid function, which is defined as:
    * \f[
    * f(\mathbf{x}) = ax^2 + by^2 \text{.}
    * \f]
    * The function has global minima at \f$\mathbf{x} = (0, 0)\f$, with \f$f(\mathbf{x}) = 0\f$.
    * The initial guesses are generated on the square \f$x_i \in \left[-100, 100\right]\f$.
    */
    class EllipticParaboloid : public CostFunction<2, EllipticParaboloid>
    {
      Real m_a{1.0}; /**< Coefficient \f$ a \f$. */
      Real m_b{1.0}; /**< Coefficient \f$ b \f$. */

    public:
      using Vector    = typename CostFunction<2, EllipticParaboloid>::Vector;
      using RowVector = typename CostFunction<2, EllipticParaboloid>::RowVector;
      using Matrix    = typename CostFunction<2, EllipticParaboloid>::Matrix;

      /**
      * Class constructor for the paraboloid function.
      */
      EllipticParaboloid()
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
      std::string name_impl() const {return "EllipticParaboloid";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(const Vector & x, Real & out) const
      {
        out = this->m_a*x(0)*x(0) + this->m_b*x(1)*x(1);
      }

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(const Vector & x, RowVector & out) const
      {
        out(0) = 2.0*this->m_a*x(0);
        out(1) = 2.0*this->m_b*x(1);
      }

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(const Vector & /*x*/, Matrix & out) const
      {
        out(0, 0) = 2.0*this->m_a;
        out(0, 1) = 0.0;
        out(1, 0) = 0.0;
        out(1, 1) = 2.0*this->m_b;
      }

    }; // class EllipticParaboloid

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_COST_FUNCTION_ELLIPTIC_PARABOLOID_HXX
