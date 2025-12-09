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

#ifndef OPTIMIST_TESTSET_ELLIPTICPARABOLOID_HH
#define OPTIMIST_TESTSET_ELLIPTICPARABOLOID_HH

#include "Optimist/Function.hh"

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
     * \tparam Scalar Floating-point number type.
     */
    template <typename Vector>
    requires TypeTrait<Vector>::IsEigen &&
      (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == 2)
    class EllipticParaboloid : public Function<Vector, typename Vector::Scalar, EllipticParaboloid<Vector>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector, Scalar, EllipticParaboloid<Scalar>>::Vector;
      using typename Function<Vector, Scalar, EllipticParaboloid<Scalar>>::RowVector;
      using typename Function<Vector, Scalar, EllipticParaboloid<Scalar>>::Matrix;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

    private:
      Scalar m_a{1.0}; /**< Coefficient \f$ a \f$. */
      Scalar m_b{1.0}; /**< Coefficient \f$ b \f$. */

    public:
      /**
       * Class constructor for the paraboloid function.
       */
      EllipticParaboloid()
      {
        this->m_solutions.emplace_back(0.0, 0.0);
        for (Scalar x{-100}; x < 100 + EPSILON; x += 100/25.0) {
          for (Scalar y{-100}; y < 100 + EPSILON; y += 100/25.0) {
            this->m_guesses.emplace_back(x, y);
          }
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "EllipticParaboloid";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       */
      bool evaluate_impl(Vector const & x, Scalar & out) const
      {
        out = this->m_a*x(0)*x(0) + this->m_b*x(1)*x(1);
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       */
      bool first_derivative_impl(Vector const & x, RowVector & out) const
      {
        out << 2.0*this->m_a*x(0), 2.0*this->m_b*x(1);
        return out.allFinite();
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       */
      bool second_derivative_impl(Vector const & /*x*/, Matrix & out) const
      {
        out << 2.0*this->m_a, 0.0,
               0.0, 2.0*this->m_b;
        return out.allFinite();
      }

    }; // class EllipticParaboloid

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_ELLIPTICPARABOLOID_HH
