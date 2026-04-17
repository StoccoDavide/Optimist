/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_ELLIPTIC_PARABOLOID_HH
#define OPTIMIST_TESTSET_ELLIPTIC_PARABOLOID_HH

#include "Optimist/Function.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the paraboloid function.
     *
     * Class container for the paraboloid function, which is defined as:
     * \f[
     * f(\mathbf{x}) = ax^2 + by^2 \text{.}
     * \f]
     * The function has global minima at \f$\mathbf{x} = (0, 0)\f$, with
     * \f$f(\mathbf{x}) = 0\f$. The initial guesses are generated on the square
     * \f$x_i \in \left[-100, 100\right]\f$.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen && (!TypeTrait<Vector>::IsFixed ||
                                              TypeTrait<Vector>::Dimension == 2)
    class EllipticParaboloid : public Function<Vector,
                                               typename Vector::Scalar,
                                               EllipticParaboloid<Vector>> {
     public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector, Scalar, EllipticParaboloid<Vector>>::
          FirstDerivative;
      using typename Function<Vector, Scalar, EllipticParaboloid<Vector>>::
          SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_a{1.0}; /**< Coefficient \f$ a \f$. */
      Scalar m_b{1.0}; /**< Coefficient \f$ b \f$. */

     public:
      /**
       * Class constructor for the paraboloid function.
       */
      EllipticParaboloid() {
        this->m_solutions.resize(1);
        if constexpr (VectorTrait::IsFixed) {
          this->m_solutions[0] << 0.0, 0.0;
          for (Scalar x{-10}; x < 10 + EPSILON; x += 10 / 2.5) {
            for (Scalar y{-10}; y < 10 + EPSILON; y += 10 / 2.5) {
              this->m_guesses.emplace_back(x, y);
            }
          }
        } else if constexpr (VectorTrait::IsDynamic) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0] << 0.0, 0.0;
          for (Scalar x{-10}; x < 10 + EPSILON; x += 10 / 2.5) {
            for (Scalar y{-10}; y < 10 + EPSILON; y += 10 / 2.5) {
              Vector guess(2);
              guess << x, y;
              this->m_guesses.push_back(guess);
            }
          }
        } else if constexpr (VectorTrait::IsSparse) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0].reserve(2);
          this->m_solutions[0].coeffRef(0) = 0.0;
          this->m_solutions[0].coeffRef(1) = 0.0;
          for (Scalar x{-10}; x < 10 + EPSILON; x += 10 / 2.5) {
            for (Scalar y{-10}; y < 10 + EPSILON; y += 10 / 2.5) {
              Vector guess(2);
              guess.reserve(2);
              guess.coeffRef(0) = x;
              guess.coeffRef(1) = y;
              this->m_guesses.push_back(guess);
            }
          }
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "EllipticParaboloid";
      }

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       */
      bool evaluate_impl(const Vector &x, Scalar &out) const {
        Scalar x_0, x_1;
        if constexpr (VectorTrait::IsSparse) {
          x_0 = x.coeff(0);
          x_1 = x.coeff(1);
        } else {
          x_0 = x(0);
          x_1 = x(1);
        }
        out = this->m_a * x_0 * x_0 + this->m_b * x_1 * x_1;
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       */
      bool first_derivative_impl(const Vector &x, FirstDerivative &out) const {
        Scalar x_0, x_1;
        if constexpr (VectorTrait::IsSparse) {
          x_0 = x.coeff(0);
          x_1 = x.coeff(1);
        } else {
          x_0 = x(0);
          x_1 = x(1);
        }

        if constexpr (VectorTrait::IsFixed) {
          out << 2.0 * this->m_a * x_0, 2.0 * this->m_b * x_1;
          return out.allFinite();
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2);
          out << 2.0 * this->m_a * x_0, 2.0 * this->m_b * x_1;
          return out.allFinite();
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2);
          out.reserve(2);
          out.coeffRef(0) = 2.0 * this->m_a * x_0;
          out.coeffRef(1) = 2.0 * this->m_b * x_1;
          for (typename FirstDerivative::InnerIterator it(out); it; ++it) {
            if (!std::isfinite(it.value())) {
              return false;
            }
          }
          return true;
        }
        return false;
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       */
      bool second_derivative_impl(const Vector & /*x*/,
                                  SecondDerivative &out) const {
        if constexpr (VectorTrait::IsFixed) {
          out << 2.0 * this->m_a, 0.0, 0.0, 2.0 * this->m_b;
          return out.allFinite();
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2, 2);
          out << 2.0 * this->m_a, 0.0, 0.0, 2.0 * this->m_b;
          return out.allFinite();
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2, 2);
          out.reserve(2);
          out.coeffRef(0, 0) = 2.0 * this->m_a;
          out.coeffRef(1, 1) = 2.0 * this->m_b;
          for (Integer k{0}; k < out.outerSize(); ++k) {
            for (typename SecondDerivative::InnerIterator it(out, k); it;
                 ++it) {
              if (!std::isfinite(it.value())) {
                return false;
              }
            }
          }
          return true;
        }
        return false;
      }

    };  // class EllipticParaboloid

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_ELLIPTIC_PARABOLOID_HH
