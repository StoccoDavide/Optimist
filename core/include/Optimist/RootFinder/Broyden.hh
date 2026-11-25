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

#ifndef OPTIMIST_ROOTFINDER_BROYDEN_HH
#define OPTIMIST_ROOTFINDER_BROYDEN_HH

#include "Optimist/RootFinder/QuasiNewton.hh"

namespace Optimist {
  namespace RootFinder {

    /**
     * \brief Class container for the Broyden's method.
     *
     * \includedoc docs/markdown/RootFinder/Broyden.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class Broyden : public QuasiNewton<Vector, Broyden<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename TypeTrait<Vector>::Scalar;
      using typename QuasiNewton<Vector, Broyden<Vector>>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Broyden solver type enumeration.
       */
      using Method = enum class Method : Integer {
        GOOD     = 0,
        BAD      = 1,
        COMBINED = 2
      };

     private:
      Method m_method{Method::COMBINED}; /**< Broyden solver method. */

     public:
      /**
       * Class constructor for the Broyden solver.
       */
      Broyden() {}

      /**
       * Get the Broyden solver name.
       * \return The Broyden solver name.
       */
      constexpr std::string name_impl() const {
        return "Broyden";
      }

      /**
       * Get the enumeration type of the Broyden solver method.
       * \return The Broyden solver enumeration type.
       */
      Method method() const {
        return this->m_method;
      }

      /**
       * Set the enumeration type of the Broyden solver method.
       * \param[in] t_method The Broyden solver enumeration type.
       */
      void method(Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Enable the \em good Broyden solver method.
       */
      void enable_good_method() {
        this->m_method = Method::GOOD;
      }

      /**
       * Enable the \em bad Broyden solver method.
       */
      void enable_bad_method() {
        this->m_method = Method::BAD;
      }

      /**
       * Enable the \em combined Broyden solver method.
       */
      void enable_combined_method() {
        this->m_method = Method::COMBINED;
      }


      /**
       * Jacobian approximation update rule for the Broyden's method.
       * \param[in] delta_x_old Old difference between points.
       * \param[in] delta_function_old Old difference between function values.
       * \param[in] jacobian_old Old jacobian approximation.
       * \param[in] delta_x_new New difference between points.
       * \param[in] delta_function_new New difference between function values.
       * \param[in] function_new New function value.
       * \param[out] jacobian_new New jacobian approximation.
       */
      void update_impl(const Vector &delta_x_old,
                       const Vector &delta_function_old,
                       const FirstDerivative &jacobian_old,
                       const Vector &delta_x_new,
                       const Vector &delta_function_new,
                       [[maybe_unused]] const Vector &function_new,
                       FirstDerivative &jacobian_new) {
        Vector tmp_1(jacobian_old * delta_function_new);
        Scalar tmp_2{delta_function_new.squaredNorm()};

        // Select the update formula. The combined method uses the criterion of
        // Martinez and Ochi: use the good update if
        // |(dx_new'*dx_old) / (dx_new'*J_old*dF_new)| <
        // |(dF_new'*dF_old) / (dF_new'*dF_new)|
        bool use_good{true};
        if (this->m_method == Method::BAD) {
          use_good = false;
        } else if (this->m_method == Method::COMBINED &&
                   this->iterations() >= Integer(2)) {
          use_good = std::abs(delta_x_new.dot(delta_x_old)) /
                         std::abs(delta_x_new.dot(tmp_1)) <
                     std::abs(delta_function_new.dot(delta_function_old)) /
                         tmp_2;
        }

        if (use_good) {
          // Broyden's good method: J_new = J_old - (J_old*dF_new-dx_new) /
          // (C'*dF_new)*C', where C = J_old'*dx_new
          Vector C_g(jacobian_old.transpose() * delta_x_new);
          Scalar den{C_g.dot(delta_function_new)};
          if (std::abs(den) <= EPSILON || !std::isfinite(den)) {
            jacobian_new = jacobian_old;
            return;
          }
          jacobian_new =
              jacobian_old - (tmp_1 - delta_x_new) / den * C_g.transpose();
        } else {
          // Broyden's bad method: J_new = J_old - (J_old*dF_new-dx_new) /
          // (dF_new'*dF_new)*dF_new'
          if (tmp_2 <= EPSILON || !std::isfinite(tmp_2)) {
            jacobian_new = jacobian_old;
            return;
          }
          jacobian_new = jacobian_old - (tmp_1 - delta_x_new) / tmp_2 *
                                            delta_function_new.transpose();
        }
      }

    };  // class Broyden

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_BROYDEN_HH
