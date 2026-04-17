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
       * Set the Broyden solver type.
       * \param[in] t_method The Broyden solver type enumeration.
       */
      void set_method(Method t_method) {
        this->m_method = t_method;
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
                       const Vector & /*function_new*/,
                       FirstDerivative &jacobian_new) {
        Vector tmp_1(jacobian_old * delta_function_new);
        Scalar tmp_2{delta_function_new.squaredNorm()};
        // Selection criteria: |(dx_new'*dx_old) / (dx_new'*J_old*dF_new)| <
        // |(dF_new'*dF_old) / (dF_new'*dF_new)|
        if (this->m_method == Method::COMBINED ||
            this->m_method == Method::GOOD || this->iterations() < Integer(2) ||
            std::abs(delta_x_new.dot(delta_x_old)) /
                    std::abs(delta_x_new.dot(tmp_1)) <
                std::abs(delta_function_new.dot(delta_function_old)) / tmp_2) {
          // Broyden's Good solver: J_new = J_old - (J_old*dF_new-dx_new) /
          // (C'*dF_new)*C', where C = J_old'*dx_new;
          Vector C_g(jacobian_old.transpose() * delta_x_new);
          jacobian_new = jacobian_old - (tmp_1 - delta_x_new) /
                                            (C_g.dot(delta_function_new)) *
                                            C_g.transpose();
        } else {
          // Broyden's Bad solver: J_new = J_old - (J_old*dF_old-dx_new) /
          // (C'*dF_old)*C', where C = dF_old;
          jacobian_new = jacobian_old - (tmp_1 - delta_x_new) / tmp_2 *
                                            delta_function_old.transpose();
        }
      }

    };  // class Broyden

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_BROYDEN_HH
