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

#ifndef OPTIMIST_ROOTFINDER_GREENSTADT_HH
#define OPTIMIST_ROOTFINDER_GREENSTADT_HH

#include "Optimist/RootFinder/QuasiNewton.hh"

namespace Optimist {
  namespace RootFinder {

    /*\
     |    ____                         _            _ _
     |   / ___|_ __ ___  ___ _ __  ___| |_ __ _  __| | |_
     |  | |  _| '__/ _ \/ _ \ '_ \/ __| __/ _` |/ _` | __|
     |  | |_| | | |  __/  __/ | | \__ \ || (_| | (_| | |_
     |   \____|_|  \___|\___|_| |_|___/\__\__,_|\__,_|\__|
     |
    \*/

    /**
     * \brief Class container for the Greenstadt's method.
     *
     * \includedoc docs/markdown/RootFinder/Greenstadt.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class Greenstadt : public QuasiNewton<Vector, Greenstadt<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename TypeTrait<Vector>::Scalar;
      using typename QuasiNewton<Vector, Greenstadt<Vector>>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Greenstadt solver type enumeration.
       */
      using Method = enum class Method : Integer {
        ONE = 1,
        TWO = 2
      };

     private:
      Method m_method{Method::ONE}; /**< Greenstadt solver type. */

     public:
      /**
       * Class constructor for the Greenstadt solver.
       */
      Greenstadt() {}

      /**
       * Get the Greenstadt solver name.
       * \return The Greenstadt solver name.
       */
      constexpr std::string name_impl() const {
        return "Greenstadt";
      }

      /**
       * Get the enumeration type of the Greenstadt solver method.
       * \return The Greenstadt solver enumeration type.
       */
      Method method() const {
        return this->m_method;
      }

      /**
       * Set the enumeration type of the Greenstadt solver method.
       * \param[in] t_method The Greenstadt solver method enumeration type.
       */
      void method(Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Enable the \em Greenstadt1 solver method.
       */
      void enable_one_method() {
        this->m_method = Method::ONE;
      }

      /**
       * Enable the \em Greenstadt2 solver method.
       */
      void enable_two_method() {
        this->m_method = Method::TWO;
      }

      /**
       * Set the Greenstadt solver method.
       * \param[in] t_method The Greenstadt solver method enumeration type.
       */
      void set_method(Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Jacobian approximation update rule for the Greenstadt's method.
       * \param[in] delta_x_old Old difference between points.
       * \param[in] delta_function_old Old difference between function values.
       * \param[in] jacobian_old Old jacobian approximation.
       * \param[in] delta_x_new New difference between points.
       * \param[in] delta_function_new New difference between function values.
       * \param[in] function_new New function value.
       * \param[out] jacobian_new New jacobian approximation.
       */
      void update_impl(const Vector & /*delta_x_old*/,
                       const Vector & /*delta_function_old*/,
                       const FirstDerivative &jacobian_old,
                       const Vector &delta_x_new,
                       const Vector &delta_function_new,
                       const Vector &function_new,
                       FirstDerivative &jacobian_new) {
        if (this->m_method == Method::ONE) {
          // Greenstadt's 1st method
          // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C = F1;
          jacobian_new =
              jacobian_old - (jacobian_old * delta_function_new - delta_x_new) /
                                 (function_new.dot(delta_function_new)) *
                                 function_new.transpose();
        } else if (this->m_method == Method::TWO) {
          // Greenstadt's 2nd method
          // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C  = J0'*J0*DF1;
          Vector C(jacobian_old.transpose() * jacobian_old *
                   delta_function_new);
          jacobian_new =
              jacobian_old - (jacobian_old * delta_function_new - delta_x_new) /
                                 (C.dot(delta_function_new)) * C.transpose();
        }
      }

    };  // class Greenstadt

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_GREENSTADT_HH
