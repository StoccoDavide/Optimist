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

#ifndef OPTIMIST_TESTSET_SCHAFFER2_HH
#define OPTIMIST_TESTSET_SCHAFFER2_HH

#include "Optimist/Function.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the Schaffer2 function.
     *
     * Class container for the Schaffer2 function, which is defined as:
     * \f[
     * f(\mathbf{x}) = 0.5 + \displaystyle\frac{\sin^{2}(x_1^2 - x_2^2) -
     * 0.5}{(1 + 0.001(x_1^2 + x_2^2))^2} \text{.}
     * \f]
     * The function has global minima at \f$\mathbf{x} = (0, 0)\f$, with
     * \f$f(\mathbf{x}) = 0\f$. The initial guesses are generated on the square
     * \f$x_i \in \left[-10, 10\right]\f$.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen && (!TypeTrait<Vector>::IsFixed ||
                                              TypeTrait<Vector>::Dimension == 2)
    class Schaffer2
        : public Function<Vector, typename Vector::Scalar, Schaffer2<Vector>> {
     public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector,
                              typename Vector::Scalar,
                              Schaffer2<Vector>>::FirstDerivative;
      using typename Function<Vector,
                              typename Vector::Scalar,
                              Schaffer2<Vector>>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Schaffer2 function.
       */
      Schaffer2() {
        this->m_solutions.resize(1);
        this->m_guesses.resize(2);
        this->m_solutions.resize(1);
        this->m_guesses.resize(2);
        if constexpr (VectorTrait::IsFixed) {
          this->m_solutions[0] << 0.0, 0.0;
          this->m_guesses[0] << -10, -10;
          this->m_guesses[1] << 10, 10;
        } else if constexpr (!VectorTrait::IsSparse) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0] << 0.0, 0.0;
          this->m_guesses[0].resize(2);
          this->m_guesses[0] << -10, -10;
          this->m_guesses[1].resize(2);
          this->m_guesses[1] << 10, 10;
        } else if constexpr (VectorTrait::IsSparse) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0].reserve(2);
          this->m_solutions[0].coeffRef(0) = 0.0;
          this->m_solutions[0].coeffRef(1) = 0.0;
          this->m_guesses[0].resize(2);
          this->m_guesses[0].reserve(2);
          this->m_guesses[0].coeffRef(0) = -10;
          this->m_guesses[0].coeffRef(1) = -10;
          this->m_guesses[1].resize(2);
          this->m_guesses[1].reserve(2);
          this->m_guesses[1].coeffRef(0) = 10;
          this->m_guesses[1].coeffRef(1) = 10;
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Schaffer2";
      }

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
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
        out = 0.5e0 + (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) *
                          pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -2);
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
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
          out << (-0.4e-2 * x_0 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
                  (0.4e-2 * pow(x_0, 3) + (0.4e1 + 0.4e-2 * x_1 * x_1) * x_0) *
                      cos(x_0 * x_0 - x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
                  0.20e-2 * x_0) *
                     pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -3),
              (-0.4e-2 * x_1 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (-0.4e-2 * pow(x_1, 3) + (-0.4e1 - 0.4e-2 * x_0 * x_0) * x_1) *
                   cos(-x_0 * x_0 + x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               0.20e-2 * x_1) *
                  pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -3);
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2);
          out << (-0.4e-2 * x_0 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
                  (0.4e-2 * pow(x_0, 3) + (0.4e1 + 0.4e-2 * x_1 * x_1) * x_0) *
                      cos(x_0 * x_0 - x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
                  0.20e-2 * x_0) *
                     pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -3),
              (-0.4e-2 * x_1 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (-0.4e-2 * pow(x_1, 3) + (-0.4e1 - 0.4e-2 * x_0 * x_0) * x_1) *
                   cos(-x_0 * x_0 + x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               0.20e-2 * x_1) *
                  pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -3);

        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2);
          out.reserve(2);
          out.coeffRef(0) =
              (-0.4e-2 * x_0 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (0.4e-2 * pow(x_0, 3) + (0.4e1 + 0.4e-2 * x_1 * x_1) * x_0) *
                   cos(x_0 * x_0 - x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               0.20e-2 * x_0) *
              pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -3);
          out.coeffRef(1) =
              (-0.4e-2 * x_1 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (-0.4e-2 * pow(x_1, 3) + (-0.4e1 - 0.4e-2 * x_0 * x_0) * x_1) *
                   cos(-x_0 * x_0 + x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               0.20e-2 * x_1) *
              pow(1 + 0.1e-2 * x_0 * x_0 + 0.1e-2 * x_1 * x_1, -3);
        }

        return out.allFinite();
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(const Vector &x,
                                  SecondDerivative &out) const {
        Scalar x_0, x_1;
        if constexpr (VectorTrait::IsSparse) {
          x_0 = x.coeff(0);
          x_1 = x.coeff(1);
        } else {
          x_0 = x(0);
          x_1 = x(1);
        }
        if constexpr (VectorTrait::IsFixed) {
          out.resize(2, 2);
          out(0, 0) =
              ((-0.8000000e7 * pow(x_0, 6) +
                (-0.1600000000e11 - 0.16000000e8 * x_1 * x_1) * pow(x_0, 4) -
                0.8000000e7 * (x_1 * x_1 + 0.100158113883042665e4) *
                    (x_1 * x_1 + 0.998418861169573461e3) * x_0 * x_0 -
                0.4000000e7 * x_1 * x_1 - 0.4000000000e10) *
                   pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (-0.28000000e8 * pow(x_0, 4) +
                (-0.2400000000e11 - 0.24000000e8 * x_1 * x_1) * x_0 * x_0 +
                0.4000000e7 * pow(x_1 * x_1 + 0.9999999999e3, 2)) *
                   cos(x_0 * x_0 - x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               (0.8000000e7 * pow(x_0, 6) +
                (0.1600000000e11 + 0.16000000e8 * x_1 * x_1) * pow(x_0, 4) +
                0.8000000e7 * pow(x_1 * x_1 + 0.9999999999e3, 2) * x_0 * x_0) *
                   pow(cos(x_0 * x_0 - x_1 * x_1), 2) -
               0.100000000e8 * x_0 * x_0 + 0.20000000e7 * x_1 * x_1 +
               0.2000000000e10) *
              pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out(0, 1) = -8000000 * x_1 * cos(-x_0 * x_0 + x_1 * x_1) * x_0 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          sin(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) *
                          x_0 * x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out(1, 0) = -8000000 * x_1 * cos(-x_0 * x_0 + x_1 * x_1) * x_0 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          sin(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) *
                          x_0 * x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out(1, 1) =
              8000000 * x_1 * x_1 * pow(cos(-x_0 * x_0 + x_1 * x_1), 2) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              4000000 * sin(x_0 * x_0 - x_1 * x_1) *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) +
              8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_1 * x_1 *
                  sin(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) +
              32000000 * sin(x_0 * x_0 - x_1 * x_1) * x_1 * x_1 *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) * x_1 *
                  x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4) +
              (-4000000 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) + 0.20000000e7) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3);
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2, 2);
          out(0, 0) =
              ((-0.8000000e7 * pow(x_0, 6) +
                (-0.1600000000e11 - 0.16000000e8 * x_1 * x_1) * pow(x_0, 4) -
                0.8000000e7 * (x_1 * x_1 + 0.100158113883042665e4) *
                    (x_1 * x_1 + 0.998418861169573461e3) * x_0 * x_0 -
                0.4000000e7 * x_1 * x_1 - 0.4000000000e10) *
                   pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (-0.28000000e8 * pow(x_0, 4) +
                (-0.2400000000e11 - 0.24000000e8 * x_1 * x_1) * x_0 * x_0 +
                0.4000000e7 * pow(x_1 * x_1 + 0.9999999999e3, 2)) *
                   cos(x_0 * x_0 - x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               (0.8000000e7 * pow(x_0, 6) +
                (0.1600000000e11 + 0.16000000e8 * x_1 * x_1) * pow(x_0, 4) +
                0.8000000e7 * pow(x_1 * x_1 + 0.9999999999e3, 2) * x_0 * x_0) *
                   pow(cos(x_0 * x_0 - x_1 * x_1), 2) -
               0.100000000e8 * x_0 * x_0 + 0.20000000e7 * x_1 * x_1 +
               0.2000000000e10) *
              pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out(0, 1) = -8000000 * x_1 * cos(-x_0 * x_0 + x_1 * x_1) * x_0 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          sin(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) *
                          x_0 * x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out(1, 0) = -8000000 * x_1 * cos(-x_0 * x_0 + x_1 * x_1) * x_0 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          sin(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(x_0 * x_0 - x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                          cos(-x_0 * x_0 + x_1 * x_1) *
                          pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
                      24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) *
                          x_0 * x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out(1, 1) =
              8000000 * x_1 * x_1 * pow(cos(-x_0 * x_0 + x_1 * x_1), 2) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              4000000 * sin(x_0 * x_0 - x_1 * x_1) *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) +
              8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_1 * x_1 *
                  sin(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) +
              32000000 * sin(x_0 * x_0 - x_1 * x_1) * x_1 * x_1 *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) * x_1 *
                  x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4) +
              (-4000000 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) + 0.20000000e7) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3);
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2, 2);
          out.reserve(4);
          out.coeffRef(0, 0) =
              ((-0.8000000e7 * pow(x_0, 6) +
                (-0.1600000000e11 - 0.16000000e8 * x_1 * x_1) * pow(x_0, 4) -
                0.8000000e7 * (x_1 * x_1 + 0.100158113883042665e4) *
                    (x_1 * x_1 + 0.998418861169573461e3) * x_0 * x_0 -
                0.4000000e7 * x_1 * x_1 - 0.4000000000e10) *
                   pow(sin(x_0 * x_0 - x_1 * x_1), 2) +
               (-0.28000000e8 * pow(x_0, 4) +
                (-0.2400000000e11 - 0.24000000e8 * x_1 * x_1) * x_0 * x_0 +
                0.4000000e7 * pow(x_1 * x_1 + 0.9999999999e3, 2)) *
                   cos(x_0 * x_0 - x_1 * x_1) * sin(x_0 * x_0 - x_1 * x_1) +
               (0.8000000e7 * pow(x_0, 6) +
                (0.1600000000e11 + 0.16000000e8 * x_1 * x_1) * pow(x_0, 4) +
                0.8000000e7 * pow(x_1 * x_1 + 0.9999999999e3, 2) * x_0 * x_0) *
                   pow(cos(x_0 * x_0 - x_1 * x_1), 2) -
               0.100000000e8 * x_0 * x_0 + 0.20000000e7 * x_1 * x_1 +
               0.2000000000e10) *
              pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out.coeffRef(0, 1) =
              -8000000 * x_1 * cos(-x_0 * x_0 + x_1 * x_1) * x_0 *
                  cos(x_0 * x_0 - x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                  sin(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                  cos(x_0 * x_0 - x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) * x_0 *
                  x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out.coeffRef(1, 0) =
              -8000000 * x_1 * cos(-x_0 * x_0 + x_1 * x_1) * x_0 *
                  cos(x_0 * x_0 - x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                  sin(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                  cos(x_0 * x_0 - x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              16000000 * sin(x_0 * x_0 - x_1 * x_1) * x_0 * x_1 *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) * x_0 *
                  x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4);
          out.coeffRef(1, 1) =
              8000000 * x_1 * x_1 * pow(cos(-x_0 * x_0 + x_1 * x_1), 2) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) -
              4000000 * sin(x_0 * x_0 - x_1 * x_1) *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) +
              8000000 * sin(x_0 * x_0 - x_1 * x_1) * x_1 * x_1 *
                  sin(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -2) +
              32000000 * sin(x_0 * x_0 - x_1 * x_1) * x_1 * x_1 *
                  cos(-x_0 * x_0 + x_1 * x_1) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3) +
              24000000 * (pow(sin(x_0 * x_0 - x_1 * x_1), 2) - 0.5e0) * x_1 *
                  x_1 * pow(x_0 * x_0 + x_1 * x_1 + 1000, -4) +
              (-4000000 * pow(sin(x_0 * x_0 - x_1 * x_1), 2) + 0.20000000e7) *
                  pow(x_0 * x_0 + x_1 * x_1 + 1000, -3);
        }

        return out.allFinite();
      }

    };  // class Schaffer2

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_SCHAFFER2_HH
