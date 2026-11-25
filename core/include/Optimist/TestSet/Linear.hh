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

#ifndef OPTIMIST_TESTSET_LINEAR_HH
#define OPTIMIST_TESTSET_LINEAR_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the linear function.
     *
     * Class container for the linear function, which is defined as:
     * \f[
     * f(x) = mx + q \text{.}
     * \f]
     * The function has roots at
     * \f[
     * x = -\displaystyle\frac{q}{m} \text{.}
     * \f]
     * The default coefficients are \f$m = 1\f$ and \f$q = 1\f$, and the
     * function initial guess is
     * \f$x = 0\f$.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class Linear : public Function<Scalar, Scalar, Linear<Scalar>> {
      Scalar m_m{-1.0}; /**< Coefficient \f$ m \f$. */
      Scalar m_q{1.0};  /**< Coefficient \f$ q \f$. */

     public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the linear function.
       */
      Linear() {
        this->m_solutions.emplace_back(-this->m_q / this->m_m);
        this->m_guesses.emplace_back(0.0);
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Linear";
      }

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(const Scalar x, Scalar &out) const {
        out = this->m_m * x + this->m_q;
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl([[maybe_unused]] const Scalar x, Scalar &out) const {
        out = this->m_m;
        return std::isfinite(out);
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl([[maybe_unused]] const Scalar x, Scalar &out) const {
        out = 0.0;
        return std::isfinite(out);
      }

    };  // class Linear

    /**
     * \brief Class container for the one-dimensional linear system.
     *
     * Vector counterpart of Linear, \f$ f(x) = mx + q \f$ with \f$ m = q = 1
     * \f$, used to test the vector solvers in the smallest dimension.
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class Linear1 : public DenseSystemFunction<Vector, 1, Linear1<Vector>> {
     public:
      using Base = DenseSystemFunction<Vector, 1, Linear1<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

     private:
      Scalar m_m{1.0}; /**< Coefficient \f$ m \f$. */
      Scalar m_q{1.0}; /**< Coefficient \f$ q \f$. */

     public:
      /**
       * Class constructor for the one-dimensional linear system.
       */
      Linear1() {
        this->add_solution(this->make({-this->m_q / this->m_m}));
        this->add_guess(this->make({0.0}));
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Linear1";
      }

      void residual(const Dense &x, Dense &F) const {
        F(0) = this->m_m * x(0) + this->m_q;
      }

      void jacobian([[maybe_unused]] const Dense &x, DenseMatrix &J) const {
        J(0, 0) = this->m_m;
      }

      void hessian([[maybe_unused]] const Dense &x,
                   [[maybe_unused]] const Integer i,
                   DenseMatrix &H) const {
        H.setZero();
      }

    };  // class Linear1

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_LINEAR_HH
