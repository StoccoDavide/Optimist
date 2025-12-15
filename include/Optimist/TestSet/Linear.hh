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

#ifndef OPTIMIST_TESTSET_LINEAR_HH
#define OPTIMIST_TESTSET_LINEAR_HH

#include "Optimist/Function.hh"

namespace Optimist
{

  namespace TestSet
  {

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
     * The default coefficients are \f$m = 1\f$ and \f$q = 1\f$, and the function initial guess is
     * \f$x = 0\f$.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    requires TypeTrait<Scalar>::IsScalar
    class Linear : public Function<Scalar, Scalar, Linear<Scalar>>
    {
      Scalar m_m{-1.0}; /**< Coefficient \f$ m \f$. */
      Scalar m_q{1.0}; /**< Coefficient \f$ q \f$. */

    public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the linear function.
       */
      Linear()
      {
        this->m_solutions.emplace_back(-this->m_q/this->m_m);
        this->m_guesses.emplace_back(0.0);
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "Linear";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Scalar const x, Scalar & out) const
      {
        out = this->m_m*x + this->m_q;
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Scalar const /*x*/, Scalar & out) const
      {
        out = this->m_m;
        return std::isfinite(out);
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Scalar const /*x*/, Scalar & out) const
      {
        out = 0.0;
        return std::isfinite(out);
      }

    }; // class Linear


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
     * The default coefficients are \f$m = 1\f$ and \f$q = 1\f$, and the function initial guess is
     * \f$x = 0\f$.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    requires TypeTrait<Vector>::IsEigen &&
      (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == 1)
    class Linear1 : public Function<Vector, Vector, Linear1<Vector>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector, Vector, Linear1<Vector>>::FirstDerivative;
      using typename Function<Vector, Vector, Linear1<Vector>>::SecondDerivative;

    private:
      Scalar m_m{1.0}; /**< Coefficient \f$ m \f$. */
      Scalar m_q{1.0}; /**< Coefficient \f$ q \f$. */

    public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the linear function.
       */
      Linear1()
      {
        this->m_solutions.resize(1);
        this->m_guesses.resize(1);
        Scalar tmp{-this->m_q/this->m_m};
        if constexpr (VectorTrait::IsFixed) {
          this->m_solutions[0] << tmp;
          this->m_guesses[0] << 0.0;
        } else if constexpr (VectorTrait::IsDynamic) {
          this->m_solutions[0].resize(1);
          this->m_solutions[0] << tmp;
          this->m_guesses[0].resize(1); this->m_guesses[0] << 0.0;
        } else if constexpr (VectorTrait::IsSparse) {
          this->m_solutions[0].resize(1); this->m_solutions[0].reserve(1);
          this->m_solutions[0].coeffRef(0) = tmp;
          this->m_guesses[0].resize(1); this->m_guesses[0].reserve(1);
          this->m_guesses[0].coeffRef(0) = 0.0;
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "Linear1";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Vector const & x, Vector & out) const
      {
        #define CMD "Optimist::TestSet::Linear1::evaluate_impl(...): "

        if constexpr (VectorTrait::IsFixed ) {
          out << this->m_m*x(0) + this->m_q;
          return std::isfinite(out(0));
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(1);
          out << this->m_m*x(0) + this->m_q;
          return std::isfinite(out(0));
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(1); out.reserve(1);
          out.coeffRef(0) = this->m_m*x.coeff(0) + this->m_q;
          return std::isfinite(out.coeff(0));
        } else {
          OPTIMIST_ERROR(CMD "input type not supported.");
          return false;
        }

        #undef CMD
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Vector const & /*x*/, FirstDerivative & out) const
      {
        #define CMD "Optimist::TestSet::Linear1::first_derivative_impl(...): "

        if constexpr (VectorTrait::IsFixed) {
          out << this->m_m;
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(1, 1);
          out << this->m_m;
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(1, 1); out.reserve(1);
          out.coeffRef(0, 0) = this->m_m;
        } else {
          OPTIMIST_ERROR(CMD "input type not supported.");
          return false;
        }
        return std::isfinite(this->m_m);

        #undef CMD
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Vector const & /*x*/, SecondDerivative & out) const
      {
        #define CMD "Optimist::TestSet::Linear1::second_derivative_impl(...): "

        out.resize(1);
        if constexpr (VectorTrait::IsFixed) {
          out[0] << 0.0;
        } else if constexpr (VectorTrait::IsDynamic) {
          out[0].resize(1, 1);
          out[0] << 0.0;
        } else if constexpr (VectorTrait::IsSparse) {
          out[0].resize(1, 1); out[0].reserve(1);
          out[0].coeffRef(0, 0) = 0.0;
        } else {
          OPTIMIST_ERROR(CMD "input type not supported.");
          return false;
        }
        return true;

        #undef CMD
      }

    }; // class Linear1

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_LINEAR_HH
