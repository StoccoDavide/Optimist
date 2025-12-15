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

#ifndef OPTIMIST_TESTSET_SCHAFFER2_HH
#define OPTIMIST_TESTSET_SCHAFFER2_HH

#include "Optimist/Function.hh"

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
     * The function has global minima at \f$\mathbf{x} = (0, 0)\f$, with \f$f(\mathbf{x}) = 0\f$.
     * The initial guesses are generated on the square \f$x_i \in \left[-10, 10\right]\f$.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    requires TypeTrait<Vector>::IsEigen &&
      (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == 2)
    class Schaffer2 : public Function<Vector, typename Vector::Scalar, Schaffer2<Vector>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector, typename Vector::Scalar, Schaffer2<Vector>>::FirstDerivative;
      using typename Function<Vector, typename Vector::Scalar, Schaffer2<Vector>>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Schaffer2 function.
       */
      Schaffer2()
      {
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
          this->m_guesses[0].resize(2); this->m_guesses[0] << -10, -10;
          this->m_guesses[1].resize(2); this->m_guesses[1] << 10, 10;
        } else if constexpr (VectorTrait::IsSparse) {
          this->m_solutions[0].resize(2); this->m_solutions[0].reserve(2);
          this->m_solutions[0].coeffRef(0) = 0.0;
          this->m_solutions[0].coeffRef(1) = 0.0;
          this->m_guesses[0].resize(2); this->m_guesses[0].reserve(2);
          this->m_guesses[0].coeffRef(0) = -10;
          this->m_guesses[0].coeffRef(1) = -10;
          this->m_guesses[1].resize(2); this->m_guesses[1].reserve(2);
          this->m_guesses[1].coeffRef(0) = 10;
          this->m_guesses[1].coeffRef(1) = 10;
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "Schaffer2";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Vector const & x, Scalar & out) const
      {
        Scalar xx_0, xx_1;
        if constexpr (VectorTrait::IsSparse) {
          xx_0 = x.coeff(0)*x.coeff(0);
          xx_1 = x.coeff(1)*x.coeff(1);
        } else {
          xx_0 = x(0)*x(0);
          xx_1 = x(1)*x(1);
        }
        Scalar xx_0_m_xx_1{xx_0 - xx_1};
        Scalar xx_0_p_xx_1{xx_0 + xx_1};
        out = 0.5 + (std::sin(xx_0_m_xx_1)*std::sin(xx_0_m_xx_1) - 0.5) /
          ((1.0 + 0.001*(xx_0_p_xx_1))*(1.0 + 0.001*(xx_0_p_xx_1)));
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Vector const & x, FirstDerivative & out) const
      {
        Scalar xx_0, xx_1;
        if constexpr (VectorTrait::IsSparse) {
          xx_0 = x.coeff(0)*x.coeff(0);
          xx_1 = x.coeff(1)*x.coeff(1);
        } else {
          xx_0 = x(0)*x(0);
          xx_1 = x(1)*x(1);
        }
        Scalar xx_0_m_xx_1{xx_0 - xx_1};
        Scalar xx_0_p_xx_1{xx_0 + xx_1};
        Scalar tmp{static_cast<Scalar>(1.0) + static_cast<Scalar>(0.001)*xx_0_p_xx_1};
        Scalar tmp2{tmp*tmp}, tmp3{tmp2*tmp};
        if constexpr (VectorTrait::IsFixed) {
          out << 2.0*x(0)*std::sin(xx_0_m_xx_1) / tmp2 - 2.0*0.001*x(0)*std::cos(xx_0_m_xx_1) / tmp3,
            -2.0*x(1)*std::sin(xx_0_m_xx_1) / tmp2 + 2.0*0.001*x(1)*std::cos(xx_0_m_xx_1) / tmp3;
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2);
          out << 2.0*x(0)*std::sin(xx_0_m_xx_1) / tmp2 - 2.0*0.001*x(0)*std::cos(xx_0_m_xx_1) / tmp3,
            -2.0*x(1)*std::sin(xx_0_m_xx_1) / tmp2 + 2.0*0.001*x(1)*std::cos(xx_0_m_xx_1) / tmp3;
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2); out.reserve(2);
          out.coeffRef(0) = +2.0*x(0)*std::sin(xx_0_m_xx_1) / tmp2 - 2.0*0.001*x(0)*std::cos(xx_0_m_xx_1) / tmp3;
          out.coeffRef(1) = -2.0*x(1)*std::sin(xx_0_m_xx_1) / tmp2 + 2.0*0.001*x(1)*std::cos(xx_0_m_xx_1) / tmp3;
        }

        return out.allFinite();
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Vector const & x, SecondDerivative & out) const
      {
        Scalar xx_0, xx_1;
        if constexpr (VectorTrait::IsSparse) {
          xx_0 = x.coeff(0)*x.coeff(0);
          xx_1 = x.coeff(1)*x.coeff(1);
        } else {
          xx_0 = x(0)*x(0);
          xx_1 = x(1)*x(1);
        }
        Scalar xx_0_m_xx_1{xx_0 - xx_1};
        Scalar xx_0_p_xx_1{xx_0 + xx_1};
        Scalar tmp{static_cast<Scalar>(1.0) + static_cast<Scalar>(0.001)*(xx_0_p_xx_1)};
        Scalar tmp2{tmp*tmp}, tmp3{tmp2*tmp}, tmp4{tmp2*tmp2};
        if constexpr (VectorTrait::IsFixed) {
          out.resize(2, 2);
          out(0, 0) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
            4.0*x(0)*x(0)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(0)*x(0)*std::cos(xx_0_m_xx_1) / tmp4;
          out(0, 1) = -2.0*x(0)*x(1)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(0)*x(1)*std::cos(xx_0_m_xx_1) / tmp4;
          out(1, 0) = out(0, 1);
          out(1, 1) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
            4.0*x(1)*x(1)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(1)*x(1)*std::cos(xx_0_m_xx_1) / tmp4;
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2, 2);
          out(0, 0) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
            4.0*x(0)*x(0)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(0)*x(0)*std::cos(xx_0_m_xx_1) / tmp4;
          out(0, 1) = -2.0*x(0)*x(1)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(0)*x(1)*std::cos(xx_0_m_xx_1) / tmp4;
          out(1, 0) = out(0, 1);
          out(1, 1) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
            4.0*x(1)*x(1)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(1)*x(1)*std::cos(xx_0_m_xx_1) / tmp4;
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2, 2); out.reserve(4);
          out.coeffRef(0, 0) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
            4.0*x(0)*x(0)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x(0)*x(0)*std::cos(xx_0_m_xx_1) / tmp4;
          out.coeffRef(0, 1) = -2.0*x.coeff(0)*x.coeff(1)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x.coeff(0)*x.coeff(1)*std::cos(xx_0_m_xx_1) / tmp4;
          out.coeffRef(1, 0) = out.coeff(0, 1);
          out.coeffRef(1, 1) = 2.0*std::sin(xx_0_m_xx_1) / tmp2 -
            4.0*x.coeff(1)*x.coeff(1)*std::sin(xx_0_m_xx_1) / tmp3 +
            6.0*0.001*x.coeff(1)*x.coeff(1)*std::cos(xx_0_m_xx_1) / tmp4;
        }

        return out.allFinite();
      }

    }; // class Schaffer2

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_SCHAFFER2_HH
