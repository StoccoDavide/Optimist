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

#ifndef OPTIMIST_TESTSET_BROWN_HH
#define OPTIMIST_TESTSET_BROWN_HH

#include "Optimist/Function.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the Brown badly scaled function.
     *
     * Class container for the Brown badly scaled function, which is defined as:
     * \f[
     * \mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 - a \\ x_2 - 2a \\ x_1x_2 -
     * 2
     * \end{bmatrix} \text{,}
     * \f]
     * where \f$a = 10^{-6}\f$. The function has one solution at \f$\mathbf{x} =
     * [a, 2a]^\top\f$, with \f$f(\mathbf{x}) = 0\f$. The initial guess is
     * generated at
     * \f$\mathbf{x} = [1, 1]^\top\f$.
     * \tparam Input Input vector type.
     * \tparam Output Output vector type.
     */
    template <typename Input, typename Output>
      requires TypeTrait<Input>::IsEigen && TypeTrait<Output>::IsEigen &&
               (!TypeTrait<Input>::IsFixed ||
                TypeTrait<Input>::Dimension == 2) &&
               (!TypeTrait<Output>::IsFixed ||
                TypeTrait<Output>::Dimension == 3)
    class Brown : public Function<Input, Output, Brown<Input, Output>> {
     public:
      using InputTrait  = TypeTrait<Input>;
      using OutputTrait = TypeTrait<Output>;
      using Scalar      = typename Input::Scalar;
      using typename Function<Input, Output, Brown<Input, Output>>::
          FirstDerivative;
      using typename Function<Input, Output, Brown<Input, Output>>::
          SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_a{
        1.0e-6}; /**< Scaling value (keep it low to guarantee bad scaling). */

     public:
      /**
       * Class constructor for the Brown function.
       */
      Brown() {
        this->m_solutions.resize(1);
        this->m_guesses.resize(1);
        if constexpr (InputTrait::IsFixed) {
          this->m_solutions[0] << this->m_a, 2.0 * this->m_a;
          this->m_guesses[0] << 1.0, 1.0;
        } else if constexpr (InputTrait::IsDynamic) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0] << this->m_a, 2.0 * this->m_a;
          this->m_guesses[0].resize(2);
          this->m_guesses[0] << 1.0, 1.0;
        } else if constexpr (InputTrait::IsSparse) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0].reserve(2);
          this->m_solutions[0].coeffRef(0) = this->m_a;
          this->m_solutions[0].coeffRef(1) = 2.0 * this->m_a;
          this->m_guesses[0].resize(2);
          this->m_guesses[0].reserve(2);
          this->m_guesses[0].coeffRef(0) = 1.0;
          this->m_guesses[0].coeffRef(1) = 1.0;
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Brown";
      }

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(const Input &x, Output &out) const {
        Scalar x_0, x_1;
        if constexpr (InputTrait::IsSparse) {
          x_0 = x.coeff(0);
          x_1 = x.coeff(1);
        } else {
          x_0 = x(0);
          x_1 = x(1);
        }

        if constexpr (OutputTrait::IsFixed) {
          out << x_0 - this->m_a, x_1 - 2.0 * this->m_a, x_0 * x_1 - 2.0;
          return out.allFinite();
        } else if constexpr (OutputTrait::IsDynamic) {
          out.resize(3);
          out << x_0 - this->m_a, x_1 - 2.0 * this->m_a, x_0 * x_1 - 2.0;
          return out.allFinite();
        } else if constexpr (OutputTrait::IsSparse) {
          out.resize(3);
          out.reserve(3);
          out.coeffRef(0) = x_0 - this->m_a;
          out.coeffRef(1) = x_1 - 2.0 * this->m_a;
          out.coeffRef(2) = x_0 * x_1 - 2.0;
          for (typename Output::InnerIterator it(out); it; ++it) {
            if (!std::isfinite(it.value())) {
              return false;
            }
          }
          return true;
        }
        return false;
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(const Input &x, FirstDerivative &out) const {
        Scalar x_0, x_1;
        if constexpr (InputTrait::IsSparse) {
          x_0 = x.coeff(0);
          x_1 = x.coeff(1);
        } else {
          x_0 = x(0);
          x_1 = x(1);
        }

        if constexpr (InputTrait::IsFixed) {
          out << 1.0, 0.0, x_1, 0.0, 1.0, x_0;
          return out.allFinite();
        } else if constexpr (InputTrait::IsDynamic) {
          out.resize(3, 2);
          out << 1.0, 0.0, x_1, 0.0, 1.0, x_0;
          return out.allFinite();
        } else if constexpr (InputTrait::IsSparse) {
          out.resize(3, 2);
          out.reserve(4);
          out.coeffRef(0, 0) = 1.0;
          out.coeffRef(1, 1) = 1.0;
          out.coeffRef(2, 0) = x_1;
          out.coeffRef(2, 1) = x_0;
          for (Integer k{0}; k < out.outerSize(); ++k) {
            for (typename FirstDerivative::InnerIterator it(out, k); it; ++it) {
              if (!std::isfinite(it.value())) {
                return false;
              }
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
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(const Input & /*x*/,
                                  SecondDerivative &out) const {
        out.resize(2);
        if constexpr (InputTrait::IsFixed) {
          out[0].setZero();
          out[1].setZero();
          out[0](2, 1) = 1.0;
          out[1](2, 0) = 1.0;
          return out[0].allFinite() && out[1].allFinite();
        } else if constexpr (InputTrait::IsDynamic) {
          for (auto &matrix : out) {
            matrix.resize(3, 2);
            matrix.setZero();
          }
          out[0](2, 1) = 1.0;
          out[1](2, 0) = 1.0;
          return out[0].allFinite() && out[1].allFinite();
        } else if constexpr (InputTrait::IsSparse) {
          for (auto &matrix : out) {
            matrix.resize(3, 2);
            matrix.reserve(1);
          }
          out[0].coeffRef(2, 1) = 1.0;
          out[1].coeffRef(2, 0) = 1.0;
          for (const auto &matrix : out) {
            for (Integer k{0}; k < matrix.outerSize(); ++k) {
              for (typename SecondDerivative::value_type::InnerIterator it(
                       matrix,
                       k);
                   it;
                   ++it) {
                if (!std::isfinite(it.value())) {
                  return false;
                }
              }
            }
          }
          return true;
        }
        return false;
      }

    };  // class Brown

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_BROWN_HH
