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

#ifndef OPTIMIST_TESTSET_BOOTH_HH
#define OPTIMIST_TESTSET_BOOTH_HH

#include "Optimist/Function.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
     * \brief Class container for the Booth function.
     *
     * Class container for the Booth function, which is defined as:
     * \f[
     * \mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 + 2x_2 - 7 \\ 2x_1 + x_2 - 5 \end{bmatrix} \text{.}
     * \f]
     * The function has one solution at \f$\mathbf{x} = [1, 3]^\top\f$, with \f$f(\mathbf{x}) = 0\f$.
     * The initial guesses are generated on the square \f$x_i \in [-10, 10]\f$, for all \f$x_i = 1, 2\f$.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Vector>
    requires TypeTrait<Vector>::IsEigen &&
      (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == 2)
    class Booth : public Function<Vector, Vector, Booth<Vector>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector, Vector, Booth<Vector>>::FirstDerivative;
      using typename Function<Vector, Vector, Booth<Vector>>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Booth function.
       */
      Booth()
      {
        this->m_solutions.resize(1);
        this->m_guesses.resize(2);
        if constexpr (VectorTrait::IsFixed) {
          this->m_solutions[0] << 1.0, 3.0;
          this->m_guesses[0] << -10, -10;
          this->m_guesses[1] << 10, 10;
        } else if constexpr (VectorTrait::IsDynamic) {
          this->m_solutions[0].resize(2);
          this->m_solutions[0] << 1.0, 3.0;
          this->m_guesses[0].resize(2); this->m_guesses[0] << -10, -10;
          this->m_guesses[1].resize(2); this->m_guesses[1] << 10, 10;
        } else if constexpr (VectorTrait::IsSparse) {
          this->m_solutions[0].resize(2); this->m_solutions[0].reserve(2);
          this->m_solutions[0].coeffRef(0) = 1.0;
          this->m_solutions[0].coeffRef(1) = 3.0;
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
      constexpr std::string name_impl() const {return "Booth";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Vector const & x, Vector & out) const
      {
        #define CMD "Optimist::TestSet::Booth::evaluate_impl(...): "

        if constexpr (VectorTrait::IsFixed) {
          out << x(0) + 2.0*x(1) - 7.0, 2.0*x(0) + x(1) - 5.0;
          return out.allFinite();
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2);
          out << x(0) + 2.0*x(1) - 7.0, 2.0*x(0) + x(1) - 5.0;
          return out.allFinite();
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2); out.reserve(2);
          out.coeffRef(0) = x.coeff(0) + 2.0*x.coeff(1) - 7.0;
          out.coeffRef(1) = 2.0*x.coeff(0) + x.coeff(1) - 5.0;
          for (typename Vector::InnerIterator it(out); it; ++it) {
              if (!std::isfinite(it.value())) {return false;}
          }
          return true;
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
        #define CMD "Optimist::TestSet::Booth::first_derivative_impl(...): "

        if constexpr (VectorTrait::IsFixed) {
          out << 1.0, 2.0, 2.0, 1.0;
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(2, 2);
          out << 1.0, 2.0, 2.0, 1.0;
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2, 2); out.reserve(4);
          std::vector<Eigen::Triplet<Scalar>> triplets; triplets.reserve(4);
          triplets.emplace_back(0, 0, 1.0);
          triplets.emplace_back(0, 1, 2.0);
          triplets.emplace_back(1, 0, 2.0);
          triplets.emplace_back(1, 1, 1.0);
          out.setFromTriplets(triplets.begin(), triplets.end());
        } else {
          OPTIMIST_ERROR(CMD "input type not supported.");
          return false;
        }
        return true;

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
        #define CMD "Optimist::TestSet::Booth::second_derivative_impl(...): "

        out.resize(2);
        std::for_each(out.begin(), out.end(), [] (auto & m) {
          if constexpr (VectorTrait::IsFixed) {
            m.setZero();}
          if constexpr (VectorTrait::IsDynamic) {
            m.resize(2, 2);
            m.setZero();}
          if constexpr (VectorTrait::IsSparse) {
            m.resize(2, 2);
            m.reserve(4);
            std::vector<Eigen::Triplet<Scalar>> triplets; triplets.reserve(4);
            triplets.emplace_back(0, 0, 0.0);
            triplets.emplace_back(0, 1, 0.0);
            triplets.emplace_back(1, 0, 0.0);
            triplets.emplace_back(1, 1, 0.0);
            m.setFromTriplets(triplets.begin(), triplets.end());
          } else {
            OPTIMIST_ERROR(CMD "input type not supported.");
          }
        });
        return true;

        #undef CMD
      }

    }; // class Booth

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_BOOTH_HH
