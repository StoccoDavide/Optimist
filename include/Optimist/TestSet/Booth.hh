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
    requires TypeTrait<Vector>::IsEigen && (TypeTrait<Vector>::IsDynamicSize ||
      (TypeTrait<Vector>::IsFixedSize && Vector::RowsAtCompileTime == 2))
    class Booth : public Function<Vector, Vector, Booth<Vector>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar = typename Vector::Scalar;
      using typename Function<Vector, Vector, Booth<Vector>>::Matrix;
      using typename Function<Vector, Vector, Booth<Vector>>::Tensor;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Booth function.
       */
      Booth()
      {
        this->m_solutions.emplace_back(1.0, 3.0);
        for (Scalar x{-10.0}; x < 10.0 + EPSILON; x += 5.0) {
          for (Scalar y{-10.0}; y < 10.0 + EPSILON; y += 5.0) {
            this->m_guesses.emplace_back(x, y);
          }
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

        if constexpr (VectorTrait::IsFixedSize) {
          out << x(0) + 2.0*x(1) - 7.0, 2.0*x(0) + x(1) - 5.0;
        } else if constexpr (VectorTrait::IsDynamicSize) {
          out.resize(2);
          out << x(0) + 2.0*x(1) - 7.0, 2.0*x(0) + x(1) - 5.0;
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(2); out.reserve(2);
          std::vector<Eigen::Triplet<Scalar>> triplets; triplets.reserve(2);
          triplets.emplace_back(0, 0, x.coeff(0) + 2.0*x.coeff(1) - 7.0);
          triplets.emplace_back(1, 0, 2.0*x.coeff(0) + x.coeff(1) - 5.0);
          out.setFromTriplets(triplets.begin(), triplets.end());
        } else {
          static_assert(VectorTrait::IsEigen, CMD "input type not supported.");
        }
        return out.allFinite();

        #undef CMD
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Vector const & /*x*/, Matrix & out) const
      {
        #define CMD "Optimist::TestSet::Booth::first_derivative_impl(...): "

        if constexpr (VectorTrait::IsFixedSize) {
          out << 1.0, 2.0, 2.0, 1.0;
        } else if constexpr (VectorTrait::IsDynamicSize) {
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
          static_assert(VectorTrait::IsEigen, CMD "input type not supported.");
        }
        return out.allFinite();

        #undef CMD
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Vector const & /*x*/, Tensor & out) const
      {
        out.resize(2);
        std::for_each(out.begin(), out.end(), [] (Matrix & m) {
          if constexpr (VectorTrait::IsDynamicSize) {m.resize(2, 2);}
          m.setZero();
        });
        return true;
      }

    }; // class Booth

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_BOOTH_HH
