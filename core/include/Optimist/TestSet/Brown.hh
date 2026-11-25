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

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Brown's badly scaled function.
     *
     * Overdetermined system of three equations in two unknowns (Moré, Garbow
     * and Hillstrom, problem 4):
     * \f[
     *   \mathbf{F}(\mathbf{x}) = \begin{bmatrix}
     *     x_1 - 10^6 \\
     *     x_2 - 2 \cdot 10^{-6} \\
     *     x_1 x_2 - 2
     *   \end{bmatrix} \text{,}
     * \f]
     * with zero residual at \f$ \mathbf{x} = (10^6, 2 \cdot 10^{-6}) \f$. The
     * twelve orders of magnitude between the components of the solution make
     * it a test for the scaling invariance of least-squares solvers.
     *
     * \tparam Input Eigen input vector type (dimension 2).
     * \tparam Output Eigen output vector type (dimension 3).
     */
    template <typename Input, typename Output>
    class Brown : public DenseSystemFunction<Input,
                                             2,
                                             Brown<Input, Output>,
                                             Output,
                                             3> {
     public:
      using Base =
          DenseSystemFunction<Input, 2, Brown<Input, Output>, Output, 3>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;
      using InputTrait  = TypeTrait<Input>;
      using OutputTrait = TypeTrait<Output>;

      /**
       * Class constructor for Brown's badly scaled function.
       */
      Brown() {
        this->add_solution(this->make({1.0e6, 2.0e-6}));
        this->add_guess(this->make({1.0, 1.0}));
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Brown";
      }

      void residual(const Dense &x, Dense &F) const {
        F(0) = x(0) - Scalar(1.0e6);
        F(1) = x(1) - Scalar(2.0e-6);
        F(2) = x(0) * x(1) - 2;
      }

      void jacobian(const Dense &x, DenseMatrix &J) const {
        J << 1, 0, 0, 1, x(1), x(0);
      }

      void hessian([[maybe_unused]] const Dense &x, const Integer i, DenseMatrix &H) const {
        H.setZero();
        if (i == 2) {
          H(0, 1) = 1;
          H(1, 0) = 1;
        }
      }

    };  // class Brown

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_BROWN_HH
