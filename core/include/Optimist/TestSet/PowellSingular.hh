/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_POWELL_SINGULAR_HH
#define OPTIMIST_TESTSET_POWELL_SINGULAR_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Powell's singular function.
     *
     * \f[
     *   \mathbf{F}(\mathbf{x}) = \begin{bmatrix}
     *     x_1 + 10x_2 \\
     *     \sqrt{5}(x_3 - x_4) \\
     *     (x_2 - 2x_3)^2 \\
     *     \sqrt{10}(x_1 - x_4)^2
     *   \end{bmatrix} \text{,}
     * \f]
     * whose root \f$ \mathbf{x} = \mathbf{0} \f$ has a singular Jacobian, so
     * that Newton-like methods only converge linearly. It is a standard
     * robustness test (Moré, Garbow and Hillstrom, problem 13).
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class PowellSingular
        : public DenseSystemFunction<Vector, 4, PowellSingular<Vector>> {
     public:
      using Base = DenseSystemFunction<Vector, 4, PowellSingular<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for Powell's singular function.
       */
      PowellSingular() {
        this->add_solution(Dense::Zero(4));
        this->add_guess(this->make({3.0, -1.0, 0.0, 1.0}));
      }

      constexpr std::string name_impl() const {
        return "PowellSingular";
      }

      void residual(const Dense &x, Dense &F) const {
        const Scalar s5{std::sqrt(Scalar(5))}, s10{std::sqrt(Scalar(10))};
        F(0) = x(0) + 10 * x(1);
        F(1) = s5 * (x(2) - x(3));
        F(2) = (x(1) - 2 * x(2)) * (x(1) - 2 * x(2));
        F(3) = s10 * (x(0) - x(3)) * (x(0) - x(3));
      }

      void jacobian(const Dense &x, DenseMatrix &J) const {
        const Scalar s5{std::sqrt(Scalar(5))}, s10{std::sqrt(Scalar(10))};
        const Scalar a{x(1) - 2 * x(2)}, b{x(0) - x(3)};
        J.setZero();
        J(0, 0) = 1;
        J(0, 1) = 10;
        J(1, 2) = s5;
        J(1, 3) = -s5;
        J(2, 1) = 2 * a;
        J(2, 2) = -4 * a;
        J(3, 0) = 2 * s10 * b;
        J(3, 3) = -2 * s10 * b;
      }

      void hessian([[maybe_unused]] const Dense &x, const Integer i, DenseMatrix &H) const {
        const Scalar s10{std::sqrt(Scalar(10))};
        H.setZero();
        if (i == 2) {
          H(1, 1) = 2;
          H(1, 2) = -4;
          H(2, 1) = -4;
          H(2, 2) = 8;
        } else if (i == 3) {
          H(0, 0) = 2 * s10;
          H(0, 3) = -2 * s10;
          H(3, 0) = -2 * s10;
          H(3, 3) = 2 * s10;
        }
      }

    };  // class PowellSingular

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_POWELL_SINGULAR_HH
