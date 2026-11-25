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

#include "Optimist/TestSet/DenseFunction.hh"

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
     * The function has a global minimum at \f$\mathbf{x} = (0, 0)\f$, with
     * \f$f(\mathbf{x}) = 0\f$, surrounded by concentric rings of local
     * minima. Local methods started far from the origin are expected to stop
     * at one of these rings, so the function is meant for global optimizers.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class Schaffer2 : public DenseCostFunction<Vector, 2, Schaffer2<Vector>> {
     public:
      using Base = DenseCostFunction<Vector, 2, Schaffer2<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for the Schaffer2 function.
       */
      Schaffer2() {
        this->add_solution(this->make({0.0, 0.0}));
        this->add_guess(this->make({-50.0, -50.0}));
        this->add_guess(this->make({50.0, 50.0}));
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Schaffer2";
      }

      /**
       * Usual half-width of the search domain.
       */
      static constexpr Scalar domain() {
        return 100.0;
      }

      // f = 1/2 + A * B, with A = sin^2(u) - 1/2, u = x^2 - y^2, and
      // B = D^(-2), D = 1 + (x^2 + y^2)/1000

      Scalar value(const Dense &x) const {
        const Scalar u{x(0) * x(0) - x(1) * x(1)};
        const Scalar D{1 + (x(0) * x(0) + x(1) * x(1)) / 1000};
        const Scalar s{std::sin(u)};
        return Scalar(0.5) + (s * s - Scalar(0.5)) / (D * D);
      }

      void gradient(const Dense &x, Dense &g) const {
        Scalar A, B;
        Dense dA(2), dB(2);
        DenseMatrix ddA(2, 2), ddB(2, 2);
        this->terms(x, A, B, dA, dB, ddA, ddB);
        g = dA * B + A * dB;
      }

      void hessian(const Dense &x, DenseMatrix &H) const {
        Scalar A, B;
        Dense dA(2), dB(2);
        DenseMatrix ddA(2, 2), ddB(2, 2);
        this->terms(x, A, B, dA, dB, ddA, ddB);
        H = ddA * B + dA * dB.transpose() + dB * dA.transpose() + A * ddB;
      }

     private:
      void terms(const Dense &x,
                 Scalar &A,
                 Scalar &B,
                 Dense &dA,
                 Dense &dB,
                 DenseMatrix &ddA,
                 DenseMatrix &ddB) const {
        const Scalar u{x(0) * x(0) - x(1) * x(1)};
        const Scalar D{1 + (x(0) * x(0) + x(1) * x(1)) / 1000};
        const Scalar s{std::sin(u)}, s2{std::sin(2 * u)}, c2{std::cos(2 * u)};
        A = s * s - Scalar(0.5);
        B = 1 / (D * D);
        // dA = sin(2u) * du, du = (2x, -2y)
        dA(0) = 2 * x(0) * s2;
        dA(1) = -2 * x(1) * s2;
        ddA(0, 0) = 2 * s2 + 8 * x(0) * x(0) * c2;
        ddA(1, 1) = -2 * s2 + 8 * x(1) * x(1) * c2;
        ddA(0, 1) = -8 * x(0) * x(1) * c2;
        ddA(1, 0) = ddA(0, 1);
        // dB = -2 D^(-3) dD, dD = (x, y)/500
        const Scalar D3{D * D * D}, D4{D3 * D};
        dB(0) = -x(0) / (250 * D3);
        dB(1) = -x(1) / (250 * D3);
        ddB(0, 0) = -1 / (250 * D3) + 3 * x(0) * x(0) / (125000 * D4);
        ddB(1, 1) = -1 / (250 * D3) + 3 * x(1) * x(1) / (125000 * D4);
        ddB(0, 1) = 3 * x(0) * x(1) / (125000 * D4);
        ddB(1, 0) = ddB(0, 1);
      }

    };  // class Schaffer2

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_SCHAFFER2_HH
