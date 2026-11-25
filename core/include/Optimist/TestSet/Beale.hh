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

#ifndef OPTIMIST_TESTSET_BEALE_HH
#define OPTIMIST_TESTSET_BEALE_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Beale's function.
     *
     * \f[
     *   f(x, y) = (1.5 - x + xy)^2 + (2.25 - x + xy^2)^2 + (2.625 - x + xy^3)^2
     *   \text{,}
     * \f]
     * with global minimum \f$ f(3, 0.5) = 0 \f$ at the bottom of a narrow
     * curved valley surrounded by very steep walls.
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class Beale : public DenseCostFunction<Vector, 2, Beale<Vector>> {
     public:
      using Base = DenseCostFunction<Vector, 2, Beale<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for Beale's function.
       */
      Beale() {
        this->add_solution(this->make({3.0, 0.5}));
        this->add_guess(this->make({1.0, 1.0}));
        this->add_guess(this->make({0.0, 0.0}));
        this->add_guess(this->make({2.0, 0.2}));
      }

      constexpr std::string name_impl() const {
        return "Beale";
      }

      Scalar value(const Dense &x) const {
        Scalar f{0.0};
        for (Integer i{1}; i <= 3; ++i) {
          const Scalar t{c(i) - x(0) + x(0) * std::pow(x(1), Scalar(i))};
          f += t * t;
        }
        return f;
      }

      void gradient(const Dense &x, Dense &g) const {
        g.setZero();
        for (Integer i{1}; i <= 3; ++i) {
          const Scalar yi{std::pow(x(1), Scalar(i))};
          const Scalar t{c(i) - x(0) + x(0) * yi};
          g(0) += 2 * t * (yi - 1);
          g(1) += 2 * t * (i * x(0) * std::pow(x(1), Scalar(i - 1)));
        }
      }

      void hessian(const Dense &x, DenseMatrix &H) const {
        H.setZero();
        for (Integer i{1}; i <= 3; ++i) {
          const Scalar si{static_cast<Scalar>(i)};
          const Scalar yi{std::pow(x(1), si)};
          const Scalar yi1{std::pow(x(1), si - 1)};
          const Scalar yi2{i >= 2 ? std::pow(x(1), si - 2) : Scalar(0)};
          const Scalar t{c(i) - x(0) + x(0) * yi};
          const Scalar tx{yi - 1}, ty{si * x(0) * yi1};
          const Scalar txy{si * yi1}, tyy{si * (si - 1) * x(0) * yi2};
          H(0, 0) += 2 * tx * tx;
          H(0, 1) += 2 * (tx * ty + t * txy);
          H(1, 1) += 2 * (ty * ty + t * tyy);
        }
        H(1, 0) = H(0, 1);
      }

     private:
      static Scalar c(const Integer i) {
        constexpr Scalar values[3]{1.5, 2.25, 2.625};
        return values[i - 1];
      }

    };  // class Beale

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_BEALE_HH
