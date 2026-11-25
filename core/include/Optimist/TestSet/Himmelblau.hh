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

#ifndef OPTIMIST_TESTSET_HIMMELBLAU_HH
#define OPTIMIST_TESTSET_HIMMELBLAU_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Himmelblau's function.
     *
     * \f[
     *   f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2 \text{,}
     * \f]
     * which has four identical local minima \f$ f = 0 \f$ (and therefore
     * tests whether a solver lands on \em one of several solutions), a local
     * maximum and four saddle points.
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class Himmelblau
        : public DenseCostFunction<Vector, 2, Himmelblau<Vector>> {
     public:
      using Base = DenseCostFunction<Vector, 2, Himmelblau<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for Himmelblau's function.
       */
      Himmelblau() {
        this->add_solution(this->make({3.0, 2.0}));
        this->add_solution(this->make({-2.805118086952745, 3.131312518250573}));
        this->add_solution(
            this->make({-3.779310253377747, -3.283185991286170}));
        this->add_solution(this->make({3.584428340330492, -1.848126526964404}));
        this->add_guess(this->make({2.0, 2.0}));
        this->add_guess(this->make({-2.0, 2.0}));
        this->add_guess(this->make({-2.0, -2.0}));
        this->add_guess(this->make({2.0, -2.0}));
      }

      constexpr std::string name_impl() const {
        return "Himmelblau";
      }

      Scalar value(const Dense &x) const {
        const Scalar a{x(0) * x(0) + x(1) - 11}, b{x(0) + x(1) * x(1) - 7};
        return a * a + b * b;
      }

      void gradient(const Dense &x, Dense &g) const {
        const Scalar a{x(0) * x(0) + x(1) - 11}, b{x(0) + x(1) * x(1) - 7};
        g(0) = 4 * x(0) * a + 2 * b;
        g(1) = 2 * a + 4 * x(1) * b;
      }

      void hessian(const Dense &x, DenseMatrix &H) const {
        H(0, 0) = 12 * x(0) * x(0) + 4 * x(1) - 42;
        H(0, 1) = 4 * (x(0) + x(1));
        H(1, 0) = H(0, 1);
        H(1, 1) = 4 * x(0) + 12 * x(1) * x(1) - 26;
      }

    };  // class Himmelblau

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_HIMMELBLAU_HH
