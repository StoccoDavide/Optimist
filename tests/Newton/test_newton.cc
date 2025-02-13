/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Optimist.hh"
#include "Optimist/TestSet.hh"
#include <gtest/gtest.h>

using namespace Optimist;
using namespace Optimist::RootFinder;
using namespace Optimist::TestSet;

// Funzione di Booth.
TEST(Newton, Booth) {
  // Function
  Booth nlfunction;
  // Non-linear solver
  Newton<2> nlsolver;
  nlsolver.task(nlfunction.name());
  // Starting entries
  Vector2 x_ini = Vector2::Zero();
  Vector2 x_out = Vector2::Zero();
  auto fun = [&nlfunction](Vector2 const & X, Vector2 & F) {nlfunction.evaluate(X, F);};
  auto jac = [&nlfunction](Vector2 const & X, Matrix2 & J) {nlfunction.jacobian(X, J);};
  for (Integer i{0}; i < static_cast<Integer>(nlfunction.guesses().size()); ++i) {
    x_ini = nlfunction.guess(i);
    // Solve without damping
    nlsolver.disable_damped_mode();
    nlsolver.rootfind(nlfunction, x_ini, x_out);
    nlsolver.solve(fun, jac, x_ini, x_out);
    EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
    EXPECT_TRUE(nlsolver.converged());
    EXPECT_TRUE(nlfunction.is_solution(x_out));
    // Solve with damping
    nlsolver.enable_damped_mode();
    nlsolver.solve(fun, jac, x_ini, x_out);
    EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
    EXPECT_TRUE(nlsolver.converged());
    EXPECT_TRUE(nlfunction.is_solution(x_out));
  }
}

// 2D Rosenbrock function
TEST(Newton, Rosenbrock2D) {
  // Non-linear solver
  Newton<2> nlsolver;
  nlsolver.task("Rosenbrock2D");
  for (Real a = 1.0; a <= 5.0; a += 1.0) {
    for (Real b = 1.0; b <= 5.0; b += 1.0) {
      // Starting entries
      Vector2 x_ini = Vector2::Zero();
      Vector2 x_out = Vector2::Zero();
      auto fun = [a, b](Vector2 const & X, Vector2 & F) {
        F << a*(1 - X(0)), b*(X(1) - X(0)*X(0));
      };
      auto jac = [a, b](Vector2 const & /*X*/, Matrix2 & JF) {
        JF << -a, Real(0.0), Real(0.0), b;
      };
      // Solve without damping
      nlsolver.disable_damped_mode();
      nlsolver.solve(fun, jac, x_ini, x_out);
      EXPECT_LE((x_out - Vector2::Ones()).maxCoeff(), EPSILON_LOW);
      EXPECT_TRUE(nlsolver.converged());
      // Solve with damping
      nlsolver.enable_damped_mode();
      nlsolver.solve(fun, jac, x_ini, x_out);
      EXPECT_LE((x_out - Vector2::Ones()).maxCoeff(), EPSILON_LOW);
      EXPECT_TRUE(nlsolver.converged());
    }
  }
}

// 3D Rosenbrock function
TEST(Newton, Rosenbrock3D) {
  // Non-linear solver
  Newton<3> nlsolver;
  nlsolver.task("Rosenbrock3D");
  for (Real a = 1.0; a <= 10.0; a += 1.0) {
    for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector3 x_ini = Vector3::Zero();
        Vector3 x_out = Vector3::Zero();
        auto fun = [a, b](Vector3 const & X, Vector3 & F) {
          F << a*(1 - X(0)), b*(X(1) - X(0)*X(0)), b*(X(2) - X(1)*X(1));
        };
        auto jac = [a, b](Vector3 const & X, Matrix3 & JF) {
          JF << -a, Real(0.0), Real(0.0), -2*b*X(0), b, Real(0.0), Real(0.0), -2*b*X(1), b;
        };
        // Solve without damping
        nlsolver.disable_damped_mode();
        nlsolver.solve(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - Vector3::Ones()).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
        // Solve with damping
        nlsolver.enable_damped_mode();
        nlsolver.solve(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - Vector3::Ones()).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
    }
  }
}

template <typename F, std::size_t... Is>
void repeat_unrolled(F&& f, std::index_sequence<Is...>) {
  (f(std::integral_constant<std::size_t, Is>{}), ...);
}

template <std::size_t... Iterations, typename F>
void for_unrolled(F&& f) {
  repeat_unrolled(std::forward<F>(f), std::index_sequence<Iterations...>{});
}

// ND Rosenbrock function
TEST(Newton, RosenbrockND) {
  for_unrolled<2, 3>([](auto D) {
    using Vector = typename Newton<D>::Vector;
    using Matrix = typename Newton<D>::Matrix;
    // Non-linear solver
    Newton<D> nlsolver;
    nlsolver.task("RosenbrockND(" + std::to_string(D) + ")");
    for (Real a = 1.0; a <= 10.0; a += 1.0) {
      for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector x_ini = Vector::Zero();
        Vector x_out = Vector::Zero();
        auto fun = [a, b, D](Vector const & X, Vector & F) {
          F(0) = a*(1 - X(0));
          for (std::size_t i = 1; i < D; ++i) {
            F(i) = b*(X(i) - X(i-1)*X(i-1));
          }
        };
        auto jac = [a, b, D](Vector const & X, Matrix & JF) {
          JF.setZero();
          JF(0, 0) = -a;
          for (std::size_t i = 1; i < D; ++i) {
            JF(i, i)   = b;
            JF(i, i-1) = -2*b*X(i-1);
          }
        };
        // Solve without damping
        nlsolver.disable_damped_mode();
        nlsolver.solve(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - Vector::Ones()).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
        // Solve with damping
        nlsolver.enable_damped_mode();
        nlsolver.solve(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - Vector::Ones()).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
      }
    }
  });
}

// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
