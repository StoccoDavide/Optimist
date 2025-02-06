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
#include <gtest/gtest.h>

using namespace Optimist;
using namespace Optimist::RootFinder;

template <typename F, std::size_t... Is>
void repeat_unrolled(F&& f, std::index_sequence<Is...>) {
  (f(std::integral_constant<std::size_t, Is>{}), ...);
}

template <std::size_t... Iterations, typename F>
void for_unrolled(F&& f) {
  repeat_unrolled(std::forward<F>(f), std::index_sequence<Iterations...>{});
}

/*\
 |   ____            _
 |  | __ )  __ _  __| |
 |  |  _ \ / _` |/ _` |
 |  | |_) | (_| | (_| |
 |  |____/ \__,_|\__,_|
 |
\*/

// Funzione di Booth.
TEST(BroydenBad, Booth) {
  // Non-linear solver
  Broyden<2> nlsolver;
  nlsolver.task("Booth");
  nlsolver.enable_bad_method();
  // Starting entries
  Vector2 x_ini = Vector2::Zero();
  Vector2 x_out = Vector2::Zero();
  std::function<void(const Vector2&, Vector2&)> fun = [](const Vector2 & X, Vector2 & F) {
    F << X(0) + 2.0*X(1) - 7.0, 2.0*X(0) + X(1) - 5.0;
  };
  std::function<void(const Vector2&, Matrix2&)> jac = [](const Vector2 & /*X*/, Matrix2 & JF) {
    JF << 1.0, 2.0, 2.0, 1.0;
  };
  // Solve without damping
  nlsolver.disable_damped_mode();
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
  // Solve with damping
  nlsolver.enable_damped_mode();
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
}

// 2D Rosenbrock function
TEST(BroydenBad, Rosenbrock2D) {
  // Non-linear solver
  Broyden<2> nlsolver;
  nlsolver.task("Rosenbrock2D");
  nlsolver.enable_bad_method();
  for (Real a = 1.0; a <= 5.0; a += 1.0) {
    for (Real b = 1.0; b <= 5.0; b += 1.0) {
      // Starting entries
      Vector2 x_ini = Vector2::Zero();
      Vector2 x_out = Vector2::Zero();
      auto fun = [a, b](const Vector2 & X, Vector2 & F) {
        F << a*(1 - X(0)), b*(X(1) - X(0)*X(0));
      };
      auto jac = [a, b](const Vector2 & /*X*/, Matrix2 & JF) {
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
TEST(BroydenBad, Rosenbrock3D) {
  // Non-linear solver
  Broyden<3> nlsolver;
  nlsolver.task("Rosenbrock3D");
  nlsolver.enable_bad_method();
  for (Real a = 1.0; a <= 10.0; a += 1.0) {
    for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector3 x_ini = Vector3::Zero();
        Vector3 x_out = Vector3::Zero();
        auto fun = [a, b](const Vector3 & X, Vector3 & F) {
          F << a*(1 - X(0)), b*(X(1) - X(0)*X(0)), b*(X(2) - X(1)*X(1));
        };
        auto jac = [a, b](const Vector3 & X, Matrix3 & JF) {
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

// ND Rosenbrock function
TEST(BroydenBad, RosenbrockND) {
  for_unrolled<2, 3>([](auto D) {
    using Vector = typename Broyden<D>::Vector;
    using Matrix = typename Broyden<D>::Matrix;
    // Non-linear solver
    Broyden<D> nlsolver;
    nlsolver.task("RosenbrockND(" + std::to_string(D) + ")");
    nlsolver.enable_bad_method();
    for (Real a = 1.0; a <= 10.0; a += 1.0) {
      for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector x_ini = Vector::Zero();
        Vector x_out = Vector::Zero();
        auto fun = [a, b, D](const Vector & X, Vector & F) {
          F(0) = a*(1 - X(0));
          for (std::size_t i = 1; i < D; ++i) {
            F(i) = b*(X(i) - X(i-1)*X(i-1));
          }
        };
        auto jac = [a, b, D](const Vector & X, Matrix & JF) {
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

/*\
 |    ____                 _
 |   / ___| ___   ___   __| |
 |  | |  _ / _ \ / _ \ / _` |
 |  | |_| | (_) | (_) | (_| |
 |   \____|\___/ \___/ \__,_|
 |
\*/

// Funzione di Booth.
TEST(BroydenGood, Booth) {
  // Non-linear solver
  Broyden<2> nlsolver;
  nlsolver.task("Booth");
  nlsolver.enable_good_method();
  // Starting entries
  Vector2 x_ini = Vector2::Zero();
  Vector2 x_out = Vector2::Zero();
  auto fun = [](const Vector2 & X, Vector2 & F) {
    F << X(0) + 2.0*X(1) - 7.0, 2.0*X(0) + X(1) - 5.0;
  };
  auto jac = [](const Vector2 & /*X*/, Matrix2 & JF) {
    JF << 1.0, 2.0, 2.0, 1.0;
  };
  // Solve without damping
  nlsolver.disable_damped_mode();
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
  // Solve with damping
  nlsolver.enable_damped_mode();
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
}

// 2D Rosenbrock function
TEST(BroydenGood, Rosenbrock2D) {
  // Non-linear solver
  Broyden<2> nlsolver;
  nlsolver.task("Rosenbrock2D");
  nlsolver.enable_good_method();
  for (Real a = 1.0; a <= 5.0; a += 1.0) {
    for (Real b = 1.0; b <= 5.0; b += 1.0) {
      // Starting entries
      Vector2 x_ini = Vector2::Zero();
      Vector2 x_out = Vector2::Zero();
      auto fun = [a, b](const Vector2 & X, Vector2 & F) {
        F << a*(1 - X(0)), b*(X(1) - X(0)*X(0));
      };
      auto jac = [a, b](const Vector2 & /*X*/, Matrix2 & JF) {
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
TEST(BroydenGood, Rosenbrock3D) {
  // Non-linear solver
  Broyden<3> nlsolver;
  nlsolver.task("Rosenbrock3D");
  nlsolver.enable_good_method();
  for (Real a = 1.0; a <= 10.0; a += 1.0) {
    for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector3 x_ini = Vector3::Zero();
        Vector3 x_out = Vector3::Zero();
        auto fun = [a, b](const Vector3 & X, Vector3 & F) {
          F << a*(1 - X(0)), b*(X(1) - X(0)*X(0)), b*(X(2) - X(1)*X(1));
        };
        auto jac = [a, b](const Vector3 & X, Matrix3 & JF) {
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

// ND Rosenbrock function
TEST(BroydenGood, RosenbrockND) {
  for_unrolled<2, 3>([](auto D) {
    using Vector = typename Broyden<D>::Vector;
    using Matrix = typename Broyden<D>::Matrix;
    // Non-linear solver
    Broyden<D> nlsolver;
    nlsolver.task("RosenbrockND(" + std::to_string(D) + ")");
    nlsolver.enable_good_method();
    for (Real a = 1.0; a <= 10.0; a += 1.0) {
      for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector x_ini = Vector::Zero();
        Vector x_out = Vector::Zero();
        auto fun = [a, b, D](const Vector & X, Vector & F) {
          F(0) = a*(1 - X(0));
          for (std::size_t i = 1; i < D; ++i) {
            F(i) = b*(X(i) - X(i-1)*X(i-1));
          }
        };
        auto jac = [a, b, D](const Vector & X, Matrix & JF) {
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

/*\
 |    ____                _     _                _
 |   / ___|___  _ __ ___ | |__ (_)_ __   ___  __| |
 |  | |   / _ \| '_ ` _ \| '_ \| | '_ \ / _ \/ _` |
 |  | |__| (_) | | | | | | |_) | | | | |  __/ (_| |
 |   \____\___/|_| |_| |_|_.__/|_|_| |_|\___|\__,_|
 |
\*/

// Funzione di Booth.
TEST(BroydenCombined, Booth) {
  // Non-linear solver
  Broyden<2> nlsolver;
  nlsolver.task("Booth");
  nlsolver.enable_combined_method();
  // Starting entries
  Vector2 x_ini = Vector2::Zero();
  Vector2 x_out = Vector2::Zero();
  auto fun = [](const Vector2 & X, Vector2 & F) {
    F << X(0) + 2.0*X(1) - 7.0, 2.0*X(0) + X(1) - 5.0;
  };
  auto jac = [](const Vector2 & /*X*/, Matrix2 & JF) {
    JF << 1.0, 2.0, 2.0, 1.0;
  };
  // Solve without damping
  nlsolver.disable_damped_mode();
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
  // Solve with damping
  nlsolver.enable_damped_mode();
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - Vector2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
}

// 2D Rosenbrock function
TEST(BroydenCombined, Rosenbrock2D) {
  // Non-linear solver
  Broyden<2> nlsolver;
  nlsolver.task("Rosenbrock2D");
  nlsolver.enable_combined_method();
  for (Real a = 1.0; a <= 5.0; a += 1.0) {
    for (Real b = 1.0; b <= 5.0; b += 1.0) {
      // Starting entries
      Vector2 x_ini = Vector2::Zero();
      Vector2 x_out = Vector2::Zero();
      auto fun = [a, b](const Vector2 & X, Vector2 & F) {
        F << a*(1 - X(0)), b*(X(1) - X(0)*X(0));
      };
      auto jac = [a, b](const Vector2 & /*X*/, Matrix2 & JF) {
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
TEST(BroydenCombined, Rosenbrock3D) {
  // Non-linear solver
  Broyden<3> nlsolver;
  nlsolver.task("Rosenbrock3D");
  nlsolver.enable_combined_method();
  for (Real a = 1.0; a <= 10.0; a += 1.0) {
    for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector3 x_ini = Vector3::Zero();
        Vector3 x_out = Vector3::Zero();
        auto fun = [a, b](const Vector3 & X, Vector3 & F) {
          F << a*(1 - X(0)), b*(X(1) - X(0)*X(0)), b*(X(2) - X(1)*X(1));
        };
        auto jac = [a, b](const Vector3 & X, Matrix3 & JF) {
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

// ND Rosenbrock function
TEST(BroydenCombined, RosenbrockND) {
  for_unrolled<2, 3>([](auto D) {
    using Vector = typename Broyden<D>::Vector;
    using Matrix = typename Broyden<D>::Matrix;
    // Non-linear solver
    Broyden<D> nlsolver;
    nlsolver.task("RosenbrockND(" + std::to_string(D) + ")");
    nlsolver.enable_combined_method();
    for (Real a = 1.0; a <= 10.0; a += 1.0) {
      for (Real b = 1.0; b <= 10.0; b += 1.0) {
        // Starting entries
        Vector x_ini = Vector::Zero();
        Vector x_out = Vector::Zero();
        auto fun = [a, b, D](const Vector & X, Vector & F) {
          F(0) = a*(1 - X(0));
          for (std::size_t i = 1; i < D; ++i) {
            F(i) = b*(X(i) - X(i-1)*X(i-1));
          }
        };
        auto jac = [a, b, D](const Vector & X, Matrix & JF) {
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
