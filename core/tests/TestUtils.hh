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

#ifndef OPTIMIST_TESTS_TESTUTILS_HH
#define OPTIMIST_TESTS_TESTUTILS_HH

// Optimist library
#include "Optimist.hh"

// Google Test framework
#include <gtest/gtest.h>

// C++ standard libraries
#include <initializer_list>
#include <sstream>

/**
 * \brief Declare a type-parameterized test fixture and register its type list.
 * \param NAME Fixture name.
 * \param TYPES A \c testing::Types<...> list.
 */
#define OPTIMIST_TYPED_TEST_SUITE(NAME, TYPES) \
  template <typename T>                        \
  struct NAME : public testing::Test {};       \
  TYPED_TEST_SUITE(NAME, TYPES)

namespace OptimistTest {

  using Optimist::Integer;
  using Optimist::TypeTrait;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Type lists
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * \brief Concatenate several \c testing::Types<...> lists into one.
   */
  template <typename... Lists>
  struct Concat;

  template <typename... A>
  struct Concat<testing::Types<A...>> {
    using type = testing::Types<A...>;
  };

  template <typename... A, typename... B, typename... Rest>
  struct Concat<testing::Types<A...>, testing::Types<B...>, Rest...> {
    using type = typename Concat<testing::Types<A..., B...>, Rest...>::type;
  };

  template <typename... Lists>
  using ConcatT = typename Concat<Lists...>::type;

  /**
   * \brief Instantiate a scalar function family for \c float and \c double.
   */
  template <template <typename> class F>
  using ScalarVariants = testing::Types<F<float>, F<double>>;

  /**
   * \brief Instantiate a vector function family for fixed-size, dynamic-size
   * and sparse vectors, both in single and double precision.
   */
  template <template <typename> class F, int N>
  using VectorVariants = testing::Types<F<Eigen::Vector<float, N>>,
                                        F<Eigen::Vector<float, Eigen::Dynamic>>,
                                        F<Eigen::SparseVector<float>>,
                                        F<Eigen::Vector<double, N>>,
                                        F<Eigen::Vector<double, Eigen::Dynamic>>,
                                        F<Eigen::SparseVector<double>>>;

  /**
   * \brief Same as VectorVariants but without sparse vectors.
   */
  template <template <typename> class F, int N>
  using DenseVectorVariants =
      testing::Types<F<Eigen::Vector<float, N>>,
                     F<Eigen::Vector<float, Eigen::Dynamic>>,
                     F<Eigen::Vector<double, N>>,
                     F<Eigen::Vector<double, Eigen::Dynamic>>>;

  /**
   * \brief Same as DenseVectorVariants but in double precision only.
   */
  template <template <typename> class F, int N>
  using DoubleVectorVariants =
      testing::Types<F<Eigen::Vector<double, N>>,
                     F<Eigen::Vector<double, Eigen::Dynamic>>,
                     F<Eigen::SparseVector<double>>>;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Vector helpers
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * \brief Build a vector of any supported type from a list of values.
   * \tparam Vector Fixed-size, dynamic-size or sparse Eigen vector, or scalar.
   * \param values Coefficients.
   * \return The vector.
   */
  template <typename Vector>
  Vector make(std::initializer_list<typename TypeTrait<Vector>::Scalar> values) {
    using Trait = TypeTrait<Vector>;
    if constexpr (Trait::IsScalar) {
      return *values.begin();
    } else {
      const Integer n{static_cast<Integer>(values.size())};
      Vector out;
      if constexpr (Trait::IsFixed) {
        EXPECT_EQ(n, Trait::Dimension);
      } else {
        out.resize(n);
      }
      if constexpr (Trait::IsSparse) {
        out.reserve(n);
      }
      Integer i{0};
      for (const auto v : values) {
        if constexpr (Trait::IsSparse) {
          out.coeffRef(i++) = v;
        } else {
          out(i++) = v;
        }
      }
      return out;
    }
  }

  /**
   * \brief Build a constant vector of any supported type.
   */
  template <typename Vector>
  Vector constant(const Integer n, const typename TypeTrait<Vector>::Scalar v) {
    using Trait = TypeTrait<Vector>;
    if constexpr (Trait::IsScalar) {
      return v;
    } else {
      Vector out;
      if constexpr (!Trait::IsFixed) {
        out.resize(n);
      }
      if constexpr (Trait::IsSparse) {
        out.reserve(n);
        for (Integer i{0}; i < n; ++i) {
          out.coeffRef(i) = v;
        }
      } else {
        out.setConstant(v);
      }
      return out;
    }
  }

  /**
   * \brief Pretty-print any supported input type.
   */
  template <typename Vector>
  std::string to_string(const Vector &x) {
    std::ostringstream os;
    os << std::setprecision(8);
    if constexpr (TypeTrait<Vector>::IsScalar) {
      os << x;
    } else {
      using Dense = Eigen::Matrix<typename TypeTrait<Vector>::Scalar,
                                  Eigen::Dynamic,
                                  1>;
      Dense d(x);
      os << "[";
      for (Eigen::Index i{0}; i < d.size(); ++i) {
        os << (i ? ", " : "") << d(i);
      }
      os << "]";
    }
    return os.str();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Solution checks
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * \brief Problem type for the solution checks.
   */
  enum class Problem { ROOTFIND, OPTIMIZE };

  /**
   * \brief Convergence tolerance given to the solvers: 1e-12 in double
   * precision, \f$ \sqrt{\epsilon} \f$ in single precision (where 1e-12 is
   * far below the machine epsilon).
   */
  template <typename Scalar>
  Scalar solver_tolerance() {
    if constexpr (std::is_same_v<Scalar, float>) {
      return std::sqrt(std::numeric_limits<Scalar>::epsilon());
    } else {
      return static_cast<Scalar>(1.0e-12);
    }
  }

  /**
   * \brief Acceptance tolerance on the distance from a known solution: 1e-10
   * in double precision, \f$ \sqrt[3]{\epsilon} \f$ in single precision.
   */
  template <typename Scalar>
  Scalar default_tolerance() {
    if constexpr (std::is_same_v<Scalar, float>) {
      return std::cbrt(std::numeric_limits<Scalar>::epsilon());
    } else {
      return static_cast<Scalar>(1.0e-10);
    }
  }

  /**
   * \brief Tolerance on the optimality gap \f$ f(\mathbf{x}) - f^\star \f$
   * (relative to \f$ 1 + |f^\star| \f$) for methods that only compare
   * function values: 1e-12 in double precision, \f$ \sqrt{\epsilon} \f$ in
   * single precision.
   */
  template <typename Scalar>
  Scalar value_gap_tolerance() {
    return solver_tolerance<Scalar>();
  }

  /**
   * \brief Acceptance tolerance on the distance from the minimizer for methods
   * that only compare function values. Near a minimizer \f$ f(\mathbf{x}^\star
   * + \boldsymbol{\delta}) - f^\star \approx \frac{1}{2}\boldsymbol{\delta}^\top
   * \nabla^2 f\boldsymbol{\delta} \f$, so a gap at rounding level only pins
   * the minimizer down to about its square root: 1e-5 = sqrt(1e-10) in double
   * precision, \f$ \sqrt[3]{\epsilon} \f$ in single precision.
   */
  template <typename Scalar>
  Scalar value_tolerance() {
    if constexpr (std::is_same_v<Scalar, float>) {
      return std::cbrt(std::numeric_limits<Scalar>::epsilon());
    } else {
      return static_cast<Scalar>(1.0e-5);
    }
  }

  /**
   * \brief Check the bookkeeping that every solver must respect after a solve.
   */
  template <typename Solver>
  testing::AssertionResult counters_are_consistent(const Solver &sol) {
    if (sol.iterations() < 0 || sol.iterations() > sol.max_iterations()) {
      return testing::AssertionFailure()
             << sol.name() << ": iterations " << sol.iterations()
             << " outside [0, " << sol.max_iterations() << "]";
    }
    if (sol.function_evaluations() <= 0 ||
        sol.function_evaluations() > sol.max_function_evaluations()) {
      return testing::AssertionFailure()
             << sol.name() << ": function evaluations "
             << sol.function_evaluations() << " outside (0, "
             << sol.max_function_evaluations() << "]";
    }
    return testing::AssertionSuccess();
  }

  /**
   * \brief Solve a problem from a given guess and check that the solver
   * converged to one of the known solutions of the function.
   * \param sol Solver instance (already configured).
   * \param fun Test function.
   * \param x_ini Initial guess.
   * \tparam P Root-finding or optimization problem.
   * \param tol Acceptance tolerance on the distance from a known solution.
   * \param x_out Optional output solution.
   */
  template <Problem P, typename Solver, typename Function, typename Input>
  testing::AssertionResult solves(
      Solver &sol,
      const Function &fun,
      const Input &x_ini,
      const typename TypeTrait<Input>::Scalar tol,
      Input *x_out = nullptr) {
    Input x_sol{x_ini};
    try {
      if constexpr (P == Problem::ROOTFIND) {
        sol.rootfind(fun, x_ini, x_sol);
      } else {
        sol.optimize(fun, x_ini, x_sol);
      }
    } catch (const std::exception &e) {
      return testing::AssertionFailure()
             << sol.name() << " on " << fun.name() << " from "
             << to_string(x_ini) << " threw: " << e.what();
    }
    if (x_out != nullptr) {
      *x_out = x_sol;
    }
    auto counters{counters_are_consistent(sol)};
    if (!counters) {
      return counters;
    }
    if (!sol.converged()) {
      return testing::AssertionFailure()
             << sol.name() << " on " << fun.name() << " from "
             << to_string(x_ini) << " did not converge after "
             << sol.iterations() << " iterations (stopped at "
             << to_string(x_sol) << ")";
    }
    if (!fun.is_solution(x_sol, tol)) {
      return testing::AssertionFailure()
             << sol.name() << " on " << fun.name() << " from "
             << to_string(x_ini) << " converged to " << to_string(x_sol)
             << ", which is not a known solution (tol = " << tol << ")";
    }
    return testing::AssertionSuccess();
  }

  /**
   * \brief Minimize from a given guess with a method that only compares
   * function values, and check both the optimality gap with respect to the
   * closest known solution and the distance from it.
   */
  template <typename Solver, typename Function, typename Input>
  testing::AssertionResult minimizes(Solver &sol,
                                     const Function &fun,
                                     const Input &x_ini,
                                     Input *x_out = nullptr) {
    using Scalar = typename TypeTrait<Input>::Scalar;
    Input x_sol{x_ini};
    auto result{solves<Problem::OPTIMIZE>(sol, fun, x_ini,
                                          value_tolerance<Scalar>(), &x_sol)};
    if (x_out != nullptr) {
      *x_out = x_sol;
    }
    if (!result) {
      return result;
    }
    // Closest known solution
    const Input *closest{nullptr};
    Scalar distance{std::numeric_limits<Scalar>::infinity()};
    for (const Input &s : fun.solutions()) {
      Scalar d;
      if constexpr (TypeTrait<Input>::IsScalar) {
        d = std::abs(x_sol - s);
      } else {
        d = (x_sol - s).norm();
      }
      if (d < distance) {
        distance = d;
        closest  = &s;
      }
    }
    Scalar f_sol, f_star;
    if (!fun.evaluate(x_sol, f_sol) || !fun.evaluate(*closest, f_star)) {
      return testing::AssertionFailure() << "evaluation failed";
    }
    const Scalar gap{f_sol - f_star};
    if (gap > value_gap_tolerance<Scalar>() * (1 + std::abs(f_star))) {
      return testing::AssertionFailure()
             << sol.name() << " on " << fun.name() << " from "
             << to_string(x_ini) << ": optimality gap " << gap << " at "
             << to_string(x_sol);
    }
    return testing::AssertionSuccess();
  }

  /**
   * \brief Run minimizes(...) from every suggested guess of the function.
   */
  template <typename Solver, typename Function>
  void expect_minimizes_all(Solver &sol, const Function &fun) {
    for (std::size_t i{0}; i < fun.guesses().size(); ++i) {
      EXPECT_TRUE(minimizes(sol, fun, fun.guess(i))) << "(guess #" << i << ")";
    }
  }

  /**
   * \brief Run solves(...) from every suggested guess of the function.
   */
  template <Problem P, typename Solver, typename Function>
  void expect_solves_all(Solver &sol,
                         const Function &fun,
                         const typename Function::Scalar tol) {
    for (std::size_t i{0}; i < fun.guesses().size(); ++i) {
      EXPECT_TRUE(solves<P>(sol, fun, fun.guess(i), tol))
          << "(guess #" << i << ")";
    }
  }

}  // namespace OptimistTest

#endif  // OPTIMIST_TESTS_TESTUTILS_HH
