/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimist library
#include "Optimist/RootFinder/Newton.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<VectorVariants<Linear1, 1>,
                          VectorVariants<Booth, 2>,
                          VectorVariants<Rosenbrock2, 2>,
                          VectorVariants<Rosenbrock4, 4>,
                          VectorVariants<Rosenbrock6, 6>,
                          VectorVariants<Rosenbrock8, 8>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, RootFind) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::Newton<Vector> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());

  for (const bool damped : {false, true}) {
    sol.damped_mode(damped);
    expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
  }

}

// Root finders handed an optimization problem solve grad(phi) = 0
using CostTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Trid2, 2>,
                          DoubleVectorVariants<Himmelblau, 2>>;

OPTIMIST_TYPED_TEST_SUITE(CostFunctions, CostTypes);

TYPED_TEST(CostFunctions, StationaryPoint) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::Newton<Vector> sol;
  sol.tolerance(solver_tolerance<Scalar>());
  for (std::size_t i{0}; i < fun.guesses().size(); ++i) {
    Vector x_out;
    sol.optimize(fun, fun.guess(i), x_out);
    EXPECT_TRUE(sol.converged()) << fun.name() << " guess #" << i;
    typename Function::FirstDerivative g;
    ASSERT_TRUE(fun.first_derivative(x_out, g));
    EXPECT_LT(g.norm(), 10 * std::sqrt(Function::EPSILON))
        << fun.name() << " guess #" << i;
    static_cast<void>(default_tolerance<Scalar>);
  }
}
