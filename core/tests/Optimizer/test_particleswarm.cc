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
#include "Optimist/Optimizer/ParticleSwarm.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

namespace {

  // Configure a solver on the symmetric box [-h, h]^n
  template <typename Solver, typename Function>
  void setup(Solver &sol, const Function &fun, const double h) {
    using Vector = typename Function::VectorTrait::Type;
    using Scalar = typename Function::Scalar;
    const Integer n{static_cast<Integer>(fun.guess(0).size())};
    sol.bounds(constant<Vector>(n, static_cast<Scalar>(-h)),
               constant<Vector>(n, static_cast<Scalar>(h)));
    sol.task(fun.name());
    sol.tolerance(solver_tolerance<Scalar>());
  sol.particles(30);
  sol.max_iterations(5000);
  sol.max_function_evaluations(200000);
  }

}  // namespace

// Multimodal functions: local methods stop at the nearest local minimum
using TestTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Rastrigin2, 2>,
                          VectorVariants<Ackley2, 2>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, GlobalMinimum) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  Function fun;
  Optimist::Optimizer::ParticleSwarm<Vector> sol;
  setup(sol, fun, 5.0);
  expect_minimizes_all(sol, fun);
}

// The minimum of Schaffer2 is extremely flat (quartic along x = y): the gap
// test is strict, the distance test only resolves the minimizer to about the
// fourth root of the rounding level, so it is tested in double precision only
using SchafferTypes = ConcatT<DoubleVectorVariants<Schaffer2, 2>>;

OPTIMIST_TYPED_TEST_SUITE(Schaffer, SchafferTypes);

TYPED_TEST(Schaffer, GlobalMinimum) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  Function fun;
  Optimist::Optimizer::ParticleSwarm<Vector> sol;
  setup(sol, fun, Function::domain());
  expect_minimizes_all(sol, fun);
}

TEST(ParticleSwarm, SameSeedSameResult) {
  Ackley2<Eigen::Vector2d> fun;
  Eigen::Vector2d x_1, x_2;
  for (Eigen::Vector2d *x : {&x_1, &x_2}) {
    Optimist::Optimizer::ParticleSwarm<Eigen::Vector2d> sol;
    setup(sol, fun, 5.0);
    sol.seed(7);
    sol.optimize(fun, fun.guess(0), *x);
  }
  EXPECT_EQ(x_1, x_2);
}

TEST(ParticleSwarm, RespectsBounds) {
  // The unconstrained minimum (1, 1) lies outside the box [-0.5, 0.5]^2
  RosenbrockCost2<Eigen::Vector2d> fun;
  Optimist::Optimizer::ParticleSwarm<Eigen::Vector2d> sol;
  setup(sol, fun, 0.5);
  Eigen::Vector2d x_out;
  sol.optimize(fun, Eigen::Vector2d(0.0, 0.0), x_out);
  EXPECT_LE(x_out.cwiseAbs().maxCoeff(), 0.5);
}

TEST(ParticleSwarm, InertiaWeightStrategies) {
  using Vector = Eigen::Vector2d;
  using namespace Optimist::Optimizer;
  EllipticParaboloid<Vector> fun;
  auto check = [&fun](auto &sol) {
    setup(sol, fun, 10.0);
    EXPECT_TRUE(minimizes(sol, fun, fun.guess(0)));
  };
  ParticleSwarm<Vector, ConstantWeight<double>> s_1;
  ParticleSwarm<Vector, LinearDecrease<double>> s_2;
  ParticleSwarm<Vector, ExponentialDecrease1<double>> s_3;
  ParticleSwarm<Vector, ExponentialDecrease2<double>> s_4;
  ParticleSwarm<Vector, ExponentialDecrease3<double>> s_5;
  check(s_1);
  check(s_2);
  check(s_3);
  check(s_4);
  check(s_5);
}

TEST(ParticleSwarm, LinearDecreaseIsDecreasing) {
  Optimist::Optimizer::LinearDecrease<double> w;
  EXPECT_DOUBLE_EQ(w(0, 100), 0.9);
  EXPECT_DOUBLE_EQ(w(100, 100), 0.4);
  EXPECT_GT(w(10, 100), w(90, 100));
}
