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
#include "Optimist/Optimizer/ConjugateGradient.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Himmelblau, 2>,
                          VectorVariants<Trid2, 2>,
                          DoubleVectorVariants<Trid6, 6>>;

template <typename Function>
void check_all_methods(const typename Function::Scalar gradient_tolerance) {
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;
  using Solver   = Optimist::Optimizer::ConjugateGradient<Vector>;
  using Alpha    = typename Solver::AlphaMethod;
  using Beta     = typename Solver::BetaMethod;

  Function fun;
  Solver sol;
  sol.task(fun.name());
  sol.tolerance(gradient_tolerance);
  sol.max_iterations(5000);
  sol.max_function_evaluations(50000);
  sol.max_gradient_evaluations(50000);

  for (const auto search : {Optimist::LineSearch::Method::STRONG_WOLFE,
                            Optimist::LineSearch::Method::MORE_THUENTE}) {
  sol.line_search(search);
  for (const auto alpha : {Alpha::ACCEPTED_STEP, Alpha::BARZILAI_BORWEIN}) {
    for (const auto beta : {Beta::FLETCHER_REEVES,
                            Beta::POLAK_RIBIERE,
                            Beta::POLAK_RIBIERE_PLUS,
                            Beta::HESTENES_STIEFEL,
                            Beta::CONJUGATE_DESCENT,
                            Beta::LIU_STOREY,
                            Beta::DAI_YUAN,
                            Beta::HAGER_ZHANG,
                            Beta::HAGER_ZHANG_PLUS}) {
      SCOPED_TRACE(testing::Message() << Optimist::LineSearch::name(search)
                                      << ", alpha method "
                                      << static_cast<int>(alpha)
                                      << ", beta method "
                                      << static_cast<int>(beta));
      sol.alpha_method(alpha);
      sol.beta_method(beta);
      expect_solves_all<Problem::OPTIMIZE>(sol, fun,
                                           default_tolerance<Scalar>());
    }
  }
  }
}

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, LocalMinimum) {
  check_all_methods<TypeParam>(solver_tolerance<typename TypeParam::Scalar>());
}

// Nonlinear CG needs many iterations in the curved valley of Rosenbrock's
// function: it is tested in double precision only
using RosenbrockTypes = DoubleVectorVariants<RosenbrockCost2, 2>;

OPTIMIST_TYPED_TEST_SUITE(RosenbrockValley, RosenbrockTypes);

TYPED_TEST(RosenbrockValley, LocalMinimum) {
  check_all_methods<TypeParam>(solver_tolerance<double>());
}
