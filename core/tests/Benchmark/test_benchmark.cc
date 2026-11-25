/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Cross-solver comparison on the cost functions of the test set. The test
// always checks that every solver handles the convex problems; configure with
// -DOPTIMIST_BENCHMARK=ON to also print the full comparison table.

// Optimist library
#include "Optimist/Solvers.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;
using Vector = Eigen::VectorXd;

namespace {

  struct Record {
    std::string solver;
    std::string function;
    bool solved;
    Integer iterations;
    Integer evaluations;
  };

  // Run a solver from the first suggested guess of a function
  template <typename Solver, typename Function>
  Record run(Solver &sol, const Function &fun, const double box) {
    const Integer n{static_cast<Integer>(fun.guess(0).size())};
    sol.bounds(Vector::Constant(n, -box), Vector::Constant(n, box));
    sol.max_iterations(5000);
    sol.max_function_evaluations(100000);
    sol.max_gradient_evaluations(100000);
    sol.max_hessian_evaluations(100000);
    Vector x_out(fun.guess(0));
    bool ok{false};
    try {
      sol.optimize(fun, fun.guess(0), x_out);
      ok = sol.converged() && fun.is_solution(x_out, 1.0e-4);
    } catch (const std::exception &) {
      ok = false;
    }
    return {sol.name(), fun.name(), ok, sol.iterations(),
            sol.function_evaluations()};
  }

  // Run every solver on a function
  template <typename Function>
  std::vector<Record> run_all(const double box) {
    using namespace Optimist::Optimizer;
    Function fun;
    std::vector<Record> out;
    auto add = [&](auto &&sol) { out.push_back(run(sol, fun, box)); };

    add(ConjugateGradient<Vector>());
    {
      BFGS<Vector> s;
      add(s);
      s.method(BFGS<Vector>::Method::LBFGS);
      add(s);
    }
    {
      TrustRegion<Vector> s;
      add(s);
      s.method(TrustRegion<Vector>::Method::DOGLEG);
      add(s);
    }
    add(NelderMead<Vector>());
    add(PatternSearch<Vector>());
    {
      ParticleSwarm<Vector> s;
      s.particles(30);
      s.tolerance(1.0e-6);
      add(s);
    }
    {
      CMAES<Vector> s;
      s.restarts(6);
      s.tolerance(1.0e-6);
      add(s);
    }
    {
      Adam<Vector> s;
      s.learning_rate(0.05);
      s.tolerance(1.0e-5);
      add(s);
    }
    return out;
  }

  [[maybe_unused]] void print(const std::vector<Record> &records) {
    std::cout << std::left << std::setw(24) << "Solver" << std::setw(22)
              << "Function" << std::setw(8) << "Solved" << std::right
              << std::setw(8) << "Iter" << std::setw(10) << "#f" << '\n';
    for (const auto &r : records) {
      std::cout << std::left << std::setw(24) << r.solver << std::setw(22)
                << r.function << std::setw(8) << (r.solved ? "yes" : "no")
                << std::right << std::setw(8) << r.iterations << std::setw(10)
                << r.evaluations << '\n';
    }
  }

}  // namespace

TEST(Benchmark, EverySolverHandlesConvexProblems) {
  for (const auto &records : {run_all<EllipticParaboloid<Vector>>(10.0),
                              run_all<Trid2<Vector>>(10.0)}) {
    for (const auto &r : records) {
      EXPECT_TRUE(r.solved) << r.solver << " on " << r.function;
    }
  }
}

TEST(Benchmark, ComparisonTable) {
#ifndef OPTIMIST_BENCHMARK
  GTEST_SKIP() << "configure with -DOPTIMIST_BENCHMARK=ON to print the table";
#else
  std::vector<Record> all;
  auto append = [&all](std::vector<Record> &&r) {
    all.insert(all.end(), r.begin(), r.end());
  };
  append(run_all<EllipticParaboloid<Vector>>(10.0));
  append(run_all<Himmelblau<Vector>>(5.0));
  append(run_all<Beale<Vector>>(4.5));
  append(run_all<Trid6<Vector>>(40.0));
  append(run_all<RosenbrockCost2<Vector>>(5.0));
  append(run_all<RosenbrockCost4<Vector>>(5.0));
  append(run_all<Rastrigin2<Vector>>(5.12));
  append(run_all<Ackley2<Vector>>(5.0));
  append(run_all<Schaffer2<Vector>>(100.0));
  print(all);
#endif
}
