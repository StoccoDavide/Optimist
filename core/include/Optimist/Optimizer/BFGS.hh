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

#ifndef OPTIMIST_OPTIMIZER_BFGS_HH
#define OPTIMIST_OPTIMIZER_BFGS_HH

#include <deque>

#include "Optimist/LineSearch.hh"
#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the BFGS and limited-memory BFGS methods.
     *
     * \includedoc docs/markdown/Optimizer/BFGS.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class BFGS : public Optimizer<Vector, BFGS<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;
      using typename Optimizer<Vector, BFGS<Vector>>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Inverse Hessian approximation strategy.
       */
      using Method = enum class Method : Integer {
        BFGS  = 0, /**< Dense inverse Hessian approximation. */
        LBFGS = 1  /**< Limited-memory two-loop recursion. */
      };

     private:
      Method m_method{Method::BFGS}; /**< Inverse Hessian strategy. */
      LineSearch::Method m_line_search{
        LineSearch::Method::STRONG_WOLFE}; /**< Line search strategy. */
      LineSearch::Settings<Scalar> m_settings; /**< Line search settings. */
      Integer m_memory{8}; /**< Number of pairs stored by L-BFGS. */

     public:
      /**
       * Class constructor for the BFGS solver.
       */
      BFGS() {}

      /**
       * Get the BFGS solver name.
       */
      constexpr std::string name_impl() const {
        return this->m_method == Method::LBFGS ? "L-BFGS" : "BFGS";
      }

      /**
       * Get the inverse Hessian strategy.
       */
      Method method() const {
        return this->m_method;
      }

      /**
       * Set the inverse Hessian strategy.
       */
      void method(const Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Get the line search strategy.
       */
      LineSearch::Method line_search() const {
        return this->m_line_search;
      }

      /**
       * Set the line search strategy. Strong Wolfe is recommended, since it
       * guarantees the curvature condition needed by the BFGS update.
       */
      void line_search(const LineSearch::Method t_line_search) {
        this->m_line_search = t_line_search;
      }

      /**
       * Access the line search settings.
       */
      LineSearch::Settings<Scalar> &line_search_settings() {
        return this->m_settings;
      }

      /**
       * Get the number of correction pairs stored by L-BFGS.
       */
      Integer memory() const {
        return this->m_memory;
      }

      /**
       * Set the number of correction pairs stored by L-BFGS.
       */
      void memory(const Integer t_memory) {
        OPTIMIST_ASSERT(t_memory > 0,
                        "Optimist::Optimizer::BFGS::memory(...): invalid "
                        "input detected.");
        this->m_memory = t_memory;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] gradient Gradient of the cost function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool solve_impl(FunctionLambda &&function,
                      GradientLambda &&gradient,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::Optimizer::BFGS::solve_impl(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        const Integer n{static_cast<Integer>(x_ini.size())};
        Dense x(Utils::to_dense(x_ini)), g(n), d(n), s(n), y(n);
        Scalar f;
        DenseMatrix H(DenseMatrix::Identity(n, n));
        std::deque<std::pair<Dense, Dense>> pairs;  // L-BFGS (s, y) pairs
        Scalar gamma{1.0};                          // L-BFGS scaling

        // Evaluate the objective and the gradient at a dense point
        auto eval_f = [this, &function](const Dense &z, Scalar &out) -> bool {
          return this->evaluate_function(function,
                                         Utils::from_dense<Vector>(z),
                                         out) &&
                 std::isfinite(out);
        };
        auto eval_g = [this, &gradient](const Dense &z, Dense &out) -> bool {
          FirstDerivative tmp;
          bool ok{this->evaluate_gradient(gradient,
                                          Utils::from_dense<Vector>(z),
                                          tmp)};
          if (!ok) {
            return false;
          }
          out = Utils::to_dense(tmp);
          return out.allFinite();
        };

        OPTIMIST_ASSERT(eval_f(x, f) && eval_g(x, g),
                        CMD "evaluation failed at the initial point.");

        std::string notes{"Initialize"};
        bool fresh{true};  // Current direction is steepest descent
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          const Scalar g_norm{g.norm()};
          if (this->m_verbose) {
            this->info(g_norm, notes);
          }
          if (g_norm < this->m_tolerance) {
            this->m_converged = true;
            break;
          }

          // Search direction
          if (this->m_method == Method::BFGS) {
            d.noalias() = -H * g;
          } else {
            d = -this->two_loop(pairs, gamma, g);
          }
          if (!(g.dot(d) < -EPSILON * g_norm * d.norm())) {
            d     = -g;
            fresh = true;
            notes = "Restart";
            H.setIdentity();
            pairs.clear();
          }

          // Line search along d
          Dense x_trial(n), g_trial(n);
          Scalar f_trial{f};
          bool g_valid{false};
          auto phi = [&](const Scalar a, Scalar &v, Scalar *dv) -> bool {
            x_trial = x + a * d;
            g_valid = false;
            if (!eval_f(x_trial, v)) {
              return false;
            }
            f_trial = v;
            if (dv != nullptr) {
              if (!eval_g(x_trial, g_trial)) {
                return false;
              }
              g_valid = true;
              *dv     = g_trial.dot(d);
            }
            return true;
          };
          Scalar alpha{fresh ? std::min(Scalar(1), 1 / g_norm) : Scalar(1)};
          Scalar phi_alpha, dphi_alpha;
          bool accepted{false};
          try {
            accepted = LineSearch::search(this->m_line_search, phi, f,
                                          g.dot(d), alpha, phi_alpha,
                                          dphi_alpha, this->m_settings);
          } catch (const std::runtime_error &) {
            notes = "Budget exhausted";
            break;
          }
          if (!accepted) {
            if (fresh) {
              notes = "Line search failed";
              break;
            }
            // Retry from steepest descent with a fresh approximation
            H.setIdentity();
            pairs.clear();
            fresh = true;
            notes = "Restart";
            continue;
          }
          // Make sure the stored trial point is the accepted one
          if ((x_trial - (x + alpha * d)).norm() > 0) {
            x_trial = x + alpha * d;
            OPTIMIST_ASSERT(eval_f(x_trial, f_trial),
                            CMD "function evaluation failed.");
            g_valid = false;
          }
          if (!g_valid) {
            OPTIMIST_ASSERT(eval_g(x_trial, g_trial),
                            CMD "gradient evaluation failed.");
          }

          s = x_trial - x;
          y = g_trial - g;
          x = x_trial;
          f = f_trial;
          g = g_trial;

          // Stop on stagnation
          if (s.norm() <= EPSILON * (1 + x.norm())) {
            this->m_converged = g.norm() < this->m_tolerance;
            notes             = "Stagnation";
            break;
          }

          // Inverse Hessian update (skipped if the curvature condition fails)
          const Scalar sy{s.dot(y)};
          if (sy > SQRT_EPSILON * s.norm() * y.norm()) {
            const Scalar rho{1 / sy};
            if (this->m_method == Method::BFGS) {
              if (fresh) {
                // Shanno-Phua scaling of the initial approximation
                H.setIdentity();
                H *= sy / y.squaredNorm();
              }
              const Dense Hy{H * y};
              const Scalar yHy{y.dot(Hy)};
              H.noalias() -= rho * (s * Hy.transpose() + Hy * s.transpose());
              H.noalias() += (rho * rho * yHy + rho) * (s * s.transpose());
            } else {
              pairs.emplace_back(s, y);
              if (static_cast<Integer>(pairs.size()) > this->m_memory) {
                pairs.pop_front();
              }
              gamma = sy / y.squaredNorm();
            }
            notes = "Update";
            fresh = false;
          } else {
            notes = "Skip update";
          }
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = Utils::from_dense<Vector>(x);
        return this->m_converged;

#undef CMD
      }

     private:
      /**
       * L-BFGS two-loop recursion: compute \f$ \mathbf{H}\mathbf{g} \f$.
       */
      static Dense two_loop(const std::deque<std::pair<Dense, Dense>> &pairs,
                            const Scalar gamma,
                            const Dense &g) {
        Dense q(g);
        std::vector<Scalar> a(pairs.size());
        for (Integer i{static_cast<Integer>(pairs.size()) - 1}; i >= 0; --i) {
          const auto &[s, y] = pairs[i];
          a[i]               = s.dot(q) / y.dot(s);
          q -= a[i] * y;
        }
        Dense r(gamma * q);
        for (Integer i{0}; i < static_cast<Integer>(pairs.size()); ++i) {
          const auto &[s, y] = pairs[i];
          const Scalar b{y.dot(r) / y.dot(s)};
          r += (a[i] - b) * s;
        }
        return r;
      }

    };  // class BFGS

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_BFGS_HH
