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

#ifndef OPTIMIST_OPTIMIZER_TRUST_REGION_HH
#define OPTIMIST_OPTIMIZER_TRUST_REGION_HH

#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the trust-region Newton method.
     *
     * \includedoc docs/markdown/Optimizer/TrustRegion.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class TrustRegion : public Optimizer<Vector, TrustRegion<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{true};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;
      using typename Optimizer<Vector, TrustRegion<Vector>>::FirstDerivative;
      using typename Optimizer<Vector, TrustRegion<Vector>>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Trust-region subproblem solver.
       */
      using Method = enum class Method : Integer {
        STEIHAUG = 0, /**< Truncated conjugate gradient (Steihaug-Toint). */
        DOGLEG   = 1  /**< Powell's dogleg (falls back to Steihaug-Toint if
                         the Hessian is not positive definite). */
      };

     private:
      Method m_method{Method::STEIHAUG}; /**< Subproblem solver. */
      Scalar m_radius{1.0};              /**< Initial trust-region radius. */
      Scalar m_max_radius{1.0e3};        /**< Maximum trust-region radius. */
      Scalar m_eta{1.0e-4}; /**< Minimum ratio to accept a step. */

     public:
      /**
       * Class constructor for the trust-region solver.
       */
      TrustRegion() {}

      /**
       * Get the trust-region solver name.
       */
      constexpr std::string name_impl() const {
        return this->m_method == Method::DOGLEG ? "TrustRegion(Dogleg)"
                                                : "TrustRegion(Steihaug)";
      }

      /**
       * Get the subproblem solver.
       */
      Method method() const {
        return this->m_method;
      }

      /**
       * Set the subproblem solver.
       */
      void method(const Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Get the initial trust-region radius.
       */
      Scalar radius() const {
        return this->m_radius;
      }

      /**
       * Set the initial trust-region radius.
       */
      void radius(const Scalar t_radius) {
        OPTIMIST_ASSERT(t_radius > 0 && t_radius <= this->m_max_radius,
                        "Optimist::Optimizer::TrustRegion::radius(...): "
                        "invalid input detected.");
        this->m_radius = t_radius;
      }

      /**
       * Get the maximum trust-region radius.
       */
      Scalar max_radius() const {
        return this->m_max_radius;
      }

      /**
       * Set the maximum trust-region radius.
       */
      void max_radius(const Scalar t_max_radius) {
        OPTIMIST_ASSERT(t_max_radius >= this->m_radius,
                        "Optimist::Optimizer::TrustRegion::max_radius(...): "
                        "invalid input detected.");
        this->m_max_radius = t_max_radius;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] gradient Gradient of the cost function.
       * \param[in] hessian Hessian of the cost function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda,
                typename GradientLambda,
                typename HessianLambda>
      bool solve_impl(FunctionLambda &&function,
                      GradientLambda &&gradient,
                      HessianLambda &&hessian,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::Optimizer::TrustRegion::solve_impl(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        const Integer n{static_cast<Integer>(x_ini.size())};
        Dense x(Utils::to_dense(x_ini)), g(n), p(n);
        DenseMatrix B(n, n);
        Scalar f, delta{this->m_radius};

        auto eval_f = [this, &function](const Dense &z, Scalar &out) -> bool {
          return this->evaluate_function(function,
                                         Utils::from_dense<Vector>(z),
                                         out) &&
                 std::isfinite(out);
        };
        auto eval_derivatives = [this, &gradient, &hessian](
                                    const Dense &z,
                                    Dense &g_out,
                                    DenseMatrix &B_out) -> bool {
          const Vector zv{Utils::from_dense<Vector>(z)};
          FirstDerivative gv;
          SecondDerivative Bv;
          bool ok{this->evaluate_gradient(gradient, zv, gv) &&
                  this->evaluate_hessian(hessian, zv, Bv)};
          if (!ok) {
            return false;
          }
          g_out = Utils::to_dense(gv);
          B_out = DenseMatrix(Bv);
          return g_out.allFinite() && B_out.allFinite();
        };

        OPTIMIST_ASSERT(eval_f(x, f) && eval_derivatives(x, g, B),
                        CMD "evaluation failed at the initial point.");

        std::string notes{"Initialize"};
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
          if (delta < EPSILON * (1 + x.norm())) {
            notes = "Radius collapsed";
            break;
          }

          // Approximately solve min m(p) = g'p + p'Bp/2, |p| <= delta
          if (this->m_method == Method::DOGLEG) {
            this->dogleg(g, B, delta, p);
          } else {
            this->steihaug(g, B, delta, p);
          }
          const Scalar predicted{-(g.dot(p) + Scalar(0.5) * p.dot(B * p))};

          // Actual vs predicted reduction
          Dense x_trial{x + p};
          Scalar f_trial;
          Scalar rho{-1};
          bool ok{false};
          try {
            ok = eval_f(x_trial, f_trial);
          } catch (const std::runtime_error &) {
            notes = "Budget exhausted";
            break;
          }
          if (ok && predicted > 0) {
            rho = (f - f_trial) / predicted;
          }

          // Radius update
          const Scalar p_norm{p.norm()};
          if (rho < 0.25) {
            delta = 0.25 * p_norm;
            notes = "Shrink";
          } else if (rho > 0.75 && p_norm >= 0.99 * delta) {
            delta = std::min(2 * delta, this->m_max_radius);
            notes = "Expand";
          } else {
            notes = "Keep";
          }

          // Step acceptance
          if (rho > this->m_eta) {
            x = x_trial;
            f = f_trial;
            try {
              OPTIMIST_ASSERT(eval_derivatives(x, g, B),
                              CMD "derivative evaluation failed.");
            } catch (const std::runtime_error &) {
              notes = "Budget exhausted";
              break;
            }
          } else {
            notes = "Reject";
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
       * Positive root \f$ \tau \f$ of \f$ \|\mathbf{z} + \tau\mathbf{d}\| =
       * \Delta \f$.
       */
      static Scalar to_boundary(const Dense &z,
                                const Dense &d,
                                const Scalar delta) {
        const Scalar a{d.squaredNorm()}, b{2 * z.dot(d)},
            c{z.squaredNorm() - delta * delta};
        const Scalar disc{std::max(Scalar(0), b * b - 4 * a * c)};
        return (-b + std::sqrt(disc)) / (2 * a);
      }

      /**
       * Steihaug-Toint truncated conjugate gradient.
       */
      static void steihaug(const Dense &g,
                           const DenseMatrix &B,
                           const Scalar delta,
                           Dense &p) {
        const Integer n{static_cast<Integer>(g.size())};
        const Scalar g_norm{g.norm()};
        const Scalar eps{std::min(Scalar(0.5), std::sqrt(g_norm)) * g_norm};
        Dense z(Dense::Zero(n)), r(g), d(-g);
        for (Integer j{0}; j < 2 * n; ++j) {
          const Dense Bd{B * d};
          const Scalar dBd{d.dot(Bd)};
          if (dBd <= 0) {
            // Negative curvature: go to the boundary
            p = z + to_boundary(z, d, delta) * d;
            return;
          }
          const Scalar alpha{r.squaredNorm() / dBd};
          const Dense z_new{z + alpha * d};
          if (z_new.norm() >= delta) {
            p = z + to_boundary(z, d, delta) * d;
            return;
          }
          const Dense r_new{r + alpha * Bd};
          if (r_new.norm() < eps) {
            p = z_new;
            return;
          }
          const Scalar beta{r_new.squaredNorm() / r.squaredNorm()};
          d = -r_new + beta * d;
          r = r_new;
          z = z_new;
        }
        p = z;
      }

      /**
       * Powell's dogleg step.
       */
      static void dogleg(const Dense &g,
                         const DenseMatrix &B,
                         const Scalar delta,
                         Dense &p) {
        const Eigen::LLT<DenseMatrix> llt(B);
        if (llt.info() != Eigen::Success) {
          steihaug(g, B, delta, p);
          return;
        }
        const Dense p_newton{-llt.solve(g)};
        if (p_newton.norm() <= delta) {
          p = p_newton;
          return;
        }
        const Dense p_cauchy{-(g.squaredNorm() / g.dot(B * g)) * g};
        if (p_cauchy.norm() >= delta) {
          p = (delta / p_cauchy.norm()) * p_cauchy;
          return;
        }
        p = p_cauchy + to_boundary(p_cauchy, p_newton - p_cauchy, delta) *
                           (p_newton - p_cauchy);
      }

    };  // class TrustRegion

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_TRUST_REGION_HH
