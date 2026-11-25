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

#ifndef OPTIMIST_ROOTFINDER_LEVENBERG_MARQUARDT_HH
#define OPTIMIST_ROOTFINDER_LEVENBERG_MARQUARDT_HH

#include "Optimist/RootFinder.hh"

namespace Optimist {

  namespace RootFinder {

    /**
     * \brief Class container for the Levenberg-Marquardt method.
     *
     * \includedoc docs/markdown/RootFinder/LevenbergMarquardt.md
     *
     * \tparam Vector Eigen vector type of the unknowns.
     * \tparam OutVector Eigen vector type of the equations (by default the
     * system is square). With more equations than unknowns the system is
     * solved in the least-squares sense.
     */
    template <typename Vector, typename OutVector = Vector>
      requires TypeTrait<Vector>::IsEigen && TypeTrait<OutVector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class LevenbergMarquardt
        : public RootFinder<Vector,
                            LevenbergMarquardt<Vector, OutVector>,
                            OutVector> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;
      using typename RootFinder<Vector,
                                LevenbergMarquardt<Vector, OutVector>,
                                OutVector>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_tau{1.0e-3}; /**< Initial damping relative to diag(J'J). */
      bool m_scaling{true}; /**< Marquardt's diagonal scaling flag. */

     public:
      /**
       * Class constructor for the Levenberg-Marquardt solver.
       */
      LevenbergMarquardt() {}

      /**
       * Get the Levenberg-Marquardt solver name.
       */
      constexpr std::string name_impl() const {
        return "LevenbergMarquardt";
      }

      /**
       * Get the initial damping factor \f$ \tau \f$.
       */
      Scalar tau() const {
        return this->m_tau;
      }

      /**
       * Set the initial damping factor \f$ \tau \f$: the initial damping is
       * \f$ \lambda_0 = \tau \max_i (\mathbf{J}^\top\mathbf{J})_{ii} \f$.
       */
      void tau(const Scalar t_tau) {
        OPTIMIST_ASSERT(t_tau > 0,
                        "Optimist::RootFinder::LevenbergMarquardt::tau(...): "
                        "invalid input detected.");
        this->m_tau = t_tau;
      }

      /**
       * Enable or disable Marquardt's scaling, i.e., damping with \f$
       * \lambda\,\mathrm{diag}(\mathbf{J}^\top\mathbf{J}) \f$ instead of \f$
       * \lambda\mathbf{I} \f$.
       */
      void scaling_mode(const bool t_scaling) {
        this->m_scaling = t_scaling;
      }

      /**
       * Get Marquardt's scaling flag.
       */
      bool scaling_mode() const {
        return this->m_scaling;
      }

      /**
       * Solve the nonlinear system of equations.
       * \param[in] function Function.
       * \param[in] jacobian Jacobian of the function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename JacobianLambda>
      bool solve_impl(FunctionLambda &&function,
                      JacobianLambda &&jacobian,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::RootFinder::LevenbergMarquardt::solve_impl(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        const Integer n{static_cast<Integer>(x_ini.size())};
        Dense x(Utils::to_dense(x_ini)), F, g(n), h(n), D(n);
        DenseMatrix J, A(n, n);

        auto eval_F = [this, &function](const Dense &z, Dense &out) -> bool {
          OutVector tmp;
          bool ok{this->evaluate_function(function,
                                          Utils::from_dense<Vector>(z),
                                          tmp)};
          if (!ok) {
            return false;
          }
          out = Utils::to_dense(tmp);
          return out.allFinite();
        };
        auto eval_J = [this, &jacobian](const Dense &z,
                                        DenseMatrix &out) -> bool {
          FirstDerivative tmp;
          bool ok{this->evaluate_jacobian(jacobian,
                                          Utils::from_dense<Vector>(z),
                                          tmp)};
          if (!ok) {
            return false;
          }
          out = DenseMatrix(tmp);
          return out.allFinite();
        };

        OPTIMIST_ASSERT(eval_F(x, F) && eval_J(x, J),
                        CMD "evaluation failed at the initial point.");
        A = J.transpose() * J;
        g = J.transpose() * F;
        Scalar lambda{this->m_tau * A.diagonal().maxCoeff()};
        if (!(lambda > 0)) {
          lambda = this->m_tau;
        }
        Scalar nu{2};

        std::string notes{"Initialize"};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          const Scalar residuals{F.norm()};
          if (this->m_verbose) {
            this->info(residuals, notes);
          }
          if (residuals < this->m_tolerance) {
            this->m_converged = true;
            break;
          }
          // Overdetermined systems may have no root: stop at a stationary
          // point of |F|^2/2, i.e., where J'F vanishes relative to |F|
          if (F.size() > n &&
              g.norm() < this->m_tolerance * std::max(Scalar(1), residuals)) {
            this->m_converged = true;
            notes             = "Least-squares solution";
            break;
          }

          // Damped normal equations: (J'J + lambda D) h = -J'F
          if (this->m_scaling) {
            D = A.diagonal().cwiseMax(SQRT_EPSILON);
          } else {
            D.setOnes();
          }
          DenseMatrix M{A};
          M.diagonal() += lambda * D;
          const Eigen::LDLT<DenseMatrix> ldlt(M);
          h = -ldlt.solve(g);
          if (ldlt.info() != Eigen::Success || !h.allFinite()) {
            lambda *= nu;
            nu *= 2;
            notes = "Singular system";
            continue;
          }
          if (h.norm() <= EPSILON * (x.norm() + EPSILON)) {
            notes = "Stagnation";
            break;
          }

          // Gain ratio between actual and predicted reduction of |F|^2/2
          Dense x_new{x + h}, F_new;
          bool ok{false};
          try {
            ok = eval_F(x_new, F_new);
          } catch (const std::runtime_error &) {
            notes = "Budget exhausted";
            break;
          }
          const Scalar predicted{Scalar(0.5) * h.dot(lambda * D.cwiseProduct(h) - g)};
          const Scalar actual{ok ? Scalar(0.5) * (F.squaredNorm() - F_new.squaredNorm())
                                 : -INFTY};
          const Scalar rho{predicted > 0 ? actual / predicted : Scalar(-1)};

          if (rho > 0) {
            // Accept the step and relax the damping (Nielsen's update)
            x = x_new;
            F = F_new;
            try {
              OPTIMIST_ASSERT(eval_J(x, J), CMD "jacobian evaluation failed.");
            } catch (const std::runtime_error &) {
              notes = "Budget exhausted";
              break;
            }
            A = J.transpose() * J;
            g = J.transpose() * F;
            const Scalar t{2 * rho - 1};
            lambda *= std::max(Scalar(1) / 3, 1 - t * t * t);
            nu    = 2;
            notes = "Accept";
          } else {
            lambda *= nu;
            nu *= 2;
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

    };  // class LevenbergMarquardt

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_LEVENBERG_MARQUARDT_HH
