/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_ROOTFINDER_QUASINEWTON_HH
#define OPTIMIST_ROOTFINDER_QUASINEWTON_HH

#include "Optimist/RootFinder.hh"

namespace Optimist {
  namespace RootFinder {

    /**
     * \brief Class container for the QuasiNewton's method.
     *
     * \includedoc docs/markdown/RootFinder/QuasiNewton.md
     *
     * \tparam Vector Eigen vector type.
     * \tparam DerivedSolver Derived solver type.
     */
    template <typename Vector, typename DerivedSolver>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class QuasiNewton : public RootFinder<Vector, DerivedSolver> {
     public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename TypeTrait<Vector>::Scalar;
      using typename RootFinder<Vector, DerivedSolver>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the QuasiNewton solver.
       */
      QuasiNewton() {}

      /**
       * Get the QuasiNewton solver name.
       * \return The QuasiNewton solver name.
       */
      constexpr std::string name_impl() const {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

      /**
       * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0
       * \f$, with \f$
       * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
       * \tparam FunctionLambda Function lambda type.
       * \tparam JacobianLambda Jacobian lambda type.
       * \param[in] function Function lambda.
       * \param[in] jacobian Jacobian lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename JacobianLambda>
      bool solve_impl(FunctionLambda &&function,
                      JacobianLambda &&jacobian,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::RootFinder::QuasiNewton::solve(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        // H approximates the inverse of the Jacobian: it starts from the
        // inverse of the exact Jacobian (or the identity if it is singular)
        Vector x(x_ini), f, step, x_new, f_new, step_new, dx, df, dx_old,
            df_old;
        FirstDerivative H, H_new;
        OPTIMIST_ASSERT(this->evaluate_function(function, x, f),
                        CMD "function evaluation failed at the initial point.");
        OPTIMIST_ASSERT(this->evaluate_jacobian(jacobian, x, H),
                        CMD "jacobian evaluation failed at the initial point.");
        this->invert_jacobian(H);

        const Scalar tolerance_residuals{this->m_tolerance};
        const Scalar tolerance_step{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          // Quasi-Newton step
          step = -H * f;

          // Convergence
          if (this->m_verbose) {
            this->info(f.norm());
          }
          if (f.norm() < tolerance_residuals || step.norm() < tolerance_step) {
            this->m_converged = true;
            break;
          }

          // Update
          if (this->m_damped) {
            this->damp(function, x, f, step, x_new, f_new, step_new);
          } else {
            x_new = x + step;
            OPTIMIST_ASSERT(this->evaluate_function(function, x_new, f_new),
                            CMD "function evaluation failed at iteration "
                                << this->m_iterations << ".");
          }

          // Secant update of the inverse Jacobian
          dx = x_new - x;
          df = f_new - f;
          this->update(dx_old, df_old, H, dx, df, f_new, H_new);
          x      = x_new;
          f      = f_new;
          dx_old = dx;
          df_old = df;
          H      = H_new;
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = x;
        return this->m_converged;

#undef CMD
      }

      /**
       * Invert the Jacobian matrix in place. If the matrix is singular, it is
       * replaced by the identity matrix.
       * \param[in,out] jacobian Jacobian matrix to be inverted.
       */
      static void invert_jacobian(FirstDerivative &jacobian) {
        using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        const Dense J(jacobian);
        const Eigen::FullPivLU<Dense> lu(J);
        Dense J_inv(J.rows(), J.cols());
        if (lu.isInvertible()) {
          J_inv = lu.inverse();
        } else {
          J_inv.setIdentity();
        }
        if constexpr (VectorTrait::IsSparse) {
          jacobian = J_inv.sparseView();
        } else {
          jacobian = J_inv;
        }
      }

      /**
       * Jacobian approximation update rule for the QuasiNewton's method.
       * \param[in] delta_x_old Old difference between points.
       * \param[in] delta_function_old Old difference between function values.
       * \param[in] jacobian_old Old jacobian approximation.
       * \param[in] delta_x_new New difference between points.
       * \param[in] delta_function_new New difference between function values.
       * \param[in] function_new New function value.
       * \param[out] jacobian_new New jacobian approximation.
       */
      void update(const Vector &delta_x_old,
                  const Vector &delta_function_old,
                  const FirstDerivative &jacobian_old,
                  const Vector &delta_x_new,
                  const Vector &delta_function_new,
                  const Vector &function_new,
                  FirstDerivative &jacobian_new) {
        static_cast<DerivedSolver *>(this)->update_impl(
            delta_x_old,
            delta_function_old,
            jacobian_old,  // Old step data
            delta_x_new,
            delta_function_new,
            function_new,
            jacobian_new  // New step data
        );
      }

    };  // class QuasiNewton

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_QUASINEWTON_HH
