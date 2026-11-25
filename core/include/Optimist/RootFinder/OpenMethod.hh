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

#ifndef OPTIMIST_ROOTFINDER_OPENMETHOD_HH
#define OPTIMIST_ROOTFINDER_OPENMETHOD_HH

#include "Optimist/RootFinder.hh"

namespace Optimist {

  namespace RootFinder {

    /**
     * \brief Base class for the open (i.e., non-bracketing) one-dimensional
     * root finders of Newton type.
     *
     * \includedoc docs/markdown/ScalarRootFinder.md
     *
     * The iteration \f$ x_{k+1} = x_k + h_k \f$ is shared by all methods,
     * which only differ in the computation of the step \f$ h_k \f$. The
     * derived class provides
     * \code
     * template <typename FunctionLambda>
     * Scalar step_impl(FunctionLambda &function, Scalar x, Scalar f,
     *                  Scalar df, Scalar ddf);
     * \endcode
     * where \c ddf is zero for methods that do not require the second
     * derivative. The iteration stops when \f$ |f(x_k)| \f$ or \f$ |h_k| \f$
     * are below the tolerance and its square, respectively. In damped mode the
     * step is relaxed by the affine-invariant test of SolverBase::damp().
     *
     * \tparam Scalar Floating-point number type.
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename Scalar, typename DerivedSolver>
      requires TypeTrait<Scalar>::IsScalar
    class OpenMethod : public RootFinder<Scalar, DerivedSolver> {
     public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Solve the root-finding problem with first-derivative methods.
       * \param[in] function Function.
       * \param[in] first_derivative First derivative of the function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename FirstDerivativeLambda>
      bool solve_impl(FunctionLambda &&function,
                      FirstDerivativeLambda &&first_derivative,
                      Scalar x_ini,
                      Scalar &x_sol) {
        return this->iterate(function, first_derivative, nullptr, x_ini, x_sol);
      }

      /**
       * Solve the root-finding problem with second-derivative methods.
       * \param[in] function Function.
       * \param[in] first_derivative First derivative of the function.
       * \param[in] second_derivative Second derivative of the function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda,
                typename FirstDerivativeLambda,
                typename SecondDerivativeLambda>
      bool solve_impl(FunctionLambda &&function,
                      FirstDerivativeLambda &&first_derivative,
                      SecondDerivativeLambda &&second_derivative,
                      Scalar x_ini,
                      Scalar &x_sol) {
        return this->iterate(function,
                             first_derivative,
                             &second_derivative,
                             x_ini,
                             x_sol);
      }

     private:
      /**
       * Shared iteration. The second derivative is evaluated only when a
       * pointer to it is given.
       */
      template <typename FunctionLambda,
                typename FirstDerivativeLambda,
                typename SecondDerivativePointer>
      bool iterate(FunctionLambda &function,
                   FirstDerivativeLambda &first_derivative,
                   SecondDerivativePointer second_derivative,
                   const Scalar x_ini,
                   Scalar &x_sol) {
#define CMD "Optimist::RootFinder::" << this->name() << "::solve(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        Scalar x{x_ini}, f{0}, df{0}, ddf{0}, step{0};
        Scalar x_new{0}, f_new{0}, step_new{0};
        OPTIMIST_ASSERT(this->evaluate_function(function, x, f),
                        CMD "function evaluation failed at the initial point.");

        const Scalar tolerance_residuals{this->m_tolerance};
        const Scalar tolerance_step{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          // Derivatives
          OPTIMIST_ASSERT(this->evaluate_first_derivative(first_derivative,
                                                          x,
                                                          df),
                          CMD "first derivative evaluation failed at iteration "
                              << this->m_iterations << ".");
          if constexpr (!std::is_same_v<SecondDerivativePointer,
                                        std::nullptr_t>) {
            OPTIMIST_ASSERT(
                this->evaluate_second_derivative(*second_derivative, x, ddf),
                CMD "second derivative evaluation failed at iteration "
                    << this->m_iterations << ".");
          }

          // Step, keeping the first derivative away from zero
          std::string notes{"-"};
          if (std::abs(df) < SQRT_EPSILON) {
            df    = std::copysign(SQRT_EPSILON, df);
            notes = "Singular derivative";
          }
          step = static_cast<DerivedSolver *>(this)->step_impl(function,
                                                               x,
                                                               f,
                                                               df,
                                                               ddf);
          OPTIMIST_ASSERT(std::isfinite(step),
                          CMD "step " << this->m_iterations
                                      << " is not finite.");

          // Convergence
          if (this->m_verbose) {
            this->info(std::abs(f), notes);
          }
          if (std::abs(f) < tolerance_residuals ||
              std::abs(step) < tolerance_step) {
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
          x = x_new;
          f = f_new;
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = x;
        return this->m_converged;

#undef CMD
      }

    };  // class OpenMethod

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_OPENMETHOD_HH
