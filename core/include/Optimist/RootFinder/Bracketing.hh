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

#ifndef OPTIMIST_ROOTFINDER_BRACKTING_HH
#define OPTIMIST_ROOTFINDER_BRACKTING_HH

#include "Optimist/RootFinder.hh"

namespace Optimist {
  namespace RootFinder {

    /**
     * \brief Class container for the Bracketing algorithms.
     *
     * The Bracketing algorithms allow to find the roots of a scalar function
     * \f$f(x)\f$ in a given interval \f$[a, b]\f$.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar, typename DerivedSolver>
    class Bracketing : public RootFinder<Scalar, DerivedSolver> {
     public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

     protected:
      Scalar m_tolerance_bracketing{
        100 * EPSILON}; /**< Tolerance for the Bracketing solver. */
      Scalar m_mu{0.5}; /**< Parameter \f$ \mu \f$. */
      Scalar m_interval_shink{0.025}; /**< Interval shrinking factor. */

      Scalar m_a{0.0}, m_fa{0.0};
      Scalar m_b{0.0}, m_fb{0.0};
      Scalar m_c{0.0}, m_fc{0.0};
      Scalar m_d{0.0}, m_fd{0.0};
      Scalar m_e{0.0}, m_fe{0.0};

     public:
      /**
       * Class constructor for the Bracketing.
       */
      Bracketing() {}

      /**
       * Get the Bracketing solver name.
       * \return The Bracketing solver name.
       */
      constexpr std::string name_impl() const {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

      /**
       * Set the tolerance for the Bracketing solver.
       * \param[in] t_tolerance The input value at which the tolerance is
       * computed.
       * \note To accurately find polynomial roots, the tolerance should be set
       * to
       * \f$ 100\epsilon(0) \f$.
       */
      void tolerance_bracketing(Scalar t_tolerance) {
        this->m_tolerance_bracketing = t_tolerance;
      }

      /**
       * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R}
       * \rightarrow \mathbb{R} \f$.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \param[in] x_ini Initialization point (not used).
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda &&function,
                      [[maybe_unused]] Scalar x_ini,
                      Scalar &x_sol) {
#define CMD "Optimist::RootFinder::Bracketing::solve_impl(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        // Function values at the bounds
        this->m_a = this->m_lower_bound;
        this->m_b = this->m_upper_bound;
        OPTIMIST_ASSERT(this->evaluate_function(function, this->m_a, this->m_fa),
                        CMD "function evaluation failed at the lower bound.");
        OPTIMIST_ASSERT(this->evaluate_function(function, this->m_b, this->m_fb),
                        CMD "function evaluation failed at the upper bound.");

        if (this->m_fa * this->m_fb > 0) {
          // Same sign at both bounds: no root is bracketed
          if (this->m_verbose) {
            this->info(std::min(std::abs(this->m_fa), std::abs(this->m_fb)),
                       "No sign change");
          }
          x_sol = std::abs(this->m_fa) < std::abs(this->m_fb) ? this->m_a
                                                              : this->m_b;
        } else {
          x_sol = this->find_root(function);
        }

        if (this->m_verbose) {
          this->bottom();
        }
        return this->m_converged;

#undef CMD
      }

      /**
       * Finds either an exact solution or an approximate solution of the
       * equation
       * \f$f(x) = 0\f$ in the interval \f$[a, b]\f$.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \return The approximate root.
       */
      template <typename FunctionLambda>
      Scalar find_root(FunctionLambda &&function) {
        return static_cast<DerivedSolver *>(this)->find_root_impl(
            function);
      }

    };  // class Bracketing

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_BRACKTING_HH
