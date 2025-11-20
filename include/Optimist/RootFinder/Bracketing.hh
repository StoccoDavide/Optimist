/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_ROOTFINDER_BRACKTING_HH
#define OPTIMIST_ROOTFINDER_BRACKTING_HH

#include "Optimist/RootFinder.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |   ____                 _        _   _
     |  | __ ) _ __ __ _  ___| | _____| |_(_)_ __   __ _
     |  |  _ \| '__/ _` |/ __| |/ / _ \ __| | '_ \ / _` |
     |  | |_) | | | (_| | (__|   <  __/ |_| | | | | (_| |
     |  |____/|_|  \__,_|\___|_|\_\___|\__|_|_| |_|\__, |
     |                                             |___/
    \*/

    /**
     * \brief Class container for the Bracketing algorithms.
     *
     * The Bracketing algorithms allow to find the roots of a scalar function \f$f(x)\f$ in a given
     * interval \f$[a, b]\f$.
     * \tparam Real Scalar number type.
     */
    template <typename Real, typename DerivedSolver>
    class Bracketing : public RootFinder<Real, 1, DerivedSolver>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Real)

    protected:
      Real m_tolerance_bracketing{100*EPSILON}; /**< Tolerance for the Bracketing solver. */
      Real m_mu{0.5};                           /**< Parameter \f$ \mu \f$. */
      Real m_interval_shink{0.025};             /**< Interval shrinking factor. */

      Real m_a{0.0}, m_fa{0.0};
      Real m_b{0.0}, m_fb{0.0};
      Real m_c{0.0}, m_fc{0.0};
      Real m_d{0.0}, m_fd{0.0};
      Real m_e{0.0}, m_fe{0.0};

    public:
      /**
       * Class constructor for the Bracketing.
       */
      Bracketing() {}

      /**
       * Get the Bracketing solver name.
       * \return The Bracketing solver name.
       */
      std::string name_impl() const
      {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

      /**
       * Set the tolerance for the Bracketing solver.
       * \param[in] t_tolerance The input value at which the tolerance is computed.
       * \note To accurately find polynomial roots, the tolerance should be set to \f$ 100\epsilon(0) \f$.
       */
      void tolerance_bracketing(Real t_tolerance) {this->m_tolerance_bracketing = t_tolerance;}

      /**
       * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R} \rightarrow \mathbb{R} \f$.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \param[in] x_ini Initialization point (not used).
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda && function, Real /*x_ini*/, Real & x_sol)
      {
        #define CMD "Optimist::Rootfinder::Bracketing::solve_impl(...): "

        // Setup internal variables
        bool success;
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        this->m_a = this->m_lower_bound;
        success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_a, this->m_fa);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the lower bound.");

        this->m_b = this->m_upper_bound;
        success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_b, this->m_fb);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the upper bound.");

        // Check if the solution exists
        if (this->m_fa*this->m_fb > 0) {return false;}
        else {x_sol = this->find_root(std::forward<FunctionLambda>(function));}

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        return this->m_converged;

        #undef CMD
      }

      /**
       * Finds either an exact solution or an approximate solution of the equation \f$f(x) = 0\f$ in
       * the interval \f$[a, b]\f$.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \return The approximate root.
       */
      template <typename FunctionLambda>
      Real find_root(FunctionLambda && function)
      {
        return static_cast<DerivedSolver *>(this)->find_root_impl(std::forward<FunctionLambda>(function));
      }

    }; // class Bracketing

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_BRACKTING_HH
