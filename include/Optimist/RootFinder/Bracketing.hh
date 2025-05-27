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

#ifndef OPTIMIST_ROOT_FINDER_BRACKTING_HH
#define OPTIMIST_ROOT_FINDER_BRACKTING_HH

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
      static constexpr bool requires_function{DerivedSolver::requires_function};
      static constexpr bool requires_first_derivative{DerivedSolver::requires_first_derivative};
      static constexpr bool requires_second_derivative{DerivedSolver::requires_second_derivative};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      // Function types
      using FunctionWrapper         = typename RootFinder<Real, 1, DerivedSolver>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename RootFinder<Real, 1, DerivedSolver>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename RootFinder<Real, 1, DerivedSolver>::SecondDerivativeWrapper;

    protected:
      Real m_tolerance_bracketing{100*EPSILON}; /**< Tolerance for the Algorithm 748 solver. */
      Real m_mu{0.5}; /**< Parameter \f$ \mu \f$. */
      Real m_interval_shink{0.025}; /**< Interval shrinking factor. */

      Real m_a{0.0}, m_fa{0.0};
      Real m_b{0.0}, m_fb{0.0};
      Real m_c{0.0}, m_fc{0.0};
      Real m_d{0.0}, m_fd{0.0};
      Real m_e{0.0}, m_fe{0.0};

    public:
      /**
      * Class constructor for the Algorithm 748.
      */
      Bracketing() {}

      /**
      * Get the Algorithm 748 solver name.
      * \return The Algorithm 748 solver name.
      */
      std::string name_impl() const
      {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

      /**
      * Set the tolerance for the Algorithm 748 solver.
      * \param[in] t_tolerance The input value at which the tolerance is computed.
      * \note To accurately find polynomial roots, the tolerance should be set to \f$ 100\epsilon(0) \f$.
      */
      void tolerance_bracketing(Real t_tolerance) {this->m_tolerance_bracketing = t_tolerance;}

      /**
      * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R} \rightarrow \mathbb{R} \f$.
      * \param[in] function Function wrapper.
      * \param[in] x_ini Initialization point (not used).
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(FunctionWrapper function, Real /*x_ini*/, Real & x_sol)
      {
        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        this->m_a = this->m_lower_bound; this->evaluate_function(function, this->m_a, this->m_fa);
        this->m_b = this->m_upper_bound; this->evaluate_function(function, this->m_b, this->m_fb);

        // Check if the solution exists
        if (this->m_fa*this->m_fb > 0.0) {return false;}
        else {x_sol = this->find_root(function);}

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        return this->m_converged;
      }

      /**
      * Finds either an exact solution or an approximate solution of the equation \f$f(x) = 0\f$ in
      * the interval \f$[a, b]\f$.
      * \param[in] function Function wrapper.
      * \return The approximate root.
      */
      Real find_root(FunctionWrapper function)
      {
        return static_cast<DerivedSolver *>(this)->find_root_impl(function);
      }

    }; // class Bracketing

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOT_FINDER_BRACKTING_HH
