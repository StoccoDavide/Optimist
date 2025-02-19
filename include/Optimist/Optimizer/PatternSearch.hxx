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

#ifndef OPTIMIST_OPTIMIZER_HJ_PATTERN_SEARCH_HXX
#define OPTIMIST_OPTIMIZER_HJ_PATTERN_SEARCH_HXX

namespace Optimist
{
  namespace Optimizer
  {


    /*\
     |   ____       _   _                  ____                      _
     |  |  _ \ __ _| |_| |_ ___ _ __ _ __ / ___|  ___  __ _ _ __ ___| |__
     |  | |_) / _` | __| __/ _ \ '__| '_ \\___ \ / _ \/ _` | '__/ __| '_ \
     |  |  __/ (_| | |_| ||  __/ |  | | | |___) |  __/ (_| | | | (__| | | |
     |  |_|   \__,_|\__|\__\___|_|  |_| |_|____/ \___|\__,_|_|  \___|_| |_|
     |
    \*/

    /**
    * \brief Class container for the Hooke and Jeeves Pattern Search algorithm.
    *
    * \includedoc docs/markdown/Optimizer/PatternSearch.md
    *
    * \tparam N Dimension of the root-finding problem.
    */
    template <Integer N>
    class PatternSearch : public Optimizer<N, PatternSearch<N>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{false};
      static constexpr bool requires_second_derivative{false};

      using Vector = typename Optimizer<N, PatternSearch<N>>::Vector;
      using Matrix = typename Optimizer<N, PatternSearch<N>>::Matrix;
      using FunctionWrapper = typename Optimizer<N, PatternSearch<N>>::FunctionWrapper;
      using Optimizer<N, PatternSearch<N>>::solve;

      using Simplex = std::vector<Vector>; /**< Simplex type. */

    private:
      Real   m_rho{Real(0.9)}; // stencil step decreasing factor (must be 0 < rho < 1)
      Real   m_h{Real(0.1)};
      bool   m_stencil_failure{false}; // stencil failure flag - used to shrink h,
                                      // stencil_failure = true means failure

    public:
      /**
      * Class constructor for the PatternSearch solver.
      */
      PatternSearch() {}

      /**
      * Get the Nelder-Mead's solver name.
      * \return The Nelder-Mead's solver name.
      */
      std::string name_impl() const {return "PatternSearch";}

      void set_max_num_stagnation( integer nstg ) {
        UTILS_ASSERT(
          nstg > 0,
          "set_max_num_stagnation({}) argument must be >0\n", nstg
        );
        m_max_num_stagnation = nstg;
      }

      // SEARCH This method call the explore method on the first
      // iteration and then continue to call explore until a stencil
      // fails. In the case of a stencil failure, it tries once to go
      // back of half a step along the search direction by setting x_center
      // equal to the base point x_best.
      // If the stencil fails again, it exits the while loop and stencil_failure
      // is set to zero in order to signal that a reduction of h is necessary.
      void
      search() {
        ++m_iteration_count; // augment counter

        this->best_nearby();

        while ( m_stencil_failure ) {
          // reduce the scale
          m_h *= m_rho;
          this->best_nearby();
          if ( m_h <= m_tolerance ) return;
          if ( m_fun_evaluation_count >= m_max_fun_evaluation ) return;
        }

        search_dir.noalias() = x_best - x_old; // Compute search direction

        // Continue exploring until stencil failure or exceed of
        Real lambda  = 1;
        Real max_der = 0;
        while ( m_fun_evaluation_count < m_max_fun_evaluation && lambda > 0.1 ) {
          x_new.noalias() = x_best + lambda*search_dir;
          Real new_f = eval_function(x_new);
          if ( new_f < f_best-(0.25*lambda)*max_der ) {
            Real der = (f_best-new_f)/lambda;
            if ( der > max_der ) max_der = der;
            x_best.noalias() = x_new;
            f_best = new_f;
            //lambda *= 2;
          }
          lambda *= 0.5;
        }
      }

      // This method explore all points on the stencil center at
      // x_temporary = x_center and updates the current iteration x to the current
      // best point x_current_best. If the current best point x_current_best is worse than the
      // base point x_best, the current iteration x will remain constant
      // (x = x_best) and stencil failure flag stencil_failure will be set to zero.
      void best_nearby()
      {
        /*
        */
        // Initialize
        m_stencil_failure = true;

        // ----------------------------------------------------------------------------------------
        // Cycle on all stencil directions

        for ( integer j = 0; j < N; ++j ) {
          Real s_dirh = search_sign(j) * m_h;
          m_p.noalias() = x_best; m_p(j) += s_dirh;
          Real fp = eval_function( m_p );
          if ( fp >= f_best ) {
            m_p1.noalias() = x_best; m_p1(j) -= s_dirh; // try the opposite direction
            Real fp1 = eval_function( m_p1 );
            if ( fp1 < fp ) {
              m_p.noalias() = m_p1; fp = fp1;
              // change priority of search direction to the opposite verse
              search_sign(j) = -search_sign(j);
            }
          }
          // Update temporary and current best point before checking
          // the remaining directions j
          if ( fp < f_best ) {
            x_best.noalias()   = m_p;   // move temporary point
            f_best             = fp;    // new current best point
            m_stencil_failure  = false; // update stencil failure flag
          }
        }
      }

      /**
      * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0 \f$, with \f$
      * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
      * \param[in] function Function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(FunctionWrapper function, Vector const &x_ini, Vector &x_sol)
      {
        // Setup internal variables
        this->reset();
        this->m_stencil_failure = false;

        // Print header
        if (this->m_verbose) {this->header();}

        // Set initial iteration
        Vector x_best(x_ini);
        Real f_best;
        this->evaluate_function(function, x_best, f_best);
        search_sign.setOnes(); // First direction will be positive for each direction
        this->m_h = h;

        // Algorithm iterations
        integer stagnations{0};
        for (this->m_iterations = Integer(0); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {

          x_old = x_best;
          f_old = f_best;

          this->search();
          if (this->m_stencil_failure) {break;}

          // Check convergence
          if (this->m_verbose) {this->info(residuals);}
          if (this->m_h < this->m_tolerance) {
            this->m_converged = true;
            break;
          }

          // Reduce the scale
          this->m_h *= this->m_rho;

          // Check for stagnation
          if (f_old <= f_best) {
            ++stagnations;
            if (stagnations > this->m_max_stagnations) {break;}
          } else {
            stagnations = 0;
          }
        }

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        x_sol =
        return this->m_converged;
      }

    }; // class PatternSearch

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_HJ_PATTERN_SEARCH_HXX
