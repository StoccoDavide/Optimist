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

#ifndef OPTIMIST_OPTIMIZER_NELDER_MEAD_HXX
#define OPTIMIST_OPTIMIZER_NELDER_MEAD_HXX

namespace Optimist
{
  namespace Optimizer
  {

    /*\
     |   ____                       _ _
     |  / ___| _ __   ___ _ __   __| | | ___ _   _
     |  \___ \| '_ \ / _ \ '_ \ / _` | |/ _ \ | | |
     |   ___) | |_) |  __/ | | | (_| | |  __/ |_| |
     |  |____/| .__/ \___|_| |_|\__,_|_|\___|\__, |
     |        |_|                            |___/
    \*/

    /**
    * \brief Class container for the Spendley's method.
    *
    * \includedoc docs/markdown/Optimizer/Spendley.md
    Spendley, W., Hext, G. R., & Himsworth, F. R. (1962). Sequential Application of Simplex Designs in Optimisation and Evolutionary Operation. Technometrics, 4(4), 441–461. doi:10.1080/00401706.1962.10490033
    *
    * \tparam N Dimension of the root-finding problem.
    */
    template <Integer N>
    class Spendley : public Optimizer<N, Spendley<N>>
    {
    public:
      static constexpr bool requires_function          = true;
      static constexpr bool requires_first_derivative  = false;
      static constexpr bool requires_second_derivative = false;

      using Vector = typename Optimizer<N, Spendley<N>>::Vector;
      using Matrix = typename Optimizer<N, Spendley<N>>::Matrix;
      using FunctionWrapper = typename Optimizer<N, Spendley<N>>::FunctionWrapper;
      using JacobianWrapper = typename Optimizer<N, Spendley<N>>::JacobianWrapper;
      using Optimizer<N, Spendley<N>>::solve;

    private:
      using Design = Eigen::Matrix<Real, N, N + 1>; /**< Design matrix type. */

      Design m_design; /**< Design matrix. */

    public:
      /**
      * Class constructor for the Spendley solver.
      */
      Spendley() {}

      /**
      * Get the Spendley's solver name.
      * \return The Spendley's solver name.
      */
      std::string name_impl() const {return "Spendley";}

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

        // Print header
        if (this->m_verbose) {this->header();}

        // Build the design matrix
        this->m_design.setZero();
        Real t1{std::sqrt(Real(N + 1)) - Real(1)};
        Real t2{Real(N) * std::sqrt(Real(2))};
        Real p{delta * (Real(N) + t1) / t2};
        Real q{delta * t1 / t2};
        Vector x_old = x_ini;
        Matrix m_p(N, N + 1);
        Vector m_f(N + 1);

        for (Integer i = 0; i < N; ++i) {
          m_p.coeffRef(i, 0) = x_ini[i];
        }
        m_f.coeffRef(0) = function(m_p.col(0));

        if (this->m_verbose) {
          this->info(m_f.coeffRef(0));
        }

        for (Integer i = 0; i < N; ++i) {
          m_p.col(i + 1) = m_p.col(0).array() + p;
          m_p.coeffRef(i, i + 1) += q;
          m_f.coeffRef(i + 1) = function(m_p.col(i + 1));

          if (this->m_verbose) {
            this->info(m_f.coeffRef(i + 1));
          }
        }

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        x_sol = x_old;
        return this->m_converged;
      }

    }; // class Spendley

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_NELDER_MEAD_HXX
