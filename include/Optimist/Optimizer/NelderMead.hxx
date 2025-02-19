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
     |   _   _      _     _           __  __                _
     |  | \ | | ___| | __| | ___ _ __|  \/  | ___  __ _  __| |
     |  |  \| |/ _ \ |/ _` |/ _ \ '__| |\/| |/ _ \/ _` |/ _` |
     |  | |\  |  __/ | (_| |  __/ |  | |  | |  __/ (_| | (_| |
     |  |_| \_|\___|_|\__,_|\___|_|  |_|  |_|\___|\__,_|\__,_|
     |
    \*/


    /**
    * \brief Class container for the Nelder-Mead's method.
    *
    * \includedoc docs/markdown/Optimizer/NelderMead.md
    *
    * \tparam N Dimension of the root-finding problem.
    */
    template <Integer N>
    class NelderMead : public Optimizer<N, NelderMead<N>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{false};
      static constexpr bool requires_second_derivative{false};

      using Vector = typename Optimizer<N, NelderMead<N>>::Vector;
      using Matrix = typename Optimizer<N, NelderMead<N>>::Matrix;
      using FunctionWrapper = typename Optimizer<N, NelderMead<N>>::FunctionWrapper;
      using Optimizer<N, NelderMead<N>>::solve;

      using Simplex = std::vector<Vector>; /**< Simplex type. */

    private:
      Simplex m_s;          /**< Simplex representing the vertices of the simplex. */
      Vector  m_c;          /**< The centroid of the simplex (ex_cluding the worst point). */
      Real    m_alpha{1.0}; /** Nelder-Mead reflection coefficient. */
      Real    m_gamma{2.0}; /** Nelder-Mead expansion coefficient. */
      Real    m_rho{0.5};   /** Nelder-Mead contraction coefficient. */
      Real    m_sigma{0.5}; /** Nelder-Mead shrink coefficient. */

      public:
      /**
      * Class constructor for the NelderMead solver.
      */
      NelderMead() {
        this->m_s.resize(N);
      }

      /**
      * Get the Nelder-Mead's solver name.
      * \return The Nelder-Mead's solver name.
      */
      std::string name_impl() const {return "NelderMead";}

      /**
      * Initialize the simplex around a given centroid.
      * \param x_ini The initial point.
      * \param delta The initial step size.
      */
      void initialize(const Vector & x_ini, Real delta)
      {
        this->m_s.clear();
        this->m_s.push_back(x_ini);
        Vector x;
        for (Integer i{0}; i < N; ++i) {
          x = x_ini;
          x(i) += delta;
          this->m_s.push_back(x);
        }
        this->m_c = this->centroid();
      }

      /**
      * Compute the centroid of the simplex.
      * \return The centroid of the simplex.
      */
      Vector centroid() const
      {
        return 1.0/Real(N) * std::accumulate(this->m_s.begin(), this->m_s.end(), this->m_c);
      }

      /**
      * \brief Sort the points of the simplex based on the function values and updates the centroid.
      */
      void sort(FunctionWrapper function)
      {
        std::sort(this->m_s.begin(), this->m_s.end(),
          [function, this] (const Vector & a, const Vector & b) {
            Real function_a, function_b;
            this->evaluate_function(function, a, function_a);
            this->evaluate_function(function, b, function_b);
            return function_a < function_b;
          }
        );
        this->m_c = this->centroid();
      }

      /**
      * \brief Computes the reflection point.
      * \return The reflection point.
      */
      Vector reflection() const
      {
        return this->m_c + this->m_alpha * (this->m_c - this->m_s.back());
      }

      /**
      * \brief Computes the expansion point.
      * \param x The reflection point.
      * \return The expansion point.
      */
      Vector expansion(const Vector & x) const
      {
        return this->m_c + this->m_gamma * (x - this->m_c);
      }

      /**
      * Computes the outward contraction point.
      * \param x The reflection point.
      * \return The contraction point.
      */
      Vector outward_contraction(const Vector & x) const
      {
        return this->m_c + this->m_rho * (x - this->m_c);
      }

      /**
      * Computes the inward contraction point.
      * \return The contraction point.
      */
      Vector inward_contraction() const
      {
        return this->m_c + this->m_rho * (this->m_s.back() - this->m_c);
      }

      /**
      * Shrinks the simplex towards the best point.
      */
      void shrink()
      {
        for (Integer i{1}; i < static_cast<Integer>(this->m_s.size()); ++i) {
          this->m_s.at(i) = this->m_s.at(1) + this->m_sigma * (this->m_s.at(i) - this->m_s.at(1));
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

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        Vector x_r, x_c, x_e, simplex_size;
        Real function_x_r, function_x_1, function_x_last, function_x_e, function_x_c;

        // Set initial iteration
        this->initialize(x_ini, 100.0);

        // Algorithm iterations
        for (this->m_iterations = Integer(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Sort the simplex points based on the function values
          this->sort(function);

          // Store trace
          this->store_trace(this->m_s.back());

          // Compute the reflection point
          x_r = this->reflection();
          this->evaluate_function(function, x_r, function_x_r);
          this->evaluate_function(function, this->m_s.at(1), function_x_1);
          this->evaluate_function(function, this->m_s.back(), function_x_last);

          if (function_x_1 <= function_x_r && function_x_r < function_x_last)
          {
            this->m_s.back() = x_r;
            if (this->m_verbose) {this->info(function_x_1, "Reflection");}
          }
          else if (function_x_r < function_x_1)
          {
            // Expansion
            x_e = this->expansion(x_r);
            this->evaluate_function(function, x_e, function_x_e);
            if (function_x_e < function_x_r) {
              this->m_s.back() = x_e;
            } else {
              this->m_s.back() = x_r;
            }
            if (this->m_verbose) {this->info(function_x_1, "Expansion");}
          }
          else if (function_x_r < function_x_last)
          {
            // Outward contraction
            x_c = this->outward_contraction(x_r);
            this->evaluate_function(function, x_c, function_x_c);
            if (function_x_c < function_x_r) {
              this->m_s.back() = x_c;
            } else {
              this->shrink();
            }
            if (this->m_verbose) {this->info(function_x_1, "Outward contraction");}
          }
          else
          {
            // Inward contraction
            x_c = this->inward_contraction();
            this->evaluate_function(function, x_c, function_x_c);
            if (function_x_c < function_x_last) {
              this->m_s.back() = x_c;
            } else {
              this->shrink();
            }
            if (this->m_verbose) {this->info(function_x_1, "Inward contraction");}
          }
          simplex_size = this->m_s.back() - this->m_s.front();
          if (simplex_size.norm() < this->m_tolerance) {
            this->m_converged = true;
            break;
          }
        }

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        x_sol = this->m_s.at(0);
        return this->m_converged;
      }

    }; // class NelderMead

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_NELDER_MEAD_HXX
