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

#ifndef OPTIMIST_OPTIMIZER_HXX
#define OPTIMIST_OPTIMIZER_HXX

namespace Optimist
{

  /**
  * \brief Namespace for multi-dimensional optimization algorithms.
  */
  namespace Optimizer
  {

    /*\
     |    ___        _   _           _
     |   / _ \ _ __ | |_(_)_ __ ___ (_)_______ _ __
     |  | | | | '_ \| __| | '_ ` _ \| |_  / _ \ '__|
     |  | |_| | |_) | |_| | | | | | | |/ /  __/ |
     |   \___/| .__/ \__|_|_| |_| |_|_/___\___|_|
     |        |_|
    \*/

    /**
    * \brief Class container for the multi-dimensional optimizer.
    *
    * \includedoc docs/markdown/Optimizer.md
    *
    * \tparam N Dimension of the optimization problem.
    * \tparam DerivedSolver Derived solver class.
    */
    template <Integer N, typename DerivedSolver>
    class Optimizer : public Solver<N, 1, DerivedSolver>
    {
    public:
      friend Solver<N, 1, Optimizer<N, DerivedSolver>>;

      static constexpr bool is_rootfinder{false};
      static constexpr bool is_optimizer{true};

      static constexpr bool requires_function{DerivedSolver::requires_function};
      static constexpr bool requires_first_derivative{DerivedSolver::requires_first_derivative};
      static constexpr bool requires_second_derivative{DerivedSolver::requires_second_derivative};

      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != Integer(0),
        "Have you ever heard of a zero-dimensional optimization problem?");
      static_assert(N != Integer(1),
        "Good try, but you should use a scalar solver for a one-dimensional optimization problem.");

      // I/O types
      using Vector = typename Solver<N, 1, DerivedSolver>::InputType; /**< Vector type. */

      // Derivative types
      using RowVector = typename Solver<N, 1, DerivedSolver>::FirstDerivativeType;  /**< Gradient (row) vector type. */
      using Matrix    = typename Solver<N, 1, DerivedSolver>::SecondDerivativeType; /**< Hessian matrix type. */

      // Function types
      using FunctionWrapper = typename Solver<N, 1, DerivedSolver>::FunctionWrapper;
      using GradientWrapper = typename Solver<N, 1, DerivedSolver>::FirstDerivativeWrapper;
      using HessianWrapper  = typename Solver<N, 1, DerivedSolver>::SecondDerivativeWrapper;

      /**
      * Class constructor for the multi-dimensional optimizer.
      */
      Optimizer() {}

      /**
      * Get the solver name.
      * \return The solver name.
      */
      std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

      /**
      * Get the number of gradient evaluations.
      * \return The number of gradient evaluations.
      */
      Integer gradient_evaluations() const {return this->first_derivative_evaluations();}

      /**
      * Get the number of maximum allowed gradient evaluations.
      * \return The number of maximum allowed gradient evaluations.
      */
      Integer max_gradient_evaluations() const {return this->max_first_derivative_evaluations();}

      /**
      * Set the number of maximum allowed gradient evaluations.
      * \param[in] t_gradient_evaluations The number of maximum allowed gradient evaluations.
      */
      void max_gradient_evaluations(Integer t_gradient_evaluations)
      {
        this->max_first_derivative_evaluations(t_gradient_evaluations);
      }

      /**
      * Get the number of hessian evaluations.
      * \return The number of hessian evaluations.
      */
      Integer hessian_evaluations() const {return this->second_derivative_evaluations();}

      /**
      * Get the number of maximum allowed hessian evaluations.
      * \return The number of maximum allowed hessian evaluations.
      */
      Integer max_hessian_evaluations() const {return this->max_second_derivative_evaluations();}

      /**
      * Set the number of maximum allowed hessian evaluations.
      * \param[in] t_hessian_evaluations The number of maximum allowed hessian evaluations.
      */
      void max_hessian_evaluations(Integer t_hessian_evaluations)
      {
        this->max_second_derivative_evaluations(t_hessian_evaluations);
      }

    protected:
      /**
      * Evaluate the gradient function.
      * \param[in] gradient Gradient function wrapper.
      * \param[in] x Input point.
      * \param[out] out Gradient value.
      */
      void evaluate_gradient(GradientWrapper gradient, const Vector & x, Matrix & out)
      {
        this->evaluate_first_derivative_impl(gradient, x, out);
      }

      /**
      * Evaluate the hessian function.
      * \param[in] hessian Hessian function wrapper.
      * \param[in] x Input point.
      * \param[out] out Hessian value.
      */
      void evaluate_hessian(HessianWrapper hessian, const Vector & x, Matrix & out)
      {
        this->evaluate_second_derivative_impl(hessian, x, out);
      }

      /**
      * Solve the root-finding problem given the function, and without derivatives.
      * \param[in] function Function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, Vector const &x_ini, Vector &x_sol)
      {
        #define CMD "Optimist::Optimizer::solve(...): "

        static_assert(DerivedSolver::requires_function,
          CMD "the solver requires the function.");
        return static_cast<DerivedSolver *>(this)->solve_impl(function, x_ini, x_sol);

        #undef CMD
      }

      /**
      * Solve the root-finding problem given the function, and its gradient.
      * \param[in] function Function wrapper.
      * \param[in] gradient Gradient function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, GradientWrapper gradient, Vector const &x_ini, Vector &x_sol)
      {
        #define CMD "Optimist::Optimizer::solve(...): "

        static_assert(DerivedSolver::requires_function,
          CMD "the solver requires the function.");
        static_assert(DerivedSolver::requires_first_derivative,
          CMD "the solver requires the first derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(function, gradient, x_ini, x_sol);

        #undef CMD
      }

      /**
      * Solve the root-finding problem given the function, and its gradient and Hessian.
      * \param[in] function Function wrapper.
      * \param[in] gradient Gradient function wrapper.
      * \param[in] hessian Hessian function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, GradientWrapper gradient, HessianWrapper hessian, Vector
        const &x_ini, Vector &x_sol)
      {
        #define CMD "Optimist::Optimizer::solve(...): "

        static_assert(DerivedSolver::requires_function,
          CMD "the solver requires the function.");
        static_assert(DerivedSolver::requires_first_derivative,
          CMD "the solver requires the first derivative.");
        static_assert(DerivedSolver::requires_second_derivative,
          CMD "the solver requires the second derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(function, gradient, hessian, x_ini, x_sol);

        #undef CMD
      }

    }; // class Optimizer

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_HXX
