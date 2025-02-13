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

#ifndef OPTIMIST_ROOT_FINDER_HXX
#define OPTIMIST_ROOT_FINDER_HXX

namespace Optimist
{

  /**
  * \brief Namespace for multi-dimensional root-finding algorithms.
  */
  namespace RootFinder
  {

    /*\
     |   ____             _   _____ _           _
     |  |  _ \ ___   ___ | |_|  ___(_)_ __   __| | ___ _ __
     |  | |_) / _ \ / _ \| __| |_  | | '_ \ / _` |/ _ \ '__|
     |  |  _ < (_) | (_) | |_|  _| | | | | | (_| |  __/ |
     |  |_| \_\___/ \___/ \__|_|   |_|_| |_|\__,_|\___|_|
     |
    \*/

    /**
    * \brief Class container for the multi-dimensional root finder.
    *
    * \includedoc docs/markdown/RootFinder.md
    *
    * \tparam N The dimension of the root-finding problem.
    * \tparam DerivedSolver Derived solver class.
    */
    template <Integer N, typename DerivedSolver>
    class RootFinder : public Solver<N, N, DerivedSolver>
    {
    public:
      friend Solver<N, N, RootFinder<N, DerivedSolver>>;

      static constexpr bool requires_function          = DerivedSolver::requires_function;
      static constexpr bool requires_first_derivative  = DerivedSolver::requires_first_derivative;
      static constexpr bool requires_second_derivative = DerivedSolver::requires_second_derivative;

      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != Integer(0),
        "Are you sure you want to solve a zero-dimensional system of equations?");
      static_assert(N != Integer(1),
        "C'mon, let's not kid ourselves. Use a scalar solver...");

      using Solver<N, N, DerivedSolver>::Solver;

      // I/O types
      using Vector = typename Solver<N, N, DerivedSolver>::InputType; /**< Vector type. */

      // Derivative types
      using Matrix = typename Solver<N, N, DerivedSolver>::FirstDerivativeType;  /**< Jacobian matrix type. */
      using Tensor = typename Solver<N, N, DerivedSolver>::SecondDerivativeType; /**< Hessian tensor type. */

      // Function types
      using FunctionWrapper = typename Solver<N, N, DerivedSolver>::FunctionWrapper;
      using JacobianWrapper = typename Solver<N, N, DerivedSolver>::FirstDerivativeWrapper;
      using HessianWrapper  = typename Solver<N, N, DerivedSolver>::SecondDerivativeWrapper;

      /**
      * Class constructor for the multi-dimensional root finder.
      */
      RootFinder() {}

      /**
      * Get the solver name.
      * \return The solver name.
      */
      std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

      /**
      * Get the number of Jacobian evaluations.
      * \return The number of Jacobian evaluations.
      */
      Integer jacobian_evaluations() const {return this->first_derivative_evaluations();}

      /**
      * Get the number of maximum allowed Jacobian evaluations.
      * \return The number of maximum allowed Jacobian evaluations.
      */
      Integer max_jacobian_evaluations() const {return this->max_first_derivative_evaluations();}

      /**
      * Set the number of maximum allowed Jacobian evaluations.
      * \param[in] t_jacobian_evaluations The number of maximum allowed Jacobian evaluations.
      */
      void max_jacobian_evaluations(Integer t_jacobian_evaluations)
      {
        this->first_derivative_evaluations(t_jacobian_evaluations);
      }

      /**
      * Get the number of Hessian evaluations.
      * \return The number of Hessian evaluations.
      */
      Integer hessian_evaluations() const {return this->first_derivative_evaluations();}

      /**
      * Get the number of maximum allowed Hessian evaluations.
      * \return The number of maximum allowed Hessian evaluations.
      */
      Integer max_hessian_evaluations() const {return this->max_first_derivative_evaluations();}

      /**
      * Set the number of maximum allowed Hessian evaluations.
      * \param[in] t_hessian_evaluations The number of maximum allowed Hessian evaluations.
      */
      void max_hessian_evaluations(Integer t_hessian_evaluations)
      {
        this->first_derivative_evaluations(t_hessian_evaluations);
      }

    protected:
      /**
      * Evaluate the Jacobian function.
      * \param[in] jacobian Jacobian function wrapper.
      * \param[in] x Input point.
      * \param[out] out Jacobian value.
      */
      void evaluate_jacobian(JacobianWrapper jacobian, const Vector & x, Matrix & out)
      {
        this->evaluate_first_derivative(jacobian, x, out);
      }

      /**
      * Evaluate the Hessian function.
      * \param[in] hessian Hessian function wrapper.
      * \param[in] x Input point.
      * \param[out] out Hessian value.
      */
      void evaluate_hessian(HessianWrapper hessian, const Vector & x, Matrix & out)
      {
        this->evaluate_first_derivative(hessian, x, out);
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
        #define CMD "Optimist::RootFinder::solve(...): "

        static_assert(DerivedSolver::requires_function,
          CMD "the solver requires a function.");
        return static_cast<DerivedSolver *>(this)->solve_impl(function, x_ini, x_sol);

        #undef CMD
      }

      /**
      * Solve the root-finding problem given the function, and its Jacobian.
      * \param[in] function Function wrapper.
      * \param[in] jacobian The Jacobian function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, JacobianWrapper jacobian, Vector const &x_ini, Vector &x_sol)
      {
        #define CMD "Optimist::RootFinder::solve(...): "

        static_assert(DerivedSolver::requires_function,
          CMD "the solver requires a function.");
        static_assert(DerivedSolver::requires_first_derivative,
          CMD "the solver requires the first derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(function, jacobian, x_ini, x_sol);

        #undef CMD
      }

      /**
      * Solve the root-finding problem given the function, and its Jacobian and Hessian.
      * \param[in] function Function wrapper.
      * \param[in] jacobian The Jacobian function wrapper.
      * \param[in] hessian The Hessian function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, JacobianWrapper jacobian, HessianWrapper hessian, Vector
        const &x_ini, Vector &x_sol)
      {
        #define CMD "Optimist::RootFinder::solve(...): "

        static_assert(DerivedSolver::requires_function,
          CMD "the solver requires the function.");
        static_assert(DerivedSolver::requires_first_derivative,
          CMD "the solver requires the first derivative.");
        static_assert(DerivedSolver::requires_second_derivative,
          CMD "the solver requires the second derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(function, jacobian, hessian, x_ini, x_sol);

        #undef CMD
      }

    }; // class RootFinder

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOT_FINDER_HXX
