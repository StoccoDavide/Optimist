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
      friend Solver<N, N, RootFinder<N, DerivedSolver>>;

      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != Integer(0),
        "Are you sure you want to solve a zero-dimensional system of equations?");
      static_assert(N != Integer(1),
        "C'mon, let's not kid ourselves. Use a scalar solver...");

    public:
      using Solver<N, N, DerivedSolver>::Solver;

      // I/O types
      using Vector = typename Solver<N, N, DerivedSolver>::InputType; /**< Vector type. */

      // Derivative types
      using Matrix = typename Solver<N, N, DerivedSolver>::FirstDerivativeType; /**< Jacobian matrix type. */

      // Function types
      using Function = typename Solver<N, N, DerivedSolver>::Function;         /**< Function type. */
      using Jacobian = typename Solver<N, N, DerivedSolver>::FirstDerivative;  /**< Jacobian function type. */
      using Hessian  = typename Solver<N, N, DerivedSolver>::SecondDerivative; /**< Hessian function type. */

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

    protected:
      /**
      * Evaluate the Jacobian function.
      * \param[in] jacobian Jacobian function pointer.
      * \param[in] in Input point.
      * \param[out] out Jacobian value.
      */
      void evaluate_jacobian(Jacobian jacobian, const Vector & in, Matrix & out)
      {
        this->evaluate_first_derivative(jacobian, in, out);
      }

      /**
      * Solve the root-finding problem given the function, and without derivatives.
      * \param[in] function Function pointer.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Function function, Vector const &x_ini, Vector &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(function, x_ini, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first derivative.
      * \param[in] function Function pointer.
      * \param[in] jacobian The Jacobian function pointer.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Function function, Jacobian jacobian, Vector const &x_ini, Vector &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(function, jacobian, x_ini, x_sol);
      }

    }; // class RootFinder

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOT_FINDER_HXX
