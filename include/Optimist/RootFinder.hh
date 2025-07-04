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

#ifndef OPTIMIST_ROOTFINDER_HH
#define OPTIMIST_ROOTFINDER_HH

#include "Optimist/SolverBase.hh"

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
    * \tparam Real Scalar number type.
    * \tparam N The dimension of the root-finding problem.
    * \tparam DerivedSolver Derived solver class.
    * \tparam ForceEigen Force the use of Eigen types for input and output.
    */
    template <typename Real, Integer N, typename DerivedSolver, bool ForceEigen = false>
    class RootFinder : public SolverBase<Real, N, N, DerivedSolver, ForceEigen>
    {
    public:
      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != static_cast<Integer>(0),
        "Are you sure you want to solve a zero-dimensional system of equations?");

      friend class SolverBase<Real, N, N, RootFinder<Real, N, DerivedSolver, ForceEigen>>;

      static constexpr bool is_rootfinder{true};
      static constexpr bool is_optimizer{false};

      static constexpr bool requires_function{DerivedSolver::requires_function};
      static constexpr bool requires_first_derivative{DerivedSolver::requires_first_derivative};
      static constexpr bool requires_second_derivative{DerivedSolver::requires_second_derivative};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      using SolverBase<Real, N, N, DerivedSolver, ForceEigen>::solve;

      // I/O types
      using Vector = typename SolverBase<Real, N, N, DerivedSolver, ForceEigen>::InputType; /**< Vector type. */

      // Derivative types
      using Matrix = typename SolverBase<Real, N, N, DerivedSolver, ForceEigen>::FirstDerivativeType;  /**< Jacobian matrix type. */
      using Tensor = typename SolverBase<Real, N, N, DerivedSolver, ForceEigen>::SecondDerivativeType; /**< Hessian tensor type. */

      // Function types
      using FunctionWrapper = typename SolverBase<Real, N, N, DerivedSolver, ForceEigen>::FunctionWrapper;
      using JacobianWrapper = typename SolverBase<Real, N, N, DerivedSolver, ForceEigen>::FirstDerivativeWrapper;
      using HessianWrapper  = typename SolverBase<Real, N, N, DerivedSolver, ForceEigen>::SecondDerivativeWrapper;

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
        this->max_first_derivative_evaluations(t_jacobian_evaluations);
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
        this->max_first_derivative_evaluations(t_hessian_evaluations);
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
        this->evaluate_second_derivative(hessian, x, out);
      }

    public:
      /**
      * Solve the root-finding problem given the function, and without derivatives.
      * \param[in] function Function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, Vector const & x_ini, Vector & x_sol)
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
      bool solve(FunctionWrapper function, JacobianWrapper jacobian, Vector const & x_ini, Vector & x_sol)
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
        const & x_ini, Vector & x_sol)
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
    * \brief Class container for the scalar scalar root-finder.
    *
    * \includedoc docs/markdown/RootFinder.md
    *
    * \tparam Real Scalar number type.
    * \tparam DerivedSolver Derived solver class.
    */
    template <typename Real, typename DerivedSolver>
    class RootFinder<Real, 1, DerivedSolver> : public SolverBase<Real, 1, 1, DerivedSolver>
    {
    public:
      friend class SolverBase<Real, 1, 1, RootFinder<Real, 1, DerivedSolver>>;

      static constexpr bool is_rootfinder{true};
      static constexpr bool is_optimizer{false};

      static constexpr bool requires_function{DerivedSolver::requires_function};
      static constexpr bool requires_first_derivative{DerivedSolver::requires_first_derivative};
      static constexpr bool requires_second_derivative{DerivedSolver::requires_second_derivative};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      using FunctionWrapper         = typename SolverBase<Real, 1, 1, DerivedSolver>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename SolverBase<Real, 1, 1, DerivedSolver>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename SolverBase<Real, 1, 1, DerivedSolver>::SecondDerivativeWrapper;

      /**
      * Class constructor for the scalar root-finder.
      */
      RootFinder<Real, 1, DerivedSolver>() {}

      /**
      * Get the solver name.
      * \return The solver name.
      */
      std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

      /**
      * Solve the root-finding problem given the function, and without derivatives.
      * \param[in] function Function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, Real x_ini, Real & x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(function, x_ini, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first derivative.
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, FirstDerivativeWrapper first_derivative, Real x_ini,
        Real & x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(function, first_derivative, x_ini, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first and second derivatives.
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] second_derivate Second derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, FirstDerivativeWrapper first_derivative, SecondDerivativeWrapper
        second_derivate, Real x_ini, Real & x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(function, first_derivative,
          second_derivate, x_ini, x_sol);
      }

    }; // class RootFinder

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_HH
