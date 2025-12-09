/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
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
     * \tparam T Input and output type (scalar or Eigen vector).
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename T, typename DerivedSolver>
    requires (TypeTrait<T>::IsScalar || TypeTrait<T>::IsEigen) &&
      (!TypeTrait<T>::IsFixed || TypeTrait<T>::Dimension > 0)
    class RootFinder : public SolverBase<T, T, DerivedSolver>
    {
    public:
      friend class SolverBase<T, T, DerivedSolver>;

      static constexpr bool IsRootFinder{true};
      static constexpr bool IsOptimizer{false};

      // Input and output types
      using Scalar = typename TypeTrait<T>::Scalar;
      using Input  = T;
      using Output = T;

      // Derivative types
      using typename SolverBase<T, T, DerivedSolver>::FirstDerivative;
      using typename SolverBase<T, T, DerivedSolver>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the multi-dimensional root finder.
       */
      RootFinder() {}

      /**
       * Get the solver name.
       * \return The solver name.
       */
      constexpr std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

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
      void max_jacobian_evaluations(Integer const t_jacobian_evaluations)
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
      void max_hessian_evaluations(Integer const t_hessian_evaluations)
      {
        this->max_first_derivative_evaluations(t_hessian_evaluations);
      }

    protected:
      /**
       * Evaluate the Jacobian function.
       * \tparam JacobianLambda The Jacobian lambda function type.
       * \param[in] jacobian Jacobian lambda function.
       * \param[in] x Input point.
       * \param[out] out Jacobian value.
       * \return The boolean flag for successful evaluation.
       */
      template <typename JacobianLambda>
      bool evaluate_jacobian(JacobianLambda && jacobian, Input const & x, FirstDerivative & out)
      {
        return this->evaluate_first_derivative(std::forward<JacobianLambda>(jacobian), x, out);
      }

      /**
       * Evaluate the Hessian function.
       * \tparam HessianLambda The Hessian lambda function type.
       * \param[in] hessian Hessian lambda function.
       * \param[in] x Input point.
       * \param[out] out Hessian value.
       * \return The boolean flag for successful evaluation.
       */
      template <typename HessianLambda>
      bool evaluate_hessian(HessianLambda && hessian, Input const & x, SecondDerivative & out)
      {
        return this->evaluate_second_derivative(std::forward<HessianLambda>(hessian), x, out);
      }

    public:
      /**
       * Solve the root-finding problem given the function, and without derivatives.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve(FunctionLambda && function, Input const & x_ini, Output & x_sol)
      {
        #define CMD "Optimist::RootFinder::solve(...): "

        static_assert(DerivedSolver::RequiresFunction,
          CMD "the solver requires a function.");
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          x_ini, x_sol);

        #undef CMD
      }

      /**
       * Solve the root-finding problem given the function, and its Jacobian.
       * \tparam FunctionLambda The lambda function type.
       * \tparam JacobianLambda The Jacobian lambda function type.
       * \param[in] function Function lambda.
       * \param[in] jacobian The Jacobian lambda function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename JacobianLambda>
      bool solve(FunctionLambda && function, JacobianLambda && jacobian, Input const & x_ini,
        Output & x_sol)
      {
        #define CMD "Optimist::RootFinder::solve(...): "

        static_assert(DerivedSolver::RequiresFunction,
          CMD "the solver requires a function.");
        static_assert(DerivedSolver::RequiresFirstDerivative,
          CMD "the solver requires the first derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<JacobianLambda>(jacobian),
          x_ini, x_sol);

        #undef CMD
      }

      /**
       * Solve the root-finding problem given the function, and its Jacobian and Hessian.
       * \tparam FunctionLambda The lambda function type.
       * \tparam JacobianLambda The Jacobian lambda function type.
       * \tparam HessianLambda The Hessian lambda function type.
       * \param[in] function Function lambda.
       * \param[in] jacobian The Jacobian lambda function.
       * \param[in] hessian The Hessian lambda function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename JacobianLambda, typename HessianLambda>
      bool solve(FunctionLambda && function, JacobianLambda && jacobian, HessianLambda && hessian,
        Input const & x_ini, Output & x_sol)
      {
        #define CMD "Optimist::RootFinder::solve(...): "

        static_assert(DerivedSolver::RequiresFunction,
          CMD "the solver requires the function.");
        static_assert(DerivedSolver::RequiresFirstDerivative,
          CMD "the solver requires the first derivative.");
        static_assert(DerivedSolver::RequiresSecondDerivative,
          CMD "the solver requires the second derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<JacobianLambda>(jacobian),
          std::forward<HessianLambda>(hessian),
          x_ini, x_sol);

        #undef CMD
      }

    }; // class RootFinder

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_HH
