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

#ifndef OPTIMIST_OPTIMIZER_HH
#define OPTIMIST_OPTIMIZER_HH

#include "Optimist/SolverBase.hh"

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
     * \tparam ForceEigen Force the use of Eigen types for input and output.
     */
    template <typename Real, Integer N, typename DerivedSolver, bool ForceEigen = false>
    class Optimizer : public SolverBase<Real, N, 1, DerivedSolver, ForceEigen>
    {
    public:
      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != 0, "Have you ever heard of a zero-dimensional optimization problem?");

      friend class SolverBase<Real, N, 1, Optimizer<Real, N, DerivedSolver, ForceEigen>>;

      static constexpr bool IsRootFinder{false};
      static constexpr bool IsOptimizer{true};

      OPTIMIST_BASIC_CONSTANTS(Real)

      // Input and output types
      using Vector = typename SolverBase<Real, N, 1, DerivedSolver, ForceEigen>::InputType;

      // Derivative types
      using RowVector = typename SolverBase<Real, N, 1, DerivedSolver, ForceEigen>::FirstDerivativeType;
      using Matrix    = typename SolverBase<Real, N, 1, DerivedSolver, ForceEigen>::SecondDerivativeType;

      /**
       * Class constructor for the multi-dimensional optimizer.
       */
      Optimizer() {}

      /**
       * Get the solver name.
       * \return The solver name.
       */
      constexpr std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

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
       * \tparam GradientLambda The gradient lambda function type.
       * \param[in] gradient Gradient lambda function.
       * \param[in] x Input point.
       * \param[out] out Gradient value.
       * \return The boolean flag for successful evaluation.
       */
      template <typename GradientLambda>
      bool evaluate_gradient(GradientLambda && gradient, Vector const & x, Matrix & out)
      {
        return this->evaluate_first_derivative(std::forward<GradientLambda>(gradient), x, out);
      }

      /**
       * Evaluate the hessian function.
       * \tparam HessianLambda The hessian lambda function type.
       * \param[in] hessian Hessian lambda function.
       * \param[in] x Input point.
       * \param[out] out Hessian value.
       * \return The boolean flag for successful evaluation.
       */
      template <typename HessianLambda>
      bool evaluate_hessian(HessianLambda && hessian, Vector const & x, Matrix & out)
      {
        return this->evaluate_second_derivative(std::forward<HessianLambda>(hessian), x, out);
      }

      /**
       * Solve the root-finding problem given the function, and without derivatives.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve(FunctionLambda && function, Vector const & x_ini, Vector & x_sol)
      {
        #define CMD "Optimist::Optimizer::solve(...): "

        static_assert(DerivedSolver::RequiresFunction,
          CMD "the solver requires the function.");
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          x_ini, x_sol);

        #undef CMD
      }

      /**
       * Solve the root-finding problem given the function, and its gradient.
       * \tparam FunctionLambda The lambda function type.
       * \tparam GradientLambda The gradient lambda function type.
       * \param[in] function Function lambda.
       * \param[in] gradient Gradient lambda function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool solve(FunctionLambda && function, GradientLambda && gradient, Vector const & x_ini, Vector & x_sol)
      {
        #define CMD "Optimist::Optimizer::solve(...): "

        static_assert(DerivedSolver::RequiresFunction,
          CMD "the solver requires the function.");
        static_assert(DerivedSolver::RequiresFirstDerivative,
          CMD "the solver requires the first derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<GradientLambda>(gradient),
          x_ini, x_sol);

        #undef CMD
      }

      /**
       * Solve the root-finding problem given the function, and its gradient and Hessian.
       * \tparam FunctionLambda The lambda function type.
       * \tparam GradientLambda The gradient lambda function type.
       * \tparam HessianLambda The hessian lambda function type.
       * \param[in] function Function lambda.
       * \param[in] gradient Gradient lambda function.
       * \param[in] hessian Hessian lambda function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename GradientLambda, typename HessianLambda>
      bool solve(FunctionLambda && function, GradientLambda && gradient, HessianLambda && hessian,
         Vector const & x_ini, Vector & x_sol)
      {
        #define CMD "Optimist::Optimizer::solve(...): "

        static_assert(DerivedSolver::RequiresFunction,
          CMD "the solver requires the function.");
        static_assert(DerivedSolver::RequiresFirstDerivative,
          CMD "the solver requires the first derivative.");
        static_assert(DerivedSolver::RequiresSecondDerivative,
          CMD "the solver requires the second derivative.");
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<GradientLambda>(gradient),
          std::forward<HessianLambda>(hessian),
          x_ini, x_sol);

        #undef CMD
      }

    }; // class Optimizer

    /**
     * \brief Class container for the scalar optimizer.
     *
     * \includedoc docs/markdown/ScalarOptimizer.md
     *
     * \tparam Real Scalar number type.
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename Real, typename DerivedSolver>
    class Optimizer<Real, 1, DerivedSolver> : public SolverBase<Real, 1, 1, DerivedSolver>
    {
    public:
      friend class SolverBase<Real, 1, 1, Optimizer<Real, 1, DerivedSolver>>;

      static constexpr bool IsRootFinder{false};
      static constexpr bool IsOptimizer{true};

      OPTIMIST_BASIC_CONSTANTS(Real)

      /**
       * Class constructor for the scalar optimizer.
       */
      Optimizer<Real, 1, DerivedSolver>() {}

      /**
       * Get the solver name.
       * \return The solver name.
       */
      constexpr std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

    }; // class Optimizer

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_HH
