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
     * \tparam Input Solver input type.
     * \tparam Output Solver output type.
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename Input, typename Output, typename DerivedSolver>
    requires (TypeTrait<Input>::IsScalar && TypeTrait<Output>::IsScalar) ||
      ((TypeTrait<Input>::IsEigen && TypeTrait<Output>::IsEigen) &&
      ((TypeTrait<Input>::IsFixed && TypeTrait<Output>::IsFixed) ||
      (TypeTrait<Input>::IsDynamic && TypeTrait<Output>::IsDynamic) ||
      (TypeTrait<Input>::IsSparse && TypeTrait<Output>::IsSparse)))
    class RootFinder : public SolverBase<Input, Output, DerivedSolver>
    {
    public:
      friend class SolverBase<Input, Output, DerivedSolver>;

      static constexpr bool IsRootFinder{true};
      static constexpr bool IsOptimizer{false};

      // Input and output types
      using Scalar = typename TypeTrait<Input>::Scalar;
      using Vector = Input;

      // Derivative types
      using Matrix = typename SolverBase<Input, Output, DerivedSolver>::FirstDerivative;
      using Tensor = typename SolverBase<Input, Output, DerivedSolver>::SecondDerivative;

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
       * \tparam JacobianLambda The Jacobian lambda function type.
       * \param[in] jacobian Jacobian lambda function.
       * \param[in] x Input point.
       * \param[out] out Jacobian value.
       * \return The boolean flag for successful evaluation.
       */
      template <typename JacobianLambda>
      bool evaluate_jacobian(JacobianLambda && jacobian, Vector const & x, Matrix & out)
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
      bool evaluate_hessian(HessianLambda && hessian, Vector const & x, Matrix & out)
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
      bool solve(FunctionLambda && function, Vector const & x_ini, Vector & x_sol)
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
      bool solve(FunctionLambda && function, JacobianLambda && jacobian, Vector const & x_ini,
        Vector & x_sol)
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
        Vector const & x_ini, Vector & x_sol)
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * \brief Class container for the scalar scalar root-finder.
     *
     * \includedoc docs/markdown/RootFinder.md
     *
     * \tparam Scalar Scalar number type.
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename Scalar, typename DerivedSolver>
    requires TypeTrait<Scalar>::IsScalar
    class RootFinder<Scalar, Scalar, DerivedSolver> : public SolverBase<Scalar, Scalar, DerivedSolver>
    {
    public:
      friend class SolverBase<Scalar, Scalar, DerivedSolver>;

      static constexpr bool IsRootFinder{true};
      static constexpr bool IsOptimizer{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the scalar root-finder.
       */
      RootFinder<Scalar, 1, DerivedSolver>() {}

      /**
       * Get the solver name.
       * \return The solver name.
       */
      constexpr std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

      /**
       * Solve the root-finding problem given the function, and without derivatives.
       * \tparam FunctionLambda The lambda function type.
       * \param[in] function Function lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve(FunctionLambda && function, Scalar x_ini, Scalar & x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          x_ini, x_sol);
      }

      /**
       * Solve the root-finding problem given the function, and its first derivative.
       * \tparam FunctionLambda The lambda function type.
       * \tparam FirstDerivativeLambda The first derivative lambda type.
       * \param[in] function Function lambda.
       * \param[in] first_derivative First derivative lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
       template <typename FunctionLambda, typename FirstDerivativeLambda>
      bool solve(FunctionLambda && function, FirstDerivativeLambda && first_derivative, Scalar x_ini,
        Scalar & x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<FirstDerivativeLambda>(first_derivative),
          x_ini, x_sol);
      }

      /**
       * Solve the root-finding problem given the function, and its first and second derivatives.
       * \tparam FunctionLambda The lambda function type.
       * \tparam FirstDerivativeLambda The first derivative lambda type.
       * \tparam SecondDerivativeLambda The second derivative lambda type.
       * \param[in] function Function lambda.
       * \param[in] first_derivative First derivative lambda.
       * \param[in] second_derivate Second derivative lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename FirstDerivativeLambda, typename SecondDerivativeLambda>
      bool solve(FunctionLambda && function, FirstDerivativeLambda && first_derivative, SecondDerivativeLambda
        && second_derivate, Scalar const x_ini, Scalar & x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<FirstDerivativeLambda>(first_derivative),
          std::forward<SecondDerivativeLambda>(second_derivate),
          x_ini, x_sol);
      }

    }; // class RootFinder

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_HH
