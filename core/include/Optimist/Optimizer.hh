/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_OPTIMIZER_HH
#define OPTIMIST_OPTIMIZER_HH

#include "Optimist/SolverBase.hh"

namespace Optimist {

  /**
   * \brief Namespace for multi-dimensional optimization algorithms.
   */
  namespace Optimizer {

    /**
     * \brief Class container for the multi-dimensional optimizer.
     *
     * \includedoc docs/markdown/Optimizer.md
     *
     * \tparam T Input type (scalar or Eigen vector).
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename T, typename DerivedSolver>
      requires(TypeTrait<T>::IsScalar || TypeTrait<T>::IsEigen) &&
              (!TypeTrait<T>::IsFixed || TypeTrait<T>::Dimension > 0)
    class Optimizer
        : public SolverBase<T, typename TypeTrait<T>::Scalar, DerivedSolver> {
     public:
      friend class SolverBase<T, typename TypeTrait<T>::Scalar, DerivedSolver>;

      static constexpr bool IsRootFinder{false};
      static constexpr bool IsOptimizer{true};

      // Input and output types
      using Scalar = typename TypeTrait<T>::Scalar;
      using Input  = T;
      using Output = Scalar;

      // Derivative types
      using typename SolverBase<T, Scalar, DerivedSolver>::FirstDerivative;
      using typename SolverBase<T, Scalar, DerivedSolver>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the multi-dimensional optimizer.
       */
      Optimizer() {}

      /**
       * Get the solver name.
       * \return The solver name.
       */
      constexpr std::string name() const {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

      /**
       * Get the number of gradient evaluations.
       * \return The number of gradient evaluations.
       */
      Integer gradient_evaluations() const {
        return this->first_derivative_evaluations();
      }

      /**
       * Get the number of maximum allowed gradient evaluations.
       * \return The number of maximum allowed gradient evaluations.
       */
      Integer max_gradient_evaluations() const {
        return this->max_first_derivative_evaluations();
      }

      /**
       * Set the number of maximum allowed gradient evaluations.
       * \param[in] t_gradient_evaluations The number of maximum allowed
       * gradient evaluations.
       */
      void max_gradient_evaluations(const Integer t_gradient_evaluations) {
        this->max_first_derivative_evaluations(t_gradient_evaluations);
      }

      /**
       * Get the number of hessian evaluations.
       * \return The number of hessian evaluations.
       */
      Integer hessian_evaluations() const {
        return this->second_derivative_evaluations();
      }

      /**
       * Get the number of maximum allowed hessian evaluations.
       * \return The number of maximum allowed hessian evaluations.
       */
      Integer max_hessian_evaluations() const {
        return this->max_second_derivative_evaluations();
      }

      /**
       * Set the number of maximum allowed hessian evaluations.
       * \param[in] t_hessian_evaluations The number of maximum allowed hessian
       * evaluations.
       */
      void max_hessian_evaluations(const Integer t_hessian_evaluations) {
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
      bool evaluate_gradient(GradientLambda &&gradient,
                             const Input &x,
                             FirstDerivative &out) {
        return this->evaluate_first_derivative(
            std::forward<GradientLambda>(gradient),
            x,
            out);
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
      bool evaluate_hessian(HessianLambda &&hessian,
                            const Input &x,
                            SecondDerivative &out) {
        return this->evaluate_second_derivative(
            std::forward<HessianLambda>(hessian),
            x,
            out);
      }

    };  // class Optimizer

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_HH
