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

#ifndef OPTIMIST_ROOTFINDER_HH
#define OPTIMIST_ROOTFINDER_HH

#include "Optimist/SolverBase.hh"

namespace Optimist {

  /**
   * \brief Namespace for multi-dimensional root-finding algorithms.
   */
  namespace RootFinder {

    /**
     * \brief Class container for the multi-dimensional root finder.
     *
     * \includedoc docs/markdown/RootFinder.md
     *
     * \tparam T Input type (scalar or Eigen vector).
     * \tparam DerivedSolver Derived solver class.
     * \tparam U Output type, equal to the input type for square systems.
     * Solvers supporting non-square systems (in the least-squares sense) may
     * use a different output dimension.
     */
    template <typename T, typename DerivedSolver, typename U = T>
      requires(TypeTrait<T>::IsScalar || TypeTrait<T>::IsEigen) &&
              (!TypeTrait<T>::IsFixed || TypeTrait<T>::Dimension > 0)
    class RootFinder : public SolverBase<T, U, DerivedSolver> {
     public:
      friend class SolverBase<T, U, DerivedSolver>;

      static constexpr bool IsRootFinder{true};
      static constexpr bool IsOptimizer{false};

      // Input and output types
      using Scalar = typename TypeTrait<T>::Scalar;
      using Input  = T;
      using Output = U;

      // Derivative types
      using typename SolverBase<T, U, DerivedSolver>::FirstDerivative;
      using typename SolverBase<T, U, DerivedSolver>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the multi-dimensional root finder.
       */
      RootFinder() {}

      /**
       * Get the solver name.
       * \return The solver name.
       */
      constexpr std::string name() const {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

      /**
       * Get the number of Jacobian evaluations.
       * \return The number of Jacobian evaluations.
       */
      Integer jacobian_evaluations() const {
        return this->first_derivative_evaluations();
      }

      /**
       * Get the number of maximum allowed Jacobian evaluations.
       * \return The number of maximum allowed Jacobian evaluations.
       */
      Integer max_jacobian_evaluations() const {
        return this->max_first_derivative_evaluations();
      }

      /**
       * Set the number of maximum allowed Jacobian evaluations.
       * \param[in] t_jacobian_evaluations The number of maximum allowed
       * Jacobian evaluations.
       */
      void max_jacobian_evaluations(const Integer t_jacobian_evaluations) {
        this->max_first_derivative_evaluations(t_jacobian_evaluations);
      }

      /**
       * Get the number of Hessian evaluations.
       * \return The number of Hessian evaluations.
       */
      Integer hessian_evaluations() const {
        return this->second_derivative_evaluations();
      }

      /**
       * Get the number of maximum allowed Hessian evaluations.
       * \return The number of maximum allowed Hessian evaluations.
       */
      Integer max_hessian_evaluations() const {
        return this->max_second_derivative_evaluations();
      }

      /**
       * Set the number of maximum allowed Hessian evaluations.
       * \param[in] t_hessian_evaluations The number of maximum allowed Hessian
       * evaluations.
       */
      void max_hessian_evaluations(const Integer t_hessian_evaluations) {
        this->max_second_derivative_evaluations(t_hessian_evaluations);
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
      bool evaluate_jacobian(JacobianLambda &&jacobian,
                             const Input &x,
                             FirstDerivative &out) {
        return this->evaluate_first_derivative(
            std::forward<JacobianLambda>(jacobian),
            x,
            out);
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
      bool evaluate_hessian(HessianLambda &&hessian,
                            const Input &x,
                            SecondDerivative &out) {
        return this->evaluate_second_derivative(
            std::forward<HessianLambda>(hessian),
            x,
            out);
      }

    };  // class RootFinder

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_HH
