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

#ifndef OPTIMIST_SOLVER_HH
#define OPTIMIST_SOLVER_HH

#include "Optimist.hh"
#include "Optimist/Function.hh"
#include "Optimist/Utils.hh"

namespace Optimist {

  /**
   * \brief Class container for the generic root-finding/optimization problem
   * solver.
   *
   * \includedoc docs/markdown/SolverBase.md
   *
   * \tparam Input Solver input type.
   * \tparam Output Solver output type.
   * \tparam DerivedSolver Derived solver class.
   */
  template <typename Input, typename Output, typename DerivedSolver>
    requires std::is_same<typename TypeTrait<Input>::Scalar,
                          typename TypeTrait<Output>::Scalar>::value &&
             (TypeTrait<Input>::IsScalar || TypeTrait<Input>::IsEigen) &&
             (TypeTrait<Output>::IsScalar || TypeTrait<Output>::IsEigen) &&
             (!TypeTrait<Input>::IsFixed || TypeTrait<Input>::Dimension > 0) &&
             (!TypeTrait<Output>::IsFixed ||
              TypeTrait<Output>::Dimension > 0) &&
             (!(TypeTrait<Input>::IsEigen && TypeTrait<Output>::IsEigen) ||
              (TypeTrait<Input>::IsFixed && TypeTrait<Output>::IsFixed) ||
              (TypeTrait<Input>::IsDynamic && TypeTrait<Output>::IsDynamic) ||
              (TypeTrait<Input>::IsSparse && TypeTrait<Output>::IsSparse))
  class SolverBase {
   public:
    // Input and output types
    using InputTrait  = TypeTrait<Input>;
    using OutputTrait = TypeTrait<Output>;
    using Scalar      = typename InputTrait::Scalar;

    // Input and output must be non-zero dimensional
    static_assert(InputTrait::Dimension != 0,
                  "Input must be non-zero dimensional");
    static_assert(OutputTrait::Dimension != 0,
                  "Output must be non-zero dimensional");

    // Input and output must have the same scalar type
    static_assert(std::is_same<typename InputTrait::Scalar,
                               typename OutputTrait::Scalar>::value,
                  "Input and output scalar types must be the same.");

    // If both input and output are eigen types they be both fixed-size,
    // dynamic-size, or sparse
    static_assert(!(InputTrait::IsEigen && OutputTrait::IsEigen) ||
                      (InputTrait::IsFixed && OutputTrait::IsFixed) ||
                      (InputTrait::IsDynamic && OutputTrait::IsDynamic) ||
                      (InputTrait::IsSparse && OutputTrait::IsSparse),
                  "Input and output Eigen types must be both fixed-size, "
                  "dynamic-size, or sparse.");

    // Derivative types
    using FirstDerivative =
        FunctionBase<Input, Output, DerivedSolver>::FirstDerivative;
    using SecondDerivative =
        FunctionBase<Input, Output, DerivedSolver>::SecondDerivative;

    OPTIMIST_BASIC_CONSTANTS(Scalar)

   protected:
    // Bounds (may not be used)
    Input m_lower_bound; /**< Lower bound. */
    Input m_upper_bound; /**< Upper bound. */

    // Evaluations
    Integer m_function_evaluations{0}; /**< Function evaluations. */
    Integer m_first_derivative_evaluations{
      0};                              /**< First derivative evaluations. */
    Integer m_second_derivative_evaluations{
      0};                              /**< Second derivative evaluations. */

    // Maximum allowed evaluations
    Integer m_max_function_evaluations{
      1000}; /**< Maximum allowed function evaluations. */
    Integer m_max_first_derivative_evaluations{
      1000}; /**< Maximum allowed first derivative evaluations. */
    Integer m_max_second_derivative_evaluations{
      1000}; /**< Maximum allowed second derivative evaluations. */

    // Iterations and relaxations
    Integer m_iterations{0};       /**< Algorithm iterations. */
    Integer m_max_iterations{100}; /**< Maximum allowed algorithm iterations. */
    Scalar m_alpha{0.8};           /**< Relaxation factor \f$ \alpha \f$. */
    Integer m_relaxations{0};      /**< Algorithm relaxations. */
    Integer m_max_relaxations{
      10}; /**< Maximum allowed algorithm relaxations. */

    // Settings
    Scalar m_tolerance{
      SQRT_EPSILON}; /**< Solver tolerance \f$ \epsilon \f$ for convergence. */
    bool m_verbose{false};               /**< Verbose mode boolean flag. */
    bool m_damped{true};                 /**< Damped mode boolean flag. */
    std::ostream *m_ostream{&std::cout}; /**< Output stream for verbose mode. */

    // Convergence output flag and trace
    std::string m_task{"Undefined"}; /**< Task name. */
    bool m_converged{false};         /**< Convergence boolean flag. */

   public:
    /**
     * Class constructor for the nonlinear solver.
     */
    SolverBase() {
      this->reset_bounds();
    }

    /**
     * Reset lower and upper bounds to default values.
     * \param[in] n Input dimension for dynamic-size types.
     */
    void reset_bounds(const Integer n = InputTrait::IsFixed
                                            ? InputTrait::Dimension
                                            : 0) {
#define CMD "Optimist::Solver::reset_bounds(...): "

      if constexpr (InputTrait::IsScalar) {
        this->m_lower_bound = -INFTY;
        this->m_upper_bound = +INFTY;
      } else if constexpr (InputTrait::IsFixed) {
        this->m_lower_bound.setConstant(-INFTY);
        this->m_upper_bound.setConstant(+INFTY);
      } else if constexpr (InputTrait::IsDynamic) {
        this->m_lower_bound.resize(n);
        this->m_lower_bound.setConstant(-INFTY);
        this->m_upper_bound.resize(n);
        this->m_upper_bound.setConstant(+INFTY);
      } else if constexpr (InputTrait::IsSparse) {
        this->m_lower_bound.resize(n);
        this->m_lower_bound.reserve(n);
        this->m_upper_bound.resize(n);
        this->m_upper_bound.reserve(n);
        for (Integer i{0}; i < n; ++i) {
          this->m_lower_bound.coeffRef(i) = -INFTY;
          this->m_upper_bound.coeffRef(i) = +INFTY;
        }
      } else {
        OPTIMIST_ERROR(CMD "unsupported input type for bounds reset.");
      }

#undef CMD
    }

    /**
     * Get the lower bound.
     * \return The lower bound.
     */
    const Input &lower_bound() const {
      return this->m_lower_bound;
    }

    /**
     * Set the lower bound.
     * \param[in] t_lower_bound The lower bound.
     */
    void lower_bound(const Input &t_lower_bound) {
      OPTIMIST_ASSERT(valid_bounds(t_lower_bound, this->m_upper_bound),
                      "Optimist::Solver::lower_bound(...): invalid or "
                      "degenerate bounds detected.");
      this->m_lower_bound = t_lower_bound;
    }

    /**
     * Get the upper bound.
     * \return The upper bound.
     */
    const Input &upper_bound() const {
      return this->m_upper_bound;
    }

    /**
     * Set the upper bound.
     * \param[in] t_upper_bound The upper bound.
     */
    void upper_bound(const Input &t_upper_bound) {
      OPTIMIST_ASSERT(valid_bounds(this->m_lower_bound, t_upper_bound),
                      "Optimist::Solver::upper_bound(...): invalid or "
                      "degenerate bounds detected.");
      this->m_upper_bound = t_upper_bound;
    }

    /**
     * Set the bounds.
     * \param[in] t_lower_bound The lower bound.
     * \param[in] t_upper_bound The upper bound.
     */
    void bounds(const Input &t_lower_bound, const Input &t_upper_bound) {
      OPTIMIST_ASSERT(
          valid_bounds(t_lower_bound, t_upper_bound) &&
              (InputTrait::IsScalar ||
               Integer(size_of(t_lower_bound)) ==
                   Integer(size_of(t_upper_bound))),
          "Optimist::Solver::bounds(...): invalid or degenerate bounds "
          "detected.");
      this->m_lower_bound = t_lower_bound;
      this->m_upper_bound = t_upper_bound;
    }

    /**
     * Check whether the solver has finite bounds on every component.
     * 
eturn True if the bounds are set and finite.
     */
    bool has_finite_bounds() const {
      if constexpr (InputTrait::IsScalar) {
        return std::isfinite(this->m_lower_bound) &&
               std::isfinite(this->m_upper_bound);
      } else {
        using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        if (this->m_lower_bound.size() == 0 ||
            this->m_lower_bound.size() != this->m_upper_bound.size()) {
          return false;
        }
        return Dense(this->m_lower_bound).allFinite() &&
               Dense(this->m_upper_bound).allFinite();
      }
    }

   protected:
    /**
     * Size of an input (1 for scalars).
     */
    static Integer size_of(const Input &x) {
      if constexpr (InputTrait::IsScalar) {
        return 1;
      } else {
        return static_cast<Integer>(x.size());
      }
    }

    /**
     * Check that the lower bound is strictly below the upper bound on every
     * component. Bounds of different sizes (e.g., not yet initialized
     * dynamic-size bounds) are not compared.
     */
    static bool valid_bounds(const Input &t_lower, const Input &t_upper) {
      if constexpr (InputTrait::IsScalar) {
        return !(t_upper <= t_lower);
      } else {
        using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        if (t_lower.size() != t_upper.size() || t_lower.size() == 0) {
          return true;
        }
        const Dense lower(t_lower), upper(t_upper);
        for (Eigen::Index i{0}; i < lower.size(); ++i) {
          if (!(upper(i) > lower(i))) {
            return false;
          }
        }
        return true;
      }
    }

   public:

    /**
     * Get the input dimension of the function.
     * \return The input dimension of the function.
     */
    constexpr Integer input_dimension() const {
      return InputTrait::Dimension;
    }

    /**
     * Get the output dimension of the function.
     * \return The output dimension of the function.
     */
    constexpr Integer output_dimension() const {
      return OutputTrait::Dimension;
    }

    /**
     * Get the number of function evaluations.
     * \return The number of function evaluations.
     */
    Integer function_evaluations() const {
      return this->m_function_evaluations;
    }

    /**
     * Set the number of maximum allowed function evaluations.
     * \param[in] t_max_function_evaluations The number of maximum allowed
     * function evaluations.
     */
    void max_function_evaluations(const Integer t_max_function_evaluations) {
      OPTIMIST_ASSERT(
          t_max_function_evaluations > 0,
          "Optimist::Solver::max_function_evaluations(...): invalid "
          "input detected.");
      this->m_max_function_evaluations = t_max_function_evaluations;
    }

    /**
     * Get the number of maximum allowed function evaluations.
     * \return The number of maximum allowed function evaluations.
     */
    Integer max_function_evaluations() const {
      return this->m_max_function_evaluations;
    }

   protected:
    /**
     * Get the number of function derivative evaluations.
     * \return The number of function derivative evaluations.
     */
    Integer first_derivative_evaluations() const {
      return this->m_first_derivative_evaluations;
    }

    /**
     * Get the number of maximum allowed first derivative evaluations.
     * \return The number of maximum allowed first derivative evaluations.
     */
    Integer max_first_derivative_evaluations() const {
      return this->m_max_first_derivative_evaluations;
    }

    /**
     * Set the number of maximum allowed first derivative evaluations.
     * \param[in] t_first_derivative_evaluations The number of maximum allowed
     * first derivative evaluations.
     */
    void max_first_derivative_evaluations(
        const Integer t_first_derivative_evaluations) {
      OPTIMIST_ASSERT(
          t_first_derivative_evaluations > 0,
          "Optimist::Solver::max_first_derivative_evaluations(...): "
          "invalid input detected.");
      this->m_max_first_derivative_evaluations = t_first_derivative_evaluations;
    }

    /**
     * Get the number of second derivative evaluations.
     * \return The number of second derivative evaluations.
     */
    Integer second_derivative_evaluations() const {
      return this->m_second_derivative_evaluations;
    }

    /**
     * Get the number of maximum allowed second derivative evaluations.
     * \return The number of maximum allowed second derivative evaluations.
     */
    Integer max_second_derivative_evaluations() const {
      return this->m_max_second_derivative_evaluations;
    }

    /**
     * Set the number of maximum allowed second derivative evaluations.
     * \param[in] t_second_derivative_evaluations The number of maximum allowed
     * second derivative evaluations.
     */
    void max_second_derivative_evaluations(
        const Integer t_second_derivative_evaluations) {
      OPTIMIST_ASSERT(
          t_second_derivative_evaluations > 0,
          "Optimist::Solver::max_second_derivative_evaluations(...): "
          "invalid input detected.");
      this->m_max_second_derivative_evaluations =
          t_second_derivative_evaluations;
    }

   public:
    /**
     * Get the number of algorithm iterations.
     */
    Integer iterations() const {
      return this->m_iterations;
    }

    /**
     * Get the number of maximum allowed iterations.
     * \return The number of maximum allowed iterations.
     */
    Integer max_iterations() const {
      return this->m_max_iterations;
    }

    /**
     * Set the number of maximum allowed iterations.
     * \param[in] t_max_iterations The number of maximum allowed iterations.
     */
    void max_iterations(const Integer t_max_iterations) {
      OPTIMIST_ASSERT(
          t_max_iterations > 0,
          "Optimist::Solver::max_iterations(...): invalid input detected.");
      this->m_max_iterations = t_max_iterations;
    }

    /**
     * Get relaxation factor \f$ \alpha \f$.
     * \return The relaxation factor \f$ \alpha \f$.
     */
    Scalar alpha() const {
      return this->m_alpha;
    }

    /**
     * Set relaxation factor \f$ \alpha \f$.
     * \param[in] t_alpha The relaxation factor \f$ \alpha \f$.
     */
    void alpha(Scalar t_alpha) {
      OPTIMIST_ASSERT(!std::isnan(t_alpha) && std::isfinite(t_alpha) &&
                          0.0 <= t_alpha && t_alpha <= 1.0,
                      "Optimist::Solver::alpha(...): invalid input detected.");
      this->m_alpha = t_alpha;
    }

    /**
     * Get the number of algorithm relaxations.
     * \return The number of algorithm relaxations.
     */
    Integer relaxations() const {
      return this->m_relaxations;
    }

    /**
     * Get the number of maximum allowed relaxations.
     * \return The number of maximum allowed relaxations.
     */
    Integer max_relaxations() const {
      return this->m_max_relaxations;
    }

    /**
     * Set the number of maximum allowed relaxations.
     * \param[in] t_max_relaxations The number of maximum allowed relaxations.
     */
    void max_relaxations(const Integer t_max_relaxations) {
      OPTIMIST_ASSERT(
          t_max_relaxations > 0,
          "Optimist::Solver::max_relaxations(...): invalid input detected.");
      this->m_max_relaxations = t_max_relaxations;
    }

    /**
     * Get the tolerance \f$ \epsilon \f$.
     * \return The tolerance \f$ \epsilon \f$.
     */
    Scalar tolerance() const {
      return this->m_tolerance;
    }

    /**
     * Set the tolerance \f$ \epsilon \f$ for which the nonlinear solver stops,
     * i.e., \f$ \left\|
     * \mathbf{F}(\mathbf{x}) \right\|_{2} < \epsilon \f$.
     * \param[in] t_tolerance The tolerance \f$ \epsilon \f$.
     */
    void tolerance(Scalar t_tolerance) {
      OPTIMIST_ASSERT(
          !std::isnan(t_tolerance) && std::isfinite(t_tolerance) &&
              t_tolerance > 0,
          "Optimist::Solver::tolerance(...): invalid input detected.");
      this->m_tolerance = t_tolerance;
    }

    /**
     * Set the verbose mode boolean flag.
     * \param[in] t_verbose The verbose mode boolean flag.
     */
    void verbose_mode(bool t_verbose) {
      this->m_verbose = t_verbose;
    }

    /**
     * Get the verbose mode boolean flag.
     * \return The verbose mode boolean flag.
     */
    bool verbose_mode() const {
      return this->m_verbose;
    }

    /**
     * Enable solver's verbose mode.
     */
    void enable_verbose_mode() {
      this->m_verbose = true;
    }

    /**
     * Disable solver's verbose mode.
     */
    void disable_verbose_mode() {
      this->m_verbose = false;
    }

    /**
     * Set the damped mode boolean flag.
     * \param[in] t_damped The damped mode boolean flag.
     */
    void damped_mode(bool t_damped) {
      this->m_damped = t_damped;
    }

    /**
     * Get the damped mode boolean flag.
     * \return The damped mode boolean flag.
     */
    bool damped_mode() const {
      return this->m_damped;
    }

    /**
     * Enable solver's damped mode.
     */
    void enable_damped_mode() {
      this->m_damped = true;
    }

    /**
     * Disable solver's damped mode.
     */
    void disable_damped_mode() {
      this->m_damped = false;
    }

    /**
     * Get the task name.
     * \return The task name.
     */
    std::string task() const {
      return this->m_task;
    }

    /**
     * Set the task name.
     * \param[in] t_task The task name.
     */
    void task(std::string t_task) {
      this->m_task = t_task;
    }

    /**
     * Get the convergence boolean flag.
     * \return The convergence boolean flag.
     */
    bool converged() const {
      return this->m_converged;
    }

    /**
     * Get the output stream for verbose mode.
     * \return The output stream for verbose mode.
     */
    std::ostream &ostream() const {
      return *this->m_ostream;
    }

    /**
     * Set the output stream for verbose mode.
     * \param[in] t_ostream The output stream for verbose mode.
     */
    void ostream(std::ostream &t_ostream) {
      this->m_ostream = &t_ostream;
    }

    /**
     * Solve the root-finding/optimization problem given the function, and
     * without derivatives.
     * \tparam FunctionLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return True if the problem is solved, false otherwise.
     */
    template <typename FunctionLambda>
    bool solve(FunctionLambda &&function, const Input &x_ini, Input &x_sol) {
#define CMD "Optimist::Solver::solve(...): "

      static_assert(!DerivedSolver::RequiresFirstDerivative,
                    CMD "the solver requires the first derivative.");
      return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          x_ini,
          x_sol);

#undef CMD
    }

    /**
     * Solve the root-finding/optimization problem given the function, and its
     * first derivative.
     * \tparam FunctionLambda First derivative lambda type.
     * \tparam FirstDerivativeLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] first_derivative First derivative lambda
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return True if the problem is solved, false otherwise.
     */
    template <typename FunctionLambda, typename FirstDerivativeLambda>
    bool solve(FunctionLambda &&function,
               FirstDerivativeLambda &&first_derivative,
               const Input &x_ini,
               Input &x_sol) {
#define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::RequiresFunction,
                    CMD "the solver requires a function.");
      static_assert(DerivedSolver::RequiresFirstDerivative,
                    CMD "the solver requires the first derivative.");
      static_assert(!DerivedSolver::RequiresSecondDerivative,
                    CMD "the solver requires the second derivative.");
      return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<FirstDerivativeLambda>(first_derivative),
          x_ini,
          x_sol);

#undef CMD
    }

    /**
     * Solve the root-finding/optimization problem given the function, and its
     * first and second derivatives.
     * \tparam FunctionLambda Function lambda type.
     * \tparam FirstDerivativeLambda First derivative lambda type.
     * \tparam SecondDerivativeLambda Second derivative lambda type.
     * \param[in] function Function lambda.
     * \param[in] first_derivative First derivative lambda.
     * \param[in] second_derivative The second derivative lambda.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return True if the problem is solved, false otherwise.
     */
    template <typename FunctionLambda,
              typename FirstDerivativeLambda,
              typename SecondDerivativeLambda>
    bool solve(FunctionLambda &&function,
               FirstDerivativeLambda &&first_derivative,
               SecondDerivativeLambda &&second_derivative,
               const Input &x_ini,
               Input &x_sol) {
#define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::RequiresFunction,
                    CMD "the solver requires the function.");
      static_assert(DerivedSolver::RequiresFirstDerivative,
                    CMD "the solver requires the first derivative.");
      static_assert(DerivedSolver::RequiresSecondDerivative,
                    CMD "the solver requires the second derivative.");
      return static_cast<DerivedSolver *>(this)->solve_impl(
          std::forward<FunctionLambda>(function),
          std::forward<FirstDerivativeLambda>(first_derivative),
          std::forward<SecondDerivativeLambda>(second_derivative),
          x_ini,
          x_sol);

#undef CMD
    }

    /**
     * Solve the root-finding problem \f$ \mathbf{f}(\mathbf{x}) = \mathbf{0}
     * \f$ given a function class.
     *
     * - Root finders solve the system directly.
     * - Optimizers minimize the least-squares merit function \f$
     *   \frac{1}{2}\|\mathbf{f}(\mathbf{x})\|_2^2 \f$, whose global minima are
     *   the roots of \f$ \mathbf{f} \f$.
     *
     * \tparam FunctionInput The function input type.
     * \tparam FunctionOutput The function output type.
     * \tparam DerivedFunction Derived function class.
     * \param[in] function Function class.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return The convergence boolean flag.
     */
    template <typename FunctionInput,
              typename FunctionOutput,
              typename DerivedFunction>
      requires(TypeTrait<FunctionInput>::Dimension == InputTrait::Dimension)
    bool rootfind(
        const FunctionBase<FunctionInput, FunctionOutput, DerivedFunction>
            &function,
        const Input &x_ini,
        Input &x_sol) {
#define CMD "Optimist::Solver::rootfind(...): "

      using FunctionOutputTrait = TypeTrait<FunctionOutput>;

      static_assert(
          !(InputTrait::IsScalar && DerivedSolver::IsOptimizer),
          CMD
          "one-dimensional optimizers do not support root-finding problems.");

      if constexpr (DerivedSolver::IsOptimizer) {
        return this->template solve<Problem::LEAST_SQUARES>(function,
                                                            x_ini,
                                                            x_sol);
      } else {
        static_assert(OutputTrait::Dimension == FunctionOutputTrait::Dimension,
                      CMD
                      "root finders require a function whose output dimension "
                      "equals the input dimension.");
        return this->template solve<Problem::ROOT>(function, x_ini, x_sol);
      }

#undef CMD
    }

    /**
     * Solve the optimization problem given a function class.
     *
     * The objective \f$ \phi(\mathbf{x}) \f$ is selected as follows:
     * - cost functions \f$ f : \mathbb{R}^n \to \mathbb{R} \f$ handed to an
     *   optimizer are minimized directly, \f$ \phi = f \f$;
     * - vector-valued functions are turned into the least-squares merit
     *   function \f$ \phi = \frac{1}{2}\|\mathbf{f}\|_2^2 \f$, and so are
     *   scalar functions handed to a one-dimensional root finder;
     * - root finders locate the stationary points of \f$ \phi \f$ by solving
     *   \f$ \nabla\phi(\mathbf{x}) = \mathbf{0} \f$, using \f$
     *   \nabla^2\phi \f$ as the Jacobian.
     *
     * \tparam FunctionInput The function input type.
     * \tparam FunctionOutput The function output type.
     * \tparam DerivedFunction Derived function class.
     * \param[in] function Function class.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return The convergence boolean flag.
     */
    template <typename FunctionInput,
              typename FunctionOutput,
              typename DerivedFunction>
      requires(TypeTrait<FunctionInput>::Dimension == InputTrait::Dimension)
    bool optimize(
        const FunctionBase<FunctionInput, FunctionOutput, DerivedFunction>
            &function,
        const Input &x_ini,
        Input &x_sol) {
#define CMD "Optimist::Solver::optimize(...): "

      using FunctionInputTrait  = TypeTrait<FunctionInput>;
      using FunctionOutputTrait = TypeTrait<FunctionOutput>;

      if constexpr (DerivedSolver::IsOptimizer) {
        if constexpr (FunctionOutputTrait::IsScalar) {
          return this->template solve<Problem::COST>(function, x_ini, x_sol);
        } else {
          return this->template solve<Problem::LEAST_SQUARES>(function,
                                                              x_ini,
                                                              x_sol);
        }
      } else {
        static_assert(!DerivedSolver::RequiresSecondDerivative,
                      CMD
                      "root finders that require second derivatives cannot "
                      "solve optimization problems (third derivatives of the "
                      "function would be needed).");
        if constexpr (FunctionInputTrait::IsEigen &&
                      FunctionOutputTrait::IsScalar) {
          return this->template solve<Problem::STATIONARY_COST>(function,
                                                                x_ini,
                                                                x_sol);
        } else {
          static_assert(
              OutputTrait::Dimension == FunctionOutputTrait::Dimension,
              CMD
              "root finders require a function whose output dimension equals "
              "the input dimension.");
          return this->template solve<Problem::STATIONARY_LEAST_SQUARES>(
              function,
              x_ini,
              x_sol);
        }
      }

#undef CMD
    }

    /**
     * Get the solver name.
     * \return The solver name.
     */
    constexpr std::string name() const {
      return static_cast<const DerivedSolver *>(this)->name_impl();
    };

   protected:
    /**
     * Kind of problem handed to the solver implementation.
     */
    using Problem = enum class Problem : Integer {
      ROOT          = 0, /**< Solve \f$ \mathbf{f} = \mathbf{0} \f$. */
      COST          = 1, /**< Minimize \f$ f \f$. */
      LEAST_SQUARES = 2, /**< Minimize \f$ \frac{1}{2}\|\mathbf{f}\|^2 \f$. */
      STATIONARY_COST = 3, /**< Solve \f$ \nabla f = \mathbf{0} \f$. */
      STATIONARY_LEAST_SQUARES =
          4 /**< Solve \f$ \mathbf{J}^\top\mathbf{f} = \mathbf{0} \f$. */
    };

    /**
     * Compute \f$ \mathbf{J}^\top\mathbf{f} \f$, the gradient of the
     * least-squares merit function \f$ \frac{1}{2}\|\mathbf{f}\|^2 \f$.
     */
    template <typename FunctionValue, typename FunctionJacobian, typename Out>
    static void least_squares_gradient(const FunctionValue &f,
                                       const FunctionJacobian &J,
                                       Out &out) {
      if constexpr (TypeTrait<FunctionValue>::IsScalar &&
                    TypeTrait<FunctionJacobian>::IsScalar) {
        out = J * f;
      } else if constexpr (TypeTrait<FunctionValue>::IsScalar) {
        out = f * J;
      } else {
        out = J.transpose() * f;
      }
    }

    /**
     * Compute \f$ \mathbf{J}^\top\mathbf{J} + \sum_i f_i \nabla^2 f_i \f$, the
     * Hessian of the least-squares merit function \f$
     * \frac{1}{2}\|\mathbf{f}\|^2 \f$.
     */
    template <typename FunctionValue,
              typename FunctionJacobian,
              typename FunctionHessian,
              typename Out>
    static void least_squares_hessian(const FunctionValue &f,
                                      const FunctionJacobian &J,
                                      const FunctionHessian &H,
                                      Out &out) {
      if constexpr (TypeTrait<FunctionValue>::IsScalar &&
                    TypeTrait<FunctionJacobian>::IsScalar) {
        out = J * J + f * H;
      } else if constexpr (TypeTrait<FunctionValue>::IsScalar) {
        out = J * J.transpose() + f * H;
      } else {
        out = J.transpose() * J;
        for (Integer i{0}; i < static_cast<Integer>(H.size()); ++i) {
          Scalar f_i;
          if constexpr (TypeTrait<FunctionValue>::IsSparse) {
            f_i = f.coeff(i);
          } else {
            f_i = f(i);
          }
          if (f_i != 0.0) {
            out += f_i * H[i];
          }
        }
      }
    }

    /**
     * Solve the root-finding/optimization problem given a function class.
     * \tparam Kind Kind of problem to be solved.
     * \tparam FunctionInput Function input type.
     * \tparam FunctionOutput Function output type.
     * \tparam DerivedFunction Derived function class.
     * \param[in] function Function class.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return The convergence boolean flag.
     */
    template <Problem Kind,
              typename FunctionInput,
              typename FunctionOutput,
              typename DerivedFunction>
    bool solve(
        const FunctionBase<FunctionInput, FunctionOutput, DerivedFunction>
            &function,
        const Input &x_ini,
        Input &x_sol) {
#define CMD "Optimist::Solver::solve(...): "

      using FunctionType =
          FunctionBase<FunctionInput, FunctionOutput, DerivedFunction>;
      using FunctionFirstDerivative  = typename FunctionType::FirstDerivative;
      using FunctionSecondDerivative = typename FunctionType::SecondDerivative;

      // Lambda generators for the function and its derivatives as seen by the
      // solver implementation
      auto function_lambda = [&function](const Input &x, Output &out) -> bool {
        if constexpr (Kind == Problem::ROOT || Kind == Problem::COST) {
          return function.evaluate(x, out);
        } else if constexpr (Kind == Problem::LEAST_SQUARES) {
          FunctionOutput f;
          if (!function.evaluate(x, f)) {
            return false;
          }
          if constexpr (TypeTrait<FunctionOutput>::IsScalar) {
            out = 0.5 * f * f;
          } else {
            out = 0.5 * f.squaredNorm();
          }
          return std::isfinite(out);
        } else if constexpr (Kind == Problem::STATIONARY_COST) {
          return function.first_derivative(x, out);
        } else {
          FunctionOutput f;
          FunctionFirstDerivative J;
          if (!function.evaluate(x, f) || !function.first_derivative(x, J)) {
            return false;
          }
          least_squares_gradient(f, J, out);
          return true;
        }
      };

      auto first_derivative_lambda = [&function](const Input &x,
                                                 FirstDerivative &out) -> bool {
        if constexpr (Kind == Problem::ROOT || Kind == Problem::COST) {
          return function.first_derivative(x, out);
        } else if constexpr (Kind == Problem::LEAST_SQUARES) {
          FunctionOutput f;
          FunctionFirstDerivative J;
          if (!function.evaluate(x, f) || !function.first_derivative(x, J)) {
            return false;
          }
          least_squares_gradient(f, J, out);
          return true;
        } else if constexpr (Kind == Problem::STATIONARY_COST) {
          return function.second_derivative(x, out);
        } else {
          FunctionOutput f;
          FunctionFirstDerivative J;
          FunctionSecondDerivative H;
          if (!function.evaluate(x, f) || !function.first_derivative(x, J) ||
              !function.second_derivative(x, H)) {
            return false;
          }
          least_squares_hessian(f, J, H, out);
          return true;
        }
      };

      auto second_derivative_lambda = [&function](
                                          const Input &x,
                                          SecondDerivative &out) -> bool {
        if constexpr (Kind == Problem::ROOT || Kind == Problem::COST) {
          return function.second_derivative(x, out);
        } else if constexpr (Kind == Problem::LEAST_SQUARES) {
          FunctionOutput f;
          FunctionFirstDerivative J;
          FunctionSecondDerivative H;
          if (!function.evaluate(x, f) || !function.first_derivative(x, J) ||
              !function.second_derivative(x, H)) {
            return false;
          }
          least_squares_hessian(f, J, H, out);
          return true;
        } else {
          static_cast<void>(function);
          static_cast<void>(x);
          static_cast<void>(out);
          OPTIMIST_ERROR(CMD "second derivative not available for "
                             "stationary-point problems.");
          return false;
        }
      };

      // Select solver method based on derivative requirements
      if constexpr (DerivedSolver::RequiresFunction &&
                    !DerivedSolver::RequiresFirstDerivative &&
                    !DerivedSolver::RequiresSecondDerivative) {
        return static_cast<DerivedSolver *>(this)->solve_impl(function_lambda,
                                                              x_ini,
                                                              x_sol);
      } else if constexpr (DerivedSolver::RequiresFunction &&
                           DerivedSolver::RequiresFirstDerivative &&
                           !DerivedSolver::RequiresSecondDerivative) {
        return static_cast<DerivedSolver *>(this)->solve_impl(
            function_lambda,
            first_derivative_lambda,
            x_ini,
            x_sol);
      } else if constexpr (DerivedSolver::RequiresFunction &&
                           DerivedSolver::RequiresFirstDerivative &&
                           DerivedSolver::RequiresSecondDerivative) {
        return static_cast<DerivedSolver *>(this)->solve_impl(
            function_lambda,
            first_derivative_lambda,
            second_derivative_lambda,
            x_ini,
            x_sol);
      } else {
        static_assert(DerivedSolver::RequiresFunction,
                      CMD "no matching function signature for the solver.");
        return false;
      }

#undef CMD
    }

    /**
     * Reset internal counters and flags.
     */
    void reset_counters() {
      this->m_function_evaluations          = 0;
      this->m_first_derivative_evaluations  = 0;
      this->m_second_derivative_evaluations = 0;
      this->m_iterations                    = 0;
      this->m_relaxations                   = 0;
      this->m_converged                     = false;
    }

    /**
     * Evaluate the function.
     * \tparam FunctionLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] x Input point.
     * \param[out] out Function value.
     * \return The boolean flag for successful evaluation.
     */
    template <typename FunctionLambda>
    bool evaluate_function(FunctionLambda &&function,
                           const Input &x,
                           Output &out) {
      OPTIMIST_ASSERT(
          this->m_function_evaluations < this->m_max_function_evaluations,
          "Optimist::" << this->name()
                       << "::evaluate_function(...): maximum allowed "
                          "function evaluations reached.");
      ++this->m_function_evaluations;
      return function(x, out);
    }

    /**
     * Evaluate the first derivative.
     * \tparam FirstDerivativeLambda First derivative lambda type.
     * \param[in] function First derivative lambda.
     * \param[in] x Input point.
     * \param[out] out First derivative value.
     * \return The boolean flag for successful evaluation.
     */
    template <typename FirstDerivativeLambda>
    bool evaluate_first_derivative(FirstDerivativeLambda &&function,
                                   const Input &x,
                                   FirstDerivative &out) {
      OPTIMIST_ASSERT(
          this->m_first_derivative_evaluations <
              this->m_max_first_derivative_evaluations,
          "Optimist::" << this->name()
                       << "::evaluate_first_derivative(...): maximum allowed "
                          "first derivative evaluations reached.");
      ++this->m_first_derivative_evaluations;
      return function(x, out);
    }

    /**
     * Evaluate the second derivative.
     * \param[in] function Second derivative lambda.
     * \param[in] x Input point.
     * \param[out] out Second derivative value.
     * \return The boolean flag for successful evaluation.
     */
    template <typename SecondDerivativeLambda>
    bool evaluate_second_derivative(SecondDerivativeLambda &&function,
                                    const Input &x,
                                    SecondDerivative &out) {
      OPTIMIST_ASSERT(
          this->m_second_derivative_evaluations <
              this->m_max_second_derivative_evaluations,
          "Optimist::" << this->name()
                       << "::evaluate_second_derivative(...): maximum allowed "
                          "second derivative evaluations reached.");
      ++this->m_second_derivative_evaluations;
      return function(x, out);
    }

    /**
     * Damp the step using the affine invariant criterion.
     * \tparam FunctionLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] x_old Old point.
     * \param[in] function_old Old function value.
     * \param[in] step_old Old step.
     * \param[in] x_new New point.
     * \param[in] function_new New function value.
     * \param[out] step_new New step.
     * \return The damping boolean flag, true if the damping is successful,
     * false otherwise.
     */
    template <typename FunctionLambda>
    bool damp(FunctionLambda &&function,
              const Input &x_old,
              const Input &function_old,
              const Input &step_old,
              Input &x_new,
              Input &function_new,
              Input &step_new) {
#define CMD "Optimist::Solver::damp(...): "

      Scalar step_norm_old, step_norm_new, residuals_old, residuals_new,
          tau{1.0};
      for (this->m_relaxations = 0;
           this->m_relaxations < this->m_max_relaxations;
           ++this->m_relaxations) {
        // Update point
        step_new = tau * step_old;
        x_new    = x_old + step_new;
        bool success{
          this->evaluate_function(function,
                                  x_new,
                                  function_new)};
        OPTIMIST_ASSERT(success,
                        CMD "function evaluation failed during damping.");

        // Compute step norms
        if constexpr (InputTrait::IsEigen) {
          step_norm_old = step_old.norm();
          step_norm_new = step_new.norm();
        } else {
          step_norm_old = std::abs(step_old);
          step_norm_new = std::abs(step_new);
        }

        // Compute residual norms
        if constexpr (OutputTrait::IsEigen) {
          residuals_old = function_old.norm();
          residuals_new = function_new.norm();
        } else {
          residuals_old = std::abs(function_old);
          residuals_new = std::abs(function_new);
        }

        // Check relaxation
        if (residuals_new < residuals_old ||
            step_norm_new < (Scalar(1.0) - tau / Scalar(2.0)) * step_norm_old) {
          return true;
        } else {
          tau *= this->m_alpha;
        }
      }
      return false;

#undef CMD
    }

    /**
     * Print the table header solver information.
     * \note This has to be properly placed in the derived classes.
     */
    void header() {
      static constexpr std::string_view c_tl{table_top_left_corner()};
      static constexpr std::string_view c_tr{table_top_right_corner()};
      static std::string h_7{table_horizontal_line<7>()};
      static std::string h_14{table_horizontal_line<14>()};
      static std::string h_23{table_horizontal_line<23>()};
      static std::string h_78{table_horizontal_line<78>()};
      static constexpr std::string vert_line{table_vertical_line()};
      static constexpr std::string v_ll{vert_line + " "};
      static constexpr std::string v_rr{" " + vert_line};
      static constexpr std::string v_lc{" " + vert_line + " "};
      static constexpr std::string_view j_tt{table_top_junction()};
      static constexpr std::string_view j_cc{table_center_cross()};
      static constexpr std::string_view j_ll{table_left_junction()};
      static constexpr std::string_view j_rr{table_right_junction()};

      *this->m_ostream << c_tl << h_78 << c_tr << std::endl
                       << v_ll << "Solver:" << std::setw(69) << this->name()
                       << v_rr << std::endl
                       << v_ll << "Task:  " << std::setw(69) << this->task()
                       << v_rr << std::endl
                       << j_ll << h_7 << j_tt << h_7 << j_tt << h_7 << j_tt
                       << h_7 << j_tt << h_7 << j_tt << h_14 << j_tt << h_23
                       << j_rr << std::endl
                       << v_ll << "#Iter" << v_lc << "   #f" << v_lc << "  #Df"
                       << v_lc << " #DDf" << v_lc << " #Rlx" << v_lc
                       << "   ║f(x)║₂  " << v_lc << "Additional notes"
                       << std::setw(9) << v_rr << std::endl
                       << j_ll << h_7 << j_cc << h_7 << j_cc << h_7 << j_cc
                       << h_7 << j_cc << h_7 << j_cc << h_14 << j_cc << h_23
                       << j_rr << std::endl;
    }

    /**
     * Print the table bottom solver information.
     * \note This has to be properly placed in the derived classes.
     */
    void bottom() {
      static constexpr std::string_view c_bl{table_bottom_left_corner()};
      static constexpr std::string_view c_br{table_bottom_right_corner()};
      static std::string h_7{table_horizontal_line<7>()};
      static std::string h_14{table_horizontal_line<14>()};
      static std::string h_23{table_horizontal_line<23>()};
      static std::string h_78{table_horizontal_line<78>()};
      static constexpr std::string vert_line{table_vertical_line()};
      static constexpr std::string v_ll{vert_line + " "};
      static constexpr std::string v_rr{" " + vert_line};
      static constexpr std::string_view j_ll{table_left_junction()};
      static constexpr std::string_view j_rr{table_right_junction()};
      static constexpr std::string_view j_bb{table_bottom_junction()};

      *this->m_ostream << j_ll << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb
                       << h_7 << j_bb << h_7 << j_bb << h_14 << j_bb << h_23
                       << j_rr << std::endl
                       << v_ll << std::setw(40)
                       << (this->m_converged ? "CONVERGED" : "NOT CONVERGED")
                       << std::setw(40) << v_rr << std::endl
                       << c_bl << h_78 << c_br << std::endl;
    }

    /**
     * Print the solver information during the algorithm iterations.
     * \note This has to be properly placed in the derived classes.
     */
    void info(Scalar residuals, const std::string &notes = "-") {
      static constexpr std::string vert_line{table_vertical_line()};
      static constexpr std::string v_rr{" " + vert_line};
      static constexpr std::string v_ll{vert_line + " "};
      static constexpr std::string v_lc{" " + vert_line + " "};

      *this->m_ostream << v_ll << std::setw(5) << this->m_iterations << v_lc
                       << std::setw(5) << this->m_function_evaluations << v_lc
                       << std::setw(5) << this->m_first_derivative_evaluations
                       << v_lc << std::setw(5)
                       << this->m_second_derivative_evaluations << v_lc
                       << std::setw(5) << this->m_relaxations << v_lc
                       << std::setw(12) << std::scientific
                       << std::setprecision(6) << residuals << v_lc << notes
                       << std::setw(25 - notes.length()) << v_rr << std::endl;
    }

  };  // class Solver

}  // namespace Optimist

#endif  // OPTIMIST_SOLVER_HH
