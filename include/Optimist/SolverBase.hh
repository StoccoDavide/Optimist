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

#ifndef OPTIMIST_SOLVER_HH
#define OPTIMIST_SOLVER_HH

#include "Optimist.hh"
#include "Optimist/Function.hh"

namespace Optimist
{

  /*\
   |   ____        _                ____
   |  / ___|  ___ | |_   _____ _ __| __ )  __ _ ___  ___
   |  \___ \ / _ \| \ \ / / _ \ '__|  _ \ / _` / __|/ _ \
   |   ___) | (_) | |\ V /  __/ |  | |_) | (_| \__ \  __/
   |  |____/ \___/|_| \_/ \___|_|  |____/ \__,_|___/\___|
   |
  \*/

  /**
  * \brief Class container for the generic root-finding/optimization problem solver.
  *
  * \includedoc docs/markdown/SolverBase.md
  *
  * \tparam Input Solver input type.
  * \tparam Output Solver output type.
  * \tparam DerivedSolver Derived solver class.
  */
  template <typename Input, typename Output, typename DerivedSolver>
  class SolverBase
  {
  public:
    // Input and output types
    using InputTrait  = TypeTrait<Input>;
    using OutputTrait = TypeTrait<Output>;
    using Scalar = typename InputTrait::Scalar;

    // Input and output must be non-zero dimensional
    static_assert(InputTrait::Dimension != 0,
      "Input must be non-zero dimensional");
    static_assert(OutputTrait::Dimension != 0,
      "Output must be non-zero dimensional");

    // Input and output must have the same scalar type
    static_assert(std::is_same<typename InputTrait::Scalar, typename OutputTrait::Scalar>::value,
      "Input and output scalar types must be the same.");

    // If both input and output are eigen types they be both fixed-size, dynamic-size, or sparse
    static_assert(!(InputTrait::IsEigen && OutputTrait::IsEigen) ||
      (InputTrait::IsFixed && OutputTrait::IsFixed) ||
      (InputTrait::IsDynamic && OutputTrait::IsDynamic) ||
      (InputTrait::IsSparse && OutputTrait::IsSparse),
      "Input and output Eigen types must be both fixed-size, dynamic-size, or sparse.");

    // Derivative types
    using FirstDerivative = std::conditional_t<InputTrait::IsEigen || OutputTrait::IsEigen,
      std::conditional_t<InputTrait::IsSparse || OutputTrait::IsSparse,
        Eigen::SparseMatrix<Scalar>,
        Eigen::Matrix<Scalar, OutputTrait::Dimension, InputTrait::Dimension>>,
      Scalar>;
    using SecondDerivative = std::conditional_t<InputTrait::IsEigen || OutputTrait::IsEigen,
      std::conditional_t<(InputTrait::Dimension == 1) || (OutputTrait::Dimension == 1),
        std::conditional_t<InputTrait::IsSparse || OutputTrait::IsSparse,
          Eigen::SparseMatrix<Scalar>,
          Eigen::Matrix<Scalar, OutputTrait::Dimension, InputTrait::Dimension>>,
        std::conditional_t<InputTrait::IsSparse || OutputTrait::IsSparse,
          std::vector<Eigen::SparseMatrix<Scalar>>,
          std::vector<Eigen::Matrix<Scalar, OutputTrait::Dimension, InputTrait::Dimension>>>>,
      Scalar>;

    OPTIMIST_BASIC_CONSTANTS(Scalar)

  protected:
    // Bounds (may not be used)
    Input m_lower_bound; /**< Lower bound. */
    Input m_upper_bound; /**< Upper bound. */

    // Evaluations
    Integer m_function_evaluations{0};          /**< Function evaluations. */
    Integer m_first_derivative_evaluations{0};  /**< First derivative evaluations. */
    Integer m_second_derivative_evaluations{0}; /**< Second derivative evaluations. */

    // Maximum allowed evaluations
    Integer m_max_function_evaluations{1000};          /**< Maximum allowed function evaluations. */
    Integer m_max_first_derivative_evaluations{1000};  /**< Maximum allowed first derivative evaluations. */
    Integer m_max_second_derivative_evaluations{1000}; /**< Maximum allowed second derivative evaluations. */

    // Iterations and relaxations
    Integer m_iterations{0};       /**< Algorithm iterations. */
    Integer m_max_iterations{100}; /**< Maximum allowed algorithm iterations. */
    Scalar  m_alpha{0.8};          /**< Relaxation factor \f$ \alpha \f$. */
    Integer m_relaxations{0};      /**< Algorithm relaxations. */
    Integer m_max_relaxations{10}; /**< Maximum allowed algorithm relaxations. */

    // Settings
    Scalar         m_tolerance{EPSILON_LOW}; /**< Solver tolerance \f$ \epsilon \f$ for convergence. */
    bool           m_verbose{false};         /**< Verbose mode boolean flag. */
    bool           m_damped{true};           /**< Damped mode boolean flag. */
    std::ostream * m_ostream{&std::cout};    /**< Output stream for verbose mode. */

    // Convergence output flag and trace
    std::string m_task{"Undefined"}; /**< Task name. */
    bool        m_converged{false};  /**< Convergence boolean flag. */

  public:
    /**
     * Class constructor for the nonlinear solver.
     */
    SolverBase() {}

    /**
     * Class constructor for the nonlinear solver.
     * \tparam FunctionLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     */
    template <typename FunctionLambda>
    SolverBase(FunctionLambda && function, Input const & x_ini, Input & x_sol) : SolverBase() {
      static_cast<const DerivedSolver *>(this)->solve_impl(
        std::forward<FunctionLambda>(function),
        x_ini, x_sol);
    }

    /**
     * Class constructor for the nonlinear solver.
     * \tparam FunctionLambda Function lambda type.
     * \tparam FirstDerivativeLambda First derivative lambda type.
     * \param[in] function Function lambda.
     * \param[in] first_derivative First derivative lambda.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     */
    template <typename FunctionLambda, typename FirstDerivativeLambda>
    SolverBase(FunctionLambda && function, FirstDerivativeLambda && first_derivative, Input const & x_ini,
      Input & x_sol) : SolverBase()
    {
      static_cast<const DerivedSolver *>(this)->solve_impl(
        std::forward<FunctionLambda>(function),
        std::forward<FirstDerivativeLambda>(first_derivative),
        x_ini, x_sol);
    }

    /**
     * Class constructor for the nonlinear solver.
     * \tparam FunctionLambda Function lambda type.
     * \tparam FirstDerivativeLambda First derivative lambda type.
     * \tparam SecondDerivativeLambda Second derivative lambda type.
     * \param[in] function Function lambda.
     * \param[in] first_derivative First derivative lambda.
     * \param[in] second_derivative The second derivative lambda.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     */
    template <typename FunctionLambda, typename FirstDerivativeLambda, typename SecondDerivativeLambda>
    SolverBase(FunctionLambda && function, FirstDerivativeLambda && first_derivative, SecondDerivativeLambda
      && second_derivative, Input const & x_ini, Input & x_sol) : SolverBase()
    {
      static_cast<const DerivedSolver *>(this)->solve_impl(
        std::forward<FunctionLambda>(function),
        std::forward<FirstDerivativeLambda>(first_derivative),
        std::forward<SecondDerivativeLambda>(second_derivative),
        x_ini, x_sol);
    }

    /**
     * Reset lower and upper bounds to default values.
     * \param[in] n Input dimension for dynamic-size types.
     */
    void reset_bounds(Integer n = InputTrait::IsDynamic ? 0 : InputTrait::Dimension)
    {
      #define CMD "Optimist::Solver::reset_bounds(...): "

      if constexpr (InputTrait::IsScalar) {
        this->m_lower_bound = -INFTY;
        this->m_upper_bound = +INFTY;
      } else if constexpr (InputTrait::IsFixed) {
        this->m_lower_bound.setConstant(-INFTY);
        this->m_upper_bound.setConstant(+INFTY);
      } else if constexpr (InputTrait::IsDynamic) {
        this->m_lower_bound.resize(n); this->m_lower_bound.setConstant(-INFTY);
        this->m_upper_bound.resize(n); this->m_upper_bound.setConstant(+INFTY);
      } else if constexpr (InputTrait::IsSparse) {
        this->m_lower_bound.resize(n); this->m_lower_bound.reserve(n);
        this->m_upper_bound.resize(n); this->m_upper_bound.reserve(n);
        std::vector<Eigen::Triplet<Scalar>> triplets; triplets.reserve(n);
        for (Integer i{0}; i < n; ++i) {triplets.emplace_back(i, 0, -INFTY);}
        this->m_lower_bound.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear(); triplets.reserve(n);
        for (Integer i{0}; i < n; ++i) {triplets.emplace_back(i, 0, +INFTY);}
        this->m_upper_bound.setFromTriplets(triplets.begin(), triplets.end());
      } else {
        OPTIMIST_ERROR(CMD "unsupported input type for bounds reset.");
      }

      #undef CMD
    }

    /**
     * Get the lower bound.
     * \return The lower bound.
     */
    Input const & lower_bound() const {return this->m_lower_bound;}

    /**
     * Set the lower bound.
     * \param[in] t_lower_bound The lower bound.
     */
    void lower_bound(Input const & t_lower_bound) {
      #define CMD "Optimist::Solver::bounds(...): "

      if constexpr (InputTrait::IsEigen) {
        OPTIMIST_ASSERT((this->m_upper_bound - t_lower_bound).minCoeff() <= 0.0,
          CMD "invalid or degenarate bounds detected.");
      } else {
        OPTIMIST_ASSERT(this->m_upper_bound > t_lower_bound,
          CMD "invalid or degenarate bounds detected.");
      }
      this->m_lower_bound = t_lower_bound;

      #undef CMD
    }

    /**
     * Get the upper bound.
     * \return The upper bound.
     */
    Input const & upper_bound() const {return this->m_upper_bound;}

    /**
     * Set the upper bound.
     * \param[in] t_upper_bound The upper bound.
     */
    void upper_bound(Input const & t_upper_bound) {
      #define CMD "Optimist::Solver::bounds(...): "

      if constexpr (InputTrait::IsEigen) {
        OPTIMIST_ASSERT((t_upper_bound - this->m_lower_bound).minCoeff() <= 0.0,
          CMD "invalid or degenarate bounds detected.");
      } else {
        OPTIMIST_ASSERT(t_upper_bound > this->m_lower_bound,
          CMD "invalid or degenarate bounds detected.");
      }
      this->m_upper_bound = t_upper_bound;

      #undef CMD
    }

    /**
     * Set the bounds.
     * \param[in] t_lower_bound The lower bound.
     * \param[in] t_upper_bound The upper bound.
     */
    void bounds(Input const & t_lower_bound, Input const & t_upper_bound)
    {
      #define CMD "Optimist::Solver::bounds(...): "

      if constexpr (InputTrait::IsEigen) {
        OPTIMIST_ASSERT((t_upper_bound - t_lower_bound).minCoeff() <= 0.0,
          CMD "invalid or degenarate bounds detected.");
      } else {
        OPTIMIST_ASSERT(t_upper_bound > t_lower_bound,
          CMD "invalid or degenarate bounds detected.");
      }
      this->m_lower_bound = t_lower_bound;
      this->m_upper_bound = t_upper_bound;

      #undef CMD
    }

    /**
     * Get the input dimension of the function.
     * \return The input dimension of the function.
     */
    constexpr Integer input_dimension() const {return InputTrait::Dimension;}

    /**
     * Get the output dimension of the function.
     * \return The output dimension of the function.
     */
    constexpr Integer output_dimension() const {return OutputTrait::Dimension;}

    /**
     * Get the number of function evaluations.
     * \return The number of function evaluations.
     */
    Integer function_evaluations() const {return this->m_function_evaluations;}

    /**
     * Set the number of maximum allowed function evaluations.
     * \param[in] t_max_function_evaluations The number of maximum allowed function evaluations.
     */
    void max_function_evaluations(Integer t_max_function_evaluations)
    {
      OPTIMIST_ASSERT(t_max_function_evaluations > 0,
        "Optimist::Solver::max_function_evaluations(...): invalid input detected.");
      this->m_max_function_evaluations = t_max_function_evaluations;
    }

    /**
     * Get the number of maximum allowed function evaluations.
     * \return The number of maximum allowed function evaluations.
     */
    Integer max_function_evaluations() const {return this->m_max_function_evaluations;}

  protected:
    /**
     * Get the number of function derivative evaluations.
     * \return The number of function derivative evaluations.
     */
    Integer first_derivative_evaluations() const {return this->m_first_derivative_evaluations;}

    /**
     * Get the number of maximum allowed first derivative evaluations.
     * \return The number of maximum allowed first derivative evaluations.
     */
    Integer max_first_derivative_evaluations() const {return this->m_max_first_derivative_evaluations;}

    /**
     * Set the number of maximum allowed first derivative evaluations.
     * \param[in] first_derivative_evaluations The number of maximum allowed first derivative evaluations.
     */
    void max_first_derivative_evaluations(Integer first_derivative_evaluations)
    {
      OPTIMIST_ASSERT(first_derivative_evaluations > 0,
        "Optimist::Solver::max_first_derivative_evaluations(...): invalid input detected.");
      this->m_max_first_derivative_evaluations = first_derivative_evaluations;
    }

    /**
     * Get the number of second derivative evaluations.
     * \return The number of second derivative evaluations.
     */
    Integer second_derivative_evaluations() const {return this->m_second_derivative_evaluations;}

    /**
     * Get the number of maximum allowed second derivative evaluations.
     * \return The number of maximum allowed second derivative evaluations.
     */
    Integer max_second_derivative_evaluations() const {return this->m_max_second_derivative_evaluations;}

    /**
     * Set the number of maximum allowed second derivative evaluations.
     * \param[in] second_derivative_evaluations The number of maximum allowed second derivative evaluations.
     */
    void max_second_derivative_evaluations(Integer second_derivative_evaluations)
    {
      OPTIMIST_ASSERT(second_derivative_evaluations > 0,
        "Optimist::Solver::max_second_derivative_evaluations(...): invalid input detected.");
      this->m_max_second_derivative_evaluations = second_derivative_evaluations;
    }

  public:
    /**
     * Get the number of algorithm iterations.
     */
    Integer iterations() const {return this->m_iterations;}

    /**
     * Get the number of maximum allowed iterations.
     * \return The number of maximum allowed iterations.
     */
    Integer max_iterations() const {return this->max_iterations;}

    /**
     * Set the number of maximum allowed iterations.
     * \param[in] t_max_iterations The number of maximum allowed iterations.
     */
    void max_iterations(Integer t_max_iterations) {
      OPTIMIST_ASSERT(t_max_iterations > 0,
        "Optimist::Solver::max_iterations(...): invalid input detected.");
      this->m_max_iterations = t_max_iterations;
    }

    /**
     * Get relaxation factor \f$ \alpha \f$.
     * \return The relaxation factor \f$ \alpha \f$.
     */
    Scalar alpha() const {return this->m_alpha;}

    /**
     * Set relaxation factor \f$ \alpha \f$.
     * \param[in] t_alpha The relaxation factor \f$ \alpha \f$.
     */
    void alpha(Scalar t_alpha)
    {
      OPTIMIST_ASSERT(
        !std::isnan(t_alpha) && std::isfinite(t_alpha) && 0.0 <= t_alpha && t_alpha <= 1.0,
        "Optimist::Solver::alpha(...): invalid input detected.");
      this->m_alpha = t_alpha;
    }

    /**
     * Get the number of algorithm relaxations.
     * \return The number of algorithm relaxations.
     */
    Integer relaxations() const {return this->m_relaxations;}

    /**
     * Get the number of maximum allowed relaxations.
     * \return The number of maximum allowed relaxations.
     */
    Integer max_relaxations() const {return this->max_relaxations;}

    /**
     * Set the number of maximum allowed relaxations.
     * \param[in] t_max_relaxations The number of maximum allowed relaxations.
     */
    void max_relaxations(Integer t_max_relaxations)
    {
      OPTIMIST_ASSERT(t_max_relaxations > 0,
        "Optimist::Solver::max_relaxations(...): invalid input detected.");
      this->m_max_relaxations = t_max_relaxations;
    }

    /**
     * Get the tolerance \f$ \epsilon \f$.
     * \return The tolerance \f$ \epsilon \f$.
     */
    Scalar tolerance() const {return this->m_tolerance;}

    /**
     * Set the tolerance \f$ \epsilon \f$ for which the nonlinear solver stops, i.e., \f$ \left\|
     * \mathbf{F}(\mathbf{x}) \right\|_{2} < \epsilon \f$.
     * \param[in] t_tolerance The tolerance \f$ \epsilon \f$.
     */
    void tolerance(Scalar t_tolerance) {
      OPTIMIST_ASSERT(!std::isnan(t_tolerance) && std::isfinite(t_tolerance) && t_tolerance > 0,
        "Optimist::Solver::tolerance(...): invalid input detected.");
      this->m_tolerance = t_tolerance;
    }

    /**
     * Set the verbose mode boolean flag.
     * \param[in] t_verbose The verbose mode boolean flag.
     */
    void verbose_mode(bool t_verbose) {this->m_verbose = t_verbose;}

    /**
     * Get the verbose mode boolean flag.
     * \return The verbose mode boolean flag.
     */
    bool verbose_mode() const {return this->m_verbose;}

    /**
     * Enable solver's verbose mode.
     */
    void enable_verbose_mode() {this->m_verbose = true;}

    /**
     * Disable solver's verbose mode.
     */
    void disable_verbose_mode() {this->m_verbose = false;}

    /**
     * Set the damped mode boolean flag.
     * \param[in] t_damped The damped mode boolean flag.
     */
    void damped_mode(bool t_damped) {this->m_damped = t_damped;}

    /**
     * Get the damped mode boolean flag.
     * \return The damped mode boolean flag.
     */
    bool damped_mode() const {return this->m_damped;}

    /**
     * Enable solver's damped mode.
     */
    void enable_damped_mode() {this->m_damped = true;}

    /**
     * Disable solver's damped mode.
     */
    void disable_damped_mode() {this->m_damped = false;}

    /**
     * Get the task name.
     * \return The task name.
     */
    std::string task() const {return this->m_task;}

    /**
     * Set the task name.
     * \param[in] t_task The task name.
     */
    void task(std::string t_task) {this->m_task = t_task;}

    /**
     * Get the convergence boolean flag.
     * \return The convergence boolean flag.
     */
    bool converged() const {return this->m_converged;}

    /**
     * Get the output stream for verbose mode.
     * \return The output stream for verbose mode.
     */
    std::ostream & ostream() const {return *this->m_ostream;}

    /**
     * Set the output stream for verbose mode.
     * \param[in] t_ostream The output stream for verbose mode.
     */
    void ostream(std::ostream & t_ostream) {this->m_ostream = &t_ostream;}

    /**
     * Solve the root-finding/optimization problem given the function, and without derivatives.
     * \tparam FunctionLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return True if the problem is solved, false otherwise.
     */
    template <typename FunctionLambda>
    bool solve(FunctionLambda && function, Input const & x_ini, Input & x_sol)
    {
      #define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::RequiresFunction,
        CMD "the solver requires a function.");
      return static_cast<DerivedSolver *>(this)->solve_impl(
        std::forward<FunctionLambda>(function),
        nullptr, nullptr, x_ini, x_sol
      );

      #undef CMD
    }

    /**
     * Solve the root-finding/optimization problem given the function, and its first derivative.
     * \tparam FunctionLambda First derivative lambda type.
     * \tparam FirstDerivativeLambda Function lambda type.
     * \param[in] function Function lambda.
     * \param[in] first_derivative First derivative lambda
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \return True if the problem is solved, false otherwise.
     */
    template <typename FunctionLambda, typename FirstDerivativeLambda>
    bool solve(FunctionLambda && function, FirstDerivativeLambda && first_derivative, Input const & x_ini,
      Input & x_sol)
    {
      #define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::RequiresFunction,
        CMD "the solver requires a function.");
      static_assert(DerivedSolver::RequiresFirstDerivative,
        CMD "the solver requires the first derivative.");
      return static_cast<DerivedSolver *>(this)->solve_impl(
        std::forward<FunctionLambda>(function),
        std::forward<FirstDerivativeLambda>(first_derivative),
        nullptr, x_ini, x_sol
      );

      #undef CMD
    }

    /**
     * Solve the root-finding/optimization problem given the function, and its first and second derivatives.
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
    template <typename FunctionLambda, typename FirstDerivativeLambda, typename SecondDerivativeLambda>
    bool solve(FunctionLambda && function, FirstDerivativeLambda && first_derivative, SecondDerivativeLambda
      && second_derivative, Input const & x_ini, Input & x_sol)
      {
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
        x_ini, x_sol
      );

      #undef CMD
    }

    /**
     * Solve the root-finding problem given a function class.
     * \tparam FunctionInput The function input dimension.
     * \tparam FunctionOutput The function output dimension.
     * \tparam DerivedFunction Derived function class.
     * \param[in] function Function class.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     */
    template <typename FunctionInput, typename FunctionOutput, typename DerivedFunction>
    requires (InputTrait::Dimension == TypeTrait<FunctionInput>::Dimension) &&
      (OutputTrait::Dimension == TypeTrait<FunctionOutput>::Dimension || OutputTrait::Dimension == 1) &&
      (!(InputTrait::Dimension == 1 && DerivedSolver::IsOptimizer))
    bool rootfind(FunctionBase<FunctionInput, FunctionOutput, DerivedFunction> const & function,
      Input const & x_ini, Input & x_sol)
    {
      using FunctionOutputTrait = TypeTrait<FunctionOutput>;
      return this->solve(function, x_ini, x_sol,
        (OutputTrait::Dimension != FunctionOutputTrait::Dimension) || InputTrait::IsEigen);
    }

    /**
     * Solve the optimization problem given a function class.
     * \tparam FunctionInput The function input dimension.
     * \tparam FunctionOutput The function output dimension.
     * \tparam DerivedFunction Derived function class.
     * \param[in] function Function class.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     */
    template <typename FunctionInput, typename FunctionOutput, typename DerivedFunction>
    bool optimize(FunctionBase<FunctionInput, FunctionOutput, DerivedFunction> const & function,
      Input const & x_ini, Input & x_sol)
    {
      #define CMD "Optimist::Solver::optimize(...): "

      using FunctionInputTrait = TypeTrait<FunctionInput>;

      static_assert(InputTrait::Dimension == FunctionInputTrait::Dimension,
        CMD "solver input dimension must be equal to the function input dimension.");
      static_assert(OutputTrait::Dimension == 1,
        CMD "solver output dimension must be equal to the function output dimension or 1.");
      static_assert(!(InputTrait::Dimension == 1 && DerivedSolver::IsRootFinder),
        CMD "one-dimensional root-finders do not support optimization problems.");
      return this->solve(function, x_ini, x_sol, true);

      #undef CMD
    }

    /**
     * Get the solver name.
     * \return The solver name.
     */
    constexpr std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();};

  protected:
    /**
     * Solve the root-finding/optimization problem given a function class.
     * \tparam Input Solver input type.
     * \tparam Output Solver output type.
     * \tparam DerivedFunction Derived function class.
     * \param[in] function Function class.
     * \param[in] x_ini Initialization point.
     * \param[out] x_sol Solution point.
     * \param[in] is_optimization Boolean flag for optimization.
     */
    template <typename FunctionInput, typename FunctionOutput, typename DerivedFunction>
    requires (InputTrait::Dimension == TypeTrait<FunctionInput>::Dimension) &&
      (OutputTrait::Dimension == TypeTrait<FunctionOutput>::Dimension || OutputTrait::Dimension == 1)
    bool solve(FunctionBase<FunctionInput, FunctionOutput, DerivedFunction>
      const & function, Input const & x_ini, Input & x_sol, bool is_optimization)
    {
      #define CMD "Optimist::Solver::solve(...): "

      using FunctionType = FunctionBase<FunctionInput, FunctionOutput, DerivedFunction>;
      using FunctionInputTrait  = TypeTrait<FunctionInput>;
      using FunctionOutputTrait = TypeTrait<FunctionOutput>;

      // Lambda generators for function and derivatives
      auto function_lambda = [&function, is_optimization] (Input const & x, Output & out) -> bool
      {
        bool success{false};
        FunctionOutput f; success = function.evaluate(x, f);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed during function computation.");

        if (is_optimization) {
          if constexpr (FunctionOutputTrait::Dimension == 1) {out = 0.5*f*f;}
          else if constexpr (OutputTrait::Dimension != FunctionOutputTrait::Dimension) {out = 0.5*f.squaredNorm();}
          else {OPTIMIST_ERROR(CMD "optimization problem with inconsistent output in function.");}
        } else {
          if constexpr (OutputTrait::Dimension == FunctionOutputTrait::Dimension) {out = f;}
          else {OPTIMIST_ERROR(CMD "root-finding problem with inconsistent output in function.");}
        }
        return success;
      };

      auto first_derivative_lambda = [&function, is_optimization, this] (Input const & x,
        FirstDerivative & out) -> bool
      {
        bool success{false};
        typename FunctionType::FirstDerivative J; success = function.first_derivative(x, J);
        OPTIMIST_ASSERT(success,
          CMD "first derivative evaluation failed during first derivative computation.");

        if (is_optimization) {
          bool success{false};
          FunctionOutput f; success = function.evaluate(x, f);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed during first derivative computation.");
          this->m_function_evaluations++;
          if constexpr (FunctionInputTrait::Dimension == 1 && FunctionOutputTrait::Dimension == 1) {out = J*f;}
          else if constexpr (OutputTrait::Dimension != FunctionOutputTrait::Dimension) {out = J.transpose()*f;}
          else {OPTIMIST_ERROR(CMD "optimization problem inconsistent output in first derivative.");}
        } else {
          if constexpr (OutputTrait::Dimension == FunctionOutputTrait::Dimension) {out = J;}
          else {OPTIMIST_ERROR(CMD "root-finding problem with inconsistent output in first derivative.");}
        }
        return success;
      };

      auto second_derivative_lambda = [&function, is_optimization, this] (Input const & x,
        SecondDerivative & out) -> bool
      {
        bool success{false};
        typename FunctionType::SecondDerivative H; success = function.second_derivative(x, H);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed during second derivative computation.");

        if (is_optimization) {
          FunctionOutput f; success = function.evaluate(x, f);
          this->m_function_evaluations++;
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed during second derivative computation.");
          typename FunctionType::FirstDerivative J; success = function.first_derivative(x, J);
          this->m_first_derivative_evaluations++;
          OPTIMIST_ASSERT(success,
            CMD "first derivative evaluation failed during second derivative computation.");
          if constexpr (FunctionInputTrait::Dimension == 1 && FunctionOutputTrait::Dimension == 1) {out = J*J + f*H;}
          else if constexpr (OutputTrait::Dimension != FunctionOutputTrait::Dimension) {
            out = J.transpose()*J;
            for (Integer i{0}; i < static_cast<Integer>(H.size()); ++i) {out += f(i)*H[i];}
          }
          else {OPTIMIST_ERROR(CMD "optimization problem with inconsistent output in second derivative.");}
        } else {
          if constexpr (OutputTrait::Dimension == FunctionOutputTrait::Dimension) {out = H;}
          else {OPTIMIST_ERROR(CMD "root-finding problem with inconsistent output in second derivative.");}
        }
        return success;
      };

      // Select solver method based on derivative requirements
      if constexpr (DerivedSolver::RequiresFunction && !DerivedSolver::RequiresFirstDerivative &&
        !DerivedSolver::RequiresSecondDerivative) {
        return static_cast<DerivedSolver *>(this)->solve_impl(function_lambda, x_ini, x_sol);
      } else if constexpr (DerivedSolver::RequiresFunction && DerivedSolver::RequiresFirstDerivative &&
        !DerivedSolver::RequiresSecondDerivative) {
        return static_cast<DerivedSolver *>(this)->solve_impl(function_lambda, first_derivative_lambda,
          x_ini, x_sol);
      } else if constexpr (DerivedSolver::RequiresFunction && DerivedSolver::RequiresFirstDerivative &&
        DerivedSolver::RequiresSecondDerivative) {
        return static_cast<DerivedSolver *>(this)->solve_impl(function_lambda, first_derivative_lambda,
          second_derivative_lambda, x_ini, x_sol);
      } else {
        OPTIMIST_ERROR(CMD "no matching function signature found for the solver.");
      }

      #undef CMD
    }

    /**
     * Reset internal counters and flags.
     */
    void reset_counters()
    {
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
    bool evaluate_function(FunctionLambda && function, Input const & x, Output & out)
    {
      OPTIMIST_ASSERT(this->m_function_evaluations < this->m_max_function_evaluations,
        "Optimist::" << this-> name() << "::evaluate_function(...): maximum allowed function evaluations reached.");
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
    bool evaluate_first_derivative(FirstDerivativeLambda && function, Input const & x, FirstDerivative & out)
    {
      OPTIMIST_ASSERT(this->m_first_derivative_evaluations < this->m_max_first_derivative_evaluations,
        "Optimist::" << this-> name() << "::evaluate_first_derivative(...): maximum allowed first derivative evaluations reached.");
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
    bool evaluate_second_derivative(SecondDerivativeLambda && function, Input const & x, SecondDerivative & out)
    {
      OPTIMIST_ASSERT(this->m_second_derivative_evaluations < this->m_max_second_derivative_evaluations,
        "Optimist::" << this-> name() << "::evaluate_second_derivative(...): maximum allowed second derivative evaluations reached.");
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
     * \return The damping boolean flag, true if the damping is successful, false otherwise.
     */
    template <typename FunctionLambda>
    bool damp(FunctionLambda && function, Input const & x_old, Input const & function_old,
      Input const & step_old, Input & x_new, Input & function_new, Input & step_new)
    {
      #define CMD "Optimist::Solver::damp(...): "

      Scalar step_norm_old, step_norm_new, residuals_old, residuals_new, tau{1.0};
      for (this->m_relaxations = 0; this->m_relaxations < this->m_max_relaxations; ++this->m_relaxations)
      {
        // Update point
        step_new = tau * step_old;
        x_new = x_old + step_new;
        bool success{this->evaluate_function(std::forward<FunctionLambda>(function), x_new, function_new)};
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
        if (residuals_new < residuals_old || step_norm_new < (Scalar(1.0)-tau/Scalar(2.0))*step_norm_old) {
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
    void header()
    {
      static constexpr std::string c_tl{table_top_left_corner()};
      static constexpr std::string c_tr{table_top_right_corner()};
      static constexpr std::string h_7{table_horizontal_line<7>()};
      static std::string h_14{table_horizontal_line<14>()};
      static std::string h_23{table_horizontal_line<23>()};
      static std::string h_78{table_horizontal_line<78>()};
      static constexpr std::string v_ll{table_vertical_line() + " "};
      static constexpr std::string v_rr{" " + table_vertical_line()};
      static constexpr std::string v_lc{" " + table_vertical_line() + " "};
      static constexpr std::string j_tt{table_top_junction()};
      static constexpr std::string j_cc{table_center_cross()};
      static constexpr std::string j_ll{table_left_junction()};
      static constexpr std::string j_rr{table_right_junction()};

      *this->m_ostream
        << c_tl << h_78 << c_tr << std::endl
        << v_ll << "Solver:" << std::setw(69) << this->name() << v_rr << std::endl
        << v_ll << "Task:  " << std::setw(69) << this->task() << v_rr << std::endl
        << j_ll << h_7 << j_tt << h_7 << j_tt << h_7 << j_tt << h_7 << j_tt << h_7 << j_tt << h_14 << j_tt << h_23 << j_rr << std::endl
        << v_ll << "#Iter" << v_lc << "   #f" << v_lc << "  #Df" << v_lc << " #DDf" << v_lc << " #Rlx" << v_lc
        << "   ║f(x)║₂  " << v_lc << "Additional notes" << std::setw(9) << v_rr << std::endl
        << j_ll << h_7 << j_cc << h_7 << j_cc << h_7 << j_cc << h_7 << j_cc << h_7 << j_cc << h_14 << j_cc << h_23 << j_rr << std::endl;
    }

    /**
     * Print the table bottom solver information.
     * \note This has to be properly placed in the derived classes.
     */
    void bottom()
    {
      static constexpr std::string c_bl{table_bottom_left_corner()};
      static constexpr std::string c_br{table_bottom_right_corner()};
      static std::string h_7{table_horizontal_line<7>()};
      static std::string h_14{table_horizontal_line<14>()};
      static std::string h_23{table_horizontal_line<23>()};
      static std::string h_78{table_horizontal_line<78>()};
      static constexpr std::string v_ll{table_vertical_line() + " "};
      static constexpr std::string v_rr{" " + table_vertical_line()};
      static constexpr std::string j_ll{table_left_junction()};
      static constexpr std::string j_rr{table_right_junction()};
      static constexpr std::string j_bb{table_bottom_junction()};

      *this->m_ostream
        << j_ll << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb << h_14 << j_bb << h_23 << j_rr << std::endl
        << v_ll << std::setw(40) << (this->m_converged ? "CONVERGED" : "NOT CONVERGED") << std::setw(40) << v_rr << std::endl
        << c_bl << h_78 << c_br << std::endl;
    }

    /**
     * Print the solver information during the algorithm iterations.
     * \note This has to be properly placed in the derived classes.
     */
    void info(Scalar residuals, std::string const & notes = "-")
    {
      static constexpr std::string v_rr{" " + table_vertical_line()};
      static constexpr std::string v_ll{table_vertical_line() + " "};
      static constexpr std::string v_lc{" " + table_vertical_line() + " "};

      *this->m_ostream << v_ll
        << std::setw(5) << this->m_iterations << v_lc
        << std::setw(5) << this->m_function_evaluations << v_lc
        << std::setw(5) << this->m_first_derivative_evaluations << v_lc
        << std::setw(5) << this->m_second_derivative_evaluations << v_lc
        << std::setw(5) << this->m_relaxations << v_lc
        << std::setw(12) << std::scientific << std::setprecision(6) << residuals << v_lc
        << notes << std::setw(25-notes.length()) << v_rr << std::endl;
    }

  }; // class Solver

} // namespace Optimist

#endif // OPTIMIST_SOLVER_HH
