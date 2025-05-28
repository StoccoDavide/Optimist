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
  * \tparam Real Scalar number type.
  * \tparam SolInDim The root-finding/optimization problem input dimension.
  * \tparam SolOutDim The root-finding/optimization problem output dimension.
  * \tparam DerivedSolver Derived solver class.
  * \tparam ForceEigen Force the use of Eigen types for input and output.
  */
  template <typename Real, Integer SolInDim, Integer SolOutDim, typename DerivedSolver, bool ForceEigen = false>
  class SolverBase
  {
    public:
    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(SolInDim > static_cast<Integer>(0) && SolOutDim > static_cast<Integer>(0),
      "Negative-dimensional optimization problem? Are you serious?");

    OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

  protected:
    // I/O types
    using InputType  = typename std::conditional_t<ForceEigen || (SolInDim > 1),  Eigen::Vector<Real, SolInDim>, Real>;  /**< Input type. */
    using OutputType = typename std::conditional_t<ForceEigen || (SolOutDim > 1), Eigen::Vector<Real, SolOutDim>, Real>; /**< Output type. */

    // Trace types
    using TraceType = typename std::vector<InputType>; /**< Input trace type. */

    // Derivative types
    using FirstDerivativeType  = std::conditional_t<ForceEigen || (SolInDim > 1) || (SolOutDim > 1), Eigen::Matrix<Real, SolOutDim, SolInDim>, Real>; /**< First derivative type. */
    using SecondDerivativeType = std::conditional_t<ForceEigen || (SolInDim > 1) || (SolOutDim > 1),
      std::conditional_t<SolInDim == 1 || SolOutDim == 1, Eigen::Matrix<Real, SolInDim, SolInDim>,
      std::vector<Eigen::Matrix<Real, SolInDim, SolInDim>>>, Real>;  /**< Second derivative type. */

    // Function types
    using FunctionWrapper         = typename std::function<void(const InputType &, OutputType &)>;           /**< Function type. */
    using FirstDerivativeWrapper  = typename std::function<void(const InputType &, FirstDerivativeType &)>;  /**< First derivative type. */
    using SecondDerivativeWrapper = typename std::function<void(const InputType &, SecondDerivativeType &)>; /**< Second derivative type. */

    // Bounds (may not be used)
    InputType m_lower_bound; /**< Lower bound. */
    InputType m_upper_bound; /**< Upper bound. */

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
    Real    m_alpha{0.8};          /**< Relaxation factor \f$ \alpha \f$. */
    Integer m_relaxations{0};      /**< Algorithm relaxations. */
    Integer m_max_relaxations{10}; /**< Maximum allowed algorithm relaxations. */

    // Settings
    Real m_tolerance{EPSILON_LOW};        /**< Solver tolerance \f$ \epsilon \f$ for convergence. */
    bool m_verbose{false};                /**< Verbose mode boolean flag. */
    bool m_damped{true};                  /**< Damped mode boolean flag. */
    std::ostream * m_ostream{&std::cout}; /**< Output stream for verbose mode. */

    // Convergence output flag and trace
    std::string m_task{"Undefined"}; /**< Task name. */
    bool        m_converged{false};  /**< Convergence boolean flag. */
    TraceType   m_trace;             /**< Trace for points \f$ \mathbf{x} \f$ values. */

  public:
    /**
    * Class constructor for the nonlinear solver.
    */
    SolverBase() {
      if constexpr (ForceEigen || SolInDim > 1) {
        this->m_lower_bound.setConstant(-INFTY);
        this->m_upper_bound.setConstant(INFTY);
      } else {
        this->m_lower_bound = -INFTY;
        this->m_upper_bound = INFTY;
      }
      this->m_trace.reserve(this->m_max_iterations * this->m_max_relaxations);
    }

    /**
    * Class constructor for the nonlinear solver.
    * \param[in] function Function wrapper.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    */
    SolverBase(FunctionWrapper function, const InputType & x_ini, InputType & x_sol) : SolverBase() {
      static_cast<const DerivedSolver *>(this)->solve_impl(function, x_ini, x_sol);
    }

    /**
    * Class constructor for the nonlinear solver.
    * \param[in] function Function wrapper.
    * \param[in] first_derivative First derivative wrapper.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    */
    SolverBase(FunctionWrapper function, FirstDerivativeWrapper first_derivative, const InputType & x_ini,
      InputType & x_sol) : SolverBase()
    {
      static_cast<const DerivedSolver *>(this)->solve_impl(function, first_derivative, x_ini, x_sol);
    }

    /**
    * Class constructor for the nonlinear solver.
    * \param[in] function Function wrapper.
    * \param[in] first_derivative First derivative wrapper.
    * \param[in] second_derivative The second derivative wrapper.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    */
    SolverBase(FunctionWrapper function, FirstDerivativeWrapper first_derivative, SecondDerivativeWrapper
      second_derivative, const InputType & x_ini, InputType & x_sol) : SolverBase()
    {
      static_cast<const DerivedSolver *>(this)->solve_impl(function, first_derivative, second_derivative, x_ini, x_sol);
    }

    /**
    * Get the lower bound.
    * \return The lower bound.
    */
    const InputType & lower_bound() const {return this->m_lower_bound;}

    /**
    * Set the lower bound.
    * \param[in] t_lower_bound The lower bound.
    */
    void lower_bound(const InputType & t_lower_bound) {
      #define CMD "Optimist::Solver::bounds(...): "

      if constexpr (ForceEigen || SolInDim > 1) {
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
    const InputType & upper_bound() const {return this->m_upper_bound;}

    /**
    * Set the upper bound.
    * \param[in] t_upper_bound The upper bound.
    */
    void upper_bound(const InputType & t_upper_bound) {
      #define CMD "Optimist::Solver::bounds(...): "

      if constexpr (ForceEigen || SolInDim > 1) {
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
    void bounds(const InputType & t_lower_bound, const InputType & t_upper_bound)
    {
      #define CMD "Optimist::Solver::bounds(...): "

      if constexpr (ForceEigen || SolInDim > 1) {
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
    constexpr Integer input_dimension() const {return SolInDim;}

    /**
    * Get the output dimension of the function.
    * \return The output dimension of the function.
    */
    constexpr Integer output_dimension() const {return SolOutDim;}

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
      OPTIMIST_ASSERT(!std::isnan(t_max_function_evaluations) && std::isfinite(t_max_function_evaluations),
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
      OPTIMIST_ASSERT(
        !std::isnan(first_derivative_evaluations) && std::isfinite(first_derivative_evaluations),
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
      OPTIMIST_ASSERT(
        !std::isnan(second_derivative_evaluations) && std::isfinite(second_derivative_evaluations),
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
      OPTIMIST_ASSERT(!std::isnan(t_max_iterations) && std::isfinite(t_max_iterations),
        "Optimist::Solver::max_iterations(...): invalid input detected.");
      this->m_max_iterations = t_max_iterations;
    }

    /**
    * Get relaxation factor \f$ \alpha \f$.
    * \return The relaxation factor \f$ \alpha \f$.
    */
    Real alpha() const {return this->m_alpha;}

    /**
    * Set relaxation factor \f$ \alpha \f$.
    * \param[in] t_alpha The relaxation factor \f$ \alpha \f$.
    */
    void alpha(Real t_alpha)
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
      OPTIMIST_ASSERT(!std::isnan(t_max_relaxations) && std::isfinite(t_max_relaxations),
        "Optimist::Solver::max_relaxations(...): invalid input detected.");
      this->m_max_relaxations = t_max_relaxations;
    }

    /**
    * Get the tolerance \f$ \epsilon \f$.
    * \return The tolerance \f$ \epsilon \f$.
    */
    Real tolerance() const {return this->m_tolerance;}

    /**
    * Set the tolerance \f$ \epsilon \f$ for which the nonlinear solver stops, i.e., \f$ \left\|
    * \mathbf{F}(\mathbf{x}) \right\|_{2} < \epsilon \f$.
    * \param[in] t_tolerance The tolerance \f$ \epsilon \f$.
    */
    void tolerance(Real t_tolerance) {
      OPTIMIST_ASSERT(
        !std::isnan(t_tolerance) && std::isfinite(t_tolerance) && t_tolerance > 0.0,
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
    * Get the trace of input values during the algorithm iterations.
    * \return The trace of input values.
    */
    const TraceType & trace() const {return this->m_trace;}

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
    * \param[in] function Function wrapper.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    * \return True if the problem is solved, false otherwise.
    */
    bool solve(FunctionWrapper function, const InputType & x_ini, InputType & x_sol)
    {
      #define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::requires_function,
        CMD "the solver requires a function.");
      return static_cast<DerivedSolver *>(this)->solve(function, nullptr, nullptr, x_ini, x_sol);

      #undef CMD
    }

    /**
    * Solve the root-finding/optimization problem given the function, and its first derivative.
    * \param[in] function Function wrapper.
    * \param[in] first_derivative First derivative wrapper
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    * \return True if the problem is solved, false otherwise.
    */
    bool solve(FunctionWrapper function, FirstDerivativeWrapper first_derivative, const InputType & x_ini,
      InputType & x_sol)
    {
      #define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::requires_function,
        CMD "the solver requires a function.");
      static_assert(DerivedSolver::requires_first_derivative,
        CMD "the solver requires the first derivative.");
      return static_cast<DerivedSolver *>(this)->solve(function, first_derivative, nullptr, x_ini, x_sol);

      #undef CMD
    }

    /**
    * Solve the root-finding/optimization problem given the function, and its first and second derivatives.
    * \param[in] function Function wrapper.
    * \param[in] first_derivative First derivative wrapper.
    * \param[in] second_derivative The second derivative wrapper.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    * \return True if the problem is solved, false otherwise.
    */
    bool solve(FunctionWrapper function, FirstDerivativeWrapper first_derivative, SecondDerivativeWrapper
      second_derivative, const InputType & x_ini, InputType & x_sol)
      {
      #define CMD "Optimist::Solver::solve(...): "

      static_assert(DerivedSolver::requires_function,
        CMD "the solver requires the function.");
      static_assert(DerivedSolver::requires_first_derivative,
        CMD "the solver requires the first derivative.");
      static_assert(DerivedSolver::requires_second_derivative,
        CMD "the solver requires the second derivative.");
      return static_cast<DerivedSolver *>(this)->solve(function, first_derivative, second_derivative,
        x_ini, x_sol);

      #undef CMD
    }

    /**
    * Solve the root-finding problem given a function class.
    * \param[in] function Function class.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    * \tparam FunInDim The function output dimension.
    * \tparam FunOutDim The function output dimension.
    * \tparam DerivedFunction Derived function class.
    * \return True if the problem is solved, false otherwise.
    */
    template <Integer FunInDim, Integer FunOutDim, typename DerivedFunction>
    bool rootfind(FunctionBase<Real, FunInDim, FunOutDim, DerivedFunction, ForceEigen && FunOutDim == 1> const & function,
      const InputType & x_ini, InputType & x_sol)
    {
      #define CMD "Optimist::Solver::rootfind(...): "

      static_assert(SolInDim == FunInDim,
        CMD "solver input dimension must be equal to the function input dimension.");
      static_assert(SolOutDim == FunOutDim || SolOutDim == 1,
        CMD "solver output dimension must be equal to the function output dimension or 1.");
      static_assert(!(SolInDim == 1 && DerivedSolver::is_optimizer),
        CMD "one-dimensional optimizers do not support root-finding problems.");
      return this->solve(function, x_ini, x_sol, SolOutDim != FunOutDim || (SolInDim > 1 && SolOutDim == 1));

      #undef CMD
    }

    /**
    * Solve the optimization problem given a function class.
    * \param[in] function Function class.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    * \tparam FunInDim The function output dimension.
    * \tparam FunOutDim The function output dimension.
    * \tparam DerivedFunction Derived function class.
    * \return True if the problem is solved, false otherwise.
    */
    template <Integer FunInDim, Integer FunOutDim, typename DerivedFunction>
    bool optimize(FunctionBase<Real, FunInDim, FunOutDim, DerivedFunction, ForceEigen && FunOutDim == 1> const & function,
      const InputType & x_ini, InputType & x_sol)
    {
      #define CMD "Optimist::Solver::optimize(...): "

      static_assert(SolInDim == FunInDim,
        CMD "solver input dimension must be equal to the function input dimension.");
      static_assert(SolOutDim == 1,
        CMD "solver output dimension must be equal to the function output dimension or 1.");
      static_assert(!(SolInDim == 1 && DerivedSolver::is_rootfinder),
        CMD "one-dimensional root-finders do not support optimization problems.");
      return this->solve(function, x_ini, x_sol, true);

      #undef CMD
    }

    /**
    * Get the solver name.
    * \return The solver name.
    */
    std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();};

  protected:
    /**
    * Solve the root-finding/optimization problem given a function class.
    * \param[in] function Function class.
    * \param[in] x_ini Initialization point.
    * \param[out] x_sol Solution point.
    * \param[in] is_optimization Boolean flag for optimization.
    * \tparam FunInDim The function output dimension.
    * \tparam FunOutDim The function output dimension.
    * \tparam DerivedFunction Derived function class.
    * \return True if the problem is solved, false otherwise.
    */
    template <Integer FunInDim, Integer FunOutDim, typename DerivedFunction>
    bool solve(FunctionBase<Real, FunInDim, FunOutDim, DerivedFunction, ForceEigen && FunOutDim == 1> const & function,
      const InputType & x_ini, InputType & x_sol, bool is_optimization)
    {
      #define CMD "Optimist::Solver::solve(...): "

      using FunctionType = FunctionBase<Real, FunInDim, FunOutDim, DerivedFunction, ForceEigen && FunOutDim == 1>;

      static_assert(SolInDim == FunInDim,
        CMD "solver input dimension must be equal to the function input dimension.");
      static_assert(SolOutDim == FunOutDim || SolOutDim == 1,
        CMD "solver output dimension must be equal to the function output dimension or 1.");

      // Lambda generators for function and derivatives
      auto function_wrapper = [&function, is_optimization](const InputType & x, OutputType & out)
      {
        typename FunctionType::OutputType f;
        function.evaluate(x, f);
        if (is_optimization) {
          if constexpr (FunOutDim == 1) {out = 0.5*f*f;}
          else if constexpr (SolOutDim != FunOutDim) {out = 0.5*f.squaredNorm();}
          else {OPTIMIST_ERROR(CMD "optimization problem with inconsistent output in function.");}
        } else {
          if constexpr (SolOutDim == FunOutDim) {
            out = f;
          }
          else {OPTIMIST_ERROR(CMD "root-finding problem with inconsistent output in function.");}
        }
      };

      auto first_derivative_wrapper = [&function, is_optimization, this](const InputType & x, FirstDerivativeType & out)
      {
        typename FunctionType::FirstDerivativeType J; function.first_derivative(x, J);

        if (is_optimization) {
          typename FunctionType::OutputType f; function.evaluate(x, f);
          this->m_function_evaluations++;
          if constexpr (FunInDim == 1 && FunOutDim == 1) {out = J*f;}
          else if constexpr (SolOutDim != FunOutDim) {out = J.transpose()*f;}
          else {OPTIMIST_ERROR(CMD "optimization problem inconsistent output in first derivative.");}
        } else {
          if constexpr (SolOutDim == FunOutDim) {
            out = J;
          }
          else {OPTIMIST_ERROR(CMD "root-finding problem with inconsistent output in first derivative.");}
        }
      };

      auto second_derivative_wrapper = [&function, is_optimization, this](const InputType & x, SecondDerivativeType & out)
      {
        typename FunctionType::SecondDerivativeType H; function.second_derivative(x, H);

        if (is_optimization) {
          typename FunctionType::OutputType f;          function.evaluate(x, f);
          typename FunctionType::FirstDerivativeType J; function.first_derivative(x, J);
          this->m_function_evaluations++;               this->m_first_derivative_evaluations++;
          if constexpr (FunInDim == 1 && FunOutDim == 1) {out = J*J + f*H;}
          else if constexpr (SolOutDim != FunOutDim) {
            out = J.transpose()*J;
            for (Integer i = 0; i < static_cast<Integer>(H.size()); ++i) {out += f(i)*H[i];}
          }
          else {OPTIMIST_ERROR(CMD "optimization problem with inconsistent output in second derivative.");}
        } else {
          if constexpr (SolOutDim == FunOutDim) {out = H;}
          else {OPTIMIST_ERROR(CMD "root-finding problem with inconsistent output in second derivative.");}
        }
      };

      // Select solver method based on derivative requirements
      if constexpr (DerivedSolver::requires_function && !DerivedSolver::requires_first_derivative &&
        !DerivedSolver::requires_second_derivative) {
        return static_cast<DerivedSolver *>(this)->solve(function_wrapper, x_ini, x_sol);
      } else if constexpr (DerivedSolver::requires_function && DerivedSolver::requires_first_derivative &&
        !DerivedSolver::requires_second_derivative) {
        return static_cast<DerivedSolver *>(this)->solve(function_wrapper, first_derivative_wrapper,
          x_ini, x_sol);
      } else if constexpr (DerivedSolver::requires_function && DerivedSolver::requires_first_derivative &&
        DerivedSolver::requires_second_derivative) {
        return static_cast<DerivedSolver *>(this)->solve(function_wrapper, first_derivative_wrapper,
          second_derivative_wrapper, x_ini, x_sol);
      } else {
        OPTIMIST_ERROR(CMD "no matching function signature found for the solver.");
      }

      #undef CMD
    }

    /**
    * Reset solver internal counters and variables.
    */
    void reset()
    {
      this->m_function_evaluations          = static_cast<Integer>(0);
      this->m_first_derivative_evaluations  = static_cast<Integer>(0);
      this->m_second_derivative_evaluations = static_cast<Integer>(0);
      this->m_iterations                    = static_cast<Integer>(0);
      this->m_relaxations                   = static_cast<Integer>(0);
      this->m_converged                     = false;
      this->m_trace.clear();
    }

    /**
    * Evaluate the function.
    * \param[in] function Function wrapper.
    * \param[in] x Input point.
    * \param[out] out Function value.
    */
    void evaluate_function(FunctionWrapper function, const InputType & x, OutputType & out)
    {
      OPTIMIST_ASSERT(this->m_function_evaluations < this->m_max_function_evaluations,
        "Optimist::" << this-> name() << "::evaluate_function(...): maximum allowed function evaluations reached.");
      ++this->m_function_evaluations;
      function(x, out);
    }

    /**
    * Evaluate the first derivative.
    * \param[in] function First derivative wrapper.
    * \param[in] x Input point.
    * \param[out] out First derivative value.
    */
    void evaluate_first_derivative(FirstDerivativeWrapper function, const InputType & x, FirstDerivativeType & out)
    {
      OPTIMIST_ASSERT(this->m_first_derivative_evaluations < this->m_max_first_derivative_evaluations,
        "Optimist::" << this-> name() << "::evaluate_first_derivative(...): maximum allowed first derivative evaluations reached.");
      ++this->m_first_derivative_evaluations;
      function(x, out);
    }

    /**
    * Evaluate the second derivative.
    * \param[in] function Second derivative wrapper.
    * \param[in] x Input point.
    * \param[out] out Second derivative value.
    */
    void evaluate_second_derivative(SecondDerivativeWrapper function, const InputType & x, SecondDerivativeType & out)
    {
      OPTIMIST_ASSERT(this->m_second_derivative_evaluations < this->m_max_second_derivative_evaluations,
        "Optimist::" << this-> name() << "::evaluate_second_derivative(...): maximum allowed second derivative evaluations reached.");
      ++this->m_second_derivative_evaluations;
      function(x, out);
    }

    /**
    * Update the history of the solver with the current point and function value.
    * \param[in] x The point \f$ \mathbf{x} \f$.
    */
    void store_trace(const InputType & x) {this->m_trace.push_back(x);}

    /**
    * Damp the step using the affine invariant criterion.
    * \param[in] function Function wrapper.
    * \param[in] x_old Old point.
    * \param[in] function_old Old function value.
    * \param[in] step_old Old step.
    * \param[in] x_new New point.
    * \param[in] function_new New function value.
    * \param[out] step_new New step.
    * \return The damping boolean flag, true if the damping is successful, false otherwise.
    */
    bool damp(FunctionWrapper function, InputType const & x_old, InputType const & function_old,
      InputType const & step_old, InputType & x_new, InputType & function_new, InputType & step_new)
    {
      Real step_norm_old, step_norm_new, residuals_old, residuals_new, tau{1.0};
      for (this->m_relaxations = static_cast<Integer>(0); this->m_relaxations < this->m_max_relaxations; ++this->m_relaxations)
      {
        // Update point
        step_new = tau * step_old;
        x_new = x_old + step_new;
        this->evaluate_function(function, x_new, function_new);

        // Check relaxation
        if constexpr (ForceEigen || (SolInDim > 1 && SolOutDim > 1)) {
          residuals_old = function_old.norm();
          residuals_new = function_new.norm();
          step_norm_old = step_old.norm();
          step_norm_new = step_new.norm();
        } else {
          residuals_old = std::abs(function_old);
          residuals_new = std::abs(function_new);
          step_norm_old = std::abs(step_old);
          step_norm_new = std::abs(step_new);
        }
        if (residuals_new < residuals_old || step_norm_new < (Real(1.0)-tau/Real(2.0))*step_norm_old) {
          return true;
        } else {
          tau *= this->m_alpha;
        }
      }
      return false;
    }

    /**
    * Print the table header solver information.
    * \note This has to be properly placed in the derived classes.
    */
    void header()
    {
      std::string c_tl{table_top_left_corner()};
      std::string c_tr{table_top_right_corner()};
      std::string h_7{table_horizontal_line<7>()};
      std::string h_14{table_horizontal_line<14>()};
      std::string h_23{table_horizontal_line<23>()};
      std::string h_78{table_horizontal_line<78>()};
      std::string v_ll{table_vertical_line() + " "};
      std::string v_rr{" " + table_vertical_line()};
      std::string v_lc{" " + table_vertical_line() + " "};
      std::string j_tt{table_top_junction()};
      std::string j_cc{table_center_cross()};
      std::string j_ll{table_left_junction()};
      std::string j_rr{table_right_junction()};

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
      std::string c_bl{table_bottom_left_corner()};
      std::string c_br{table_bottom_right_corner()};
      std::string h_7{table_horizontal_line<7>()};
      std::string h_14{table_horizontal_line<14>()};
      std::string h_23{table_horizontal_line<23>()};
      std::string h_78{table_horizontal_line<78>()};
      std::string v_ll{table_vertical_line() + " "};
      std::string v_rr{" " + table_vertical_line()};
      std::string j_ll{table_left_junction()};
      std::string j_rr{table_right_junction()};
      std::string j_bb{table_bottom_junction()};

      *this->m_ostream
        << j_ll << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb << h_7 << j_bb << h_14 << j_bb << h_23 << j_rr << std::endl
        << v_ll << std::setw(40) << (this->m_converged ? "CONVERGED" : "NOT CONVERGED") << std::setw(40) << v_rr << std::endl
        << c_bl << h_78 << c_br << std::endl;
    }

    /**
    * Print the solver information during the algorithm iterations.
    * \note This has to be properly placed in the derived classes.
    */
    void info(Real residuals, std::string const & notes = "-")
    {
      std::string v_rr{" " + table_vertical_line()};
      std::string v_ll{table_vertical_line() + " "};
      std::string v_lc{" " + table_vertical_line() + " "};

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
