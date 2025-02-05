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

#ifndef OPTIMIST_SOLVER_HXX
#define OPTIMIST_SOLVER_HXX

namespace Optimist
{

  /*\
   |   ____        _
   |  / ___|  ___ | |_   _____ _ __
   |  \___ \ / _ \| \ \ / / _ \ '__|
   |   ___) | (_) | |\ V /  __/ |
   |  |____/ \___/|_| \_/ \___|_|
   |
  \*/

  /**
  * \brief Class container for the generic root-finding/optimization problem solver.
  *
  * \includedoc docs/markdown/Solver.md
  *
  * \tparam InputDim The root-finding/optimization problem input dimension.
  * \tparam OutputDim The root-finding/optimization problem output dimension.
  */
  template <Integer InputDim, Integer OutputDim>
  class Solver
  {
    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(InputDim > Integer(0) && OutputDim > Integer(0),
      "Negative-dimensional optimization problem? Are you serious?");

  protected:
    // I/O types
    using InputType  = typename std::conditional_t<InputDim == 1,  Real, Eigen::Vector<Real, InputDim>>;  /**< Input type. */
    using OutputType = typename std::conditional_t<OutputDim == 1, Real, Eigen::Vector<Real, OutputDim>>; /**< Output type. */

    // Trace types
    using TraceType = typename std::vector<InputType>; /**< Input trace type. */

    // Derivative types
    using FirstDerivativeType  = std::conditional_t<InputDim == 1 && OutputDim == 1, Real, Eigen::Matrix<Real, OutputDim, InputDim>>; /**< First derivative type. */
    using SecondDerivativeType = std::conditional_t<InputDim == 1 && OutputDim == 1, Real, Eigen::Matrix<Real, InputDim, InputDim>>;  /**< Second derivative type. */

    // Function types
    using Function         = typename std::function<void(const InputType &, OutputType &)>;           /**< Function type. */
    using FirstDerivative  = typename std::function<void(const InputType &, FirstDerivativeType &)>;  /**< First derivative type. */
    using SecondDerivative = typename std::function<void(const InputType &, SecondDerivativeType &)>; /**< Second derivative type. */

    // Function pointers
    Function         m_function{nullptr};          /**< Function pointer. */
    FirstDerivative  m_first_derivative{nullptr};  /**< First derivative pointer. */
    SecondDerivative m_second_derivative{nullptr}; /**< Second derivative pointer. */

    // Evaluations
    Integer m_function_evaluations{0};          /**< Function evaluations. */
    Integer m_first_derivative_evaluations{0};  /**< First derivative evaluations. */
    Integer m_second_derivative_evaluations{0}; /**< Second derivative evaluations. */

    // Maximum allowed evaluations
    Integer m_max_function_evaluations{100};          /**< Maximum allowed function evaluations. */
    Integer m_max_first_derivative_evaluations{100};  /**< Maximum allowed first derivative evaluations. */
    Integer m_max_second_derivative_evaluations{100}; /**< Maximum allowed second derivative evaluations. */

    // Iterations and relaxations
    Integer m_iterations{0};       /**< Algorithm iterations. */
    Integer m_max_iterations{100}; /**< Maximum allowed algorithm iterations. */
    Real    m_alpha{0.8};          /**< Relaxation factor \f$ \alpha \f$. */
    Integer m_relaxations{0};      /**< Algorithm relaxations. */
    Integer m_max_relaxations{10}; /**< Maximum allowed algorithm relaxations. */

    // Settings
    Real m_tolerance{EPSILON_HIGH};       /**< Solver tolerance \f$ \epsilon \f$ for convergence. */
    bool m_verbose{true};                 /**< Verbose mode boolean flag. */
    bool m_damped{true};                  /**< Damped mode boolean flag. */
    std::ostream * m_ostream{&std::cout}; /**< Output stream for verbose mode. */

    // Convergence output flag and trace
    bool      m_converged{false}; /**< Convergence boolean flag. */
    TraceType m_trace;            /**< Trace for points \f$ \mathbf{x} \f$ values. */

    // Ouput stream for verbose mode

  public:
    /**
    * Class constructor for the nonlinear solver.
    */
    Solver() {this->m_trace.reserve(this->m_max_iterations * this->m_max_relaxations);}

    /**
    * Class constructor for the nonlinear solver.
    * \param[in] t_function The function.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    Solver(Function t_function, const InputType & x_ini, OutputType & x_sol) : Solver() {
      this->solve(t_function, x_ini, x_sol);
    }

    /**
    * Class constructor for the nonlinear solver.
    * \param[in] t_function The function.
    * \param[in] t_first_derivative The first derivative of the function.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    Solver(Function t_function, FirstDerivative t_first_derivative, const InputType & x_ini,
      OutputType & x_sol) : Solver()
    {
      this->solve(t_function, t_first_derivative, x_ini, x_sol);
    }

    /**
    * Class constructor for the nonlinear solver.
    * \param[in] t_function The function.
    * \param[in] t_first_derivative The first derivative of the function.
    * \param[in] t_second_derivative The second derivative of the function.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    Solver(Function t_function, FirstDerivative t_first_derivative, SecondDerivative t_second_derivative,
      const InputType & x_ini, OutputType & x_sol) : Solver()
    {
      this->solve(t_function, t_first_derivative, t_second_derivative, x_ini, x_sol);
    }

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
    * Get the number of function first derivative evaluations.
    * \return The number of function first derivative evaluations.
    */
    Integer first_derivative_evaluations() const {return this->m_first_derivative_evaluations;}

    /**
    * Get the number of maximum allowed first derivative evaluations.
    * \return The number of maximum allowed first derivative evaluations.
    */
    Integer max_first_derivative_evaluations() const {return this->m_max_first_derivative_evaluations;}

    /**
    * Set the number of maximum allowed first derivative evaluations.
    * \param[in] t_first_derivative_evaluations The number of maximum allowed first derivative evaluations.
    */
    void max_first_derivative_evaluations(Integer t_first_derivative_evaluations)
    {
      OPTIMIST_ASSERT(
        !std::isnan(t_first_derivative_evaluations) && std::isfinite(t_first_derivative_evaluations),
        "Optimist::Solver::max_first_derivative_evaluations(...): invalid input detected.");
      this->m_max_first_derivative_evaluations = t_first_derivative_evaluations;
    }

    /**
    * Get the number of function second derivative evaluations.
    * \return The number of function second derivative evaluations.
    */
    Integer second_derivative_evaluations() const {return this->m_second_derivative_evaluations;}

    /**
    * Get the number of maximum allowed second derivative evaluations.
    * \return The number of maximum allowed second derivative evaluations.
    */
    Integer max_second_derivative_evaluations() const {return this->m_max_second_derivative_evaluations;}

    /**
    * Set the number of maximum allowed second derivative evaluations.
    * \param[in] t_second_derivative_evaluations The number of maximum allowed second derivative evaluations.
    */
    void max_second_derivative_evaluations(Integer t_second_derivative_evaluations)
    {
      OPTIMIST_ASSERT(
        !std::isnan(t_second_derivative_evaluations) && std::isfinite(t_second_derivative_evaluations),
        "Optimist::Solver::max_second_derivative_evaluations(...): invalid input detected.");
      this->m_max_second_derivative_evaluations = t_second_derivative_evaluations;
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
      OPTIMIST_ASSERT(!std::isnan(t_alpha) && std::isfinite(t_alpha) && Real(0.0) <= t_alpha && t_alpha <= Real(1.0),
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
        !std::isnan(t_tolerance) && std::isfinite(t_tolerance) && t_tolerance > Real(0.0),
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
    * Solve the root-finding/optimization problem without derivatives given the function.
    * \param[in] t_function The function.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    bool solve(Function t_function, const InputType & x_ini, OutputType & x_sol)
    {
      this->m_function          = t_function;
      this->m_first_derivative  = nullptr;
      this->m_second_derivative = nullptr;
      OPTIMIST_ASSERT(this->check(), "Optimist::Solver::solve(...): in solver " << this->name() <<
        ", insufficient information is provided.");
      return this->solve(x_ini, x_sol);
    }

    /**
    * Solve the root-finding/optimization problem given the function and its first derivative.
    * \param[in] t_function The function.
    * \param[in] t_first_derivative The first derivative of the function.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    bool solve(Function t_function, FirstDerivative t_first_derivative, const InputType & x_ini,
      OutputType & x_sol)
    {
      this->m_function          = t_function;
      this->m_first_derivative  = t_first_derivative;
      this->m_second_derivative = nullptr;
      OPTIMIST_ASSERT(this->check(), "Optimist::Solver::solve(...): in solver " << this->name() <<
        ", insufficient information is provided.");
      return this->solve(x_ini, x_sol);
    }

    /**
    * Solve the root-finding/optimization problem given the function and its first and second
    * derivatives.
    * \param[in] t_function The function.
    * \param[in] t_first_derivative The first derivative of the function.
    * \param[in] t_second_derivative The second derivative of the function.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    bool solve(Function t_function, FirstDerivative t_first_derivative, SecondDerivative t_second_derivative,
      const InputType & x_ini, OutputType & x_sol)
    {
      this->m_function          = t_function;
      this->m_first_derivative  = t_first_derivative;
      this->m_second_derivative = t_second_derivative;
      OPTIMIST_ASSERT(this->check(), "Optimist::Solver::solve(...): in solver " << this->name() <<
        ", insufficient information is provided.");
      return this->solve(x_ini, x_sol);
    }

    /**
    * Get the solver name.
    * \return The solver name.
    */
    virtual std::string name() const = 0;

    /**
    * Check if the solver is able to solve the problem with the given input.
    * \return The check boolean flag.
    */
    virtual bool check() const = 0;

  protected:
    /**
    * Reset solver internal counters and variables.
    */
    void reset()
    {
      this->m_function_evaluations          = Integer(0);
      this->m_first_derivative_evaluations  = Integer(0);
      this->m_second_derivative_evaluations = Integer(0);
      this->m_iterations                    = Integer(0);
      this->m_relaxations                   = Integer(0);
      this->m_converged                     = false;
      this->m_trace.clear();
    }

    /**
    * Evaluate the function.
    * \param[in] x Input point.
    * \param[out] function Function value.
    */
    void evaluate_function(const InputType & x, OutputType & function)
    {
      ++this->m_function_evaluations;
      this->m_function(x, function);
    }

    /**
    * Evaluate the first derivative.
    * \param[in] x Input point.
    * \param[out] first_derivative First derivative value.
    */
    void evaluate_first_derivative(const InputType & x, FirstDerivativeType & first_derivative)
    {
      ++this->m_first_derivative_evaluations;
      this->m_first_derivative(x, first_derivative);
    }

    /**
    * Evaluate the second derivative.
    * \param[in] x Input point.
    * \param[out] second_derivative Second derivative value.
    */
    void evaluate_second_derivative(const InputType & x, SecondDerivativeType & second_derivative)
    {
      ++this->m_second_derivative_evaluations;
      this->m_second_derivative(x, second_derivative);
    }

    /**
    * Update the history of the solver with the current point and function value.
    * \param[in] x The point \f$ \mathbf{x} \f$.
    */
    void store_trace(const InputType & x) {this->m_trace.push_back(x);}

    /**
    * Damp the step using the affine invariant criterion.
    * \param[in] x_old Old point.
    * \param[in] function_old Old function value.
    * \param[in] step_old Old step.
    * \param[out] x_new New point.
    * \param[out] function_new New function value.
    * \param[out] step_new New step.
    * \return The damping boolean flag, true if the damping is successful, false otherwise.
    */
    bool damp(InputType const & x_old, InputType const & function_old, InputType const & step_old,
      InputType & x_new, InputType & function_new, InputType & step_new)
    {
      Real step_norm_old, step_norm_new, residuals_old, residuals_new, tau{1.0};
      for (this->m_relaxations = Integer(0); this->m_relaxations < this->m_max_relaxations; ++this->m_relaxations)
      {
        // Update point
        step_new = tau * step_old;
        x_new = x_old + step_new;
        this->evaluate_function(x_new, function_new);

        // Check relaxation
        if constexpr (InputDim == 1 && OutputDim == 1) {
          residuals_old = std::abs(function_old);
          residuals_new = std::abs(function_new);
          step_norm_old = std::abs(step_old);
          step_norm_new = std::abs(step_new);
        } else {
          residuals_old = function_old.norm();
          residuals_new = function_new.norm();
          step_norm_old = step_old.norm();
          step_norm_new = step_new.norm();
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
      *this->m_ostream
        << "Solver Name: " << this->name() << std::endl
        << CTL << H14 << TU << H7 << TU << H7 << TU << H7 << TU << H7 << CTR << std::endl
        << VL << "   ║f(x)║₂  " << VC << "#Iter" << VC << "   #f" << VC << "  #Df" << VC << " #DDf" << VC << "Additional notes" << std::endl
        << TL << H14 << C << H7 << C << H7 << C << H7 << C << H7 << TR << std::endl;
    }

    /**
    * Print the table bottom solver information.
    * \note This has to be properly placed in the derived classes.
    */
    void bottom()
    {
      *this->m_ostream
        << CBL << H14 << TD << H7 << TD << H7 << TD << H7 << TD << H7 << CBR << std::endl
        << this->name() << ": " << (this->m_converged ? "CONVERGED" : "NOT CONVERGED") << std::endl;
    }

    /**
    * Print the solver information during the algorithm iterations.
    * \note This has to be properly placed in the derived classes.
    */
    void info(Real residuals, std::string notes = "")
    {
      if (this->m_verbose)
      {
        *this->m_ostream << VL
          << std::setw(12) << std::scientific << std::setprecision(6) << residuals << VC
          << std::setw(5) << this->m_iterations << VC
          << std::setw(5) << this->m_function_evaluations << VC
          << std::setw(5) << this->m_first_derivative_evaluations << VC
          << std::setw(5) << this->m_second_derivative_evaluations << VC
          << notes << std::endl;
      }
    }

    /**
    * Solve root-finding/optimization problem according to the solver.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    * \return The convergence boolean flag.
    */
    virtual bool solve(const InputType & x_ini, OutputType &x_sol) = 0;

  }; // class Solver

} // namespace Optimist

#endif // OPTIMIST_SOLVER_HXX
