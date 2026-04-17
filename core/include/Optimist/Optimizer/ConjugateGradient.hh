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

#ifndef OPTIMIST_OPTIMIZER_CONJUGATE_GRADIENT_HH
#define OPTIMIST_OPTIMIZER_CONJUGATE_GRADIENT_HH

#include "Optimist/Optimizer.hh"

namespace Optimist {
  namespace Optimizer {

    /**
     * \brief Class container for nonlinear conjugate gradient methods.
     *
     * \includedoc docs/markdown/Optimizer/ConjugateGradient.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class ConjugateGradient
        : public Optimizer<Vector, ConjugateGradient<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename TypeTrait<Vector>::Scalar;
      using typename Optimizer<Vector,
                               ConjugateGradient<Vector>>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Conjugate gradient update method enumeration.
       */
      using BetaMethod = enum class BetaMethod : Integer {
        FLETCHER_REEVES    = 0,
        POLAK_RIBIERE      = 1,
        POLAK_RIBIERE_PLUS = 2,
        HESTENES_STIEFEL   = 3,
        CONJUGATE_DESCENT  = 4,
        LIU_STOREY         = 5,
        DAI_YUAN           = 6,
        HAGER_ZHANG        = 7,
        HAGER_ZHANG_PLUS   = 8
      };

      /**
       * Step-length initialization strategy enumeration.
       */
      using AlphaMethod = enum class AlphaMethod : Integer {
        ACCEPTED_STEP    = 0,
        BARZILAI_BORWEIN = 1
      };

     private:
      BetaMethod m_beta_method{
        BetaMethod::HAGER_ZHANG};      /**< Conjugate gradient method. */
      AlphaMethod m_alpha_method{
        AlphaMethod::ACCEPTED_STEP};   /**< Trial step initialization method. */

      Scalar m_initial_step{1.0};      /**< Initial line-search step. */
      Scalar m_min_step{SQRT_EPSILON}; /**< Minimum admissible step. */
      Scalar m_max_step{static_cast<Scalar>(1.0) /
                        SQRT_EPSILON}; /**< Maximum admissible step. */
      Scalar m_step_expand{2.0};       /**< Line-search expansion factor. */
      Scalar m_step_shrink{0.5};       /**< Line-search shrink factor. */
      Scalar m_armijo_parameter{
        1.0e-4}; /**< Armijo sufficient decrease constant. */
      Scalar m_descent_parameter{1.0e-2}; /**< Sufficient descent constant. */
      Scalar m_curvature_parameter{0.1};  /**< Wolfe curvature constant. */
      Scalar m_truncation_parameter{
        1.0e-2}; /**< Truncation parameter for Hager-Zhang+. */
      Integer m_max_line_search_iterations{
        25};     /**< Maximum line-search iterations. */

      bool m_powell_restart{true};     /**< Powell restart flag. */
      Scalar m_restart_tolerance{0.9}; /**< Powell restart tolerance. */

     private:
      /**
       * Get the lower truncation bound used by the Hager-Zhang+ method.
       * \param[in] gradient_new Current gradient.
       * \param[in] direction_old Previous search direction.
       * \return The lower truncation bound.
       */
      Scalar hager_zhang_plus_bound(const Vector &gradient_new,
                                    const Vector &direction_old) const {
        Scalar direction_norm{direction_old.norm()};
        Scalar gradient_norm{gradient_new.norm()};
        Scalar denominator{
          direction_norm *
          std::min(this->m_truncation_parameter, gradient_norm)};
        if (!std::isfinite(denominator) || denominator <= SQRT_EPSILON) {
          return 0.0;
        }
        return -1.0 / denominator;
      }

      /**
       * Clamp the line-search step inside the configured interval.
       * \param[in] t_step Trial step.
       * \return Clamped step.
       */
      Scalar clamp_step(Scalar t_step) const {
        if (!std::isfinite(t_step) || t_step <= 0.0) {
          t_step = this->m_initial_step;
        }
        return std::min(this->m_max_step, std::max(this->m_min_step, t_step));
      }

      /**
       * Compute the next trial step using the Barzilai-Borwein formula.
       *
       * The update uses the accepted step along the previous search direction
       * and falls back to the accepted step when the curvature estimate is not
       * positive or is numerically unreliable.
       *
       * \param[in] accepted_step Accepted line-search step.
       * \param[in] direction_old Previous search direction.
       * \param[in] gradient_old Previous gradient.
       * \param[in] gradient_new Current gradient.
       * \return The next trial step.
       */
      Scalar barzilai_borwein_step(const Scalar accepted_step,
                                   const FirstDerivative &direction_old,
                                   const FirstDerivative &gradient_old,
                                   const FirstDerivative &gradient_new) const {
        const Scalar fallback_step{this->clamp_step(accepted_step)};
        const FirstDerivative y{gradient_new - gradient_old};
        const Scalar denominator{direction_old.dot(y)};
        if (!std::isfinite(denominator) || denominator <= SQRT_EPSILON) {
          return fallback_step;
        }

        const Scalar numerator{accepted_step * direction_old.squaredNorm()};
        if (!std::isfinite(numerator) || numerator <= SQRT_EPSILON) {
          return fallback_step;
        }

        return this->clamp_step(numerator / denominator);
      }

      /**
       * Compute the next trial step passed to the line search.
       * \param[in] accepted_step Accepted line-search step.
       * \param[in] direction_old Previous search direction.
       * \param[in] gradient_old Previous gradient.
       * \param[in] gradient_new Current gradient.
       * \return The next trial step.
       */
      Scalar next_trial_step(const Scalar accepted_step,
                             const FirstDerivative &direction_old,
                             const FirstDerivative &gradient_old,
                             const FirstDerivative &gradient_new) const {
        switch (this->m_alpha_method) {
          case AlphaMethod::BARZILAI_BORWEIN:
            return this->barzilai_borwein_step(accepted_step,
                                               direction_old,
                                               gradient_old,
                                               gradient_new);

          case AlphaMethod::ACCEPTED_STEP:
          default:
            return this->clamp_step(accepted_step);
        }
      }

      /**
       * Check whether the current direction is a descent direction.
       * \param[in] gradient Current gradient.
       * \param[in] direction Current search direction.
       * \return True if the direction is a descent direction, false otherwise.
       */
      bool is_descent_direction(const Vector &gradient,
                                const Vector &direction) const {
        const Scalar directional_derivative{gradient.dot(direction)};
        const Scalar threshold{
          -SQRT_EPSILON *
          std::max(static_cast<Scalar>(1.0), gradient.squaredNorm())};
        return std::isfinite(directional_derivative) &&
               directional_derivative < threshold;
      }

      /**
       * Check whether the current direction satisfies a sufficient descent
       * condition.
       * \param[in] gradient Current gradient.
       * \param[in] direction Current search direction.
       * \return True if the direction satisfies a sufficient descent bound,
       * false otherwise.
       */
      bool has_sufficient_descent(const Vector &gradient,
                                  const Vector &direction) const {
        const Scalar directional_derivative{gradient.dot(direction)};
        const Scalar threshold{-this->m_descent_parameter *
                               gradient.squaredNorm()};
        return std::isfinite(directional_derivative) &&
               directional_derivative <= threshold;
      }

      /**
       * Check the Powell restart condition.
       * \param[in] gradient_old Previous gradient.
       * \param[in] gradient_new Current gradient.
       * \return True if the search direction must be restarted.
       */
      bool should_restart(const Vector &gradient_old,
                          const Vector &gradient_new) const {
        if (!this->m_powell_restart) {
          return false;
        }
        const Scalar gradient_old_norm{gradient_old.norm()};
        const Scalar gradient_new_norm{gradient_new.norm()};
        if (gradient_old_norm <= SQRT_EPSILON ||
            gradient_new_norm <= SQRT_EPSILON) {
          return false;
        }
        return gradient_new.dot(gradient_old) <= -this->m_restart_tolerance *
                                                     gradient_new_norm *
                                                     gradient_old_norm;
      }

      /**
       * Compute the conjugate gradient update parameter.
       * \param[in] gradient_old Previous gradient.
       * \param[in] gradient_new Current gradient.
       * \param[in] direction_old Previous search direction.
       * \param[in] y Gradient difference.
       * \return The update parameter \f$\beta_k\f$.
       */
      Scalar compute_beta(const Vector &gradient_old,
                          const Vector &gradient_new,
                          const Vector &direction_old,
                          const Vector &y) const {
        const Scalar gradient_old_norm_sq{gradient_old.squaredNorm()};
        const Scalar gradient_new_norm_sq{gradient_new.squaredNorm()};
        const Scalar direction_gradient{-direction_old.dot(gradient_old)};
        const Scalar direction_y{direction_old.dot(y)};
        const Scalar gradient_y{gradient_new.dot(y)};

        switch (this->m_beta_method) {
          case BetaMethod::FLETCHER_REEVES:
            if (gradient_old_norm_sq <= SQRT_EPSILON) {
              return 0.0;
            }
            return gradient_new_norm_sq / gradient_old_norm_sq;

          case BetaMethod::POLAK_RIBIERE:
            if (gradient_old_norm_sq <= SQRT_EPSILON) {
              return 0.0;
            }
            return gradient_y / gradient_old_norm_sq;

          case BetaMethod::POLAK_RIBIERE_PLUS:
            if (gradient_old_norm_sq <= SQRT_EPSILON) {
              return 0.0;
            }
            return std::max(gradient_y / gradient_old_norm_sq,
                            static_cast<Scalar>(0.0));

          case BetaMethod::HESTENES_STIEFEL:
            if (std::abs(direction_y) <= SQRT_EPSILON) {
              return 0.0;
            }
            return gradient_y / direction_y;

          case BetaMethod::CONJUGATE_DESCENT:
            if (direction_gradient <= SQRT_EPSILON) {
              return 0.0;
            }
            return gradient_new_norm_sq / direction_gradient;

          case BetaMethod::LIU_STOREY:
            if (direction_gradient <= SQRT_EPSILON) {
              return 0.0;
            }
            return gradient_y / direction_gradient;

          case BetaMethod::DAI_YUAN:
            if (std::abs(direction_y) <= SQRT_EPSILON) {
              return 0.0;
            }
            return gradient_new_norm_sq / direction_y;

          case BetaMethod::HAGER_ZHANG: {
            if (std::abs(direction_y) <= SQRT_EPSILON) {
              return 0.0;
            }
            Vector tmp(direction_old);
            tmp *= 2.0 * y.squaredNorm() / direction_y;
            tmp = y - tmp;
            return gradient_new.dot(tmp) / direction_y;
          }

          case BetaMethod::HAGER_ZHANG_PLUS: {
            if (std::abs(direction_y) <= SQRT_EPSILON) {
              return 0.0;
            }
            Vector tmp(direction_old);
            tmp *= 2.0 * y.squaredNorm() / direction_y;
            tmp = y - tmp;
            Scalar beta_hz{gradient_new.dot(tmp) / direction_y};
            return std::max(
                beta_hz,
                this->hager_zhang_plus_bound(gradient_new, direction_old));
          }

          default:
            return 0.0;
        }
      }

      /**
       * Perform a standard Wolfe line search.
       * \tparam FunctionLambda Function lambda type.
       * \tparam GradientLambda Gradient lambda type.
       * \param[in] function Function lambda.
       * \param[in] gradient Gradient lambda.
       * \param[in] x_old Current point.
       * \param[in] function_old Current objective value.
       * \param[in] gradient_old Current gradient.
       * \param[in] direction_old Current search direction.
       * \param[in] trial_step Initial line-search step.
       * \param[out] x_new Updated point.
       * \param[out] function_new Updated objective value.
       * \param[out] gradient_new Updated gradient.
       * \param[out] accepted_step Accepted line-search step.
       * \return The boolean flag for successful line search.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool line_search(FunctionLambda &function,
                       GradientLambda &gradient,
                       const Vector &x_old,
                       const Scalar function_old,
                       const FirstDerivative &gradient_old,
                       const Vector &direction_old,
                       const Scalar trial_step,
                       Vector &x_new,
                       Scalar &function_new,
                       FirstDerivative &gradient_new,
                       Scalar &accepted_step) {
#define CMD "Optimist::Optimizer::ConjugateGradient::line_search(...): "

        Scalar directional_derivative{gradient_old.dot(direction_old)};
        OPTIMIST_ASSERT(directional_derivative < 0.0,
                        CMD "non-descent search direction detected.");

        accepted_step = this->clamp_step(trial_step);
        Scalar lower_step{0.0};
        Scalar upper_step{INFTY};
        for (this->m_relaxations = 0;
             this->m_relaxations < this->m_max_line_search_iterations;
             ++this->m_relaxations) {
          x_new = x_old + accepted_step * direction_old;

          bool success{this->evaluate_function(function, x_new, function_new)};
          OPTIMIST_ASSERT(success,
                          CMD "function evaluation failed during line search.");

          if (function_new > function_old + this->m_armijo_parameter *
                                                accepted_step *
                                                directional_derivative) {
            upper_step = accepted_step;
          } else {
            success = this->evaluate_gradient(gradient, x_new, gradient_new);
            OPTIMIST_ASSERT(success,
                            CMD
                            "gradient evaluation failed during line search.");

            Scalar directional_derivative_new{gradient_new.dot(direction_old)};
            if (std::abs(directional_derivative_new) <=
                -this->m_curvature_parameter * directional_derivative) {
              return true;
            }

            if (directional_derivative_new >= 0.0) {
              upper_step = accepted_step;
            } else {
              lower_step = accepted_step;
            }
          }

          if (std::isfinite(upper_step)) {
            if (lower_step <= this->m_min_step) {
              accepted_step *= this->m_step_shrink;
            } else {
              accepted_step =
                  static_cast<Scalar>(0.5) * (lower_step + upper_step);
            }
          } else {
            accepted_step =
                this->clamp_step(accepted_step * this->m_step_expand);
          }
          if (accepted_step < this->m_min_step ||
              (std::isfinite(upper_step) &&
               upper_step - lower_step < this->m_min_step)) {
            break;
          }
        }
        return false;

#undef CMD
      }

     public:
      /**
       * Class constructor for the ConjugateGradient solver.
       */
      ConjugateGradient() {}

      /**
       * Get the ConjugateGradient solver name.
       * \return The ConjugateGradient solver name.
       */
      constexpr std::string name_impl() const {
        return "ConjugateGradient";
      }

      /**
       * Get the enumeration type of the conjugate gradient mbeta ethod.
       * \return The conjugate gradient method.
       */
      BetaMethod beta_method() const {
        return this->m_beta_method;
      }

      /**
       * Set the enumeration type of the conjugate gradient beta method.
       * \param[in] t_beta_method Conjugate gradient method.
       */
      void beta_method(const BetaMethod t_beta_method) {
        this->m_beta_method = t_beta_method;
      }

      /**
       * Enable the Fletcher-Reeves method.
       */
      void enable_fletcher_reeves_method() {
        this->m_beta_method = BetaMethod::FLETCHER_REEVES;
      }

      /**
       * Enable the Polak-Ribiere method.
       */
      void enable_polak_ribiere_method() {
        this->m_beta_method = BetaMethod::POLAK_RIBIERE;
      }

      /**
       * Enable the Polak-Ribiere+ method.
       */
      void enable_polak_ribiere_plus_method() {
        this->m_beta_method = BetaMethod::POLAK_RIBIERE_PLUS;
      }

      /**
       * Enable the Hestenes-Stiefel method.
       */
      void enable_hestenes_stiefel_method() {
        this->m_beta_method = BetaMethod::HESTENES_STIEFEL;
      }

      /**
       * Enable the conjugate descent method.
       */
      void enable_conjugate_descent_method() {
        this->m_beta_method = BetaMethod::CONJUGATE_DESCENT;
      }

      /**
       * Enable the Liu-Storey method.
       */
      void enable_liu_storey_method() {
        this->m_beta_method = BetaMethod::LIU_STOREY;
      }

      /**
       * Enable the Dai-Yuan method.
       */
      void enable_dai_yuan_method() {
        this->m_beta_method = BetaMethod::DAI_YUAN;
      }

      /**
       * Enable the Hager-Zhang method.
       */
      void enable_hager_zhang_method() {
        this->m_beta_method = BetaMethod::HAGER_ZHANG;
      }

      /**
       * Enable the Hager-Zhang+ method.
       */
      void enable_hager_zhang_plus_method() {
        this->m_beta_method = BetaMethod::HAGER_ZHANG_PLUS;
      }

      /**
       * Set the conjugate gradient beta method.
       * \param[in] t_method Conjugate gradient beta method.
       */
      void set_method(const BetaMethod t_method) {
        this->m_beta_method = t_method;
      }

      /**
       * Get the enumeration type of the alpha computation method.
       * \return The alpha computation method.
       */
      AlphaMethod alpha_method() const {
        return this->m_alpha_method;
      }

      /**
       * Set the enumeration type of the alpha computation method.
       * \param[in] t_alpha_method Alpha computation method.
       */
      void alpha_method(const AlphaMethod t_alpha_method) {
        this->m_alpha_method = t_alpha_method;
      }

      /**
       * Enable the reuse of the last accepted line-search step.
       */
      void enable_accepted_step_alpha_method() {
        this->m_alpha_method = AlphaMethod::ACCEPTED_STEP;
      }

      /**
       * Enable the Barzilai-Borwein trial-step update.
       */
      void enable_barzilai_borwein_alpha_method() {
        this->m_alpha_method = AlphaMethod::BARZILAI_BORWEIN;
      }

      /**
       * Set the alpha computation method.
       * \param[in] t_alpha_method Alpha computation method.
       */
      void set_alpha_method(const AlphaMethod t_alpha_method) {
        this->m_alpha_method = t_alpha_method;
      }

      /**
       * Get the initial line-search step.
       * \return The initial line-search step.
       */
      Scalar initial_step() const {
        return this->m_initial_step;
      }

      /**
       * Set the initial line-search step.
       * \param[in] t_initial_step Initial line-search step.
       */
      void initial_step(const Scalar t_initial_step) {
        OPTIMIST_ASSERT(
            !std::isnan(t_initial_step) && std::isfinite(t_initial_step) &&
                t_initial_step > 0.0,
            "Optimist::Optimizer::ConjugateGradient::initial_step(...): "
            "invalid input detected.");
        this->m_initial_step = t_initial_step;
      }

      /**
       * Get the minimum admissible line-search step.
       * \return The minimum admissible line-search step.
       */
      Scalar min_step() const {
        return this->m_min_step;
      }

      /**
       * Set the minimum admissible line-search step.
       * \param[in] t_min_step Minimum admissible line-search step.
       */
      void min_step(const Scalar t_min_step) {
        OPTIMIST_ASSERT(
            !std::isnan(t_min_step) && std::isfinite(t_min_step) &&
                t_min_step > 0.0 && t_min_step <= this->m_max_step,
            "Optimist::Optimizer::ConjugateGradient::min_step(...): "
            "invalid input detected.");
        this->m_min_step = t_min_step;
      }

      /**
       * Get the maximum admissible line-search step.
       * \return The maximum admissible line-search step.
       */
      Scalar max_step() const {
        return this->m_max_step;
      }

      /**
       * Set the maximum admissible line-search step.
       * \param[in] t_max_step Maximum admissible line-search step.
       */
      void max_step(const Scalar t_max_step) {
        OPTIMIST_ASSERT(
            !std::isnan(t_max_step) && std::isfinite(t_max_step) &&
                t_max_step >= this->m_min_step,
            "Optimist::Optimizer::ConjugateGradient::max_step(...): "
            "invalid input detected.");
        this->m_max_step = t_max_step;
      }

      /**
       * Get the backtracking shrink factor.
       * \return The backtracking shrink factor.
       */
      Scalar step_shrink() const {
        return this->m_step_shrink;
      }

      /**
       * Get the line-search expansion factor.
       * \return The line-search expansion factor.
       */
      Scalar step_expand() const {
        return this->m_step_expand;
      }

      /**
       * Set the line-search expansion factor.
       * \param[in] t_step_expand Line-search expansion factor.
       */
      void step_expand(const Scalar t_step_expand) {
        OPTIMIST_ASSERT(
            !std::isnan(t_step_expand) && std::isfinite(t_step_expand) &&
                t_step_expand > 1.0,
            "Optimist::Optimizer::ConjugateGradient::step_expand(...): "
            "invalid input detected.");
        this->m_step_expand = t_step_expand;
      }

      /**
       * Set the backtracking shrink factor.
       * \param[in] t_step_shrink Backtracking shrink factor.
       */
      void step_shrink(const Scalar t_step_shrink) {
        OPTIMIST_ASSERT(
            !std::isnan(t_step_shrink) && std::isfinite(t_step_shrink) &&
                t_step_shrink > 0.0 && t_step_shrink < 1.0,
            "Optimist::Optimizer::ConjugateGradient::step_shrink(...): "
            "invalid input detected.");
        this->m_step_shrink = t_step_shrink;
      }

      /**
       * Get the Armijo sufficient decrease parameter.
       * \return The Armijo sufficient decrease parameter.
       */
      Scalar armijo_parameter() const {
        return this->m_armijo_parameter;
      }

      /**
       * Set the Armijo sufficient decrease parameter.
       * \param[in] t_armijo_parameter Armijo sufficient decrease parameter.
       */
      void armijo_parameter(const Scalar t_armijo_parameter) {
        OPTIMIST_ASSERT(
            !std::isnan(t_armijo_parameter) &&
                std::isfinite(t_armijo_parameter) && t_armijo_parameter > 0.0 &&
                t_armijo_parameter < 1.0,
            "Optimist::Optimizer::ConjugateGradient::armijo_parameter(...): "
            "invalid input detected.");
        this->m_armijo_parameter = t_armijo_parameter;
      }

      /**
       * Get the Wolfe curvature parameter.
       * \return The Wolfe curvature parameter.
       */
      Scalar curvature_parameter() const {
        return this->m_curvature_parameter;
      }

      /**
       * Set the Wolfe curvature parameter.
       * \param[in] t_curvature_parameter Wolfe curvature parameter.
       */
      void curvature_parameter(const Scalar t_curvature_parameter) {
        OPTIMIST_ASSERT(
            !std::isnan(t_curvature_parameter) &&
                std::isfinite(t_curvature_parameter) &&
                t_curvature_parameter > this->m_armijo_parameter &&
                t_curvature_parameter < 1.0,
            "Optimist::Optimizer::ConjugateGradient::curvature_parameter(...): "
            "invalid input detected.");
        this->m_curvature_parameter = t_curvature_parameter;
      }

      /**
       * Get the truncation parameter used by the Hager-Zhang+ method.
       * \return The truncation parameter.
       */
      Scalar truncation_parameter() const {
        return this->m_truncation_parameter;
      }

      /**
       * Set the truncation parameter used by the Hager-Zhang+ method.
       * \param[in] t_truncation_parameter Truncation parameter.
       */
      void truncation_parameter(const Scalar t_truncation_parameter) {
        OPTIMIST_ASSERT(!std::isnan(t_truncation_parameter) &&
                            std::isfinite(t_truncation_parameter) &&
                            t_truncation_parameter > 0.0,
                        "Optimist::Optimizer::ConjugateGradient::truncation_"
                        "parameter(...): "
                        "invalid input detected.");
        this->m_truncation_parameter = t_truncation_parameter;
      }

      /**
       * Get the sufficient descent parameter.
       * \return The sufficient descent parameter.
       */
      Scalar descent_parameter() const {
        return this->m_descent_parameter;
      }

      /**
       * Set the sufficient descent parameter.
       * \param[in] t_descent_parameter Sufficient descent parameter.
       */
      void descent_parameter(const Scalar t_descent_parameter) {
        OPTIMIST_ASSERT(
            !std::isnan(t_descent_parameter) &&
                std::isfinite(t_descent_parameter) &&
                t_descent_parameter > 0.0 && t_descent_parameter < 1.0,
            "Optimist::Optimizer::ConjugateGradient::descent_parameter(...): "
            "invalid input detected.");
        this->m_descent_parameter = t_descent_parameter;
      }

      /**
       * Get the maximum number of line-search iterations.
       * \return The maximum number of line-search iterations.
       */
      Integer max_line_search_iterations() const {
        return this->m_max_line_search_iterations;
      }

      /**
       * Set the maximum number of line-search iterations.
       * \param[in] t_max_line_search_iterations Maximum line-search iterations.
       */
      void max_line_search_iterations(
          const Integer t_max_line_search_iterations) {
        OPTIMIST_ASSERT(t_max_line_search_iterations > 0,
                        "Optimist::Optimizer::ConjugateGradient::max_line_"
                        "search_iterations(...): "
                        "invalid input detected.");
        this->m_max_line_search_iterations = t_max_line_search_iterations;
      }

      /**
       * Get the Powell restart flag.
       * \return The Powell restart flag.
       */
      bool powell_restart() const {
        return this->m_powell_restart;
      }

      /**
       * Set the Powell restart flag.
       * \param[in] t_powell_restart Powell restart flag.
       */
      void powell_restart(const bool t_powell_restart) {
        this->m_powell_restart = t_powell_restart;
      }

      /**
       * Enable the Powell restart safeguard.
       */
      void enable_powell_restart() {
        this->m_powell_restart = true;
      }

      /**
       * Disable the Powell restart safeguard.
       */
      void disable_powell_restart() {
        this->m_powell_restart = false;
      }

      /**
       * Get the Powell restart tolerance.
       * \return The Powell restart tolerance.
       */
      Scalar restart_tolerance() const {
        return this->m_restart_tolerance;
      }

      /**
       * Set the Powell restart tolerance.
       * \param[in] t_restart_tolerance Powell restart tolerance.
       */
      void restart_tolerance(const Scalar t_restart_tolerance) {
        OPTIMIST_ASSERT(
            !std::isnan(t_restart_tolerance) &&
                std::isfinite(t_restart_tolerance) &&
                t_restart_tolerance > 0.0 && t_restart_tolerance < 1.0,
            "Optimist::Optimizer::ConjugateGradient::restart_tolerance(...): "
            "invalid input detected.");
        this->m_restart_tolerance = t_restart_tolerance;
      }

      /**
       * Solve the unconstrained optimization problem
       * \f$ \min_{\mathbf{x}} f(\mathbf{x}) \f$, with
       * \f$ f: \mathbb{R}^n \rightarrow \mathbb{R} \f$.
       * \tparam FunctionLambda Function lambda type.
       * \tparam GradientLambda Gradient lambda type.
       * \param[in] function Function lambda.
       * \param[in] gradient Gradient lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool solve_impl(FunctionLambda &&function,
                      GradientLambda &&gradient,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::Optimizer::ConjugateGradient::solve_impl(...): "

        auto &function_ref{function};
        auto &gradient_ref{gradient};

        // Reset internal variables
        this->reset_counters();

        // Print header
        if (this->m_verbose) {
          this->header();
        }

        // Initialize variables
        bool success{false};
        Scalar function_old{0.0}, function_new{0.0},
            step_length{this->clamp_step(this->m_initial_step)}, beta{0.0},
            accepted_step{0.0}, gradient_norm{0.0};
        const Integer size{static_cast<Integer>(x_ini.size())};
        Vector x_old(x_ini), x_new(x_ini), y(size), step(size), direction_new;
        FirstDerivative gradient_old(size), gradient_new(size),
            direction_old(size);
        gradient_old.setZero();
        gradient_new.setZero();
        direction_old.setZero();
        direction_new.setZero();
        y.setZero();
        step.setZero();
        std::string notes{"Initialize"};

        // Set initial iteration
        success = this->evaluate_function(function_ref, x_old, function_old);
        OPTIMIST_ASSERT(success,
                        CMD "function evaluation failed at the initial point.");
        success = this->evaluate_gradient(gradient_ref, x_old, gradient_old);
        OPTIMIST_ASSERT(success,
                        CMD "gradient evaluation failed at the initial point.");
        direction_old = -gradient_old;

        // Algorithm iterations
        Scalar tolerance_gradient{this->m_tolerance};
        Scalar tolerance_step{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          gradient_norm = gradient_old.norm();
          if (this->m_verbose) {
            this->info(gradient_norm, notes);
          }
          if (gradient_norm < tolerance_gradient) {
            this->m_converged = true;
            break;
          }

          if (!this->has_sufficient_descent(gradient_old, direction_old)) {
            direction_old = -gradient_old;
            notes         = "Restart";
          }

          bool accepted{this->line_search(function_ref,
                                          gradient_ref,
                                          x_old,
                                          function_old,
                                          gradient_old,
                                          direction_old,
                                          step_length,
                                          x_new,
                                          function_new,
                                          gradient_new,
                                          accepted_step)};
          if (!accepted) {
            notes = "Line search failed";
            break;
          }

          step = x_new - x_old;
          if (step.norm() < tolerance_step ||
              gradient_new.norm() < tolerance_gradient) {
            x_old             = x_new;
            function_old      = function_new;
            gradient_old      = gradient_new;
            this->m_converged = true;
            break;
          }

          y = gradient_new - gradient_old;
          beta =
              this->compute_beta(gradient_old, gradient_new, direction_old, y);

          direction_new = -gradient_new + beta * direction_old;
          if (!std::isfinite(beta) ||
              this->should_restart(gradient_old, gradient_new) ||
              !this->has_sufficient_descent(gradient_new, direction_new)) {
            direction_new = -gradient_new;
            notes         = "Restart";
          }

          step_length   = this->next_trial_step(accepted_step,
                                              direction_old,
                                              gradient_old,
                                              gradient_new);
          x_old         = x_new;
          function_old  = function_new;
          gradient_old  = gradient_new;
          direction_old = direction_new;
        }

        // Print bottom
        if (this->m_verbose) {
          this->bottom();
        }

        // Convergence data
        x_sol = x_old;
        return this->m_converged;

#undef CMD
      }

    };  // class ConjugateGradient

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_CONJUGATE_GRADIENT_HH
