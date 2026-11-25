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

#ifndef OPTIMIST_OPTIMIZER_PARTICLE_SWARM_HH
#define OPTIMIST_OPTIMIZER_PARTICLE_SWARM_HH

#include <random>

#include "Optimist/Optimizer.hh"

namespace Optimist {
  namespace Optimizer {

    /**
     * \brief Inertia weight functor, which returns a constant weight.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    class ConstantWeight {
      Scalar m_weight{0.7298}; /**< Constant inertia weight (Clerc). */

     public:
      /**
       * Class constructor with default constant inertia weight equal to 1.0.
       */
      ConstantWeight() {}

      /**
       * Class constructor that accepts the constant inertia weight.
       * \param[in] weight Constant inertia weight.
       */
      ConstantWeight(const Scalar weight) : m_weight(weight) {}

      /**
       * Get the inertia weight.
       * \param[in] iteration Current iteration.
       * \param[in] max_iterations Maximum number of iterations.
       * \return The inertia weight.
       */
      Scalar operator()([[maybe_unused]] const Integer iteration,
                        [[maybe_unused]] const Integer max_iterations) const {
        return this->m_weight;
      }

    };  // class ConstantWeight

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * \brief Inertia weight functor, which decreases linearly over time.
     *
     * \tparam Scalar Floating-point number type.
     * The inertia weight is calculated by the following formula:
     * \f[
     *   w = w_{\min} + (w_{\max} - w_{\min}) * \displaystyle\frac{t}{t_{\max}}
     * \text{.}
     * \f]
     */
    template <typename Scalar>
    class LinearDecrease {
      Scalar m_min_weight{0.4}; /**< Minimum inertia weight. */
      Scalar m_max_weight{0.9}; /**< Maximum inertia weight. */

     public:
      /**
       * Class constructor with default minimum and maximum weights equal
       * to 0.4 and 0.9, respectively.
       */
      LinearDecrease() {}

      /**
       * Class constructor that accepts the minimum and maximum weight of the
       * linear decrease.
       * \param[in] min_weight Lower bound of the inertia weight.
       * \param[in] max_weight Upper bound of the inertia weight.
       */
      LinearDecrease(const Scalar min_weight, const Scalar max_weight)
          : m_min_weight(min_weight), m_max_weight(max_weight) {}

      /**
       * Get the inertia weight.
       * \param[in] iteration Current iteration.
       * \param[in] max_iterations Maximum number of iterations.
       * \return The inertia weight.
       */
      Scalar operator()(const Integer iteration,
                        const Integer max_iterations) const {
        Scalar factor{static_cast<Scalar>(iteration) /
                      static_cast<Scalar>(max_iterations)};
        return this->m_max_weight -
               (this->m_max_weight - this->m_min_weight) * factor;
      }

    };  // class LinearDecrease

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * \brief Inertia weight functor, which decreases exponentially over time.
     *
     * \tparam Scalar Floating-point number type.
     * The inertia weight is calculated by the following formula:
     * \f[
     *   w = wMin + (wMax - wMin) * \exp\left(-\frac{t}{\frac{tMax}{10}}\right)
     * \text{.}
     * \f]
     */
    template <typename Scalar>
    class ExponentialDecrease1 {
      Scalar m_min_weight{0.4}; /**< Minimum inertia weight. */
      Scalar m_max_weight{0.9}; /**< Maximum inertia weight. */

     public:
      /**
       * Class constructor with minimum and maximum weights equal to 0.4 and
       * 0.9, respectively.
       */
      ExponentialDecrease1() {}

      /**
       * Class constructor that accepts the minimum and maximum weight of the
       * exponential decrease.
       * \param[in] min_weight Lower bound of the inertia weight
       * \param[in] max_weight Upper bound of the inertia weight */
      ExponentialDecrease1(const Scalar min_weight, const Scalar max_weight)
          : m_min_weight(min_weight), m_max_weight(max_weight) {}

      /**
       * Get the inertia weight.
       * \param[in] iteration Current iteration.
       * \param[in] max_iterations Maximum number of iterations.
       * \return The inertia weight.
       */
      Scalar operator()(const Integer iteration,
                        const Integer max_iterations) const {
        Scalar exponent{static_cast<Scalar>(iteration) /
                        (static_cast<Scalar>(max_iterations) / 10.0)};
        return this->m_min_weight +
               (this->m_max_weight - this->m_min_weight) * std::exp(-exponent);
      }

    };  // class ExponentialDecrease1

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /** \brief Inertia weight functor, which decreases exponentially over time.
     *
     * \tparam Scalar Floating-point number type.
     * The inertia weight is calculated by the following formula:
     * \f[
     *   w = w_{\min} + (w_{\max} - w_{\min}) *
     * \exp\left(-\left(\frac{t}{\frac{t_{\max}}{4}}\right)^2
     *   \right) \text{.}
     * \f]
     */
    template <typename Scalar>
    class ExponentialDecrease2 {
      Scalar m_min_weight{0.4}; /**< Minimum inertia weight */
      Scalar m_max_weight{0.9}; /**< Maximum inertia weight */

     public:
      /**
       * Class constructor with minimum and maximum weights equal to 0.4 and
       * 0.9, respectively.
       */
      ExponentialDecrease2() {}

      /**
       * Class constructor that accepts the minimum and maximum weight of the
       * exponential decrease.
       * \param[in] min_weight Lower bound of the inertia weight
       * \param[in] max_weight Upper bound of the inertia weight
       */
      ExponentialDecrease2(const Scalar min_weight, const Scalar max_weight)
          : m_min_weight(min_weight), m_max_weight(max_weight) {}

      /**
       * Get the inertia weight.
       * \param[in] iteration Current iteration.
       * \param[in] max_iterations Maximum number of iterations.
       * \return The inertia weight.
       */
      Scalar operator()(const Integer iteration,
                        const Integer max_iterations) const {
        Scalar exponent{static_cast<Scalar>(iteration) /
                        (static_cast<Scalar>(max_iterations) / 4.0)};
        exponent *= exponent;
        return this->m_min_weight +
               (this->m_max_weight - this->m_min_weight) * std::exp(-exponent);
      }

    };  // class ExponentialDecrease2

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /** \brief Inertia weight functor, which decreases exponentially over time.
     *
     * \tparam Scalar Floating-point number type.
     * The inertia weight is calculated by the following formula:
     * \f[
     *   w = (w_{\max} - w_{\min} - d1) * \exp\left(\displaystyle\frac{1}{1 + d2
     * \frac{t}{t_{\max}}}
     *   \right) \text{.}
     * \f]
     */
    template <typename Scalar>
    class ExponentialDecrease3 {
      Scalar m_min_weight{0.4};  /**< Minimum inertia weight */
      Scalar m_max_weight{0.95}; /**< Maximum inertia weight */
      Scalar m_d_1{0.2};         /**< First control factor */
      Scalar m_d_2{7.0};         /**< Second control factor */

     public:
      /**
       * Class constructor with minimum and maximum weights equal to 0.4 and
       * 0.95, respectively, and control factors equal to 0.2 and 7.0.
       */
      ExponentialDecrease3() {}

      /**
       * Class constructor that accepts the minimum and maximum weight of the
       * exponential decrease and the two control factors.
       * \param[in] min_weight Lower bound of the inertia weight
       * \param[in] max_weight Upper bound of the inertia weight
       * \param[in] d_1 First control factor
       * \param[in] d_2 Second control factor
       */
      ExponentialDecrease3(const Scalar min_weight,
                           const Scalar max_weight,
                           const Scalar d_1,
                           const Scalar d_2)
          : m_min_weight(min_weight),
            m_max_weight(max_weight),
            m_d_1(d_1),
            m_d_2(d_2) {}

      /**
       * Get the inertia weight.
       * \param[in] iteration Current iteration.
       * \param[in] max_iterations Maximum number of iterations.
       * \return The inertia weight.
       */
      Scalar operator()(const Integer iteration,
                        const Integer max_iterations) const {
        const Scalar factor{static_cast<Scalar>(iteration) /
                            static_cast<Scalar>(max_iterations)};
        const Scalar exponent{1.0 / (1.0 + this->m_d_2 * factor)};
        return (this->m_max_weight - this->m_min_weight - this->m_d_1) *
               std::exp(exponent);
      }

    };  // class ExponentialDecrease3

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * \brief Class container for the particle swarm optimizer.
     *
     * \includedoc docs/markdown/Optimizer/ParticleSwarm.md
     *
     * \tparam Vector Eigen vector type.
     * \tparam InertiaWeightStrategy Inertia weight strategy functor.
     */
    template <typename Vector,
              typename InertiaWeightStrategy =
                  LinearDecrease<typename Vector::Scalar>>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class ParticleSwarm
        : public Optimizer<Vector,
                           ParticleSwarm<Vector, InertiaWeightStrategy>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      InertiaWeightStrategy m_weight_strategy; /**< Inertia weight strategy. */
      Integer m_particles{0}; /**< Swarm size (0 means automatic). */
      Scalar m_phi_p{1.49445}; /**< Cognitive coefficient. */
      Scalar m_phi_g{1.49445}; /**< Social coefficient. */
      Scalar m_max_velocity{0.2}; /**< Maximum velocity (fraction of box). */
      Scalar m_radius{1.0}; /**< Search radius around the initial guess when
                               the bounds are not finite. */
      std::mt19937 m_generator{42}; /**< Random number generator. */

     public:
      /**
       * Class constructor for the particle swarm solver.
       */
      ParticleSwarm() {}

      /**
       * Get the particle swarm solver name.
       */
      constexpr std::string name_impl() const {
        return "ParticleSwarm";
      }

      /**
       * Set the seed of the random number generator.
       */
      void seed(const unsigned t_seed) {
        this->m_generator.seed(t_seed);
      }

      /**
       * Get the swarm size (0 means \f$ 10 + 2\sqrt{n} \f$).
       */
      Integer particles() const {
        return this->m_particles;
      }

      /**
       * Set the swarm size (0 means \f$ 10 + 2\sqrt{n} \f$).
       */
      void particles(const Integer t_particles) {
        OPTIMIST_ASSERT(t_particles >= 0,
                        "Optimist::Optimizer::ParticleSwarm::particles(...): "
                        "invalid input detected.");
        this->m_particles = t_particles;
      }

      /**
       * Set the cognitive and social coefficients.
       */
      void coefficients(const Scalar t_phi_p, const Scalar t_phi_g) {
        OPTIMIST_ASSERT(t_phi_p >= 0 && t_phi_g >= 0,
                        "Optimist::Optimizer::ParticleSwarm::coefficients(...)"
                        ": invalid input detected.");
        this->m_phi_p = t_phi_p;
        this->m_phi_g = t_phi_g;
      }

      /**
       * Set the maximum velocity as a fraction of the search box width.
       */
      void max_velocity(const Scalar t_max_velocity) {
        OPTIMIST_ASSERT(t_max_velocity > 0,
                        "Optimist::Optimizer::ParticleSwarm::max_velocity(...)"
                        ": invalid input detected.");
        this->m_max_velocity = t_max_velocity;
      }

      /**
       * Set the search radius around the initial guess, used when the bounds
       * are not finite.
       */
      void radius(const Scalar t_radius) {
        OPTIMIST_ASSERT(t_radius > 0,
                        "Optimist::Optimizer::ParticleSwarm::radius(...): "
                        "invalid input detected.");
        this->m_radius = t_radius;
      }

      /**
       * Set the inertia weight strategy.
       */
      void inertia_weight_strategy(const InertiaWeightStrategy &t_strategy) {
        this->m_weight_strategy = t_strategy;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] x_ini Initialization point (first particle).
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda &&function,
                      const Vector &x_ini,
                      Vector &x_sol) {
        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        const Integer n{static_cast<Integer>(x_ini.size())};
        const Integer m{this->m_particles > 0
                            ? this->m_particles
                            : 10 + static_cast<Integer>(
                                       2 * std::sqrt(static_cast<Scalar>(n)))};
        const Dense x0(Utils::to_dense(x_ini));

        // Search box
        Dense lb(n), ub(n);
        if (this->has_finite_bounds()) {
          lb = Utils::to_dense(this->m_lower_bound);
          ub = Utils::to_dense(this->m_upper_bound);
        } else {
          lb = x0.array() - this->m_radius;
          ub = x0.array() + this->m_radius;
        }
        const Dense width{ub - lb};
        const Dense v_max{this->m_max_velocity * width};

        // Swarm initialization (the first particle is the initial guess)
        std::uniform_real_distribution<Scalar> unif(0, 1);
        DenseMatrix x(n, m), v(n, m), p(n, m);
        Dense f(m), f_p(m);
        for (Integer i{0}; i < m; ++i) {
          for (Integer j{0}; j < n; ++j) {
            x(j, i) = lb(j) + unif(this->m_generator) * width(j);
            v(j, i) = (2 * unif(this->m_generator) - 1) * v_max(j);
          }
        }
        x.col(0) = x0.cwiseMax(lb).cwiseMin(ub);
        for (Integer i{0}; i < m; ++i) {
          f(i) = this->eval(function, x.col(i));
        }
        p   = x;
        f_p = f;
        Integer best;
        f_p.minCoeff(&best);

        std::string notes{"Initialize"};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          // Convergence: personal bests collapsed onto the global best
          Scalar radius{0};
          for (Integer i{0}; i < m; ++i) {
            radius = std::max(radius, (p.col(i) - p.col(best)).norm());
          }
          if (this->m_verbose) {
            this->info(radius, notes);
          }
          // Converged when the personal bests collapsed onto the global best,
          // or when they are all equal to rounding level in value (the swarm
          // cannot resolve the minimizer any further)
          const Scalar scale{1 + p.col(best).norm()};
          const Scalar spread{f_p.maxCoeff() - f_p(best)};
          if (radius < this->m_tolerance * scale ||
              (spread <= 100 * EPSILON * std::max(Scalar(1), std::abs(f_p(best))) &&
               radius < std::sqrt(this->m_tolerance) * scale)) {
            this->m_converged = true;
            break;
          }
          if (this->m_function_evaluations + m >
              this->m_max_function_evaluations) {
            notes = "Budget exhausted";
            break;
          }

          // Velocity and position update
          const Scalar w{this->m_weight_strategy(this->m_iterations,
                                                 this->m_max_iterations)};
          for (Integer i{0}; i < m; ++i) {
            for (Integer j{0}; j < n; ++j) {
              Scalar vel{w * v(j, i) +
                         this->m_phi_p * unif(this->m_generator) *
                             (p(j, i) - x(j, i)) +
                         this->m_phi_g * unif(this->m_generator) *
                             (p(j, best) - x(j, i))};
              vel     = std::clamp(vel, -v_max(j), v_max(j));
              v(j, i) = vel;
              x(j, i) = x(j, i) + vel;
              // Reflect on the box walls
              if (x(j, i) < lb(j)) {
                x(j, i) = lb(j);
                v(j, i) = -v(j, i) / 2;
              } else if (x(j, i) > ub(j)) {
                x(j, i) = ub(j);
                v(j, i) = -v(j, i) / 2;
              }
            }
            f(i) = this->eval(function, x.col(i));
            if (f(i) < f_p(i)) {
              f_p(i)   = f(i);
              p.col(i) = x.col(i);
            }
          }
          const Integer best_old{best};
          f_p.minCoeff(&best);
          notes = best != best_old ? "New best" : "-";
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = Utils::from_dense<Vector>(p.col(best));
        return this->m_converged;
      }

     private:
      template <typename FunctionLambda, typename Derived>
      Scalar eval(FunctionLambda &function,
                  const Eigen::MatrixBase<Derived> &x) {
        Scalar f;
        if (!this->evaluate_function(function, Utils::from_dense<Vector>(x), f) ||
            !std::isfinite(f)) {
          return INFTY;
        }
        return f;
      }

    };  // class ParticleSwarm

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_PARTICLE_SWARM_HH
