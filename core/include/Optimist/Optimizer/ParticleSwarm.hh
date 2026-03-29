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

#include "Optimist/Optimizer.hh"

namespace Optimist {
  namespace Optimizer {

    /**
     * \brief Inertia weight functor, which returns a constant weight.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    class ConstantWeight {
      Scalar m_weight{1.0}; /**< Constant inertia weight. */

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
      Scalar operator()(const Integer /*iteration*/,
                        const Integer /*max_iterations*/) const {
        return this->m_weight;
      }

    };  // class ConstantWeight

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - -
    // - - - - - - - -

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
       * Class class onstructor with default minimum and maximum weights equal
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
        return this->m_min_weight +
               (this->m_max_weight - this->m_min_weight) * factor;
      }

    };  // class LinearDecrease

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - -
    // - - - - - - - -

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
    // - -
    // - - - - - - - -

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
    // - -
    // - - - - - - - -

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
    // - -
    // - - - - - - - -

    /*\
     |   ____            _   _      _      ____
     |  |  _ \ __ _ _ __| |_(_) ___| | ___/ ___|_      ____ _ _ __ _ __ ___
     |  | |_) / _` | '__| __| |/ __| |/ _ \___ \ \ /\ / / _` | '__| '_ ` _ \
     |  |  __/ (_| | |  | |_| | (__| |  __/___) \ V  V / (_| | |  | | | | | |
     |  |_|   \__,_|_|   \__|_|\___|_|\___|____/ \_/\_/ \__,_|_|  |_| |_| |_|
     |
    \*/

    /**
     * \brief Class container for the particle swarm optimizer.
     *
     * The optimization process can be configured by providing an inertia weight
     * strategy functor. The inertia weight functor determines the amount of
     * velocity, which is maintained from the previous iterations. It has a huge
     * effect on convergence speed and stability of the optimization.
     *
     * \tparam Vector Eigen vector type.
     * \tparam InertiaWeightStrategy Inertia weight strategy functor type.
     */
    template <typename Vector,
              typename InertiaWeightStrategy =
                  ConstantWeight<typename Vector::Scalar>>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class ParticleSwarm
        : public Optimizer<Vector, ParticleSwarm<Vector>, true> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;

      using Particles = std::vector<Vector>; /**< Particles container type. */

     private:
      InertiaWeightStrategy
          m_weight_strategy;        /**< Inertia weight strategy functor. */

      Scalar m_phi_p{SQRT_EPSILON}; /**< Cognitive coefficient. */
      Scalar m_phi_g{SQRT_EPSILON}; /**<
      Scalar m_max_velocity;        /**< Maximum velocity for each particle. */

      Scalar xeps_;
      Scalar feps_;
      std::function<Scalar()> this->m_dice;

     public:
      /**
       * Class constructor for the ParticleSwarm solver.
       */
      ParticleSwarm() {}

      /**
       * Get the ParticleSwarm's solver name.
       * \return The ParticleSwarm's solver name.
       */
      constexpr std::string name_impl() const {
        return "ParticleSwarm";
      }

      /**
       * Randomize the particles within the given bounds.
       * \param[out] particles Matrix containing the particles to be randomized.
       */
      void randomize_particles(Particles &particles) const {
        for (Integer i{0}; i < particles.size(); ++i) {
          for (Integer j{0}; j < particles.rows(); ++j) {
            Scalar lb{this->m_lower_bound(j)};
            Scalar ub{this->m_upper_bound(j)};
            particles(j, i) = lb + (this->m_dice() * (ub - lb));
          }
        }
      }

      /**
       * Randomize the velocities within the given bounds.
       * \param[out] velocities Matrix containing the velocities to be
       * randomized.
       */
      void randomize_velocities(Particles &velocities) const {
        for (Integer i{0}; i < velocities.cols(); ++i) {
          for (Integer j{0}; j < velocities.rows(); ++j) {
            Scalar lb{this->m_lower_bound(j)};
            Scalar ub{this->m_upper_bound(j)};
            Scalar diff{ub - lb};
            Scalar vel{-diff + (this->m_dice() * 2.0 * diff)};
            velocities(j, i) =
                std::min(m_max_velocity, std::max(-m_max_velocity, vel));
          }
        }
      }

      void evaluate_objective(const Particles &particles, Vector &fvals) {
        for (Integer i{0}; i < particles.size(); ++i) {
          fvals(i) = objective_(particles.col(i));
        }
      }

      void maintainBounds(Matrix &particles) const {
        for (Integer i{0}; i < particles.size(); ++i) {
          for (Integer j{0}; j < particles.rows(); ++j) {
            Scalar minval{bounds(0, j)};
            Scalar maxval{bounds(1, j)};
            particles[j](i) =
                std::min(maxval, std::max(minval, particles[j](i)));
          }
        }
      }

      void calculateVelocities(const Matrix &particles,
                               const Matrix &best_particles,
                               const Integer gbest,
                               const Integer iteration,
                               Matrix &velocities) {
        assert(velocities.rows() == particles.rows());
        assert(velocities.cols() == particles.size());
        assert(velocities.rows() == best_particles.rows());
        assert(velocities.cols() == best_particles.size());
        assert(gbest < best_particles.size());

        Scalar weight = m_weight_strategy(iteration, this->m_max_iterations);

        for (Integer i{0}; i < velocities.cols(); ++i) {
          for (Integer j{0}; j < velocities.rows(); ++j) {
            Scalar velp =
                this->m_dice() * (best_particles(j, i) - particles(j, i));
            Scalar velg =
                this->m_dice() * (best_particles(j, gbest) - particles(j, i));
            Scalar vel =
                weight * velocities(j, i) + m_phi_p * velp + m_phi_g * velg;

            if (m_max_velocity > 0)
              vel = std::min(m_max_velocity, std::max(-m_max_velocity, vel));

            velocities(j, i) = vel;
          }
        }
      }

      Result _minimize(const Matrix &bounds, Matrix &particles) {
        Matrix velocities(particles.rows(), particles.size());

        Vector fvals(particles.size());

        Matrix best_particles = particles;
        Vector bestFvals(particles.size());

        Matrix prevParticles(particles.rows(), particles.size());
        Vector prevFvals(particles.size());

        Vector diff(particles.rows());

        Index gbest = 0;

        // initialize velocities randomly
        randomize_velocities(bounds, velocities);

        // evaluate objective function for the initial particles
        evaluate_objective(particles, fvals);
        bestFvals = fvals;
        bestFvals.minCoeff(&gbest);

        // init stop conditions
        Index iterations = 0;
        Scalar fchange   = feps_ + 1;
        Scalar xchange   = xeps_ + 1;

        while ((this->m_iterations < this->m_max_iterations) &&
               fchange > feps_ && xchange > xeps_) {
          // calculate new velocities
          calculateVelocities(particles,
                              best_particles,
                              gbest,
                              iterations,
                              velocities);

          // move particles by velocity and stay within bounds
          particles += velocities;
          maintainBounds(bounds, particles);

          // evaluate objective for moved particles
          evaluate_objective(particles, fvals);

          prevParticles = best_particles;
          prevFvals     = bestFvals;

          for (Integer i{0}; i < fvals.size(); ++i) {
            // check if there was an improvement and update best vals
            if (fvals(i) < bestFvals(i)) {
              bestFvals(i)          = fvals(i);
              best_particles.col(i) = particles.col(i);
            }
          }
          bestFvals.minCoeff(&gbest);

          // calculate new diffs
          xchange = (best_particles - prevParticles).colwise().norm().sum();
          fchange = (bestFvals - prevFvals).array().abs().sum();

          xchange /= best_particles.size();
          fchange /= bestFvals.size();

          if (this->m_verbose) {
            this->info(gbest);
          }

          ++iterations;
        }

        result.iterations = iterations;
        result.converged  = fchange <= feps_ || xchange <= xeps_;
        result.fval       = bestFvals(gbest);
        result.xval       = best_particles.col(gbest);

        return result;
      }

     public:
      ParticleSwarmOptimization()
          : m_weight_strategy(),
            xeps_(static_cast<Scalar>(1e-6)),
            feps_(static_cast<Scalar>(1e-6)),
            m_phi_p(static_cast<Scalar>(2.0)),
            m_phi_g(static_cast<Scalar>(2.0)),
            m_max_velocity(static_cast<Scalar>(0.0)),
            this->m_verbose(0),
            this->m_dice() {
        std::default_random_engine gen(std::time(0));
        std::uniform_real_distribution<Scalar> distrib(0.0, 1.0);
        this->m_dice = std::bind(distrib, gen);
      }

      /**
       * Set the minimum average change of particles per iteration.
       * If the average change of particles (input parameters) falls below
       * this value, the optimization terminates.
       * \param change minimum change of input paramaters
       */
      void setMinParticleChange(
          const Scalar change = ParticleSwarm::SQRT_EPSILON) {
        xeps_ = change;
      }

      /** Set the minimum average change of function values per iteration.
       * If the average change of functions values falls below
       * this value, the optimization terminates.
       * \param change minimum change of function values */
      void setMinFunctionChange(
          const Scalar change = ParticleSwarm::SQRT_EPSILON) {
        feps_ = change;
      }

      /** Set the tendency of particles to move towards their local optimum
       * found so far.
       * Each particle individually maintains a memory of where it has
       * visited the lowest function value so far.
       * Increasing this value increases the particles' tendency to move
       * towards that point.
       * \param phip tendency to move towards individual optimum */
      void setPhiParticles(const Scalar phip) {
        m_phi_p = phip;
      }

      /** Set the tendency of particles to move towards the global optimum
       * found so far.
       * The swarm maintains a collective memory of where it has visited the
       * lowest function value so far.
       * Increasing this value increases the particles' tendency to move
       * towards that point.
       * \param phig tendency to move towards collective optimum */
      void setPhiGlobal(const Scalar phig) {
        m_phi_g = phig;
      }

      /**
       * Set an upper bound for the velocity of particles.
       * A particle cannot move faster than this value, which may prevent
       * divergence.
       * \param[in] max_velocity maximum velocity of a particle */
      void setMaxVelocity(const Scalar max_velocity) {
        this->m_max_velocity = max_velocity;
      }

      void inertia_weight_strategy(
          const InertiaWeightStrategy &t_weightStrategy) {
        this->m_weight_strategy = t_weightStrategy;
      }

      /** Perform minimization with the given bounds and number of particels.
       *
       * The swarm of particles will be drawn uniform randomly within the
       * given bounds.
       *
       * The bounds matrix has to have 2 rows and one column per dimension
       * of particle. The first row holds the minimum value of the respective
       * dimension and the second row holds the maximum value.
       *
       * \param bounds 2xM matrix for bounds of M-dimensional particles
       * \param cnt number of particles used for optimization */
      Result minimize(const Integer cnt) {
        if (cnt == 0)
          throw std::runtime_error("particle count cannot be 0");
        if (bounds.rows() != 2)
          throw std::runtime_error("bounds has not exactly 2 rows (min, max)");
        for (Integer i{0}; i < bounds.cols(); ++i) {
          if (bounds(0, i) >= bounds(1, i))
            throw std::runtime_error("bounds min is greater than max");
        }

        Matrix particles(bounds.cols(), cnt);
        randomize_particles(bounds, particles);

        return _minimize(bounds, particles);
      }

      /** Perform minimization with the given bounds, number of particels and
       * initial guess.
       *
       * The swarm of particles will be drawn uniform randomly within the
       * given bounds.
       *
       * The bounds matrix has to have 2 rows and one column per dimension
       * of particle. The first row holds the minimum value of the respective
       * dimension and the second row holds the maximum value.
       *
       * The initial guess vector has to have the same length as the number
       * of columns of the bounds. It will be included as one particle of
       * the swarm.
       *
       * \param bounds 2xM matrix for bounds of M-dimensional particles
       * \param cnt number of particles used for optimization
       * \param initGuess initial guess for a particle */
      Result minimize(const Integer cnt, const Vector &initGuess) {
        if (cnt == 0)
          throw std::runtime_error("particle count cannot be 0");
        if (bounds.rows() != 2)
          throw std::runtime_error("bounds has not exactly 2 rows (min, max)");
        for (Integer i{0}; i < bounds.cols(); ++i) {
          if (bounds(0, i) >= bounds(1, i))
            throw std::runtime_error("bounds min is greater than max");
        }
        if (bounds.cols() != initGuess.size())
          throw std::runtime_error(
              "init guess and bounds have different dimensions");

        Matrix particles(bounds.cols(), cnt);
        randomize_particles(bounds, particles);
        particles.col(0) = initGuess;
        maintainBounds(bounds, particles);

        return _minimize(bounds, particles);
      }

      /** Perform minimization with the given bounds and a pre-computed
       * swarm of particles.
       *
       * The bounds matrix has to have 2 rows and one column per dimension
       * of particle. The first row holds the minimum value of the respective
       * dimension and the second row holds the maximum value.
       *
       * \param bounds 2xM matrix for bounds of M-dimensional particles
       * \param particles initial swarm used for optimization */
      Result minimize(Matrix &particles) {
        if (bounds.rows() != 2)
          throw std::runtime_error("bounds has not exactly 2 rows (min, max)");
        if (bounds.cols() != particles.rows())
          throw std::runtime_error(
              "columns of bounds and rows of "
              "particles do not match");
        for (Integer i{0}; i < bounds.cols(); ++i) {
          if (bounds(0, i) >= bounds(1, i))
            throw std::runtime_error("bounds min is greater than max");
        }

        maintainBounds(bounds, particles);

        return _minimize(bounds, particles);
      }

      void getRandomParticles(const Integer cnt, Matrix &particles) {
        particles.resize(bounds.cols(), cnt);
        randomize_particles(bounds, particles);
      }

    };  // class ParticleSwarmOptimization

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_PARTICLE_SWARM_HH
