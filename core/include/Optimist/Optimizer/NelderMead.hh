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

#ifndef OPTIMIST_OPTIMIZER_NELDER_MEAD_HH
#define OPTIMIST_OPTIMIZER_NELDER_MEAD_HH

#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the Nelder-Mead's method.
     *
     * \includedoc docs/markdown/Optimizer/NelderMead.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class NelderMead : public Optimizer<Vector, NelderMead<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Simplex moves (used for verbose output).
       */
      using Move = enum class Move : Integer {
        INITIALIZE = 0,
        REFLECT    = 1,
        EXPAND_E   = 2,
        EXPAND_R   = 3,
        CONTRACT_O = 4,
        CONTRACT_I = 5,
        SHRINK     = 6,
        RESTART    = 7
      };

      /**
       * Initial simplex shape.
       */
      using Method = enum class Method : Integer {
        STANDARD = 0, /**< Axis-aligned simplex (best of \f$ \pm\delta \f$). */
        SPENDLEY = 1  /**< Regular simplex of Spendley, Hext and Himsworth. */
      };

     private:
      Method m_method{Method::STANDARD}; /**< Initial simplex shape. */
      Scalar m_delta{1.0};               /**< Initial simplex size. */
      Scalar m_rho{1.0};                 /**< Reflection coefficient. */
      Scalar m_chi{2.0};                 /**< Expansion coefficient. */
      Scalar m_gamma{0.5};               /**< Contraction coefficient. */
      Scalar m_sigma{0.5};               /**< Shrink coefficient. */
      bool m_adaptive{false}; /**< Dimension-dependent coefficients flag. */
      Scalar m_volume_tolerance{
        1.0e-4}; /**< Normalized volume below which the simplex restarts. */

      // Internal variables
      Integer m_dim{0};       /**< Problem dimension. */
      DenseMatrix m_p;        /**< Simplex points (one per column). */
      Dense m_f;              /**< Function values at simplex points. */
      std::vector<Integer> m_order; /**< Simplex points sorted by value. */
      Scalar m_diameter{0};   /**< Simplex diameter. */

     public:
      /**
       * Class constructor for the Nelder-Mead solver.
       */
      NelderMead() {}

      /**
       * Class constructor for the Nelder-Mead solver.
       * \param[in] t_delta Initial simplex size.
       */
      explicit NelderMead(const Scalar t_delta) : NelderMead() {
        this->delta(t_delta);
      }

      /**
       * Get the Nelder-Mead solver name.
       * \return The Nelder-Mead solver name.
       */
      constexpr std::string name_impl() const {
        return "NelderMead";
      }

      /**
       * Get the initial simplex shape.
       */
      Method method() const {
        return this->m_method;
      }

      /**
       * Set the initial simplex shape.
       */
      void method(const Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Get the initial simplex size \f$ \delta \f$.
       */
      Scalar delta() const {
        return this->m_delta;
      }

      /**
       * Set the initial simplex size \f$ \delta \f$.
       */
      void delta(const Scalar t_delta) {
        OPTIMIST_ASSERT(t_delta > 0,
                        "Optimist::Optimizer::NelderMead::delta(...): invalid "
                        "input detected.");
        this->m_delta = t_delta;
      }

      /**
       * Get the reflection coefficient \f$ \rho \f$.
       */
      Scalar rho() const {
        return this->m_rho;
      }

      /**
       * Set the reflection coefficient \f$ \rho > 0 \f$.
       */
      void rho(const Scalar t_rho) {
        OPTIMIST_ASSERT(t_rho > 0,
                        "Optimist::Optimizer::NelderMead::rho(...): invalid "
                        "input detected.");
        this->m_rho = t_rho;
      }

      /**
       * Get the expansion coefficient \f$ \chi \f$.
       */
      Scalar chi() const {
        return this->m_chi;
      }

      /**
       * Set the expansion coefficient \f$ \chi > 1 \f$.
       */
      void chi(const Scalar t_chi) {
        OPTIMIST_ASSERT(t_chi > 1,
                        "Optimist::Optimizer::NelderMead::chi(...): invalid "
                        "input detected.");
        this->m_chi = t_chi;
      }

      /**
       * Get the contraction coefficient \f$ \gamma \f$.
       */
      Scalar gamma() const {
        return this->m_gamma;
      }

      /**
       * Set the contraction coefficient \f$ 0 < \gamma < 1 \f$.
       */
      void gamma(const Scalar t_gamma) {
        OPTIMIST_ASSERT(t_gamma > 0 && t_gamma < 1,
                        "Optimist::Optimizer::NelderMead::gamma(...): invalid "
                        "input detected.");
        this->m_gamma = t_gamma;
      }

      /**
       * Get the shrink coefficient \f$ \sigma \f$.
       */
      Scalar sigma() const {
        return this->m_sigma;
      }

      /**
       * Set the shrink coefficient \f$ 0 < \sigma < 1 \f$.
       */
      void sigma(const Scalar t_sigma) {
        OPTIMIST_ASSERT(t_sigma > 0 && t_sigma < 1,
                        "Optimist::Optimizer::NelderMead::sigma(...): invalid "
                        "input detected.");
        this->m_sigma = t_sigma;
      }

      /**
       * Enable or disable the dimension-dependent coefficients of Gao and Han
       * (2012), which improve the behaviour in high dimension. When enabled,
       * the user-defined coefficients are ignored.
       */
      void adaptive_mode(const bool t_adaptive) {
        this->m_adaptive = t_adaptive;
      }

      /**
       * Get the adaptive coefficients flag.
       */
      bool adaptive_mode() const {
        return this->m_adaptive;
      }

      /**
       * Get the normalized volume tolerance that triggers a restart.
       */
      Scalar volume_tolerance() const {
        return this->m_volume_tolerance;
      }

      /**
       * Set the normalized volume tolerance that triggers a restart.
       */
      void volume_tolerance(const Scalar t_volume_tolerance) {
        OPTIMIST_ASSERT(
            t_volume_tolerance > 0 && t_volume_tolerance < 1,
            "Optimist::Optimizer::NelderMead::volume_tolerance(...): "
            "invalid input detected.");
        this->m_volume_tolerance = t_volume_tolerance;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda &&function,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::Optimizer::NelderMead::solve_impl(...): "

        // Reset internal variables
        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        // Coefficients
        Scalar rho{this->m_rho}, chi{this->m_chi}, gamma{this->m_gamma},
            sigma{this->m_sigma};
        this->m_dim = static_cast<Integer>(x_ini.size());
        if (this->m_adaptive) {
          const Scalar n{static_cast<Scalar>(this->m_dim)};
          rho   = 1.0;
          chi   = 1.0 + 2.0 / n;
          gamma = 0.75 - 1.0 / (2.0 * n);
          sigma = this->m_dim > 1 ? 1.0 - 1.0 / n : 0.5;
        }

        // Initialize simplex
        auto eval = [this, &function](const Dense &x) -> Scalar {
          Scalar f;
          bool success{this->evaluate_function(function,
                                               Utils::from_dense<Vector>(x),
                                               f)};
          OPTIMIST_ASSERT(success, CMD "function evaluation failed.");
          return std::isfinite(f) ? f : INFTY;
        };
        this->initialize_simplex(eval, Utils::to_dense(x_ini), this->m_delta);

        // Algorithm iterations
        const Scalar tol{this->m_tolerance};
        Move move{Move::INITIALIZE};
        Dense centroid(this->m_dim), x_r(this->m_dim), x_e(this->m_dim),
            x_c(this->m_dim);
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          // Sort the simplex points: best (low), second worst (mid) and worst
          // (high)
          this->sort();
          const Integer low{this->m_order.front()};
          const Integer high{this->m_order.back()};
          const Integer mid{this->m_order[this->m_dim - 1]};
          const Scalar f_low{this->m_f(low)}, f_mid{this->m_f(mid)},
              f_high{this->m_f(high)};

          // Check convergence on both the function values and the simplex size
          this->m_diameter = 0.0;
          for (Integer i{0}; i <= this->m_dim; ++i) {
            this->m_diameter =
                std::max(this->m_diameter,
                         (this->m_p.col(i) - this->m_p.col(low)).norm());
          }
          const Scalar f_range{f_high - f_low};
          if (this->m_verbose) {
            this->info(f_range, move_name(move));
          }
          const Scalar x_scale{
              static_cast<Scalar>(1.0) +
              this->m_p.col(low).template lpNorm<Eigen::Infinity>()};
          if (f_range <= tol * (1.0 + std::abs(f_low)) &&
              this->m_diameter <= tol * x_scale) {
            this->m_converged = true;
            break;
          }

          // Restart the simplex if it degenerated (flat simplex)
          if (this->m_dim > 1 && this->m_diameter > 0 &&
              this->normalized_volume(low) < this->m_volume_tolerance) {
            const Dense x_best(this->m_p.col(low));
            this->initialize_simplex(eval, x_best, this->m_diameter);
            move = Move::RESTART;
            continue;
          }

          // Centroid of all points but the worst
          centroid = (this->m_p.rowwise().sum() - this->m_p.col(high)) /
                     static_cast<Scalar>(this->m_dim);
          const Dense direction{centroid - this->m_p.col(high)};

          // Reflection
          x_r = centroid + rho * direction;
          const Scalar f_r{eval(x_r)};
          if (f_r < f_low) {
            // Expansion
            x_e = centroid + rho * chi * direction;
            const Scalar f_e{eval(x_e)};
            if (f_e < f_r) {
              this->replace(high, x_e, f_e);
              move = Move::EXPAND_E;
            } else {
              this->replace(high, x_r, f_r);
              move = Move::EXPAND_R;
            }
          } else if (f_r < f_mid) {
            this->replace(high, x_r, f_r);
            move = Move::REFLECT;
          } else if (f_r < f_high) {
            // Outside contraction
            x_c = centroid + gamma * rho * direction;
            const Scalar f_c{eval(x_c)};
            if (f_c <= f_r) {
              this->replace(high, x_c, f_c);
              move = Move::CONTRACT_O;
            } else {
              this->shrink(eval, low, sigma);
              move = Move::SHRINK;
            }
          } else {
            // Inside contraction
            x_c = centroid - gamma * direction;
            const Scalar f_c{eval(x_c)};
            if (f_c < f_high) {
              this->replace(high, x_c, f_c);
              move = Move::CONTRACT_I;
            } else {
              this->shrink(eval, low, sigma);
              move = Move::SHRINK;
            }
          }
        }

        // Print bottom
        if (this->m_verbose) {
          this->bottom();
        }

        // Convergence data
        this->sort();
        x_sol = Utils::from_dense<Vector>(this->m_p.col(this->m_order.front()));
        return this->m_converged;

#undef CMD
      }

     private:
      /**
       * Get the name of a simplex move.
       */
      static std::string move_name(const Move move) {
        switch (move) {
          case Move::INITIALIZE:
            return "Initialize";
          case Move::REFLECT:
            return "Reflect";
          case Move::EXPAND_E:
            return "Expand";
          case Move::EXPAND_R:
            return "Expand (reflection)";
          case Move::CONTRACT_O:
            return "Contract (outside)";
          case Move::CONTRACT_I:
            return "Contract (inside)";
          case Move::SHRINK:
            return "Shrink";
          case Move::RESTART:
            return "Restart";
          default:
            return "Unknown";
        }
      }

      /**
       * Build the initial simplex around a point.
       * \param[in] eval Function evaluator.
       * \param[in] x Simplex base point.
       * \param[in] delta Simplex size.
       */
      template <typename Evaluator>
      void initialize_simplex(Evaluator &&eval,
                              const Dense &x,
                              const Scalar delta) {
        const Integer n{this->m_dim};
        this->m_p.resize(n, n + 1);
        this->m_f.resize(n + 1);
        this->m_p.col(0) = x;
        this->m_f(0)     = eval(x);
        if (this->m_method == Method::SPENDLEY) {
          // Regular simplex with edge length delta
          const Scalar sn{static_cast<Scalar>(n)};
          const Scalar t{std::sqrt(sn + 1)};
          const Scalar s2{std::sqrt(static_cast<Scalar>(2.0))};
          const Scalar p{delta * (sn - 1 + t) / (sn * s2)};
          const Scalar q{delta * (t - 1) / (sn * s2)};
          for (Integer i{0}; i < n; ++i) {
            this->m_p.col(i + 1) = x.array() + q;
            this->m_p(i, i + 1)  = x(i) + p;
            this->m_f(i + 1)     = eval(this->m_p.col(i + 1));
          }
        } else {
          // Axis-aligned simplex, choosing the best of +delta and -delta
          for (Integer i{0}; i < n; ++i) {
            Dense x_p(x), x_m(x);
            x_p(i) += delta;
            x_m(i) -= delta;
            const Scalar f_p{eval(x_p)}, f_m{eval(x_m)};
            if (f_p <= f_m) {
              this->m_p.col(i + 1) = x_p;
              this->m_f(i + 1)     = f_p;
            } else {
              this->m_p.col(i + 1) = x_m;
              this->m_f(i + 1)     = f_m;
            }
          }
        }
        this->m_order.resize(n + 1);
        std::iota(this->m_order.begin(), this->m_order.end(), Integer(0));
      }

      /**
       * Sort the simplex points by increasing function value.
       */
      void sort() {
        std::stable_sort(this->m_order.begin(),
                         this->m_order.end(),
                         [this](const Integer a, const Integer b) {
                           return this->m_f(a) < this->m_f(b);
                         });
      }

      /**
       * Replace a simplex point.
       */
      void replace(const Integer j, const Dense &x, const Scalar f) {
        this->m_p.col(j) = x;
        this->m_f(j)     = f;
      }

      /**
       * Shrink the simplex towards the best point.
       */
      template <typename Evaluator>
      void shrink(Evaluator &&eval, const Integer low, const Scalar sigma) {
        for (Integer i{0}; i <= this->m_dim; ++i) {
          if (i != low) {
            this->m_p.col(i) =
                this->m_p.col(low) +
                sigma * (this->m_p.col(i) - this->m_p.col(low));
            this->m_f(i) = eval(this->m_p.col(i));
          }
        }
      }

      /**
       * Volume of the simplex divided by the volume of a regular simplex
       * whose edge equals the current diameter. It is 1 for a regular
       * simplex and tends to 0 as the simplex flattens.
       */
      Scalar normalized_volume(const Integer low) const {
        const Integer n{this->m_dim};
        DenseMatrix edges(n, n);
        Integer j{0};
        for (Integer i{0}; i <= n; ++i) {
          if (i != low) {
            edges.col(j++) = (this->m_p.col(i) - this->m_p.col(low)) /
                             this->m_diameter;
          }
        }
        // Regular simplex with unit edge: sqrt(n+1) / (n! * 2^(n/2)), the
        // factorial cancels out in the ratio
        const Scalar sn{static_cast<Scalar>(n)};
        const Scalar regular{std::sqrt(sn + 1) /
                             std::pow(static_cast<Scalar>(2.0), sn / 2)};
        return std::abs(edges.determinant()) / regular;
      }

    };  // class NelderMead

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_NELDER_MEAD_HH
