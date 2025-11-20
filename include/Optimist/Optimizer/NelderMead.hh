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

#ifndef OPTIMIST_OPTIMIZER_NELDER_MEAD_HH
#define OPTIMIST_OPTIMIZER_NELDER_MEAD_HH

#include "Optimist/Optimizer.hh"

namespace Optimist
{
  namespace Optimizer
  {

    /*\
     |   _   _      _     _           __  __                _
     |  | \ | | ___| | __| | ___ _ __|  \/  | ___  __ _  __| |
     |  |  \| |/ _ \ |/ _` |/ _ \ '__| |\/| |/ _ \/ _` |/ _` |
     |  | |\  |  __/ | (_| |  __/ |  | |  | |  __/ (_| | (_| |
     |  |_| \_|\___|_|\__,_|\___|_|  |_|  |_|\___|\__,_|\__,_|
     |
    \*/


    /**
     * \brief Class container for the Nelder-Mead's method.
     *
     * \includedoc docs/markdown/Optimizer/NelderMead.md
     *
     * \tparam Real Scalar number type.
     * \tparam N Dimension of the root-finding problem.
     */
    template <typename Real, Integer N>
    class NelderMead : public Optimizer<Real, N, NelderMead<Real, N>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{false};
      static constexpr bool requires_second_derivative{false};

      using typename Optimizer<Real, N, NelderMead<Real, N>>::Vector;
      using typename Optimizer<Real, N, NelderMead<Real, N>>::Matrix;

      using VectorF = Eigen::Vector<Real, N+1>;
      using MatrixP = Eigen::Matrix<Real, N, N+1>;
      using MatrixD = Eigen::Matrix<Real, N+1, N+1>;

      /**
       * NelderMead move type enumeration.
       */
      using Move = enum class Move : Integer {INITIALIZE = 0, REFLECT = 1, EXPAND_E = 2, EXPAND_R = 3,
        CONTRACT_O = 4, CONTRACT_I = 5, SHRINK = 6, RESTART = 7};

      /**
       * NelderMead solver type enumeration.
       */
      using Method = enum class Method : Integer {STANDARD = 0, SPENDLEY = 1};

    private:
      // Solver parameters
      Method m_method{Method::STANDARD}; /**< NelderMead solver method. */
      Real   m_delta{1.0};               /**< Initial simplex size. */
      Real   m_rho{1.0};                 /**< Reflection coefficient. */
      Real   m_chi{2.0};                 /**< Expansion coefficient. */
      Real   m_gamma{0.5};               /**< Contraction coefficient. */
      Real   m_sigma{0.25};              /**< Shrink coefficient. */
      Real   m_volume_tolerance{std::sqrt(this->m_tolerance)}; /**< Simplex volume tolerance. */

      // Internal variables
      Real m_dim_factorial{0};          /**< Factorial of the problem dimension. */
      Real m_regular_simplex_volume{0}; /**< Volume of a regular simplex. */

      Integer m_low{0};  /**< Index of the lowest (best) point. */
      Integer m_mid{0};  /**< Index of the middle point. */
      Integer m_high{0}; /**< Index of the highest (worst) point. */

      VectorF m_f;      /**< Function values at simplex points. */
      MatrixP m_p;      /**< Simplex points. */
      MatrixD m_dist;   /**< Distance matrix between simplex points. */
      Vector  m_f_work; /**< Work function values. */
      Matrix  m_p_work; /**< Work simplex points. */
      Vector  m_p_sum;  /**< Sum of simplex points. */

      Eigen::PartialPivLU<Matrix> m_lu; /**< LU decomposition. */

      Vector m_p_r; /**< Reflection point. */
      Vector m_p_e; /**< Expansion point. */
      Vector m_p_c; /**< Contraction point. */

      Real m_diameter{0};       /**< Simplex diameter. */
      Real m_simplex_volume{0}; /**< Simplex volume. */

    public:
      /**
       * Class constructor for the NelderMead solver.
       */
      NelderMead()
      {
        // Precompute constants
        this->m_dim_factorial = 1;
        for (Integer i{2}; i <= N+1; ++i) {this->m_dim_factorial *= i;}
        this->m_regular_simplex_volume = (std::sqrt(static_cast<Real>(N+1)) /
          static_cast<Real>(this->m_dim_factorial)) / std::pow(2.0, static_cast<Real>(N)/2.0);

        // Initialize internal variables
        this->m_f.setZero();
        this->m_p.setZero();
        this->m_dist.setZero();
        this->m_f_work.setZero();
        this->m_p_work.setZero();
        this->m_p_sum.setZero();
        this->m_p_r.setZero();
        this->m_p_e.setZero();
        this->m_p_c.setZero();
      }

      /**
       * Class constructor for the NelderMead solver.
       * \param[in] delta Initial simplex size.
       */
       NelderMead(Real const delta) : NelderMead()
      {
        this->delta(delta);
      }

      /**
       * Get the NelderMead solver name.
       * \return The NelderMead solver name.
       */
      std::string name_impl() const {
        std::ostringstream os;
        os << "NelderMead";
        if (this->m_method == Method::STANDARD) {
          os << "Standard";
        } else if (this->m_method == Method::SPENDLEY) {
          os << "Spendley";
        }
        return os.str();
      }

      /**
       * Get the enumeration type of the NelderMead solver method.
       * \return The NelderMead solver enumeration type.
       */
      Method method() const {return this->m_method;}

      /**
       * Set the enumeration type of the NelderMead solver method.
       * \param[in] t_method The NelderMead solver enumeration type.
       */
      void method(Method t_method) {this->m_method = t_method;}

      /**
       * Get the initial simplex size.
       * \return The initial simplex size.
       */
      Real delta() const {return this->m_delta;}

      /**
       * Set the initial simplex size.
       * \param[in] t_delta The initial simplex size.
       */
      void delta(Real t_delta)
      {
        OPTIMIST_ASSERT(t_delta > 0,
          "Optimist::Optimizer::NelderMead::delta(...): invalid input detected.");
        this->m_delta = t_delta;
      }

      /**
       * Get the reflection coefficient.
       * \return The reflection coefficient.
       */
      Real rho() const {return this->m_rho;}

      /**
       * Set the reflection coefficient.
       * \param[in] t_rho The reflection coefficient.
       */
      void rho(Real t_rho)
      {
        OPTIMIST_ASSERT(t_rho > 0,
          "Optimist::Optimizer::NelderMead::rho(...): invalid input detected.");
        this->m_rho = t_rho;
      }

      /**
       * Get the expansion coefficient.
       * \return The expansion coefficient.
       */
      Real chi() const {return this->m_chi;}

      /**
       * Set the expansion coefficient.
       * \param[in] t_chi The expansion coefficient.
       */
      void chi(Real t_chi)
      {
        OPTIMIST_ASSERT(t_chi > 1,
          "Optimist::Optimizer::NelderMead::chi(...): invalid input detected.");
        this->m_chi = t_chi;
      }

      /**
       * Get the contraction coefficient.
       * \return The contraction coefficient.
       */
      Real gamma() const {return this->m_gamma;}

      /**
       * Set the contraction coefficient.
       * \param[in] t_gamma The contraction coefficient.
       */
      void gamma(Real t_gamma)
      {
        OPTIMIST_ASSERT(t_gamma > 0 && t_gamma < 1,
          "Optimist::Optimizer::NelderMead::gamma(...): invalid input detected.");
        this->m_gamma = t_gamma;
      }

      /**
       * Get the shrink coefficient.
       * \return The shrink coefficient.
       */
      Real sigma() const {return this->m_sigma;}

      /**
       * Set the shrink coefficient.
       * \param[in] t_sigma The shrink coefficient.
       */
      void sigma(Real t_sigma)
      {
        OPTIMIST_ASSERT(t_sigma > 0 && t_sigma < 1,
          "Optimist::Optimizer::NelderMead::sigma(...): invalid input detected.");
        this->m_sigma = t_sigma;
      }

      /**
       * Get the simplex volume tolerance.
       * \return The simplex volume tolerance.
       */
      Real volume_tolerance() const {return this->m_volume_tolerance;}

      /**
       * Set the simplex volume tolerance.
       * \param[in] t_volume_tolerance The simplex volume tolerance.
       */
      void volume_tolerance(Real t_volume_tolerance)
      {
        OPTIMIST_ASSERT(t_volume_tolerance > 0,
          "Optimist::Optimizer::NelderMead::volume_tolerance(...): invalid input detected.");
        this->m_volume_tolerance = t_volume_tolerance;
      }

      /**
       * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0 \f$, with \f$
       * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[in] jacobian Jacobian lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda && function, Vector const & x_ini, Vector & x_sol)
      {
        #define CMD "Optimist::Optimizer::NelderMead::solve_impl(...): "

        // Initialize simplex
        if (this->m_method == Method::STANDARD) {
          this->diamond(std::forward<FunctionLambda>(function), x_ini, this->m_delta);
        } else if (this->m_method == Method::SPENDLEY) {
          this->spendley(std::forward<FunctionLambda>(function), x_ini, this->m_delta);
        } else {
          OPTIMIST_ERROR(CMD "invalid method detected.");
        }

        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        this->dist_init();
        this->m_p_sum.noalias() = this->m_p.col(N);
        for (Integer i{0}; i < N; ++i) {this->m_p_sum.noalias() += this->m_p.col(i);}

        // Algorithm iterations
        Move move{Move::INITIALIZE};
        for (this->m_iterations = 1; this->m_iterations <= this->m_max_iterations; ++this->m_iterations)
        {
          // Determine which point is the highest (worst), next-highest, and lowest (best).
          this->m_low  = 0;
          this->m_mid  = 1;
          this->m_high = 2;
          if (this->m_f(this->m_low) > this->m_f(this->m_mid))  {std::swap(this->m_low, this->m_mid);}
          if (this->m_f(this->m_low) > this->m_f(this->m_high)) {std::swap(this->m_low, this->m_high);}
          if (this->m_f(this->m_mid) > this->m_f(this->m_high)) {std::swap(this->m_mid, this->m_high);}
          for (Integer i{3}; i <= N; ++i)
          {
            if (this->m_f(i) < this->m_f(this->m_low) ) {
              this->m_low = i; // New minima
            } else if (this->m_f(i) > this->m_f(this->m_high) ) {
              this->m_mid  = this->m_high;
              this->m_high = i;
            } else if (this->m_f(i) > this->m_f(this->m_mid) ) {
              this->m_mid = i;
            }
          }

          Real f_low{this->m_f(this->m_low)};
          Real f_mid{this->m_f(this->m_mid)};
          Real f_high{this->m_f(this->m_high)};

          // Compute the fractional range from highest to lowest and return if satisfactory
          Real rtol{std::abs(f_high - f_low) / (std::abs(f_low) + this->m_tolerance)};

          if (this->m_verbose) {
            std::string move_str;
            switch (move) {
              case Move::INITIALIZE: move_str = "Initialize"; break;
              case Move::REFLECT:    move_str = "Reflect"; break;
              case Move::EXPAND_E:   move_str = "Expand"; break;
              case Move::EXPAND_R:   move_str = "Expand (reflection)"; break;
              case Move::CONTRACT_O: move_str = "Contract (outside)"; break;
              case Move::CONTRACT_I: move_str = "Contract (inside)"; break;
              case Move::SHRINK:     move_str = "Shrink"; break;
              case Move::RESTART:    move_str = "Restart"; break;
              default:               move_str = "Unknown"; break;
            }
            this->info(rtol, move_str);
          }
          if (rtol < this->m_tolerance || this->m_diameter < this->m_tolerance) {
            this->m_converged = true;
            break;
          }

          Real ratio{std::pow(this->m_simplex_volume / this->m_regular_simplex_volume, 1.0 / N) /
            this->m_diameter};
          if (ratio <= this->m_volume_tolerance)
          {
            if (this->m_method == Method::STANDARD) {
              this->diamond(std::forward<FunctionLambda>(function), this->m_p.col(this->m_low),
                this->m_diameter*this->m_volume_tolerance);
            } else if (this->m_method == Method::SPENDLEY) {
              this->spendley(std::forward<FunctionLambda>(function), this->m_p.col(this->m_low),
                this->m_diameter*this->m_volume_tolerance);
            } else {
              OPTIMIST_ERROR(CMD "invalid method detected.");
            }
            this->dist_init();
            this->m_p_sum.noalias() = this->m_p.col(N);
            for (Integer i{0}; i < N; ++i) {this->m_p_sum.noalias() += this->m_p.col(i);}
            move = Move::RESTART;
            continue;
          }

          // Begin a new iteration: first extrapolate by a factor alpha (default: −1 a reflection)
          // through the face of the simplex across from the high point, i.e., reflect the simplex
          // from the high point
          Real f_r{this->reflect(std::forward<FunctionLambda>(function), this->m_p_r)};
          if (f_r < f_low) {
            Real f_e{this->expand(std::forward<FunctionLambda>(function), this->m_p_e)};
            if (f_e < f_r) {
              move = Move::EXPAND_E;
              this->replace_point(f_e, this->m_p_e, this->m_high);
            } else {
              move = Move::EXPAND_R;
              this->replace_point(f_r, this->m_p_r, this->m_high);
            }
          } else if (f_r < f_mid) {
            move = Move::REFLECT;
            this->replace_point(f_r, this->m_p_r, this->m_high);
          } else if (f_r < f_high) {
            Real f_c{this->outside(std::forward<FunctionLambda>(function), this->m_p_c)};
            if (f_c < f_r) {
              move = Move::CONTRACT_O;
              this->replace_point(f_c, this->m_p_c, this->m_high);
            } else {
              move = Move::REFLECT;
              this->replace_point(f_r, this->m_p_r, this->m_high);
            }
          } else {
            Real f_c{this->inside(std::forward<FunctionLambda>(function), this->m_p_c)};
            if (f_c < f_high) {
              move = Move::CONTRACT_I;
              this->replace_point(f_c, this->m_p_c, this->m_high);
            } else {
              move = Move::SHRINK;
              this->shrink(std::forward<FunctionLambda>(function));
            }
          }
        }

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        x_sol = this->m_p.col(this->m_low);
        return this->m_converged;

        #undef CMD
      }

    private:
      /**
       * Replace simplex point.
       * \param[in] f Function value at the new point.
       * \param[in] p New simplex point.
       * \param[in] j Index of the point to be replaced.
       */
      void replace_point(Real const f, Vector const & p, Integer const j)
      {
        this->m_f(j)               = f;
        this->m_p_sum             += p - this->m_p.col(j);
        this->m_p.col(j).noalias() = p;
        this->dist_update(j);
      }

      /**
       * Initialize simplex with Spendley method.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[in] x Initial point.
       * \param[in] delta Initial simplex size.
       */
      template <typename FunctionLambda>
      void spendley(FunctionLambda && function, Vector const x, Real const delta)
      {
        #define CMD "Optimist::Optimizer::NelderMead::spendley(...): "

        Real t1{std::sqrt(static_cast<Real>(N + 1)) - 1.0};
        Real t2{static_cast<Real>(N)*std::sqrt(static_cast<Real>(2))};
        Real p{delta * (N+t1)/t2};
        Real q{delta * t1/t2};
        bool success;
        for (Integer i{0}; i < N; ++i) {this->m_p(i,0) = x(i);}
        success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_p.col(0), this->m_f(0));
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the initial point.");
        for (Integer i{0}; i < N; ++i) {
          this->m_p.col(i+1) = this->m_p.col(0).array()+p;
          this->m_p(i, i+1) += q;
          success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_p.col(i+1), this->m_f(i+1));
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed at simplex point " << i+1 << ".");
        }

        #undef CMD
      }

      /**
       * Initialize simplex with Diamond method.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[in] x Initial point.
       * \param[in] delta Initial simplex size.
       */
      template <typename FunctionLambda>
      void diamond(FunctionLambda && function, Vector const x, Real const delta)
      {
        #define CMD "Optimist::Optimizer::NelderMead::diamond(...): "

        bool success;
        for (Integer i{0}; i < N; ++i) {this->m_p.row(i).setConstant(x(i));}
        success = this->evaluate_function(std::forward<FunctionLambda>(function), x, this->m_f(0));
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the initial point.");
        for (Integer i{0}; i < N; ++i)
        {
          Real f_p, f_m;
          this->m_p(i, i+1) = x(i) + delta;
          success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_p.col(i+1), f_p);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed at simplex point " << i+1 << ".");
          this->m_p(i, i+1) = x(i) - delta;
          success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_p.col(i+1), f_m);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed at simplex point " << i+1 << ".");
          if (f_p < f_m) {
            this->m_f(i+1)    = f_p;
            this->m_p(i, i+1) = x(i) + delta;
          } else {
            this->m_f(i+1) = f_m;
          }
        }

        #undef CMD
      }

      /**
       * Compute simplex volume.
       * \param[in] k Index of the point to be used as reference.
       */
      void simplex_volume(Integer k)
      {
        Integer j{0};
        for (Integer i{0}; i <= N; ++i)
        {
          if (i != k) {
            this->m_p_work.col(j).noalias() = this->m_p.col(i) - this->m_p.col(k);
            this->m_f_work(j)               = this->m_f(i) - this->m_f(k);
            ++j;
          }
        }
        this->m_lu.compute(this->m_p_work.transpose());
        this->m_simplex_volume = std::abs(this->m_lu.determinant()) / this->m_dim_factorial;
      }

      /**
       * Initialize distance matrix.
       */
      void dist_init()
      {
        for (Integer i{0}; i < N; ++i) {
          this->m_dist(i, i) = 0;
          for (Integer j{i+1}; j <= N; ++j) {
            this->m_dist(i, j) = (this->m_p.col(i) - this->m_p.col(j)).norm();
            this->m_dist(j, i) = this->m_dist(i, j);
          }
        }
        this->m_diameter = this->m_dist.maxCoeff();
        this->simplex_volume(0);
      }

      /**
       * Update distance matrix.
       * \param[in] j Index of the updated point.
       */
      void dist_update(Integer j)
      {
        for (Integer i{0}; i <= N; ++i) {
          if (i != j) {
            this->m_dist(i, j) = (this->m_p.col(i) - this->m_p.col(j) ).norm();
            this->m_dist(j, i) = this->m_dist(i, j);
          }
        }
        this->m_diameter = this->m_dist.maxCoeff();
        this->simplex_volume(this->m_low);
      }

      /**
       * Shrink the simplex towards the lowest (best) point.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       */
      template <typename FunctionLambda>
      void shrink(FunctionLambda && function)
      {
        #define CMD "Optimist::Optimizer::NelderMead::shrink(...): "

        bool success;
        Real const c_1{1.0 - this->m_sigma}, c_2{this->m_sigma};
        for (Integer i{0}; i <= N; ++i)
        {
          if (i != this->m_low) {
            this->m_p.col(i) = c_1 * this->m_p.col(i) + c_2 * this->m_p.col(this->m_low);
            success = this->evaluate_function(std::forward<FunctionLambda>(function), this->m_p.col(i), this->m_f(i));
            OPTIMIST_ASSERT(success,
              CMD "function evaluation failed at simplex point " << i << ".");
          }
        }
        this->m_p_sum.noalias() = this->m_p.col(N);
        for (Integer i{0}; i < N; ++i) {this->m_p_sum.noalias() += this->m_p.col(i);}
        this->m_dist           *= c_1;
        this->m_diameter       *= c_1;
        this->m_simplex_volume *= std::pow(c_1, static_cast<Real>(N));
        #undef CMD
      }

      /**
       * Extrapolate simplex point by a factor alpha through the face of the simplex.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[in] alpha Extrapolation factor.
       * \param[in] j Index of the point to be extrapolated.
       * \param[out] p Extrapolated simplex point.
       * \return Function value at the extrapolated point.
       */
      template <typename FunctionLambda>
      Real extrapolate(FunctionLambda && function, Real const alpha, Integer const j, Vector & p)
      {
        // Extrapolates by a factor alpha through the face of the simplex
        // across from the high point.
        //
        // pface   = (sum(p)-pfrom)/N
        // pcenter = sum(p)/(N+1);
        // pe      = pface + alpha*(pface-pfrom)
        //         = pface*(1+alpha)-pfrom*alpha
        //         = (sum(p)-pfrom)/N*(1+alpha)-pfrom*alpha
        //         = sum(p)/N*(1+alpha)-pfrom*(alpha+(1+alpha)/N)
        //         = sum(p)/(N+1) * (1+alpha)*((N+1)/N) - pfrom*[1/N+alpha*(N+1)/N]
        //         = pcenter * (1+alpha)*((N+1)/N) - pfrom*(1/N-(N+1)/N+(1+alpha)*(N+1)/N)
        //         = pcenter * (1+alpha)*((N+1)/N) - pfrom*(-1+(1+alpha)*(N+1)/N)
        //
        // let beta = (1+alpha)*((N+1)/N);
        //
        // pe = pcenter * beta + pfrom * (1-beta)

        #define CMD "Optimist::Optimizer::NelderMead::extrapolate(...): "

        Real const beta_1{(1.0 + alpha)/static_cast<Real>(N)};
        Real const beta{(static_cast<Real>(N) + 1.0)*beta_1};
        p.noalias() = beta_1*this->m_p_sum + (1.0 - beta) * this->m_p.col(j);
        Real f;
        bool success{this->evaluate_function(std::forward<FunctionLambda>(function), p, f)};
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed during extrapolation.");
        return f;

        #undef CMD
      }

      /**
       * Reflect simplex point.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[out] p Reflected simplex point.
       * \return Function value at the reflected point.
       */
      template <typename FunctionLambda>
      Real reflect(FunctionLambda && function, Vector & p)
      {
        return this->extrapolate(std::forward<FunctionLambda>(function), this->m_rho, this->m_high, p);
      }

      /**
       * Expand simplex point.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[out] p Expanded simplex point.
       * \return Function value at the expanded point.
       */
      template <typename FunctionLambda>
      Real expand(FunctionLambda && function, Vector & p)
      {
        return this->extrapolate(std::forward<FunctionLambda>(function), this->m_rho*this->m_chi, this->m_high, p);
      }

      /**
       * Outside contraction simplex point.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[out] p Outside contracted simplex point.
       * \return Function value at the outside contracted point.
       */
      template <typename FunctionLambda>
      Real outside(FunctionLambda && function, Vector & p)
      {
        return this->extrapolate(std::forward<FunctionLambda>(function), this->m_rho*this->m_gamma, this->m_high, p);
      }

      /**
       * Inside contraction simplex point.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \param[out] p Inside contracted simplex point.
       * \return Function value at the inside contracted point.
       */
      template <typename FunctionLambda>
      Real inside(FunctionLambda && function, Vector & p)
      {
        return this->extrapolate(std::forward<FunctionLambda>(function), -this->m_gamma, this->m_high, p);
      }

    }; // class NelderMead

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_NELDER_MEAD_HH
