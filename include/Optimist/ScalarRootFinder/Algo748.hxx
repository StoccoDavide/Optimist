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

#ifndef OPTIMIST_SCALAR_ROOT_FINDER_ALGO748_HXX
#define OPTIMIST_SCALAR_ROOT_FINDER_ALGO748_HXX

namespace Optimist
{
  namespace ScalarRootFinder
  {

    /*\
     |      _    _           _____ _  _    ___
     |     / \  | | __ _  __|___  | || |  ( _ )
     |    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
     |   / ___ \| | (_| | (_) / / |__   _| (_) |
     |  /_/   \_\_|\__, |\___/_/     |_|  \___/
     |             |___/
    \*/

    /**
    * \brief Class container for the Algorithm 748.
    *
    * The algorithm 748 allows to find the roots of a scalar function \f$f(x)\f$ in a given interval
    * \f$[a, b]\f$. The algorithm is based on the work of G. Alefeld, F. Potra, and Y. Shi, *Algorithm
    * 748: Enclosing Zeros of Continuous Functions*, ACM Transactions on Mathematical Software, 21
    * (1995), pp. 327-344, 10.1145/210089.210111.
    */
    class Algo748 : public Bracketing<Algo748>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{false};
      static constexpr bool requires_second_derivative{false};

      // Function types
      using FunctionWrapper         = typename Bracketing<Algo748>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename Bracketing<Algo748>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename Bracketing<Algo748>::SecondDerivativeWrapper;

      /**
      * Class constructor for the Algorithm 748.
      */
      Algo748() {}

      /**
      * Get the Algorithm 748 solver name.
      * \return The Algorithm 748 solver name.
      */
      std::string name_impl() const {return "Algo748";}

    private:
      /**
      * Check if the input values are all different.
      * \param[in] a The first value.
      * \param[in] b The second value.
      * \param[in] c The third value.
      * \param[in] d The fourth value.
      * \return The boolean flag, true if all the values are different, false otherwise.
      */
      bool all_different(Real a, Real b, Real c, Real d) const {
        return a != b && a != c && a != d && b != c && b != d && c != d;
      }

      /**
      * Given current enclosing interval \f$[a, b]\f$ and \f$a\f$ number \f$c\f$ in \f$(a, b)\f$:
      *   1. if \f$f(c) = 0\f$ then sets the output \f$a = c\f$.
      *   2. Otherwise determines the new enclosing interval \f$[a, b] = [a, c]\f$ or \f$[a, b] = [c, b]\f$.
      *      also updates the termination criterion corresponding
      *      to the new enclosing interval.
      * Adjust \f$c\f$ if \f$(b-a)\f$ is very small or if \f$c\f$ is very close to \f$a\f$ or \f$b\f$.
      * \param[in] function Function wrapper.
      * \return The boolean flag, true if the root is found, false otherwise.
      */
      bool bracketing(FunctionWrapper function) {

        {
          Real tolerance{Real(0.7)*this->m_tolerance_bracketing};
          Real hba{(this->m_b - this->m_a)/2.0};
          if (hba <= tolerance) {this->m_c = this->m_a + hba;}
          else if (this->m_c <= m_a + tolerance) {m_c = m_a + tolerance;}
          else if (this->m_c >= m_b - tolerance) {m_c = m_b - tolerance;}
        }

        // If f(c) = 0 set a = c and return true
        this->evaluate_function(function, this->m_c, this->m_fc);
        if (this->m_fc == 0) {
          this->m_a = m_c; this->m_fa = 0.0;
          this->m_d = 0.0; this->m_fd = 0.0;
          return true;
        } else {
          // If f(c) is not zero, then determine the new enclosing interval
          if (this->m_fa*this->m_fc < 0) {
            // D <- B <- C
            this->m_d = this->m_b; this->m_fd = this->m_fb;
            this->m_b = this->m_c; this->m_fb = this->m_fc;
          } else {
            // D <- A <- C
            this->m_d = this->m_a; this->m_fd = this->m_fa;
            this->m_a = this->m_c; this->m_fa = this->m_fc;
          }
          return false;
        }
      }

      /**
      * Perform the cubic inverse interpolation of \f$f(x)\f$ at \f$a\f$, \f$b\f$, \f$d\f$, and
      * \f$e\f$ to get an approximate root \f$r\f$ (\f$f(r) = 0\f$).
      * \return The approximate root \f$r\f$.
      */
      Real cubic_interpolation() {

        // Compute the coefficients for the inverse cubic interpolation
        Real d_1{this->m_b - this->m_a};
        Real d_2{this->m_d - this->m_a};
        Real d_3{this->m_e - this->m_a};

        Real dd_0{d_1/(this->m_fb - this->m_fa)};
        Real dd_1{(d_1-d_2)/(this->m_fb - this->m_fd)};
        Real dd_2{(d_2-d_3)/(this->m_fd - this->m_fe)};

        Real ddd_0 {(dd_0 - dd_1)/(this->m_fa - this->m_fd)};
        Real ddd_1 {(dd_1 - dd_2)/(this->m_fb - this->m_fe)};

        Real dddd_0{(ddd_0-ddd_1)/(this->m_fa-this->m_fe)};

        // Compute the approximate root
        Real root{this->m_a - this->m_fa*(dd_0 - this->m_fb*(ddd_0 - this->m_fd*dddd_0))};
        Real tolerance{0.7*this->m_tolerance_bracketing};
        if (root <= this->m_a + tolerance || root >= this->m_b - tolerance) {
          root = (this->m_a + this->m_b)/2.0;
        }
        return root;
      }

      /**
      * Perform \f$n\f$ Newton steps to approximate the zero in \f$(a, b)\f$ of the quadratic
      * polynomial interpolating \f$f(x)\f$ at \f$a\f$, \f$b\f$, and \f$d\f$. Safeguard is used to
      * avoid overflow.
      * \param[in] n The number of steps.
      * \return The approximate root.
      */
      Real newton_quadratic(Integer n) {

        // Compute the coefficients for the quadratic interpolation
        Real a_0{this->m_fa};
        Real a_1{(this->m_fb - this->m_fa)/(this->m_b - this->m_a)};
        Real a_2{((this->m_fd - this->m_fb)/(this->m_d - this->m_b)-a_1)/(this->m_d - this->m_a)};

        // Safeguard to avoid overflow
        if (a_2 == 0) {return this->m_a - a_0/a_1;}

        // Determine the starting point of newton steps
        Real c{a_2*this->m_fa > 0.0 ? this->m_a : this->m_b};

        // Start the safeguarded newton steps
        bool is_nonzero{true};
        Real pc{0.0}, pdc{0.0};
        for (Integer i{0}; i < n && is_nonzero ; ++i ) {
          pc  = a_0 + (a_1 + a_2*(c - this->m_b))*(c - this->m_a);
          pdc = a_1 + a_2*((2.0 * c) - (this->m_a + this->m_b));
          is_nonzero = pdc != 0;
          if (is_nonzero) {c -= pc/pdc;}
        }

        if (is_nonzero) {return c;}
        else {return this->m_a - a_0/a_1;}
      }

    public:
      /**
      * Finds either an exact solution or an approximate solution of the equation \f$f(x) = 0\f$ in
      * the interval \f$[a,b]\f$. At the beginning of each iteration, the current enclosing interval
      * is recorded as \f$[a_0, b_0]\f$. The first iteration is simply a secant step. Starting with
      * the second iteration, three steps are taken in each iteration.
      *   1. The first two steps are either quadratic interpolation or cubic inverse interpolation.
      *   2. The third step is a double-size secant step.
      * If the diameter of the enclosing interval obtained after these three steps is larger than
      * \f$\mu*(b_0-a_0)\f$, an additional bisection step will be used to shrink the enclosing interval.
      * \param[in] function Function wrapper.
      * \return The approximate root.
      */
      Real find_root_impl(FunctionWrapper function) {

        Real tolerance_step{this->m_tolerance_bracketing};
        Real tolerance_function{this->m_tolerance_bracketing};

        // Check for trivial solution
        this->m_converged = this->m_fa == 0; if (this->m_converged) {return this->m_a;}
        this->m_converged = this->m_fb == 0; if (this->m_converged) {return this->m_b;}

        // Initialize to dumb values
        this->m_e  = QUIET_NAN;
        this->m_fe = QUIET_NAN;

        // While f(left) or f(right) are infinite perform bisection
        while (!(std::isfinite(m_fa) && std::isfinite(m_fb))) {
          ++this->m_iterations;
          this->m_c  = (this->m_a + this->m_b)/2.0;
          this->evaluate_function(function, this->m_c, this->m_fc);
          this->m_converged = this->m_fc == 0;
          if (this->m_converged) {return this->m_c;}
          if (this->m_fa*this->m_fc < 0.0) {
            // -> [a, c]
            this->m_b = this->m_c; this->m_fb = this->m_fc;
          } else {
            // -> [c, b]
            this->m_a = this->m_c; this->m_fa = this->m_fc;
          }

          // Check for convergence
          Real abs_fa{std::abs(this->m_fa)};
          Real abs_fb{std::abs(this->m_fb)};
          this->m_converged = (this->m_b - this->m_a) < tolerance_step ||
                               abs_fa < tolerance_function ||abs_fb < tolerance_function;
          if (this->m_converged) {return abs_fb < abs_fa ? this->m_b : this->m_a;}
        }

        {
          // Impede the step c to be too close to the bounds a or b.
          Real b_a{this->m_b  - this->m_a};
          Real r{b_a/(this->m_fb - this->m_fa)};
          this->m_c = std::abs(this->m_fb) < std::abs(this->m_fa) ?
            this->m_b + this->m_fb*r : this->m_a - this->m_fa*r;
          Real delta{this->m_interval_shink*b_a}, a_min, b_max;
          if (this->m_c < (a_min = this->m_a + delta)) {this->m_c = a_min;}
          else if (this->m_c > (b_max = this->m_b - delta)) {this->m_c = b_max;}
        }

        // Call "bracketing" to get a shrinked enclosing interval and to update the termination criterion
        this->m_converged = this->bracketing(function);
        if (this->m_converged ) {return this->m_a;}
        this->m_converged = false;

        // The enclosing interval is recorded as [a0, b0] before executing the iteration steps
        while (++this->m_iterations < this->m_max_iterations) {
          ++this->m_iterations;
          Real a0_b0{this->m_b - this->m_a};

          // Check the termination criterion
          {
            Real abs_fa {std::abs(this->m_fa)};
            Real abs_fb {std::abs(this->m_fb)};
            this->m_converged = abs_fa < tolerance_function || abs_fb < tolerance_function;
            if (this->m_converged) {return abs_fa < abs_fb ? this->m_a : this->m_b;}
          }

          // Starting with the second iteration: newton quadratic or cubic interpolation
          bool do_newton_quadratic{false};
          if (!std::isfinite(this->m_fe)) {
            do_newton_quadratic = true;
          } else if (!this->all_different( this->m_fa, this->m_fb, this->m_fd, this->m_fe)) {
            do_newton_quadratic = true;
          } else {
            this->m_c = this->cubic_interpolation();
            do_newton_quadratic = (this->m_c-this->m_a)*(this->m_c-this->m_b) >= 0.0;
          }
          if (do_newton_quadratic) {this->m_c = this->newton_quadratic(2);}

          this->m_e  = this->m_d;
          this->m_fe = this->m_fd;

          // Call "bracketing" to get a shrinked enclosing interval and to update the termination criterion
          this->m_converged = this->bracketing(function) || (this->m_b-this->m_a) <= this->m_tolerance;
          if (this->m_converged) {
            return std::abs(this->m_fa) < std::abs(this->m_fb) ? this->m_a : this->m_b;
          }
          if (!this->all_different(this->m_fa, this->m_fb, this->m_fd, this->m_fe)) {
            do_newton_quadratic = true;
          } else {
            this->m_c = this->cubic_interpolation();
            do_newton_quadratic = (this->m_c - this->m_a)*(this->m_c - this->m_b) >= 0.0;
          }
          if (do_newton_quadratic) {this->m_c = this->newton_quadratic(3);}

          // Call "bracketing" to get a shrinked enclosing interval and to update the termination criterion
          {
            this->m_converged = this->bracketing(function) || (this->m_b-this->m_a) <= this->m_tolerance;
            Real abs_fa{std::abs(this->m_fa)};
            Real abs_fb{std::abs(this->m_fb)};
            if (this->m_converged) {return abs_fa < abs_fb ? this->m_a : this->m_b;}

            this->m_e  = this->m_d;
            this->m_fe = this->m_fd;

            // Takes the double-size secant step
            Real u, fu;
            if (abs_fa < abs_fb) {u = this->m_a; fu = this->m_fa;}
            else {u = this->m_b; fu = this->m_fb;}
            Real hba{(this->m_b - this->m_a)/2.0};
            this->m_c = u - 4.0*(fu/(this->m_fb - this->m_fa))*hba;
            if (std::abs(this->m_c - u) > hba) {this->m_c = this->m_a + hba;}
          }

          // Call "bracketing" to get a shrinked enclosing interval and to update the termination criterion
          this->m_converged = this->bracketing(function) || (this->m_b-this->m_a) <= this->m_tolerance_bracketing;
          if (this->m_converged) {
            return std::abs(this->m_fa) < std::abs(this->m_fb) ? this->m_a : this->m_b;
          }

          // Determine whether an additional bisection step is needed
          if ((this->m_b-this->m_a) < this->m_mu*a0_b0) {continue;}

          this->m_e  = this->m_d;
          this->m_fe = this->m_fd;

          // Call "bracketing" to get a shrinked enclosing interval and to update the termination criterion
          {
            Real b_a{this->m_b-this->m_a};
            this->m_c = this->m_a + b_a/2.0;
            this->m_converged = this->bracketing(function) || b_a <= this->m_tolerance_bracketing;
          }
        }
        // Return the approximate root
        return std::abs(this->m_fa) < std::abs(this->m_fb) ? this->m_a : this->m_b;
      }

    }; // class Algo748

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_ALGO748_HXX
