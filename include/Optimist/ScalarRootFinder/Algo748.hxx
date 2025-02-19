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
    * \includedoc docs/markdown/ScalarRootFinder/Algo748.md
    */
    class Algo748 : public ScalarRootFinder<Algo748>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{false};
      static constexpr bool requires_second_derivative{false};

      // Function types
      using FunctionWrapper         = typename ScalarRootFinder<Algo748>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename ScalarRootFinder<Algo748>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename ScalarRootFinder<Algo748>::SecondDerivativeWrapper;

    private:
      Real m_mu{0.5};
      Real m_interval_shink{0.025};

      Real m_a{0.0}, m_fa{0.0};
      Real m_b{0.0}, m_fb{0.0};
      Real m_c{0.0}, m_fc{0.0};
      Real m_d{0.0}, m_fd{0.0};
      Real m_e{0.0}, m_fe{0.0};

    public:
      /**
      * Class constructor for the Algorithm 748.
      */
      Algo748() {}

      /**
      * Get the Algorithm 748 solver name.
      * \return The Algorithm 748 solver name.
      */
      std::string name_impl() const {return "Algo748";}

      void set_tolerance( Real B ) {
        this->tolerance(2.0*EPSILON + 4.0*std::abs(B)*EPSILON);
      }

      bool bracketing(FunctionWrapper function) {

        /*
        Given current enclosing interval [a,b] and a number c in (a,b):
         a) if f(c)=0 then sets the output a=c.
         b) Otherwise determines the new enclosing interval:
            [a,b]=[a,c] or [a,b]=[c,b].
            also updates the termination criterion corresponding
            to the new enclosing interval.
        Adjust c if (b-a) is very small or if c is very close to a or b.
        */
        {
          Real tol{0.7*this->m_tolerance};
          Real hba{(this->m_b - this->m_a)/2.0};
          if (hba <= tol) {this->m_c = this->m_a + hba;}
          else if (this->m_c <= this->m_a+tol) {this->m_c = this->m_a + tol;}
          else if (this->m_c >= this->m_b-tol) {this->m_c = this->m_b - tol;}
        }

        // Evaluate f(c)
        this->evaluate_function(function, this->m_c, this->m_fc);

        // If f(c)=0, then set a=c and return, this will terminate the procedure.
        if (this->m_fc == 0) {
          this->m_a = this->m_c; this->m_fa = 0;
          this->m_d = 0;   this->m_fd = 0;
          return true;
        } else {
          // If f(c) is not zero, then determine the new enclosing interval
          if (this->m_fa*this->m_fc < 0.0) {
            // D <-- B <-- C
            this->m_d = this->m_b; this->m_fd = this->m_fb;
            this->m_b = this->m_c; this->m_fb = this->m_fc;
          } else {
            // D <-- A <-- C
            this->m_d = this->m_a; this->m_fd = this->m_fa;
            this->m_a = this->m_c; this->m_fa = this->m_fc;
          }
          // Update the termination criterion according to the new enclosing interval
          if (std::abs(this->m_fb) <= std::abs(this->m_fa)) {this->set_tolerance(this->m_b);}
          else {this->set_tolerance(this->m_a);}
          return false;
        }
      }

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      // Uses cubic inverse interpolation of f(x) at a, b, d, and e to get an approximate root of f(x).
      Real pzero() {
        Real D1 = m_b-m_a;
        Real D2 = m_d-m_a;
        Real D3 = m_e-m_a;

        Real DD0 = D1/(m_fb-m_fa);
        Real DD1 = (D1-D2)/(m_fb-m_fd);
        Real DD2 = (D2-D3)/(m_fd-m_fe);

        Real DDD0  = (DD0-DD1)/(m_fa-m_fd);
        Real DDD1  = (DD1-DD2)/(m_fb-m_fe);

        Real DDDD0 = (DDD0-DDD1)/(m_fa-m_fe);

        Real c = m_a - m_fa*(DD0-m_fb*(DDD0-m_fd*DDDD0));

        Real tol = Real(0.7)*m_tolerance;
        if ( c <= m_a+tol || c >= m_b-tol ) c = (m_a+m_b)/2;

        return c;
      }

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Real newton_quadratic( Integer niter ) {
        // Uses `niter` newton steps to approximate the zero in (a,b) of the
        // quadratic polynomial interpolating f(x) at a, b, and d.
        // Safeguard is used to avoid overflow.

        Real A0 = m_fa;
        Real A1 = (m_fb-m_fa)/(m_b-m_a);
        Real A2 = ((m_fd-m_fb)/(m_d-m_b)-A1)/(m_d-m_a);

        // Safeguard to avoid overflow.
        if ( A2 == 0 ) return m_a-A0/A1;

        // Determine the starting point of newton steps.
        Real c = A2*m_fa > 0 ? m_a : m_b;

        // Start the safeguarded newton steps.
        bool ok = true;
        for ( Integer i=0; i < niter && ok ; ++i ) {
          Real PC  = A0+(A1+A2*(c-m_b))*(c-m_a);
          Real PDC = A1+A2*((2*c)-(m_a+m_b));
          ok = PDC != 0;
          if ( ok ) c -= PC/PDC;
        }

        if ( ok ) return c;
        else      return m_a-A0/A1;
      }

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      bool solve_impl(FunctionWrapper function, Real x_ini, Real &x_sol, Real a = -INFTY, Real b = INFTY)
      {
        // Setup internal variables
        this->reset();

        #define CMD "Optimist::ScalarRootFinder::Algo748::solve(...): "
          OPTIMIST_ASSERT(a < b, CMD "invalid interval.");
          OPTIMIST_ASSERT(a < x_ini && x_ini < b, CMD "invalid initial point.");
        #undef CMD

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        this->m_a = a; this->evaluate_function(function, this->m_a, this->m_fa);
        this->m_b = b; this->evaluate_function(function, this->m_b, this->m_fb);

        // Check if the solution exists
        if (this->m_fa*this->m_fb > 0.0) {return this->m_a;}
        else {x_sol = this->eval(function);}

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        return this->m_converged;
      }


      bool all_different( Real a, Real b, Real c, Real d ) const {
        return a != b &&
               a != c &&
               a != d &&
               b != c &&
               b != d &&
               c != d;
      }


      Real eval(FunctionWrapper function) {

    // check for trivial solution
    m_converged = m_fa == 0; if ( m_converged ) return m_a;
    m_converged = m_fb == 0; if ( m_converged ) return m_b;

    // Finds either an exact solution or an approximate solution
    // of the equation f(x)=0 in the interval [a,b].
    //
    // At the beginning of each iteration, the current enclosing interval
    // is recorded as [a0,b0].
    // The first iteration is simply a secant step.
    // Starting with the second iteration, three steps are taken in each iteration.
    // First two steps are either quadratic interpolation
    // or cubic inverse interpolation.
    // The third step is a double-size secant step.
    // If the diameter of the enclosing interval obtained after
    // those three steps is larger than 0.5*(b0-a0),
    // then an additional bisection step will be taken.

    m_e  = QUIET_NAN;
    m_fe = QUIET_NAN; // Dumb values

    //
    // While f(left) or f(right) are infinite perform bisection
    //
    while ( !( std::isfinite(m_fa) && std::isfinite(m_fb) ) ) {
      ++m_iterations;
      m_c  = (m_a+m_b)/2;
      this->evaluate_function(function, m_c, m_fc);
      m_converged = m_fc == 0;
      if ( m_converged ) return m_c;
      if ( m_fa*m_fc < 0 ) {
        // --> [a,c]
        m_b = m_c; m_fb = m_fc;
      } else {
        // --> [c,b]
        m_a = m_c; m_fa = m_fc;
      }
      m_converged = (m_b-m_a) <= m_tolerance;
      if ( std::abs(m_fb) <= std::abs(m_fa) ) {
        this->set_tolerance(m_b);
        if ( m_converged ) return m_b;
      } else {
        this->set_tolerance(m_a);
        if ( m_converged ) return m_a;
      }
    }
    {
      Real ba  = m_b  - m_a;
      Real R   = ba/(m_fb - m_fa);
      m_c = std::abs(m_fb) < std::abs(m_fa) ? m_b+m_fb*R : m_a-m_fa*R;
      // impedisce m_c troppo vicino ad m_a o m_b
      Real delta = m_interval_shink*ba, amin, bmax;
      if      ( m_c < (amin=m_a+delta) ) m_c = amin;
      else if ( m_c > (bmax=m_b-delta) ) m_c = bmax;
    }
    //
    // Call "bracketing" to get a shrinked enclosing interval as
    // well as to update the termination criterion.
    // Stop the procedure if the criterion is satisfied or the
    // exact solution is obtained.
    //
    m_converged = bracketing(function);
    if ( m_converged ) return m_a;
    // Iteration starts.
    // The enclosing interval before executing the iteration is recorded as [a0, b0].
    m_converged = false;

    // ITERATION STARTS. THE ENCLOSING INTERVAL BEFORE EXECUTING THE
    // ITERATION IS RECORDED AS [A0, B0].
    while ( !m_converged ) {
      ++m_iterations;
      Real BA0 = m_b-m_a;

      // Calculates the termination criterion.
      // Stops the procedure if the criterion is satisfied.

      {
        Real abs_fa = std::abs(m_fa);
        Real abs_fb = std::abs(m_fb);
        Real c      = abs_fb <= abs_fa ? m_b : m_a;
        this->set_tolerance( c );
        m_converged = BA0 <= m_tolerance;
        if ( m_converged ) return c;
      }

      //
      // Starting with the second iteration, in the first two steps, either
      // quadratic interpolation is used by calling the subroutine "newtonquadratic"
      // or the cubic inverse interpolation is used by calling the subroutine
      // "pzero". in the following, if "prof" is not equal to 0, then the
      // four function values "fa", "fb", "fd", and "fe" are distinct, and
      // hence "pzero" will be called.
      //

      bool do_newton_quadratic = false;
      if ( !std::isfinite(m_fe) ) {
        do_newton_quadratic = true;
      } else if ( !this->all_different( m_fa, m_fb, m_fd, m_fe ) ) {
        do_newton_quadratic = true;
      } else {
        m_c = this->pzero();
        do_newton_quadratic = (m_c-m_a)*(m_c-m_b) >= 0;
      }
      if ( do_newton_quadratic ) m_c = this->newton_quadratic(2);

      m_e  = m_d;
      m_fe = m_fd;

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      m_converged = this->bracketing(function) || (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return std::abs(m_fa) < std::abs(m_fb) ? m_a : m_b;
      if ( !this->all_different( m_fa, m_fb, m_fd, m_fe ) ) {
        do_newton_quadratic = true;
      } else {
        m_c = this->pzero();
        do_newton_quadratic = (m_c-m_a)*(m_c-m_b) >= 0;
      }
      if ( do_newton_quadratic ) m_c = this->newton_quadratic(3);

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //

      {
        m_converged = this->bracketing(function) || (m_b-m_a) <= m_tolerance;
        Real abs_fa = std::abs(m_fa);
        Real abs_fb = std::abs(m_fb);
        if ( m_converged ) return abs_fa < abs_fb ? m_a : m_b;

        m_e  = m_d;
        m_fe = m_fd;
        // Takes the double-size secant step.
        Real u, fu;
        if ( abs_fa < abs_fb ) { u = m_a; fu = m_fa; }
        else                   { u = m_b; fu = m_fb; }
        Real hba = (m_b-m_a)/2;
        m_c = u-4*(fu/(m_fb-m_fa))*hba;
        if ( std::abs(m_c-u) > hba ) m_c = m_a + hba;
      }

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      m_converged = this->bracketing(function) || (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return std::abs(m_fa) < std::abs(m_fb) ? m_a : m_b;
      //
      // Determines whether an additional bisection step is needed.
      // Takes it if necessary.
      //
      if ( (m_b-m_a) < m_mu*BA0 ) continue;

      m_e  = m_d;
      m_fe = m_fd;
      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      {
        Real ba = m_b-m_a;
        m_c = m_a + ba/2;
        m_converged = this->bracketing(function) || ba <= m_tolerance;
      }
    }
    // TERMINATES THE PROCEDURE AND RETURN THE "ROOT".
    return std::abs(m_fa) < std::abs(m_fb) ? m_a : m_b;
      }

    }; // class Algo748

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_ALGO748_HXX
