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

#ifndef OPTIMIST_COS_HXX
#define OPTIMIST_COS_HXX

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the cosine function.
    *
    * Class container for the cosine function, which is defined as:
    * \f[
    * f(x) = \cos(x) \text{.}
    * \f]
    * The function has roots at \f$x = i\pi\f$, with \f$f(x) = 0\f$, and minima at \f$x = -\pi + 2i\pi\f$,
    * with \f$f(x) = -1\f$ and \f$i = 0, 1, \ldots, n\f$. The initial guesses are generated on the
    * range \f$x \in \left[-\pi, \pi\right]\f$.
    */
    class Cos : public ScalarFunction<Cos>
    {
    public:
      /**
      * Class constructor for the cosine function.
      */
      Cos()
      {
        this->m_solutions.emplace_back(PI);
        this->m_solutions.emplace_back(0.0);
        this->m_solutions.emplace_back(-PI);
        for (Real x{-PI}; x < PI + EPSILON; x += PI/2.0) {this->m_guesses.emplace_back(x);}
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Cos";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = std::cos(x);}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real x, Real & out) const {out = -std::sin(x);}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real x, Real & out) const {out = -std::cos(x);}

    }; // class Cos

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_COS_HXX
