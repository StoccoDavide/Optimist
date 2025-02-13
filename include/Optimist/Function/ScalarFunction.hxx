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

#ifndef OPTIMIST_SCALAR_FUNCTION_HXX
#define OPTIMIST_SCALAR_FUNCTION_HXX

namespace Optimist
{

  /**
  * \brief Class container for the scalar function.
  *
  * \includedoc docs/markdown/Function.md
  *
  * \tparam DerivedFunction Derived scalar function class.
  */
  template <typename DerivedFunction>
  class ScalarFunction : public Function<1, 1, DerivedFunction>
  {
    friend Function<1, 1, ScalarFunction<DerivedFunction>>;

  public:
    /**
    * Class constructor for the function.
    */
    ScalarFunction() {}

    /**
    * Get the function name.
    * \return The function name.
    */
    std::string name() const {return static_cast<const DerivedFunction *>(this)->name_impl();}

    /**
    * Compute the function value at the input point.
    * \param[in] x Input point.
    * \param[out] out The function value.
    */
    void evaluate(Real x, Real & out) const
    {
      static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
    * Compute the function first derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function first derivative.
    */
    void first_derivative(Real x, Real & out) const
    {
      static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
    * Compute the function second derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function second derivative.
    */
    void second_derivative(Real x, Real & out) const
    {
      static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class ScalarFunction

} // namespace Optimist

#endif // OPTIMIST_SCALAR_FUNCTION_HXX
