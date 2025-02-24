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

#ifndef OPTIMIST_COST_FUNCTION_HXX
#define OPTIMIST_COST_FUNCTION_HXX

namespace Optimist
{

  /**
  * \brief Class container for the cost function.
  *
  * \includedoc docs/markdown/Function.md
  *
  * \tparam N The dimension of the cost function input.
  * \tparam DerivedFunction Derived cost function class.
  */
  template <typename Real, Integer N, typename DerivedFunction>
  class CostFunction : public Function<Real, N, 1, DerivedFunction>
  {
  public:
    friend Function<Real, N, 1, CostFunction<Real, N, DerivedFunction>>;

    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(N != Integer(0),
      "Are you sure you want to a zero-dimensional system of equations?");
    static_assert(N != Integer(1),
      "C'mon, let's not kid ourselves. Use a scalar function...");

    // I/O types
    using Vector = typename Function<Real, N, 1, DerivedFunction>::InputType; /**< Vector type. */

    // Derivative types
    using RowVector = typename Function<Real, N, 1, DerivedFunction>::FirstDerivativeType;  /**< Gradient (row) vector type. */
    using Matrix    = typename Function<Real, N, 1, DerivedFunction>::SecondDerivativeType; /**< Hessian matrix type. */

    /**
    * Class constructor for the function.
    */
    CostFunction() {}

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
    void evaluate(const Vector & x, Vector & out) const
    {
      static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
    * Compute the function first derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function first derivative.
    */
    void gradient(const Vector & x, RowVector & out) const
    {
      static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
    * Compute the function second derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function second derivative.
    */
    void hessian(const Vector & x, Matrix & out) const
    {
      static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class CostFunction

} // namespace Optimist

#endif // OPTIMIST_COST_FUNCTION_HXX
