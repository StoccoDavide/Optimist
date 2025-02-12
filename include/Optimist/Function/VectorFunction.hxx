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

#ifndef OPTIMIST_VECTOR_FUNCTION_HXX
#define OPTIMIST_VECTOR_FUNCTION_HXX

namespace Optimist
{

  /**
  * \brief Class container for the vector-valued function.
  *
  * \includedoc docs/markdown/Function.md
  *
  * \tparam N The dimension of the vector-valued function.
  * \tparam DerivedFunction Derived vector-valued function class.
  */
  template <Integer N, typename DerivedFunction>
  class VectorFunction : public Function<N, N, DerivedFunction>
  {
    friend Function<N, N, VectorFunction<N, DerivedFunction>>;

    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(N != Integer(0),
      "Are you sure you want to a zero-dimensional system of equations?");
    static_assert(N != Integer(1),
      "C'mon, let's not kid ourselves. Use a scalar function...");

  public:
    // I/O types
    using Vector = typename Function<N, N, DerivedFunction>::InputType; /**< Vector type. */

    // Derivative types
    using Matrix = typename Function<N, N, DerivedFunction>::FirstDerivativeType;  /**< Jacobian matrix type. */
    using Tensor = typename Function<N, N, DerivedFunction>::SecondDerivativeType; /**< Hessian tensor type. */

    /**
    * Class constructor for the function.
    * \param[in] solutions Number of known solutions.
    * \param[in] guesses Number of initial guesses.
    */
    VectorFunction(const Integer solutions, const Integer guesses)
    : Function<N, N, DerivedFunction>(solutions, guesses) {}

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
    void jacobian(const Vector & x, Matrix & out) const
    {
      static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
    * Compute the function second derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function second derivative.
    */
    void hessian(const Vector & x, Tensor & out) const
    {
      static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class VectorFunction

} // namespace Optimist

#endif // OPTIMIST_VECTOR_FUNCTION_HXX
