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
  * \tparam N The input dimension of the vector-valued function.
  * \tparam M The output dimension of the vector-valued function.
  * \tparam DerivedFunction Derived vector-valued function class.
  */
  template <typename Real, Integer N, Integer M, typename DerivedFunction>
  class VectorFunction : public Function<Real, N, M, DerivedFunction>
  {
  public:
    friend Function<Real, N, M, VectorFunction<Real, N, M, DerivedFunction>>;

    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(N != Integer(0) && M != Integer(0),
      "Are you sure you want to a zero-dimensional system of equations?");
    static_assert(N != Integer(1) && M != Integer(1),
      "C'mon, let's not kid ourselves. Use a scalar function...");
    static_assert(M != Integer(1),
      "Good try, but you're looking for an objective function, not a vector-valued function.");

    // I/O types
    using InputVector = typename Function<Real, N, M, DerivedFunction>::InputType; /**< Input vector type. */
    using OutputVector = typename Function<Real, N, M, DerivedFunction>::OutputType; /**< Output vector type. */

    // Derivative types
    using Matrix = typename Function<Real, N, M, DerivedFunction>::FirstDerivativeType;  /**< Jacobian matrix type. */
    using Tensor = typename Function<Real, N, M, DerivedFunction>::SecondDerivativeType; /**< Hessian tensor type. */

    /**
    * Class constructor for the vector-valued function.
    */
    VectorFunction() {}

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
    void evaluate(const InputVector & x, OutputVector & out) const
    {
      static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
    * Compute the function first derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function first derivative.
    */
    void jacobian(const InputVector & x, Matrix & out) const
    {
      static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
    * Compute the function second derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function second derivative.
    */
    void hessian(const InputVector & x, Tensor & out) const
    {
      static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class VectorFunction

} // namespace Optimist

#endif // OPTIMIST_VECTOR_FUNCTION_HXX
