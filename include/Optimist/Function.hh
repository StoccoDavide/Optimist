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

#ifndef OPTIMIST_FUNCTION_HH
#define OPTIMIST_FUNCTION_HH

#include "Optimist.hh"

namespace Optimist
{

  /*\
   |   _____                 _   _             ____
   |  |  ___|   _ _ __   ___| |_(_) ___  _ __ | __ )  __ _ ___  ___
   |  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \|  _ \ / _` / __|/ _ \
   |  |  _|| |_| | | | | (__| |_| | (_) | | | | |_) | (_| \__ \  __/
   |  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|____/ \__,_|___/\___|
   |
  \*/

  /**
  * \brief Class container for the generic function.
  *
  * \includedoc docs/markdown/FunctionBase.md
  *
  * \tparam Real Scalar number type.
  * \tparam FunInDim The function problem input dimension.
  * \tparam FunOutDim The function problem output dimension.
  * \tparam DerivedFunction Derived function class.
  */
  template <typename Real, Integer FunInDim, Integer FunOutDim, typename DerivedFunction>
  class FunctionBase
  {
  public:
    // Fancy static assertios (just for fun, don't take it too seriously)
    static_assert(FunInDim > static_cast<Integer>(0) && FunOutDim > static_cast<Integer>(0),
      "Negative-dimensional function? Are you serious?");

    OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

    // I/O types
    using InputType  = typename std::conditional_t<FunInDim == 1,  Real, Eigen::Vector<Real, FunInDim>>;  /**< Input type. */
    using OutputType = typename std::conditional_t<FunOutDim == 1, Real, Eigen::Vector<Real, FunOutDim>>; /**< Output type. */

    // Derivative types
    using FirstDerivativeType  = std::conditional_t<FunInDim == 1 && FunOutDim == 1, Real, Eigen::Matrix<Real, FunOutDim, FunInDim>>; /**< First derivative type. */
    using SecondDerivativeType = std::conditional_t<FunInDim == 1 && FunOutDim == 1, Real,
      std::conditional_t<FunInDim == 1 || FunOutDim == 1, Eigen::Matrix<Real, FunInDim, FunInDim>,
      std::vector<Eigen::Matrix<Real, FunInDim, FunInDim>>>>;  /**< Second derivative type. */

  protected:
    std::vector<InputType> m_solutions; /** Known solutions used for test purposes. */
    std::vector<InputType> m_guesses;   /** Suggested initial guess used for testing. */

  public:
    /**
    * Class constructor for the function.
    */
    FunctionBase() {}

    /**
    * Get the function name.
    * \return The function name.
    */
    std::string name() const {return static_cast<const DerivedFunction *>(this)->name();};

    /**
    * Compute the function value at the input point.
    * \param[in] x Input point.
    * \param[out] out The function value.
    */
    void evaluate(const InputType & x, OutputType & out) const
    {
      static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
    * Compute the function first derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function first derivative.
    */
    void first_derivative(const InputType & x, FirstDerivativeType & out) const
    {
      static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
    * Compute the function second derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function second derivative.
    */
    void second_derivative(const InputType & x, SecondDerivativeType & out) const
    {
      static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

    /**
    * Get the input dimension of the function.
    * \return The input dimension of the function.
    */
    constexpr Integer input_dimension() const {return FunInDim;}

    /**
    * Get the output dimension of the function.
    * \return The output dimension of the function.
    */
    constexpr Integer output_dimension() const {return FunOutDim;}

    /**
    * Get the vector of known solutions.
    * \return The vector of known solutions.
    */
    const std::vector<InputType> & solutions() const {return this->m_solutions;}

    /**
    * Get the vector of initial guesses.
    * \return The vector of initial guesses.
    */
    const std::vector<InputType> & guesses() const {return this->m_guesses;}

    /**
    * Retrieve the known solution at the index.
    * \param[in] i The index of the known solution.
    * \return The known solution.
    */
    const InputType & solution(const Integer i) const {return this->m_solutions.at(i);}

    /**
    * Retrieve the initial guess at the index.
    * \param[in] i The index of the initial guess.
    * \return The initial guess.
    */
    const InputType & guess(const Integer i) const {return this->m_guesses.at(i);}

    /**
    * Check if the input point is a known solution.
    * \param[in] x Input point.
    * \param[in] tol Tolerance.
    * \return True if the input point is a known solution, false otherwise.
    */
    bool is_solution(const InputType & x, const Real tol = EPSILON_LOW) const
    {
      for (const auto & s : this->m_solutions) {
        if constexpr (FunInDim == 1) {
          if (std::abs(x - s) < tol) {return true;}
        } else if constexpr (FunInDim > 1) {
          if((x - s).norm() < tol) {return true;}
        } else {
          OPTIMIST_ERROR("Optimist::FunctionBase::is_solution(...): invalid input dimension.");
          return false;
        }
      }
      return false;
    }

  }; // class FunctionBase

  /*\
   |   _____                 _   _
   |  |  ___|   _ _ __   ___| |_(_) ___  _ __
   |  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \
   |  |  _|| |_| | | | | (__| |_| | (_) | | | |
   |  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
   |
  \*/

  /**
  * \brief Class container for the vector-valued function.
  *
  * \tparam N The input dimension of the vector-valued function.
  * \tparam M The output dimension of the vector-valued function.
  * \tparam DerivedFunction Derived vector-valued function class.
  */
  template <typename Real, Integer N, Integer M, typename DerivedFunction>
  class Function : public FunctionBase<Real, N, M, DerivedFunction>
  {
  public:
    friend class FunctionBase<Real, N, M, Function<Real, N, M, DerivedFunction>>;

    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(N != static_cast<Integer>(0) && M != static_cast<Integer>(0),
      "Are you sure you want to a zero-dimensional system of equations?");

    // I/O types
    using InputVector = typename FunctionBase<Real, N, M, DerivedFunction>::InputType; /**< Input vector type. */
    using OutputVector = typename FunctionBase<Real, N, M, DerivedFunction>::OutputType; /**< Output vector type. */

    // Derivative types
    using Matrix = typename FunctionBase<Real, N, M, DerivedFunction>::FirstDerivativeType;  /**< Jacobian matrix type. */
    using Tensor = typename FunctionBase<Real, N, M, DerivedFunction>::SecondDerivativeType; /**< Hessian tensor type. */

    /**
    * Class constructor for the vector-valued function.
    */
    Function() {}

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

  }; // class Function

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
  * \brief Class container for the cost function.
  *
  * \tparam N The dimension of the cost function input.
  * \tparam DerivedFunction Derived cost function class.
  */
  template <typename Real, Integer N, typename DerivedFunction>
  class Function<Real, N, 1, DerivedFunction> : public FunctionBase<Real, N, 1, DerivedFunction>
  {
  public:
    friend class FunctionBase<Real, N, 1, Function<Real, N, 1, DerivedFunction>>;

    // Fancy static assertions (just for fun, don't take it too seriously)
    static_assert(N != static_cast<Integer>(0),
      "Are you sure you want to a zero-dimensional system of equations?");

    // I/O types
    using Vector = typename FunctionBase<Real, N, 1, DerivedFunction>::InputType; /**< Vector type. */

    // Derivative types
    using RowVector = typename FunctionBase<Real, N, 1, DerivedFunction>::FirstDerivativeType;  /**< Gradient (row) vector type. */
    using Matrix    = typename FunctionBase<Real, N, 1, DerivedFunction>::SecondDerivativeType; /**< Hessian matrix type. */

    /**
    * Class constructor for the function.
    */
    Function<Real, N, 1, DerivedFunction>() {}

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

  }; // class Function

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
  * \brief Class container for the scalar function.
  *
  * \tparam Real Scalar number type.
  * \tparam DerivedFunction Derived scalar function class.
  */
  template <typename Real, typename DerivedFunction>
  class Function<Real, 1, 1, DerivedFunction> : public FunctionBase<Real, 1, 1, DerivedFunction>
  {
  public:
    friend class FunctionBase<Real, 1, 1, Function<Real, 1, 1, DerivedFunction>>;

    /**
    * Class constructor for the function.
    */
    Function<Real, 1, 1, DerivedFunction>() {}

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

  }; // class Function

} // namespace Optimist

#endif // OPTIMIST_FUNCTION_HH
