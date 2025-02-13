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

#ifndef OPTIMIST_FUNCTION_HXX
#define OPTIMIST_FUNCTION_HXX

namespace Optimist
{

  /**
  * \brief Class container for the generic function.
  *
  * \includedoc docs/markdown/Function.md
  *
  * \tparam FunInDim The function problem input dimension.
  * \tparam FunOutDim The function problem output dimension.
  * \tparam DerivedFunction Derived function class.
  */
  template <Integer FunInDim, Integer FunOutDim, typename DerivedFunction>
  class Function
  {

    // Fancy static assertios (just for fun, don't take it too seriously)
    static_assert(FunInDim > Integer(0) && FunOutDim > Integer(0),
      "Negative-dimensional function? Are you serious?");

  public:
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
    Function() {}

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
        if (x.isApprox(s, tol)) {return true;}
      }
      return false;
    }

  }; // class Function

} // namespace Optimist

#endif // OPTIMIST_FUNCTION_HXX
