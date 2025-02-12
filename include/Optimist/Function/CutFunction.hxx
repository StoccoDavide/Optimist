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

#ifndef OPTIMIST_CUT_FUNCTION_HXX
#define OPTIMIST_CUT_FUNCTION_HXX

namespace Optimist
{

  /**
  * \brief Class container for the scalar function defined on a "cut" of a vector-valued function.
  *
  * \includedoc docs/markdown/Function.md
  *
  * \tparam N The dimension of the vector-valued function.
  * \tparam DerivedFunction Derived scalar function class.
  */
  template <Integer N, typename DerivedFunction>
  class CutFunction : public CostFunction<N, DerivedFunction>
  {
    Integer m_d{-1}; /**< The cut dimension. */

  public:
    /**
    * Class constructor for the function.
    */
    CutFunction() : CostFunction<N, DerivedFunction>() {}

    /**
    * Class constructor for the function.
    */
    CutFunction(const Integer d) : CutFunction() {this->cut(d);}

    /**
    * Set the cut dimension.
    * \param[in] d The cut dimension.
    */
    void cut(const Integer d)
    {
      #define CMD "Optimist::CutFunction::cut(d) "

      OPTIMIST_ASSERT(d >= Integer(0) && d < N, CMD "cut dimension out of bounds.");
      this->m_d = d;

      #undef CMD
    }

    /**
    * Get the cut dimension.
    * \return The cut dimension.
    */
    Integer cut() const {return this->m_d;}

    /**
    * Get the function name.
    * \return The function name.
    */
    std::string name() const
    {
      return this->m_function.name() + "(" +
        this->m_d > Integer(0) ? std::to_string(this->cut()) : std::string("undefined")
        + "-dimension cut)";
    };

    /**
    * Compute the function value at the input point.
    * \param[in] x Input point.
    * \param[out] out The function value.
    */
    void evaluate(Real x, Real & out) const
    {
      this->m_function.evaluate(x, out)(this->cut());
    }

    /**
    * Compute the function first derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function first derivative.
    */
    void first_derivative(Real x, Real & out) const
    {
      this->m_function.first_derivative(x, out)(this->cut());
    }

    /**
    * Compute the function second derivative at the input point.
    * \param[in] x Input point.
    * \param[out] out The function second derivative.
    */
    void second_derivative(Real x, Real & out) const
    {
      this->m_function.second_derivative(x, out)(this->cut());
    }

    /**
    * Retrieve the known solution at the index.
    * \param[in] i The index of the known solution.
    * \return The known solution.
    */
    Real solution(const Integer i) const
    {
      return this->m_function.solution(i)(this->cut());
    }

    /**
    * Retrieve the initial guess at the index.
    * \param[in] i The index of the initial guess.
    * \return The initial guess.
    */
    Real guess(const Integer i) const
    {
      return this->m_function.guess(i)(this->cut());
    }

  }; // class CutFunction

} // namespace Optimist

#endif // OPTIMIST_CUT_FUNCTION_HXX
