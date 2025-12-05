/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
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
   * \tparam Input Function input type.
   * \tparam Output Function output type.
   * \tparam DerivedFunction Derived function class.
   */
  template <typename Input, typename Output, typename DerivedFunction>
  class FunctionBase
  {
  public:

    // Input and output types
    using InputTrait  = TypeTrait<Input>;
    using OutputTrait = TypeTrait<Output>;
    using Scalar      = typename InputTrait::Scalar;

    // Input and output must have the same scalar type
    static_assert(std::is_same<Scalar, typename OutputTrait::Scalar>::value,
      "Input and output scalar types must be the same.");

    // If both input and output are eigen types they be both fixed-size, dynamic-size, or sparse
    static_assert(!(InputTrait::IsEigen && OutputTrait::IsEigen) ||
      (InputTrait::IsFixed && OutputTrait::IsFixed) ||
      (InputTrait::IsDynamic && OutputTrait::IsDynamic) ||
      (InputTrait::IsSparse && OutputTrait::IsSparse),
      "Input and output Eigen types must be both fixed-size, dynamic-size, or sparse.");

    // Derivative types
    using FirstDerivative = std::conditional_t<InputTrait::IsEigen || OutputTrait::IsEigen,
      std::conditional_t<InputTrait::IsSparse || OutputTrait::IsSparse,
        Eigen::SparseMatrix<Scalar>,
        Eigen::Matrix<Scalar, OutputTrait::Dimension, InputTrait::Dimension>>,
      Scalar>;
    using SecondDerivative = std::conditional_t<InputTrait::IsEigen || OutputTrait::IsEigen,
      std::conditional_t<(InputTrait::Dimension == 1) || (OutputTrait::Dimension == 1),
        std::conditional_t<InputTrait::IsSparse || OutputTrait::IsSparse,
          Eigen::SparseMatrix<Scalar>,
          Eigen::Matrix<Scalar, OutputTrait::Dimension, InputTrait::Dimension>>,
        std::conditional_t<InputTrait::IsSparse || OutputTrait::IsSparse,
          std::vector<Eigen::SparseMatrix<Scalar>>,
          std::vector<Eigen::Matrix<Scalar, OutputTrait::Dimension, InputTrait::Dimension>>>>,
      Scalar>;

    OPTIMIST_BASIC_CONSTANTS(Scalar)

  protected:
    std::vector<Input> m_solutions; /**< Known solutions used for test purposes. */
    std::vector<Input> m_guesses;   /**< Suggested initial guess used for testing. */

  public:
    /**
     * Class constructor for the function.
     */
    FunctionBase() {}

    /**
     * Get the function name.
     * \return The function name.
     */
    constexpr std::string name() const {return static_cast<const DerivedFunction *>(this)->name();};

    /**
     * Compute the function value at the input point.
     * \param[in] x Input point.
     * \param[out] out The function value.
     * \return The boolean flag for successful evaluation.
     */
    bool evaluate(const Input & x, Output & out) const
    {
      return static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
     * Compute the function first derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function first derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool first_derivative(const Input & x, FirstDerivative & out) const
    {
      return static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
     * Compute the function second derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function second derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool second_derivative(const Input & x, SecondDerivative & out) const
    {
      return static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

    /**
     * Get the input dimension of the function.
     * \return The input dimension of the function.
     */
    constexpr Integer input_dimension() const {return InputTrait::Dimension;}

    /**
     * Get the output dimension of the function.
     * \return The output dimension of the function.
     */
    constexpr Integer output_dimension() const {return OutputTrait::Dimension;}

    /**
     * Get the vector of known solutions.
     * \return The vector of known solutions.
     */
    const std::vector<Input> & solutions() const {return this->m_solutions;}

    /**
     * Get the vector of initial guesses.
     * \return The vector of initial guesses.
     */
    const std::vector<Input> & guesses() const {return this->m_guesses;}

    /**
     * Retrieve the known solution at the index.
     * \param[in] i The index of the known solution.
     * \return The known solution.
     */
    const Input & solution(Integer const i) const {return this->m_solutions.at(i);}

    /**
     * Retrieve the initial guess at the index.
     * \param[in] i The index of the initial guess.
     * \return The initial guess.
     */
    const Input & guess(Integer const i) const {return this->m_guesses.at(i);}

    /**
     * Check if the input point is a known solution.
     * \param[in] x Input point.
     * \param[in] tol Tolerance.
     * \return True if the input point is a known solution, false otherwise.
     */
    bool is_solution(const Input & x, Scalar const tol = EPSILON_LOW) const
    {
      for (const auto & s : this->m_solutions) {
        if constexpr (InputTrait::IsEigen) {
          if((x - s).norm() < tol) {return true;}
        } else {
          if (std::abs(x - s) < tol) {return true;}
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
   * \brief Class container for the vector-valued function (both input and output are vectors).
   *
   * \tparam Input Function input type.
   * \tparam Output Function output type.
   * \tparam DerivedFunction Derived function class.
   */
  template <typename Input, typename Output, typename DerivedFunction>
  class Function : public FunctionBase<Input, Output, DerivedFunction>
  {
  public:
    friend class FunctionBase<Input, Output, DerivedFunction>;

    // Derivative types
    using Matrix = typename FunctionBase<Input, Output, DerivedFunction>::FirstDerivative;
    using Tensor = typename FunctionBase<Input, Output, DerivedFunction>::SecondDerivative;

    /**
     * Class constructor for the vector-valued function.
     */
    Function() {}

    /**
     * Get the function name.
     * \return The function name.
     */
    constexpr std::string name() const {return static_cast<const DerivedFunction *>(this)->name_impl();}

    /**
     * Compute the function value at the input point.
     * \param[in] x Input point.
     * \param[out] out The function value.
     * \return The boolean flag for successful evaluation.
     */
    bool evaluate(const Input & x, Output & out) const
    {
      return static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
     * Compute the function first derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function first derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool jacobian(const Input & x, Matrix & out) const
    {
      return static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
     * Compute the function second derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function second derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool hessian(const Input & x, Tensor & out) const
    {
      return static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class Function

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
   * \brief Class container for the cost function.
   *
   * \tparam Input Function input type.
   * \tparam DerivedFunction Derived cost function class.
   */
  template <typename Input, typename DerivedFunction>
  requires TypeTrait<Input>::IsEigen
  class Function<Input, typename Input::Scalar, DerivedFunction> : public FunctionBase<Input, typename Input::Scalar, DerivedFunction>
  {
  public:
    friend class FunctionBase<Input, typename Input::Scalar, DerivedFunction>;

    // Input and output types
    using Vector = typename FunctionBase<Input, typename Input::Scalar, DerivedFunction>::Input;

    // Derivative types
    using RowVector = typename FunctionBase<Input, typename Input::Scalar, DerivedFunction>::FirstDerivative;
    using Matrix    = typename FunctionBase<Input, typename Input::Scalar, DerivedFunction>::SecondDerivative;
    /**
     * Class constructor for the function.
     */
    Function<Input, typename Input::Scalar, DerivedFunction>() {}
    /**
     * Get the function name.
     * \return The function name.
     */
    constexpr std::string name() const {return static_cast<const DerivedFunction *>(this)->name_impl();}

    /**
     * Compute the function value at the input point.
     * \param[in] x Input point.
     * \param[out] out The function value.
     * \return The boolean flag for successful evaluation.
     */
    bool evaluate(Vector const & x, Vector & out) const
    {
      return static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
     * Compute the function first derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function first derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool gradient(Vector const & x, RowVector & out) const
    {
      return static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
     * Compute the function second derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function second derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool hessian(Vector const & x, Matrix & out) const
    {
      return static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class Function

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /**
  * \brief Class container for the scalar function.
  *
  * \tparam Scalar Floating-point number type.
  * \tparam DerivedFunction Derived scalar function class.
  */
  template <typename Scalar, typename DerivedFunction>
  requires std::is_floating_point_v<Scalar>
  class Function<Scalar, Scalar, DerivedFunction> : public FunctionBase<Scalar, Scalar, DerivedFunction>
  {
  public:
    friend class FunctionBase<Scalar, Scalar, DerivedFunction>;

    /**
     * Class constructor for the function.
     */
    Function<Scalar, Scalar, DerivedFunction>() {}

    /**
     * Get the function name.
     * \return The function name.
     */
    constexpr std::string name() const {return static_cast<const DerivedFunction *>(this)->name_impl();}

    /**
     * Compute the function value at the input point.
     * \param[in] x Input point.
     * \param[out] out The function value.
     * \return The boolean flag for successful evaluation.
     */
    bool evaluate(Scalar x, Scalar & out) const
    {
      return static_cast<const DerivedFunction *>(this)->evaluate_impl(x, out);
    }

    /**
     * Compute the function first derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function first derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool first_derivative(Scalar x, Scalar & out) const
    {
      return static_cast<const DerivedFunction *>(this)->first_derivative_impl(x, out);
    }

    /**
     * Compute the function second derivative at the input point.
     * \param[in] x Input point.
     * \param[out] out The function second derivative.
     * \return The boolean flag for successful evaluation.
     */
    bool second_derivative(Scalar x, Scalar & out) const
    {
      return static_cast<const DerivedFunction *>(this)->second_derivative_impl(x, out);
    }

  }; // class Function

} // namespace Optimist

#endif // OPTIMIST_FUNCTION_HH
