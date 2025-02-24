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

#ifndef OPTIMIST_ROOTFINDER_GREENSTADT_HXX
#define OPTIMIST_ROOTFINDER_GREENSTADT_HXX

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |    ____                         _            _ _
     |   / ___|_ __ ___  ___ _ __  ___| |_ __ _  __| | |_
     |  | |  _| '__/ _ \/ _ \ '_ \/ __| __/ _` |/ _` | __|
     |  | |_| | | |  __/  __/ | | \__ \ || (_| | (_| | |_
     |   \____|_|  \___|\___|_| |_|___/\__\__,_|\__,_|\__|
     |
    \*/

    /**
    * \brief Class container for the Greenstadt's method.
    *
    * \includedoc docs/markdown/RootFinder/Greenstadt.md
    *
    * \tparam Real Scalar number type.
    * \tparam N Dimension of the root-finding problem.
    */
    template <typename Real, Integer N>
    class Greenstadt : public QuasiNewton<Real, N, Greenstadt<Real, N>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{false};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      using Method = enum class Method : Integer {ONE = 1, TWO = 2}; /**< Greenstadt solver type. */
      using Vector = typename QuasiNewton<Real, N, Greenstadt<Real, N>>::Vector;
      using Matrix = typename QuasiNewton<Real, N, Greenstadt<Real, N>>::Matrix;
      using FunctionWrapper = typename QuasiNewton<Real, N, Greenstadt<Real, N>>::FunctionWrapper;
      using JacobianWrapper = typename QuasiNewton<Real, N, Greenstadt<Real, N>>::JacobianWrapper;
      using QuasiNewton<Real, N, Greenstadt<Real, N>>::solve;

    private:
      Method m_method{Method::ONE}; /**< Greenstadt solver type. */

    public:
      /**
      * Class constructor for the Greenstadt solver.
      */
      Greenstadt() {}

      /**
      * Get the Greenstadt solver name.
      * \return The Greenstadt solver name.
      */
      std::string name_impl() const
      {
        std::ostringstream os;
        os << "Greenstadt";
        if (this->m_method == Method::ONE) {
          os << "1";
        } else if (this->m_method == Method::TWO) {
          os << "2";
        }
        return os.str();
      }

      /**
      * Get the enumeration type of the Greenstadt solver method.
      * \return The Greenstadt solver enumeration type.
      */
      Method method() const {return this->m_method;}

      /**
      * Set the enumeration type of the Greenstadt solver method.
      * \param[in] t_method The Greenstadt solver method enumeration type.
      */
      void method(Method t_method) {this->m_method = t_method;}

      /**
      * Enable the \em Greenstadt1 solver method.
      */
      void enable_one_method() {this->m_method = Method::ONE;}

      /**
      * Enable the \em Greenstadt2 solver method.
      */
      void enable_two_method() {this->m_method = Method::TWO;}

      /**
      * Set the Greenstadt solver method.
      * \param[in] t_method The Greenstadt solver method enumeration type.
      */
      void set_method(Method t_method) {this->m_method = t_method;}

      /**
      * Jacobian approximation update rule for the Greenstadt's method.
      * \param[in] delta_x_old Old difference between points.
      * \param[in] delta_function_old Old difference between function values.
      * \param[in] jacobian_old Old jacobian approximation.
      * \param[in] delta_x_new New difference between points.
      * \param[in] delta_function_new New difference between function values.
      * \param[in] function_new New function value.
      * \param[out] jacobian_new New jacobian approximation.
      */
      void update_impl(
        Vector const & /*delta_x_old*/, Vector const & /*delta_function_old*/, Matrix const & jacobian_old,
        Vector const & delta_x_new,     Vector const & delta_function_new,     Vector const & function_new,
        Matrix       & jacobian_new
      ) {
        if (this->m_method == Method::ONE) {
          // Greenstadt's 1st method
          // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C = F1;
          jacobian_new = jacobian_old - (jacobian_old*delta_function_new-delta_x_new)/(function_new.transpose()*delta_function_new)*function_new.transpose();
        } else if (this->m_method == Method::TWO) {
          // Greenstadt's 2nd method
          // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C  = J0'*J0*DF1;
          Vector C(jacobian_old.transpose()*jacobian_old*delta_function_new);
          jacobian_new = jacobian_old - (jacobian_old*delta_function_new-delta_x_new)/(C.transpose()*delta_function_new)*C.transpose();
        }
      }

    }; // class Greenstadt

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_GREENSTADT_HXX
