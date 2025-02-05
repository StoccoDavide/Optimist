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

#ifndef OPTIMIST_OPTIMIZER_HXX
#define OPTIMIST_OPTIMIZER_HXX

namespace Optimist
{

  /**
  * \brief Namespace for multi-dimensional optimization algorithms.
  */
  namespace Optimizer
  {

    /*\
     |    ___        _   _           _
     |   / _ \ _ __ | |_(_)_ __ ___ (_)_______ _ __
     |  | | | | '_ \| __| | '_ ` _ \| |_  / _ \ '__|
     |  | |_| | |_) | |_| | | | | | | |/ /  __/ |
     |   \___/| .__/ \__|_|_| |_| |_|_/___\___|_|
     |        |_|
    \*/

    /**
    * \brief Class container for the multi-dimensional optimizer.
    *
    * \includedoc docs/markdown/Optimizer.md
    *
    * \tparam N The dimension of the optimization problem.
    */
    template <Integer N>
    class Optimizer : public Solver<N, 1>
    {
      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != Integer(0),
        "Have you ever heard of a zero-dimensional optimization problem?");
      static_assert(N != Integer(1),
        "Good try, but you should use a scalar solver for a one-dimensional optimization problem.");

    public:
      // I/O types
      using Vector = typename Solver<N, 1>::InputType; /**< Vector type. */

      // Derivative types
      using RowVector = typename Solver<N, 1>::FirstDerivativeType;  /**< Gradient (row) vector type. */
      using Matrix    = typename Solver<N, 1>::SecondDerivativeType; /**< Hessian matrix type. */

      // Function types
      using Function = typename Solver<N, 1>::Function;         /**< Function type. */
      using Gradient = typename Solver<N, 1>::FirstDerivative;  /**< Gradient function type. */
      using Hessian  = typename Solver<N, 1>::SecondDerivative; /**< Hessian function type. */

      /**
      * Class constructor for the multi-dimensional optimizer.
      */
      Optimizer() {}

      /**
      * Get the number of function first derivative evaluations.
      * \return The number of function first derivative evaluations.
      */
      Integer gradient_evaluations() const {return this->first_derivative_evaluations();}

      /**
      * Get the number of maximum allowed gradient evaluations.
      * \return The number of maximum allowed gradient evaluations.
      */
      Integer max_gradient_evaluations() const {return this->max_first_derivative_evaluations();}

      /**
      * Set the number of maximum allowed gradient evaluations.
      * \param[in] t_gradient_evaluations The number of maximum allowed gradient evaluations.
      */
      void max_gradient_evaluations(Integer t_gradient_evaluations)
      {
        this->first_derivative_evaluations(t_gradient_evaluations);
      }

      /**
      * Get the number of function first derivative evaluations.
      * \return The number of function first derivative evaluations.
      */
      Integer hessian_evaluations() const {return this->first_derivative_evaluations();}

      /**
      * Get the number of maximum allowed hessian evaluations.
      * \return The number of maximum allowed hessian evaluations.
      */
      Integer max_hessian_evaluations() const {return this->max_first_derivative_evaluations();}

      /**
      * Set the number of maximum allowed hessian evaluations.
      * \param[in] t_hessian_evaluations The number of maximum allowed hessian evaluations.
      */
      void max_hessian_evaluations(Integer t_hessian_evaluations)
      {
        this->first_derivative_evaluations(t_hessian_evaluations);
      }

    protected:
      /**
      * Evaluate the gradient.
      * \param[in] x Input point.
      * \param[out] gradient gradient value.
      */
      void evaluate_gradient(const Vector & x, RowVector & gradient)
      {
        this->evaluate_first_derivative(x, gradient);
      }

      /**
      * Evaluate the hessian.
      * \param[in] x Input point.
      * \param[out] hessian hessian value.
      */
      void evaluate_hessian(const Vector & x, Matrix & hessian)
      {
        this->evaluate_first_derivative(x, hessian);
      }

    }; // class Optimizer

  } // namespace Optimizer

} // namespace Optimist

#endif // OPTIMIST_OPTIMIZER_HXX
