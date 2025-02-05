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

#ifndef OPTIMIST_ROOT_FINDER_HXX
#define OPTIMIST_ROOT_FINDER_HXX

namespace Optimist
{

  /**
  * \brief Namespace for multi-dimensional root-finding algorithms.
  */
  namespace RootFinder
  {

    /*\
     |   ____             _   _____ _           _
     |  |  _ \ ___   ___ | |_|  ___(_)_ __   __| | ___ _ __
     |  | |_) / _ \ / _ \| __| |_  | | '_ \ / _` |/ _ \ '__|
     |  |  _ < (_) | (_) | |_|  _| | | | | | (_| |  __/ |
     |  |_| \_\___/ \___/ \__|_|   |_|_| |_|\__,_|\___|_|
     |
    \*/

    /**
    * \brief Class container for the multi-dimensional root finder.
    *
    * \includedoc docs/markdown/RootFinder.md
    *
    * \tparam N The dimension of the root-finding problem.
    */
    template <Integer N>
    class RootFinder : public Solver<N, N>
    {
      // Fancy static assertions (just for fun, don't take it too seriously)
      static_assert(N != Integer(0),
        "Are you sure you want to solve a zero-dimensional system of equations?");
      static_assert(N != Integer(1),
        "C'mon, let's not kid ourselves. Use a scalar solver...");

    public:
      using Solver<N, N>::Solver;

      // I/O types
      using Vector = typename Solver<N, N>::InputType; /**< Vector type. */

      // Derivative types
      using Matrix = typename Solver<N, N>::FirstDerivativeType; /**< Jacobian matrix type. */

      // Function types
      using Function = typename Solver<N, N>::Function;        /**< Function type. */
      using Jacobian = typename Solver<N, N>::FirstDerivative; /**< Jacobian function type. */

      /**
      * Class constructor for the multi-dimensional root finder.
      */
      RootFinder() {}

      /**
      * Get the number of function first derivative evaluations.
      * \return The number of function first derivative evaluations.
      */
      Integer jacobian_evaluations() const {return this->first_derivative_evaluations();}

      /**
      * Get the number of maximum allowed Jacobian evaluations.
      * \return The number of maximum allowed Jacobian evaluations.
      */
      Integer max_jacobian_evaluations() const {return this->max_first_derivative_evaluations();}

      /**
      * Set the number of maximum allowed Jacobian evaluations.
      * \param[in] t_jacobian_evaluations The number of maximum allowed Jacobian evaluations.
      */
      void max_jacobian_evaluations(Integer t_jacobian_evaluations)
      {
        this->first_derivative_evaluations(t_jacobian_evaluations);
      }

    protected:
      /**
      * Evaluate the Jacobian function.
      * \param[in] x Input point.
      * \param[out] jacobian Jacobian value.
      */
      void evaluate_jacobian(const Vector & x, Matrix & jacobian)
      {
        this->evaluate_first_derivative(x, jacobian);
      }

    }; // class RootFinder

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOT_FINDER_HXX
