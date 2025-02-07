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

#ifndef OPTIMIST_SCALAR_ROOT_FINDER_HXX
#define OPTIMIST_SCALAR_ROOT_FINDER_HXX

namespace Optimist
{

  /**
  * \brief Namespace for scalar root-finding algorithms.
  */
  namespace ScalarRootFinder
  {

    /*\
     |   ____            _            ____             _   _____ _           _
     |  / ___|  ___ __ _| | __ _ _ __|  _ \ ___   ___ | |_|  ___(_)_ __   __| | ___ _ __
     |  \___ \ / __/ _` | |/ _` | '__| |_) / _ \ / _ \| __| |_  | | '_ \ / _` |/ _ \ '__|
     |   ___) | (_| (_| | | (_| | |  |  _ < (_) | (_) | |_|  _| | | | | | (_| |  __/ |
     |  |____/ \___\__,_|_|\__,_|_|  |_| \_\___/ \___/ \__|_|   |_|_| |_|\__,_|\___|_|
     |
    \*/

    /**
    * \brief Class container for the scalar scalar root-finder.
    *
    * \includedoc docs/markdown/ScalarRootFinder.md
    *
    * \tparam DerivedSolver Derived solver class.
    */
    template <typename DerivedSolver>
    class ScalarRootFinder : public Solver<1, 1, DerivedSolver>
    {
      friend Solver<1, 1, ScalarRootFinder<DerivedSolver>>;

    public:
      // Function types
      using Function         = typename Solver<1, 1, DerivedSolver>::Function;         /**< Function type. */
      using FirstDerivative  = typename Solver<1, 1, DerivedSolver>::FirstDerivative;  /**<  First derivative type. */
      using SecondDerivative = typename Solver<1, 1, DerivedSolver>::SecondDerivative; /**< Second derivative type. */

      /**
      * Class constructor for the scalar root-finder.
      */
      ScalarRootFinder() {}

      /**
      * Get the solver name.
      * \return The solver name.
      */
      std::string name() const {return static_cast<const DerivedSolver *>(this)->name_impl();}

    protected:
      /**
      * Solve the root-finding problem given the function, and without derivatives.
      * \param[in] function Function pointer.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Function function, Real const &x_ini, Real &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(x_ini, function, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first derivative.
      * \param[in] function Function pointer.
      * \param[in] first_derivative First derivative of the function.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Function function, FirstDerivative first_derivative, Real const &x_ini, Real &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(x_ini, function, first_derivative, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first and second derivatives.
      * \param[in] function Function pointer.
      * \param[in] first_derivative First derivative of the function.
      * \param[in] second_derivative Second derivative of the function.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Function function, FirstDerivative first_derivative, SecondDerivative second_derivate,
        Real const &x_ini, Real &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(x_ini, function, first_derivative,
          second_derivate, x_sol);
      }

    }; // class ScalarRootFinder

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_HXX
