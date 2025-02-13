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
    public:
      friend Solver<1, 1, ScalarRootFinder<DerivedSolver>>;

      static constexpr bool requires_function          = DerivedSolver::requires_function;
      static constexpr bool requires_first_derivative  = DerivedSolver::requires_first_derivative;
      static constexpr bool requires_second_derivative = DerivedSolver::requires_second_derivative;

      using FunctionWrapper         = typename Solver<1, 1, DerivedSolver>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename Solver<1, 1, DerivedSolver>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename Solver<1, 1, DerivedSolver>::SecondDerivativeWrapper;

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
      * \param[in] function Function wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, Real x_ini, Real &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(x_ini, function, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first derivative.
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, FirstDerivativeWrapper first_derivative, Real x_ini,
        Real &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(x_ini, function, first_derivative, x_sol);
      }

      /**
      * Solve the root-finding problem given the function, and its first and second derivatives.
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] second_derivative Second derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve(FunctionWrapper function, FirstDerivativeWrapper first_derivative, SecondDerivativeWrapper
        second_derivate, Real x_ini, Real &x_sol)
      {
        return static_cast<DerivedSolver *>(this)->solve_impl(x_ini, function, first_derivative,
          second_derivate, x_sol);
      }

    }; // class ScalarRootFinder

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_HXX
