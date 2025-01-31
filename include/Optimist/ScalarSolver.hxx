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

#ifndef OPTIMIST_SCALAR_SOLVER_HXX
#define OPTIMIST_SCALAR_SOLVER_HXX

namespace Optimist
{

  /**
  * Namespace for the solution of scalar problems.
  */
  namespace ScalarSolver
  {

    /*\
     |   ____            _            ____        _
     |  / ___|  ___ __ _| | __ _ _ __/ ___|  ___ | |_   _____ _ __
     |  \___ \ / __/ _` | |/ _` | '__\___ \ / _ \| \ \ / / _ \ '__|
     |   ___) | (_| (_| | | (_| | |   ___) | (_) | |\ V /  __/ |
     |  |____/ \___\__,_|_|\__,_|_|  |____/ \___/|_| \_/ \___|_|
     |
    \*/

    /**
    * \brief Class container for the scalar solver.
    *
    * \includedoc docs/markdown/ScalarSolver.md
    */
    class ScalarSolver : public Solver<1, 1>
    {
    public:
      // Function types
      using Function         = typename Solver<1, 1>::Function;         /**< Function type. */
      using FirstDerivative  = typename Solver<1, 1>::FirstDerivative;  /**<  First derivative type. */
      using SecondDerivative = typename Solver<1, 1>::SecondDerivative; /**< Second derivative type. */

      /**
      * Class constructor for the scalar solver.
      */
      ScalarSolver() : Solver<1, 1>() {}

    }; // class ScalarSolver

  } // namespace ScalarSolver

} // namespace Optimist

#endif // OPTIMIST_SCALAR_SOLVER_HXX
