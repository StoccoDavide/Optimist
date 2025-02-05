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

#ifndef OPTIMIST_SCALAR_OPTIMIZER_HXX
#define OPTIMIST_SCALAR_OPTIMIZER_HXX

namespace Optimist
{

  /**
  * \brief Namespace for scalar optimization algorithms.
  */
  namespace ScalarOptimizer
  {

    /*\
     |   ____            _             ___        _   _           _
     |  / ___|  ___ __ _| | __ _ _ __ / _ \ _ __ | |_(_)_ __ ___ (_)_______ _ __
     |  \___ \ / __/ _` | |/ _` | '__| | | | '_ \| __| | '_ ` _ \| |_  / _ \ '__|
     |   ___) | (_| (_| | | (_| | |  | |_| | |_) | |_| | | | | | | |/ /  __/ |
     |  |____/ \___\__,_|_|\__,_|_|   \___/| .__/ \__|_|_| |_| |_|_/___\___|_|
     |                                     |_|
    \*/

    /**
    * \brief Class container for the scalar optimizer.
    *
    * \includedoc docs/markdown/ScalarOptimizer.md
    */
    class ScalarOptimizer : public Solver<1, 1>
    {
    public:
      // Function types
      using Function         = typename Solver<1, 1>::Function;         /**< Function type. */
      using FirstDerivative  = typename Solver<1, 1>::FirstDerivative;  /**<  First derivative type. */
      using SecondDerivative = typename Solver<1, 1>::SecondDerivative; /**< Second derivative type. */

      /**
      * Class constructor for the scalar optimizer.
      */
      ScalarOptimizer() {}

    }; // class ScalarOptimizer

  } // namespace ScalarOptimizer

} // namespace Optimist

#endif // OPTIMIST_SCALAR_OPTIMIZER_HXX
