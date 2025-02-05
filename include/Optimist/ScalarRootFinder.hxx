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
    */
    class ScalarRootFinder : public Solver<1, 1>
    {
    public:
      // Function types
      using Function         = typename Solver<1, 1>::Function;         /**< Function type. */
      using FirstDerivative  = typename Solver<1, 1>::FirstDerivative;  /**<  First derivative type. */
      using SecondDerivative = typename Solver<1, 1>::SecondDerivative; /**< Second derivative type. */

      /**
      * Class constructor for the scalar root-finder.
      */
      ScalarRootFinder() {}

    protected:
      /**
      * Damp the step using the affine invariant criterion.
      * \param[in] x_old Old point.
      * \param[in] function_old Old function value.
      * \param[in] step_old Old step.
      * \param[out] x_new New point.
      * \param[out] function_new New function value.
      * \param[out] step_new New step.
      * \return The damping boolean flag, true if the damping is successful, false otherwise.
      */
      bool damp(Real const & x_old, Real const & function_old, Real const & step_old,
        Real & x_new, Real & function_new, Real & step_new)
      {
        Real step_norm_old, step_norm_new, residuals_old, residuals_new, tau{1.0};
        for (this->m_relaxations = Integer(0); this->m_relaxations < this->m_max_relaxations; ++this->m_relaxations)
        {
          // Update point
          step_new = tau * step_old;
          x_new = x_old + step_new;
          this->evaluate_function(x_new, function_new);

          // Check relaxation
          residuals_old = std::abs(function_old);
          residuals_new = std::abs(function_new);
          step_norm_old = std::abs(step_old);
          step_norm_new = std::abs(step_new);
          if (residuals_new < residuals_old || step_norm_new < (Real(1.0)-tau/Real(2.0))*step_norm_old) {
            return true;
          } else {
            tau *= this->m_alpha;
          }
        }
        return false;
      }

    }; // class ScalarRootFinder

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_HXX
