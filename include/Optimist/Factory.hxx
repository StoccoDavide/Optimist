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

#ifndef OPTIMIST_FACTORY_HXX
#define OPTIMIST_FACTORY_HXX

namespace Optimist
{
  /**
  * Root-finding with jacobian function signature.
  * \tparam N The dimension of the nonlinear system of equations.
  */
  template <Integer N>
  using RootfindFunction = std::function<bool(
      typename RootFinder::RootFinder<N>::Function,
      typename RootFinder::RootFinder<N>::Jacobian,
      const typename RootFinder::RootFinder<N>::Vector&,
      typename RootFinder::RootFinder<N>::Vector&)>;


  /**
  * Newton's method wrapper.
  * \param[in] t_function The function.
  * \param[in] t_jacobian The Jacobian function.
  * \param[in] x_ini The initialization point.
  * \param[in] x_sol The solution point.
  * \tparam N The dimension of the nonlinear system of equations.
  */
  template <Integer N>
  bool NewtonMethod(
      typename RootFinder::RootFinder<N>::Function t_function,
      typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
      const typename RootFinder::RootFinder<N>::Vector& x_ini,
      typename RootFinder::RootFinder<N>::Vector& x_sol)
  {
    RootFinder::Newton<N> solver;
    return solver.solve(t_function, t_jacobian, x_ini, x_sol);
  }

  /**
  * Broyden's good method wrapper.
  * \param[in] t_function The function.
  * \param[in] t_jacobian The Jacobian function.
  * \param[in] x_ini The initialization point.
  * \param[in] x_sol The solution point.
  */
  template <Integer N>
  bool BroydenGoodMethod(
      typename RootFinder::RootFinder<N>::Function t_function,
      typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
      const typename RootFinder::RootFinder<N>::Vector& x_ini,
      typename RootFinder::RootFinder<N>::Vector& x_sol)
  {
    RootFinder::Broyden<N> solver;
    solver.enable_good_method();
    return solver.solve(t_function, t_jacobian, x_ini, x_sol);
  }

  /**
  * Broyden's bad method wrapper.
  * \param[in] t_function The function.
  * \param[in] t_jacobian The Jacobian function.
  * \param[in] x_ini The initialization point.
  * \param[in] x_sol The solution point.
  */
  template <Integer N>
  bool BroydenBadMethod(
      typename RootFinder::RootFinder<N>::Function t_function,
      typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
      const typename RootFinder::RootFinder<N>::Vector& x_ini,
      typename RootFinder::RootFinder<N>::Vector& x_sol)
  {
    RootFinder::Broyden<N> solver;
    solver.enable_bad_method();
    return solver.solve(t_function, t_jacobian, x_ini, x_sol);
  }

  /**
  * Broyden's combined method wrapper.
  * \param[in] t_function The function.
  * \param[in] t_jacobian The Jacobian function.
  * \param[in] x_ini The initialization point.
  * \param[in] x_sol The solution point.
  */
  template <Integer N>
  bool BroydenCombinedMethod(
      typename RootFinder::RootFinder<N>::Function t_function,
      typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
      const typename RootFinder::RootFinder<N>::Vector& x_ini,
      typename RootFinder::RootFinder<N>::Vector& x_sol)
  {
    RootFinder::Broyden<N> solver;
    solver.enable_combined_method();
    return solver.solve(t_function, t_jacobian, x_ini, x_sol);
  }

  /**
  * Greenstadt1's method wrapper.
  * \param[in] t_function The function.
  * \param[in] t_jacobian The Jacobian function.
  * \param[in] x_ini The initialization point.
  * \param[in] x_sol The solution point.
  */
  template <Integer N>
  bool Greenstadt1Method(
      typename RootFinder::RootFinder<N>::Function t_function,
      typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
      const typename RootFinder::RootFinder<N>::Vector& x_ini,
      typename RootFinder::RootFinder<N>::Vector& x_sol)
  {
    RootFinder::Greenstadt<N> solver;
    solver.enable_one_method();
    return solver.solve(t_function, t_jacobian, x_ini, x_sol);
  }

  /**
  * Greenstadt2's method wrapper.
  * \param[in] t_function The function.
  * \param[in] t_jacobian The Jacobian function.
  * \param[in] x_ini The initialization point.
  * \param[in] x_sol The solution point.
  */
  template <Integer N>
  bool Greenstadt2Method(
      typename RootFinder::RootFinder<N>::Function t_function,
      typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
      const typename RootFinder::RootFinder<N>::Vector& x_ini,
      typename RootFinder::RootFinder<N>::Vector& x_sol)
  {
    RootFinder::Greenstadt<N> solver;
    solver.enable_two_method();
    return solver.solve(t_function, t_jacobian, x_ini, x_sol);
  }

  /**
  * Factory function to create the root-finding solvers map.
  * \tparam N The dimension of the nonlinear system of equations.
  */
  template <Integer N>
  std::unordered_map<std::string, RootfindFunction<N>>& get_rootfind_map() {
    static std::unordered_map<std::string, RootfindFunction<N>> rootfind_map = {
      {"Newton",          &NewtonMethod<N>},
      {"BroydenGood",     &BroydenGoodMethod<N>},
      {"BroydenBad",      &BroydenBadMethod<N>},
      {"BroydenCombined", &BroydenCombinedMethod<N>},
      {"Greenstadt1",     &Greenstadt1Method<N>},
      {"Greenstadt2",     &Greenstadt2Method<N>}
    };
    return rootfind_map;
  }

  /**
  * Root-finding function.
  */
  template <Integer N>
  bool Rootfind(
    const std::string & solver_name,
    typename RootFinder::RootFinder<N>::Function t_function,
    typename RootFinder::RootFinder<N>::Jacobian t_jacobian,
    const typename RootFinder::RootFinder<N>::Vector & x_ini,
    typename RootFinder::RootFinder<N>::Vector & x_sol
  ) {
    #define CMD "Optimist::Rootfind<" << N << ">(...): "

    // Get the root-finding map
    auto & rootfind_map = get_rootfind_map<N>();

    auto it = rootfind_map.find(solver_name);
    if (it == rootfind_map.end()) {
      OPTIMIST_ERROR(CMD << N << "solver '" << solver_name << "' not found.");
      return false;
    }

    // Call the solver function from the map
    return it->second(t_function, t_jacobian, x_ini, x_sol);

    #undef CMD
  }

} // namespace Optimist

#endif // OPTIMIST_FACTORY_HXX
