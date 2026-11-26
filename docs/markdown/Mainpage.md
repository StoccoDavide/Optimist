# Optimist

`Optimist` is a library for the numerical solution of optimization problems and zero-finding problems. It is based on the [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) library and written in C++20.

Are you looking for the online documentation? Visit [this link](https://stoccodavide.github.io/Optimist/)!

## Installation

### Quick and dirty

`Optimist` is a header-only library that depends only on [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) (version >= 5.0.1), so the quick and dirty way of installing it is by simply copying the `include` directory to your project and make sure to have [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) available however you see fit. Alternatively, you can do things properly and use `CMake` (version >= 3.14).

### CMake

If you are using CMake, you can add the library as a subdirectory in your project.

```cmake
add_subdirectory(path/to/Optimist)
target_link_libraries(your_target PRIVATE Optimist::Optimist)
```

If you already have `Optimist` somewhere on your system, you can use `find_package` directly.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  Optimist
  ${YOUR_DESIRED_OPTIMIST_VERSION}
  NO_MODULE
)

target_link_libraries(your_target PRIVATE Optimist::Optimist)
```

## Solvers

All solvers work with scalars or with fixed-size, dynamic-size and sparse `Eigen` vectors, in single or double precision.

| Root finders | Kind |
|---|---|
| `NewtonRaphson`, `Halley`, `Chebyshev`, `Varona` | Scalar, open (derivatives) |
| `Algo748`, `Chandrupatla`, `ITP` | Scalar, bracketing |
| `Newton`, `Broyden`, `Greenstadt` | Vector, (quasi-)Newton |
| `LevenbergMarquardt` | Vector, damped Gauss-Newton (also overdetermined) |

| Optimizers | Kind |
|---|---|
| `Brent`, `GoldenSection` | Scalar, derivative-free |
| `ConjugateGradient`, `BFGS` (and L-BFGS), `Adam` | Vector, gradient-based |
| `TrustRegion` (Steihaug-Toint or dogleg) | Vector, Hessian-based |
| `NelderMead`, `PatternSearch` | Vector, derivative-free, local |
| `ParticleSwarm`, `CMAES` | Vector, derivative-free, global |

The descent optimizers share the line searches (Armijo, Goldstein, strong Wolfe, Moré-Thuente) in `Optimist/LineSearch.hh`. Include `Optimist/Solvers.hh` to get every solver at once.

Given a function class, every solver offers two entry points:

- `rootfind(fun, x_ini, x_sol)` solves \f$\mathbf{f}(\mathbf{x}) = \mathbf{0}\f$: root finders solve the system directly, optimizers minimize \f$\frac{1}{2}\|\mathbf{f}\|^2\f$;
- `optimize(fun, x_ini, x_sol)` minimizes: optimizers minimize a cost function \f$f : \mathbb{R}^n \to \mathbb{R}\f$ directly (and vector-valued functions in the least-squares sense), root finders look for the stationary points by solving \f$\nabla\phi = \mathbf{0}\f$.

Plain lambdas can be used through `solve(function, [derivatives...], x_ini, x_sol)`.

## Structure

The library is header-only and lives in `core/include`:

| Header | Content |
|---|---|
| `Optimist.hh` | Type traits, error macros and common definitions |
| `Optimist/Function.hh` | `FunctionBase` and `Function`, the interface of the functions handed to the solvers |
| `Optimist/SolverBase.hh` | Shared solver machinery: problem semantics, counters, limits, bounds, damping, verbose output |
| `Optimist/RootFinder.hh`, `Optimist/Optimizer.hh` | Bases of the two solver families |
| `Optimist/RootFinder/*.hh` | Root finders; `OpenMethod` and `Bracketing` are the bases of the scalar ones, `QuasiNewton` of Broyden and Greenstadt |
| `Optimist/Optimizer/*.hh` | Optimizers; `ScalarMinimizer` is the base of Brent and golden-section search |
| `Optimist/LineSearch.hh` | Line searches shared by the descent optimizers |
| `Optimist/FiniteDifferences.hh` | Finite-difference gradients, Jacobians and Hessians |
| `Optimist/Utils.hh` | Conversions between fixed-size, dynamic-size and sparse vectors |
| `Optimist/Solvers.hh` | Convenience header including every solver |
| `Optimist/TestSet/*.hh` | Test problems with known solutions; `DenseFunction` lets a problem be written once for every vector type |

A new solver derives from `RootFinder` or `Optimizer` (or from one of the family bases above), declares which derivatives it needs (`RequiresFirstDerivative`, `RequiresSecondDerivative`) and implements `solve_impl(...)`; `SolverBase` takes care of the rest.

## Examples

The `examples` folder holds small programs built with `-DOPTIMIST_BUILD_EXAMPLES=ON` (and run by CTest).

`kepler_orbits` follows four satellites around Earth for three days of two-body motion. At every three-minute step it solves Kepler's equation \f$M = E - e\sin E\f$ for the eccentric anomaly with the `NewtonRaphson` and `Halley` root finders, then writes the positions as CSV files and draws the two figures below. More eccentric orbits need more iterations per solve, and Halley's method needs fewer than Newton's.

| Satellite | Eccentricity | Period | Newton iterations | Halley iterations |
|---|---|---|---|---|
| Low Earth orbit | 0.001 | 1.54 h | 2.00 | 2.00 |
| Eccentric orbit | 0.300 | 5.08 h | 3.86 | 2.99 |
| Molniya orbit | 0.740 | 11.99 h | 4.56 | 3.65 |
| Transfer orbit | 0.730 | 10.54 h | 4.53 | 3.63 |

![Satellite orbits around Earth](kepler_orbits.svg)

![Satellite altitude over three days](kepler_altitude.svg)

## Tests

The tests use [`GoogleTest`](https://github.com/google/googletest), which is included as a submodule together with `Eigen`.

```sh
git submodule update --init
cmake -S . -B build -DOPTIMIST_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Configure with `-DOPTIMIST_BENCHMARK=ON` to make `test_benchmark` print a comparison of all optimizers on the cost functions of the test set.

## Author

- Davide Stocco <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: davide.stocco@unitn.it

## License

The `Optimist` project is distributed under the BSD 2-Clause License - see the [LICENSE](https://StoccoDavide.github.io/Optimist/LICENSE) file for details.

Here's what the license entails:

1. Anyone can copy, modify and distribute this software.
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. This software is provided without warranty.
6. The software author or license can not be held liable for any damages inflicted by the software.
