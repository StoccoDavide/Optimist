# Optimist

`Optimist` is a library for the numerical solution of optimization problems and zero-finding problems. It is based on the [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) library and written in C++17.

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
