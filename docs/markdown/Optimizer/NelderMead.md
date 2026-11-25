# Nelder-Mead method

The `NelderMead` class implements the downhill simplex method in the form of Lagarias et al. (1998). A simplex of \f$n+1\f$ points is moved by reflection (\f$\rho\f$), expansion (\f$\chi\f$), outside and inside contraction (\f$\gamma\f$) and shrinking (\f$\sigma\f$), using function values only.

- The initial simplex is either axis-aligned around the initial guess (`STANDARD`, choosing the best of \f$\pm\delta\f$ on each axis) or regular (`SPENDLEY`).
- When the simplex degenerates, i.e., its volume normalized by that of a regular simplex with the same diameter falls below a tolerance, it is rebuilt around the best point.
- The adaptive coefficients of Gao and Han (2012), \f$\rho = 1\f$, \f$\chi = 1 + 2/n\f$, \f$\gamma = 3/4 - 1/(2n)\f$, \f$\sigma = 1 - 1/n\f$, can be enabled with `adaptive_mode(true)` to improve the behaviour in higher dimension.

The iteration stops when both the spread of the function values and the simplex diameter are below the tolerance. The simplex is stored densely, so that fixed-size, dynamic-size and sparse vectors are all supported.
