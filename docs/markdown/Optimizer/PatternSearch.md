# Hooke-Jeeves pattern search

The `PatternSearch` class implements the method of Hooke and Jeeves (1961). An exploratory move perturbs each coordinate of the base point by \f$\pm h\f$, keeping every change that lowers the cost (and remembering the successful sign for the next exploration). If the exploration succeeds, a pattern move jumps to \f$2\mathbf{x}_{\text{new}} - \mathbf{x}_{\text{base}}\f$ and explores again around it, accelerating along the successful direction; otherwise the step is reduced, \f$h \leftarrow \rho h\f$. The iteration stops when \f$h\f$ falls below the tolerance.

The method uses function values only and is robust, but it moves along the coordinate directions and therefore progresses slowly in narrow valleys that are not aligned with the axes.
