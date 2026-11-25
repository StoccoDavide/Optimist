# Particle swarm optimization

The `ParticleSwarm` class implements a global-best particle swarm. Each particle \f$i\f$ has a position \f$\mathbf{x}_i\f$, a velocity \f$\mathbf{v}_i\f$ and a personal best \f$\mathbf{p}_i\f$, and the swarm keeps the global best \f$\mathbf{g}\f$:
\f[
  \mathbf{v}_i \leftarrow w\mathbf{v}_i + \phi_p r_p (\mathbf{p}_i - \mathbf{x}_i) + \phi_g r_g (\mathbf{g} - \mathbf{x}_i) \text{,}
  \qquad \mathbf{x}_i \leftarrow \mathbf{x}_i + \mathbf{v}_i \text{,}
\f]
with \f$r_p, r_g \sim U(0, 1)\f$ drawn per component. The inertia weight \f$w\f$ follows a strategy functor (constant, linear decrease from 0.9 to 0.4 by default, or one of three exponential decreases), and the velocities are clamped to a fraction of the search box.

The search box is given by the solver bounds or, if they are not finite, by a cube of given radius around the initial guess, which is also the first particle. The swarm has converged when all personal bests have collapsed onto the global best. The random generator is seeded deterministically, see `seed(...)`.
