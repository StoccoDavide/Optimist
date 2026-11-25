/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_SOLVERS_HH
#define OPTIMIST_SOLVERS_HH

/**
 * \file Solvers.hh
 * \brief Convenience header including every solver of the library.
 */

// Line searches
#include "Optimist/LineSearch.hh"

// Root finders
#include "Optimist/RootFinder/Algo748.hh"
#include "Optimist/RootFinder/Broyden.hh"
#include "Optimist/RootFinder/Chandrupatla.hh"
#include "Optimist/RootFinder/Chebyshev.hh"
#include "Optimist/RootFinder/Greenstadt.hh"
#include "Optimist/RootFinder/Halley.hh"
#include "Optimist/RootFinder/ITP.hh"
#include "Optimist/RootFinder/LevenbergMarquardt.hh"
#include "Optimist/RootFinder/Newton.hh"
#include "Optimist/RootFinder/NewtonRaphson.hh"
#include "Optimist/RootFinder/Varona.hh"

// Optimizers
#include "Optimist/Optimizer/Adam.hh"
#include "Optimist/Optimizer/BFGS.hh"
#include "Optimist/Optimizer/Brent.hh"
#include "Optimist/Optimizer/CMAES.hh"
#include "Optimist/Optimizer/ConjugateGradient.hh"
#include "Optimist/Optimizer/GoldenSection.hh"
#include "Optimist/Optimizer/NelderMead.hh"
#include "Optimist/Optimizer/ParticleSwarm.hh"
#include "Optimist/Optimizer/PatternSearch.hh"
#include "Optimist/Optimizer/TrustRegion.hh"

#endif  // OPTIMIST_SOLVERS_HH
