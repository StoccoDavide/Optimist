#include "Optimist/RootFinder/Newton.hh"
#include "Optimist/RootFinder/NewtonRaphson.hh"
#include "Optimist/TestSet/Quadratic.hh"

int main()
{
  Optimist::RootFinder::NewtonRaphson<double> newton_raphson;
  Optimist::RootFinder::Newton<double, 1> newton_1;
  Optimist::RootFinder::Newton<double, 2> newton_2;
  Optimist::TestSet::Quadratic<double> fun;
  return 0;
}
