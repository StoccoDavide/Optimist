#include "Optimist/ScalarRootFinder/Chandrupatla.hh"
#include "Optimist/TestSet/Quadratic.hh"

int main()
{
  Optimist::ScalarRootFinder::Chandrupatla<double> sol;
  Optimist::TestSet::Quadratic<double> fun;
  return 0;
}
