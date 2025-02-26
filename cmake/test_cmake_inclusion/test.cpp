#include "Optimist.hh"
#include "Optimist/TestSet.hh"
#include <iostream>

int main()
{
  Optimist::ScalarRootFinder::Chandrupatla<double> sol;
  Optimist::TestSet::Quadratic<double> fun;
  return 0;
}
