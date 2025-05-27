#include "Optimist/RootFinder/Newton.hh"
#include "Optimist/TestSet/ScalarFunction/Quadratic.hh"

int main()
{
  Optimist::RootFinder::Newton<double> newton_1;
  Optimist::RootFinder::Newton<double, 2> newton_2;
  Optimist::TestSet::Quadratic<double> fun;
  return 0;
}
