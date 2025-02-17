#include "Optimist.hh"
#include "Optimist/TestSet.hh"
#include <iostream>

int main()
{
  Optimist::Info(std::cout);
  Optimist::TestSet::TestSetInfo(std::cout);
  return 0;
}
