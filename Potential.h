#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <iostream>
#include <string>
#include <stdlib.h>
#include <limits>
#include <math.h>

using namespace std;

class Potential
{
  public:
      Potential();
      double HardSphere(float, float);
      double LennardJones(float,float,float);
};

#endif // POTENTIAL_H
