#ifndef CLOSURE_H
#define CLOSURE_H
#include <iostream>
#include <string>
#include <stdlib.h>
#include <limits>
#include <math.h>

using namespace std;

class Closure
{
  public:
      Closure();
      float PercusYevick(float, double, float);
      float HyperNettedChain(float,double,float);
      float MeanSphericalApproximation(float,double,float);
};

#endif // POTENTIAL_H
