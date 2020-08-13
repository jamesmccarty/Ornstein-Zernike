#include "Closure.h"

Closure::Closure()
{

}


float Closure::PercusYevick(float r, double Ur, float gamma)
{
  return (exp(-Ur)-1.0)*(1.0 + gamma);
}

float Closure::HyperNettedChain(float r, double Ur, float gamma)
{
  return exp(gamma-Ur)-1.0 - gamma;
}

float Closure::MeanSphericalApproximation(float r, double Ur, float gamma)
{
  double inf = numeric_limits<double>::infinity();
  if(Ur==inf){
    // apply hard core condition manually
    return -1.0 - gamma;
  }
  else{
  return -Ur;
  }
}
