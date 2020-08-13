#include "Potential.h"

Potential::Potential()
{

}

double Potential::HardSphere(float r, float dmol)
{
  double inf = numeric_limits<double>::infinity();
  if(r < dmol){
    return inf;
  }
  else{
    return 0.0;
  }
}
double Potential::LennardJones(float r, float epsilon, float sigma)
{
  float rinv = sigma/r;
  return 4.0*epsilon*(pow(rinv,12)-pow(rinv,6));
}
