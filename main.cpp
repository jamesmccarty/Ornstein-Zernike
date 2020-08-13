#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Closure.h"
#include "Potential.h"
#include "fftpack.h"

using namespace std;

void sinft (float *inputarray, int N, int direc, float deli);

int main()
{

  // input parameters here
  float rho = 0.7;   // number density
  float dmol = 1.0;  // hard sphere diameter
  float amix = 0.1;  // mixing paramter
  //  float kBT = 1.0; // Thermal temperature in reduced units of kBT
  //  float beta = 1.0/kBT // inverse tempearture

  // some other needed parameters here
  int NM=14;
  int N=pow(2,NM);    // number of points
  float dr=0.015;    // r-spacing
  float Tolp=1.E-6;  // tolerance
  float dk = M_PI/(N*dr); // k-spacking
  int niter = 10000; // max number of iterations

  // declare variables and arrays
  float Hr[N], Cr[N], gamma[N], Hk[N], Ck[N], Crn[N];
  float r, Cki, error;
  double Ur;
  int i, iter;

  // initiate potential class and closure class
  Potential potential;
  Closure closure;

  // initialize c(r) to zero
  for(i=0;i<N;i++){
    Cr[i]=0.0;
    gamma[i]=0.0;
    r = (i+1)*dr;
  }

  cout << "Start Picard iteration" << endl;
  for(iter=1;iter<niter;iter++){
    // Calculate c(k) from c(r) using FFT
    for(i=0;i<N;i++){
      Ck[i] = Cr[i];
    }
    sinft(Ck,N,1,dr);
    // Calculate h(k) from Ornstein -Zernike equation
    for(i=0;i<N;i++){
      Cki = Ck[i];
      Hk[i] = Cki/(1.0-rho*Cki);   // Ornstein - Zernike equation
      // A change of variables is necessary in order to use potentials with hard cores
      gamma[i] = (Cki/(1.0-rho*Cki)) - Cki;   // define function gamma(r) = h(r) - c(r)
    }
    // calculate gamma(r) from gamma(k) using FFT
    sinft(gamma,N,-1,dk);
    // calculate new c(r) from the closure relation
    for(i=0;i<N;i++){
      r = (i+1)*dr;
      Ur = potential.HardSphere(r,dmol); // evaluate the potential energy at distance r
    //  Ur = potential.LennardJones(r, 0.25,dmol);
      Crn[i] = closure.PercusYevick(r,Ur, gamma[i]); // apply closure relation to get new c(r)
//      Crn[i] = closure.HyperNettedChain(r, Ur, gamma[i]); // apply closure relation to get new c(r)
//      Crn[i] = closure.MeanSphericalApproximation(r, Ur, gamma[i]); // apply closure relation to get new c(r)
    }
    // calculate errors
    error = 0.0;
    for(i=0;i<N;i++){
      error = error + pow((dr*(Crn[i]-Cr[i])),2);
    }
    error = sqrt(error);
  if( iter % 100 == 0){
    cout << "iteration " <<  iter << "  error " << error << endl;
   }
   if(error < Tolp){
     // stop iteration if error is within the specified tolerance
     break;
   }
   for(i=0;i<N;i++){
     // update c(r) with new c(r). Here we mix the old and new c(r) for numerical stability
     Cr[i] = amix*Crn[i] + (1.0 - amix)*Cr[i];
   }
  }
  if(iter == niter){
    cout << "Picard iteration does not converge" << endl;
  }
  else{
    cout << "End Picard iteration" << endl;
  }

  // Fourier transform h(k) to h(r) using FFT
  for(i=0;i<N;i++)
  {
     Hr[i]=Hk[i];
  }
  sinft(Hr,N,-1,dk);

  // Now write output files
  ofstream hrfile;
  ofstream hkfile;
  ofstream crfile;
  ofstream ckfile;
  hrfile.open("hr.dat");
  hkfile.open("hk.dat");
  crfile.open("cr.dat");
  ckfile.open("ck.dat");
  for(i=0;i<N;i++)
  {
    r=(i+1)*dr;
    hrfile << r << " " << Hr[i] << "\n";
    crfile << r << " " << Cr[i] << "\n";
    hkfile << (i+1)*dk << " " << Hk[i] << "\n";
    ckfile << (i+1)*dk << " " << Ck[i] << "\n";
  }

  hrfile.close();
  crfile.close();
  hkfile.close();
  ckfile.close();
  return 0;
}


void sinft (float *inputarray, int N, int direc, float deli){
  // function to use the netlib fast Fourier transform
  int i;
  float finnorm[N];
  float wsave[3*N+15];
  float delf = M_PI/(N*deli);
  if(direc==1){
    for(i=0;i<N;i++){
      finnorm[i] = (2.0*pow(M_PI,2))/(delf*delf*(i+1)*N);
    }
  }
  else if(direc==-1){
    for(i=0;i<N;i++){
      finnorm[i] = (1.0)/(4.0*M_PI*delf*delf*(i+1)*N);
    }
  }
  else{
    std::cerr << "ERROR: Invalid transform direction paramter\n";
    exit (EXIT_FAILURE);
  }
  for(i=0;i<N;i++){
    inputarray[i]=(deli*(i+1))*inputarray[i];
  }
  sinti(N,wsave);
  sint(N,inputarray,wsave);

  // normalize the transform
  for(i=0;i<N;i++){
    inputarray[i] = inputarray[i]*(finnorm[i]);
  }
  return;
}
