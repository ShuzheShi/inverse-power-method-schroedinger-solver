#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "power_inverse.h"

using namespace std;
int main(int argv, char* argc[])
{
  EMVb ipm;  
  ipm.callMemory();

  int N_iteration = 200;
  ofstream fout;
  double En_list[3] = {6.79506, 8.58193, 10.1978};//these are obtained from imp without b field.
  for (int n=0; n<3; n++){
    double En = En_list[n];//5 + (1.5 + 2*n);
    double En0 = En;
    double En1 = En;
    fout.open("./results/En_"+to_string(n)+".dat");
    for (int ib=0; ib<51; ib++){
      double eB = 0.02*ib;
      double corr = 0;
      En = 2*En0 - En1;
      for(int i=0; i<N_iteration; i++){
        ipm.InitializingWaveFunction();
        En = En + corr / (2+exp(0.02*i)); // factor smaller than unity to ensure stability
        corr = ipm.Calculate(En, 0, eB);
        if (fabs(corr) < 2e-5) {break;} // stop if converges
      }
      En1 = En0; // second last eB
      if(ib==0){
        En1=En;
        ipm.PrintWaveFunction("./results/WFvac_"+to_string(n));
      }
      En0 = En;  // very last eB
      fout << eB <<"\t"<< En << endl;;
    }
    fout.close();
  }

  ipm.freeMemory();
  return 1;
}
