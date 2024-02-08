#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <chrono>
using namespace std;

class EMVb
{
public:
  void callMemory();
  void freeMemory();
  void InitializingWaveFunction();
  void set_sign(int sign_in){sign = sign_in;}
  void set_jz(double jz_in){jz = jz_in;}
  double Calculate(double En, double lambda, double eB);
  void PrintWaveFunction(string head);

private:
  double mass = 5;
  double omega_square = 1;
  double h = 0.001;
  int N = 2502;
  int M = 10; // number of angular states
  int N_ipm = 50;  // number of ipm iterations
  int PrintStep;
  int PrintOder;
  string FileNameHead;
  double jz = -0.5;
  int sign = -1; // sign of the j = 1/2 state.
  
  double averageR;
  double averageRR;
  float **Psi;
  float **Phi;

  double V(double r);
  double V_prime(double r);
  gsl_matrix* GetInverse(gsl_matrix *A);
  void gsl_matrix_mul(gsl_matrix *a, gsl_matrix *b, gsl_matrix *c);
};

//-----------------------------------------------
double EMVb::V(double r){
  return 0.5*mass*omega_square*r*r;
}
double EMVb::V_prime(double r){
  return mass*omega_square*r;
}
inline double ajjz(double j, double jz){
  if(j>=fabs(jz)){
    return sqrt(j*j-jz*jz)/j/2.;
  }else{
    return 0;
  }
}
inline double bjjz(double j, double jz){
  return jz*ajjz(j,jz)/(j*j-1.);
}
inline double cjjz(double j, double jz){
  return 0.5*(1.+jz*jz/(j*(j+1.)));
}

gsl_matrix* EMVb::GetInverse(gsl_matrix *A){
  int n = A->size1;   
  int sign = 0;
  gsl_matrix *inverse;
  gsl_matrix *tmpA = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(tmpA, A);
  gsl_permutation *p = gsl_permutation_alloc(n);
  gsl_linalg_LU_decomp(tmpA, p, &sign);
  inverse = gsl_matrix_alloc(n, n);
  gsl_linalg_LU_invert(tmpA, p, inverse);
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
  return inverse;
}
void EMVb::gsl_matrix_mul(gsl_matrix *a, gsl_matrix *b, gsl_matrix *c){
  assert(a->size1 == c->size1);
  assert(a->size2 == b->size1);
  assert(b->size2 == c->size2);
  for (size_t i=0; i<a->size1; i++){
    for (size_t j=0; j<b->size2; j++){
      double sum=0.0;
      for (size_t k=0;k<b->size1;k++){
        sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,sum);
    }
  }
}
void EMVb::callMemory(){
  Psi=(float **)malloc(N*sizeof(float *));
  for(int i=0; i<N; i++){
    Psi[i]=(float *)malloc(M*sizeof(float));
  }
  Phi=(float **)malloc(N*sizeof(float *));
  for(int i=0; i<N; i++){
    Phi[i]=(float *)malloc(M*sizeof(float));
  }
}
void EMVb::freeMemory(){
  for(int i=0; i<N; i++){
    free(Psi[i]);
  }
  free(Psi);
  for(int i=0; i<N; i++){
    free(Phi[i]);
  }
  free(Phi);
}
void EMVb::InitializingWaveFunction(){
  for(int i=0; i<N; i++){
    double r = i*h;
    Psi[i][0] = exp(-r*r/8.) *r;
  }
}
double EMVb::Calculate(double En, double lambda, double qB){
  float V_gyro[M][M];
  float V_dple[M][M];
  float V_tnsr[M][M];
  auto time_0 = chrono::system_clock::now();
  /*
    Wavefunctions are stored at Psi[k][l], where k labels the radial coordinate, 
     and l for the angular momentum state.
    When sign = -1, l = 0, 2, 4, ... stand for $|1/2>_-$, $|5/2>_-$, $|9/2>_-$, ..., 
     whereas l = 1, 3, ... stand for $|3/2>_+$, $|7/2>_+$, ....
  */
  // boundaries must be zero
  for(int j=0; j<M; j++){ // wavefunctions
    Psi[0][j]   = 0.0;
    Psi[N-1][j] = 0.0;
    Phi[0][j]   = 0.0;
    Phi[N-1][j] = 0.0;
  }
  //-------initialize------------
  for(int j=0; j<M; j++) for(int l=0; l<M; l++){ // angular dependent matrices
    V_gyro[j][l] = 0;
    V_dple[j][l] = 0;
    V_tnsr[j][l] = 0;
  }
  //--------------------------
  for(int j=0; j<M; j++){
    int sgn_j = 1-2*(j%2); //(-1)^j
    sgn_j *= sign;
    V_gyro[j][j] = cjjz(j+0.5, jz);
    if (j<M-1){
      V_gyro[j][j+1] = bjjz(j+1.5, jz);
      V_gyro[j+1][j] = V_gyro[j][j+1];
    }
    if (j<M-2){
      V_gyro[j][j+2] = -ajjz(j+2.5, jz)*ajjz(j+1.5,jz);
      V_gyro[j+2][j] = V_gyro[j][j+2];
    }
    V_tnsr[j][j] = 0.5*jz/(j+0.5)/(j+1.5);
    if (j<M-1){
      V_tnsr[j][j+1] = -ajjz(j+1.5, jz);
      V_tnsr[j+1][j] = V_tnsr[j][j+1];
    }
    double ll = j+0.5;
    if (sgn_j==1){ ll += 1;}
    V_dple[j][j] = - sgn_j * jz / ll;
    if ((sgn_j==1) and (j<M-1)) {
      V_dple[j][j+1] = -sqrt(ll*ll - jz*jz) / ll;
      V_dple[j+1][j] = V_dple[j][j+1];
    }
  }

  double two_m = 2*mass;//1;
  float global_C[N][M][M];
  float global_E[N][M][M];
  float global_invDbar[N][M][M];
  for(int k=0; k<N; k++) for(int j=0; j<M; j++) for(int l=0; l<M; l++){
    global_C[k][j][l] = 0; global_E[k][j][l] = 0; global_invDbar[k][j][l] = 0;
  }

  // I. Preparation: assign C, E, and inverse of D-bar at each [k].
  gsl_matrix * local_C = gsl_matrix_alloc (M, M);
  gsl_matrix * local_E = gsl_matrix_alloc (M, M);
  gsl_matrix * local_Dbar = gsl_matrix_alloc (M, M);
  gsl_matrix * local_invDbar = gsl_matrix_alloc (M, M);
  gsl_matrix * local_buffer1 = gsl_matrix_alloc (M, M);
  gsl_matrix * local_buffer2 = gsl_matrix_alloc (M, M);

  for(int j=0; j<M; j++) for(int l=0; l<M; l++){
    gsl_matrix_set(local_invDbar, j, l, 0);
    global_invDbar[0][j][l] = 0;
  }
  for(int k=1; k<N; k++){
    double r = k*h;
    double V_crct = V_prime(r)/(En+mass-V(r));
    for(int j=0; j<M; j++){
      global_C[k][j][j] = -1/two_m/h/h + V_crct/2./h;
      global_E[k][j][j] = -1/two_m/h/h - V_crct/2./h;
    }

    for(int j=0; j<M; j++) for(int l=0; l<M; l++){
      int sgn_j = 1-2*(j%2); //(-1)^j
      sgn_j *= sign;
      double D_res = 0.5*qB*(0.5*qB*V_gyro[j][l] + (r*V_crct-1)*V_dple[j][l] + r*V_crct*V_tnsr[j][l] ) / two_m;
      if (l==j){
        D_res += 2/two_m/h/h - lambda + 
          ( (j+0.5*(1+sgn_j))*(j+1+0.5*(1+sgn_j))/r/r
            -qB*jz - sgn_j*(j+1)*V_crct/r 
            - En*En + mass*mass - V(r)*V(r)
            + 2*En*V(r)) / two_m;
      }
      gsl_matrix_set(local_Dbar, j, l, D_res); // set Dbar to be D
      gsl_matrix_set(local_C, j, l, global_C[k][j][l]);
      gsl_matrix_set(local_E, j, l, global_E[k-1][j][l]);
    }
    gsl_matrix_mul(local_C, local_invDbar, local_buffer1);
    gsl_matrix_mul(local_buffer1, local_E, local_buffer2);
    gsl_matrix_sub(local_Dbar, local_buffer2); // subtract C[k].(bar D[k-1])^{-1}.E[k-1]
    local_invDbar = GetInverse(local_Dbar);    // compute (bar D[k])^{-1} 
    
    for(int j=0; j<M; j++) for(int l=0; l<M; l++){
      double elem = gsl_matrix_get(local_invDbar, j, l);
      global_invDbar[k][j][l] = elem;
    }
  }
  gsl_matrix_free(local_Dbar);

  gsl_matrix * vec_X1   = gsl_matrix_alloc (M, 1);
  gsl_matrix * vec_X2   = gsl_matrix_alloc (M, 1);
  gsl_matrix * vec_Phi  = gsl_matrix_alloc (M, 1);
  gsl_matrix * vec_Psi_ = gsl_matrix_alloc (M, 1);
  float Psi_bar[N][M];
  double lambda_final = lambda;

  //-------------------------------------------------
  // II. inverse power method iteration
  int l_ipm;
  for(l_ipm=1; l_ipm<N_ipm; l_ipm++){
    for(int j=0; j<M; j++){
      gsl_matrix_set(vec_Psi_, j, 0, 0);
    }
    // II.d. forward iteration to compute Psi-bar
    for(int k=1; k<N; k++){
      for(int j=0; j<M; j++) for(int l=0; l<M; l++){
        gsl_matrix_set(local_invDbar, j, l, global_invDbar[k-1][j][l]);
        gsl_matrix_set(local_C, j, l, global_C[k][j][l]);
      }
      gsl_matrix_mul(local_invDbar, vec_Psi_, vec_X1);
      gsl_matrix_mul(local_C, vec_X1, vec_X2);
      for(int j=0; j<M; j++){
        Psi_bar[k][j] = Psi[k][j] - gsl_matrix_get(vec_X2, j, 0);
        gsl_matrix_set(vec_Psi_, j, 0, Psi_bar[k][j]);
      }
    }
    // II.e backward iteration to compute Phi
    for(int j=0; j<M; j++){
      gsl_matrix_set(vec_Phi, j, 0, 0);
    }
    for(int k=N-1; k>0; k--){
      for(int j=0; j<M; j++) for(int l=0; l<M; l++){
        gsl_matrix_set(local_invDbar, j, l, global_invDbar[k][j][l]);
        gsl_matrix_set(local_E, j, l, global_E[k][j][l]);
      }
      gsl_matrix_mul(local_E, vec_Phi, vec_X1); // X1 = E[k].Phi[k+1]
      for(int j=0; j<M; j++){
        gsl_matrix_set(vec_X2, j, 0 , Psi_bar[k][j] - gsl_matrix_get(vec_X1, j, 0));
        // X2 = Psi_var[k] - X1
      }
      gsl_matrix_mul(local_invDbar, vec_X2, vec_Phi); // Phi[k+1] = invDbar[k].X2
      for(int j=0; j<M; j++){
        Phi[k][j] = gsl_matrix_get(vec_Phi, j, 0);
      }
    }
    // II.f renormalize, compute correction
    double sum3 = 0.0;
    double sum4 = 0.0;
    for(int k=1; k<N-1; k++){
      for(int j=0;j<M;j++){
        sum3 += (Phi[k][j]*Phi[k][j]);
        sum4 += (Psi[k][j]*Phi[k][j]);
      }
    }
    for(int k=1; k<N-1; k++){
      for(int j=0; j<M; j++){
        Psi[k][j] = Phi[k][j]/sqrt(sum3*h);
      }
    }
    //cout <<"\t"<< l_ipm <<"\t"<< (sum4*h) << endl;
    if (fabs(lambda + 1./(sum4*h) - lambda_final) < 1e-8){break;} // stop if converges
    lambda_final = lambda + 1./(sum4*h);
  }// end of ipm step loop

  gsl_matrix_free(local_C);
  gsl_matrix_free(local_E);
  gsl_matrix_free(local_buffer1);
  gsl_matrix_free(local_buffer2);
  gsl_matrix_free(local_invDbar);

  gsl_matrix_free(vec_X1);
  gsl_matrix_free(vec_X2);
  gsl_matrix_free(vec_Phi);
  gsl_matrix_free(vec_Psi_);
  //==================================================
  double averageR = 0.0;
  double averageRR= 0.0;
  double normal = 0.0;
  for(int k=1; k<N-1; k++){
    double r = k*h;
    for(int j=0; j<M; j++){
      normal += (Psi[k][j]*Psi[k][j])*h;
      averageR += r*(Psi[k][j]*Psi[k][j])*h;
      averageRR += r*r*(Psi[k][j]*Psi[k][j])*h;
    }
  }
  auto time_1 = chrono::system_clock::now();
  chrono::duration<double> d1 = time_1 - time_0;
  cout << "[" << ((int)(d1.count()*1000))<<"ms/"<< l_ipm <<"step]" << setprecision(10) 
  << qB <<"\t"<< En <<"\t"<< lambda_final <<"\t"<< averageR<< "\t" << sqrt(averageRR) << endl;
  return lambda_final;
}

void EMVb::PrintWaveFunction(string head){
  for(int j=0; j<M; j++){
    string number_lm = to_string(j);
    string wavename = head+"_"+number_lm+".dat";
    ofstream waveOut;
    waveOut.open(wavename);
    for(int k=1; k<N-1; k++){
      double r = k*h;
      waveOut << r << "\t" << Psi[k][j]/r << endl;
    }
    waveOut.close();
  }
}
