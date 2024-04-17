#ifndef CELL_H
#define CELL_H

#include "BCTool.h"
#define Omega 100 // size of dicrete epigenetic space
struct Cellpara
{
  double m1,m2,m3;// proliferation
  double gamma0,gamma1;// apoptosis in mitosis or quiescence
  double v0,kappa;// average sensitivity to stress-adaptation
  double Sigma,yth; // threshold of x to v
  double beta0,beta1; //proliferation
  double K; // capacity of TME
  double zeta,eta,alpha1,alpha2,psi,n; // gene partition 

};
class CCell
{
private:
  int _cellid;
  int _ProfQ;// _ProfQ=1 cell in proliferation phase;_ProfQ=0 cell in quiescent phase
  double Y[2];// Y[0]:CTP(x);Y[1]: pheotype (y)
  
  Cellpara _par;
public:
  CCell();
  ~CCell();
  Dcell CellFateDecision(long N0,float dt,float t,float c);
  bool Update_y(double dt,float c);
  bool Initialized(int k,char fn[]);
  void ReadPar(char fpar[]);
  void Outputparameter();
  double fbeta(long N,float dt);  // proliferation rate
  double fgamma(float dt,float c);//  death rate in the quiescenct phase
  
  double GetnextEpi(double x,float c,float t);      // Get the epigenetic state of the daughter cell given the value y after DNA duplication.
  
  friend class System;
};


#endif