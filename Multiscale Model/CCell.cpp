#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "time.h"
#include "iostream"

using namespace std;

#include "CCell.h"
#include "System.h"
#include "BCTool.h"
#include "Random.h"

extern CRandom Rand;
extern struct IMD _MD;
CCell::CCell()
{
}
CCell::~CCell()
{
}
bool CCell::Update_y(double dt, float c) // Update phentype s
{
  double eta = Rand.WienerGen() * sqrt(dt);
  double v;
  if (c > 0) // (y<0.8:reflection boundary)
  {
    v = _par.v0 * c * Y[0];
    v = v * (0.8 - Y[1]);
  }
  else
  {
    v = _par.v0 * Y[0] * (0.2 - Y[1]);
  }
  Y[1] = Y[1] + v * dt + sqrt(2 * _par.Sigma) * eta;
  return 1;
}
bool CCell::Initialized(int k, char fn[]) //
{
  double x1;
  double x2;
  ReadPar(fn);
  _cellid = k;
  _ProfQ = 0;
  Y[0] = Rand.Value(0, 1);
  Y[1] = Rand.NormalDistribution(0.2, 0.01);
  return 1;
}
Dcell CCell::CellFateDecision(long N0, float dt, float t, float c) // update
{
  double beta, mu, gamma, k;
  double _rand;
  Dcell nextcell;
  beta = fbeta(N0, dt);  // proliferation rate
  gamma = fgamma(dt, c); // normal death rate
  nextcell.type = 0;
  nextcell._ProfQ = _ProfQ;
  nextcell.Y1[0] = Y[0];
  nextcell.Y1[1] = Y[1];
  if (_ProfQ == 0)
  {
    _rand = Rand();
    // type = 0: unchanged(at the resting phase of during the proliferating phase);
    // type = 1: mitosis, generate two daughter cells;
    // type = 2: cell death, removed from the resting phase
    // type = 3: enter the proliferating state

    if (_rand < gamma) // cell death
    {
      nextcell.type = 2;
    }
    else if (_rand < gamma + beta) // cell enter proliferation phase
    {
      nextcell.type = 3;
      nextcell.Y1[0] = Y[0];
      nextcell.Y1[1] = Y[1];
      nextcell._ProfQ = 1;
    }
    else if (Update_y(dt, c))
    {
      nextcell.Y1[1] = Y[1];
    }
  }
  else if (_ProfQ == 1)
  {
    // mitosis two daughter cells
    nextcell.type = 1;
    nextcell._ProfQ = 0;
    nextcell.Y1[0] = GetnextEpi(Y[0], c, t);
    nextcell.Y2[0] = GetnextEpi(Y[0], c, t);
    nextcell.Y1[1] = Y[1];
    nextcell.Y2[1] = Y[1];
  }
  return (nextcell);
}

double CCell::fbeta(long N, float dt)
{
  double beta;
  double beta1;
  // theta=_MD.k5*pow(M,1)/(1+_MD.k5*_MD.h4*pow(M,1));
  beta1 = _par.beta1 * (_par.m1*Y[0]+pow(_par.m2*Y[0],6))/(1+pow(_par.m3*Y[0],6));
  beta1 = (_par.beta0+beta1) * (1 - N * 1.0 / _par.K);
  // beta1=_par.beta0;
  beta = beta1 * dt;
  return (beta);
}
double CCell::fgamma(float dt, float c)
{
  double gamma;
  if (Y[1] >= _par.yth)
  {
    // gamma=_par.gamma0+_MD.k2*CTL;
    gamma = _par.gamma0;
  }
  else
  {
    // gamma=_par.gamma0+_MD.k2*CTL+_par.gamma1*c;
    gamma = _par.gamma0 + _par.gamma1 * c;
  }
  // gamma=_par.gamma0+_MD.k2*CTL+_par.gamma1*c*(1-Y[1]);
  gamma = gamma * dt; // abnormal cell is cleared
  return (gamma);
}

double CCell::GetnextEpi(double x, float c, float t)
{
  double phi;
  double a, b;
  double z;
  double f;
  double eta, zeta, alpha, psi;
  eta = _par.eta;
  zeta = _par.zeta;
  psi = _par.psi;
  if (c == 0 || t < _MD.T0)
  {
    alpha = _par.alpha1;
  }
  else
  {
    alpha = _par.alpha2;
  }
  phi = zeta + psi * pow(alpha * x, _par.n) / (1 + pow(alpha * x, _par.n));
  a = eta * phi;
  b = eta * (1 - phi);
  z = Rand.BetaDistribution(a, b);

  return (z);
}
void CCell::ReadPar(char fpar[])
{
  FILE *fp;
  char str[Strlength], *pst;
  if ((fp = fopen(fpar, "r")) == NULL)
  {
    cout << "Cannot open the cell parameter file" << endl;
    exit(0);
  }
  rewind(fp);
  while (!feof(fp))
  {
    fgets(str, Strlength, fp);
    if (str[0] == '#')
    {
      continue;
    }

    if ((pst = strstr(str, "m1=")) != NULL)
    {
      _par.m1 = atof(pst + 3);
    }
    if ((pst = strstr(str, "m2=")) != NULL)
    {
      _par.m2 = atof(pst + 3);
    }
    if ((pst = strstr(str, "m3=")) != NULL)
    {
      _par.m3 = atof(pst + 3);
    }
    if ((pst = strstr(str, "gamma0=")) != NULL)
    {
      _par.gamma0 = atof(pst + 7);
    }
    if ((pst = strstr(str, "gamma1=")) != NULL)
    {
      _par.gamma1 = atof(pst + 7);
    }
    if ((pst = strstr(str, "beta0=")) != NULL)
    {
      _par.beta0 = atof(pst + 6);
    }
    if ((pst = strstr(str, "beta1=")) != NULL)
    {
      _par.beta1 = atof(pst + 6);
    }
    if ((pst = strstr(str, "Sigma=")) != NULL)
    {
      _par.Sigma = atof(pst + 6);
    }
    if ((pst = strstr(str, "yth=")) != NULL)
    {
      _par.yth = atof(pst + 4);
    }
    if ((pst = strstr(str, "v0=")) != NULL)
    {
      _par.v0 = atof(pst + 3);
    }
    if ((pst = strstr(str, "K=")) != NULL)
    {
      _par.K = atof(pst + 2);
    }
    if ((pst = strstr(str, "n=")) != NULL)
    {
      _par.n = atof(pst + 2);
    }
    if ((pst = strstr(str, "zeta=")) != NULL)
    {
      _par.zeta = atof(pst + 5);
    }
    if ((pst = strstr(str, "psi=")) != NULL)
    {
      _par.psi = atof(pst + 4);
    }
    if ((pst = strstr(str, "eta=")) != NULL)
    {
      _par.eta = atof(pst + 4);
    }
    if ((pst = strstr(str, "alpha1=")) != NULL)
    {
      _par.alpha1 = atof(pst + 7);
    }
    if ((pst = strstr(str, "alpha2=")) != NULL)
    {
      _par.alpha2 = atof(pst + 7);
    }
  }
  fclose(fp);
}
void CCell::Outputparameter()
{
  printf("m1=%f m2=%f m3=%f alpha1=%f alpha2=%f\n", _par.m1, _par.m2, _par.m3, _par.alpha1, _par.alpha2);
  printf("gamma0=%f gamma1=%f \n", _par.gamma0, _par.gamma1);
  printf("v0=%f yth=%f K=%f\n", _par.v0, _par.yth, _par.K);
  printf("zeta=%f eta=%f \n", _par.zeta, _par.eta);
  printf("beta0=%f Sigma=%f\n", _par.beta0, _par.Sigma);
}
