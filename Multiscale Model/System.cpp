#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "time.h"
#include "iostream"
using namespace std;

#include "BCTool.h"
#include "CCell.h"
#include "Random.h"
#include "System.h"

extern struct IMD _MD;
extern CRandom Rand;
System::System()
{
  NumCell = 0;
}
System::~System()
{
  delete _cell;
}
System::System(int N0)
{
  _Prolif=1.0;
  NumPoolCell = N0;      // Number of cells in the simulation pool
  NumCell = N0;          // Number of total cell numbers
  MaxNumCell = MaxCellNumber;  // Number of maxium cell number in the simulation pool
  TME.CTL=0;  // intial CTL
  TME.MDSC=0;  // intial MDSC
  _cell = new CCell[MaxNumCell];
}
bool System::Initialized()
{
  int m;
  m=1;
  do
  {
    if((*this)(m).Initialized(m,_MD.cellpar))
    {
      m++;
    }
    else
    {
    }
  } while (m<=MaxNumCell);
  //(*this)(1).Outputparameter();
  return 1;
}
void System::Run()
{
  double t;
  int step;
  int m;
  FILE *fmdsysm;
  char fn[Comlength];
  sprintf(fn,"%s.dat",_MD.mdcrd);
  if((fmdsysm=fopen(fn,"w"))==NULL)
  {
    cout<<"Cannot open the input file"<<endl;
    exit(0);
  }
  step=0;
  OutPutSys(step);
  _N0=NumPoolCell;_N1=0;_N2=0;_N3=0;_N4=0;
  for(t=0;t<=_MD.T1;t=t+_MD.dt)
  {
    step=step+1;
    fprintf(fmdsysm,"%f %f %f %d %d %f %f\n",t,_Prolif, NumCell, NumPoolCell,_N1,TME.CTL,TME.MDSC); //output system result
    if(SystemUpdate(t))
    {
      //Update_TME(NumPoolCell,_MD.dt); // Update Tumormicroenvironment
      if((_MD.ntpx>0)&&(step%_MD.ntpx==0))
      {
        OutPutSys(step);// output cell CTP x and phenotype y
        //fprintf(fmdsysm,"%f %f %f %f %d %d %d\n",t,nCoV,AI.IL6,AI.T1/100000.0,_N0,_N1,_N2); //output system result
      }
    }
  }
  fclose(fmdsysm);
  return ; 
}

bool System::SystemUpdate(double t)
{
    int m;
    int Ntemp;
    float c;
    double *X1, *S1;      //  X1 stores the CTP x, S1 stores y for all cells
    bool *PQ;     //  store the state of _ProlifQ for all cells
    
    X1 = new double[2*NumPoolCell+1];
    S1 = new double[2*NumPoolCell+1];
    PQ = new bool[2*NumPoolCell+1];
    
    Ntemp = 0;        // Number of cells after a cycle.
    // type = 0: unchanged(at the resting phase of during the proliferating phase);
    // type = 1: mitosis, generate two daughter cells;
    // type = 2: cell death, removed from the resting phase
    // type = 3: enter the proliferating state
    _N0=0;           //_N0:Number of cells in the pool after cell fate decision.
    _N1=0;           //_N1:Number of cells removed from resting phase by death;
    _N2=0;           //_N2:Number of cells undergoing mitosis;
    _N3=0;           //_N3:Number of cells removed from the proliferating phase;
    _N4=0;           //_N4:Number of cells undergoing the proliferating phase
    int n=(int)(t/_MD.dt);
   for (m=1; m<=NumPoolCell; m++)
    {
      c=_MD.c[n];
      nextcell=(*this)(m).CellFateDecision(NumCell,_MD.dt,t,c);
      switch(nextcell.type)
      {
        case 2: // cells  death
        _N1++;
        break;
      case 0:   // cells in resting phasse
      case 3: // Enter the proliferating state
        Ntemp++;
        X1[Ntemp]=nextcell.Y1[0];
        S1[Ntemp]=nextcell.Y1[1];
        PQ[Ntemp]=nextcell._ProfQ;
        
        if(nextcell._ProfQ==0)
        {
          _N0++;
        }
        else
        {
          _N4++;
        }
        break;
      case 1:
        Ntemp++;
        X1[Ntemp] = nextcell.Y1[0];
        S1[Ntemp] = nextcell.Y1[1];
        PQ[Ntemp] = nextcell._ProfQ;
                
        Ntemp++;
        X1[Ntemp] = nextcell.Y2[0];
        S1[Ntemp] = nextcell.Y2[1];
        PQ[Ntemp] = nextcell._ProfQ;
        
        _N0 = _N0+2;
        _N2 = _N2+1;
                
        break;
      
      }
    }
    
    _Prolif = 1.0*Ntemp/NumPoolCell;
   NumCell =NumCell * _Prolif;

    if(Ntemp==0)
    {
        return(0);
    }

    int i;
    double p0;

    p0 = 1.0*MaxNumCell/Ntemp;
    m=0;
    for (i=1; i<=Ntemp; i++)
    {
        if(Rand()<p0 && m<MaxNumCell)
        {
            m=m+1;
            (*this)(m).Y[0] = X1[m];
            (*this)(m).Y[1] = S1[m];
            (*this)(m)._ProfQ = PQ[m];
        }
        if (m==MaxNumCell)
        {
            break;
        }
    }
    NumPoolCell = m;
    
    delete X1;
    delete S1;
    delete PQ;
    return(1);
}
/*
TumorMicroenvironment System::Update_TME(float N00,double dt) //ODE23 
{
  double *y,*yp,*yc;
  
  y=new double [2];yp=new double [2];yc=new double [2];
  y[1]=TME.MDSC;y[0]=TME.CTL;

  yp[0]=y[0]+dt*(_MD.s1+_MD.k1*pow(N00,2)*y[0]/(1+_MD.k2*_MD.h1*pow(N00,2))-_MD.k3*y[0]*y[1]-_MD.muT*y[0]);
  yp[1]=y[1]+dt*(_MD.s2*_MD.k4*N00/(1+_MD.h3*_MD.k4*N00)-_MD.muM*y[1]);

  yc[0]=y[0]+dt*(_MD.s1+_MD.k1*pow(N00,2)*yp[0]/(1+_MD.k2*_MD.h1*pow(N00,2))-_MD.k3*yp[0]*yp[1]-_MD.muT*yp[0]);
  yc[1]=y[1]+dt*(_MD.s2*_MD.k4*N00/(1+_MD.h3*_MD.k4*N00)-_MD.muM*yp[1]);

  TME.CTL=0.5*(yc[0]+yp[0]);
  TME.MDSC=0.5*(yc[1]+yp[1]);
  delete y;
  delete yp;
  delete yc;
  return TME;
}
*/
void System::OutPutSys(int step)
{
  FILE *fp;
  int m;
  char fnc[Strlength];
  sprintf(fnc,"%s-%d.dat",_MD.mdcrd,step);

  if((fp=fopen(fnc,"w"))==NULL)
  {
    cout<<"Cannot open the fmdy file"<<endl;
    exit(0);
  }
  for(int m=1;m<=NumPoolCell;m++)
  {
    fprintf(fp,"%d %6.3f %6.3f\n",m, (*this)(m).Y[0],(*this)(m).Y[1]);
  }
  fclose(fp);
}

void System::OutPutCells(double t)
{
	
}

// operator 
CCell&System ::operator()(int index)
{
  if(index>=0&&index<=MaxCellNumber)
  {
    return *(_cell+(index-1));
  }
  else
  {
    cout<<"Err<< System () >>Dimensional error"<<endl;
	  exit(0);
  }
  
}
