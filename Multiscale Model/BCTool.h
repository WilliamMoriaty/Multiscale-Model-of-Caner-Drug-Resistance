#ifndef BCTOOL_H
#define BCTOOL_H

#define Comlength 100
#define Strlength 1024
#define UnitMAX 2147483574
#define MaxCellNumber 1000000
#define omega 401
struct IMD
{
  char mdcrd[Comlength];
  char cellpar[Comlength];
  char dose[Comlength];
  int N0;                 // Number of cells at the initial state.
	int seed;				// Seed of the random numbers.
  float c[7200];       // drug dose
	double dt;              // Dynamical parameters.
  double T0;              // The time of relaxation from the initial state. We need the run the program for t < T0, then change the system parameter to study the system dynamics
	double T1;              // Dynamical parameters, total time of simulation.
  int ntpr;               // Output the result, if ntpr = 0, then no output.
	int ntpx;                // Output the path result for every ntpx steps. if ntpx = 0, then no result is outputed. 

};
struct Dcell
{
  int type;         // type = 0: unchanged(at the resting phase of during the proliferating phase);
                    // type = 1: mitosis, generate two daughter cells;
                    // type = 2: cell death, removed from the resting phase
                    // type = 3: enter the proliferating state
                    
    double Y1[2];      // State of the daughter cell 1.
    double Y2[2];      // State of the daughter cell 2.
    int _ProfQ;             // States of the cell at the proliferating phase
    
};
struct TumorMicroenvironment
{
  float CTL;
  float MDSC;
};

#endif