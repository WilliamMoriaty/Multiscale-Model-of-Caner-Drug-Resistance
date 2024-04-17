#ifndef SYSTEM_H
#define SYSTEM_H
#include "CCell.h"
#include "BCTool.h"


class System
{
private:
  int NumPoolCell; // Number of cells in the simulation pool.
  double NumCell; // Number of total cells in the animal (in unit 10^6).
  int MaxNumCell;    // Maximal cell number.
  int _N0,_N1,_N2,_N3,_N4;
  //_N0:Number of cells in the pool after cell fate decision.
  //_N1:Number of cells removed from resting phase by cell death;
  //_N2:Number of cells undergoing mitosis;
  //_N3: Number of cells undergoing the proliferating phase
  //_N4:Number of cells remain unchanged at the resting phase
  
  TumorMicroenvironment TME; //Adapated Immune to virus
  CCell *_cell;
  double _Prolif;       // Proliferation rate
public:
  System();
  ~System();
  System(int N0); //Initialized system
  CCell& operator()(int indx); // Edit to return the type of your own class.
  bool Initialized();// Initialize parameter
  void OutPutCells(double t); // output cell result
	void OutPutSys(int step); // output system result
  bool SystemUpdate(double t);
	void Run();
  TumorMicroenvironment Update_TME(float N00,double dt);// Update T cell , MDSC
 
  Dcell nextcell;
  
  friend class CCell;

};




#endif