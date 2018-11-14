///////////////////////////////////////////////////////////////////////////////////////////////
// Author: Sarah Mechbal, smechbal@ucsc.edu
// SCIPP, Department of Physics, UC Santa Cruz, February 1st 2018
// This file contains all utility functions to be used in the data processing from AESOP-Lite
///////////////////////////////////////////////////////////////////////////////////////////////


#include "headers.h"

//Functions from MakeRawEventMC


float Discretize(int L,vector<float> x, vector<float> y,vector<float> z,vector<float> cz,vector<float>type,vector<float> Edep,int*chip,int* fstrip,int* fstripID,int*nstrip,float offsetLL, float offsetRL,bool MCflag);
float StriptoCoord(int strip,float OffsetLL,float OffsetRL,bool MCflag);
int CoordtoStrip(float Coord,float SecCoord,float OffsetLL,float OffsetRL,bool MCflag);

//functions from ALEvent

