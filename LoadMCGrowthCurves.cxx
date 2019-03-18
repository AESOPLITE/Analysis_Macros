
////////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Sarah Mechbal, smechbal@ucsc.edu
///    Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, December 21st 2018
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "LoadMCGrowthCurves.h"


void LoadMCGrowthCurves(string filename, float*Depth, float**eMC, float**pMC) {
// Read the file filename and load results from MC atmospheric simulation
//of the Daniel and Stephens growth curve fit functions for 7 energy bins

//Number of lines in the file
int n=24; 
cout << "called LoadMCGrowthCurves " << endl;	
//Create stream from filename 
ifstream file;
file.open(filename, ios_base::in); // open file
int j=0;
for(string line; getline(file, line); )   //read stream line by line
   {
  //  cout << line << endl;
    istringstream in(line);      //make a stream for the line itself
    float depth;
    in >> depth; //and read the first whitespace-separated token
    float ebin0;
    in >> ebin0;
    float pbin0;
    in >> pbin0;	
    float ebin1;
    in >> ebin1;
    float pbin1;
    in >> pbin1;
    float ebin2;
    in >> ebin2;
    float pbin2;
    in >> pbin2;
    float ebin3;
    in >> ebin3;
    float pbin3;
    in >> pbin3;
    float ebin4;
    in >> ebin4;
    float pbin4;
    in >> pbin4;
    float ebin5;
    in >> ebin5;
    float pbin5;
    in >> pbin5;
    float ebin6;
    in >> ebin6;
    float pbin6;
    in >> pbin6;
    float ebin7;
    in >> ebin7;
    float pbin7;
    in >> pbin7;
    float ebin8;
    in >> ebin8;
    float pbin8;
    in >> pbin8;
    float ebin9;
    in >> ebin9;
    float pbin9;
    in >> pbin9;
    float ebin10;
    in >> ebin10;
    float pbin10;
    in >> pbin10;
    float ebin11;
    in >> ebin11;
    float pbin11;
    in >> pbin11;
    float ebin12;
    in >> ebin12;
    float pbin12;
    in >> pbin12;
    float ebin13;
    in >> ebin13;
    float pbin13;
    in >> pbin13;
    float ebin14;
    in >> ebin14;
    float pbin14;
    in >> pbin14;
    float ebin15;
    in >> ebin15;
    float pbin15;
    in >> pbin15;
    if(j>0) {
	   Depth[j-1]=depth;
	   eMC[0][j-1]=ebin0;
	   pMC[0][j-1]=pbin0;
	   eMC[1][j-1]=ebin1;
	   pMC[1][j-1]=pbin1;
	   eMC[2][j-1]=ebin2;
	   pMC[2][j-1]=pbin2;
	   eMC[3][j-1]=ebin3;
	   pMC[3][j-1]=pbin3;
	   eMC[4][j-1]=ebin4;
	   pMC[4][j-1]=pbin4;
	   eMC[5][j-1]=ebin5;
	   pMC[5][j-1]=pbin5;
	   eMC[6][j-1]=ebin6;
	   pMC[6][j-1]=pbin6;
	   eMC[7][j-1]=ebin7;
	   pMC[7][j-1]=pbin7;
	   eMC[8][j-1]=ebin8;
	   pMC[8][j-1]=pbin8;
	   eMC[9][j-1]=ebin9;
   	   pMC[9][j-1]=pbin9;		
	   eMC[10][j-1]=ebin10;
   	   pMC[10][j-1]=pbin10;	
	   eMC[11][j-1]=ebin11;
   	   pMC[11][j-1]=pbin11;	
	   eMC[12][j-1]=ebin12;
   	   pMC[12][j-1]=pbin12;	
	   eMC[13][j-1]=ebin13;
   	   pMC[13][j-1]=pbin13;	
	   eMC[14][j-1]=ebin14;
   	   pMC[14][j-1]=pbin14;	
	   eMC[15][j-1]=ebin15;
   	   pMC[15][j-1]=pbin15;	
	}
   j++;	   
}	
 

} //end function
 