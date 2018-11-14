////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 11, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "LoadDataparameters.h"

void LoadDataparameters(string filename,float*zL,float*OffsetLL,float*OffsetRL,float*TrigThresh)
{
//Read the file filename and load the detector geometry information
//For each of the  7 layers of the tracker 
//1st column: layers ID L0 to L6
//2nd column: board ID beteen BA and BI depending on the board used  
//3rd column: Z position of the layer with respect to the bottom of scintillator T3 (in cm)
//4th column: Offset of the left ladder (in cm)
//5th column: Offset of the right ladder (in cm)
    
 //Number of lines in the file
 int n=7; 
 //Prefix 2 characters
 string pre[12]={"L0","L1","L2","L3","L4","L5","L6","T1","T2","T3","T4","Gu"}; 
  
 //Create stream from filename 
 ifstream file;
 file.open(filename, ios_base::in); // open file
 int j=0;
 for(string line; getline(file, line); )   //read stream line by line
   {
    cout << line << endl;
    istringstream in(line);      //make a stream for the line itself
    string prefix;
    in >> prefix;                  //and read the first whitespace-separated token
    string board;//not used
    float ztmp;
    float oLLtmp;
    float oRLtmp;
    float valthres;
   
    //check prefix to load the appropriate region variable 
    for(int i=0;i<n;i++)
      {
       if(prefix.compare(pre[i]) == 0)
        {
         if(i<7)
	   {
            in >> board;//not used
            in >> ztmp;
            in >> oLLtmp;    
            in >> oRLtmp;
	    zL[i]=ztmp;
            OffsetLL[i]=oLLtmp;
            OffsetRL[i]=oRLtmp;
           }
	 else
	   {
            in >> valthres;
    	    TrigThresh[i-7]=valthres;    
	   } 
	 j++;
        }
      }
   }
  if (j!=n) cout << "Error when loading the Data parameters: Wrong number of parameters in the file " << filename << endl;
}
