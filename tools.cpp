///////////////////////////////////////////////////////////////////////////////////////////////
// Author: Sarah Mechbal, smechbal@ucsc.edu
// SCIPP, Department of Physics, UC Santa Cruz, February 1st 2018
// This file contains all utility functions to be used in the data processing from AESOP-Lite
///////////////////////////////////////////////////////////////////////////////////////////////


#include "tools.h"
using namespace std;

float Discretize(int L,vector<float> x, vector<float> y,vector<float> z,vector<float> cz,vector<float>type,vector<float> Edep,int*chip,int* fstrip,int* fstripID,int*nstrip,float offsetLL, float offsetRL,bool MCflag)
{
 //cout << "In Discretize" <<endl;
 int Nn=(int )x.size(); 
 // cout << "Number of segments: "  <<Nn << " in Layer "<< L << endl;

 float* X=new float[Nn];
 float* Y=new float[Nn];
 float* Z=new float[Nn];
 float* CZ=new float[Nn];
 float* T=new float[Nn];
 float* E=new float[Nn];
 float* xx=new float[Nn];
 float* Ss=new float[Nn];
 //Output
 float Xout=-999.;
 //Strip pitch in cm
 float  strippitch=0.0228; 

 for(int i=0;i<Nn;i++)
  { 
   if(L==0||L==4||L==6)
    {
     X[i]=x.at(i);
     Y[i]=y.at(i); 
    }
   else if(L==1||L==2||L==3||L==5)
    {
     X[i]=y.at(i);
     Y[i]=x.at(i); //Inverse coordinated for bending plane	
    }
   E[i]=Edep.at(i);  
   Z[i]=z.at(i); 
   CZ[i]=cz.at(i); 
   T[i]=type.at(i); 
   Ss[i]=-1; //Set to defaults
   //cout << "X"<< i << "= " << X[i] << ", " ;
   //cout << "Y"<< i << "= " << Y[i] << ", "  ;
   //cout << "Z"<< i << "= " << Z[i] <<endl;
  }//i  

//ASSUMPTION: the four frames are perfectly aligned without space between each of them

//Get stripID for each segment of the track in the layer
for(int i=0;i<Nn;i++)
 {  
  Ss[i]=CoordtoStrip(X[i],Y[i],offsetLL,offsetRL,MCflag);
  //cout << "StripID="<< Ss[i] <<endl ;
 }//i

//Check if there is any strip touched
bool nostrip=true;
for(int i=0;i<Nn;i++)
 {  
  if(Ss[i]>-1) {nostrip=false;i=Nn;}
 } //i
 
if(nostrip)
 { 
 // cout << "No strip hitted" <<endl;
  return Xout;
 }

//Get first and last a segment in one strip
int kfirst=-1;
int klast=-1;

for(int i=0;i<Nn;i++)
  {
   //cout << "i = " << i << ", y = " << Y[i] << " , Ss[i] = " << Ss[i] << endl;
   if(kfirst==-1 && Ss[i]>-1) {kfirst=i;klast=i;}    
   if(kfirst>-1 && Ss[i]>-1) klast=i;    
  }
  
if(kfirst<0 || klast<0)
 { 
  return Xout; 
 }

 Xout= (StriptoCoord(Ss[kfirst],offsetLL,offsetRL,MCflag)+StriptoCoord(Ss[klast],offsetLL,offsetRL,MCflag))/2.;

 if(Ss[klast]>=Ss[kfirst]) *fstripID=Ss[kfirst];
 else *fstripID=Ss[klast];

 *chip=(int)Ss[kfirst]%12;
 *fstrip=(int)Ss[kfirst]%64;
 
 *nstrip=(int) abs(Ss[klast]-Ss[kfirst])+1;
 
// cout << "Xout="<< Xout <<endl ;

 
 if(*nstrip>100)
  {
   cout << "Number of segments: "  <<Nn << " in Layer "<< L << endl;
   for(int i=0;i<Nn;i++)
     { 
      cout << "X"<< i << "= " << X[i] << ", " ;
      cout << "Y"<< i << "= " << Y[i] << ", "  ;
      cout << "Z"<< i << "= " << Z[i] << ", ";
      cout << "CZ"<< i << "= " << CZ[i] << ", ";
      cout << "Type"<< i << "= " << T[i] << endl;
      //cout << "Edep"<< i << "= " << Edep[i] <<endl;
     }//i
  }
 return Xout;
}

float StriptoCoord(int strip,float OffsetLL,float OffsetRL,bool MCflag)
{
 //Determine the coordinates from the strips
 //Equations from Sarah's email of September 4 2017.
 
 //Center of X,Y from alignment PIN in cm from Robert (January 24th 2018)
 float Xo=9.74745;   
 float Yo=4.48005; 
 
 float coord=-999;
 
 //Strip pitch in cm
 float  strippitch=0.0228;

 //for MC: 
 float N=384;                                //Number of strip in one module
 float Offset1=0.1088;                       //In cm. Offset from the center of the first/last strip to the edge of the module.
 float Offset2=0.0964;                       //In cm. Offset from the start/end of the strips to the edge of the module.

 float offsetLL=OffsetLL;//Left ladder offset
 float offsetRL=OffsetRL;//Right ladder offset   

 if(MCflag)//Assume perfect alignmentof the four pads
  {
   offsetLL=Xo-Offset1-(N-1)*strippitch; 
   offsetRL=Xo+Offset1; 
  }
 
 //First 6 chips: 0 to 5; strip number 0 to 383
 if(strip>=0 &&strip<N)
  {
   coord=offsetLL-Xo+strip*strippitch; 
  }
 //Last 6 chips: 6 to 11; strip number 384 to 767
 if(strip>=N)
  {
   coord=offsetRL-Xo+(strip-N)*strippitch;
  }
 
 return coord;
}



int CoordtoStrip(float Coord,float SecCoord,float OffsetLL,float OffsetRL,bool MCflag)
{
 //Geometry parameters from schematics
 //Used for discretisation
 float N=384;                                //Number of strip in one module
 float Offset1=0.1088;                       //In cm. Offset from the center of the first/last strip to the edge of the module.
 float Offset2=0.0964;                       //In cm. Offset from the start/end of the strips to the edge of the module.
 float sw=0.0056;                            //In cm. Width of a strip
 float sl=8.7572;                            //In cm. Length of a strip in one module = Size of the frame (8.95cm) - 2*Offset2  ??????????????????????????????
 float sh=0.04;                              //In cm. Height of a strip
 float percsh=0.2;                           //fraction of sh crossed by track to trigger a strip: To be fine tuned later
 float DeltaMax=0.0114;                      //In cm. Maximum allowed distance to the center of a strip
 float strippitch=0.0228;                    //In cm. Strip pitch 
 int strip=-1;
 
 //Center of X,Y from alignment PIN in cm from Robert (January 24th 2018)
 float Xo=9.74745;   
 float Yo=4.48005; 
 float offsetLL=OffsetLL;//Left ladder offset
 float offsetRL=OffsetRL;//Right ladder offset   
 if(MCflag)//Assume perfect alignmentof the four pads
  {
   offsetLL=Xo-Offset1-(N-1)*strippitch; 
   offsetRL=Xo+Offset1; 
  }
 float strip0=offsetLL;       //In cm. Position of the center of strip 0 

 //First: Remove hits outside the frame along Y
 if(SecCoord>sl+Offset2+DeltaMax) return strip;
 if(SecCoord<-sl-Offset2-DeltaMax) return strip;

 //Second: Remove hits outside the frame along X
 if(Coord>offsetRL-Xo+(N-1)*strippitch+DeltaMax) return strip;
 if(Coord<offsetLL-Xo+strippitch-DeltaMax) return strip;

 //Third: Remove Inner Cross
 if(Coord<offsetRL-Xo-DeltaMax && Coord>offsetLL-Xo+(N-1)*strippitch+DeltaMax) return strip;
 if(abs(SecCoord)<Offset2-DeltaMax) return strip;
 
 //Fourth: Determine strip ID from 0 to 767 
 float relX=0;
 int IrelX=0; 
 float RrelX=0;

 if(Coord<0)
  {
   relX=Coord+Xo-offsetLL;   
   IrelX=(int)(relX/strippitch);
   RrelX=relX-(float)IrelX*strippitch;  
   if(RrelX>DeltaMax)IrelX++;
   else if(RrelX<-DeltaMax)IrelX--;
   if(IrelX>=0 && IrelX<N) strip=IrelX; 
  }  
  
 if(Coord>0) 
  {
   relX=Coord+Xo-offsetRL;   
   IrelX=(int)(relX/strippitch);
   RrelX=relX-(float)IrelX*strippitch;  
   if(RrelX>DeltaMax)IrelX++;
   else if(RrelX<-DeltaMax)IrelX--;
   if(IrelX>=0 && IrelX<N) strip=IrelX+N; 
  }    
 return strip;
    
}





