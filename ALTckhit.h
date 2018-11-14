
////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

#include <iostream>
#include "headers.h"


#ifndef ALTckhit_H
#define ALTckhit_H

using namespace std;

class ALTckhit:public TObject
{
 private:
 
  float mregMC; // region if MC
  float mtrackMC; // mtrack of track
  int typeMC;   //type of particle if MC
  float eMC;     //total energy if MC (entrance of the tracker layer)
  float xin;      //x coordinate of hit MC (entrance of the tracker layer)
  float yin;      //y coordinate of hit MC (entrance of the tracker layer)
  float zin;      //z coordinate of hit MC (entrance of the tracker layer)
  float xout;     //x coordinate of hit MC (exit the tracker layer)
  float yout;     //y coordinate of hit MC (exit the tracker layer)
  float zout;     //z coordinate of hit MC (exit the tracker layer)
  float age;    //Age of particle if MC  (entrance of the tracker layer)
  float cx;     //cosineX of momentum MC (Average in the tracker layer)
  float cy;     //cosineY of momentum MC (Average in the tracker layer)
  float cz;     //cosineZ of momentum MC (Average in the tracker layer)
  int flag;     //flag value: Track=1, Hit=0
  float DeltaE;     //Energy along the track or at hit location Hit=0  
  
   // Time data from the corressponding "ASI" LINE  
  int year;//Year from first ASI line of the event
  int m;//Month from first ASI line of the event
  int d;//Day from first ASI line of the event
  int h;//Hour from first ASI line of the event
  int mi;//Minute from first ASI line of the event
  int s;//Second from first ASI line of the event
  
  //Raw Data from cluster information
  int L;         //Layer from 0 to 6 top to bottom, for MC & data
  int chip;      //Chip ID: 0 to 11
  int nstrips;   //Number of strips in the cluster
  int nstripsNC; //Number of strips in the next chip if it is a boundary cluster (except boundary chips 5-6)
  int fstrip;    //First strip ID from 0 to 63 
  int fstripID;  //First strip on the layer from 0 767
  int noisy;  //is 1 if one the strip of the cluster is noisy
  int parityerr[2];//Parity error of the clusters that make the hit
  int chiperr[2];//Chip error of the clusters that make the hit
  int overflow[2];//overflow of the clusters that make the hit
  
  //Coordinates of the cluster in cm determined from the raw data
  float x;				//used for MC and data
  float y;				//used for MC and data
  float z;				//used for MC and data

  
  //Information from Pattern Recognition
 
  float xPR;       //x coordinate of "chosen" hit MC (selected by Pattern Recognition)
  float yPR;       //y coordinate of "chosen" hit MC (selected by Pattern Recognition)
  float zPR;       //z coordinate of "chosen" hit MC (selected by Pattern Recognition)
  float cxPR;	   //cosineX from PR fit
  float cyPR;	   //cosineY from PR fit
  float czPR;	   //cosineZ from PR fit
  bool fGhost;	   //flag true = "fake" hit from PR extrapolation, false = true hit	
  bool flagPR;	   // flag  if hit chosen for reconstruction, flag = 0 if not 
 
  //Reconstructed information
  float xreco;      //x coordinate of hit 
  float yreco;      //y coordinate of hit 
  float zreco;      //z coordinate of hit 
  float agereco;    //Age of particle  
  float cxreco;     //cosineX of momentum 
  float cyreco;     //cosineY of momentum 
  float czreco;     //cosineZ of momentum    
  float ereco;     //kinetic energy 
  bool fUsed;      //flag to tell if hit was accepted by reconstruction algorithm
  int k;     	   //kth hit in event  
 
 public: 
   //Constructors
   ALTckhit();// Default
   //Destructor
   ~ALTckhit(){};  
     ////////////////////////////////
   //"setting" member methods
   ////////////////////////////////
   void set_mregMC(float a){mregMC=a;}
   void set_mtrackMC(float a){mtrackMC=a;}
   void set_typeMC(int a){typeMC=a;}
   void set_eMC(float a){eMC=a;}
   void set_xin(float a){xin=a;}
   void set_yin(float a){yin=a;}
   void set_zin(float a){zin=a;}
   void set_xout(float a){xout=a;}
   void set_yout(float a){yout=a;}
   void set_zout(float a){zout=a;}
   void set_age(float a){age=a;}
   void set_cx(float a){cx=a;}
   void set_cy(float a){cy=a;}
   void set_cz(float a){cz=a;}
   void set_flag(int a){flag=a;}
   void set_DeltaE(float a){DeltaE=a;}
   ////////////////////////////////
   void set_xreco(float a){xreco=a;}
   void set_yreco(float a){yreco=a;}
   void set_zreco(float a){zreco=a;}
   void set_cxreco(float a){cxreco=a;}
   void set_cyreco(float a){cyreco=a;}
   void set_czreco(float a){czreco=a;}
   void set_ereco(float a){ereco=a;}
   void set_agereco(float a){agereco=a;}
   void set_fUsed(bool a){fUsed=a;}
   void set_k(int a){k=a;}
   ////////////////////////////////
   void set_year(int a){y=a;}
   void set_m(int a){m=a;}
   void set_d(int a){d=a;}
   void set_h(int a){h=a;}
   void set_mi(int a){mi=a;}
   void set_s(int a){s=a;}       
   void set_L(int a){L=a;}       
   void set_chip(int a){chip=a;}      
   void set_nstrips(int a){nstrips=a;}   
   void set_nstripsNC(int a){nstripsNC=a;}   
   void set_fstrip(int a){fstrip=a;}    
   void set_fstripID(int a){fstripID=a;}  
   void set_noisy(int a){noisy=a;}  
   void set_parityerr(int a, int b){parityerr[0]=a;parityerr[1]=b;}
   void set_chiperr(int a, int b){chiperr[0]=a;chiperr[1]=b;}
   void set_overflow(int a, int b){overflow[0]=a;overflow[1]=b;}
   void set_parityerr(unsigned int a, int b){parityerr[a]=b;}
   void set_chiperr(unsigned int a, int b){chiperr[a]=b;}
   void set_overflow(unsigned int a, int b){overflow[a]=b;}


   void set_x(float a){x=a;}
   void set_y(float a){y=a;}
   void set_z(float a){z=a;}    
   /////////////////////////////////
   
   void set_xPR(float a){xPR=a;}
   void set_yPR(float a){yPR=a;}
   void set_zPR(float a){zPR=a;}
   void set_cxPR(float a){cxPR=a;}
   void set_cyPR(float a){cyPR=a;}
   void set_czPR(float a){czPR=a;}
   void set_fGhost(bool a){fGhost=a;}
   void set_flagPR(bool a){flagPR=a;}
  
   ////////////////////////////////
   //"Getting" member methods
   ////////////////////////////////
   float get_mregMC( ){return mregMC;}
   float get_mtrackMC( ){return mtrackMC;}
   int get_typeMC( ){return typeMC;}
   float get_eMC( ){return eMC;}
   float get_xin( ){return xin;}
   float get_yin( ){return yin;}
   float get_zin( ){return zin;}
   float get_xout( ){return xout;}
   float get_yout( ){return yout;}
   float get_zout( ){return zout;}
   float get_cx( ){return cx;}
   float get_cy( ){return cy;}
   float get_cz( ){return cz;}
   int get_flag( ){return flag;}
   float get_DeltaE( ){return DeltaE;}
     
   ////////////////////////////////
   int get_year(){return year;}
   int get_m(){return m;}
   int get_d(){return d;}
   int get_h(){return h;}
   int get_mi(){return mi;}
   int get_s(){return s;}       
   int get_L(){return L;}         
   int get_chip(){return chip;}      
   int get_nstrips(){return nstrips;}   
   int get_nstripsNC(){return nstripsNC;}   
   int get_fstrip(){return fstrip;}    
   int get_fstripID(){return fstripID;}  
   int get_noisy(){return noisy;} 
   int* get_parityerr(){return parityerr;}//Parity error of the clusters that make the hit
   int* get_chiperr(){return chiperr;}//Chip error of the clusters that make the hit
   int* get_overflow(){return overflow;}//overflow of the clusters that make the hit
   int get_parityerr(int i){if(i<2)return parityerr[i];else return -1;}//Parity error of the clusters that make the hit
   int get_chiperr(int i){if(i<2)return chiperr[i];else return -1;}//Chip error of the clusters that make the hit
   int get_overflow(int i){if(i<2)return overflow[i];else return -1;}//overflow of the clusters that make the hit

   float get_x(){return x;}
   float get_y(){return y;}
   float get_z(){return z;}    
   ////////////////////////////////

   float get_xPR() {return xPR;}
   float get_yPR() {return yPR;}
   float get_zPR() {return zPR;}
   float get_cxPR() {return cxPR;}
   float get_cyPR() {return cyPR;}
   float get_czPR() {return czPR;}
   bool  get_fGhost() {return fGhost;}
   bool  get_flagPR() {return flagPR;}
   
   ////////////////////////////////
   float get_xreco( ){return xreco;}
   float get_yreco( ){return yreco;}
   float get_zreco( ){return zreco;}
   float get_cxreco( ){return cxreco;}
   float get_cyreco( ){return cyreco;}
   float get_czreco( ){return czreco;}
   float get_ereco( ){return ereco;}
   float get_agereco( ){return agereco;}
   bool  get_fUsed() {return fUsed;}
   int get_k( ){return k;}

   ////////////////////////////////
   ClassDef(ALTckhit,1)
};


#endif