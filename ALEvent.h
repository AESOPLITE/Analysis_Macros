////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 
#ifndef __ALEVENT__
#define __ALEVENT__
#include "headers.h"
#include "ALTckhit.h"
#include "tools.h"


class ALEvent:public TObject
{
 private:
  
   int eventnumber; //Event number
   
   //Data FROM "PHA" LINE, trigger information
   int yPHA;//Year from PHA line linked to the event
   int mPHA;//Month from PHA line linked to the event
   int dPHA;//Day from PHA line linked to the event
   int hPHA;//Hour from PHA line linked to the event
   int miPHA;//Minute from PHA line linked to the event
   int sPHA;//Second from PHA line linked to the event
   double GoPHA;//Go counter
   double tPHA;//timer
   
   //Data FROM "EVT" LINE, tracker information
   int yEVT;//Year from EVT line linked to the event
   int mEVT;//Month from EVT line linked to the event
   int dEVT;//Day from EVT line linked to the event
   int hEVT;//Hour from EVT line linked to the event
   int miEVT;//Minute from EVT line linked to the event
   int sEVT;//Second from EVT line linked to the event  
   string EVT; //Data from EVT line linked to the event  
   double GoEVT;//Go counter
   double tEVT;//timer
   
   int nHitLEVT;//nHitL from EVT line linked to the event added 12/05/2017
   int CCEVT;//CC from EVT line linked to the event   added 12/05/2017
   int PatternEVT;//CC from EVT line linked to the event   added 12/05/2017
   int Q1EVT;//Q1 from EVT line linked to the event   added 12/05/2017
   double TrigEVT;//Q1 from EVT line linked to the event   added 12/05/2017
   
   //Data FROM "ASI" LINE
   string L[7];//Data from  ASI lines of the event
   int flagL[7];//1 if ASI line was present
   
   
   //Monte Carlo information: Truth variable names finish with MC 
   int ncase; 
   int typeMC; 					//type of particle
   double EkMC;   				//kinetic energy of the particle at injection point
   double pMC;	        			//momentum of particle at injection point
   double X0MC,Y0MC,Z0MC;			//Coordinates of the partcle at the injection point 
   double CX0MC,CY0MC,CZ0MC; 			//Incidence cosines of the partcle at the injection point 
   int Nhits; 					//Number of hits in the event (MC && data)
   
   
   //Pattern Recognition info
   
   double EkPR, p0PR;				//kinetic energy and momentum of particle from least squares fit
   double chi2NBPR, chi2BPR, clNBPR, clBPR;	//chi2 of parabolic/linear fit in bending/nonbending plane
   double aPR, bPR, cPR;			//parameters of parabolic fit ( y(x) = c + bx + ax*x)
   double slopePR, interPR;			//parameters of linear fit 
   double deflecPR;                       	//deflection from layer 2 to layer 6 in the beding plane: Difference of the slope of straight line
   
  //Reconstruction information: variables finish with reco
   int typereco;                    		//type of particle
   double Ekreco, p0reco;           		//kinetic energy and momentum of the particle
   double X0reco,Y0reco,Z0reco;     		//Coordinates of the partocle at the injection point 
   double CX0reco,CY0reco,CZ0reco;  		//Incidence cosines of the particle at the injection point 
   int ndf;                         		// number of degrees of freedom
   double chi2, cl;                 		//chi2 of fit and confidence level (cl = Prob(chi2, ndf)
   double d0, phi0, cpa, dz, tanl;  		//reconstructed helical track parameters
   double phi0_init, cpa_init, tanl_init;  	//initial helical track parameters
   double d0err2, phi0err2, cpaerr2, dzerr2, tanlerr2; //err^2 of track parameters
   TMatrixD Cov_init, Cov_last;			//covariance matrix at initialization and last site
   
   //Hits information
   std::vector<ALTckhit*> hits;  
   
   //NEW: all segments informations for MC (full position & momentum & age & type)
   std::vector<float> posX;
   std::vector<float> posY;
   std::vector<float> posZ;
   std::vector<float> posType;
   std::vector<float> posAge;
   std::vector<float> posP;
   //External Triggers
   bool T1;
   bool T2;
   bool T3;
   bool T4;
   bool guard;
   //Energy for MC, PHA pulse height for data
   std::vector<double> EneT1;  
   std::vector<double> EneT2;  
   std::vector<double> EneT3;
   std::vector<double> EneT4;
   std::vector<double> Eneg; 
   std::vector<double> PHA6; //for data only 6th PHA 
   //For MC, energy deposit in insulation and shell
   std::vector<double> EneIsofoam;  
   std::vector<double> EneShell;  
     
   //Timing for MC, Not available for data
   std::vector<double> timeT1;  
   std::vector<double> timeT2;  
   std::vector<double> timeT3;
   std::vector<double> timeT4;
   std::vector<double> timeg;     
   std::vector<double> timeIsofoam;  
   std::vector<double> timeShell;  
   //Internal Triggers 
   //Tracker layer for 0 to 6. top layer is 0
   //integer value coded with 7 bits values
   //\Sum\Limits_{k=0}^{6} 2^k
   int Ti=0;
    //HOUSEKEEPING FROM COUNTERS 1 AND 3
   //Data FROM "CT1" LINE
   int yCT1;//Year from CT1 line linked to the event (last read CT1 line)
   int mCT1;//Month from CT1 line linked to the event (last read CT1 line)
   int dCT1;//Day from CT1 line linked to the event (last read CT1 line)
   int hCT1;//Hour from CT1 line linked to the event (last read CT1 line)
   int miCT1;//Minute from CT1 line linked to the event (last read CT1 line)
   int sCT1;//Second from CT1 line linked to the event (last read CT1 line)
   float TempCT1; //Temperature measured on the board CT1

   int OnTimeCT1;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
   int LastCT1;//The last command received by the payload, expressed as a decimal number (it is in HeX on the GUI display)
   int CountCT1;//Count of commands received by the payload since power on
   
   //Barometer information: NOT INTERPRETED from line CT1
   float Baro1T;//Barometer 1 Temperature
   float Baro1P;//Barometer 1 Pressure
   float Baro2T;//Barometer 2 Temperature
   float Baro2P;//Barometer 2 Pressure
   //Barometer information: INTERPRETED from line CT1
   float TempB1;//Barometer 1 Temperature
   float TempB2;//Barometer 2 Temperature
   float PressB1;//Barometer 1 Pressure
   float PressB2;//Barometer 2 Pressure
   
   float GOCT1;
   float coinCT1;
   
   //Voltages
   float Volt5VCT1;  // Positive 5V from line CT1
   float Volt15VCT1; // Positive 15V from line CT1
 
   //Data FROM "CT3" LINE
   int yCT3;//Year from CT3 line linked to the event (last read CT3 line)
   int mCT3;//Month from CT3 line linked to the event (last read CT3 line)
   int dCT3;//Day from CT3 line linked to the event (last read CT3 line)
   int hCT3;//Hour from CT3 line linked to the event (last read CT3 line)
   int miCT3;//Minute from CT3 line linked to the event (last read CT3 line)
   int sCT3;//Second from CT3 line linked to the event (last read CT3 line)
   float TempCT3; //Temperature measured on the board CT3

   int OnTimeCT3;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
   int LastCT3;//The last command received by the payload, expressed as a decimal number (it is in HeX on the GUI display)
   int CountCT3;//Count of commands received by the payload since power on
  
   //Voltages
   float Volt5VCT3=-999;  // Positive 5V from line CT3
   float Volt15VCT3=-999; // Positive 15V from line CT3
  
   //TRIGGER RATES (PHA AND LOGIC) from CT3
   float T1L;
   float T1A;
   float T2L;
   float T2A;
   float T3L;
   float T3A;
   float T4L;
   float T4A;
   float GRDL;
   float GRDA;
  
   //HOUSEKEEPING FROM POW
   int yPOW;//Year from POW line linked to the event (last read POW line)
   int mPOW;//Month from POW line linked to the event (last read POW line)
   int dPOW;//Day from POW line linked to the event (last read POW line)
   int hPOW;//Hour from POW line linked to the event (last read POW line)
   int miPOW;//Minute from POW line linked to the event (last read POW line)
   int sPOW;//Second from POW line linked to the event (last read POW line)
   int OnTimePOW;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
   float MainC;
   float MainV;
   float HeatC;
   float HeatV;
   float TrackC;
   float TrackV;  
 public: 
   //Constructors
   ALEvent();// Default
   //Destructor
   ~ALEvent(){};   
   ///////////////////////////////
   // Methods
   ///////////////////////////////
   void Copy(ALEvent*);
   ////////////////////////////////
   //"setting" member methods
   ////////////////////////////////
   void set_eventnumber(int a){eventnumber=a;}
   
   void set_yPHA(int a){yPHA=a;}
   void set_mPHA(int a){mPHA=a;}
   void set_dPHA(int a){dPHA=a;}
   void set_hPHA(int a){hPHA=a;}
   void set_miPHA(int a){miPHA=a;}
   void set_sPHA(int a){sPHA=a;}
   void set_GoPHA(double a){GoPHA=a;}
   void set_tPHA(double a){tPHA=a;}
   void set_yEVT(int a){yEVT=a;}
   void set_mEVT(int a){mEVT=a;}
   void set_dEVT(int a){dEVT=a;}
   void set_hEVT(int a){hEVT=a;}
   void set_miEVT(int a){miEVT=a;}
   void set_sEVT(int a){sEVT=a;}
   void set_EVT(string a){EVT=a;}
   void set_GoEVT(double a){GoEVT=a;}
   void set_tEVT(double a){tEVT=a;}

   void set_nHitLEVT(int a){nHitLEVT=a;}
   void set_CCEVT(int a){CCEVT=a;}
   void set_PatternEVT(int a){PatternEVT=a;}
   void set_Q1EVT(int a){Q1EVT=a;}
   void set_TrigEVT(double a){TrigEVT=a;}
 
   void set_L(int k, string a){if(k<7)L[k]=string(a);}
   void set_flagL(int k, int a){if(k<7)flagL[k]=a;}
   
   ////////////////////////////////  
   void set_ncase(int a){ncase=a;}
   void set_typeMC(int a){typeMC=a;}
   void set_EkMC(double a){EkMC=a;}
   void set_pMC(double a){pMC=a;}
   void set_X0MC(double a){X0MC=a;}
   void set_Y0MC(double a){Y0MC=a;}
   void set_Z0MC(double a){Z0MC=a;}
   void set_CX0MC(double a){CX0MC=a;}
   void set_CY0MC(double a){CY0MC=a;}
   void set_CZ0MC(double a){CZ0MC=a;}
   ////////////////////////////////
   void set_Nhits(int a){Nhits=a;}
   void add_Nhits(){Nhits++;}
   ////////////////////////////////  
   
   void set_EkPR(double a){EkPR=a;}
   void set_p0PR(double a){p0PR=a;}
   void set_aPR(double b){aPR=b;}
   void set_bPR(double a){bPR=a;}
   void set_cPR(double a){cPR=a;}
   void set_interPR(double a){interPR=a;}
   void set_slopePR(double a){slopePR=a;}
   void set_chi2BPR(double a){chi2BPR=a;}
   void set_chi2NBPR(double a){chi2NBPR=a;}
   void set_clBPR(double a){clBPR=a;}
   void set_clNBPR(double a){clNBPR=a;}
   void set_deflecPR(double a){deflecPR=a;}
   
   ////////////////////////////////
   void set_typereco(int a){typereco=a;}
   void set_Ekreco(double a){Ekreco=a;}
   void set_p0reco(double a){p0reco=a;}
   void set_X0reco(double a){X0reco=a;}
   void set_Y0reco(double a){Y0reco=a;}
   void set_Z0reco(double a){Z0reco=a;}
   void set_CX0reco(double a){CX0reco=a;}
   void set_CY0reco(double a){CY0reco=a;}
   void set_CZ0reco(double a){CZ0reco=a;}
   void set_ndf(int a){ndf=a;}
   void set_chi2(double a){chi2=a;}
   void set_cl(double a){cl=a;}
   void set_d0(double a){d0=a;}
   void set_phi0(double a){phi0=a;}
   void set_cpa(double a){cpa=a;}
   void set_dz(double a){dz=a;}
   void set_tanl(double a){tanl=a;}
   void set_phi0_init(double a){phi0_init=a;}
   void set_cpa_init(double a){cpa_init=a;}
   void set_tanl_init(double a){tanl_init=a;}
   void set_d0err2(double a){d0err2=a;}
   void set_phi0err2(double a){phi0err2=a;}
   void set_cpaerr2(double a){cpaerr2=a;}
   void set_dzerr2(double a){dzerr2=a;}
   void set_tanlerr2(double a){tanlerr2=a;}
   void set_Cov_init(TMatrixD a){Cov_init=a;}
   void set_Cov_last(TMatrixD a){Cov_last=a;}

   
   void add_hit(ALTckhit* h){hits.push_back(h);Nhits++;}
   
   //Set Pattern Reco variable at the position of the hit of index k
   void set_hxPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_xPR(a);}
   void set_hyPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_yPR(a);}
   void set_hzPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_zPR(a);}
   void set_hcxPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_cxPR(a);}
   void set_hcyPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_cyPR(a);}
   void set_hczPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_czPR(a);}
   
   //Set reconstruted variable at the position of the hit of index k
   void set_hxreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_xreco(a);}
   void set_hyreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_yreco(a);}
   void set_hzreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_zreco(a);}
   void set_hagereco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_agereco(a);}
   void set_hcxreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_cxreco(a);}
   void set_hcyreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_cyreco(a);}
   void set_hczreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_czreco(a);}
   void set_hereco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_ereco(a);}
   
   //set all positions in vector
   void add_posX(float a){posX.push_back(a);} 
   void add_posY(float a){posY.push_back(a);} 
   void add_posZ(float a){posZ.push_back(a);} 
   void add_posType(float a){posType.push_back(a);} 
   void add_posAge(float a){posAge.push_back(a);}
   void add_posP(float a){posP.push_back(a);} 
   void set_T1(bool a){T1=a;}
   void set_T2(bool a){T2=a;}
   void set_T3(bool a){T3=a;}
   void set_T4(bool a){T4=a;}
   void set_guard(bool a){guard=a;}
   void add_EneT1(double a){EneT1.push_back(a);} 
   void add_EneT2(double a){EneT2.push_back(a);}
   void add_EneT3(double a){EneT3.push_back(a);}
   void add_EneT4(double a){EneT4.push_back(a);}
   void add_Eneg(double a){Eneg.push_back(a);}
   void add_PHA6(double a){PHA6.push_back(a);}
   void add_timeT1(double a){timeT1.push_back(a);}
   void add_timeT2(double a){timeT2.push_back(a);}
   void add_timeT3(double a){timeT3.push_back(a);}
   void add_timeT4(double a){timeT4.push_back(a);}
   void add_timeg(double a){timeg.push_back(a);}   
   void set_Ti(int a){Ti=a;}
   void add_EneIsofoam(double a){EneIsofoam.push_back(a);} 
   void add_EneShell(double a){EneShell.push_back(a);}   
   void add_timeIsofoam(double a){timeIsofoam.push_back(a);}
   void add_timeShell(double a){timeShell.push_back(a);}

  //HOUSEKEEPING FROM COUNTERS 1 AND 3   
   void set_yCT1(int a){yCT1=a;}
   void set_mCT1(int a){mCT1=a;}
   void set_dCT1(int a){dCT1=a;}
   void set_hCT1(int a){hCT1=a;}
   void set_miCT1(int a){miCT1=a;}
   void set_sCT1(int a){sCT1=a;}
   void set_TempCT1(float a){TempCT1=a;}
   void set_OnTimeCT1(int a){OnTimeCT1=a;}
   void set_LastCT1(int a){LastCT1=a;}
   void set_CountCT1(int a){CountCT1=a;}
   void set_Baro1T(float a){Baro1T=a;}
   void set_Baro1P(float a){Baro1P=a;}
   void set_Baro2T(float a){Baro2T=a;}
   void set_Baro2P(float a){Baro2P=a;}
   void set_TempB1(float a){TempB1=a;} 
   void set_TempB2(float a){TempB2=a;} 
   void set_PressB1(float a){PressB1=a;}
   void set_PressB2(float a){PressB2=a;}
   void set_Volt5VCT1(float a){Volt5VCT1=a;}  
   void set_Volt15VCT1(float a){Volt15VCT1=a;} 
   void set_GOCT1(float a){GOCT1=a;} 
   void set_coinCT1(float a){coinCT1=a;} 

   //Data FROM "CT3" LINE
   void set_yCT3(int a){yCT3=a;}
   void set_mCT3(int a){mCT3=a;}
   void set_dCT3(int a){dCT3=a;}
   void set_hCT3(int a){hCT3=a;}
   void set_miCT3(int a){miCT3=a;}
   void set_sCT3(int a){sCT3=a;}
   void set_TempCT3(float a){TempCT3=a;}
   void set_OnTimeCT3(int a){OnTimeCT3=a;}
   void set_LastCT3(int a){LastCT3=a;}
   void set_CountCT3(int a){CountCT3=a;}
   void set_T1L(float a){T1L=a;}
   void set_T1A(float a){T1A=a;}
   void set_T2L(float a){T2L=a;}
   void set_T2A(float a){T2A=a;}
   void set_T3A(float a){T3A=a;}
   void set_T3L(float a){T3L=a;}
   void set_T4A(float a){T4A=a;}
   void set_T4L(float a){T4L=a;}
   void set_GRDL(float a){GRDL=a;}
   void set_GRDA(float a){GRDA=a;}
   void set_Volt5VCT3(float a){Volt5VCT3=a;}  
   void set_Volt15VCT3(float a){Volt15VCT3=a;}
   
   
   void set_yPOW(int a){yPOW=a;}
   void set_mPOW(int a){mPOW=a;}
   void set_dPOW(int a){dPOW=a;}
   void set_hPOW(int a){hPOW=a;}
   void set_miPOW(int a){miPOW=a;}
   void set_sPOW(int a){sPOW=a;}
   void set_OnTimePOW(int a){OnTimePOW=a;}

   void set_MainC(float a){MainC=a;}
   void set_MainV(float a){MainV=a;}
   void set_HeatC(float a){HeatC=a;}
   void set_HeatV(float a){HeatV=a;}
   void set_TrackC(float a){TrackC=a;}
   void set_TrackV(float a){TrackV=a;}
   ////////////////////////////////
   //"Getting" member methods
   ////////////////////////////////
   int get_eventnumber(){return eventnumber;}

   int get_yPHA(){return yPHA;}
   int get_mPHA(){return mPHA;}
   int get_dPHA(){return dPHA;}
   int get_hPHA(){return hPHA;}
   int get_miPHA(){return miPHA;}
   int get_sPHA(){return sPHA;}
   double get_GoPHA(){return GoPHA;}
   double get_tPHA(){return tPHA;}
   int get_yEVT(){return yEVT;}
   int get_mEVT(){return mEVT;}
   int get_dEVT(){return dEVT;}
   int get_hEVT(){return hEVT;}
   int get_miEVT(){return miEVT;}
   int get_sEVT(){return sEVT;}
   double get_GoEVT(){return GoEVT;}
   double get_tEVT(){return tEVT;}
   string get_EVT(){return EVT;}
   
   int get_nHitLEVT(){return nHitLEVT;}
   int get_CCEVT(){return CCEVT;}
   int get_PatternEVT(){return PatternEVT;}
   int get_Q1EVT(){return Q1EVT;}
   double get_TrigEVT(){return TrigEVT;}
  
   string get_L(int k){if(k<7)return L[k];else return "";}
   int get_flagL(int k){if(k<7)return flagL[k];else return 0;}
   
   ////////////////////////////////
  
   int get_ncase(){return ncase;}
   int get_typeMC(){return typeMC;}
   double get_EkMC(){return EkMC;}
   double get_pMC(){return pMC;}
   double get_X0MC(){return X0MC;}
   double get_Y0MC(){return Y0MC;}
   double get_Z0MC(){return Z0MC;}
   double get_CX0MC(){return CX0MC;}
   double get_CY0MC(){return CY0MC;}
   double get_CZ0MC(){return CZ0MC;}


   ////////////////////////////////
   int get_Nhits(){return Nhits;}
   ////////////////////////////////
   
   double get_EkPR(){return EkPR;}
   double get_p0PR(){return p0PR;}
   double get_aPR(){return aPR;}
   double get_bPR(){return bPR;}
   double get_cPR(){return cPR;}
   double get_interPR(){return interPR;}
   double get_slopePR(){return slopePR;}
   double get_chi2BPR(){return chi2BPR;}
   double get_chi2NBPR(){return chi2NBPR;}
   double get_clBPR(){return clBPR;}
   double get_clNBPR(){return clNBPR;}
   double get_deflecPR(){return deflecPR;}

   ////////////////////////////////
   int get_typereco(){return typereco;}
   double get_Ekreco(){return Ekreco;}
   double get_p0reco(){return p0reco;}
   double get_X0reco(){return X0reco;}
   double get_Y0reco(){return Y0reco;}
   double get_Z0reco(){return Z0reco;}
   double get_CX0reco(){return CX0reco;}
   double get_CY0reco(){return CY0reco;}
   double get_CZ0reco(){return CZ0reco;}
   int get_ndf(){return ndf;}
   double get_chi2(){return chi2;}
   double get_cl(){return cl;}
   double get_d0(){return d0;}
   double get_phi0(){return phi0;}
   double get_cpa(){return cpa;}
   double get_dz(){return dz;}
   double get_tanl(){return tanl;}
   double get_phi0_init(){return phi0_init;}
   double get_cpa_init(){return cpa_init;}
   double get_tanl_init(){return tanl_init;}
   double get_d0err2(){return d0err2;}
   double get_phi0err2(){return phi0err2;}
   double get_cpaerr2(){return cpaerr2;}
   double get_dzerr2(){return dzerr2;}
   double get_tanlerr2(){return tanlerr2;}
   TMatrixD get_Cov_init(){return Cov_init;}
   TMatrixD get_Cov_last(){return Cov_last;}
   std::vector<ALTckhit*>& get_hits(){return hits;}  
   std::vector<float>&  get_posX(){return posX;}
   std::vector<float>&  get_posY(){return posY;}
   std::vector<float>&  get_posZ(){return posZ;}
   std::vector<float>&  get_posType(){return posType;}
   std::vector<float>&  get_posAge(){return posAge;}
   std::vector<float>&  get_posP(){return posP;}



 
   bool get_T1(){return T1;}
   bool get_T2(){return T2;}
   bool get_T3(){return T3;}
   bool get_T4(){return T4;}
   bool get_guard(){return guard;}
   std::vector<double>&  get_EneT1(){return EneT1;}
   std::vector<double>&  get_EneT2(){return EneT2;}
   std::vector<double>&  get_EneT3(){return EneT3;}
   std::vector<double>&  get_EneT4(){return EneT4;}
   std::vector<double>&  get_Eneg(){return Eneg;}
   std::vector<double>&  get_PHA6(){return PHA6;}
   std::vector<double>&  get_timeT1(){return timeT1;}
   std::vector<double>&  get_timeT2(){return timeT2;}
   std::vector<double>&  get_timeT3(){return timeT3;}
   std::vector<double>&  get_timeT4(){return timeT4;}
   std::vector<double>&  get_timeg(){return timeg;}
   int get_Ti(){return Ti;}
   std::vector<double>&  get_EneIsofoam(){return EneIsofoam;}
   std::vector<double>&  get_EneShell(){return EneShell;}
   std::vector<double>&  get_timeIsofoam(){return timeIsofoam;}
   std::vector<double>&  get_timeShell(){return timeShell;}
   
    //HOUSEKEEPING FROM COUNTERS 1 AND 3   
   int get_yCT1(){return yCT1;}
   int get_mCT1(){return mCT1;}
   int get_dCT1(){return dCT1;}
   int get_hCT1(){return hCT1;}
   int get_miCT1(){return miCT1;}
   int get_sCT1(){return sCT1;}
   float get_TempCT1(){return TempCT1;}
   int get_OnTimeCT1(){return OnTimeCT1;}
   int get_LastCT1(){return LastCT1;}
   int get_CountCT1(){return CountCT1;}
   float get_Baro1T(){return Baro1T;}
   float get_Baro1P(){return Baro1P;}
   float get_Baro2T(){return Baro2T;}
   float get_Baro2P(){return Baro2P;}
   float get_TempB1(){return TempB1;} 
   float get_TempB2(){return TempB2;} 
   float get_PressB1(){return PressB1;}
   float get_PressB2(){return PressB2;}
   float get_Volt5VCT1(){return Volt5VCT1;}  
   float get_Volt15VCT1(){return Volt15VCT1;} 
   float get_GOCT1(){return GOCT1;} 
   float get_coinCT1(){return coinCT1;} 

   int get_yCT3(){return yCT3;}
   int get_mCT3(){return mCT3;}
   int get_dCT3(){return dCT3;}
   int get_hCT3(){return hCT3;}
   int get_miCT3(){return miCT3;}
   int get_sCT3(){return sCT3;}
   float get_TempCT3(){return TempCT3;}
   int get_OnTimeCT3(){return OnTimeCT3;}
   int get_LastCT3(){return LastCT3;}
   int get_CountCT3(){return CountCT3;}
   float get_T1L(){return T1L;}
   float get_T1A(){return T1A;}
   float get_T2L(){return T2L;}
   float get_T2A(){return T2A;}
   float get_T3L(){return T3L;}
   float get_T3A(){return T3A;}
   float get_T4L(){return T4L;}
   float get_T4A(){return T4A;}
   float get_GRDL(){return GRDL;}
   float get_GRDA(){return GRDA;}
   float get_Volt5VCT3(){return Volt5VCT3;}  
   float get_Volt15VCT3(){return Volt15VCT3;} 

   int get_yPOW(){return yPOW;}
   int get_mPOW(){return mPOW;}
   int get_dPOW(){return dPOW;}
   int get_hPOW(){return hPOW;}
   int get_miPOW(){return miPOW;}
   int get_sPOW(){return sPOW;}
   int get_OnTimePOW(){return OnTimePOW;}

   float get_MainC(){return MainC;}
   float get_MainV(){return MainV;}
   float get_HeatC(){return HeatC;}
   float get_HeatV(){return HeatV;}
   float get_TrackC(){return TrackC;}
   float get_TrackV(){return TrackV;}

   ////////////////////////////////
   //Methods to get number of layerS and layer with hits
   ////////////////////////////////

    int  get_NLayers();
    int get_Layer(int);
    void get_Layers(int*);

   ////////////////////////////////
   ClassDef(ALEvent,1)

};

#endif



