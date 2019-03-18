////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

#include "ALEvent.h"
ClassImp(ALEvent)
ClassImp(ALTckhit)

//Constructors
ALEvent::ALEvent()// Default
{
 eventnumber=0; //Event number
 
 yPHA=-1;//Year from PHA line linked to the event
 mPHA=-1;//Month from PHA line linked to the event
 dPHA=-1;//Day from PHA line linked to the event
 hPHA=-1;//Hour from PHA line linked to the event
 miPHA=-1;//Minute from PHA line linked to the event
 sPHA=-1;//Second from PHA line linked to the event
 GoPHA=-1;//Go counter from PHA line linked to the event
 tPHA=-1;//timer from from PHA line linked to the event
 yEVT=-1;//Year from EVT line linked to the event
 mEVT=-1;//Month from EVT line linked to the event
 dEVT=-1;//Day from EVT line linked to the event
 hEVT=-1;//Hour from EVT line linked to the event
 miEVT=-1;//Minute from EVT line linked to the event
 sEVT=-1;//Second from EVT line linked to the event  
 EVT="";//Data from EVT line linked to the event
 GoEVT=-1;//Go counter from EVT line linked to the event
 tEVT=-1;//timer from from EVT line linked to the event
 nHitLEVT=-999;//nHitL from EVT line linked to the event added 12/05/2017
 CCEVT=-999;//CC from EVT line linked to the event   added 12/05/2017
 PatternEVT=-999;//CC from EVT line linked to the event   added 12/05/2017
 Q1EVT=-999;//Q1 from EVT line linked to the event   added 12/05/2017
 TrigEVT=-999;//Q1 from EVT line linked to the event   added 12/05/2017

 for(int i=0;i<7;i++) L[i]=string();//Data from  ASI lines of the event
 for(int i=0;i<7;i++) flagL[i]=0;////1 if ASI line was present

 //////////////////////////
 ///////MC variable////////
 //////////////////////////
	
 ncase=0; 
 typeMC=-99; //type of particle
 EkMC=0;   //kinetic energy of the particle
 pMC=0;		//momentum at point of injection
 X0MC=Y0MC=Z0MC=0;//Coordinates of the partcle at the injection point 
 CX0MC=CY0MC=CZ0MC=0; //Incidence cosines of the partcle at the injection point

 typePP=-99;
 EkPP=0;
 ZenPP=0;
 AziPP=0;
 CoLatSP=0;
 CoLonSP=0;
 	
 Nhits=0; //Number of hits in the event
 typereco=-999; //type of particle
 Ekreco=-999;   //kinetic energy of the particle
 p0reco=-999;  //momentum of the particle
 X0reco=Y0reco=Z0reco=0;//Coordinates of the partcle at the injection point 
 CX0reco=CY0reco=CZ0reco=0; //Incidence cosines of the partcle at the injection point 
 ndf=0;
 chi2=cl=-1;
 d0=phi0=cpa=dz=tanl=0;
 phi0_init=cpa_init=tanl_init=0;
 d0err2=phi0err2=cpaerr2=dzerr2=tanlerr2=0;
 Cov_init=Cov_last=0;
	
 EkPR=-999; 
 p0PR =-999;
 aPR=bPR=cPR=0;
 interPR=slopePR=0;
 chi2BPR=chi2NBPR=clBPR=clNBPR=-1;
 deflecPR=0;
 //Triggers: default is false
 T1=false;
 T2=false;
 T3=false;
 T4=false;
 guard=false;
 Ti=0;
 //HOUSEKEEPING FROM COUNTERS 1 AND 3
 //Data FROM "CT1" LINE
 yCT1=-1;//Year from CT1 line linked to the event (last read CT1 line)
 mCT1=-1;//Month from CT1 line linked to the event (last read CT1 line)
 dCT1=-1;//Day from CT1 line linked to the event (last read CT1 line)
 hCT1=-1;//Hour from CT1 line linked to the event (last read CT1 line)
 miCT1=-1;//Minute from CT1 line linked to the event (last read CT1 line)
 sCT1=-1;//Second from CT1 line linked to the event (last read CT1 line)
 TempCT1=-999; //Temperature measured on the board CT1

 OnTimeCT1=-1;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
 LastCT1=-1;//The last command received by the payload, expressed as a decimal number (it is in HeX on the GUI display)
 CountCT1=-1;//Count of commands received by the payload since power on
 
 //Barometer information: NOT INTERPRETED from line CT1
 Baro1T=-999;//Barometer 1 Temperature
 Baro1P=-999;//Barometer 1 Pressure
 Baro2T=-999;//Barometer 2 Temperature
 Baro2P=-999;//Barometer 2 Pressure
 //Barometer information: INTERPRETED from line CT1
 TempB1=-999;//Barometer 1 Temperature
 TempB2=-999;//Barometer 2 Temperature
 PressB1=-999;//Barometer 1 Pressure
 PressB2=-999;//Barometer 2 Pressure
 
 GOCT1=-999;
 coinCT1=-199;
 
 //Voltages
 Volt5VCT1=-999;  // Positive 5V from line CT1
 Volt15VCT1=-999; // Positive 15V from line CT1

 //Data FROM "CT3" LINE
 yCT3=-1;//Year from CT3 line linked to the event (last read CT3 line)
 mCT3=-1;//Month from CT3 line linked to the event (last read CT3 line)
 dCT3=-1;//Day from CT3 line linked to the event (last read CT3 line)
 hCT3=-1;//Hour from CT3 line linked to the event (last read CT3 line)
 miCT3=-1;//Minute from CT3 line linked to the event (last read CT3 line)
 sCT3=-1;//Second from CT3 line linked to the event (last read CT3 line)
 TempCT3=-999; //Temperature measured on the board CT3

 OnTimeCT3=-1;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
 LastCT3=-1;//The last command received by the payload, expressed as a decimal number (it is in HeX on the GUI display)
 CountCT3=-1;//Count of commands received by the payload since power on

 Volt5VCT3=-999;  // Positive 5V from line CT3
 Volt15VCT3=-999; // Positive 15V from line CT3
 
 //TRIGGER RATES (PHA AND LOGIC) from CT3
 T1L=-1;
 T1A=-1;
 T2L=-1;
 T2A=-1;
 T3L=-1;
 T3A=-1;
 T4L=-1;
 T4A=-1;
 GRDL=-1;
 GRDA=-1;
 
 //HOUSEKEEPING FROM POW
 int yPOW=-1;//Year from POW line linked to the event (last read POW line)
 int mPOW=-1;//Month from POW line linked to the event (last read POW line)
 int dPOW=-1;//Day from POW line linked to the event (last read POW line)
 int hPOW=-1;//Hour from POW line linked to the event (last read POW line)
 int miPOW=-1;//Minute from POW line linked to the event (last read POW line)
 int sPOW=-1;//Second from POW line linked to the event (last read POW line)
 int OnTimePOW=-1;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
 float MainC=-999;
 float MainV=-999;
 float HeatC=-999;
 float HeatV=-999;
 float TrackC=-999;
 float TrackV=-999;
}


void ALEvent::Copy(ALEvent* e)
{
  //Single variables 
  
  eventnumber =e->get_eventnumber();
  yPHA=e->get_yPHA();
  mPHA=e->get_mPHA();
  dPHA=e->get_dPHA();
  hPHA=e->get_hPHA();
  miPHA=e->get_miPHA();
  sPHA=e->get_sPHA();
  GoPHA=e->get_GoPHA();
  tPHA=e->get_tPHA();
  yEVT=e->get_yEVT();
  mEVT=e->get_mEVT();
  dEVT=e->get_dEVT();
  hEVT=e->get_hEVT();
  miEVT=e->get_miEVT();
  sEVT=e->get_sEVT();
  GoEVT=e->get_GoEVT();
  tEVT=e->get_tEVT();
  EVT=e->get_EVT();
  nHitLEVT=e->get_nHitLEVT();
  CCEVT=e->get_CCEVT();
  PatternEVT=e->get_PatternEVT();
  Q1EVT=e->get_Q1EVT();
  TrigEVT=e->get_TrigEVT();
  
  for(int i=0;i<7;i++)L[i]=e->get_L(i);
  for(int i=0;i<7;i++)flagL[i]=e->get_flagL(i);

   ncase =e->get_ncase();
   typeMC =e->get_typeMC();
   EkMC =e->get_EkMC();
   pMC = e->get_pMC();
   X0MC =e->get_X0MC();
   Y0MC =e->get_Y0MC();
   Z0MC =e->get_Z0MC();
   CX0MC =e->get_CX0MC();
   CY0MC =e->get_CY0MC();
   CZ0MC =e->get_CZ0MC();
   typePP=e->get_typePP();
   EkPP=e->get_EkPP();
   ZenPP=e->get_ZenPP();
   AziPP=e->get_AziPP();
   CoLatSP=e->get_CoLatSP();
   CoLonSP=e->get_CoLonSP();
   Nhits =e->get_Nhits();
   typereco =e->get_typereco();
   Ekreco =e->get_Ekreco();
   p0reco =e->get_p0reco();
   X0reco =e->get_X0reco();
   Y0reco =e->get_Y0reco();
   Z0reco =e->get_Z0reco();
   CX0reco =e->get_CX0reco();
   CY0reco =e->get_CY0reco();
   CZ0reco =e->get_CZ0reco();
   ndf =e->get_ndf();
   chi2 =e->get_chi2();
   cl =e->get_cl();
   d0 =e->get_d0();
   phi0 =e->get_phi0();
   cpa =e->get_cpa();
   dz =e->get_dz();
   tanl =e->get_tanl();
   phi0_init=get_phi0_init();
   cpa_init=get_cpa_init();
   tanl_init=get_tanl_init();
   d0err2 =e->get_d0err2();
   phi0err2 =e->get_phi0err2();
   cpaerr2 =e->get_cpaerr2();
   dzerr2 =e->get_dzerr2();
   tanlerr2 =e->get_tanlerr2();
   Cov_init = e->get_Cov_init();
   Cov_last = e->get_Cov_last();

   EkPR = e->get_EkPR();
   p0PR = e->get_p0PR();
   aPR = e->get_aPR();
   bPR = e->get_bPR();
   cPR = e->get_cPR();
   interPR = e->get_interPR();
   slopePR = e->get_slopePR();
   chi2BPR = e->get_chi2BPR();
   chi2NBPR = e->get_chi2NBPR();
   clBPR = e->get_clBPR();
   clNBPR = e->get_clNBPR();
   deflecPR = e->get_deflecPR();
   
   T1 =e->get_T1();
   T2 =e->get_T2();
   T3 =e->get_T3();
   T4 =e->get_T4();
   guard =e->get_guard();
   Ti =e->get_Ti();

   for(int i=0;i<(int)(e->get_posX()).size();i++) posX.push_back((e->get_posX()).at(i));
   for(int i=0;i<(int)(e->get_posY()).size();i++) posY.push_back((e->get_posY()).at(i));
   for(int i=0;i<(int)(e->get_posZ()).size();i++) posZ.push_back((e->get_posZ()).at(i));
   for(int i=0;i<(int)(e->get_posType()).size();i++) posType.push_back((e->get_posType()).at(i));
   for(int i=0;i<(int)(e->get_posAge()).size();i++) posAge.push_back((e->get_posAge()).at(i));
   for(int i=0;i<(int)(e->get_posP()).size();i++) posP.push_back((e->get_posP()).at(i));


   
   //Vectors of double
   for(int i=0;i<(int)(e->get_EneT1()).size();i++) EneT1.push_back((e->get_EneT1()).at(i));
   for(int i=0;i<(int)(e->get_EneT2()).size();i++) EneT2.push_back((e->get_EneT2()).at(i));
   for(int i=0;i<(int)(e->get_EneT3()).size();i++) EneT3.push_back((e->get_EneT3()).at(i));
   for(int i=0;i<(int)(e->get_EneT4()).size();i++) EneT4.push_back((e->get_EneT4()).at(i));
   for(int i=0;i<(int)(e->get_Eneg()).size();i++) Eneg.push_back((e->get_Eneg()).at(i));
   for(int i=0;i<(int)(e->get_PHA6()).size();i++) PHA6.push_back((e->get_PHA6()).at(i));
   for(int i=0;i<(int)(e->get_timeT1()).size();i++) timeT1.push_back((e->get_timeT1()).at(i));
   for(int i=0;i<(int)(e->get_timeT2()).size();i++) timeT2.push_back((e->get_timeT2()).at(i));
   for(int i=0;i<(int)(e->get_timeT3()).size();i++) timeT3.push_back((e->get_timeT3()).at(i));
   for(int i=0;i<(int)(e->get_timeT4()).size();i++) timeT4.push_back((e->get_timeT4()).at(i));
   for(int i=0;i<(int)(e->get_timeg()).size();i++) timeg.push_back((e->get_timeg()).at(i));
   for(int i=0;i<(int)(e->get_EneIsofoam()).size();i++) EneIsofoam.push_back((e->get_EneIsofoam()).at(i));
   for(int i=0;i<(int)(e->get_EneShell()).size();i++) EneShell.push_back((e->get_EneShell()).at(i));
   for(int i=0;i<(int)(e->get_timeIsofoam()).size();i++) timeIsofoam.push_back((e->get_timeIsofoam()).at(i));
   for(int i=0;i<(int)(e->get_timeShell()).size();i++) timeShell.push_back((e->get_timeShell()).at(i));
   //Vectors of ALTckhit
   for(int i=0;i<(int)(e->get_hits()).size();i++) hits.push_back((e->get_hits()).at(i));
  //HOUSEKEEPING FROM COUNTERS 1 AND 3   
   yCT1  = e->get_yCT1();
   mCT1  = e->get_mCT1();
   dCT1  = e->get_dCT1();
   hCT1  = e->get_hCT1();
   miCT1  = e->get_miCT1();
   sCT1  = e->get_sCT1();
   TempCT1  = e->get_TempCT1();
   OnTimeCT1  = e->get_OnTimeCT1();
   LastCT1  = e->get_LastCT1();
   CountCT1  = e->get_CountCT1();
   Baro1T  = e->get_Baro1T();
   Baro1P  = e->get_Baro1P();
   Baro2T  = e->get_Baro2T();
   Baro2P  = e->get_Baro2P();
   TempB1  = e->get_TempB1();
   TempB2  = e->get_TempB2();
   PressB1  = e->get_PressB1();
   PressB2  = e->get_PressB2();
   GOCT1=  e->get_GOCT1();
   coinCT1=  e->get_coinCT1();
 
   Volt5VCT1  = e->get_Volt5VCT1();
   Volt15VCT1  = e->get_Volt15VCT1();
   yCT3  = e->get_yCT3();
   mCT3  = e->get_mCT3();
   dCT3  = e->get_dCT3();
   hCT3  = e->get_hCT3();
   miCT3  = e->get_miCT3();
   sCT3  = e->get_sCT3();
   TempCT3  = e->get_TempCT3();
   OnTimeCT3  = e->get_OnTimeCT3();
   LastCT3  = e->get_LastCT3();
   CountCT3  = e->get_CountCT3();
   T1L  = e->get_T1L();
   T1A  = e->get_T1A();
   T2L  = e->get_T2L();
   T2A  = e->get_T2A();
   T3L  = e->get_T3L();
   T3A  = e->get_T3A();
   T4L  = e->get_T4L();
   T4A  = e->get_T4A();
   GRDL  = e->get_GRDL();
   GRDA  = e->get_GRDA();
   
   Volt5VCT3  = e->get_Volt5VCT3();
   Volt15VCT3  = e->get_Volt15VCT3();

   //HOUSEKEEPING FROM POW
   yPOW= e->get_yPOW();
   mPOW= e->get_mPOW();
   dPOW= e->get_dPOW();
   hPOW= e->get_hPOW();
   miPOW= e->get_miPOW();
   sPOW= e->get_sPOW();
   OnTimePOW= e->get_OnTimePOW();
   MainC= e->get_MainC();
   MainV= e->get_MainV();
   HeatC= e->get_HeatC();
   HeatV= e->get_HeatV();
   TrackC= e->get_TrackC();
   TrackV= e->get_TrackV();


 
}
////////////////////////////////
//Methods to get number of layerS and layer with hits
////////////////////////////////

int ALEvent::get_NLayers()
 {
  uint8_t tmpTi=(uint8_t)Ti;
 // cout << "get_NLayers tmpTi = " << unsigned(tmpTi) << endl;
  int NL=0;
  for(int ij=0;ij<7;ij++) NL+=(int)((tmpTi >>ij) & 0x01);
  return NL; 
 }
int ALEvent::get_Layer(int i)
 {
  uint8_t tmpTi=(uint8_t)Ti;
  int Ni=0;
  if(i<7) Ni=(int)((tmpTi >>i) & 0x01);
 // cout << " Layer " << i << " Ni = " << Ni << endl;
  return Ni; 
 }
void ALEvent::get_Layers(int*Lay)
 {
  uint8_t tmpTi=(uint8_t)Ti;
  // cout << "get_Layers tmpTi = " << unsigned(tmpTi) << endl;
  for(int ij=0;ij<7;ij++) Lay[ij]=(int)((tmpTi >>ij) & 0x01);
 }




