////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, May , 2017
////   Modififed by Sarah Mechbal, smechbal@ucsc.edu
////   Department of Physics, University of California Santa Cruz, December 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadMCparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void Scattering(int t, int nene, int ncycles, string s) ;


void Scattering(int t, int nene, int ncycles, string s)
{
	
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<2;i++)ShellReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
  double B = 0.3;		//average magnetic field in T
  double c = TMath::C();
 string MCparamfile="./MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);
 	
 
 
 int NTConfig=1;//Number of configuration of triggers
 string TConfig[1]={"T1&T4&T3"};
 //double TrigThresh[3]={0.3, 0.29, 0.57};		//trigger thresholds in %s determined for T1/T3/T4 
 

	
	
//Input file 
string Inppath="~/AESOPLITE/MCProduction/Detector/ProcessedFiles";	
 string startfile="aesopliteNonUniB_V4";
 string endfile="_fort.99";
 string directory= "~/AESOPLITE/Analysis/MCAnalysis/MC_Synthesis/V4";
string RecoInd;
 if(t==3 || t==4)	RecoInd= "";
 if(t==11||t==10)  RecoInd="";
 string TitleRecoID ="allPRevents";
 string RecoID = s;
	
//type of particle
int type=t;
string stype;	
if(t==3)	stype="e^{#minus}";
if(t==4)	stype="e^{#plus}";
if(t==11)	stype="#mu^{#minus}";
if(t==10)	stype="#mu^{#plus}";
if(t==1)	stype="protons";
if(t==6)	stype="#alpha";

 float mass=0.000511;//electon mass in GeV
 if(t==11 || t==10)	mass=0.10566;//muon mass in GeV
//Number of energies
 int Nene=nene;
//Energies
 int*Ene=new int[Nene];
if(t==3 || t==4 || t==1)
{  
 
 Ene[0]= 20;
 Ene[1]= 40;
 Ene[2]= 100;
 Ene[3]= 300;

 }
 string UNIT="MeV";	
 int Nbins=400;
 float EneMin=0;
 float EneMax=500;
 float Dmin=-2.;
 float Dmax=2.;
  
//muons	
 if(t==11||t==10)
  {
   Ene[0]=1;
   Ene[1]=2;
   Ene[2]=3;
   UNIT="GeV";
   Nbins=500;
   EneMin=500;
   EneMax=5000;
   Dmin=-0.05;
   Dmax=0.05;
  }	
 int mcolor[5]={1,2,4,8,33}; 
//zz0
float zz0=0;	
	


//Number of cycles per energy
 int* Ncycles=new int[Nene]; 
 for(int i=0;i<Nene;i++)Ncycles[i]=ncycles;	
	

//define histograms output

 double ymin = -20.0;				//y-bounds of T1
 double ymax = +20.0;				//y-bounds of T1
 double thetamin = -0.8;
 double thetamax = 0.8;

 //true deflection insRecoIDe tracker
 TH1F **deflectionNB=new TH1F*[Nene];//L6-L0
 TH1F **deflectionB=new TH1F*[Nene];//L5-L1

 //Deflection Pattern Recognition
 TH1F **deflectionBPR=new TH1F*[Nene];//L5-L1
 TH1F **deflectionNBPR=new TH1F*[Nene];//L6-L0 
 
 //Deflection from Reco
 TH1F **deflectionBReco=new TH1F*[Nene];//L5-L1
 TH1F **deflectionNBReco=new TH1F*[Nene];//L6-L0 
 
 ////Coordinates pull distributions 
 TH1F ****LayerPullPR=new TH1F***[Nene];
 TH1F ****LayerPullReco=new TH1F***[Nene];
	
	
 //Loss of energy between injection point and L0
 TH1F **Elossfront=new TH1F*[Nene];

 //MC inverse energy/momentum at L0
 TH1F **MomMC=new TH1F*[Nene];
 //PR energy - E at L0
 TH1F **EResoPR=new TH1F*[Nene];	
 TH1F **MomPR=new TH1F*[Nene];	
 TH1F **EPR=new TH1F*[Nene];	
 
 //Reco energy - E at L0
  TH1F **EResoReco=new TH1F*[Nene];	
 TH1F **MomReco=new TH1F*[Nene];	
 TH1F **EReco=new TH1F*[Nene];
	
 //PR & KF deflection scatter 

TH2F **ScatterDeflec=new TH2F*[Nene];


 for(int i=0; i<Nene; i++) 
   {
	if((t==3 || t==4) && Ene[i]>=60)Nbins=2*Nbins;
   	if((t==3 || t==4) && Ene[i]>=100)Nbins=3*Nbins;
 	deflectionNB[i]=new TH1F(Form("Scattering (truth), NB plane L6-L0, %d%s",Ene[i],UNIT.c_str()),Form("Scattering (truth), NB plane L6-L0, %d%s",Ene[i],UNIT.c_str()),Nbins,Dmin,Dmax);
    deflectionB[i]=new TH1F(Form("Deflection and Scattering (truth), B plane L5-L1, %d%s",Ene[i],UNIT.c_str()),Form("Deflection and Scattering (truth), B plane L5-L1, %d%s",Ene[i],UNIT.c_str()),Nbins,Dmin,Dmax);
	deflectionNB[i]->GetXaxis()->SetTitle("#Delta#theta=#theta_{L6}-#theta_{L0} (rad)");   
	deflectionB[i]->GetXaxis()->SetTitle("#Delta#theta=#theta_{L5}-#theta_{L1} (rad)");   
	deflectionNB[i]->GetYaxis()->SetTitle("arbitrary unit");   
	deflectionB[i]->GetYaxis()->SetTitle("arbitrary unit");	   
	
	deflectionBPR[i]=new TH1F(Form("Deflection and Scattering (PR), B plane L5-L1, %d%s",Ene[i],UNIT.c_str()),Form("Deflection and Scattering (PR), B plane L5-L1, %d%s",Ene[i],UNIT.c_str()),Nbins,Dmin,Dmax);
	deflectionBPR[i]->GetXaxis()->SetTitle("#Delta#theta=#theta_{L5}-#theta_{L1} (rad)");   
	deflectionBPR[i]->GetYaxis()->SetTitle("arbitrary unit");	   
	
	deflectionNBPR[i]=new TH1F(Form("Deflection and Scattering (PR), NB plane L6-L0, %d%s",Ene[i],UNIT.c_str()),Form("Deflection and Scattering (PR), NB plane L6-L0, %d%s",Ene[i],UNIT.c_str()),Nbins,Dmin,Dmax);
	deflectionNBPR[i]->GetXaxis()->SetTitle("#Delta#theta=#theta_{L6}-#theta_{L0} (rad)");   
	deflectionNBPR[i]->GetYaxis()->SetTitle("arbitrary unit");	

	deflectionBReco[i]=new TH1F(Form("Deflection and Scattering (Reco), B plane L5-L1, %d%s",Ene[i],UNIT.c_str()),Form("Deflection and Scattering (Reco), B plane L5-L1, %d%s",Ene[i],UNIT.c_str()),Nbins,Dmin,Dmax);   	
	deflectionNBReco[i]=new TH1F(Form("Deflection and Scattering (Reco), NB plane L6-L0, %d%s",Ene[i],UNIT.c_str()),Form("Deflection and Scattering (Reco), NB plane L6-L0, %d%s",Ene[i],UNIT.c_str()),Nbins,Dmin,Dmax);

	Elossfront[i] =new TH1F(Form("Eloss (truth), Injection Point -L0, %d%s",Ene[i],UNIT.c_str()),Form("Eloss (truth), Injection Point -L0, %d%s",Ene[i],UNIT.c_str()),100,0.,10.); 
    Elossfront[i]->GetXaxis()->SetTitle("#DeltaE=E_{IP}-E_{L0} (MeV)");   
   	Elossfront[i]->GetYaxis()->SetTitle("arbitrary unit");	   
    EResoPR[i] =new TH1F(Form("EL0 - EPR %d%s",Ene[i],UNIT.c_str()),Form("EL0 - EPR %d%s",Ene[i],UNIT.c_str()),2000,-1000.,1000.); 
    EResoPR[i]->GetXaxis()->SetTitle("#DeltaE=E_{L0}-E_{PR} (MeV)");   
   	EResoPR[i]->GetYaxis()->SetTitle("arbitrary unit");	
	EResoReco[i] =new TH1F(Form("EL0 - EPR & EL0 - E_{Reco} %d%s",Ene[i],UNIT.c_str()),Form("EL0 - E_{Reco} %d%s",Ene[i],UNIT.c_str()),2000,-1000.,1000.); 
    MomPR[i] =new TH1F(Form(" 1/p, %d%s",Ene[i],UNIT.c_str()),Form("1/p, %d%s",Ene[i],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
    MomPR[i]->GetXaxis()->SetTitle(Form("Inverse Momentum PR & Reco (%s^{-1})",UNIT.c_str()));   
   	MomPR[i]->GetYaxis()->SetTitle("#entries");	
    MomReco[i] =new TH1F(Form("Momentum from Reco, %d%s",Ene[i],UNIT.c_str()),Form("Momentum from Reco, %d%s",Ene[i],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
    MomMC[i] =new TH1F(Form("Momentum from MC, %d%s",Ene[i],UNIT.c_str()),Form("Momentum from MC, %d%s",Ene[i],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
    EPR[i] =new TH1F(Form("Energy from PR, %d%s",Ene[i],UNIT.c_str()),Form("Energy from PR, %d%s",Ene[i],UNIT.c_str()),1000,EneMin,EneMax); 
    EPR[i]->GetXaxis()->SetTitle(Form("Energy from PR & Reco (%s)",UNIT.c_str()));   
    EPR[i]->GetYaxis()->SetTitle("arbitrary unit");	
    EReco[i] =new TH1F(Form("Energy from Reco, %d%s",Ene[i],UNIT.c_str()),Form("Energy from Reco, %d%s",Ene[i],UNIT.c_str()),1000,EneMin,EneMax); 

// TH2F
   ScatterDeflec[i] = new TH2F(Form("Scatter deflection, %d%s",Ene[i],UNIT.c_str()),Form("Scatter deflection %d%s",Ene[i],UNIT.c_str()),500,-1,1,500,-1,1);
   ScatterDeflec[i]->SetTitle(Form("Deflection in bending plane scatter %d%s",Ene[i],UNIT.c_str()));
    } // end i 
   //Pull distrubutions
   for(int i=0;i<Nene;i++) {


    LayerPullPR[i]=new TH1F**[7];
	LayerPullReco[i]=new TH1F**[7];
       for(int j=0;j<7;j++)
      {
       LayerPullPR[i][j]= new TH1F*[3];
       LayerPullReco[i][j]= new TH1F*[3];
		  for(int k=0;k<3;k++) 
		 {
		 LayerPullPR[i][j][k] = new TH1F(Form("Pulls from PR L%d,%d MeV, coord %d",j+1,Ene[i],k),Form("Pulls from PR L%d,%d MeV, coord %d",j+1,Ene[i], k),1000,-10,10); 
		 LayerPullReco[i][j][k] = new TH1F(Form("Pulls from Reco L%d,%d MeV, coord %d",j+1,Ene[i],k),Form("Pulls from reco L%d,%d MeV, coord %d",j+1,Ene[i], k),1000,-10,10); 
     		 } //end k
	  } //end j
   } //end i
	
//cout << "chaining events" << endl;	
	
//Create TChain to merge the files for a given energy
	
 TChain **chain=new TChain*[Nene]; 

//output ROOT file
 TFile*fileout=new TFile(Form("%s/%d/PR_Analysis_%s_%d_%s.root", directory.c_str(),type,startfile.c_str(), type,RecoID.c_str()),"RECREATE");

	
for(int i=0; i<Nene; i++) 
  {  	
   cout << "Energy: " << Ene[i] << UNIT <<endl;
   chain[i]=new TChain("MC");
   for(int j=0;j<Ncycles[i];j++)//Number of cycles
     { 
	chain[i]->Add(Form("%s/%d/RecoEvent_%s_%d_%d%s%s_%s.root", Inppath.c_str(),type,startfile.c_str(),type,Ene[i],UNIT.c_str(),endfile.c_str(),RecoID.c_str()));

     }

   //Define variables to read event
   ALEvent *e = new ALEvent();      
   //Set address to access event data
   chain[i]->SetBranchAddress("Revent",&e); 
 
   // Get number of event in Tree
   int nentries=chain[i]->GetEntries();
   cout << "Number  of events: " << nentries << endl;  
      
   for (int j=0;j<nentries;j++)
     {
      chain[i]->GetEntry(j); //Load the entry i in the variable e 
      bool* t=new bool[NTConfig];
      for(int k=0;k<NTConfig;k++) t[k]=true;
      if(j%1000000==0) cout << "Event: " << j <<endl;
	  uint8_t Ti=(uint8_t)e->get_Ti();
	  //Number of layers with hit(s)
         int NL= e->get_NLayers();

 
		 
      //Fill total energy deposited in the scintillators
	  //Internal trigger requirement: at least 5 hits, event chosen by PR
	  if((e->get_deflecPR()==0)) continue;
		//apply external trigger requirements
	bool T1=e->get_T1();
	bool T3=e->get_T3();
    bool T4=e->get_T4();
	 //Number of layers with hit(s) in bending/non-bending plane
	 int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	 int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
         if(NLNB==2 || NLB==3) continue;
		if(!T1 || !T4 || !T3) continue;
		 
		  //"truth" directional informations at entrance point T1 
		  double cxin = e->get_CX0MC();
		  double cyin = e->get_CY0MC();
		  double czin = e->get_CZ0MC();
		  double X0 = e->get_X0MC();
		  double Y0 = e->get_Y0MC();
		  double Z0 = e->get_Z0MC();
		  double slopeT1B = (cyin/czin);
		  double thetaT1B = TMath::ATan(slopeT1B); 
		  double slopeT1NB = (cxin/czin);
		  double thetaT1NB = TMath::ATan(slopeT1NB);
	  	  // bending plane PR variables at L1 and L5 
		  double a = e->get_aPR(); 
		  double b = e->get_bPR(); 
		  double c = e->get_cPR(); 
		  float limL0=TckZPos[0];
		  float limL6=TckZPos[6];
		  float limL5=TckZPos[5];  		//z position of last bending plane layer
   		  float limL1=TckZPos[1];		    //z position of first bending plane layer
		  double slopeL1PR = 2*a*limL1;
		  double slopeL5PR = 2*a*limL5;
		  double thetaL1PR = TMath::ATan(slopeL1PR);
		  double thetaL5PR = TMath::ATan(slopeL5PR);
		  double deflectionPR = e->get_deflecPR();
		  double deflecPRTheta = thetaL5PR-thetaL1PR;
		  float P=fabs(e->get_p0PR());   //in GeV
		  float E=e->get_EkPR();	  
		  double Ekreco = e->get_Ekreco();			//in MeV
		  double pReco = e->get_p0reco();
		  MomPR[i]->Fill(1/(1000.*P));			//in MeV
	          EPR[i]->Fill(1000.*E);
		  if(Ekreco==-999) continue;	 
		  EReco[i]->Fill(Ekreco);
		  MomReco[i]->Fill(1/pReco);
		  deflectionBPR[i]->Fill(deflecPRTheta);
		  double slopeMCL0=0; 
		  double thetaMCL0=0; 
		  double slopePRL0=0;
		  double thetaPRL0=0;
		  double slopeRecoL0=0;
		  double thetaRecoL0=0;
		  double slopeMCL1=0; 
		  double thetaMCL1=0; 
		  double slopeRecoL1=0;
		  double thetaRecoL1=0;
		  double slopeMCL5=0; 
		  double thetaMCL5=0; 
		  double slopeRecoL5=0;
		  double thetaRecoL5=0;
		  double slopeMCL6=0; 
		  double thetaMCL6=0; 
		  double slopePRL6=0;
		  double thetaPRL6=0;
	          double slopeRecoL6=0;
		  double thetaRecoL6=0;
		  double phi0 = e->get_phi0();
		  double tanl = e->get_tanl();
		  double cpa = e->get_cpa();
		  double slope_Reco = 1/(tanl);
		  double a_Reco = c*B*cpa/2;
		  double thetaNB = TMath::ATan(e->get_slopePR());
 		  double thetaReco = TMath::ACos(1/(TMath::Sqrt(1+(tanl*tanl))));
		//  cout << "aPR = " << a << "  aReco = " << a_Reco << endl;
		//  cout << "slopePR = " << e->get_slopePR() << "  slopeReco = " << slope_Reco << endl;
		//  cout << " thetaNB = " << thetaNB << "thetaReco = " << thetaReco << endl; 
  
		  for(int k=0;k<e->get_Nhits();k++) 
		    {
			 int Lindex=e->get_hits().at(k)->get_L();
			 if((int)e->get_hits().at(k)->get_flagPR()==1) {
		           double cxPR = e->get_hits().at(k)->get_cxPR();
			   double cyPR = e->get_hits().at(k)->get_cyPR();
			   double czPR = e->get_hits().at(k)->get_czPR();
			 //  cout << "cxPR = " << cxPR << " cyPR = " << cyPR << " czPR = " << czPR << endl;
			 //LAYER L0, NB plane		
			 if(Lindex==0) 
			  {								
			   //truth variables at L0
			   double cxL0 = e->get_hits().at(k)->get_cx();
			   double cyL0 = e->get_hits().at(k)->get_cy();
			   double czL0 = e->get_hits().at(k)->get_cz();
			   slopeMCL0= (cxL0/czL0);
			   thetaMCL0 = TMath::ATan(slopeMCL0);
			   Elossfront[i]->Fill(1000*(e->get_EkMC()+mass - e->get_hits().at(k)->get_eMC()));	
			   EResoPR[i]->Fill(1000.*(e->get_hits().at(k)->get_eMC()) - E);
			   double eMC = 1000.*(e->get_hits().at(k)->get_eMC());		//energy at L0 in MeV
			   double pMC = TMath::Sqrt((eMC*eMC) - (mass*mass));
			   MomMC[i]->Fill(1/pMC);
			   if(Ekreco ==-999) continue;
			    EResoReco[i]->Fill((1000.*(e->get_hits().at(k)->get_eMC())) - Ekreco);
			
			   //PR variables at L0
			   double cxL0PR = e->get_hits().at(k)->get_cxPR();
			   double cyL0PR = e->get_hits().at(k)->get_cyPR();
			   double czL0PR = e->get_hits().at(k)->get_czPR();
			   slopePRL0= (cxL0PR/czL0PR);
			   thetaPRL0 = TMath::ATan(slopePRL0); 	
			 //  cout << " cxL0PR = " << cxL0PR << " czL0PR = " << czL0PR << "  thetaPRL0= " << thetaPRL0 << endl;
			  //Reco variables at L6
			   if((int)e->get_hits().at(k)->get_fUsed()==1) {
			   double cxL0Reco = e->get_hits().at(k)->get_cxreco();
			   double cyL0Reco = e->get_hits().at(k)->get_cyreco();
			   double czL0Reco = e->get_hits().at(k)->get_czreco();
			   slopeRecoL0= (cxL0Reco/czL0Reco);
			   thetaRecoL0 = TMath::ATan(slopeRecoL0); 
			     }	
		   	  }
				
             //LAYER L1, B plane		
		     if(Lindex==1)
		      {								
			   //truth variables at L1
			   double cxL1 = e->get_hits().at(k)->get_cx();
			   double cyL1 = e->get_hits().at(k)->get_cy();
		           double czL1 = e->get_hits().at(k)->get_cz();
			   slopeMCL1= (cyL1/czL1);
			   thetaMCL1 = TMath::ATan(slopeMCL1);	
			   //Reco variables at L1
			   if((int)e->get_hits().at(k)->get_fUsed()==1){
			   double cxL1Reco = e->get_hits().at(k)->get_cxreco();
			   double cyL1Reco = e->get_hits().at(k)->get_cyreco();
			   double czL1Reco = e->get_hits().at(k)->get_czreco();
			   slopeRecoL1= (cyL1Reco/czL1Reco);
			   thetaRecoL1 = TMath::ATan(slopeRecoL1); 	
			    }
		      }	
		     //LAYER L5, B PLANE
		     if(Lindex==5)
		      {						
			   //truth variables at L5
			   double cyout = e->get_hits().at(k)->get_cy();
			   double czout = e->get_hits().at(k)->get_cz();
			   double cxL5PR = e->get_hits().at(k)->get_cxPR();
			   double czL5PR = e->get_hits().at(k)->get_czPR();
			//   cout << "cxL5PR = " << cxL5PR << "  czL5PR = " << czL5PR << endl;
			   slopeMCL5= (cyout/czout);
			   thetaMCL5 = TMath::ATan(slopeMCL5);
			   double deflectionMC = slopeMCL5 - slopeT1B;
                           double deflecMCTheta = thetaMCL5 - thetaT1B;
			   //Reco variables at L5
			   	if((int)e->get_hits().at(k)->get_fUsed()==1) {
			   double cxL5Reco = e->get_hits().at(k)->get_cxreco();
			   double cyL5Reco = e->get_hits().at(k)->get_cyreco();
			   double czL5Reco = e->get_hits().at(k)->get_czreco();
			   slopeRecoL5= (cyL5Reco/czL5Reco);
			   thetaRecoL5 = TMath::ATan(slopeRecoL5); 	
			   }
		      }
		     //NON-BENDING PLANE
		     if(Lindex==6)
		      {
			   //truth variables
			   double cxout = e->get_hits().at(k)->get_cx();
			   double czout = e->get_hits().at(k)->get_cz();
			   double cxL6PR = e->get_hits().at(k)->get_cxPR();
			   double czL6PR = e->get_hits().at(k)->get_czPR();     
			   slopeMCL6= (cxout/czout);
			   thetaMCL6 = TMath::ATan(slopeMCL6);
			   slopePRL6 = (cxL6PR/czL6PR);
			   thetaPRL6 = TMath::ATan(slopePRL6);
			   //Reco variables at L6
			//   cout << "cxL6PR = " << cxL6PR << "  czL6PR = " << czL6PR << " thetaPRL6 = " << thetaPRL6 << endl;
			   if((int)e->get_hits().at(k)->get_fUsed()==1) {
			   double cxL6Reco = e->get_hits().at(k)->get_cxreco();
			   double cyL6Reco = e->get_hits().at(k)->get_cyreco();
			   double czL6Reco = e->get_hits().at(k)->get_czreco();
			   slopeRecoL6= (cxL6Reco/czL6Reco);
			   thetaRecoL6 = TMath::ATan(slopeRecoL6); 
			   	}
			  } //end if L6
			 }// end if flagPR
		    }	//end for

		   //Fill histograms	 
		   deflectionNB[i]->Fill(thetaMCL6-thetaMCL0); 	 
		   deflectionB[i]->Fill(thetaMCL5-thetaMCL1);
		 //  cout << "filling deflectionNBPR = " << thetaPRL6 - thetaPRL0 << endl; 
		 //  deflectionNBPR[i]->Fill(thetaPRL6-thetaPRL0); 
		   if(thetaRecoL5!=0 && thetaRecoL1!=0) {
		
		   deflectionBReco[i]->Fill(thetaRecoL5-thetaRecoL1);
                   ScatterDeflec[i]->Fill(deflecPRTheta,(thetaRecoL5-thetaRecoL1));
		   }
		   if(thetaRecoL6!=0 && thetaRecoL0!=0) 
			{
		   deflectionNBReco[i]->Fill(thetaRecoL6-thetaRecoL0);
			}
		
		//Fill pull distributions histograms in mm
		for(int k=0;k<e->get_Nhits();k++)
         		{
          if((e->get_hits().at(k)->get_flagPR()) && (!e->get_hits().at(k)->get_fGhost()) )       //Check that the hit was used for PR
           {
	    float tmp=0;
	    int Lindex=(int)e->get_hits().at(k)->get_L();
	    tmp=(e->get_hits().at(k)->get_x())*10;
	    LayerPullPR[i][Lindex][0]->Fill(tmp-(e->get_hits().at(k)->get_xPR())*10);
	    tmp=(e->get_hits().at(k)->get_y())*10;
	    LayerPullPR[i][Lindex][1]->Fill(tmp-(e->get_hits().at(k)->get_yPR())*10);
	    tmp=(e->get_hits().at(k)->get_z())*10;
	    LayerPullPR[i][Lindex][2]->Fill(tmp-(e->get_hits().at(k)->get_zPR())*10);	
	    } //end if flagPR
         		
         if((e->get_hits().at(k)->get_fUsed()) && (!e->get_hits().at(k)->get_fGhost()))       //Check that the hit was used for reconstruction
           {
	    float tmp=0;
	    int Lindex=(int)e->get_hits().at(k)->get_L();
	    tmp=(e->get_hits().at(k)->get_x())*10;
	    LayerPullReco[i][Lindex][0]->Fill(tmp-e->get_hits().at(k)->get_xreco());
	    tmp=(e->get_hits().at(k)->get_y())*10;
	    LayerPullReco[i][Lindex][1]->Fill(tmp-e->get_hits().at(k)->get_yreco());
	    tmp=(e->get_hits().at(k)->get_z())*10;
	    LayerPullReco[i][Lindex][2]->Fill(tmp-e->get_hits().at(k)->get_zreco());	
	    } //end if fUsed
		
	  } //end k

     }  //end j
  } //end i

	   
//////////////////////////////   
// Display  
//////////////////////////////   
 
 gStyle->SetOptStat(1);	

//Non-Bending plane	
 TCanvas*canNB=new TCanvas("deflections in NB","deflections in NB",200,10,1500,1000);
 canNB->Divide(2,2);
 TLegend**legNB=new TLegend*[Nene];
 TF1***fDeflecNB=new TF1**[Nene];
 TPaveText**PTNB=new TPaveText*[Nene];

	
 for(int i=0;i<Nene;i++) 
   {
   
       fDeflecNB[i] = new TF1*[2];
	canNB->cd(i+1);
	double MaxYNB=0;
        double tmp=0;
	tmp = deflectionNB[i]->GetMaximum();
	double tmpReco=0;
	tmpReco=deflectionNBReco[i]->GetMaximum();
	if(tmp>tmpReco)  MaxYNB=tmp;
	else MaxYNB=tmpReco;	
	//deflectionNB[i]->SetMaximum(1.05*MaxYNB);
	if(t==11||t==10)deflectionNB[i]->GetXaxis()->SetRangeUser(-0.01,0.01);
        deflectionNB[i]->SetTitle(Form("Scattering in non bending plane, L6-L0 %d%s", Ene[i],UNIT.c_str()));
	deflectionNB[i]->SetLineColor(mcolor[0]);
        deflectionNB[i]->Scale(MaxYNB/deflectionNB[i]->GetMaximum());
        deflectionNB[i]->Draw("hist");		
	deflectionNBReco[i]->SetLineColor(mcolor[2]);
	//deflectionNBReco[i]->Scale(MaxYNB/deflectionNBReco[i]->GetMaximum());
	deflectionNBReco[i]->Draw("same");
	legNB[i] = new TLegend(0.1,0.8,0.3,0.6);
	legNB[i]->AddEntry(deflectionNB[i],Form("MC truth %s,%d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");
	legNB[i]->AddEntry(deflectionNBReco[i],Form("Reco %s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");
        legNB[i]->SetBorderSize(0);
        legNB[i]->SetFillStyle(0);	 
	legNB[i]->Draw("same");
	for(int j=0;j<2;j++) 
	{
	if(t==3 || t==4)PTNB[i]=new TPaveText(1.0,MaxYNB-(j+1)*0.3*MaxYNB,2.,MaxYNB-j*0.1*MaxYNB);
	if(t==11 || t==10)PTNB[i]=new TPaveText(0.003,MaxYNB-(j+1)*0.1*MaxYNB,0.01,MaxYNB-j*0.1*MaxYNB); 
	if(j==0){
		fDeflecNB[i][j]=new TF1(Form("fDeflecNB%d, %d",i,j),"gaus",-2,1);
		deflectionNB[i]->Fit(fDeflecNB[i][j],"0","",-0.4,0.4);
		fDeflecNB[i][j]->SetLineColor(mcolor[0]);
	}

	if(j==1) {
 	    fDeflecNB[i][j]=new TF1(Form("fDeflecNB%d, %d",i,j),"gaus",-2,1);
            deflectionNBReco[i]->Fit(fDeflecNB[i][j],"0","",-0.4,0.4);
	    fDeflecNB[i][j]->SetLineColor(mcolor[2]);
	}
	 PTNB[i]->SetFillStyle(0);
	PTNB[i]->SetBorderSize(0);
    	    PTNB[i]->AddText(Form("%s, %d%s, #mu=%5.4f mrad, #sigma= %5.4f mrad",stype.c_str(),Ene[i],UNIT.c_str(),1000.*fDeflecNB[i][j]->GetParameter(1),1000.*fDeflecNB[i][j]->GetParameter(2))); 
	//fDeflecNB[i][j]->Draw("same");
	PTNB[i]->Draw("same");
	} //end j
	fileout->cd();
	deflectionNB[i]->Write();
	deflectionNBReco[i]->Write();
	canNB->Write();

   }	   

   //  canNB->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf(", directory.c_str(),type, startfile.c_str(), type,RecoID.c_str()));
cout << "printing first canvas" << endl;
		
 //Bending plane
 
 TCanvas*canB=new TCanvas("deflections in B","deflections in B",200,10,1500,1000);
 canB->Divide(2,2);
 TLegend**legB=new TLegend*[Nene];
 TF1***fDeflecB=new TF1**[Nene];
 TPaveText**PTB=new TPaveText*[Nene];


	
 for(int i=0;i<Nene;i++) 
   {
    fDeflecB[i] = new TF1*[3];
    canB->cd(i+1);
    double tmp=deflectionB[i]->GetMaximum();
    double tmpPR=deflectionBPR[i]->GetMaximum();
    double tmpReco=deflectionBReco[i]->GetMaximum();
    double MaxYB=0;
	if(tmp>tmpPR && tmp>tmpReco)  MaxYB=tmp; 
	else if(tmpPR>tmpReco) MaxYB=tmpPR;
	else MaxYB=tmpReco;
	deflectionB[i]->SetLineColor(mcolor[0]);
	deflectionB[i]->SetMaximum(1.10*MaxYB);
	if(t==11||t==10)deflectionB[i]->GetXaxis()->SetRangeUser(-1,0.1);
        deflectionB[i]->SetTitle(Form("Scattering in bending plane, L5-L1 %d%s", Ene[i],UNIT.c_str()));
	//deflectionB[i]->Scale(MaxYB/deflectionB[i]->GetMaximum());
	deflectionB[i]->Draw("hist");	
	deflectionBPR[i]->SetLineColor(mcolor[1]);
	//deflectionBPR[i]->Scale(MaxYB/deflectionBPR[i]->GetMaximum());
	deflectionBPR[i]->Draw("sames");
	//deflectionBReco[i]->Scale(MaxYB/deflectionBReco[i]->GetMaximum());
	deflectionBReco[i]->Draw("sames");
	deflectionBReco[i]->SetLineColor(mcolor[2]);
	legB[i] = new TLegend(0.2,0.5,0.3,0.8);
	legB[i]->AddEntry(deflectionB[i],Form("MC truth %s,%d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");
	legB[i]->AddEntry(deflectionBPR[i],Form("PR %s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");
	legB[i]->AddEntry(deflectionBReco[i],Form("Reco %s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");
        legB[i]->SetBorderSize(0);
        legB[i]->SetFillStyle(0);	
	legB[i]->Draw("same");
	for(int j=0;j<3;j++) 
	{
 	//fDeflecB[i][j]=new TF1(Form("fDeflecB%d, %d",i,j),"gaus",-2,1);
	if(j==0){
 	    float xMax=deflectionB[i]->GetXaxis()->GetBinCenter(deflectionB[i]->GetMaximumBin());
	    fDeflecB[i][j]=new TF1(Form("fDeflecB%d %d",i,j),"gaus",xMax-fabs(0.2*xMax),xMax+fabs(0.2*xMax));
	   deflectionB[i]->Fit(fDeflecB[i][j],"R0");
		}
	if(j==1) {
 	    float xMax=deflectionBPR[i]->GetXaxis()->GetBinCenter(deflectionBPR[i]->GetMaximumBin());
	    fDeflecB[i][j]=new TF1(Form("fDeflecB%d %d",i,j),"gaus",xMax-fabs(0.2*xMax),xMax+fabs(0.2*xMax));
	    deflectionBPR[i]->Fit(fDeflecB[i][j],"R0");
		}
	if(j==2) {
 	    float xMax=deflectionBReco[i]->GetXaxis()->GetBinCenter(deflectionBReco[i]->GetMaximumBin());
	    fDeflecB[i][j]=new TF1(Form("fDeflecB%d %d",i,j),"gaus",xMax-fabs(0.2*xMax),xMax+fabs(0.2*xMax));
	    deflectionBReco[i]->Fit(fDeflecB[i][j],"R0");
	}
	fDeflecB[i][j]->SetLineColor(mcolor[j]);
	fDeflecB[i][j]->Draw("same"); 
	if(t==3||t==4)PTB[i]=new TPaveText(1.0,MaxYB-(j+1)*0.3*MaxYB,2.,MaxYB-j*0.1*MaxYB);
	if(t==11)PTB[i]=new TPaveText(-0.019,MaxYB-(j+1)*0.1*MaxYB,-0.008,MaxYB-j*0.1*MaxYB);
	if(t==10)PTB[i]=new TPaveText(0.008,MaxYB-(j+1)*0.1*MaxYB,0.019,MaxYB-j*0.1*MaxYB);
    PTB[i]->AddText(Form("%s, %d%s, #mu=%5.0f mrad, #sigma= %5.0f mrad",stype.c_str(),Ene[i],UNIT.c_str(),1000.*fDeflecB[i][j]->GetParameter(1),1000.*fDeflecB[i][j]->GetParameter(2)));  
    PTB[i]->SetFillStyle(0);PTB[i]->SetBorderSize(0);
	PTB[i]->Draw("same");
	} //end j 
	fileout->cd();
	deflectionB[i]->Write();
	deflectionBPR[i]->Write();
	deflectionBReco[i]->Write();
   }	//end i	   

	canB->Write();
//     canB->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf", directory.c_str(),type, startfile.c_str(), type,RecoID.c_str()));
	ofstream toTeX;
	toTeX.open("TestToTeX.txt");
	toTeX << "Hello World " << endl;
////////////////////////////
//KF and PR deflection PR//
///////////////////////////

  TCanvas *CanScatter = new TCanvas("Deflection scatter", "deflection scatter",200,10,1500,1000);
  CanScatter->Divide(2,2);
  for(int i=0;i<Nene;i++) {
	CanScatter->cd(i+1);
	ScatterDeflec[i]->Draw("colz");
	 TF1*yx=new TF1("yx","x",-0.5,0.5);
 	 yx->SetLineColor(kGray);
	 yx->SetLineStyle(6);
	yx->Draw("same");
        ScatterDeflec[i]->GetXaxis()->SetRangeUser(-0.5,0.5);
	ScatterDeflec[i]->GetYaxis()->SetRangeUser(-0.5,0.5);
	ScatterDeflec[i]->GetXaxis()->SetTitle("#Delta#theta_{PR}");
	ScatterDeflec[i]->GetYaxis()->SetTitle("#Delta#theta_{KF}");
	fileout->cd();
	ScatterDeflec[i]->Write();
 	}


 CanScatter->Update();
 //CanScatter->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf", directory.c_str(),type, startfile.c_str(), type,RecoID.c_str()));


  ////////////////////////////////////////////////////////
//////	Momentum reconstruction ////////////////////////
////////////////////////////////////////////////////////

 TCanvas*can3=new TCanvas("Momentum PR","Momentum PR",200,10,1500,1000);
 can3->Divide(2,2);	
 TLegend**leg3=new TLegend*[Nene];
 TF1***fMom=new TF1**[Nene];
 TPaveText**PTMom=new TPaveText*[Nene];
 TLine**MomLine=new TLine*[Nene];
 TMultiGraph* preco= new TMultiGraph();
 TMultiGraph* resoreco= new TMultiGraph();
 preco->SetTitle("Inverse momentum recontruction PR (red) and KF (blue) vs truth momentum");
 resoreco->SetTitle("Momentum recontruction resolution PR (red) and KF (blue) vs truth momentum");


 TGraphErrors* PRreco = new TGraphErrors();
 PRreco->SetLineColor(kRed);
 PRreco->SetMarkerStyle(4);
 PRreco->SetMarkerColor(kRed);
 TGraphErrors* KFreco = new TGraphErrors();
 KFreco->SetLineColor(kBlue);
 KFreco->SetMarkerStyle(4);
 KFreco->SetMarkerColor(kBlue);
 TGraphErrors* PRreso = new TGraphErrors();
 PRreso->SetLineColor(kRed);
 PRreso->SetMarkerStyle(4);
 PRreso->SetMarkerColor(kRed);
 TGraphErrors* KFreso = new TGraphErrors();
 KFreso->SetLineColor(kBlue);
 KFreso->SetMarkerStyle(4);
 KFreso->SetMarkerColor(kBlue);
 int sigma=4;
 gStyle->SetOptStat(01);
 for(int i=0;i<Nene;i++) {
    can3->cd(i+1);
	leg3[i]=new TLegend(0.2,0.5,0.4,0.9);
	leg3[i]->SetBorderSize(0);
	leg3[i]->SetFillStyle(0);
	fMom[i]=new TF1*[2];
	double MaxYMom=0;
        double tmp=MomPR[i]->GetMaximum();
	double tmpReco=MomReco[i]->GetMaximum();
	double tmpMC= MomMC[i]->GetMaximum();
	if(tmp> tmpReco)  MaxYMom=tmp; 
//        MaxYMom=tmpMC;
	else MaxYMom=tmpReco;
	//MomLine[i]->SetLineWRecoIDth(2);
	MomPR[i]->SetTitle(Form("Inverse momentum PR & KF, %s %d%s",stype.c_str(),Ene[i],UNIT.c_str()));  
	MomPR[i]->Draw("hist");
	MomPR[i]->GetXaxis()->SetRangeUser(MomReco[i]->GetMean()-sigma*MomReco[i]->GetRMS(),MomReco[i]->GetMean()+sigma*MomReco[i]->GetRMS());
	MomPR[i]->SetLineColor(mcolor[1]);		   
	MomReco[i]->SetLineColor(mcolor[2]);
//        MomPR[i]->GetYaxis()->SetRangeUser(0.,1.10*MaxYMom);
      
        MomLine[i]=new TLine(1./Ene[i],0.,1/Ene[i],1.10*MaxYMom);
        MomLine[i]->SetLineColor(kGreen);
        MomPR[i]->SetMaximum(1.10*MaxYMom);	
	MomReco[i]->Draw("sames");
//	MomLine[i]->Draw();
		 cout <<"creating tpavestats"<<endl; 
	gPad->Update();
	TPaveStats *tps = (TPaveStats*) MomPR[i]->FindObject("stats");
	tps->SetName("PR Stats");
	double Y1 = tps->GetY1NDC();
	double Y2 = tps->GetY2NDC();
	tps->SetTextColor(kRed);
	tps->SetLineColor(kRed);
        TPaveStats *tpreco = (TPaveStats*) MomReco[i]->FindObject("stats");
        tpreco->SetName("KF Stats");
	tpreco->SetY2NDC(Y1);
	tpreco->SetY1NDC(Y1-(Y2-Y1));
	tpreco->SetTextColor(kBlue);
	tpreco->SetLineColor(kBlue);
	gPad->Update();
	 
//	TPaveStats *tps2 = (TPaveStats*) MomReco[i]->FindObject("stats");	
	
	
	  cout<<"dRecoID it crash"<<endl;
	//leg3[i]->AddEntry(MomMC[i],Form("MC momentum %s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");
	
	for (int j=0;j<2;j++) {
		if (j==0) {
	float xMaxPR=MomPR[i]->GetXaxis()->GetBinCenter(MomPR[i]->GetMaximumBin());
	fMom[i][j] = new TF1(Form("fMomPR %d %d",i,j),"gaus",xMaxPR-fabs(0.2*xMaxPR),xMaxPR+fabs(0.2*xMaxPR));
	MomPR[i]->Fit(fMom[i][j],"R0");
	PRreco->SetPoint(i,1./Ene[i],fMom[i][j]->GetParameter(1));
	PRreco->SetPointError(i,0,fMom[i][j]->GetParameter(2));	
	PRreso->SetPoint(i,Ene[i],100*fabs((fMom[i][j]->GetParameter(2)/(fMom[i][j]->GetParameter(1)))));
		}
		if (j==1) {
	float xMaxReco=MomReco[i]->GetXaxis()->GetBinCenter(MomReco[i]->GetMaximumBin());
	fMom[i][j] = new TF1(Form("fMomReco %d %d",i,j),"gaus",xMaxReco-fabs(0.2*xMaxReco),xMaxReco+fabs(0.2*xMaxReco));
	MomReco[i]->Fit(fMom[i][j],"R0");
	KFreco->SetPoint(i,1./Ene[i],fMom[i][j]->GetParameter(1));
	KFreco->SetPointError(i,0,fMom[i][j]->GetParameter(2));	
	KFreso->SetPoint(i,Ene[i],100*fabs((fMom[i][j]->GetParameter(2)/(fMom[i][j]->GetParameter(1)))));
		}
	fMom[i][j]->SetLineColor(mcolor[j+1]);
	fMom[i][j]->Draw("same"); 
	if(t==3||t==4)PTMom[i]=new TPaveText(1.0,MaxYMom-(j+1)*0.3*MaxYMom,2.,MaxYMom-j*0.1*MaxYMom);
	if(t==11)PTMom[i]=new TPaveText(-0.019,MaxYMom-(j+1)*0.1*MaxYMom,-0.008,MaxYMom-j*0.1*MaxYMom);
	if(t==10)PTMom[i]=new TPaveText(0.008,MaxYMom-(j+1)*0.1*MaxYMom,0.019,MaxYMom-j*0.1*MaxYMom);

	 } //end j
	  cout <<"and now?"<<endl;
	leg3[i]->AddEntry(MomPR[i],Form("MomPR fit #mu = %5.4f %s^{-1}, #sigma_{reso} = %4.1f%% ",fMom[i][0]->GetParameter(1),UNIT.c_str(),(100*fabs(((fMom[i][0]->GetParameter(2))/(fMom[i][0]->GetParameter(1)))))),"l");
	leg3[i]->AddEntry(MomReco[i],Form("MomKF fit #mu = %5.4f %s^{-1}, #sigma_{reso} = %4.1f%% ",fMom[i][1]->GetParameter(1),UNIT.c_str(),(100*fabs(((fMom[i][1]->GetParameter(2))/(fMom[i][1]->GetParameter(1)))))),"l");
    leg3[i]->AddEntry(MomLine[i],"1/eMC ", "l");
	leg3[i]->Draw("same");
	tps->Draw("same");
	 /*
	 	cout << "0"<<endl;
	  gPad->Update();
	tps2->SetX1NDC(X1);
	cout << "1"<<endl;
	tps2->SetX2NDC(X2);
	 cout << "2"<<endl;
	tps2->SetY1NDC(Y1-(Y2-Y1));
	 cout << "3"<<endl;
	tps2->SetY2NDC(Y1);
	 cout << "4"<<endl;
	tps2->SetTextColor(kBlue);
	tps2->SetLineColor(kBlue);
    tps2->Draw("same");
*/
	 fileout->cd();
    MomMC[i]->Write();
    MomPR[i]->Write();
    MomReco[i]->Write();	
   }	   
 
 cout << " print third canvas" << endl;
// can3->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf", directory.c_str(),type, startfile.c_str(), type,RecoID.c_str()));
 can3->Write();
 
 
    //////////////////////////////////////////////////////////////////
TCanvas*canreso=new TCanvas("Resolution PR&KF","Resolution PR&KF",200,10,1500,1000);
canreso->Divide(2,1);
preco->Add(PRreco);
preco->Add(KFreco);
resoreco->Add(PRreso);
resoreco->Add(KFreso);
canreso->cd(1);
 TF1*yx=new TF1("yx","x",0,1./Ene[0]);	
 yx->SetLineColor(kGray);
 yx->SetLineStyle(6);
preco->Draw("AP");
yx->Draw("same");
canreso->cd(2);
resoreco->Draw("AP");
 preco->GetXaxis()->SetRangeUser(0,Ene[0]*1.10);
 preco->GetYaxis()->SetTitle(Form("1/p_{reco} (%s^{-1})",UNIT.c_str()));
 preco->GetXaxis()->SetTitle(Form("1/p_{MC} (%s^{-1})",UNIT.c_str()));
 resoreco->GetXaxis()->SetRangeUser(0,Ene[3]*1.10);
 resoreco->GetYaxis()->SetTitle("#sigma_{reco} (%)");
 resoreco->GetXaxis()->SetTitle(Form("p_{MC} (%s)",UNIT.c_str()));
 canreso->Update();
// canreso->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf)", directory.c_str(),type, startfile.c_str(), type,RecoID.c_str()));
 fileout->Close(); 	
	
 toTeX.close();

} //end function

