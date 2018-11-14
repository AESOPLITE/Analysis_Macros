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

ClassImp(ALTckhit)
ClassImp(ALEvent)

void Scattering(int t, int nene) ;


void Scattering(int t, int nene)
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

//Number of reconstruction types to study 


string RecoInd="KFone";


//Input file 
string Inppath="/home/sarah/AESOPLITE/MCProduction/Detector/ProcessedFiles";	
string startfile="aesopliteNonUniB_V4";
string endfile="_fort.99";
string source = "30cmSource";
string directory= "/home/sarah/AESOPLITE/Analysis/MCAnalysis/MC_Synthesis/V4/3";

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
if(t==11 || t==10)     mass=0.10566;//muon mass in GeV
if(t==1)     mass=0.93827;//proton mass in GeV
if(t==6)     mass=3.7273;//alpha-particle mass in GeV

//float Gsource = 2035.01;		//geometry factor of 30cm beam source, thetaMax = 28.1 deg
float Gsource = 5476.16;		//geometry factor for 50 cm beam, thetaMax 28.1 deg
//Number of energies
int Nene=nene;
//Energies
int*Ene=new int[Nene];
//Number of cycles per energy
 int* Ncycles=new int[Nene];
 for(int i=0;i<Nene;i++)Ncycles[i]=10;

if(t==3 || t==4)
{ 
 
 Ene[0]= 20;
 //Ene[1]= 40;
// Ene[2]= 100;
 //Ene[1]= 300;
  }
 string UNIT="MeV";	
 int Nbins=500;
 float EneMin=0;
 float EneMax=500;
 float Dmin=-5.;
 float Dmax=5.;
  
//muons	
 if(t==1||t==10)
  {
   Ene[0]=1;
   Ene[1]=2;
   Ene[2]=3;
   //Ene[3]=15;
   //Ene[4]=20;
   UNIT="GeV";
   Nbins=500;
   EneMin=500;
   EneMax=5000;
   Dmin=-0.05;
   Dmax=0.05;
  }	
 int mcolor[4]={2,4,6,8}; 
//zz0
float zz0=0;	
//Number of entries per reconstruction type per energy 
	

//define histograms output

 double ymin = -20.0;				//y-bounds of T1
 double ymax = +20.0;				//y-bounds of T1
 double thetamin = -0.8;
 double thetamax = 0.8;

 //true deflection insRecoInde tracker
 TH1F **deflectionNB=new TH1F*[Nene];//L6-L0
 TH1F **deflectionB=new TH1F*[Nene];//L5-L1

 //Deflection Pattern Recognition
 TH1F **deflectionBPR=new TH1F*[Nene];//L5-L1
 TH1F **deflectionNBPR=new TH1F*[Nene];//L6-L0 
 
 ////Coordinates pull distributions 
 TH1F ****LayerPull=new TH1F***[Nene];;	
	
	
 //Loss of energy between injection point and L0
 TH1F **Elossfront=new TH1F*[Nene];

 //PR energy - E at L0
 TH1F **EResoPR=new TH1F*[Nene];	
 TH1F **Mom=new TH1F*[Nene];	
 TH1F **EPR=new TH1F*[Nene];	
	

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
	
	Elossfront[i] =new TH1F(Form("Eloss (truth), Injection Point -L0, %d%s",Ene[i],UNIT.c_str()),Form("Eloss (truth), Injection Point -L0, %d%s",Ene[i],UNIT.c_str()),100,0.,10.); 
    Elossfront[i]->GetXaxis()->SetTitle("#DeltaE=E_{IP}-E_{L0} (MeV)");   
   	Elossfront[i]->GetYaxis()->SetTitle("arbitrary unit");	   
    EResoPR[i] =new TH1F(Form("EL0 - EPR %d%s",Ene[i],UNIT.c_str()),Form("EL0 - EPR %d%s",Ene[i],UNIT.c_str()),2000,-1000.,1000.); 
    EResoPR[i]->GetXaxis()->SetTitle("#DeltaE=E_{L0}-E_{PR} (MeV)");   
   	EResoPR[i]->GetYaxis()->SetTitle("arbitrary unit");	
    Mom[i] =new TH1F(Form("Momentum, %d%s",Ene[i],UNIT.c_str()),Form("Momentum, %d%s",Ene[i],UNIT.c_str()),1000,EneMin,EneMax); 
    Mom[i]->GetXaxis()->SetTitle(Form("Momentum PR (%s)",UNIT.c_str()));   
   	Mom[i]->GetYaxis()->SetTitle("arbitrary unit");	
    EPR[i] =new TH1F(Form("Energy from PR, %d%s",Ene[i],UNIT.c_str()),Form("Energy from PR, %d%s",Ene[i],UNIT.c_str()),1000,EneMin,EneMax); 
    EPR[i]->GetXaxis()->SetTitle(Form("Energy from PR (%s)",UNIT.c_str()));   
    EPR[i]->GetYaxis()->SetTitle("arbitrary unit");	
    } // end i 
   //Pull distrubutions
   for(int i=0;i<Nene;i++) {
    LayerPull[i]=new TH1F**[7];
       for(int j=0;j<7;j++)
      {
       LayerPull[i][j]= new TH1F*[3];
		  for(int k=0;k<3;k++) 
		 {
		 LayerPull[i][j][k] = new TH1F(Form("Pulls L%d,%d MeV, coord %d",j+1,Ene[i],k),Form("Pulls L%d,%d MeV, coord %d",j+1,Ene[i], k),200,-2,2); 
     		 } //end k
	  } //end j
   } //end i
	
//cout << "chaining events" << endl;	
	
//Create TChain to merge the files for a given energy
	
 TChain **chain=new TChain*[Nene]; 

//output ROOT file
 //TFile*fileout=new //TFile(Form("./BinnedDeflectionHistogram_%s_%d_%s.root",startfile.c_str(),type,TitleRecoInd.c_str()),"RECREATE");

	
for(int i=0; i<Nene; i++) 
  {  	
   cout << "Energy: " << Ene[i] << UNIT <<endl;
   chain[i]=new TChain("MC");
   for(int j=0;j<Ncycles[i];j++)//Number of cycles
     { 
        chain[i]->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene[i],UNIT.c_str(),endfile.c_str(),RecoInd.c_str()));
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
      if(j%1000000==0) cout << "Event: " << j <<endl;
	  uint8_t Ti=(uint8_t)e->get_Ti();
	  //Number of layers with hit(s)
        int NL= e->get_NLayers();
	 
	  //Internal trigger requirement: at least 5 hits, event chosen by PR
	  if((e->get_deflecPR()!=0))
	   {
		//apply external trigger requirements
		bool T1=e->get_T1();
		bool T3=e->get_T3();
		bool T4=e->get_T4();
		if(T1 && T4)
		 {
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
		  deflectionBPR[i]->Fill(deflecPRTheta);
		  //Signed curvature at the 4 points of the bending plane
	  float zzz[4]={TckZPos[5],TckZPos[3],TckZPos[2],TckZPos[1]};
          float curv[4]={0,0,0,0};
          TF1* fcurv=new TF1("fcurv","2*[0]/TMath::Power(1+TMath::Power(2*[0]*x+[1],2),3./2.)",-20,20);
          fcurv->SetParameter(0,a);
          fcurv->SetParameter(1,2*a*zz0+b);
		  float CurvMean=0; 	 
          for(int ij=0;ij<4;ij++)
            {
             curv[ij]= fcurv->Eval(zzz[ij]);   
			 CurvMean+=	curv[ij]/4;
			}   
	      float Rmean=1./curv[2];
		  if(CurvMean!=0)	 Rmean=1./CurvMean;
	 
		  // non-bending plane PR variables
		  double slopeNB = e->get_slopePR();
		  double thetaNB = TMath::ATan(slopeNB);
			 
		  //Extract a simple estimation of energy from the average of the 4 curvature a radius	 
		  float B=0.3; //in T	
		  float Pt=0.3 * B *	0.01*Rmean; //in GeV
		  float P= Pt / TMath::Cos(fabs(thetaNB));   //in GeV
		 	 
		  Mom[i]->Fill(1000.*P);	 
			 	 
		  float E=TMath::Sqrt(P*P+mass*mass);
	          EPR[i]->Fill(1000.*E);		 
			 
		  double slopeMCL0=0; 
		  double thetaMCL0=0; 
		  double slopePRL0=0;
		  double thetaPRL0=0;
		  double slopeMCL1=0; 
		  double thetaMCL1=0; 
		  double slopeMCL5=0; 
		  double thetaMCL5=0; 
		  double slopeMCL6=0; 
		  double thetaMCL6=0; 
		  double slopePRL6=0;
		  double thetaPRL6=0;
		  for(int k=0;k<e->get_Nhits();k++) 
		    {
			 int Lindex=e->get_hits().at(k)->get_L();
			 if((int)e->get_hits().at(k)->get_flagPR()!=1) continue;
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
			   EResoPR[i]->Fill(1000.*(e->get_hits().at(k)->get_eMC() - E));
			  //PR variables at L0
			   double cxL0PR = e->get_hits().at(k)->get_cxPR();
			   double cyL0PR = e->get_hits().at(k)->get_cyPR();
			   double czL0PR = e->get_hits().at(k)->get_czPR();
			   slopePRL0= (cxL0PR/czL0PR);
			   thetaPRL0 = TMath::ATan(slopePRL0); 		 
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
		      }	
		     //LAYER L5, B PLANE
		     if(Lindex==5)
		      {						
			   //truth variables at L5
			   double cyout = e->get_hits().at(k)->get_cy();
			   double czout = e->get_hits().at(k)->get_cz();
			   slopeMCL5= (cyout/czout);
			   thetaMCL5 = TMath::ATan(slopeMCL5);
			   double deflectionMC = slopeMCL5 - slopeT1B;
               double deflecMCTheta = thetaMCL5 - thetaT1B;
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
			  } //end if
		    }	//end for

		   //Fill histograms	 
		   deflectionNB[i]->Fill(thetaMCL6-thetaMCL0); 	 
		   deflectionB[i]->Fill(thetaMCL5-thetaMCL1);
		   deflectionNBPR[i]->Fill(thetaPRL6-thetaPRL0); 
		
		//Fill pull distributions histograms
		for(int k=0;k<e->get_Nhits();k++)
         {
          if(e->get_hits().at(k)->get_flagPR())       //Check that the hit was used for PR
           {
	    float tmp=0;
	    int Lindex=(int)e->get_hits().at(k)->get_L();
	    tmp=(e->get_hits().at(k)->get_x());
	    LayerPull[i][Lindex][0]->Fill(tmp-e->get_hits().at(k)->get_xPR());
	    tmp=(e->get_hits().at(k)->get_y());
	    LayerPull[i][Lindex][1]->Fill(tmp-e->get_hits().at(k)->get_yPR());
	    tmp=(e->get_hits().at(k)->get_z());
	    LayerPull[i][Lindex][2]->Fill(tmp-e->get_hits().at(k)->get_zPR());	
	    } //end if flagPR
	  } //end k
		  } //end if thresholds
	    } //end internal trigger condition 
     }  //end j
  } //end i

	   
//////////////////////////////   
// Display  
//////////////////////////////   
 
 gStyle->SetOptStat(0);	

	
 TCanvas*can=new TCanvas("MC deflection truth","MC deflection truth",200,10,1500,1000);
 TLegend*leg=new TLegend(0.2,0.5,0.3,0.9);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);
 for(int i=0;i<Nene;i++) leg->AddEntry(deflectionNB[i],Form("%s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");

 can->Divide(1,2);
 //Non-Bending plane
 can->cd(1);	
 float MaxYNB=0;
 TF1**fDeflecNB=new TF1*[Nene];
	
 for(int i=0;i<Nene;i++) 
   {	
 	fDeflecNB[i]=new TF1(Form("fDeflecNB%d",i),"gaus",-2,1);
	deflectionNB[i]->Fit(fDeflecNB[i],"0","",-0.4,0.4);
	fDeflecNB[i]->SetLineColor(mcolor[i]);
   }	
	
 for(int i=0;i<Nene;i++)
   {
    float tmp=deflectionNB[i]->GetBinContent(deflectionNB[i]->GetMaximumBin());
	if(tmp>MaxYNB)  MaxYNB=tmp; 
   } 
 deflectionNB[0]->SetLineColor(mcolor[0]);
 deflectionNB[0]->SetMaximum(1.05*MaxYNB);
 if(t==11||t==10)deflectionNB[0]->GetXaxis()->SetRangeUser(-0.01,0.01);

 deflectionNB[0]->SetTitle("MC Scattering (truth) in non bending plane, L6-L0");
 deflectionNB[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {
	deflectionNB[i]->SetLineColor(mcolor[i]);
	deflectionNB[i]->Draw("same");   
   }	   
	
 	
 leg->Draw("same");
 TPaveText**PTNB=new TPaveText*[Nene];
	
 for(int i=0;i<Nene;i++) 
   {	
	if(t==3 || t==4)PTNB[i]=new TPaveText(0.4,MaxYNB-(i+1)*0.1*MaxYNB,0.9,MaxYNB-i*0.1*MaxYNB);
	if(t==11 || t==10)PTNB[i]=new TPaveText(0.003,MaxYNB-(i+1)*0.1*MaxYNB,0.01,MaxYNB-i*0.1*MaxYNB);
    PTNB[i]->AddText(Form("%s, %d%s, #mu=%3.0f mrad, #sigma= %3.0f mrad",stype.c_str(),Ene[i],UNIT.c_str(),1000.*fDeflecNB[i]->GetParameter(1),1000.*fDeflecNB[i]->GetParameter(2)));
    PTNB[i]->SetFillStyle(0);PTNB[i]->SetBorderSize(0);
	PTNB[i]->Draw();
	fDeflecNB[i]->Draw("same");  
   }	
		
 //Bending plane
 can->cd(2);
 float MaxYB=0;
  TF1**fDeflecB=new TF1*[Nene];
	
 for(int i=0;i<Nene;i++) 
   {	
 	float xMax=deflectionB[i]->GetXaxis()->GetBinCenter(deflectionB[i]->GetMaximumBin());
	fDeflecB[i]=new TF1(Form("fDeflecB%d",i),"gaus",xMax-fabs(0.2*xMax),xMax+fabs(0.2*xMax));
	deflectionB[i]->Fit(fDeflecB[i],"R0");
	fDeflecB[i]->SetLineColor(mcolor[i]);
   }   
	
 for(int i=0;i<Nene;i++)
   {
    float tmp=deflectionB[i]->GetBinContent(deflectionB[i]->GetMaximumBin());
	if(tmp>MaxYB)  MaxYB=tmp; 
   } 
	
 deflectionB[0]->SetLineColor(mcolor[0]);
 deflectionB[0]->SetMaximum(1.05*MaxYB);
 deflectionB[0]->GetXaxis()->SetRangeUser(-1,0.1);
 if(t==11)deflectionB[0]->GetXaxis()->SetRangeUser(-0.02,0.005);
 if(t==10)deflectionB[0]->GetXaxis()->SetRangeUser(-0.005,0.02);
 deflectionB[0]->SetTitle("MC Deflection and Scattering (truth) in bending plane, L5-L1");
 deflectionB[0]->Draw("hist");
 cout << " Or here? " << endl;	
 for(int i=1;i<Nene;i++) 
   {	 
	deflectionB[i]->SetLineColor(mcolor[i]);
	deflectionB[i]->Draw("sames");   
   }		
 TPaveText**PTB=new TPaveText*[Nene];
	
 for(int i=0;i<Nene;i++) 
   {	
	if(t==3||t==4)PTB[i]=new TPaveText(-0.9,MaxYB-(i+1)*0.1*MaxYB,-0.5,MaxYB-i*0.1*MaxYB);
	if(t==11)PTB[i]=new TPaveText(-0.019,MaxYB-(i+1)*0.1*MaxYB,-0.008,MaxYB-i*0.1*MaxYB);
	if(t==10)PTB[i]=new TPaveText(0.008,MaxYB-(i+1)*0.1*MaxYB,0.019,MaxYB-i*0.1*MaxYB);
    PTB[i]->AddText(Form("%s, %d%s, #mu=%3.0f mrad, #sigma=%3.0f mrad, #sigma= %4.2f %%",stype.c_str(),Ene[i],UNIT.c_str(),1000.*fDeflecB[i]->GetParameter(1),1000.*fDeflecB[i]->GetParameter(2),100*fDeflecB[i]->GetParameter(2)/abs(fDeflecB[i]->GetParameter(1))));
    PTB[i]->SetFillStyle(0);PTB[i]->SetBorderSize(0);
	PTB[i]->Draw();
	fDeflecB[i]->Draw("same");  
   }

  can->Print(Form("PR_Analysis_%s_%d_%s.pdf(", startfile.c_str(), type,RecoInd.c_str()));
cout << " print first canvas " << endl;
								
 TCanvas*can2=new TCanvas("Energy loss truth","Energy loss truth",200,10,1500,1000);
 TLegend*leg2=new TLegend(0.8,0.5,0.9,0.9);
 leg2->SetBorderSize(0);
 leg2->SetFillStyle(0);
 for(int i=0;i<Nene;i++) leg2->AddEntry(Elossfront[i],Form("%s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");

 can2->cd(1);	
 float MaxYEl=0;
 for(int i=0;i<Nene;i++)
   {
    float tmp=Elossfront[i]->GetBinContent(Elossfront[i]->GetMaximumBin());
	if(tmp>MaxYEl)  MaxYEl=tmp; 
   } 
 Elossfront[0]->SetLineColor(mcolor[0]);
 Elossfront[0]->SetMaximum(1.05*MaxYEl);
 Elossfront[0]->SetTitle("Energy loss between injection point and  L0 (No CK gas)");
 Elossfront[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {	 
	Elossfront[i]->SetLineColor(mcolor[i]);
	Elossfront[i]->Draw("same");   
   }	   
 leg2->Draw("same");	  
  
  cout << " print second canvas" << endl;
  can2->Print(Form("PR_Analysis_%s_%d_%s.pdf", startfile.c_str(), type,RecoInd.c_str()));

 TCanvas*can3=new TCanvas("Momentum PR","Momentum PR",200,10,1500,1000);
 TLegend*leg3=new TLegend(0.8,0.5,0.9,0.9);
 leg3->SetBorderSize(0);
 leg3->SetFillStyle(0);
 for(int i=0;i<Nene;i++) leg3->AddEntry(Mom[i],Form("%s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");

 can3->cd(1);	
 float MaxYMom=0;
 for(int i=0;i<Nene;i++)
   {
    float tmp=Mom[i]->GetBinContent(Mom[i]->GetMaximumBin());
	if(tmp>MaxYMom)  MaxYMom=tmp; 
   } 
 Mom[0]->SetLineColor(mcolor[0]);
 Mom[0]->SetMaximum(1.05*MaxYMom);
 Mom[0]->SetTitle("Momentum from Pattern Recognition");
 Mom[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {	 
	Mom[i]->SetLineColor(mcolor[i]);
	Mom[i]->Draw("same");   
   }	   
 leg3->Draw("same");	
 
 cout << " print third canvas" << endl;
 can3->Print(Form("PR_Analysis_%s_%d_%s.pdf", startfile.c_str(), type,RecoInd.c_str()));
	
 TCanvas*can4=new TCanvas("Energy PR","Energy PR",200,10,1500,1000);
 TLegend*leg4=new TLegend(0.8,0.5,0.9,0.9);
 leg4->SetBorderSize(0);
 leg4->SetFillStyle(0);
 for(int i=0;i<Nene;i++) leg4->AddEntry(EPR[i],Form("%s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");

 can4->cd(1);	
 float MaxYEPR=0;
 for(int i=0;i<Nene;i++)
   {
    float tmp=EPR[i]->GetBinContent(EPR[i]->GetMaximumBin());
	if(tmp>MaxYEPR)  MaxYEPR=tmp; 
   } 
 EPR[0]->SetLineColor(mcolor[0]);
 EPR[0]->SetMaximum(1.05*MaxYEPR);
 EPR[0]->SetTitle("Energy from Pattern Recognition");
 EPR[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {	 
	EPR[i]->SetLineColor(mcolor[i]);
	EPR[i]->Draw("same");   
   }	   
 leg4->Draw("same");	
 cout << "print fourth canvas" << endl;
 can4->Print(Form("PR_Analysis_%s_%d_%s.pdf", startfile.c_str(), type, RecoInd.c_str()));

		
 TCanvas*can5=new TCanvas("Resolution Energy PR","Reolution Energy PR",200,10,1500,1000);
 TLegend*leg5=new TLegend(0.8,0.5,0.9,0.9);
 leg5->SetBorderSize(0);
 leg5->SetFillStyle(0);
 for(int i=0;i<Nene;i++) leg5->AddEntry(EResoPR[i],Form("%s, %d%s",stype.c_str(),Ene[i],UNIT.c_str()),"l");

 can5->cd(1);	
 float MaxYResoPR=0;
 for(int i=0;i<Nene;i++)
   {
    float tmp=EResoPR[i]->GetBinContent(EResoPR[i]->GetMaximumBin());
	if(tmp>MaxYResoPR)  MaxYResoPR=tmp; 
   } 
 EResoPR[0]->SetLineColor(mcolor[0]);
 EResoPR[0]->SetMaximum(1.05*MaxYResoPR);
 EResoPR[0]->SetTitle("Energy Resolution from Pattern Recognition");
 EResoPR[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {	 
	EResoPR[i]->SetLineColor(mcolor[i]);
	EResoPR[i]->Draw("same");   
   }	   
 leg5->Draw("same");	

 can5->Print(Form("PR_Analysis_%s_%d_%s.pdf", startfile.c_str(), type,RecoInd.c_str()));
								
 TCanvas*can6=new TCanvas("MC Deflection and Scattering (PR) in bending plane, L5-L1","MC Deflection and Scattering (PR) in bending plane, L5-L1",200,10,1500,1000);
 TLegend*leg6=new TLegend(0.2,0.7,0.6,0.9);
 leg6->SetBorderSize(0);
 leg6->SetFillStyle(0);
 
   TLegend*leg7=new TLegend(0.2,0.7,0.6,0.9);
 leg7->SetBorderSize(0);
 leg7->SetFillStyle(0);
 can6->cd();	
 float MaxYBPR=0;
 TF1**fDeflecBPR=new TF1*[Nene];
	
 for(int i=0;i<Nene;i++) 
   {	
 	float xMax=deflectionBPR[i]->GetXaxis()->GetBinCenter(deflectionBPR[i]->GetMaximumBin());
	fDeflecBPR[i]=new TF1(Form("fDeflecBPR%d",i),"gaus",xMax-fabs(0.2*xMax),xMax+fabs(0.2*xMax));
	deflectionBPR[i]->Fit(fDeflecBPR[i],"0");
	fDeflecBPR[i]->SetLineColor(mcolor[i]);
   }   

	
 for(int i=0;i<Nene;i++)
   {
    float tmp=deflectionBPR[i]->GetBinContent(deflectionBPR[i]->GetMaximumBin());
	if(tmp>MaxYBPR)  MaxYBPR=tmp; 
   } 
 deflectionBPR[0]->SetLineColor(mcolor[0]);
 deflectionBPR[0]->SetMaximum(1.05*MaxYBPR);
 deflectionBPR[0]->GetXaxis()->SetRangeUser(-1,0.1);
 if(t==11)deflectionBPR[0]->GetXaxis()->SetRangeUser(-0.03,0.005);
  if(t==10)deflectionBPR[0]->GetXaxis()->SetRangeUser(-0.005,0.03);
 deflectionBPR[0]->SetTitle("MC Deflection and Scattering (PR) in bending plane, L5-L1");

 TPaveText**PTBPR=new TPaveText*[Nene];
	
 deflectionBPR[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {	 
	deflectionBPR[i]->SetLineColor(mcolor[i]);
	deflectionBPR[i]->Draw("same");   
   }	   
		
	
 for(int i=0;i<Nene;i++) 
   {	
	PTBPR[i]=new TPaveText(0.2,MaxYBPR-(i+1)*0.1*MaxYBPR,1,MaxYBPR-i*0.1*MaxYBPR);
    PTBPR[i]->AddText(Form("e^{#minus}, %d, #mu=%3.0f mrad, #sigma=%3.0f mrad, #sigma= %4.2f %%",Ene[i],1000.*fDeflecBPR[i]->GetParameter(1),1000.*fDeflecBPR[i]->GetParameter(2),100*fDeflecBPR[i]->GetParameter(2)/abs(fDeflecB[i]->GetParameter(1))));
    PTBPR[i]->SetFillStyle(0);PTBPR[i]->SetBorderSize(0);
	PTBPR[i]->Draw();
	 leg6->AddEntry(deflectionBPR[i],Form("%s, %d%s, #mu=%3.0f mrad, #sigma=%3.0f mrad, #sigma= %4.2f %%",stype.c_str(),Ene[i],UNIT.c_str(),1000.*fDeflecBPR[i]->GetParameter(1),1000.*fDeflecBPR[i]->GetParameter(2),100*fDeflecBPR[i]->GetParameter(2)/abs(fDeflecBPR[i]->GetParameter(1))),"l");  
	fDeflecBPR[i]->Draw("same");  
   }
	
 leg6->Draw("same");
 //deflectionBPR[0]->GetXaxis()->SetLimits(-1,0.1);

 gPad->Update();		
//TGRAPH
	
 TGraphErrors*DefPRvsTruth=new TGraphErrors();
	
 for(int i=0;i<Nene;i++) 
   {
	DefPRvsTruth->SetPoint(DefPRvsTruth->GetN(),1000*fDeflecB[i]->GetParameter(1),1000*fDeflecBPR[i]->GetParameter(1)); 
	DefPRvsTruth->SetPointError(DefPRvsTruth->GetN()-1,1000*fDeflecB[i]->GetParameter(2),1000*fDeflecBPR[i]->GetParameter(2)); 
   }
/*
can6->cd(2);
float MaxYNBPR=0;
 TF1**fDeflecNBPR=new TF1*[Nene];
	
 for(int i=0;i<Nene;i++) 
   {	
 	float xMax=deflectionNBPR[i]->GetXaxis()->GetBinCenter(deflectionBPR[i]->GetMaximumBin());
	fDeflecNBPR[i]=new TF1(Form("fDeflecNBPR%d",i),"gaus",xMax-fabs(0.2*xMax),xMax+fabs(0.2*xMax));
	deflectionNBPR[i]->Fit(fDeflecNBPR[i],"0");
	fDeflecNBPR[i]->SetLineColor(mcolor[i]);
   }   

	
 for(int i=0;i<Nene;i++)
   {
    float tmp=deflectionNBPR[i]->GetBinContent(deflectionNBPR[i]->GetMaximumBin());
	if(tmp>MaxYNBPR)  MaxYNBPR=tmp; 
   } 
 deflectionNBPR[0]->SetLineColor(mcolor[0]);
 deflectionNBPR[0]->SetMaximum(1.05*MaxYBPR);
 deflectionNBPR[0]->GetXaxis()->SetRangeUser(-1,0.1);
 if(t==11||t==10)deflectionNBPR[0]->GetXaxis()->SetRangeUser(-0.03,0.005);
 deflectionNBPR[0]->SetTitle("MC Deflection and Scattering (PR) in non-bending plane, L6-L0");

 TPaveText**PTNBPR=new TPaveText*[Nene];
	
 deflectionNBPR[0]->Draw("hist");	
 for(int i=1;i<Nene;i++) 
   {	 
	deflectionNBPR[i]->SetLineColor(mcolor[i]);
	deflectionNBPR[i]->Draw("same");   
   }	   
		
	
 for(int i=0;i<Nene;i++) 
   {	
	PTNBPR[i]=new TPaveText(0.2,MaxYNBPR-(i+1)*0.1*MaxYNBPR,1,MaxYNBPR-i*0.1*MaxYNBPR);
    PTNBPR[i]->AddText(Form("e^{#minus}, %d, #mu=%3.0f mrad, #sigma=%3.0f mrad, #sigma= %4.2f %%",Ene[i],1000.*fDeflecNBPR[i]->GetParameter(1),1000.*fDeflecNBPR[i]->GetParameter(2),100*fDeflecNBPR[i]->GetParameter(2)/abs(fDeflecNB[i]->GetParameter(1))));
    PTNBPR[i]->SetFillStyle(0);PTNBPR[i]->SetBorderSize(0);
	PTNBPR[i]->Draw();
	 leg7->AddEntry(deflectionNBPR[i],Form("%s, %d%s, #mu=%3.0f mrad, #sigma=%3.0f mrad, #sigma= %4.2f %%",stype.c_str(),Ene[i],UNIT.c_str(),1000.*fDeflecNBPR[i]->GetParameter(1),1000.*fDeflecNBPR[i]->GetParameter(2),100*fDeflecBPR[i]->GetParameter(2)/abs(fDeflecNBPR[i]->GetParameter(1))),"l");  
	fDeflecNBPR[i]->Draw("same");  
   }
	
 leg7->Draw("same");
 //deflectionBPR[0]->GetXaxis()->SetLimits(-1,0.1);

 gPad->Update();		
*/
 can6->Print(Form("PR_Analysis_%s_%d_%s.pdf", startfile.c_str(), type,RecoInd.c_str()));

   //////////////////////////////////////////////////////////////////

/*   
   //Pull diswtributions histograms
  TCanvas ***CanXYZReso=new TCanvas**[Nene];
  for(int h=0; h<Nene;h++) { 
  CanXYZReso[h]=new TCanvas*[7];
    for(int i=0;i<7;i++) 
	 {

     CanXYZReso[h][i]=new TCanvas(Form("XYZReso Layer %d, %dMeV", i, Ene[h]),Form("XYZReso Layer %d, %dMeV",i,Ene[h]),200,10,1200,800);     
     CanXYZReso[h][i]->DivRecoInde(3,1);
	 for(int j=0;j<3;j++) 
	 {
      CanXYZReso[h][i]->cd(j+1);
      if(j==0) {
        LayerPull[h][i][j]->SetLineColor(kRed);
	LayerPull[h][i][j]->SetTitle(Form("Pull X (line fit), %d MeV, Tracker layer %d",Ene[h],i+1));
	LayerPull[h][i][j]->GetXaxis()->SetTitle("truth-reco (cm)");
   }
		    if(j==1) {
        LayerPull[h][i][j]->SetLineColor(kBlue);
	LayerPull[h][i][j]->SetTitle(Form("Pull Y (parabola fit), %d MeV, Tracker layer %d",Ene[h],i+1));
	LayerPull[h][i][j]->GetXaxis()->SetTitle("truth-reco (cm)");
   }
		    if(j==2) {
        LayerPull[h][i][j]->SetLineColor(kBlack);
	LayerPull[h][i][j]->SetTitle(Form("Pull Z, %d MeV, Tracker layer %d",Ene[h],i+1));
	LayerPull[h][i][j]->GetXaxis()->SetTitle("truth-reco (cm)");
   }
        
	float tmp2=0;
	tmp2=LayerPull[h][i][j]->GetBinContent(LayerPull[h][i][j]->GetMaximumBin());
	LayerPull[h][i][j]->SetMaximum(1.10*tmp2);
        LayerPull[h][i][j]->Draw("hist");
	LayerPull[h][i][j]->GetXaxis()->SetRangeUser(LayerPull[h][i][j]->GetMean()-3*LayerPull[h][i][j]->GetStdDev(),               LayerPull[h][i][j]->GetMean()+3*LayerPull[h][i][j]->GetStdDev());
	
	//if(j==0)leg->Draw("same");
       }//j   
 
    CanXYZReso[h][i]->Print(Form("PR_Analysis_%s_%d.pdf", startfile.c_str(), type));
    
   } //i
  } //h
*/	
	
 TCanvas*can7=new TCanvas("MC Deflection in bending plane, L5-L1","MC Deflection in bending plane, L5-L1",200,10,1200,1200);

 can7->cd(1);	
 gPad->SetLeftMargin(0.15);
 DefPRvsTruth->GetYaxis()->SetTitleOffset(1.7);
 DefPRvsTruth->GetXaxis()->SetTitle("True deflection (mrad)");
 DefPRvsTruth->GetYaxis()->SetTitle("Pattern recognition deflection (mrad)");
 DefPRvsTruth->SetTitle("MC deflection in bending plane, L5-L1");
 DefPRvsTruth->SetMarkerStyle(20);
 DefPRvsTruth->SetMarkerSize(1);
 if(t==3)
  {	 
   DefPRvsTruth->GetXaxis()->SetLimits(-900,0);
   DefPRvsTruth->SetMinimum(-900);
   DefPRvsTruth->SetMaximum(0);
  }
 if(t==11)
  {	 
   DefPRvsTruth->GetXaxis()->SetLimits(-20,0);
   DefPRvsTruth->SetMinimum(-20);
   DefPRvsTruth->SetMaximum(0);
  } 
   if(t==10)
  {	 
   DefPRvsTruth->GetXaxis()->SetLimits(0,20);
   DefPRvsTruth->SetMinimum(0);
   DefPRvsTruth->SetMaximum(20);
  } 
 DefPRvsTruth->Draw("ape");
	
 TF1*yx=new TF1("yx","x",-750,0);	
 yx->SetLineColor(kGray);
 yx->SetLineStyle(6);
 yx->Draw("same");
	
	
 can7->Print(Form("PR_Analysis_%s_%d_%s.pdf)", startfile.c_str(), type,RecoInd.c_str()));
  cout << "last canvas printed" << endl;
	
} //end function
