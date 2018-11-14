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

void Pulls(int t, int nene, int ncycles, string s) ;

void Pulls(int t, int nene, int ncycles, string s)
{
	
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
  double B = 0.3;		//average magnetic field in T
  double c = TMath::C();
 string MCparamfile="./MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh);
 	
 
 
 int NTConfig=1;//Number of configuration of triggers
 string TConfig[1]={"T1&T4&T3"};
 //double TrigThresh[3]={0.3, 0.29, 0.57};		//trigger thresholds in %s determined for T1/T3/T4 
 

	
	
//Input file 
string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V3/Disc";
 //string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V2_NotDisc";
 //string Outpath="/home/smechbal/Documents/AESOPLITE/Analysis/V2/";	
 string startfile="aesopliteNonUniB_V3";
 string endfile="_fort.99";
 string directory= "Disc/V3/NonUniB/toshow";
string RecoInd;
 if(t==3 || t==4)	RecoInd= "";
 if(t==11||t==10)  RecoInd="";
 string TitleID ="allPRevents";
 string ID = s;
	
//type of particle
 int type=t;
 string stype;	
 if(t==3)	stype="e^{#minus}";
 if(t==4)	stype="e^{#plus}";
 if(t==11)	stype="#mu^{#minus}";
 if(t==10)	stype="#mu^{#plus}";

 float mass=0.000511;//electon mass in GeV
 if(t==11 || t==10)	mass=0.10566;//muon mass in GeV
//Number of energies
 int Nene=nene;
//Energies
 int*Ene=new int[Nene];
if(t==3 || t==4)
{  
 
 Ene[0]= 20;
 Ene[1]= 40;
 Ene[2]= 100;
 Ene[3]= 300;
 
 //Ene[0]=40;
 //Ene[1]=100;
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

 
 ////Coordinates pull distributions 
 TH1F ****LayerPullPR=new TH1F***[Nene];
 TH1F ****LayerPullReco=new TH1F***[Nene];
	

	

 for(int i=0; i<Nene; i++) 
   {
	if((t==3 || t==4) && Ene[i]>=60)Nbins=2*Nbins;
   	if((t==3 || t==4) && Ene[i]>=100)Nbins=3*Nbins; 
   //Pull distrubutions
    LayerPullPR[i]=new TH1F**[7];
	LayerPullReco[i]=new TH1F**[7];
       for(int j=0;j<7;j++)
      {
       LayerPullPR[i][j]= new TH1F*[3];
       LayerPullReco[i][j]= new TH1F*[3];
		  for(int k=0;k<3;k++) 
		 {
		 LayerPullPR[i][j][k] = new TH1F(Form("Pulls from PR L%d,%d MeV, coord %d",j+1,Ene[i],k),Form("Pulls from PR L%d,%d MeV, coord %d",j+1,Ene[i], k),5000,-5,5); 
		 LayerPullReco[i][j][k] = new TH1F(Form("Pulls from Reco L%d,%d MeV, coord %d",j+1,Ene[i],k),Form("Pulls from reco L%d,%d MeV, coord %d",j+1,Ene[i], k),5000,-5,5); 
     		 } //end k
	  } //end j
   } //end i
	
//cout << "chaining events" << endl;	
	
//Create TChain to merge the files for a given energy
	
 TChain **chain=new TChain*[Nene]; 

//output ROOT file
 TFile*fileout=new TFile(Form("%s/%d/Pulls_%s_%d_%s.root", directory.c_str(),type,startfile.c_str(), type,ID.c_str()),"RECREATE");

	
for(int i=0; i<Nene; i++) 
  {  	
   cout << "Energy: " << Ene[i] << UNIT <<endl;
   chain[i]=new TChain("MC");
   for(int j=0;j<Ncycles[i];j++)//Number of cycles
     { 
if(t==3 || t==4)      chain[i]->Add(Form("%s/%d/RecoEvent_%s_%d_%d%s%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene[i],UNIT.c_str(),j+1,endfile.c_str(),ID.c_str()));
if(t==10 || t==11)        chain[i]->Add(Form("%s/%d/Reco_%s_%d_%d%s%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene[i],UNIT.c_str(),j+1,endfile.c_str(),ID.c_str()));  
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
		if(!T1 || !T4) continue;
		 
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

 gStyle->SetOptStat(1111);	


 /////////////////////////////////////////////////////////
 ////Pull distributions histograms//////////////////////
////////////////////////////////////////////////////////

  TCanvas ***CanXYZReso=new TCanvas**[Nene];
  TLegend***legPull=new TLegend**[Nene];
  TF1****fPullPR= new TF1***[Nene];
  TF1****fPullReco =new TF1***[Nene];
  TPaveText****PTPull= new TPaveText***[Nene];
  int nsig=3;
  for(int h=0; h<Nene;h++) { 
  CanXYZReso[h]=new TCanvas*[7];
  legPull[h]=new TLegend*[7];
  fPullPR[h] = new TF1**[7];
  fPullReco[h] = new TF1**[7];
  PTPull[h] = new TPaveText**[7];
    for(int i=0;i<7;i++) 
	 {     
     CanXYZReso[h][i]=new TCanvas(Form("XYZReso Layer %d, %d%s", i, Ene[h],UNIT.c_str()),Form("XYZReso Layer %d, %d%s",i,Ene[h],UNIT.c_str()),200,10,1200,800);    
     CanXYZReso[h][i]->Divide(3,1);
     legPull[h][i]= new TLegend(0.2,0.8,0.7,0.6);
     legPull[h][i]->SetBorderSize(0);
     legPull[h][i]->SetFillStyle(0);
     fPullPR[h][i] = new TF1*[2];
     fPullReco[h][i] = new TF1*[2];
     PTPull[h][i] = new TPaveText*[2];
	 for(int j=0;j<3;j++) 
	 {
      LayerPullPR[h][i][j]->SetLineColor(kRed);
      LayerPullReco[h][i][j]->SetLineColor(kBlue);
 
      CanXYZReso[h][i]->cd(j+1);
      if(j==0) {
	LayerPullPR[h][i][j]->SetTitle(Form("Pull X, %d %s, Tracker layer %d",Ene[h],UNIT.c_str(),i));
	LayerPullPR[h][i][j]->GetXaxis()->SetTitle("pull (mm)");
	float xMaxPR=LayerPullPR[h][i][j]->GetXaxis()->GetBinCenter(LayerPullPR[h][i][j]->GetMaximumBin());
	float xMaxReco=LayerPullReco[h][i][j]->GetXaxis()->GetBinCenter(LayerPullReco[h][i][j]->GetMaximumBin());
	fPullPR[h][i][j] = new TF1(Form("fPullPR%d %d %d",h,i,j),"gaus",xMaxPR-fabs(0.2*xMaxPR),xMaxPR+fabs(0.2*xMaxPR));
	fPullReco[h][i][j] = new TF1(Form("fPullReco%d %d %d",h,i,j),"gaus",xMaxReco-fabs(0.2*xMaxReco),xMaxReco+fabs(0.2*xMaxReco));
        LayerPullPR[h][i][j]->Fit(fPullPR[h][i][j],"R0");
	LayerPullReco[h][i][j]->Fit(fPullReco[h][i][j],"R0");
        fPullPR[h][i][j]->SetLineColor(kRed);
        fPullReco[h][i][j]->SetLineColor(kBlue);
	double meanReco = LayerPullReco[h][i][j]->GetMean();
	double sigmaReco = LayerPullReco[h][i][j]->GetStdDev();
	LayerPullPR[h][i][j]->GetXaxis()->SetRangeUser(meanReco -nsig*sigmaReco, meanReco + nsig*sigmaReco);
	}
		    if(j==1) {
	LayerPullPR[h][i][j]->SetTitle(Form("Pull Y , %d %s, Tracker layer %d",Ene[h],UNIT.c_str(),i));
	LayerPullPR[h][i][j]->GetXaxis()->SetTitle("pull(mm)");
	float xMaxPR=LayerPullPR[h][i][j]->GetXaxis()->GetBinCenter(LayerPullPR[h][i][j]->GetMaximumBin());
	float xMaxReco=LayerPullReco[h][i][j]->GetXaxis()->GetBinCenter(LayerPullReco[h][i][j]->GetMaximumBin());
	fPullPR[h][i][j] = new TF1(Form("fPullPR%d %d %d",h,i,j),"gaus",xMaxPR-fabs(0.2*xMaxPR),xMaxPR+fabs(0.2*xMaxPR));
	fPullReco[h][i][j] = new TF1(Form("fPullReco%d %d %d",h,i,j),"gaus",xMaxReco-fabs(0.2*xMaxReco),xMaxReco+fabs(0.2*xMaxReco));
        LayerPullPR[h][i][j]->Fit(fPullPR[h][i][j],"R0");
	LayerPullReco[h][i][j]->Fit(fPullReco[h][i][j],"R0");
	fPullPR[h][i][j]->SetLineColor(kRed);
        fPullReco[h][i][j]->SetLineColor(kBlue);
	double meanReco = LayerPullReco[h][i][j]->GetMean();
        double sigmaReco = LayerPullReco[h][i][j]->GetStdDev();
        LayerPullPR[h][i][j]->GetXaxis()->SetRangeUser(meanReco -nsig*sigmaReco, meanReco + nsig*sigmaReco);

   }
		    if(j==2) {
	LayerPullPR[h][i][j]->SetTitle(Form("Pull Z, %d %s, Tracker layer %d",Ene[h],UNIT.c_str(),i));
	LayerPullPR[h][i][j]->GetXaxis()->SetTitle("pull (mm)");
   }
     
   float tmp1=0;
   float max=0;
   tmp1=LayerPullPR[h][i][j]->GetBinContent(LayerPullPR[h][i][j]->GetMaximumBin());
   float tmp2=0;
   tmp2=LayerPullReco[h][i][j]->GetBinContent(LayerPullReco[h][i][j]->GetMaximumBin());
   if(tmp1>tmp2) max=tmp1;
   else max=tmp2;
   LayerPullPR[h][i][j]->SetMaximum(1.10*max);
   LayerPullPR[h][i][j]->Draw("hist");
   gPad->Update();
   LayerPullReco[h][i][j]->Draw("sames");
   gPad->Update();
   if(j!=2) {
   //fPullPR[h][i][j]->Draw("same");
   //fPullReco[h][i][j]->Draw("same");
   if(t==3||t==4)PTPull[h][i][j]=new TPaveText(0.6,0.7,0.9,0.9,"NDC");
   if(t==11)PTPull[h][i][j]=new TPaveText(-0.019,max-(j+1)*0.1*max,-0.008,max-j*0.1*max);
   if(t==10)PTPull[h][i][j]=new TPaveText(0.008,max-(j+1)*0.1*max,0.019,max-j*0.1*max);
    PTPull[h][i][j]->AddText(Form("PR, #mu=%5.6f mm, #sigma= %5.6f mm",fPullPR[h][i][j]->GetParameter(1),fPullPR[h][i][j]->GetParameter(2)));  
    PTPull[h][i][j]->AddText(Form("Reco, #mu=%5.6f mm, #sigma= %5.6f mm",fPullReco[h][i][j]->GetParameter(1),fPullReco[h][i][j]->GetParameter(2)));   
    PTPull[h][i][j]->SetFillStyle(0);PTPull[h][i][j]->SetBorderSize(0);
    PTPull[h][i][j]->Draw("same");
  cout<<"retrieving PR stat box"<<endl;
	gPad->Update();
         TPaveStats *sPR = (TPaveStats*)LayerPullPR[h][i][j]->FindObject("stats");
        cout<<"retrieved PR stat box"<<endl; 
        sPR->SetTextColor(kRed);
        double y1 = sPR->GetY1NDC();
        double y2 = sPR->GetY2NDC();
        cout<<"y1="<<y1<<"	, y2="<<y2<<endl;
	gPad->Update();
        TPaveStats *sKF = (TPaveStats*)LayerPullReco[h][i][j]->FindObject("stats");
        gPad->Update(); 
        sKF->SetTextColor(kBlue);
        sKF->SetY2NDC(y1);
        sKF->SetY1NDC(y1-fabs(y2-y1));
	gPad->Update();
 
  }
   fileout->cd();
   LayerPullPR[h][i][j]->Write();
   LayerPullReco[h][i][j]->Write();
 
	  if (j==0) {
  	  legPull[h][i]->AddEntry(LayerPullPR[h][i][j],"Pull from PR","l");
	  legPull[h][i]->AddEntry(LayerPullReco[h][i][j],"Pull from Reco","l");	  
	  legPull[h][i]->Draw("same");
	  }
	
       }//j   
    if(h==0 && i==0) CanXYZReso[h][i]->Print(Form("%s/%d/Pulls_%s_%d_%s.pdf(", directory.c_str(),type, startfile.c_str(), type,ID.c_str()));
    else if(h==(Nene-1) && i==6) CanXYZReso[h][i]->Print(Form("%s/%d/Pulls_%s_%d_%s.pdf)", directory.c_str(),type, startfile.c_str(), type,ID.c_str())); 
    else CanXYZReso[h][i]->Print(Form("%s/%d/Pulls_%s_%d_%s.pdf", directory.c_str(),type, startfile.c_str(), type,ID.c_str()));
    CanXYZReso[h][i]->Write();
   } //i
  } //h


} //end function
