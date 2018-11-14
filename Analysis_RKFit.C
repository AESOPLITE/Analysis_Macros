////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Sarah Mechbal, smechbal@ucsc.edu
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

void RKFit(int t, int nene, int ncycles, string s) ;


void RKFit(int t, int nene, int ncycles, string s)
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


 //MC inverse energy/momentum at L0
 TH1F **MomReco=new TH1F*[Nene];
 TH1F **MomMC=new TH1F*[Nene];
 TH1F **MomPR=new TH1F*[Nene];	
 TH1F **EPR=new TH1F*[Nene];	
 


 for(int i=0; i<Nene; i++) 
   {
	if((t==3 || t==4) && Ene[i]>=60)Nbins=2*Nbins;
   	if((t==3 || t==4) && Ene[i]>=100)Nbins=3*Nbins;

    MomPR[i] =new TH1F(Form(" 1/p, %d%s",Ene[i],UNIT.c_str()),Form("1/p, %d%s",Ene[i],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
    MomPR[i]->GetXaxis()->SetTitle(Form("Inverse Momentum PR & Reco (%s^{-1})",UNIT.c_str()));   
   	MomPR[i]->GetYaxis()->SetTitle("#entries");	
    MomReco[i] =new TH1F(Form("Momentum from Reco, %d%s",Ene[i],UNIT.c_str()),Form("Momentum from Reco, %d%s",Ene[i],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
    MomMC[i] =new TH1F(Form("Momentum from MC, %d%s",Ene[i],UNIT.c_str()),Form("Momentum from MC, %d%s",Ene[i],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
 }
	
//cout << "chaining events" << endl;	
	
//Create TChain to merge the files for a given energy
	
 TChain **chain=new TChain*[Nene]; 

//output ROOT file
 TFile*fileout=new TFile(Form("%s/%d/Analysis_%s_%d_%s.root", directory.c_str(),type,startfile.c_str(), type,ID.c_str()),"RECREATE");

	
for(int i=0; i<Nene; i++) 
  {  	
   cout << "Energy: " << Ene[i] << UNIT <<endl;
   chain[i]=new TChain("MC");
   for(int j=0;j<Ncycles[i];j++)//Number of cycles
     { 
      chain[i]->Add(Form("%s/%d/RecoEvent_%s_%d_%d%s%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene[i],UNIT.c_str(),j+1,endfile.c_str(),ID.c_str()));
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

	  if((e->get_deflecPR()==0)) continue;
		//apply external trigger requirements
	bool T1=e->get_T1();
	bool T3=e->get_T3();
    bool T4=e->get_T4();
	 //Number of layers with hit(s) in bending/non-bending plane
	int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	if(NL!=7) continue;
	if(!T1 || !T4 || !T3) continue;
		 
		  //"truth" directional informations at entrance point T1 
		  double cxin = e->get_CX0MC();
		  double cyin = e->get_CY0MC();
		  double czin = e->get_CZ0MC();
		  double X0 = e->get_X0MC();
		  double Y0 = e->get_Y0MC();
		  double Z0 = e->get_Z0MC();

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
		  if(Ekreco==-999) continue;	 
		  MomReco[i]->Fill(1/pReco);

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
			   double eMC = 1000.*(e->get_hits().at(k)->get_eMC());		//energy at L0 in MeV
			   double pMC = TMath::Sqrt((eMC*eMC) - (mass*mass));
			   MomMC[i]->Fill(1/pMC);
			 }
		}
	  } //end k

     }  //end j
  } //end i

	   
//////////////////////////////   
// Display  
//////////////////////////////   
 
gStyle->SetOptStat(1);	


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
 preco->SetTitle("Inverse momentum recontruction PR (red) and RK (blue) vs truth momentum");
 resoreco->SetTitle("Momentum recontruction resolution PR (red) and RK (blue) vs truth momentum");


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
	else MaxYMom=tmpReco;
	MomReco[i]->SetTitle(Form("Inverse momentum PR & RK, %s %d%s",stype.c_str(),Ene[i],UNIT.c_str()));  
	MomReco[i]->Draw("hist");
	MomReco[i]->GetXaxis()->SetRangeUser(MomReco[i]->GetMean()-sigma*MomReco[i]->GetRMS(),MomReco[i]->GetMean()+sigma*MomReco[i]->GetRMS());
	MomReco[i]->SetLineColor(mcolor[1]);		   
	MomReco[i]->SetLineColor(mcolor[2]);
    MomReco[i]->SetMaximum(1.10*MaxYMom);	
	MomPR[i]->Draw("sames");
	gPad->Update();
	TPaveStats *tps = (TPaveStats*) MomPR[i]->FindObject("stats");
	tps->SetName("PR Stats");
	double Y1 = tps->GetY1NDC();
	double Y2 = tps->GetY2NDC();
	tps->SetTextColor(kRed);
	tps->SetLineColor(kRed);
        TPaveStats *tpreco = (TPaveStats*) MomReco[i]->FindObject("stats");
        tpreco->SetName("RK Stats");
	tpreco->SetY2NDC(Y1);
	tpreco->SetY1NDC(Y1-(Y2-Y1));
	tpreco->SetTextColor(kBlue);
	tpreco->SetLineColor(kBlue);
	gPad->Update();
	TPaveStats *tps2 = (TPaveStats*) MomReco[i]->FindObject("stats");	
	
	
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
	leg3[i]->AddEntry(MomPR[i],Form("MomPR fit #mu = %5.4f %s^{-1}, #sigma_{reso} = %4.1f%% ",fMom[i][0]->GetParameter(1),UNIT.c_str(),(100*fabs(((fMom[i][0]->GetParameter(2))/(fMom[i][0]->GetParameter(1)))))),"l");
	leg3[i]->AddEntry(MomReco[i],Form("MomKF fit #mu = %5.4f %s^{-1}, #sigma_{reso} = %4.1f%% ",fMom[i][1]->GetParameter(1),UNIT.c_str(),(100*fabs(((fMom[i][1]->GetParameter(2))/(fMom[i][1]->GetParameter(1)))))),"l");
    leg3[i]->AddEntry(MomLine[i],"1/eMC ", "l");
	leg3[i]->Draw("same");
	tps->Draw("same");

	 fileout->cd();
    MomMC[i]->Write();
    MomPR[i]->Write();
    MomReco[i]->Write();	
   } //i	   
 
 can3->Print(Form("%s/%d/Analysis_%s_%d_%s.pdf(", directory.c_str(),type, startfile.c_str(), type,ID.c_str()));
 can3->Write();
 
 
    //////////////////////////////////////////////////////////////////
TCanvas*canreso=new TCanvas("Resolution PR&RK","Resolution PR&RK",200,10,1500,1000);
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
 canreso->Print(Form("%s/%d/Analysis_%s_%d_%s.pdf)", directory.c_str(),type, startfile.c_str(), type,ID.c_str()));
 fileout->Close(); 	
	

} //end function

