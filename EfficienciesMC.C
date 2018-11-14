#include "headers.h"
#include "ALEvent.h"
#include "LoadMCparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void EfficienciesMC(int t, int nene, string s)
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

const int NReco=1;
string RecoInd[NReco]={"KFone"};

//Number of configurations to study

const int nEfficiencies=11;
string sEff[nEfficiencies] = {"TopShell","T1","T1\\&T3","T1\\&T3\\&L0","T1\\&T3\\&toL1","T1\\&T3\\&toL2","T1\\&T3\\&toL3",
							 "T1\\&T3\\&toL4","T1\\&T3\\&toL5","T1\\&T3\\&toL6","T1\\&T3\\&toL6\\&T4"};	
//string sEff[nEfficiencies] = {"T1&T4", "T1&T4&reco", "T1&T4&T3&reco","T1&T4 & 1+ layers","T1&T4 & 2+ layers", " T1&T4 & 3+ layers", "T1&T4 & 4+ layers", "T1&T4 & 5+ layers ", "T1&T4 & 5+ layers & track","T1&T4 & 6+ layers","T1&T4 & 6+ layers & track", "T1&T4 & 7+ layers", "T1&T4 7 layers & track"};	
//const int NConfig=12;
//string TConfig[NConfig]={"T1&T4", "T1&T4&reco","T1&T4&T3&reco","T1&T4&T3&reco chi>0" ,"T1&T2&T3&reco & 0+ sites used","T1&T4&T3&reco & 1+ sites","T1&T4&T3&reco & 2+ sites","T1&T4&T3&reco & 3+ sites","T1&T4&T3&reco & 4+sites","T1&T4&T3&reco & 5+ sites","T1&T4&T3&reco & 6+ sites","T1&T4&T3&reco & 7 sites"};					   
const int NConfig=4;
string TConfig[NConfig]={"T1&T3", "T1&T3&5layers+","T1&T3&6layers+","T1&T3&7layers"};					   
const int nLayers=7;
string LayerConfig[nLayers]={"L0","L1","L2","L3","L4","L5","L6"};	

//Input file 
string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4";
string startfile="aesopliteNonUniB_V4";
string endfile="_fort.99";
string source = s;
string directory= "/home/smechbal/Documents/AESOPLITE/Analysis/MCAnalysis/Efficiencies/V4";

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
float Gsource =	8882.64;	//geometry factor for 30 cm beam source, thetaMax = 90 deg
//float Gsource = 5476.16;		//geometry factor for 50 cm beam, thetaMax 28.1 deg
//float Gsource =	24674.01;		//geometry factor for 50 cm beam source, thetaMax = 90 deg

//Number of energies
int Nene=nene;
//Energies
int*Ene=new int[Nene];
//Number of cycles per energy
 int* Ncycles=new int[Nene];
 for(int i=0;i<Nene;i++)Ncycles[i]=10;

if(t==3 || t==4)
{ 
 
 Ene[0]= 10;
 Ene[1]= 20;
 Ene[2]= 30;
 Ene[3]= 40;
 Ene[4]= 50;
 Ene[5]= 60;
 Ene[6]= 100;
 Ene[7]= 300;
 }
 string UNIT="MeV";	
 int Nbins=500;
 float EneMin=0;
 float EneMax=500;
 float Dmin=-2.;
 float Dmax=2.;
  
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
 int mcolor[8]={432-3,400-9,616-9,750,900,800,920,840}; 
//zz0
float zz0=0;	
//Number of entries per reconstruction type per energy 

int** nevents= new int*[NReco];
int** nT1T4= new int*[NReco];
int** nT1T3= new int*[NReco];

//Rate Histograms
TH1F***nEff= new TH1F**[NReco];	

//Reco efficiency histograms
TH1F***nEffReco= new TH1F**[NReco];
//Scattering histograms
TH1F ****ScatNB=new TH1F***[NReco];
TH1F ****ScatB=new TH1F***[NReco];	
//Energy Loss histograms
TH1F ****ELoss=new TH1F***[NReco];

//Summary graph geo factor
TGraphErrors *GeoFac = new TGraphErrors();
TGaxis::SetMaxDigits(2);	

 for(int i=0;i<NReco;i++) 
   {
	nEff[i]=new TH1F*[Nene];
	nEffReco[i]=new TH1F*[Nene];
	ScatNB[i]=new TH1F**[Nene];
	ScatB[i]=new TH1F**[Nene];
	ELoss[i]=new TH1F**[Nene];
	nevents[i] = new int[Nene];
	nT1T4[i] = new int[Nene];
	nT1T3[i] = new int[Nene];
	for(int j=0; j<Nene; j++) 
		{	
			nevents[i][j]=0;
			nT1T4[i][j]=0;
			nT1T3[i][j]=0;
			//Rate histogram
			nEff[i][j] = new TH1F(Form("%s, %d%s, Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s, %d%s, Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),nEfficiencies,0,nEfficiencies);
			nEff[i][j]->SetStats(0);
			nEff[i][j]->GetXaxis()->SetAlphanumeric();
			for (int l=1;l<=nEfficiencies;l++) nEff[i][j]->GetXaxis()->SetBinLabel(l,sEff[l-1].c_str());
			//Reconstruction Rate histogram
			nEffReco[i][j] = new TH1F(Form("%s, %d%s,Reconstruction Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s, %d%s, Reconstruction Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),NConfig,0,NConfig);
			//nEffReco[i][j]->GetYaxis()->SetMaxDigits(3);			
		    nEffReco[i][j]->SetStats(0);
			nEffReco[i][j]->GetXaxis()->SetAlphanumeric();
			for (int l=1;l<=NConfig;l++) nEffReco[i][j]->GetXaxis()->SetBinLabel(l,TConfig[l-1].c_str());
			ScatNB[i][j]=new TH1F*[nLayers];
			ScatB[i][j]=new TH1F*[nLayers];
			ELoss[i][j]=new TH1F*[nLayers];
		    for(int k=0; k<nLayers; k++) 
			{
			if((t==3 || t==4) && Ene[j]>=60)Nbins=2*Nbins;
			if((t==3 || t==4) && Ene[j]>=100)Nbins=3*Nbins;
			ScatNB[i][j][k]=new TH1F(Form("%s %d%s, L0 - %s, ",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),Form("%s %d%s, L0 - %s, ",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),500,-1500,1500); 
			ScatNB[i][j][k]->GetXaxis()->SetTitle("Scattering in non-bending plane (in mrad)");   
			ScatNB[i][j][k]->GetYaxis()->SetTitle("#entries");
			ScatNB[i][j][k]->Sumw2(true);	
			ScatB[i][j][k]=new TH1F(Form("%s %d%s, L0 - %s",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),Form("%s %d%s, L0 -%s",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),500,-1500,1500); 
			ScatB[i][j][k]->GetXaxis()->SetTitle("Scattering in bending plane (in mrad)");   
			ScatB[i][j][k]->GetYaxis()->SetTitle("#entries");
			ScatB[i][j][k]->Sumw2(true);	
			ELoss[i][j][k] =new TH1F(Form("%s %d%s, %s, ELoss",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),Form("%s %d%s, %s, ELoss",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),100,0,20); 
			ELoss[i][j][k]->GetXaxis()->SetTitle("Eloss (in MeV)");   
			ELoss[i][j][k]->GetYaxis()->SetTitle("#entries");
			ELoss[i][j][k]->Sumw2(true);	

 				} // end k
    		} // end j
 	} // end i
cout <<"done creating histograms"<< endl;
//output ROOT file  
TFile **fileout = new TFile*[NReco];
TChain ***chain=new TChain**[NReco]; 

for(int i=0; i<NReco; i++) 
	{
	fileout[i]=new TFile(Form("%s/%s/EfficienciesMC_%s_%d_%s.root", directory.c_str(),source.c_str(),startfile.c_str(), type,RecoInd[i].c_str()),"RECREATE");

	//Create TChain to merge the files for a given energy
	chain[i] =new TChain*[Nene]; 
  	for(int j=0; j<Nene; j++) 
  		{	  	
   		cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT <<endl;
   		chain[i][j]=new TChain("MC");
        chain[i][j]->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s0001%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene[j],UNIT.c_str(),endfile.c_str(),RecoInd[i].c_str()));
//       chain[i][j]->Add("/home/sarah/AESOPLITE/MCProduction/Detector/ProcessedFiles/3/30cmSourceNoB90Deg/RecoEvent_aesopliteNonUniB_V4_3_20MeV0001_fort.99_KFone.root");

		//Define variables to read event


   //Define variables to read event
   ALEvent *e = new ALEvent();      
   //Set address to access event data
   chain[i][j]->SetBranchAddress("Revent",&e); 
 
   // Get number of event in Tree
   int nentries=chain[i][j]->GetEntries();

   cout << "Number  of events: " << nentries << endl;  
   for (int k=0;k<nentries;k++)
   {
	nevents[i][j]++;
	chain[i][j]->GetEntry(k); //Load the entry i in the variable e 
	if(k%100000==0)   {
		cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT << "  Event: " << k <<endl;
		//cout << " eMC = " << e->get_EkMC() << endl;
		}
	uint8_t Ti=(uint8_t)e->get_Ti();
    int  nhits = e->get_Nhits();
	//check trigger requirements
	bool TShell;
	bool T1=e->get_T1();
	bool T2=e->get_T2();
	bool T3=e->get_T3();
	bool T4=e->get_T4();
	bool guard=e->get_guard();
	double EFoam=0;
	double EShell=0;
	double tFoam=0;
	double tShell=0;
	double E1=0;
	double tT1=0;
	double E2=0;
	double E3=0;
	double E4=0;
	//Insulation Foam
	int nTFoam=e->get_EneIsofoam().size();
	int nTShell=e->get_EneShell().size();
	if(nTFoam>0 && nTShell>0)TShell=true;
	else TShell=false;	   
	//Number of layers with hit(s) in bending/non-bending plane
	int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	//Number of layers with hit(s)
	int NL= e->get_NLayers();
	bool Lhit[7] = {false,false,false,false,false,false,false};
//	cout << "Event " << k << ", total number of hits = " << nhits << endl;
//	cout << " NL = " << NL << " , NLB = " << NLB << " , NLNB = " << NLNB << endl;
	//Reconstruction variables
	float pPR=1000*fabs(e->get_p0PR());   //in MeV
	float ePR=e->get_EkPR();
	double deflecPR = e->get_deflecPR();	  
	double chiNBPR = e->get_chi2NBPR();
	double chiBPR = e->get_chi2BPR();
	double Ekreco = e->get_Ekreco();			//in MeV
	double preco = e->get_p0reco();
     double ndf = e->get_ndf();
	double chireco = e->get_chi2();
	double phi0 = e->get_phi0();
	double tanl = e->get_tanl();
	double cpa_reco = e->get_cpa();
	double sigma_cpa = TMath::Sqrt(e->get_cpaerr2());	//in GeV
	double e0MC;
	double p0MC;
	double cpa_MC;
	int typeMC = e->get_typeMC();
	int QMC;
	if(typeMC==3 || typeMC==11) QMC = -1;
	if(typeMC==4 || typeMC==10) QMC = 1;
	int sitesUsed=0;
	double thetaNBL0=0;
	double thetaBL0=0;
    double thetaNB=0;
	double thetaB=0;
 //get MC informations at every layers
	for(int l=0;l<e->get_Nhits();l++) 
     {		
		 bool used = e->get_hits().at(l)->get_fUsed();
		 if(used) sitesUsed++;
	     int Lindex=e->get_hits().at(l)->get_L();
	     Lhit[Lindex] = true;
         double cx = e->get_hits().at(l)->get_cx();
         double cy = e->get_hits().at(l)->get_cy();
         double cz = e->get_hits().at(l)->get_cz();
		 if(Lindex==0) {
			thetaNBL0 = TMath::ATan(cx/cz);
			thetaBL0 =  TMath::ATan(cy/cz);
		 	}
		 thetaNB = TMath::ATan(cx/cz);
		 thetaB =  TMath::ATan(cy/cz);
		 double eMC = 1000*(e->get_hits().at(l)->get_eMC());
		 double e0MC = 1000*(e->get_EkMC());
         if(cx!=0 || cy!=0 || cz!=0 ) {	 
			ScatNB[i][j][Lindex]->Fill(1000*(thetaNBL0-thetaNB));
			ScatB[i][j][Lindex]->Fill(1000*(thetaBL0-thetaB));
		 	}
		 ELoss[i][j][Lindex]->Fill(e0MC-eMC);		 
	} //end loop on hits
	   	
	 //Fill counter histograms
	 if(TShell){nEff[i][j]->Fill(sEff[0].c_str(),1);}
	 if(TShell && T1){nEff[i][j]->Fill(sEff[1].c_str(),1);} 
	 if(TShell && T1 && T3){nEff[i][j]->Fill(sEff[2].c_str(),1);} 
	 if(TShell && T1 && T3 && Lhit[0]){nEff[i][j]->Fill(sEff[3].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1]){nEff[i][j]->Fill(sEff[4].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2]){nEff[i][j]->Fill(sEff[5].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3]){nEff[i][j]->Fill(sEff[6].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4]){nEff[i][j]->Fill(sEff[7].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5]){nEff[i][j]->Fill(sEff[8].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5] && Lhit[6]){nEff[i][j]->Fill(sEff[9].c_str(),1);}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5] && Lhit[6] && T4){nEff[i][j]->Fill(sEff[10].c_str(),1);}

   if(T1 && T3 ){
	   nEffReco[i][j]->Fill(TConfig[0].c_str(),1);
	   nT1T3[i][j]++;} 	
	 if(T1 && T3 && NL>=5){nEffReco[i][j]->Fill(TConfig[1].c_str(),1);} 
	 if(T1 && T3 && NL>=6){nEffReco[i][j]->Fill(TConfig[2].c_str(),1);}
	 if(T1 && T3 && NL>=7){nEffReco[i][j]->Fill(TConfig[3].c_str(),1);}

	    } // end k, number of entries 
	} // end j, energies reconstructed
    } // end i, reco type

//////////////////////////////
// Display histograms////////
/////////////////////////////   

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);	
gROOT->SetStyle("Modern");
//////////////////////////////
////Efficiencies/////////////
/////////////////////////////
cout << "start display canvas" << endl;
TCanvas**CanEff=new TCanvas*[NReco]; 
TLegend**LegEff=new TLegend*[NReco];
TText****txt=new TText***[NReco];
double barwidth=0.1;
double baroffset=0.1;
	
 for(int i=0;i<NReco;i++)
   {
     CanEff[i] = new TCanvas(Form("%s, efficiencies", RecoInd[i].c_str()),Form("%s, efficiencies", RecoInd[i].c_str()),200,10,1200,800);
     txt[i]=new TText**[Nene];
	 LegEff[i] = new TLegend(0.8,0.8,0.9,0.9);
     CanEff[i]->cd();
     CanEff[i]->SetBottomMargin(0.32);
    for (int j=0; j<Nene;j++) 
	{
	   gPad->SetGridy(1);
	   nEff[i][0]->SetTitle(Form("%s, Trigger Rates",source.c_str()));
   	   nEff[i][j]->SetLineWidth(0);
   	   nEff[i][j]->SetMarkerStyle(20);
       nEff[i][j]->SetFillColor(mcolor[j]);
   	   nEff[i][j]->SetMarkerSize(0);
       nEff[i][j]->Scale(100./nevents[i][j]);
       nEff[i][j]->SetBarWidth(barwidth);
	   nEff[i][j]->SetBarOffset(baroffset+j*barwidth);
       nEff[i][j]->LabelsOption("v","X"); 
       nEff[i][j]->GetYaxis()->SetRangeUser(0,20);
       nEff[i][j]->GetYaxis()->SetNdivisions(511,0);
   	   nEff[i][j]->GetYaxis()->SetTitle("Rates (normalized to total # events)");
   	   if (j==0) { 
		nEff[i][j]->Draw("bar");
	        nEff[i][j]->SetTitle(Form("%s",RecoInd[i].c_str()));
		}
	 else {
		nEff[i][j]->Draw("bar same");
	  	}
       LegEff[i]->AddEntry(nEff[i][j],Form("%d%s",Ene[j],UNIT.c_str()),"f");
	   txt[i][j]=new TText*[nEfficiencies];
	  for(int l=0;l<nEfficiencies;l++)
     		 {
	   if(l>1){txt[i][j][l]=new TText((nEff[i][j]->GetXaxis()->GetBinCenter(l+1)+(j-1)*barwidth)-0.15,nEff[i][j]->GetBinContent(l+1)+1,Form("%2.2f cm2 sr +/- %2.2f cm2 sr", Gsource*(nEff[i][j]->GetBinContent(l+1)/100), Gsource*(nEff[i][j]->GetBinError(l+1)/100))); 
       txt[i][j][l]->SetTextSize(0.010);
       txt[i][j][l]->SetTextAngle(55);    
//       txt[i][j][l]->Draw();  
			  }
      }  // end l
   } //end j
	 
	 LegEff[i]->Draw();
     fileout[i]->cd();
     gPad->Update();
     CanEff[i]->Write();
    if (i==0) CanEff[i]->Print(Form("%s/%s/Efficiencies_%d_%s_%s.pdf(", directory.c_str(),source.c_str(),type,source.c_str(),RecoInd[i].c_str()));
	} //end i

	
for(int l=0;l<nEfficiencies;l++) {
    	cout << Form("%s& ", sEff[l].c_str());
	for (int j=0; j<Nene;j++) {
		 for(int i=0;i<NReco;i++) {		
		 cout << Form("%2.2f",Gsource*(nEff[i][j]->GetBinContent(l+1)/100));
		 if(j!=Nene-1){cout << " & " ;}
		 //cout << Form("%2.2f",Gsource*(nEff[i][j]->GetBinContent(l+1)/100)) << " $\\pm$ " << Form("%2.2f",Gsource*(nEff[i][j]->GetBinError(l+1)/100)) << " $cm^2 sr$ & " ;
		 if(l==nEfficiencies-2) {
			 GeoFac->SetPoint(j+1,Ene[j],Gsource*(nEff[i][j]->GetBinContent(l+1)/100));
			 GeoFac->SetPointError(j+1,0.0,Gsource*(nEff[i][j]->GetBinError(l+1)/100));
		 }
		 } 
		}
	cout << " \\\\" << endl;
 	}	
	
//////////////////////////
//Reconstruction Rates////
/////////////////////////

TCanvas**CanEffReco=new TCanvas*[NReco];
TLegend**LegEffReco=new TLegend*[NReco];
TText****txtReco=new TText***[NReco];
 for(int i=0;i<NReco;i++)
   {
         //  CanEffReco[i]= new TCanvas(Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),200,10,1200,800);
           CanEffReco[i]= new TCanvas(Form("%s,geometric efficiencies", RecoInd[i].c_str()),Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),200,10,1200,800);
	       txtReco[i]=new TText**[Nene];
           LegEffReco[i] = new TLegend(0.8,0.8,0.9,0.9);
           CanEffReco[i]->cd();
           CanEffReco[i]->SetBottomMargin(0.32);

    for (int j=0; j<Nene;j++)
        {
           gPad->SetGridy(1);
           nEffReco[i][0]->SetTitle(Form("%s, %s efficiencies",source.c_str(),RecoInd[i].c_str()));
           nEffReco[i][j]->SetLineWidth(0);
           nEffReco[i][j]->SetMarkerStyle(20);
           nEffReco[i][j]->SetFillColor(mcolor[j]);
           nEffReco[i][j]->SetMarkerSize(0);
           nEffReco[i][j]->Scale(100./nT1T3[i][j]);
           nEffReco[i][j]->SetBarWidth(barwidth);
           nEffReco[i][j]->SetBarOffset(baroffset+j*barwidth);
           nEffReco[i][j]->LabelsOption("v","X");
           nEffReco[i][j]->GetYaxis()->SetRangeUser(0,110);
           nEffReco[i][j]->GetYaxis()->SetNdivisions(511,0);
           nEffReco[i][j]->GetYaxis()->SetTitle("Normalized to T1/T3 rates");
 if (j==0) {
                nEffReco[i][j]->Draw("bar");
                nEffReco[i][j]->SetTitle(Form("%s reconstruction efficiencies",RecoInd[i].c_str()));
                }
         else {
                nEffReco[i][j]->Draw("bar same");
                }
           LegEffReco[i]->AddEntry(nEffReco[i][j],Form("%d%s",Ene[j],UNIT.c_str()),"f");
           txtReco[i][j]=new TText*[NConfig];
          for(int l=0;l<NConfig;l++)
                 {
    //   if(nEffReco[i][j]->GetBinContent(l+1)>=100&&l==0) txtReco[i][j][l]=new TText(nEffReco[i][j]->GetXaxis()->GetBinCenter(l+1),nEffReco[i][j]->GetBinContent(l+1)+1,Form("%d",(int)nEffReco[i][j]->GetBinContent(l+1)));
        txtReco[i][j][l]=new TText((nEffReco[i][j]->GetXaxis()->GetBinCenter(l+1)+(j-1)*barwidth)-0.15,nEffReco[i][j]->GetBinContent(l+1)+1,Form("%.1f",nEffReco[i][j]->GetBinContent(l+1)));
       txtReco[i][j][l]->SetTextSize(0.015);
       txtReco[i][j][l]->SetTextAngle(55);
      // txtReco[i][j][l]->Draw();
      }  // end l
   } //end j
     fileout[i]->cd();
     LegEffReco[i]->Draw();
     gPad->Update();
     CanEffReco[i]->Write();
     CanEffReco[i]->Print(Form("%s/%s/Efficiencies_%d_%s_%s.pdf", directory.c_str(),source.c_str(),type,source.c_str(),RecoInd[i].c_str()));
        } //end i


/////////////////////////////////
/////////GeoFactor plot/////////
///////////////////////////////
 for(int i=0;i<NReco;i++)
   {
	TCanvas *canGeo=new TCanvas(Form("Geo Factor %s",source.c_str()),Form("Geo Factor %s",source.c_str()),200,10,800,1500);
	canGeo->cd(1);
	GeoFac->GetYaxis()->SetTitleOffset(1.4);
	GeoFac->GetYaxis()->SetTitle("Geometry factor (cm^{2} sr)");
	GeoFac->GetXaxis()->SetTitle(Form("Energy (%s)", UNIT.c_str()));
	GeoFac->SetTitle(Form("Geometry Factor %s",source.c_str()));
	GeoFac->GetYaxis()->SetLimits(0,25);
	GeoFac->GetXaxis()->SetLimits(0,350);
	GeoFac->SetMarkerStyle(20);
	GeoFac->SetMarkerSize(1);	 
	GeoFac->Draw("aple");
	fileout[i]->cd();
	gPad->Update();
	canGeo->Write();
        GeoFac->Write();
	GeoFac->Print(Form("%s/%s/Efficiencies_%d_%s_%s.pdf)", directory.c_str(),source.c_str(),type,source.c_str(),RecoInd[i].c_str()));
 }

/////////////////////////////////////
///////Scattering distribution/////////
///////////////////////////////////

 gStyle->SetOptStat(11111111);
 TCanvas****cancos=new TCanvas***[NReco];
 for(int i =0; i<NReco;i++) { 
 	 fileout[i]->cd();
	 cancos[i] = new TCanvas**[Nene];
        for(int j=0;j<Nene;j++) {
            cancos[i][j]=new TCanvas*[nLayers];
	         for(int k=0; k<nLayers;k++) {	 
 		      cancos[i][j][k]=new TCanvas(Form("%s %d%s,scattering", LayerConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s, scattering", LayerConfig[k].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
			  cancos[i][j][k]->SetTitle(Form("L0-%s",LayerConfig[k].c_str()));
              cancos[i][j][k]->Divide(2,1);
              cancos[i][j][k]->cd(1);
      	      ScatNB[i][j][k]->Draw("hist");
              cancos[i][j][k]->cd(2);
      	      ScatB[i][j][k]->Draw("hist");             			 
	      cancos[i][j][k]->Write();
       		 } //end k
  			} //end j
		  } //end i


/////////////////////////////////////
///////ELoss distribution/////////
///////////////////////////////////

 gStyle->SetOptStat(11111111);
 TCanvas***canLoss=new TCanvas**[NReco];
 for(int i =0; i<NReco;i++) { 
 	 fileout[i]->cd();
	 canLoss[i] = new TCanvas*[Nene];
        for(int j=0;j<Nene;j++) {
 		     canLoss[i][j]=new TCanvas(Form("%d%s, ELoss",Ene[j],UNIT.c_str()),Form("%d%s, ELoss",Ene[j],UNIT.c_str()),200,10,1200,800);
		     canLoss[i][j]->Divide(2,3);
	         for(int k=0; k<nLayers;k++) {	 
			  ELoss[i][j][k]->SetTitle(Form("ELoss from injection point to %s",LayerConfig[k].c_str()));
              canLoss[i][j]->cd(k+1);
      	      ELoss[i][j][k]->Draw("hist");			 
       		 } //end k
	      canLoss[i][j]->Write();
  			} //end j
	      fileout[i]->Close();	  
	} //end i

} //end function






