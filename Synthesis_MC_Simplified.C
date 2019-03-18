#include "headers.h"
#include "ALEvent.h"
#include "LoadMCparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void SynthesisSimplified(int t, int nene, int ncycles)
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
  double B = 0.3;               //average magnetic field in T
  double c = TMath::C();
 string MCparamfile="./MCparameters.dat";
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);


//Number of reconstruction types to study 

const int NReco=3;
string RecoInd[NReco]={"KFone","KFtwo", "RKfit"};

//Number of configurations to study
const int NConfig=1;
string TConfig[NConfig]={"T1&T3&reco"};					   
const int nEfficiencies=13; 	 		//count and record number of entries for all events, T1&T4 no hits, T1&T4 1 layer with hit, T1&T4 2 layer with hits, etc
string sEff[nEfficiencies] = {"T1&T3","T1&T3&hits","T1&T3 & 1+ layers","T1&T3 & 2+ layers", " T1&T3 & 3+ layers", "T1&T3 & 4+ layers", "T1&T3 & 5+ layers ", "T1T3 & 5+ layers & track","T1&T3 & 6+ layers","T1&T4 & 6+ layers & track", "T1&T4 & 7+ layers", "T1&T4 7 layers & track","T1&T4&reco"};	
//Input file 
string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4";	
string source="30cmSourceYesB90Deg";
string startfile="aesopliteNonUniB_V4";
string endfile="_fort.99";
string directory= "/home/smechbal/Documents/AESOPLITE/Analysis/MCAnalysis/MC_Synthesis/V4";


//type of particle
int type=t;
string stype;
if(t==3)        stype="e^{#minus}";
if(t==4)        stype="e^{#plus}";
if(t==11)       stype="#mu^{#minus}";
if(t==10)       stype="#mu^{#plus}";
if(t==1)        stype="protons";
if(t==6)        stype="#alpha";

float mass=0.000511;//electon mass in GeV
if(t==11 || t==10)     mass=0.10566;//muon mass in GeV
if(t==1)     mass=0.93827;//proton mass in GeV
if(t==6)     mass=3.7273;//alpha-particle mass in GeV

if(t==11 || t==10)      mass=0.10566;//muon mass in GeV
//Number of energies
int Nene=nene;
//Energies
int*Ene=new int[Nene];
//Number of cycles per energy
 int* Ncycles=new int[Nene];
 for(int i=0;i<Nene;i++)Ncycles[i]=ncycles;


if(t==3 || t==4)
{

 Ene[0]= 10;
 Ene[1]= 20;
 Ene[2]= 30;
 Ene[3]= 40;
 Ene[4]= 50;
 Ene[5]= 60;
 Ene[6]= 80;
 Ene[7]= 100;
 Ene[8]= 200;
 Ene[9]= 300;
 Ene[10]= 400;
 Ene[11]= 500;
 }

 string UNIT="MeV";	
 int Nbins=500;
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
 int mcolor[4]={432-3,400-9,616-9,632}; 
//zz0
float zz0=0;	
//Number of entries per reconstruction type per energy 
int** nevents= new int*[NReco];

//Rate Histograms
TH1F***nEff= new TH1F**[NReco];	

//Reco efficiency histograms
TH1F***nEffReco= new TH1F**[NReco];
//Histograms MC inverse energy/momentum at L0
TH1F ****pMC=new TH1F***[NReco];
TH1F ****p0PR=new TH1F***[NReco];	
TH1F ****p0reco=new TH1F***[NReco];	

//Momentum pull distributions
TH1F 	***ResoReco=new TH1F**[NReco];
//Chi2 distributions
TH1F ****chi2reco=new TH1F***[NReco];
TH1F ****chi2BPR=new TH1F***[NReco];
TH1F ****chi2NBPR=new TH1F***[NReco];

  //TeX output file
ofstream totex;
totex.open(Form("%s/%d/MC_Synthesis_Simplified_TeX.txt",directory.c_str(),type));
  
 for(int i=0;i<NReco;i++) 
   {
	nEff[i]=new TH1F*[Nene];
	nEffReco[i]=new TH1F*[Nene];
	p0PR[i]=new TH1F**[Nene];
	p0reco[i]=new TH1F**[Nene];
	pMC[i]=new TH1F**[Nene];
	ResoReco[i]=new TH1F*[Nene];
	chi2reco[i]=new TH1F**[Nene];
  chi2BPR[i]=new TH1F**[Nene];
  chi2NBPR[i]=new TH1F**[Nene];
	nevents[i] = new int[Nene];
	for(int j=0; j<Nene; j++) 
		{	
			nevents[i][j]=0;
//Rate histogram
      nEff[i][j] = new TH1F(Form("%s, %d%s, Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s, %d%s, Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),nEfficiencies,0,nEfficiencies);
      nEff[i][j]->SetStats(0);
      nEff[i][j]->GetXaxis()->SetAlphanumeric();
      for (int l=1;l<=nEfficiencies;l++) nEff[i][j]->GetXaxis()->SetBinLabel(l,sEff[l-1].c_str());
//Reconstruction Rate histogram
      nEffReco[i][j] = new TH1F(Form("%s, %d%s,Reconstruction Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s, %d%s, Reconstruction Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str()),NConfig,0,NConfig);
      nEffReco[i][j]->SetStats(0);
      nEffReco[i][j]->GetXaxis()->SetAlphanumeric();
      for (int l=1;l<=NConfig;l++) nEffReco[i][j]->GetXaxis()->SetBinLabel(l,TConfig[l-1].c_str());
      p0PR[i][j]=new TH1F*[NConfig];
      p0reco[i][j]=new TH1F*[NConfig];
      pMC[i][j]=new TH1F*[NConfig];
		  chi2reco[i][j]=new TH1F*[NConfig];
		  chi2BPR[i][j]=new TH1F*[NConfig];
      chi2NBPR[i][j]=new TH1F*[NConfig];
			//Momentum pull histograms
		  ResoReco[i][j]=new TH1F(Form("Resolution momentum from reconstruction, %d%s",Ene[j],UNIT.c_str()),Form("Resolution momentum from reconstruction, %d%s",Ene[j],UNIT.c_str()),10000,-1,1);
			ResoReco[i][j]->GetYaxis()->SetTitle("#entries");
   		ResoReco[i][j]->GetXaxis()->SetTitle(Form("Inverse momentum pull #frac{#kappa_{reco}-#kappa_{truth}}{#sigma_{#kappa}}"));
			ResoReco[i][j]->Sumw2(true);
			ResoReco[i][j]->SetTitle(Form("Momentum Pull: %s %d%s",RecoInd[i].c_str(), Ene[j], UNIT.c_str()));
		    for(int k=0; k<NConfig; k++) 
			{
        if((t==3 || t==4) && Ene[j]>=60)Nbins=2*Nbins;
        if((t==3 || t==4) && Ene[j]>=100)Nbins=3*Nbins;
        p0PR[i][j][k] =new TH1F(Form("%s, %s, 1/p_{0PR}, %d%s",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s, %s, 1/p_{0PR}, %d%s",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),300,1./EneMin,1./EneMax); 
        p0PR[i][j][k]->GetXaxis()->SetTitle(Form("1/p_{0PR} (%s^{-1})",UNIT.c_str()));   
        p0PR[i][j][k]->GetYaxis()->SetTitle("#entries");
        p0PR[i][j][k]->Sumw2(true);	
        p0reco[i][j][k] =new TH1F(Form("%s, %s, 1/p_{0reco}, %d%s",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s, %s, 1/p_{0reco}, %d%s",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),300,1./EneMin,1./EneMax); 
        p0reco[i][j][k]->GetXaxis()->SetTitle(Form("1/p_{0reco} (%s^{-1})",UNIT.c_str()));   
        p0reco[i][j][k]->GetYaxis()->SetTitle("#entries");
        p0reco[i][j][k]->Sumw2(true); 	
        pMC[i][j][k] =new TH1F(Form("%s, %s, 1/p_{0MC}, %d%s",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s, %s, 1/p_{0MC}, %d%s",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),1000,1./EneMin,1./EneMax); 
        pMC[i][j][k]->GetXaxis()->SetTitle(Form("1/p_{0MC} (%s^{-1})",UNIT.c_str()));   
        pMC[i][j][k]->GetYaxis()->SetTitle("#entries");
        pMC[i][j][k]->Sumw2(true); 
        chi2reco[i][j][k] =new TH1F(Form("%s, %s, %d%s chi2reco",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s, %s, %d%s chi2reco",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),60,0,30);
        chi2reco[i][j][k]->GetXaxis()->SetTitle("chi2");
        chi2reco[i][j][k]->GetYaxis()->SetTitle("#entries");
        chi2BPR[i][j][k] =new TH1F(Form("%s, %s, %d%s chi2BPR",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s, %s, %d%s chi2BPR",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),100,0,100);
        chi2BPR[i][j][k]->GetXaxis()->SetTitle("chi2BPR");
        chi2BPR[i][j][k]->GetYaxis()->SetTitle("#entries");
        chi2NBPR[i][j][k] =new TH1F(Form("%s, %s, %d%s chi2NBPR",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),Form("%s, %s, %d%s chi2NBPR",RecoInd[i].c_str(),TConfig[k].c_str(),Ene[j],UNIT.c_str()),100,0,100);
        chi2NBPR[i][j][k]->GetXaxis()->SetTitle("chi2BPR");
        chi2NBPR[i][j][k]->GetYaxis()->SetTitle("#entries");
        chi2NBPR[i][j][k]->SetTitle(Form("chi2NBPR, %s %d%s %s",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),TConfig[k].c_str()));
        chi2BPR[i][j][k]->SetTitle(Form("chi2BPR, %s %d%s %s",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),TConfig[k].c_str()));
        chi2reco[i][j][k]->SetTitle(Form("chi2reco, %s %d%s %s",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),TConfig[k].c_str()));
	 
              } // end k
                } // end j
          } // end i
  
  
cout <<"done creating histograms"<< endl;
//output ROOT file  
TFile **fileout = new TFile*[NReco];
TChain ***chain=new TChain**[NReco]; 

for(int i=0; i<NReco; i++) 
	{
	fileout[i]=new TFile(Form("%s/%d/Synthesis_Simplified%s_%d_%s.root", directory.c_str(),type,startfile.c_str(), type,RecoInd[i].c_str()),"RECREATE");
	//Create TChain to merge the files for a given energy
	chain[i] =new TChain*[Nene]; 
  	for(int j=0; j<Nene; j++) 
  		{	  	
  // 		cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT <<endl;
   		chain[i][j]=new TChain("MC");
		for(int k=0;k<Ncycles[j];k++)//Number of cycles
			{
        		chain[i][j]->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s%04d%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene[j],UNIT.c_str(),k+1,endfile.c_str(),RecoInd[i].c_str()));

			} // end k


   //Define variables to read event
   ALEvent *e = new ALEvent();      
   //Set address to access event data
   chain[i][j]->SetBranchAddress("Revent",&e); 
 
   // Get number of event in Tree
   int nentries=chain[i][j]->GetEntries();

   cout << "Number  of events: " << nentries << endl;  
   for (int k=0;k<nentries;k++)
   {
	chain[i][j]->GetEntry(k); //Load the entry i in the variable e 
	if(k%100000==0)   cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT << "  Event: " << k <<endl;	 
	uint8_t Ti=(uint8_t)e->get_Ti();
        int  nhits = e->get_Nhits();
	//Number of layers with hit(s)
	int NL= e->get_NLayers();
	//check trigger requirements
	bool T1=e->get_T1();
	bool T3=e->get_T3();
	bool T4=e->get_T4();
	//Number of layers with hit(s) in bending/non-bending plane
	int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	float pPR=1000*fabs(e->get_p0PR());   //in MeV
	float ePR=e->get_EkPR();
	double deflecPR = e->get_deflecPR();	  
	double chiNBPR = e->get_chi2NBPR();
	double chiBPR = e->get_chi2BPR();
	double Ekreco = e->get_Ekreco();			//in MeV
	double pReco = e->get_p0reco();
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
     //get MC informations at L0
        for(int l=0;l<e->get_Nhits();l++) 
	  {
             bool used = e->get_hits().at(l)->get_fUsed();
             if(used) sitesUsed++;
	     int Lindex=int(e->get_hits().at(l)->get_mregMC())%11;
	     //LAYER L0, NB plane		
	     if(Lindex==0) 
		 {								
		 //truth variables at L0
		   e0MC = e->get_hits().at(l)->get_eMC();		
		   p0MC = TMath::Sqrt((e0MC*e0MC) - (mass*mass));  //momentum at L0 in GeV
       cpa_MC = QMC/p0MC;
		   p0MC=1000*p0MC;
		 }
          }
//Fill counter histograms
         if(T1 && T3){
		nEff[i][j]->Fill(sEff[0].c_str(),1);
         	nevents[i][j]++;
		   }
        if(T1 && T3 && nhits>0){nEff[i][j]->Fill(sEff[1].c_str(),1);} 
	if(T1 && T3 && nhits>0 && NL>=1){nEff[i][j]->Fill(sEff[2].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL>=2){nEff[i][j]->Fill(sEff[3].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL>=3){nEff[i][j]->Fill(sEff[4].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL>=4){nEff[i][j]->Fill(sEff[5].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL>=5){nEff[i][j]->Fill(sEff[6].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL>=5 && chiNBPR>=0 && chiBPR>=0){nEff[i][j]->Fill(sEff[7].c_str(),1);}	
         if(T1 && T3 && nhits>0 && NL>=6){nEff[i][j]->Fill(sEff[8].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL>=6 && chiNBPR>=0 && chiBPR>=0){nEff[i][j]->Fill(sEff[9].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL==7){nEff[i][j]->Fill(sEff[10].c_str(),1);}
         if(T1 && T3 && nhits>0 && NL==7 && chiNBPR>=0 && chiBPR>=0){nEff[i][j]->Fill(sEff[11].c_str(),1);}
         if(T1 && T3 && nhits>0 && pReco>0){nEff[i][j]->Fill(sEff[12].c_str(),1);}

//Fill histograms

	if(T1 && T3 && pReco>0) 
	  {
	  nEffReco[i][j]->Fill(TConfig[0].c_str(),1);
	  p0PR[i][j][0]->Fill(1/pPR);
	  p0reco[i][j][0]->Fill(1/pReco);
	  pMC[i][j][0]->Fill(1/p0MC);
	  chi2reco[i][j][0]->Fill(chireco);
	  chi2NBPR[i][j][0]->Fill(chiNBPR);
	  chi2BPR[i][j][0]->Fill(chiBPR);		
	  ResoReco[i][j]->Fill((cpa_reco-cpa_MC)/(sigma_cpa)); 
	 }
	    } // end k, number of entries 
	} // end j, energies reconstructed
    } // end i, reco type

//////////////////////////////
// Display histograms////////
/////i////////////////////////   

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);	
gROOT->SetStyle("Modern");

cout << "done processing events " << endl;
///////////////////////////
//Inverse momentum////////
/////////////////////////

TCanvas***can=new TCanvas**[NReco];
TF1*****fgauss=new TF1****[NReco];
double nsigma = 3;
for(int i=0;i<NReco;i++) 
   {
//output in TeX format, one table for each reconstruction type
	//output in TeX 
	totex << "\\begin{table}[H]"<<endl;
	totex << "\\centering" << endl;
	totex << "\\caption{"<<RecoInd[i].c_str()<<"} \\label{tab:title} \\\\" << endl;
	totex << "\\begin{tabular}{|c|c|c|c|}"<<endl;
	totex << "\\hline"<<endl;
	totex << " Energy (in MeV) & fit $\\mu$/ $\\sigma$ (MeV)^-1  & resolution & $\\chi$^2/ndf \\\\ " << endl;
	totex << "\\hline"<<endl;

 	can[i]=new TCanvas*[Nene];
	fgauss[i]=new TF1***[Nene];
	for(int j=0;j<Nene;j++)
	   {
		can[i][j]=new TCanvas(Form("%s %d%s Inverse momentum reconstruction", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s Inverse momentum reconstruction", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
	//	can[i][j]->Divide(4,3);
		fgauss[i][j]=new TF1**[NConfig];
		for(int k=0;k<NConfig;k++)
		   { 		
         fgauss[i][j][k] = new TF1*[2];
         can[i][j]->cd(k+1);
         double Max = 0;
         int newbin = TMath::Sqrt(p0reco[i][j][k]->GetEntries());
//         if(newbin>0) {
		 //double binratio = 300/newbin;
       		 // p0reco[i][j][k]->Rebin(1.3*binratio);
       		 // p0PR[i][j][k]->Rebin(1.3*binratio);
       		 // pMC[i][j][k]->Rebin(binratio);
	//	}
         double tmpreco=p0reco[i][j][k]->GetMaximum();
         double tmpPR=p0PR[i][j][k]->GetMaximum();
         if (tmpreco > tmpPR) {Max = tmpreco;}
         else {Max = tmpPR;}
         p0reco[i][j][k]->SetTitle(Form("%s %s %s %d%s",stype.c_str(),RecoInd[i].c_str(), TConfig[k].c_str(), Ene[j],UNIT.c_str()));
         p0reco[i][j][k]->SetLineColor(kBlue);
         double muReco = p0reco[i][j][k]->GetMean();
         double sigReco = p0reco[i][j][k]->GetStdDev();
         double muPR = p0PR[i][j][k]->GetMean();
         double sigPR = p0PR[i][j][k]->GetStdDev();
    //		 cout << Form("%s %s %s %d%s",stype.c_str(),RecoInd[i].c_str(), TConfig[k].c_str(), Ene[j],UNIT.c_str()) << ", meanReco = " << muReco << "  sigmaReco = " << sigReco << endl;
         p0reco[i][j][k]->Scale(0.8*Max/p0reco[i][j][k]->GetMaximum());
         pMC[i][j][k]->Scale(0.8*Max/pMC[i][j][k]->GetMaximum());
         p0PR[i][j][k]->Scale(0.8*Max/p0PR[i][j][k]->GetMaximum());
         p0reco[i][j][k]->GetXaxis()->SetRangeUser(muReco-3*sigReco,muReco+3*sigReco);
         p0reco[i][j][k]->Draw("hist");
         p0reco[i][j][k]->SetMaximum(1.1*Max);		
        // p0PR[i][j][k]->Draw("hist sames");
         p0PR[i][j][k]->SetLineColor(kRed);
         p0PR[i][j][k]->SetMarkerColor(kRed);
         pMC[i][j][k]->Draw("hist sames");
         pMC[i][j][k]->SetLineColor(kBlack);
         gPad->Update();
         TPaveStats *tpsreco = (TPaveStats*) p0reco[i][j][k]->FindObject("stats");
         tpsreco->SetName("Reco stats");
         double y1 = tpsreco->GetY1NDC();
         double y2 = tpsreco->GetY2NDC();
         tpsreco->SetTextColor(kBlue);
         tpsreco->SetLineColor(kBlue);
         /*
	 TPaveStats *tpsPR = (TPaveStats*) p0PR[i][j][k]->FindObject("stats");
         tpsPR->SetName("PR stats");
         tpsPR->SetY2NDC(y1);
         tpsPR->SetY1NDC(y1-(y2-y1));
         tpsPR->SetTextColor(kRed);
         tpsPR->SetLineColor(kRed);
         double y1PR = tpsPR->GetY1NDC();
         double y2PR = tpsPR->GetY2NDC();
	 */
         TPaveStats *tpsMC = (TPaveStats*) pMC[i][j][k]->FindObject("stats");
         tpsMC->SetName("MC truth stats");
         tpsMC->SetY2NDC(y1);
         tpsMC->SetY1NDC(y1-(y2-y1));
         tpsMC->SetTextColor(kBlack);
         tpsMC->SetLineColor(kBlack);

           for(int l=0;l<2;l++) 
            {
    if(l==0)
      {
      fgauss[i][j][k][l]=new TF1(Form("fgauss%d%d%d%d",i,j,k,l),"gaus",muReco-nsigma*sigReco,muReco+nsigma*sigReco);
      p0reco[i][j][k]->Fit(fgauss[i][j][k][l],"R0");
      fgauss[i][j][k][l]->SetLineColor(kBlue);
      }
    if(l==1)	   
            {
                  fgauss[i][j][k][l]=new TF1(Form("fgauss%d%d%d%d",i,j,k,l),"gaus",muPR-nsigma*sigPR,muPR+nsigma*sigPR);
                  p0PR[i][j][k]->Fit(fgauss[i][j][k][l],"R0");
                  fgauss[i][j][k][l]->SetLineColor(kRed);
            } 
          } //end l
 //TeX output
	totex << Ene[j] << "&" << fgauss[i][j][k][0]->GetParameter(1)<<"/"<<fgauss[i][j][k][0]->GetParameter(2)<<"&"\
<<100*fabs(((fgauss[i][j][k][0]->GetParameter(2))/(fgauss[i][j][k][0]->GetParameter(1))))
<<"&"<<fgauss[i][j][k][0]->GetChisquare()<<"/"<<fgauss[i][j][k][0]->GetNDF()<<" //" <<endl;
		
		
	fileout[i]->cd();
	p0reco[i][j][k]->Write();
	p0PR[i][j][k]->Write();
	pMC[i][j][k]->Write();
    } //end k	   
     can[i][j]->Write();
     can[i][j]->SaveAs(Form("%s/%d/Reconstructions_Simplified_%s_%d%s.eps", directory.c_str(),type,RecoInd[i].c_str(),Ene[j],UNIT.c_str()));
    if (i==0&&j==0) can[i][j]->Print(Form("%s/%d/Reconstructions_Simplified.pdf(", directory.c_str(),type));
    else if (i==NReco-1&&j==Nene-1) can[i][j]->Print(Form("%s/%d/Reconstructions_Simplified.pdf)",directory.c_str(),type));
    else can[i][j]->Print(Form("%s/%d/Reconstructions_Simplified.pdf", directory.c_str(),type));
  } //end j
	
	totex <<" \\hline" <<endl;
	totex <<"end{tabular}\\par"<< endl;
	totex<< "\\end{table}"<<endl;
} //end i

/////////////////////////////////////
///chi2 distribution canvas/////////
////////////////////////////////////

cout << "chi2 histograms" << endl;

 gStyle->SetOptStat(11111111);
 TCanvas***canchi2=new TCanvas**[NReco];
 for(int i =0; i<NReco;i++) { 
 	canchi2[i] = new TCanvas*[Nene];
        for(int j=0;j<Nene;j++) {
            canchi2[i][j]=new TCanvas(Form("%s %d%s Chi2 distribution", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s Chi2 distribution", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
	    for(int k=0; k<NConfig;k++) {	 
            	canchi2[i][j]->cd(k+1);
      	    //  gPad->SetLogx();
      	    //  fitchi2[i] = new TF1(Form("fchi2%d",i),"ROOT::Math::chisquared_pdf(x,1)",0,1000);
	    //  fitchi2[i]->SetParameter(0,7);
      	    //  fitchi2[i]->SetLineColor(kRed);
	    //  fitchi2[i]->SetParameter(1,chi2[i]->GetEntries());
     	    //  chi2[i]->Fit(fitchi2[i],"S");
       	    // int nbin=TMath::FloorNint(10000/(TMath::Sqrt(chi2[i]->GetEntries())));
       	    // chi2[i]->Rebin(nbin);
      	       chi2reco[i][j][k]->SetTitle(Form("#chi^{2} distribution, %s %s %s %d%s",stype.c_str(),RecoInd[i].c_str(), TConfig[k].c_str(), Ene[j],UNIT.c_str()));
      	       chi2reco[i][j][k]->Draw("hist");
               fileout[i]->cd();
               chi2reco[i][j][k]->Write();
	   //    fitchi2[i]->Draw("same")
        } //end k


    canchi2[i][j]->SaveAs(Form("%s/%d/Chi2_Simplified_%s_%d%s.eps", directory.c_str(),type,RecoInd[i].c_str(),Ene[j],UNIT.c_str()));
    if (i==0&&j==0) canchi2[i][j]->Print(Form("%s/%d/Chi2_Simplified.pdf(", directory.c_str(),type));
    else if (i==NReco-1&&j==Nene-1) canchi2[i][j]->Print(Form("%s/%d/Chi2_Simplified.pdf)",directory.c_str(),type));
    else canchi2[i][j]->Print(Form("%s/%d/Chi2_Simplified.pdf", directory.c_str(),type));
  } //end j
} //end i

//////////////////////////////////////////////////////
///////// Momentum Pull Distributions ////////////////
//////////////////////////////////////////////////////

cout << "momentum pull histograms" << endl;

TCanvas**CanMomPulls=new TCanvas*[NReco];
gStyle->SetOptStat("nemruo");
for(int i=0;i<NReco;i++) {
	CanMomPulls[i]=new TCanvas(Form("%s momentum resolution", RecoInd[i].c_str()),Form("%s momentum resolution", RecoInd[i].c_str()),200,10,1500,1000);
	CanMomPulls[i]->Divide(2,2);
	for(int j=0;j<Nene;j++) {
		CanMomPulls[i]->cd(j+1);
      	        int nbin = TMath::FloorNint(10000/(TMath::Sqrt(ResoReco[i][j]->GetEntries())));
       	        ResoReco[i][j]->Rebin(1.5*nbin);
	        ResoReco[i][j]->GetXaxis()->SetRangeUser(ResoReco[i][j]->GetMean()-2*ResoReco[i][j]->GetStdDev(),ResoReco[i][j]->GetMean()+2*ResoReco[i][j]->GetStdDev());        
       	        ResoReco[i][j]->Draw("hist");
                fileout[i]->cd();
                ResoReco[i][j]->Write();
        } //end j
    CanMomPulls[i]->SaveAs(Form("%s/%d/MomReso_Simplified_%s.eps", directory.c_str(),type,RecoInd[i].c_str()));
    if (i==0) CanMomPulls[i]->Print(Form("%s/%d/MomReso_Simplified.pdf(", directory.c_str(),type));
    else if (i==NReco-1) CanMomPulls[i]->Print(Form("%s/%d/MomReso_Simplified.pdf)",directory.c_str(),type));
    else CanMomPulls[i]->Print(Form("%s/%d/MomReso_Simplified.pdf", directory.c_str(),type));
} //end i


} //end function


