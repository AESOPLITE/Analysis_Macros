#include "headers.h"
#include "ALEvent.h"
#include "LoadMCparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void Synthesis(int t, int nene)
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
const int NConfig=11;
string TConfig[NConfig]={"T1&T4&reco", "T1&T4&T3&reco","T1&T4&T3&reco chi>0" ,"T1&T2&T3&reco & 0+ sites used","T1&T4&T3&reco & 1+ sites","T1&T4&T3&reco & 2+ sites","T1&T4&T3&reco & 3+ sites","T1&T4&T3&reco & 4+sites","T1&T4&T3&reco & 5+ sites","T1&T4&T3&reco & 6+ sites","T1&T4&T3&reco & 7 sites"};					   
const int nEfficiencies=13; 	 		//count and record number of entries for all events, T1&T4 no hits, T1&T4 1 layer with hit, T1&T4 2 layer with hits, etc
string sEff[nEfficiencies] = {"T1&T4","T1&T4&hits","T1&T4 & 1+ layers","T1&T4 & 2+ layers", " T1&T4 & 3+ layers", "T1&T4 & 4+ layers", "T1&T4 & 5+ layers ", "T1T4 & 5+ layers & track","T1&T4 & 6+ layers","T1&T4 & 6+ layers & track", "T1&T4 & 7+ layers", "T1&T4 7 layers & track","T1&T4&reco"};	
//Input file 
string Inppath="/home/sarah/AESOPLITE/MCProduction/Detector/ProcessedFiles";	
string startfile="aesopliteNonUniB_V4";
string endfile="_fort.99";
string directory= "/home/sarah/AESOPLITE/Analysis/MCAnalysis/Efficiencies/V4/20cmSource";

//type of particle
int type=t;
string stype;	
if(t==3)	stype="e^{#minus}";
if(t==4)	stype="e^{#plus}";
if(t==11)	stype="#mu^{#minus}";
if(t==10)	stype="#mu^{#plus}";

float mass=0.000511;//electon mass in GeV
if(t==11 || t==10)      mass=0.10566;//muon mass in GeV
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
 Ene[1]= 40;
 Ene[2]= 100;
 Ene[3]= 300;
 
 //Ene[0]=40;
 //Ene[1]=100;
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
		        ResoReco[i][j]=new TH1F(Form("Resolution momentum from reconstruction, %d%s",Ene[j],UNIT.c_str()),Form("Resolution momentum from reconstruction, %d%s",Ene[j],UNIT.c_str()),10000,-2,2);
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
	fileout[i]=new TFile(Form("%s/Synthesis_%s_%d_%s.root", directory.c_str(),startfile.c_str(), type,RecoInd[i].c_str()),"RECREATE");

	//Create TChain to merge the files for a given energy
	chain[i] =new TChain*[Nene]; 
  	for(int j=0; j<Nene; j++) 
  		{	  	
   		cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT <<endl;
   		chain[i][j]=new TChain("MC");
        chain[i][j]->Add(Form("%s/%d/RecoEvent_%s_%d_%d%s%s_%s.root", Inppath.c_str(),type,startfile.c_str(),type,Ene[j],UNIT.c_str(),endfile.c_str(),RecoInd[i].c_str()));
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
     //get MC informations at L0
        for(int l=0;l<e->get_Nhits();l++) 
	  {
             bool used = e->get_hits().at(l)->get_fUsed();
             if(used) sitesUsed++;
	     int Lindex=e->get_hits().at(l)->get_L();
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
         if(T1 && T4){
		nEff[i][j]->Fill(sEff[0].c_str(),1);
         	nevents[i][j]++;
		   }
        if(T1 && T4 && nhits>0){nEff[i][j]->Fill(sEff[1].c_str(),1);} 
	if(T1 && T4 && nhits>0 && NL>=1){nEff[i][j]->Fill(sEff[2].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=2){nEff[i][j]->Fill(sEff[3].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=3){nEff[i][j]->Fill(sEff[4].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=4){nEff[i][j]->Fill(sEff[5].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=5){nEff[i][j]->Fill(sEff[6].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=5 && chiNBPR>=0 && chiBPR>=0){nEff[i][j]->Fill(sEff[7].c_str(),1);}	
         if(T1 && T4 && nhits>0 && NL>=6){nEff[i][j]->Fill(sEff[8].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=6 && chiNBPR>=0 && chiBPR>=0){nEff[i][j]->Fill(sEff[9].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL==7){nEff[i][j]->Fill(sEff[10].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL==7 && chiNBPR>=0 && chiBPR>=0){nEff[i][j]->Fill(sEff[11].c_str(),1);}
         if(T1 && T4 && nhits>0 && preco>0){nEff[i][j]->Fill(sEff[12].c_str(),1);}

//Fill histograms

	if(T1 && T4 && preco>0) 
	  {
	  nEffReco[i][j]->Fill(TConfig[0].c_str(),1);
	  p0PR[i][j][0]->Fill(1/pPR);
	  p0reco[i][j][0]->Fill(1/preco);
	  pMC[i][j][0]->Fill(1/p0MC);
	  chi2reco[i][j][0]->Fill(chireco);
	  chi2NBPR[i][j][0]->Fill(chiNBPR);
	  chi2BPR[i][j][0]->Fill(chiBPR);		
	  ResoReco[i][j]->Fill((cpa_reco-cpa_MC)/(sigma_cpa)); 
	 }
        if(T1 && T3 && T4 && preco>0)
          {
          nEffReco[i][j]->Fill(TConfig[1].c_str(),1);
          p0PR[i][j][1]->Fill(1/pPR);
          p0reco[i][j][1]->Fill(1/preco);
          pMC[i][j][1]->Fill(1/p0MC);
          chi2reco[i][j][1]->Fill(chireco);
          chi2NBPR[i][j][1]->Fill(chiNBPR);
          chi2BPR[i][j][1]->Fill(chiBPR);
          }
	if(T1 && T4 && T3 && preco>0 && chireco>=0)
          {
          nEffReco[i][j]->Fill(TConfig[2].c_str(),1);
          p0PR[i][j][2]->Fill(1/pPR);
          p0reco[i][j][2]->Fill(1/preco);
          pMC[i][j][2]->Fill(1/p0MC);
          chi2reco[i][j][2]->Fill(chireco);
          chi2NBPR[i][j][2]->Fill(chiNBPR);
          chi2BPR[i][j][2]->Fill(chiBPR);
          }
       
	  if(T1 && T4 && T3 && preco>0 && sitesUsed>=0 )
          {
          nEffReco[i][j]->Fill(TConfig[3].c_str(),1);
          p0PR[i][j][3]->Fill(1/pPR);
          p0reco[i][j][3]->Fill(1/preco);
          pMC[i][j][3]->Fill(1/p0MC);
          chi2reco[i][j][3]->Fill(chireco);
          chi2NBPR[i][j][3]->Fill(chiNBPR);
          chi2BPR[i][j][3]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=1 )
          {
          nEffReco[i][j]->Fill(TConfig[4].c_str(),1);
          p0PR[i][j][4]->Fill(1/pPR);
          p0reco[i][j][4]->Fill(1/preco);
          pMC[i][j][4]->Fill(1/p0MC);
          chi2reco[i][j][4]->Fill(chireco);
          chi2NBPR[i][j][4]->Fill(chiNBPR);
          chi2BPR[i][j][4]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=2 )
          {
          nEffReco[i][j]->Fill(TConfig[5].c_str(),1);
          p0PR[i][j][5]->Fill(1/pPR);
          p0reco[i][j][5]->Fill(1/preco);
          pMC[i][j][5]->Fill(1/p0MC);
          chi2reco[i][j][5]->Fill(chireco);
          chi2NBPR[i][j][5]->Fill(chiNBPR);
          chi2BPR[i][j][5]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=3 )
          {
          nEffReco[i][j]->Fill(TConfig[6].c_str(),1);
          p0PR[i][j][6]->Fill(1/pPR);
          p0reco[i][j][6]->Fill(1/preco);
          pMC[i][j][6]->Fill(1/p0MC);
          chi2reco[i][j][6]->Fill(chireco);
          chi2NBPR[i][j][6]->Fill(chiNBPR);
          chi2BPR[i][j][6]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=4 )
          {
          nEffReco[i][j]->Fill(TConfig[7].c_str(),1);
          p0PR[i][j][7]->Fill(1/pPR);
          p0reco[i][j][7]->Fill(1/preco);
          pMC[i][j][7]->Fill(1/p0MC);
          chi2reco[i][j][7]->Fill(chireco);
          chi2NBPR[i][j][7]->Fill(chiNBPR);
          chi2BPR[i][j][7]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=5 )
          {
          nEffReco[i][j]->Fill(TConfig[8].c_str(),1);
          p0PR[i][j][8]->Fill(1/pPR);
          p0reco[i][j][8]->Fill(1/preco);
          pMC[i][j][8]->Fill(1/p0MC);
          chi2reco[i][j][8]->Fill(chireco);
          chi2NBPR[i][j][8]->Fill(chiNBPR);
          chi2BPR[i][j][8]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=6 )
          {
          nEffReco[i][j]->Fill(TConfig[9].c_str(),1);
          p0PR[i][j][9]->Fill(1/pPR);
          p0reco[i][j][9]->Fill(1/preco);
          pMC[i][j][9]->Fill(1/p0MC);
          chi2reco[i][j][9]->Fill(chireco);
          chi2NBPR[i][j][9]->Fill(chiNBPR);
          chi2BPR[i][j][9]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed==7 )
          {
          nEffReco[i][j]->Fill(TConfig[10].c_str(),1);
          p0PR[i][j][10]->Fill(1/pPR);
          p0reco[i][j][10]->Fill(1/preco);
          pMC[i][j][10]->Fill(1/p0MC);
          chi2reco[i][j][10]->Fill(chireco);
          chi2NBPR[i][j][10]->Fill(chiNBPR);
          chi2BPR[i][j][10]->Fill(chiBPR);
          }

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
double barwidth=0.2;
double baroffset=0.1;
 for(int i=0;i<NReco;i++)
   {
           CanEff[i] = new TCanvas(Form("%s, efficiencies", RecoInd[i].c_str()),Form("%s, efficiencies", RecoInd[i].c_str()),200,10,1200,800);
     txt[i]=new TText**[Nene];
	   LegEff[i] = new TLegend(0.1,0.9,0.2,1.0);
           CanEff[i]->cd();
           CanEff[i]->SetBottomMargin(0.32);

    for (int j=0; j<Nene;j++) 
	{
   	   gPad->SetGridy(1);
	   nEff[i][0]->SetTitle(Form("Trigger Rates"));
   	   nEff[i][j]->SetLineWidth(0);
   	   nEff[i][j]->SetMarkerStyle(20);
           nEff[i][j]->SetFillColor(mcolor[j]);
   	   nEff[i][j]->SetMarkerSize(0);
    	   nEff[i][j]->Scale(100./nevents[i][j]);
    	   nEff[i][j]->SetBarWidth(barwidth);
	   nEff[i][j]->SetBarOffset(baroffset+j*barwidth);
           nEff[i][j]->LabelsOption("v","X"); 
    	   nEff[i][j]->GetYaxis()->SetRangeUser(0,110);
    	   nEff[i][j]->GetYaxis()->SetNdivisions(511,0);
   	   nEff[i][j]->GetYaxis()->SetTitle("Normalized rates");
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
       if(nEff[i][j]->GetBinContent(l+1)>=100&&l==0) txt[i][j][l]=new TText(nEff[i][j]->GetXaxis()->GetBinCenter(l+1),nEff[i][j]->GetBinContent(l+1)+1,Form("%d",(int)nEff[i][j]->GetBinContent(l+1))); 
       else txt[i][j][l]=new TText((nEff[i][j]->GetXaxis()->GetBinCenter(l+1)+(j-1)*barwidth)-0.15,nEff[i][j]->GetBinContent(l+1)+1,Form("%.1f",nEff[i][j]->GetBinContent(l+1))); 
       txt[i][j][l]->SetTextSize(0.015);
       txt[i][j][l]->SetTextAngle(55);    
       txt[i][j][l]->Draw();    
      }  // end l
   } //end j
     fileout[i]->cd();
     gPad->Update();
     CanEff[i]->Write();
     CanEff[i]->SaveAs(Form("%s/Efficiencies_%d.eps",directory.c_str(),type));
    if (i==0) CanEff[i]->Print(Form("%s/Efficiencies.pdf(", directory.c_str()));
	} //end i

//////////////////////////
//Reconstruction Rates////
/////////////////////////

TCanvas**CanEffReco=new TCanvas*[NReco];
TLegend**LegEffReco=new TLegend*[NReco];
TText****txtReco=new TText***[NReco];
 for(int i=0;i<NReco;i++)
   {
           CanEffReco[i]= new TCanvas(Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),200,10,1200,800);
    	   txtReco[i]=new TText**[Nene];
           LegEffReco[i] = new TLegend(0.1,0.9,0.2,1.0);
           CanEffReco[i]->cd();
           CanEffReco[i]->SetBottomMargin(0.32);

    for (int j=0; j<Nene;j++)
        {
           gPad->SetGridy(1);
           nEffReco[i][0]->SetTitle(Form("%s efficiencies",RecoInd[i].c_str()));
           nEffReco[i][j]->SetLineWidth(0);
           nEffReco[i][j]->SetMarkerStyle(20);
           nEffReco[i][j]->SetFillColor(mcolor[j]);
           nEffReco[i][j]->SetMarkerSize(0);
           nEffReco[i][j]->Scale(100./p0reco[i][j][0]->GetEntries());
           nEffReco[i][j]->SetBarWidth(barwidth);
           nEffReco[i][j]->SetBarOffset(baroffset+j*barwidth);
           nEffReco[i][j]->LabelsOption("v","X");
           nEffReco[i][j]->GetYaxis()->SetRangeUser(0,110);
           nEffReco[i][j]->GetYaxis()->SetNdivisions(511,0);
           nEffReco[i][j]->GetYaxis()->SetTitle("Normalized rates");
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
       if(nEffReco[i][j]->GetBinContent(l+1)>=100&&l==0) txtReco[i][j][l]=new TText(nEffReco[i][j]->GetXaxis()->GetBinCenter(l+1),nEffReco[i][j]->GetBinContent(l+1)+1,Form("%d",(int)nEffReco[i][j]->GetBinContent(l+1)));
       else txtReco[i][j][l]=new TText((nEffReco[i][j]->GetXaxis()->GetBinCenter(l+1)+(j-1)*barwidth)-0.15,nEffReco[i][j]->GetBinContent(l+1)+1,Form("%.1f",nEffReco[i][j]->GetBinContent(l+1)));
       txtReco[i][j][l]->SetTextSize(0.015);
       txtReco[i][j][l]->SetTextAngle(55);
       txtReco[i][j][l]->Draw();
      }  // end l
   } //end j
     fileout[i]->cd();
     LegEffReco[i]->Draw();
     gPad->Update();
     CanEffReco[i]->Write();
     CanEffReco[i]->SaveAs(Form("%s/%sEfficienciesReco.eps", directory.c_str(),RecoInd[i].c_str()));
      if (i==NReco-1)CanEffReco[i]->Print(Form("%s/Efficiencies.pdf)", directory.c_str()));
     else CanEffReco[i]->Print(Form("%s/Efficiencies.pdf", directory.c_str()));
        } //end i

/*
///////////////////////////
//Inverse momentum////////
/////////////////////////

TCanvas***can=new TCanvas**[NReco];
TF1*****fgauss=new TF1****[NReco];
double nsigma = 3;
for(int i=0;i<NReco;i++) 
   {
 	can[i]=new TCanvas*[Nene];
	fgauss[i]=new TF1***[Nene];
	for(int j=0;j<Nene;j++)
	   {
		can[i][j]=new TCanvas(Form("%s %d%s Inverse momentum reconstruction", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s Inverse momentum reconstruction", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
		can[i][j]->Divide(4,3);
		fgauss[i][j]=new TF1**[NConfig];
		for(int k=0;k<NConfig;k++)
		   { 		
       		 fgauss[i][j][k] = new TF1*[2];
		 can[i][j]->cd(k+1);
		 double Max = 0;
		 int newbin = TMath::Sqrt(p0reco[i][j][k]->GetEntries());
		 double binratio = 300/newbin;
		 p0reco[i][j][k]->Rebin(1.3*binratio);
		 p0PR[i][j][k]->Rebin(1.3*binratio);
		// pMC[i][j][k]->Rebin(binratio);
      		 double tmp0reco=p0reco[i][j][k]->GetMaximum();
		 double tmpPR=p0PR[i][j][k]->GetMaximum();
		 if (tmp0reco > tmpPR) {Max = tmp0reco;}
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
		 p0PR[i][j][k]->Draw("hist sames");
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
                 TPaveStats *tpsPR = (TPaveStats*) p0PR[i][j][k]->FindObject("stats");
                 tpsPR->SetName("PR stats");
                 tpsPR->SetY2NDC(y1);
                 tpsPR->SetY1NDC(y1-(y2-y1));
                 tpsPR->SetTextColor(kRed);
                 tpsPR->SetLineColor(kRed);
		 double y1PR = tpsPR->GetY1NDC();
		 double y2PR = tpsPR->GetY2NDC();
                 TPaveStats *tpsMC = (TPaveStats*) pMC[i][j][k]->FindObject("stats");
                 tpsMC->SetName("MC truth stats");
                 tpsMC->SetY2NDC(y1PR);
                 tpsMC->SetY1NDC(y1PR-(y2PR-y1PR));
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

	fileout[i]->cd();
	p0reco[i][j][k]->Write();
	p0PR[i][j][k]->Write();
	pMC[i][j][k]->Write();
    } //end k	   
     can[i][j]->Write();
     can[i][j]->SaveAs(Form("%s/%d/KF/Reconstructions_%s_%d%s.eps", directory.c_str(),type,RecoInd[i].c_str(),Ene[j],UNIT.c_str()));
    if (i==0&&j==0) can[i][j]->Print(Form("%s/%d/KF/Reconstructions.pdf(", directory.c_str(),type));
    else if (i==NReco-1&&j==Nene-1) can[i][j]->Print(Form("%s/%d/KF/Reconstructions.pdf)",directory.c_str(),type));
    else can[i][j]->Print(Form("%s/%d/KF/Reconstructions.pdf", directory.c_str(),type));
  } //end j
} //end i

/////////////////////////////////////
///chi2 distribution canvas/////////
////////////////////////////////////

 gStyle->SetOptStat(11111111);
 TCanvas***canchi2=new TCanvas**[NReco];
 for(int i =0; i<NReco;i++) { 
 	canchi2[i] = new TCanvas*[Nene];
        for(int j=0;j<Nene;j++) {
            canchi2[i][j]=new TCanvas(Form("%s %d%s Chi2 distribution", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s Chi2 distribution", RecoInd[i].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
            canchi2[i][j]->Divide(4,3);
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


    canchi2[i][j]->SaveAs(Form("%s/%d/KF/Chi2_%s_%d%s.eps", directory.c_str(),type,RecoInd[i].c_str(),Ene[j],UNIT.c_str()));
    if (i==0&&j==0) canchi2[i][j]->Print(Form("%s/%d/KF/Chi2.pdf(", directory.c_str(),type));
    else if (i==NReco-1&&j==Nene-1) canchi2[i][j]->Print(Form("%s/%d/KF/Chi2.pdf)",directory.c_str(),type));
    else canchi2[i][j]->Print(Form("%s/%d/KF/Chi2.pdf", directory.c_str(),type));
  } //end j
} //end i

//////////////////////////////////////////////////////
///////// Momentum Pull Distributions ////////////////
//////////////////////////////////////////////////////

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
    CanMomPulls[i]->SaveAs(Form("%s/%d/KF/MomReso_%s.eps", directory.c_str(),type,RecoInd[i].c_str()));
    if (i==0) CanMomPulls[i]->Print(Form("%s/%d/KF/MomReso.pdf(", directory.c_str(),type));
    else if (i==NReco-1) CanMomPulls[i]->Print(Form("%s/%d/KF/MomReso.pdf)",directory.c_str(),type));
    else CanMomPulls[i]->Print(Form("%s/%d/KF/MomReso.pdf", directory.c_str(),type));
} //end i


		
////////////////////////////////////////////////////////
//////	Momentum reconstruction ////////////////////////
////////////////////////////////////////////////////////


 TMultiGraph* GraphRecoall= new TMultiGraph();
 TMultiGraph* resorecoall= new TMultiGraph();
 GraphRecoall->SetTitle("Inverse momentum recontruction vs truth momentum");
 resorecoall->SetTitle("Momentum recontruction resolution  vs truth momentum");

 TGraphErrors** GraphReco = new TGraphErrors*[NReco];
 TGraphErrors** GraphReso = new TGraphErrors*[NReco];
 for(i=0;i<NReco;i++) { 	
	TGraphErrors* p0reco = new TGraphErrors();
 	p0reco->SetLineColor(mcolor[i]);
 	p0reco->SetMarkerStyle(4);
 	p0reco->SetMarkerColor(mcolor[i]);

 	TGraphErrors* GraphReso[i]= new TGraphErrors();
    GraphReso[i]->SetLineColor(mcolor[i]);
 	GraphReso[i]->SetMarkerStyle(4);
    GraphReso[i]->SetMarkerColor(mcolor[i]);

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
	double tmp0reco=MomReco[i]->GetMaximum();
	double tmpMC= MomMC[i]->GetMaximum();
	if(tmp> tmp0reco)  MaxYMom=tmp; 
//        MaxYMom=tmpMC;
	else MaxYMom=tmp0reco;
	//MomLine[i]->SetLineWidth(2);
	MomPR[i]->SetTitle(Form("Inverse momentum PR & KF, %s %d%s",stype.c_str(),Ene[j],UNIT.c_str()));  
	MomPR[i]->Draw("hist");
	MomPR[i]->GetXaxis()->SetRangeUser(MomReco[i]->GetMean()-sigma*MomReco[i]->GetRMS(),MomReco[i]->GetMean()+sigma*MomReco[i]->GetRMS());
	MomPR[i]->SetLineColor(mcolor[1]);		   
	MomReco[i]->SetLineColor(mcolor[2]);
//        MomPR[i]->GetYaxis()->SetRangeUser(0.,1.10*MaxYMom);
      
        MomLine[i]=new TLine(1./Ene[j],0.,1/Ene[j],1.10*MaxYMom);
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
        TPaveStats *tp0reco = (TPaveStats*) MomReco[i]->FindObject("stats");
        tp0reco->SetName("KF Stats");
	tp0reco->SetY2NDC(Y1);
	tp0reco->SetY1NDC(Y1-(Y2-Y1));
	tp0reco->SetTextColor(kBlue);
	tp0reco->SetLineColor(kBlue);
	gPad->Update();
	 
//	TPaveStats *tps2 = (TPaveStats*) MomReco[i]->FindObject("stats");	
	
	
	
	  cout<<"did it crash"<<endl;
	//leg3[i]->AddEntry(MomMC[i],Form("MC momentum %s, %d%s",stype.c_str(),Ene[j],UNIT.c_str()),"l");
	
	for (int j=0;j<2;j++) {
		if (j==0) {
	float xMaxPR=MomPR[i]->GetXaxis()->GetBinCenter(MomPR[i]->GetMaximumBin());
	fMom[i][j] = new TF1(Form("fMomPR %d %d",i,j),"gaus",xMaxPR-fabs(0.2*xMaxPR),xMaxPR+fabs(0.2*xMaxPR));
	MomPR[i]->Fit(fMom[i][j],"R0");
	p0reco->SetPoint(i,1./Ene[j],fMom[i][j]->GetParameter(1));
	p0reco->SetPointError(i,0,fMom[i][j]->GetParameter(2));	
	PRreso->SetPoint(i,Ene[j],100*fabs((fMom[i][j]->GetParameter(2)/(fMom[i][j]->GetParameter(1)))));
		}
		if (j==1) {
	float xMaxReco=MomReco[i]->GetXaxis()->GetBinCenter(MomReco[i]->GetMaximumBin());
	fMom[i][j] = new TF1(Form("fMomReco %d %d",i,j),"gaus",xMaxReco-fabs(0.2*xMaxReco),xMaxReco+fabs(0.2*xMaxReco));
	MomReco[i]->Fit(fMom[i][j],"R0");
	KFreco->SetPoint(i,1./Ene[j],fMom[i][j]->GetParameter(1));
	KFreco->SetPointError(i,0,fMom[i][j]->GetParameter(2));	
	KFreso->SetPoint(i,Ene[j],100*fabs((fMom[i][j]->GetParameter(2)/(fMom[i][j]->GetParameter(1)))));
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

	 fileout->cd();
    MomMC[i]->Write();
    MomPR[i]->Write();
    MomReco[i]->Write();	
   }	   
 
 cout << " print third canvas" << endl;
 can3->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf", directory.c_str(),type, startfile.c_str(), type,ID.c_str()));
 can3->Write();
 
 
    //////////////////////////////////////////////////////////////////
TCanvas*canreso=new TCanvas("Resolution PR&KF","Resolution PR&KF",200,10,1500,1000);
canreso->Divide(2,1);
p0reco->Add(p0reco);
p0reco->Add(KFreco);
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
 canreso->Print(Form("%s/%d/PR_Analysis_%s_%d_%s.pdf)", directory.c_str(),type, startfile.c_str(), type,ID.c_str()));
 fileout->Close(); 	
	
 toTeX.close();
	
	*/

} //end function


