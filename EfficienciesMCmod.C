#include "headers.h"
#include "ALEvent.h"
#include "ALTckhit.h"
#include "LoadMCparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void EfficienciesMCmod(int t, int nene, int ncycles, string s, float gsource)
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
//const int nEfficiencies=2;
//string sEff[nEfficiencies] = {"T1\\&T3","T1\\&T3\\&L6"};
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

float Gsource = gsource;

//Number of energies
int Nene=nene;
//Energiesi
int*Ene=new int[Nene];
//Number of cycles per energy
int Ncycles=ncycles;

////TeX output file
ofstream totex;
totex.open(Form("%s/%s/EfficienciesMC_%s.txt",directory.c_str(),source.c_str(),source.c_str()));

ofstream summary;
summary.open(Form("%s/AllRadius.txt",directory.c_str()),fstream::app);



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
 Ene[8]= 300;
 Ene[9]= 400;
 Ene[10]= 500;













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

float*** nevents= new float**[NReco];
float*** nT1T4= new float**[NReco];
float*** nT1T3= new float**[NReco];
float*** nT1= new float**[NReco];
float*** nT1T3sig= new float**[NReco];
//Average value over all cycles
float*** Geo= new float**[NReco];
float*** GeoSig= new float**[NReco];

//Rate Histograms
TH1F****nEff= new TH1F***[NReco];	
//Reco efficiency histograms
TH1F****nEffReco= new TH1F***[NReco];
//Scattering histograms
TH1F ****ScatNB=new TH1F***[NReco];
TH1F ****ScatB=new TH1F***[NReco];	
//Energy Loss histograms
TH1F ****ELoss=new TH1F***[NReco];
//Zenith angle distribution histograms
TH1F ****ZenDist=new TH1F***[NReco];
//Energy Deposited in T1/T2/T3 
TH1F *****EDepPHA=new TH1F****[NReco];

//2D Histogram efficiencies as a function of zenith angle
TH2F ****EffZen=new TH2F***[NReco];


//Summary graph geo factor
TGraphErrors *GeoFac = new TGraphErrors();
TGaxis::SetMaxDigits(3);	

 for(int i=0;i<NReco;i++) 
   {
	nEff[i]=new TH1F**[Nene];
	nEffReco[i]=new TH1F**[Nene];
	ScatNB[i]=new TH1F**[Nene];
	ScatB[i]=new TH1F**[Nene];
	ELoss[i]=new TH1F**[Nene];
	ZenDist[i]=new TH1F**[Nene];
	EDepPHA[i]=new TH1F***[Nene];
	EffZen[i]=new TH2F**[Nene];
	nevents[i] = new float*[Nene];
	nT1T4[i] = new float*[Nene];
	nT1T3[i] = new float*[Nene];
        nT1[i] = new float*[Nene];
        nT1T3sig[i] = new float*[Nene];
	Geo[i] = new float*[Nene];
        GeoSig[i] = new float*[Nene];
	for(int j=0; j<Nene; j++) 
		{	
			nevents[i][j]=new float[Ncycles];
			nT1T4[i][j]=new float[nEfficiencies];
			nT1T3[i][j]=new float[nEfficiencies];
		        nT1[i][j]=new float[nEfficiencies];
			nT1T3sig[i][j]=new float[nEfficiencies];
			Geo[i][j]=new float[nEfficiencies];
		        GeoSig[i][j]=new float[nEfficiencies];
			nEff[i][j]=new TH1F*[Ncycles];
			nEffReco[i][j]=new TH1F*[Ncycles];
		        ZenDist[i][j]=new TH1F*[nEfficiencies];
			EffZen[i][j]=new TH2F*[nEfficiencies];
			EDepPHA[i][j]=new TH1F**[nEfficiencies];
                        for(int k=0;k<Ncycles;k++) 
			{
				//Rate histogram
				nevents[i][j][k]=0.0;
				nEff[i][j][k] = new TH1F(Form("%s, %d%s, cycles %d Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),k+1),Form("%s, %d%s, cycles %d Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),k+1),nEfficiencies,0,nEfficiencies);
				nEff[i][j][k]->SetStats(0);
				nEff[i][j][k]->GetXaxis()->SetAlphanumeric();
				for (int l=1;l<=nEfficiencies;l++) nEff[i][j][k]->GetXaxis()->SetBinLabel(l,sEff[l-1].c_str());
				//Reconstruction Rate histogram
				nEffReco[i][j][k] = new TH1F(Form("%s, %d%s, cycles %d Reconstruction Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),k+1),Form("%s, %d%s, cycle %d Reconstruction Rate efficiencies",RecoInd[i].c_str(),Ene[j],UNIT.c_str(),k+1),NConfig,0,NConfig);
				//nEffReco[i][j]->GetYaxis()->SetMaxDigits(3);			
		        	nEffReco[i][j][k]->SetStats(0);
				nEffReco[i][j][k]->GetXaxis()->SetAlphanumeric();
				for (int l=1;l<=NConfig;l++) nEffReco[i][j][k]->GetXaxis()->SetBinLabel(l,TConfig[l-1].c_str());
			} //end k cycles
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
				ELoss[i][j][k] =new TH1F(Form("%s %d%s, %s, ELoss",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),Form("%s %d%s, %s, ELoss",source.c_str(),Ene[j],UNIT.c_str(),LayerConfig[k].c_str()),100,0,10); 
				ELoss[i][j][k]->GetXaxis()->SetTitle("Energy (in MeV)");   
				ELoss[i][j][k]->GetYaxis()->SetTitle("#entries");
				ELoss[i][j][k]->Sumw2(true);	

 			} // end k
			for(int m=0; m<nEfficiencies; m++)
                        {
				Geo[i][j][m]=0;
				GeoSig[i][j][m]=0;
		        	ZenDist[i][j][m] =new TH1F(Form("%s %d%s, %s ZenDist",source.c_str(),Ene[j],UNIT.c_str(),sEff[m].c_str()),Form("%s %d%s, %s, ZenDist",source.c_str(),Ene[j],UNIT.c_str(),sEff[m].c_str()),200,0,90);
                       		 ZenDist[i][j][m]->GetXaxis()->SetTitle("Zenith angle at the source (in degrees)");
                       		 ZenDist[i][j][m]->GetYaxis()->SetTitle("#entries");
                       		 ZenDist[i][j][m]->Sumw2(true);
			         EffZen[i][j][m]=new TH2F(Form("%s %d%s, %s EffZen",source.c_str(),Ene[j],UNIT.c_str(),sEff[m].c_str()),Form("%s %d%s, %s, EffZen",source.c_str(),Ene[j],UNIT.c_str(),sEff[m].c_str()),300,0,90,100000,0,100000);
				 EffZen[i][j][m]->GetXaxis()->SetTitle("Zenith angle at the source (in degrees)");
                                 EffZen[i][j][m]->GetYaxis()->SetTitle("nT1T3");


                        	 EDepPHA[i][j][m]=new TH1F*[3];
			for(int l=0; l<3; l++)
				{
				EDepPHA[i][j][m][l] =new TH1F(Form("%s %d%s, %s Eloss in T%d",source.c_str(),Ene[j],UNIT.c_str(),sEff[m].c_str(),l),Form("%s %d%s, %s, Eloss in T%d",source.c_str(),Ene[j],UNIT.c_str(),sEff[m].c_str(),l),200,0,10);
                        	if(l==0) EDepPHA[i][j][m][l]->GetXaxis()->SetTitle("Energy deposited in T1 (in MeV)");
                        	if(l==1) EDepPHA[i][j][m][l]->GetXaxis()->SetTitle("Energy deposited in T2 (in MeV)");
                      	 	if(l==2) EDepPHA[i][j][m][l]->GetXaxis()->SetTitle("Energy deposited in T3 (in MeV)");
                     		EDepPHA[i][j][m][l]->GetYaxis()->SetTitle("#entries");
                        	EDepPHA[i][j][m][l]->Sumw2(true);
			} // end l
		
		} //end m
   	} // end j
   } // end i
cout <<"done creating histograms"<< endl;

//output ROOT file  
TFile **fileout = new TFile*[NReco];
TChain ****chain=new TChain***[NReco];

for(int i=0; i<NReco; i++)
        {
        fileout[i]=new TFile(Form("%s/%s/EfficienciesMC_%s_%d_%s.root", directory.c_str(),source.c_str(),startfile.c_str(), type,RecoInd[i].c_str()),"RECREATE");
        //Create TChain to merge the files for a given energy
        chain[i] =new TChain**[Nene];
        for(int j=0; j<Nene; j++)
                {
                chain[i][j]=new TChain*[Ncycles];
              //  cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT <<endl;
                int nL0=0;
		int nL1=0;
                int nL2=0;
                int nL3=0;       
		int nL4=0;
		int nL5=0;
		int nL6=0;    
                int L0=0;
                int L1=0;
                int L2=0;
                int L3=0;
                int L4=0;
                int L5=0;
                int L6=0;

	    for (int l=0;l<Ncycles;l++)
                {
                 chain[i][j][l]=new TChain("MC");
       //          chain[i][j][l]->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s%04d%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene[j],UNIT.c_str(),l+1,endfile.c_str(),RecoInd[i].c_str()));
                 chain[i][j][l]->Add(Form("%s/%d/%s/RawEvent_%s_%d_%d%s%04d%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene[j],UNIT.c_str(),l+1,endfile.c_str()));


                //Define variables to read event


   //Define variables to read event
   ALEvent *e = new ALEvent();
   //Set address to access event data
//   chain[i][j][l]->SetBranchAddress("Revent",&e);
   chain[i][j][l]->SetBranchAddress("event",&e);

   // Get number of event in Tree
   int nentries=chain[i][j][l]->GetEntries();
   for (int k=0;k<nentries;k++)
   {
        nevents[i][j][l]++;
	chain[i][j][l]->GetEntry(k); //Load the entry i in the variable e 
	if(k%10000==0)   {
//		cout << "Reco Type: " << RecoInd[i] << "  Energy: " << Ene[j] << UNIT << ", Cycle: " << l+1 <<  "  Event: " << k <<endl;
	//	cout << "nL4 = " << nL4 << ", nL5 = " << nL5 << ", nL6 = " << nL6 << endl;
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
        //T1
        int nnT1=e->get_EneT1().size();
        if(nnT1!=0) {
        	for(int j=0;j<nnT1;j++) {E1+=1000.*e->get_EneT1().at(j);}
	}
        //T2
        int nT2=e->get_EneT2().size();
        if(nT2!=0) {
                for(int j=0;j<nT2;j++) {E2+=1000.*e->get_EneT2().at(j);}
	}
        //T3
        int nT3=e->get_EneT3().size();
        if(nT3!=0){
                for(int j=0;j<nT3;j++) {E3+=1000.*e->get_EneT3().at(j);}
	}

	//Number of layers with hit(s) in bending/non-bending plane
	int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	//Number of layers with hit(s)
	int NL= e->get_NLayers();
	bool Lhit[7] = {false,false,false,false,false,false,false};
//	cout << "Event " << k << ", total number of hits = " << nhits << endl;
//	cout << " NL = " << NL << " , NLB = " << NLB << " , NLNB = " << NLNB << endl;
	//Variables at the source
	double CX0 = e->get_CX0MC();
        double CY0 = e->get_CY0MC();
        double CZ0 = e->get_CZ0MC();
	double EkMC = 1000*e->get_EkMC();
	//Calculate zenith angle in degrees
	double acosCZ = TMath::ACos(CZ0);
	double zen = 180 - (180.0/TMath::Pi())*acosCZ; 
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
	     int L=e->get_hits().at(l)->get_L();
	     int mreg = e->get_hits().at(l)->get_mregMC();
             if (mreg==0)  continue;
             int Lindex = int(e->get_hits().at(l)->get_mregMC())%11;
             if (L!=Lindex) cout << "event " << k << ", L = " << L << ", mreg = " << mreg << endl;
	     if(Lindex==0) nL0++;
             if(Lindex==1) nL1++;
             if(Lindex==2) nL2++;
	     if(Lindex==3) nL3++;	   
	     if(Lindex==4) nL4++;
	     if(Lindex==5) nL5++;
	     if(Lindex==6) nL6++;
	     Lhit[L] = true;
             if(L==0) L0++;
             if(L==1) L1++;
             if(L==2) L2++;
             if(L==3) L3++;
             if(L==4) L4++;
             if(L==5) L5++;
             if(L==6) L6++;
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
	//	ELoss[i][j][Lindex]->Fill(eMC);
			}
	} //end loop on hits

/*
         if(TShell && T1 && T3){
                nEff[i][j][l]->Fill(sEff[0].c_str(),1);
                ZenDist[i][j][0]->Fill(zen);
                EDepPHA[i][j][0][0]->Fill(E1);
                EDepPHA[i][j][0][1]->Fill(E2);
                EDepPHA[i][j][0][2]->Fill(E3);
                nT1T3[i][j][l]++;
	        EffZen[i][j][0]->Fill(zen,nT1T3[i][j][l]);

        }
         if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5] && Lhit[6]){
                nEff[i][j][l]->Fill(sEff[1].c_str(),1);
                ZenDist[i][j][1]->Fill(zen);
                double eMC = 1000*(e->get_hits().at(6)->get_eMC());
                ELoss[i][j][1]->Fill(EkMC-eMC);
                EDepPHA[i][j][1][0]->Fill(E1);
                EDepPHA[i][j][1][1]->Fill(E2);
                EDepPHA[i][j][1][2]->Fill(E3);
                }
*/



	 //Fill counter histograms
	 if(TShell){
		nEff[i][j][l]->Fill(sEff[0].c_str(),1);
		ZenDist[i][j][0]->Fill(zen);
		EDepPHA[i][j][0][0]->Fill(E1);
                EDepPHA[i][j][0][1]->Fill(E2);
                EDepPHA[i][j][0][2]->Fill(E3);

		}
	 if(TShell && T1){
		nEff[i][j][l]->Fill(sEff[1].c_str(),1);
                ZenDist[i][j][1]->Fill(zen);
                EDepPHA[i][j][1][0]->Fill(E1);
                EDepPHA[i][j][1][1]->Fill(E2);
                EDepPHA[i][j][1][2]->Fill(E3);		
                nT1[i][j][l]++;
	} 
	 if(TShell && T1 && T3){
		nEff[i][j][l]->Fill(sEff[2].c_str(),1);
		ZenDist[i][j][2]->Fill(zen);
                EDepPHA[i][j][2][0]->Fill(E1);
                EDepPHA[i][j][2][1]->Fill(E2);
                EDepPHA[i][j][2][2]->Fill(E3);
                nT1T3[i][j][l]++;
	
	} 
	 if(TShell && T1 && T3 && Lhit[0]){
		 nEff[i][j][l]->Fill(sEff[3].c_str(),1);
                 ZenDist[i][j][3]->Fill(zen);               
		 double eMC = 1000*(e->get_hits().at(0)->get_eMC());      
                 ELoss[i][j][0]->Fill(EkMC-eMC);
                EDepPHA[i][j][3][0]->Fill(E1);
                EDepPHA[i][j][3][1]->Fill(E2);
                EDepPHA[i][j][3][2]->Fill(E3);

		}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1]){
		nEff[i][j][l]->Fill(sEff[4].c_str(),1);
                ZenDist[i][j][4]->Fill(zen);	
		double eMC = 1000*(e->get_hits().at(1)->get_eMC());
                ELoss[i][j][1]->Fill(EkMC-eMC);
                EDepPHA[i][j][4][0]->Fill(E1);
                EDepPHA[i][j][4][1]->Fill(E2);
                EDepPHA[i][j][4][2]->Fill(E3);

		}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2]){
                ZenDist[i][j][5]->Fill(zen);
		nEff[i][j][l]->Fill(sEff[5].c_str(),1);
	        double eMC = 1000*(e->get_hits().at(2)->get_eMC());
                ELoss[i][j][2]->Fill(EkMC-eMC);
                EDepPHA[i][j][5][0]->Fill(E1);
                EDepPHA[i][j][5][1]->Fill(E2);
                EDepPHA[i][j][5][2]->Fill(E3);

		}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3]){
                ZenDist[i][j][6]->Fill(zen);
		nEff[i][j][l]->Fill(sEff[6].c_str(),1);
                double eMC = 1000*(e->get_hits().at(3)->get_eMC());
                ELoss[i][j][3]->Fill(EkMC-eMC);
                EDepPHA[i][j][6][0]->Fill(E1);
                EDepPHA[i][j][6][1]->Fill(E2);
                EDepPHA[i][j][6][2]->Fill(E3);
		 }
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4]){
		nEff[i][j][l]->Fill(sEff[7].c_str(),1);
                ZenDist[i][j][7]->Fill(zen);
		double eMC = 1000*(e->get_hits().at(4)->get_eMC());
                ELoss[i][j][4]->Fill(EkMC-eMC);
                EDepPHA[i][j][7][0]->Fill(E1);
                EDepPHA[i][j][7][1]->Fill(E2);
                EDepPHA[i][j][7][2]->Fill(E3);
		}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5]){
		nEff[i][j][l]->Fill(sEff[8].c_str(),1);
                ZenDist[i][j][8]->Fill(zen);
                double eMC = 1000*(e->get_hits().at(5)->get_eMC());
                ELoss[i][j][5]->Fill(EkMC-eMC);
                EDepPHA[i][j][8][0]->Fill(E1);
                EDepPHA[i][j][8][1]->Fill(E2);
                EDepPHA[i][j][8][2]->Fill(E3);
		}
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5] && Lhit[6]){
		nEff[i][j][l]->Fill(sEff[9].c_str(),1);
                ZenDist[i][j][9]->Fill(zen);
                double eMC = 1000*(e->get_hits().at(6)->get_eMC());
                ELoss[i][j][6]->Fill(EkMC-eMC);
                EDepPHA[i][j][9][0]->Fill(E1);
                EDepPHA[i][j][9][1]->Fill(E2);
                EDepPHA[i][j][9][2]->Fill(E3);
		}	
	 if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && Lhit[3] && Lhit[4] && Lhit[5] && Lhit[6] && T4){
		nEff[i][j][l]->Fill(sEff[10].c_str(),1);
                ZenDist[i][j][10]->Fill(zen);
                EDepPHA[i][j][10][0]->Fill(E1);
                EDepPHA[i][j][10][1]->Fill(E2);
                EDepPHA[i][j][10][2]->Fill(E3);

		}



   if(T1 && T3 )  nEffReco[i][j][l]->Fill(TConfig[0].c_str(),1);
	 if(T1 && T3 && NL>=5){nEffReco[i][j][l]->Fill(TConfig[1].c_str(),1);} 
	 if(T1 && T3 && NL>=6){nEffReco[i][j][l]->Fill(TConfig[2].c_str(),1);}
	 if(T1 && T3 && NL>=7){nEffReco[i][j][l]->Fill(TConfig[3].c_str(),1);}

	    } // end k, number of entries 


        // cout << "Ene = " << Ene[j] << "MeV, nTotal = " << chain[i][j][l]->GetEntries() <<", nT1 = " <<  nT1[i][j][l] << endl;
         //cout << "Ene = " << Ene[j] << "MeV, nT1T3 = " <<  nT1T3[i][j][l] << endl;
	
         //summary << "Ene = " << Ene[j] << "MeV, nT1T3 = " <<  nEff[i][j][l]->GetBinContent(3) << " / = " <<  nEff[i][j][l]->GetBinError(3) << endl;
	 //cout << " T1T3 bin content = " << nEff[i][j][l]->GetBinContent(3) << " with error = " <<  nEff[i][j][l]->GetBinError(3) << endl;
//       cout << "source " << source.c_str() << ", Ene = " << Ene[j] << "MeV, cycle = " << l+1 <<",  T1T3 bin content = " << nEff[i][j][l]->GetBinContent(3) << " +/- " <<  nEff[i][j][l]->GetBinError(3) << endl;
  //       cout << "source " << source.c_str() << ", Ene = " << Ene[j] << "MeV, cycle = " << l+1 <<",  tShell bin content = " << nEff[i][j][l]->GetBinContent(1) << " +/- " <<  nEff[i][j][l]->GetBinError(1) << endl;

        chain[i][j][l]->Delete();
	
	  } //end l number of cycles
	} // end j, energies reconstruction
    } // end i, reco type


//////////////////////////////
// Display histograms////////
/////////////////////////////   

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);	
gROOT->SetStyle("Modern");

//////////////////////////////////////////////////////
//////Calculate Geo Factor average of all cycles//////
//////////////////////////////////////////////////////

TH1F ****GeoHist = new TH1F***[NReco];
for(int i=0;i<NReco;i++) {
	GeoHist[i]=new TH1F**[Nene];
	for(int j=0; j<Nene;j++) {
		GeoHist[i][j]=new TH1F*[nEfficiencies];
		for(int l=0;l<nEfficiencies;l++) {
			GeoHist[i][j][l] = new TH1F(Form("%s %d%s, %s GeoHist",source.c_str(),Ene[j],UNIT.c_str(),sEff[l].c_str()),Form("%s %d%s, %s, GeoHist",source.c_str(),Ene[j],UNIT.c_str(),sEff[l].c_str()),20000,0,1000);
		}
	}
}

cout << "checking initialization" << endl;
cout << "nEff = " << nEfficiencies << ", Nene = " << Nene << ", NReco = " << NReco << endl;
for(int l=0;l<nEfficiencies;l++) {
        for (int j=0; j<Nene;j++) {
          for(int i=0;i<NReco;i++) {
	       Geo[i][j][l]=0.0;
	       GeoSig[i][j][l]=0.0;	
	//       cout << Form(" initial %s %d%s geo = %2.5f", sEff[l].c_str(), Ene[j],UNIT.c_str(), Geo[i][j][l]) << endl;
		}
	}
}

cout << "Calculate Geo Factor average of all cycles" << endl;
for(int l=0;l<nEfficiencies;l++) {
  //  	totex << Form("%s& ", sEff[l].c_str());
	for (int j=0; j<Nene;j++) {
          for(int i=0;i<NReco;i++) {
	    for (int k=0;k<Ncycles;k++) {
                cout <<  Form("%d%s cycle %d %s %2.5f",Ene[j],UNIT.c_str(),k+1,sEff[l].c_str(),Gsource*nEff[i][j][k]->GetBinContent(l+1)/nevents[i][j][k]) << " +/- " << Form("%2.5f",Gsource*nEff[i][j][k]->GetBinError(l+1)/nevents[i][j][k]) << endl;
		 GeoHist[i][j][l]->Fill(Gsource*(nEff[i][j][k]->GetBinContent(l+1))/(nevents[i][j][k]));
		 } //end k cycles

        //cout  << Form("%2.2f",Geo[i][j][l]);
        //if(j!=Nene-1){totex << " & " ;} 
		} //end i
	       } //end j
//	totex << " \\\\" << endl;
 	 } //end l	
	

cout << "Calculate Geo Factor variance and std dev " << endl;
//////////////////////////////////////////////////////
//////Calculate Geo Factor variance and std dev///////
//////////////////////////////////////////////////////
	
        for(int l=0;l<nEfficiencies;l++) {
            totex << Form("%s& ", sEff[l].c_str());   	    
	  for (int j=0; j<Nene;j++) {
		for(int i=0;i<NReco;i++) {
//		for (int k=0;k<Ncycles;k++) {
				//variance
	//	GeoSig[i][j][l]+=(float)TMath::Power(Gsource*(nEff[i][j][k]->GetBinContent(l+1)/nevents[i][j][k]) - Geo[i][j][l],2);     			
//    	 } //end k cycles
			//std dev is sqrt of variance
 	       // GeoSig[i][j][l] = (float)TMath::Sqrt(GeoSig[i][j][l]/(Ncycles-1));
          //      GeoSig[i][j][l] = (float)TMath::Sqrt(GeoSig[i][j][l]/(Ncycles));		
	//	nT1T3sig[i][j][l] = (float)TMath::Sqrt(nT1T3sig[i][j][l]/(Ncycles-1));
        cout << Form("from hist %s %d%s sig geo = %2.5f", sEff[l].c_str(), Ene[j],UNIT.c_str(), GeoHist[i][j][l]->GetStdDev()) << endl;

			//standard error of the mean SEM = SD /sqrt(n)i
		float GeoMean = GeoHist[i][j][l]->GetMean();
		float GeoSEM = GeoHist[i][j][l]->GetStdDev()/(TMath::Sqrt(Ncycles));
                 totex << Form("%2.2f $\\pm$ %2.2f",GeoMean,GeoSEM);
        if(j!=Nene-1){totex << " & " ;}
              //     cout << Form("%s %d%s, %2.2f +/- %2.2f",sEff[l].c_str(),Ene[j],UNIT.c_str(),Geo[i][j][l],GeoSig[i][j][l]) << endl;
                   cout << Form("from hist %s %d%s, %2.2f +/- %2.2f",sEff[l].c_str(),Ene[j],UNIT.c_str(),GeoMean,GeoSEM) << endl;

	if(l==nEfficiencies-2) {
                   GeoFac->SetPoint(j+1,Ene[j],GeoMean);
                   GeoFac->SetPointError(j+1,0.0,GeoSEM);
		} // end if
	    } //end i
	} //end j
    totex << " \\\\" << endl;
  } //end l


/////////////////////////////////
/////////GeoFactor plot/////////
///////////////////////////////
cout << "GeoFactor plot " << endl;
 for(int i=0;i<NReco;i++)
   {
	TCanvas *canGeo=new TCanvas(Form("Geo Factor %s",source.c_str()),Form("Geo Factor %s",source.c_str()),200,10,800,1500);
	canGeo->cd(1);
	GeoFac->GetYaxis()->SetTitleOffset(1.4);
	GeoFac->GetYaxis()->SetTitle("Geometry factor (cm^{2} sr)");
	GeoFac->GetXaxis()->SetTitle(Form("Energy (%s)", UNIT.c_str()));
	GeoFac->SetTitle(Form("Geometry Factor %s",source.c_str()));
	GeoFac->GetYaxis()->SetLimits(0,25);
	GeoFac->GetXaxis()->SetLimits(0,550);
	GeoFac->SetMarkerStyle(20);
	GeoFac->SetMarkerSize(1);	 
	GeoFac->Draw("aple");
	cout << "here?" <<endl;
	fileout[i]->cd();
	gPad->Update();
	canGeo->Write();
        GeoFac->Write();
	canGeo->SaveAs(Form("%s/Figures/GeoFac_%d_%s.eps", directory.c_str(),type,source.c_str()));
	cout << "or here? " << endl;
 }


/////////////////////////////////////
///////Scattering distribution//////
///////////////////////////////////

cout << "scattering distribution " << endl;
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
////Zenith angle distribution///////
///////////////////////////////////
cout << "zenith angle " << endl;
 gStyle->SetOptStat(11111111);
 TCanvas****canzen=new TCanvas***[NReco];
 for(int i =0; i<NReco;i++) {
         fileout[i]->cd();
         canzen[i] = new TCanvas**[Nene];
        for(int j=0;j<Nene;j++) {
            canzen[i][j]=new TCanvas*[nEfficiencies];
            for(int k=0; k<nEfficiencies;k++) {
              canzen[i][j][k]=new TCanvas(Form("%s %d%s,zenith angle",sEff[k].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s, zenith angle",sEff[k].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
              canzen[i][j][k]->SetTitle(Form("%s initial zenith angle",sEff[k].c_str()));
              canzen[i][j][k]->cd(1);
              ZenDist[i][j][k]->Draw("hist");
	      EffZen[i][j][k]->Write(); 
              canzen[i][j][k]->Write();
                 } //end k
               } //end j
           } //end i


/////////////////////////////////////////////
////Energy Deposited in scintillators///////
///////////////////////////////////////////

 cout << " Energy Deposited in scintillators" << endl;
 gStyle->SetOptStat(11111111);
 TCanvas****candep=new TCanvas***[NReco];
 for(int i =0; i<NReco;i++) {
         fileout[i]->cd();
         candep[i] = new TCanvas**[Nene];
        for(int j=0;j<Nene;j++) {
            candep[i][j]=new TCanvas*[nEfficiencies];
            for(int k=0; k<nEfficiencies;k++) {
              candep[i][j][k]=new TCanvas(Form("%s %d%s,edep PHA",sEff[k].c_str(),Ene[j],UNIT.c_str()),Form("%s %d%s, edep PHA",sEff[k].c_str(),Ene[j],UNIT.c_str()),200,10,1200,800);
              candep[i][j][k]->SetTitle(Form("%s Energy deposited in scintillators",sEff[k].c_str()));
	      candep[i][j][k]->Divide(1,3);
              candep[i][j][k]->cd(1);
              EDepPHA[i][j][k][0]->Draw("hist");
              candep[i][j][k]->cd(2);
              EDepPHA[i][j][k][1]->Draw("hist");             \
	      candep[i][j][k]->cd(3);
              EDepPHA[i][j][k][2]->Draw("hist");
              candep[i][j][k]->Write();
                 } //end k
               } //end j
           } //end i






/////////////////////////////////////
///////ELoss distribution///////////
///////////////////////////////////
 cout << "ELoss distribution" << endl;
 gStyle->SetOptStat(11111111);
 TCanvas***canLoss=new TCanvas**[NReco];
 for(int i =0; i<NReco;i++) { 
 	 fileout[i]->cd();
	 canLoss[i] = new TCanvas*[Nene];
        for(int j=0;j<Nene;j++) {
 		     canLoss[i][j]=new TCanvas(Form("%d%s, ELoss",Ene[j],UNIT.c_str()),Form("%d%s, ELoss",Ene[j],UNIT.c_str()),200,10,1200,800);
		     canLoss[i][j]->Divide(2,4);
          //      cout << " E3 =  " << E3 << " MeV" << endl;
	         for(int k=0; k<nLayers;k++) {	 
			  ELoss[i][j][k]->SetTitle(Form("Energy loss at %s",LayerConfig[k].c_str()));
              	          canLoss[i][j]->cd(k+1);
      	      		  ELoss[i][j][k]->Draw("hist");			 
       		 } //end K
	      canLoss[i][j]->Write();
  	} //end j
	cout << "1" << endl;
	fileout[i]->Close();
	cout << "2" << endl;	  
    } //end i
	totex.close();
	cout << "3" <<endl;
	summary.close();
	cout << "4" << endl;
 
 } //end function








/*

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
        for(int k=0; k<Ncycles;k++)
          {
//         for(int l=0;l<nEfficiencies;l++) cout <<  Form("%s %d%s %d %4.0f",sEff[l].c_str(),Ene[j],UNIT.c_str(),k+1,(nEff[i][j][k]->GetBinContent(l+1))) <<endl;
           gPad->SetGridy(1);
           nEff[i][j][0]->SetTitle(Form("%s, Trigger Rates",source.c_str()));
           nEff[i][j][k]->SetLineWidth(0);
           nEff[i][j][k]->SetMarkerStyle(20);
           nEff[i][j][k]->SetFillColor(mcolor[j]);
           nEff[i][j][k]->SetMarkerSize(0);
           nEff[i][j][k]->Scale(100./nevents[i][j][k]);
           nEff[i][j][k]->SetBarWidth(barwidth);
           nEff[i][j][k]->SetBarOffset(baroffset+j*barwidth);
           nEff[i][j][k]->LabelsOption("v","X"); 
           nEff[i][j][k]->GetYaxis()->SetRangeUser(0,20);
           nEff[i][j][k]->GetYaxis()->SetNdivisions(511,0);
           nEff[i][j][k]->GetYaxis()->SetTitle("Rates (normalized to total # events)");
           if (j==0) { 
                nEff[i][j][k]->Draw("bar");
                nEff[i][j][k]->SetTitle(Form("%s",RecoInd[i].c_str()));
                }
         else {
                nEff[i][j][k]->Draw("bar same");
                }
       LegEff[i]->AddEntry(nEff[i][j][k],Form("%d%s",Ene[j],UNIT.c_str()),"f");
        //   txt[i][j]=new TText*[nEfficiencies];
          for(int l=0;l<nEfficiencies;l++)               
        {       
         //  if(l>1){txt[i][j][l]=new TText((nEff[i][j]->GetXaxis()->GetBinCenter(l+1)+(j-1)*barwidth)-0.15,nEff[i][j]->GetBinContent(l+1)+1,Form("%2.2f cm2 sr +/- %2.2f cm2 sr", Gsource*(nEff[i][j]->GetBinContent(l+1)/100), Gsource*(nEff[i][j]->GetBinError(l+1)/100))); 
     //  txt[i][j][l]->SetTextSize(0.010);
      // txt[i][j][l]->SetTextAngle(55);    
//       txt[i][j`][l]->Draw();  
//                        }
        } //end l
      }  // end k
   } //end      j
         LegEff[i]->Draw();
     fileout[i]->cd();
     gPad->Update();
     CanEff[i]->Write();
    if (i==0) CanEff[i]->Print(Form("%s/%s/Efficiencies_%d_%s_%s.pdf(", directory.c_str(),source.c_str(),type,source.c_str(),RecoInd[i].c_str()));
    else CanEff[i]->Print(Form("%s/%s/Efficiencies_%d_%s_%s.pdf", directory.c_str(),source.c_str(),type,source.c_str(),RecoInd[i].c_str()));
        } //end i
*/





