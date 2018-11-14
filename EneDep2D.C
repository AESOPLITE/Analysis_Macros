////////////////////////////////////////////////////////////////////////////////////////// 
////   Author: Sarah Mechbal, smechbal@ucsc.edu
////   Department of Physics, University of California Santa Cruz, August 2018
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
	
void EneDep2D(int t, int nene, int ncycles);

void EneDep2D(int t, int nene, int ncycles, string Reco) 
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
 	
 
 
const int NTConfig=3;//Number of configuration of triggers
 string TConfig[NTConfig]={"T4 vs T1", "T3 vs T1", "T2 vs T1"};

//Input file 
string Inppath="~/AESOPLITE/MCProduction/Detector/ProcessedFiles";	
string Outpath ="~/AESOPLITE/Analysis/MCAnalysis/Scintillators/V4";
string startfile="aesopliteNonUniB_V4"; 
string endfile="_fort.99";
string RecoID= Reco;

//type of particle
int type=t;
string stype;	
if(t==3)	stype="e^{#minus}";
if(t==4)	stype="e^{#plus}";
if(t==11)	stype="#mu^{#minus}";
if(t==10)	stype="#mu^{#plus}";
if(t==1)	stype="protons";
if(t==6)	stype="#alpha";
//Number of energies
int Nene=nene;
//Energies
int*Ene=new int[Nene];
string UNIT="GeV";		
 if(t==3||t==4)
  {
 UNIT="MeV"; 
 Ene[0]= 20;
 Ene[1]= 40;
 Ene[2]= 100;
 Ene[3]= 300;
  }	
 if(t==11||t==10||t==6)
  {
   Ene[0]=1;
   Ene[1]=2;
   Ene[2]=3;
 //  Ene[3]=60;
//   Ene[4]=20;
}
//particle mass

float mass=0.000511;//electon mass in GeV
if(t==11 || t==10)     mass=0.10566;//muon mass in GeV
if(t==1)     mass=0.93827;//proton mass in GeV
if(t==6)     mass=3.7273;//alpha-particle mass in GeV

	
//2D Histograms of energy deposited in scintillators in MeV
 TH2F***TrigEne=new TH2F**[Nene];
 
 for(int i=0;i<Nene;i++)
   {
	 TrigEne[i] = new TH2F*[NTConfig];
	 for(int j=0; j<NTConfig;j++) 
	 	{
		TrigEne[i][j] =  new TH2F(Form("%s %d%s, Energy %s",stype.c_str(), Ene[i], UNIT.c_str(),TConfig[j].c_str()),Form("%s %d%s, Energy %s",stype.c_str(), Ene[i], UNIT.c_str(),TConfig[j].c_str()),400,0.,40,400,0.,40);
		TrigEne[i][j]->SetTitle(Form("%s",TConfig[j].c_str()));
		TrigEne[i][j]->GetXaxis()->SetTitle("Energy deposited in T1 (MeV)");
		TrigEne[i][j]->GetYaxis()->SetTitleOffset(2.);
   	} //end j 
 }	//end i
	
//Number of cycles per energy
 int* Ncycles=new int[Nene]; 
 for(int i=0;i<Nene;i++)Ncycles[i]=ncycles;		
//ROOT output file
TFile *fileout = new TFile(Form("%s/ELoss2D_%d_%s.root",Outpath.c_str(),type, RecoID.c_str()),"RECREATE");

//TChain for input files
TChain **chain = new TChain*[Nene];

for(int i=0; i<Nene; i++) 
	{
	chain[i] = new TChain("MC");
	chain[i]->Add(Form("%s/%d/RecoEvent_%s_%d_%d%s%s_%s.root", Inppath.c_str(),type,startfile.c_str(),type,Ene[i],UNIT.c_str(),endfile.c_str(),RecoID.c_str()));


	
	//Define variables to read event
	ALEvent *e = new ALEvent();      
	//Set address to access event data
	chain[i]->SetBranchAddress("Revent",&e); 
	 int nentries=chain[i]->GetEntries();
	 cout << "Number  of events: " << nentries << endl;	
	
	//loop over all events
	for (int k=0;k<nentries;k++)
		 {
		chain[i]->GetEntry(k);			//load the entry i in the variable e
		double X0 = e->get_X0MC();
		double CX0 = e->get_CX0MC();
		double Z0 = e->get_Z0MC();
		double CZ0 = e->get_CZ0MC();
		int nhits = e->get_Nhits();
		bool T1=e->get_T1();
		bool T3=e->get_T3();
		bool T4=e->get_T4();

	//cout << "event " << i << endl;
		if( !T1 || !T3 || !T4 ) continue;
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
		if(nTFoam==0) continue;
		for(int j=0;j<1;j++) {
			EFoam+=1000.*e->get_EneIsofoam().at(j);
			tFoam=e->get_timeIsofoam().at(j);
		}
		//Shell
		int nTShell=e->get_EneShell().size();
		if(nTShell==0) continue;
		for(int j=0;j<1;j++) {
			EShell+=1000.*e->get_EneShell().at(j);
			tShell=e->get_timeShell().at(j);
		}
		//T1
		int nT1=e->get_EneT1().size();
		for(int j=0;j<nT1;j++) {E1+=1000.*e->get_EneT1().at(j);}
		//T2
		int nT2=e->get_EneT2().size();
		for(int j=0;j<nT2;j++) {E2+=1000.*e->get_EneT2().at(j);}
		
		//T3
		int nT3=e->get_EneT3().size();
		for(int j=0;j<nT3;j++) {E3+=1000.*e->get_EneT3().at(j);}

		//T4
		int nT4=e->get_EneT4().size();
		for(int j=0;j<nT4;j++) {E4+=1000.*e->get_EneT4().at(j);}
		//cout << "E1 = " << E1 << "MeV, E3 = " << E3 << "MeV, E4 = " << E4 << "MeV" << endl;
		TrigEne[i][0]->Fill(E1,E4);			//T1 vs T4 plot
		TrigEne[i][1]->Fill(E1,E3);			//T1 vs T3 plot
		TrigEne[i][2]->Fill(E1,E2);			//T1 vs T3 plot


	} //end loop on all events
}     //end loop on energies
	cout << "done looping the events" << endl;
	
	
//////////////////////////////////
////////Histogram Display/////////
//////////////////////////////////
	
 fileout->cd();
 TCanvas**canTrigEne=new TCanvas*[Nene];
 for (int i=0; i< Nene; i++) 	
 {
	 
	 canTrigEne[i] = new TCanvas(Form("canTrigEne %d%s", Ene[i], UNIT.c_str()), Form("canTrigEne %d%s", Ene[i], UNIT.c_str()),200,200,800,800); 
     canTrigEne[i]->Divide(1,3);
    for(int j=0;j<NTConfig;j++) 
		{
		canTrigEne[i]->cd(j+1); 
		gPad->SetLeftMargin(0.15);
		TrigEne[i][j]->GetYaxis()->SetTitleOffset(2);  
		TrigEne[i][j]->Draw("colz");
		TrigEne[i][j]->Write();
   		}  
	 	canTrigEne[i]->Write();
	
		////////////////////
		//Save plots into a pdf file
		//////////////////  
//		if(i==0) canTrigEne[i]->Print(Form("%s/Eloss2D_%d.pdf(",Outpath.c_str(),type));	
//	 	else if(i==Nene-1) canTrigEne[i]->Print(Form("%s/Eloss2D_%d.pdf)",Outpath.c_str(),type));	
//	 	else canTrigEne[i]->Print(Form("%s/Eloss2D_%d.pdf",Outpath.c_str(),type));	
 } //end loop on energies
fileout->Close();
} // end function
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

