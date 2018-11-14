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
	

void EneDep(int t, int nene, int ncycles, string s) 
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
 	
 
 
 int NTConfig=6;//Number of configuration of triggers
 string TConfig[6]={"Isofoam","Shell", "T1","T2","T3","T4"};

//Input file 
string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4";
string Outpath ="~/Documents/AESOPLITE/Analysis/MCAnalysis/Scintillators/V4";
string startfile="aesopliteNonUniB_V4"; 
string endfile="_fort.99";
string source = s;
string RecoID= "KFone";
//Input file

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
 if(t==3||t==1)
  {
 UNIT="MeV"; 
 Ene[0]= 10;
 Ene[1]= 20;
 Ene[2]= 30;
 Ene[3]= 40;
 Ene[4]= 50;
 Ene[5]= 60; 
 Ene[6]= 100;
 Ene[7]= 300;
  }	
 if(t==11||t==10||t==6)
  {
   Ene[0]=1;
   Ene[1]=2;
   Ene[2]=3;
   Ene[3]=15;
   Ene[4]=20;
}
//particle mass

float mass=0.000511;//electon mass in GeV
if(t==11 || t==10)     mass=0.10566;//muon mass in GeV
if(t==1)     mass=0.93827;//proton mass in GeV
if(t==6)     mass=3.7273;//alpha-particle mass in GeV

	
//Histograms of energy deposited in scintillators in MeV
 TH1F***TrigEne=new TH1F**[Nene];
 
 for(int i=0;i<Nene;i++)
   {
	TrigEne[i] =  new TH1F*[NTConfig];
	for(int j=0; j<NTConfig; j++) 
		{
    TrigEne[i][j]=new TH1F(Form("%s %d%s, Energy %s",stype.c_str(), Ene[i], UNIT.c_str(),TConfig[j].c_str()),Form("%s %d%s, Energy %s",stype.c_str(), Ene[i], UNIT.c_str(),TConfig[j].c_str()),400,0.,20);
	//TrigEne[i][j]->SetTitle(Form
	TrigEne[i][j]->GetXaxis()->SetTitle("Deposited energy (MeV)");
    TrigEne[i][j]->GetYaxis()->SetTitle("Events");
    TrigEne[i][j]->GetYaxis()->SetTitleOffset(2.);
    TrigEne[i][j]->SetLineColor(kBlack);     
   	} 
 }	
	
//Number of cycles per energy
 int* Ncycles=new int[Nene]; 
 for(int i=0;i<Nene;i++)Ncycles[i]=ncycles;		
//ROOT output file
TFile *fileout = new TFile(Form("%s/ELoss_%d_%s.root",Outpath.c_str(),type, source.c_str()),"RECREATE");
//TChain for input files
TChain **chain = new TChain*[Nene];
for(int i=0; i<Nene; i++) 
	{
	chain[i] = new TChain("MC");
	for(int j=0; j<Ncycles[i]; j++)
	{
//	chain[i]->Add(Form("%s/%d/RawEvent_%s_%d_%d%s%s.root", Inppath.c_str(),type,startfile.c_str(),type,Ene[i],UNIT.c_str(),endfile.c_str()));
	chain[i]->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s0001%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene[i],UNIT.c_str(),endfile.c_str(),RecoID.c_str()));


	} //end j
	

	//Define variables to read event
	ALEvent *e = new ALEvent();      
	//Set address to access event data
//	chain[i]->SetBranchAddress("event",&e); 
	chain[i]->SetBranchAddress("Revent",&e); 

	 int nentries=chain[i]->GetEntries();
	// cout << "Number  of events: " << nentries << endl;	
	
	//loop over all events
	for (int k=0;k<nentries;k++)
		 {
		
//		cout << "event " << k << endl;
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
		//if( !T1 || !T3 || !T4 ) continue;
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
           //   cout << "nTFoam = " << nTFoam << endl;  
              if(nTFoam==0) continue;
                for(int j=0;j<nTFoam;j++) {
                        EFoam+=1000.*e->get_EneIsofoam().at(j);
                        tFoam=e->get_timeIsofoam().at(j);
			//           cout << "tFoam = " << tFoam << "  EFoam " << EFoam << " MeV" << endl;


                }
                TrigEne[i][0]->Fill(EFoam);
           //          cout << "EFoam " << EFoam << " MeV" << endl;


                //Shell
                int nTShell=e->get_EneShell().size();
            //    cout << "nTShell = " << nTShell << endl;
              if(nTShell==0) continue;
                for(int j=0;j<nTShell;j++) {
                        EShell+=1000.*e->get_EneShell().at(j);
                        tShell=e->get_timeShell().at(j);
				//	    cout << "tShell = " << tShell << "  EShell " << EShell << " MeV" << endl;

                }
		     //                cout << "EShell " << EShell << " MeV" << endl;

                TrigEne[i][1]->Fill(EShell);
                //T1
                int nT1=e->get_EneT1().size();
                if(nT1==0)continue;
                for(int j=0;j<nT1;j++) {
                        E1+=1000.*e->get_EneT1().at(j);
                        tT1=e->get_timeT1().at(j);
                }
		       //cout << " E1 =  " << E1 << " MeV" << endl;
                TrigEne[i][2]->Fill(E1);
                //T2
                int nT2=e->get_EneT2().size();
                if(nT2==0)continue;
		for(int j=0;j<nT2;j++) {E2+=1000.*e->get_EneT2().at(j);}
                TrigEne[i][3]->Fill(E2);
             //    cout << " E2 =  " << E2 << " MeV" << endl;



                //T3
                int nT3=e->get_EneT3().size();
                if(nT3==0)continue;
                for(int j=0;j<nT3;j++) {E3+=1000.*e->get_EneT3().at(j);}
                TrigEne[i][4]->Fill(E3);
          //      cout << " E3 =  " << E3 << " MeV" << endl;


                //T4
                int nT4=e->get_EneT4().size();
                if(nT4==0)continue;
                for(int j=0;j<nT4;j++) {E4+=1000.*e->get_EneT4().at(j);}
                TrigEne[i][5]->Fill(E4);
       //         cout << " E4 =  " << E4 << " MeV" << endl;



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
	 
	 canTrigEne[i] = new TCanvas(Form("canTrigEne %d%s", Ene[i], UNIT.c_str()), Form("canTrigEne %d%s", Ene[i], UNIT.c_str()),200,10,1500,500); 
     canTrigEne[i]->Divide(2,3);
    for(int j=0;j<NTConfig;j++) 
		{
		canTrigEne[i]->cd(j+1); 
		gPad->SetLeftMargin(0.15);
		TrigEne[i][j]->GetYaxis()->SetTitleOffset(2);  
		TrigEne[i][j]->Draw("hists");
		TrigEne[i][j]->Write();
   		}  
	 	canTrigEne[i]->Write();
	
		////////////////////
		//Save plots into a pdf file
		//////////////////  
//		if(i==0) canTrigEne[i]->Print(Form("%s/Eloss_%s.pdf(",Outpath.c_str(),stype.c_str()));	
//	 	else if(i==Nene-1) canTrigEne[i]->Print(Form("%s/Eloss_%s.pdf)",Outpath.c_str(),stype.c_str()));	
//	 	else canTrigEne[i]->Print(Form("%s/Eloss_%s.pdf",Outpath.c_str(),stype.c_str()));	
 } //end loop on energies

} // end function
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
