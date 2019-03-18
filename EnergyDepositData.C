////////////////////////////////////////////////////////////////////////////////////////// 
////   Author: Sarah Mechbal, smechbal@ucsc.edu
////   Department of Physics, University of California Santa Cruz, August 2018
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)



void EneDepData(string RecoID) 
{

	
double Tstart = 0;
double Tend = 134;	//total time of flight in hours)
double Launch = 22 + (7/60);
double LaunchDay = 15 + (Launch*0.041667);
string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2";
string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/PHA";
TFile *fileRaw;
TFile *fileCal;
string ID = RecoID;	
int StartFile = 0;
int EndFile = 22;
int NFiles = 23;
int NTConfig=5;//Number of configuration of triggers
string TConfig[5]={"T1","T2","T3","T4", "Guard"};
string CutConfig="AllFlight";


//Histograms of energy deposited in scintillators in PHA value
 TH1F**TrigEne=new TH1F*[NTConfig];
 
 for(int i=0;i<NTConfig;i++)
   {
    TrigEne[i]=new TH1F(Form("Energy %s",TConfig[i].c_str()),Form("Energy %s",TConfig[i].c_str()),100,0.,300);
	TrigEne[i]->GetXaxis()->SetTitle("Deposited energy (PHA value)");
    TrigEne[i]->GetYaxis()->SetTitle("Events");
    TrigEne[i]->GetYaxis()->SetTitleOffset(2.);
    TrigEne[i]->SetLineColor(kBlack);     
   	} 
	
	//ROOT fileout
	 TFile *fileout = new TFile(Form("%s/FlightPHA_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
	for (int j=StartFile; j<EndFile; j++) {
	 if (j<10) {
	fileRaw=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 }
		 if (j>=10) {
	fileRaw=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 }

	//Get Tree from the Raw file
	TTree *treeRaw = (TTree*)fileRaw->Get("Data");
	ALEvent *e = new ALEvent();
	treeRaw->SetBranchAddress("event",&e);  	

	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	

	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
	{
		treeRaw->GetEntry(i);
		treeCal->GetEntry(i);
		if(i%1000==0) cout << "Event " << i << endl;
	//Trigger information and layers information
		uint8_t Ti=(uint8_t)de->get_Ti();
		int  nhits = de->get_Nhits();
		int NL= de->get_NLayers();
		bool T1=de->get_T1();
		bool T2=de->get_T2();
		bool T3=de->get_T3();
		bool T4=de->get_T4();
		bool G=de->get_guard();
	//Time of event from CT1
		double yCT1 = (double)de->get_yCT1();
		double mCT1 = (double)de->get_mCT1();
		double dCT1 = (double)de->get_dCT1();
		double hCT1= (double)de->get_hCT1();
		double miCT1 = (double)de->get_miCT1();
		double sCT1 = (double)de->get_sCT1();
		float tmptime = hCT1 + (miCT1/60) + (sCT1/3600);
		tmptime = tmptime-Launch;
		tmptime = tmptime + 24*(dCT1 - LaunchDay);
	//PHA signals (calibrated)
		double T1PHACal = de->get_EneT1().at(0);
		double T2PHACal = de->get_EneT2().at(0);
		double T3PHACal = de->get_EneT3().at(0);
		double T4PHACal = de->get_EneT4().at(0);
		double GPHACal = de->get_Eneg().at(0);
	//PHA signals (uncalibrated)
		double T1PHA = e->get_EneT1().at(0);
		double T2PHA = e->get_EneT2().at(0);
		double T3PHA = e->get_EneT3().at(0);
		double T4PHA = e->get_EneT4().at(0);
		double GPHA = e->get_Eneg().at(0);
	//Barometer data (calibrated)
		float B2Pres = de->get_PressB2();
	//Convert pressure to g/cm^-2
		B2Pres = B2Pres* 1.3595;
	//Reconstructed energy 
		float pPR=1000*fabs(de->get_p0PR());   //in MeV
		double deflecPR = de->get_deflecPR();
		float P0sign = pPR*TMath::Sign(1,deflecPR);
		double chiNBPR = de->get_chi2NBPR();
		double chiBPR = de->get_chi2BPR();

		//T1
		TrigEne[0]->Fill(T1PHACal);
		//T2
		TrigEne[1]->Fill(T2PHACal);
		//T3
		TrigEne[2]->Fill(T3PHACal);
		//T4
		TrigEne[3]->Fill(T4PHACal);
		//Guard
		TrigEne[4]->Fill(GPHACal);
	 } //end i, loop over entries
   } //end j, loop over files


	cout << "done looping the events" << endl;
	
	
//////////////////////////////////
////////Histogram Display/////////
//////////////////////////////////
	
 fileout->cd();
 TCanvas*canTrigEne = new TCanvas(Form("canTrigEne"), Form("canTrigEne"),200,10,1500,500); 
 canTrigEne->Divide(2,3);
 for(int i=0;i<NTConfig;i++) 
	{

	canTrigEne->cd(i+1); 
	gPad->SetLeftMargin(0.15);
	TrigEne[i]->GetYaxis()->SetTitleOffset(2);  
	//Rebin histograms
//	int newbin = TMath::Sqrt(TrigEne[i]->GetEntries());
//	double binratio = 200/newbin;
//	TrigEne[i]->Rebin(1.3*binratio);
	TrigEne[i]->Draw("hists");
	TrigEne[i]->Write();
	}  
	canTrigEne->Write();
	fileout->Close();
	////////////////////
	//Save plots into a pdf file
	//////////////////  
	
	canTrigEne->Print(Form("%s/FlightPHA_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));	
	

} // end function
	
	
void EneDepDataMomBin(string RecoID, double cutT2) 
{

	
double Tstart = 0;
double Tend = 134;	//total time of flight in hours)
double Launch = 22 + (7/60);
double LaunchDay = 15 + (Launch*0.041667);
string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2";
string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/PHA";
TFile *fileRaw;
TFile *fileCal;
string ID = RecoID;	
int StartFile = 0;
int EndFile = 23;
int NFiles = 23;
int NTConfig=4;//Number of configuration of triggers
string TConfig[4]={"T1","T2","T3","T4"};
string CutConfig="AllFlightT2Above120NoGHits5to15";
const int NPBins = 12;
double PBins[NPBins] = {-300,-150,-80,-50,-30,-10,10,30,50,80,150,300};
double T2Cut = cutT2;	//value of cut on T2

//Histograms of energy deposited in scintillators in PHA value, for each momentum bin
 TH1F***TrigEne=new TH1F**[NPBins];
 
 for(int i=0;i<NPBins;i++)
   {
	TrigEne[i] = new TH1F*[NTConfig];
	for(int j=0;j<NTConfig;j++) 
	 {
//more bins for T2		
	if(j==1) TrigEne[i][j]=new TH1F(Form("%s %3.0f < pMom <%3.0f ",TConfig[j].c_str(),PBins[i],PBins[i+1]),Form("%s %3.0f < pMom <%3.0f",TConfig[j].c_str(),PBins[i],PBins[i+1]),200,0.,1000);
    else TrigEne[i][j]=new TH1F(Form("%s %3.0f < pMom <%3.0f ",TConfig[j].c_str(),PBins[i],PBins[i+1]),Form("%s %3.0f < pMom <%3.0f",TConfig[j].c_str(),PBins[i],PBins[i+1]),100,0.,300);
	TrigEne[i][j]->GetXaxis()->SetTitle("Deposited energy (PHA value)");
    TrigEne[i][j]->GetYaxis()->SetTitle("Events");
    TrigEne[i][j]->GetYaxis()->SetTitleOffset(2.);
    TrigEne[i][j]->SetLineColor(kBlack);     
   	} 
 }
	//ROOT fileout
	 TFile *fileout = new TFile(Form("%s/BinnedFlightPHA_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
	for (int j=StartFile; j<EndFile; j++) {
	 if (j<10) {
	fileRaw=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 }
		 if (j>=10) {
	fileRaw=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 }

	//Get Tree from the Raw file
	TTree *treeRaw = (TTree*)fileRaw->Get("Data");
	ALEvent *e = new ALEvent();
	treeRaw->SetBranchAddress("event",&e);  	

	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	

	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
	{
		treeRaw->GetEntry(i);
		treeCal->GetEntry(i);
		if(i%1000==0) cout << "Event " << i << endl;
	//Trigger information and layers information
		uint8_t Ti=(uint8_t)de->get_Ti();
		int  nhits = de->get_Nhits();
		int NL= de->get_NLayers();
		bool T1=de->get_T1();
		bool T2=de->get_T2();
		bool T3=de->get_T3();
		bool T4=de->get_T4();
		bool G=de->get_guard();
	//Time of event from CT1
		double yCT1 = (double)de->get_yCT1();
		double mCT1 = (double)de->get_mCT1();
		double dCT1 = (double)de->get_dCT1();
		double hCT1= (double)de->get_hCT1();
		double miCT1 = (double)de->get_miCT1();
		double sCT1 = (double)de->get_sCT1();
		float tmptime = hCT1 + (miCT1/60) + (sCT1/3600);
		tmptime = tmptime-Launch;
		tmptime = tmptime + 24*(dCT1 - LaunchDay);
	//PHA signals (calibrated)
		double T1PHACal = de->get_EneT1().at(0);
		double T2PHACal = de->get_EneT2().at(0);
		double T3PHACal = de->get_EneT3().at(0);
		double T4PHACal = de->get_EneT4().at(0);
		double GPHACal = de->get_Eneg().at(0);
	//PHA signals (uncalibrated)
		double T1PHA = e->get_EneT1().at(0);
		double T2PHA = e->get_EneT2().at(0);
		double T3PHA = e->get_EneT3().at(0);
		double T4PHA = e->get_EneT4().at(0);
		double GPHA = e->get_Eneg().at(0);
	//Barometer data (calibrated)
		float B2Pres = de->get_PressB2();
	//Convert pressure to g/cm^-2
		B2Pres = B2Pres* 1.3595;
	//Reconstructed energy 
		float pPR=1000*fabs(de->get_p0PR());   //in MeV
		double deflecPR = de->get_deflecPR();
		float P0sign = pPR*TMath::Sign(1,deflecPR);
		double chiNBPR = de->get_chi2NBPR();
		double chiBPR = de->get_chi2BPR();
		if(deflecPR==0) continue;
	//Determine the correct momentum bin
		int pbin;
		for(int k=0;k<NPBins-1;k++) 
			{

			if((PBins[k] <= P0sign) & (P0sign < PBins[k+1])) 
				{
				pbin = k;
			  //  cout << PBins[k] << " < " << P0sign << " < " << PBins[k+1] << endl;
			//	cout << "pbin = " << pbin << endl;
				break;
			 }
			}		
		//	cout << "Low edge of energy bin is " << PBins[pbin] << endl;

		//Apply cuts to only get clean events
		if(T1 & T2 & T3 & T2PHACal>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
		   {

		//   cout << "made the cut" << endl;
		   TrigEne[pbin][0]->Fill(T1PHACal);
		   TrigEne[pbin][1]->Fill(T2PHACal);
		   TrigEne[pbin][2]->Fill(T3PHACal);
		   TrigEne[pbin][3]->Fill(T4PHACal);

			}	 
	 } //end i, loop over entries
   } //end j, loop over files


	cout << "done looping the events" << endl;
	
	
//////////////////////////////////
////////Histogram Display/////////
//////////////////////////////////
	
 fileout->cd();
 TCanvas**canTrigEne = new TCanvas*[NPBins];
 for(int i=0;i<NPBins-1;i++) 
	{
	 if(i==5) continue;
	 canTrigEne[i] = new TCanvas(Form("canTrigEne %d",i), Form("canTrigEne %d",i),200,10,1500,500); 
     canTrigEne[i]->Divide(2,2);
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

	if(i==0) canTrigEne[i]->Print(Form("%s/BinnedFlightPHA_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	else if(i==NPBins-2) canTrigEne[i]->Print(Form("%s/BinnedFlightPHA_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	else canTrigEne[i]->Print(Form("%s/BinnedFlightPHA_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
 } //end loop on bins
    fileout->Close();	

} // end function
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
