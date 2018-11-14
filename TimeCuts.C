////////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Sarah Mechbal, smechbal@ucsc.edu
///    Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, July 3rd, 2018
////////////////////////////////////////////////////////////////////////////////////////////////////////

//This script is used to create two datasets from flight data: one for the nighttime spectrum, the other for the daytime
#include "headers.h"
#include "ALEvent.h"
#include "TMath.h"
#include "TGaxis.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)


void TimeCuts(string RecoID) {

double Tstart = 0;
double Tend = 134;	//total time of flight in hours)
double Launch = 22 + (7/60);
double Launchday = 15;
int NFiles = 5;	
	
string ID = RecoID;	

TFile *fileCal;	
TFile *fileNight = new TFile(Form("/home/sarah/AESOPLITE/FlightData/BBR2/TimeCuts/18A1.BPD.EVENT_%s_NightTime.root",ID.c_str()),"RECREATE");
TFile *fileDay = new TFile(Form("/home/sarah/AESOPLITE/FlightData/BBR2/TimeCuts/18A1.BPD.EVENT_%s_DayTime.root",ID.c_str()),"RECREATE"); 

	
// Create a TTree for NightTime files
TTree *NTtree = new TTree("NightTime"," Data Event Night Time");
//Define variables to make the Reco event
ALEvent *nt=new ALEvent();
// Create a branch with event
NTtree->Branch("NightTime",&nt); 	

// Create a TTree for DayTime files
TTree *DTtree = new TTree("DayTime"," Data Event Day Time");
//Define variables to make the Reco event
ALEvent *dt=new ALEvent();
// Create a branch with event
DTtree->Branch("DayTime",&dt); 
	
//ROOT filein
	
TChain*chain=new TChain("Data");
chain->Add(Form("/home/sarah//AESOPLITE/FlightData/BBR2/PRonly_Calibrated/*.root"));
//Define variables to read event
//Set address to access event data
ALEvent *de = new ALEvent();
chain->SetBranchAddress("Calevent",&de); 
int nentries=chain->GetEntries();	


cout << "Number  of events: " << nentries << endl;  
	
for(int i=0; i<nentries; i++) 
//for(int i=0;i<1000000; i++)
{
	chain->GetEntry(i);
    if(i%1000==0) cout << "Event " << i << endl;

//Time of event from CT1
	double yCT1 = (double)de->get_yCT1();
	double mCT1 = (double)de->get_mCT1();
	double dCT1 = (double)de->get_dCT1();
	double hCT1= (double)de->get_hCT1();
	double miCT1 = (double)de->get_miCT1();
	double sCT1 = (double)de->get_sCT1();
	float tmptime = hCT1 + (miCT1/60) + (sCT1/3600);
	tmptime = tmptime-Launch;
	tmptime = tmptime + 24*(dCT1 - Launchday);
	
	//Make time cuts and copy to new ROOT files
	
	// Fill Daytime file
	if(tmptime > 0) {
	if(( tmptime <19.8) || (tmptime > 30 & tmptime < 42) || (tmptime > 56 & tmptime <65) || (tmptime > 82 & tmptime < 90))
		{
			//cout << "Time = " << tmptime << "h, filling DAYTIME file" << endl;
		  	dt=new ALEvent();
    		dt->Copy(de); 
			DTtree->Fill();
			delete dt;
		}
    //Fill NightTime file
	else 
		{
		//	cout << "Time = " << tmptime << "h, filling NIGHTTIME file" << endl;
		    nt=new ALEvent();
    		nt->Copy(de); 
			NTtree->Fill();
			delete nt;	
		}
	}

   }
	//end i, loop on events

cout << "end loop on events" << endl;

//Write tree in output files
cout << "opening fileDay" << endl;
fileDay->cd();
DTtree->Write();	
cout << "Write DT" << endl;
cout << "FileDay closed" << endl;


fileNight->cd();
NTtree->Write();	
cout << "Write NT" << endl;
fileNight->Close();
fileDay->Close();



} //end function	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	