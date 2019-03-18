////////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Sarah Mechbal, smechbal@ucsc.edu
///    Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, July 1st, 2018
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "ALEvent.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)
	

			
void Calibrate(string RecoID);
float CalBar1(float B1);
float CalBar2(float B2);

//polynomial functions for PHA calibrations
double fT1[7] = {-5.97251215e-19, 5.69598557e-15, -1.80443371e-11, 2.38516848e-08,-1.02964831e-05, 2.55724193e-01 ,-1.30215976e+00};
double fT1lin[2] = {3.939, 4.687};
double fT2[7] = {-2.67643335e-19,2.70130010e-15,-8.39634550e-12,9.79081638e-09,-1.71828195e-06 ,2.39692564e-01,2.86714168e-01};
double fT2lin[2] = {4.170, -1.099};
double fT3[7] = {-4.30463942e-19,2.79141311e-15,-2.09123453e-12,-8.98195046e-09,1.68429543e-05,2.40090583e-01,-1.28605349e+00};
double fT3lin[2] = {4.072, 6.693};
double fT4[7] = {-2.41073434e-19,2.15753396e-15,-4.80080031e-12,5.71441605e-10,8.21287524e-06,2.34287148e-01,-7.36243155e-01};
double fT4lin[2] = {4.211, 3.991};
double fG[7] = {-4.71822390e-19,3.88281288e-15,-8.79628332e-12,5.93250483e-09,2.98770348e-06,2.42955811e-01,-7.38266412e-01};
double fGlin[2] = {4.094, 3.353};

double Pol6(double *fParam,double x) {
	double y = fParam[0]*pow(x,6)+fParam[1]*pow(x,5)+fParam[2]*pow(x,4)+fParam[3]*pow(x,3)+fParam[4]*pow(x,2)+fParam[5]*pow(x,1) + fParam[6];
	return y;
}

double Pol1(double *fParam,double x) {
	double y = fParam[0]*x + fParam[1];
	return y;
}
double CalT1(double T1) {
	//from PHA to mV eq.
	double mVCal = Pol6(fT1,T1);
	double T1Cal = Pol1(fT1lin,mVCal);
	return T1Cal;
}

double CalT2(double T2) {
	//from PHA to mV eq.
	double mVCal = Pol6(fT2,T2);
	double T2Cal = Pol1(fT2lin,mVCal);
	return T2Cal;
}
double CalT3(double T3) {
	//from PHA to mV eq.
	double mVCal = Pol6(fT3,T3);
	double T3Cal = Pol1(fT3lin,mVCal);
	return T3Cal;
}
double CalT4(double T4) {
	//from PHA to mV eq.
	double mVCal = Pol6(fT4,T4);
	double T4Cal = Pol1(fT4lin,mVCal);
	return T4Cal;
}
double CalG(double G) {
	//from PHA to mV eq.
	double mVCal = Pol6(fG,G);
	double GCal = Pol1(fGlin,mVCal);
	return GCal;
}


//Define calibration functions from fits

float CalBar1(float B1) {
	float B1Cal;
	if(B1 < 599.80) {
		B1Cal = B1 - 0.298;
	}
	else {
		B1Cal = (0.916 * B1) + 50.159;
	}
	return B1Cal;
}

float CalBar2(float B2) {
	float B2Cal;
	if(B2 < 6.46) {
		B2Cal = (1.004*B2) - 0.335;
	}
	else {
		B2Cal = B2 - 0.312;
	}
	return B2Cal;
}


void Calibrate(string RecoID) {
	
string ID = RecoID;
string Inpath = "/data/smechbal/Data/FlightData/18A1_SplitBPD";
	
TFile *filein;
TFile *fileout;
for (int j=0; j<22; j++) {
 if (j<10) {
	//ROOT input file (uncalibrated)
filein=new TFile(Form("%s/18A1_00%d.BPD.EVENT_%s.root",Inpath.c_str(),j,ID.c_str()),"READ");
cout << "Input file is open" <<endl;
//ROOT output file
fileout=new TFile(Form("%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),j,ID.c_str()),"RECREATE");
cout << "Output file is open" << endl;
 }
	 if (j>=10) {
	//ROOT input file (uncalibrated)
filein=new TFile(Form("%s/18A1_0%d.BPD.EVENT_%s.root",Inpath.c_str(),j,ID.c_str()),"READ");
cout << "Input file is open" <<endl;
//ROOT output file
fileout=new TFile(Form("%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),j,ID.c_str()),"RECREATE");
cout << "Output file is open" << endl;
 }
//Get Tree from the input file
TTree *tree = (TTree*)filein->Get("Data");
//Define variables to read event
ALEvent *e = new ALEvent();
//Set address to access event data
tree->SetBranchAddress("event",&e);  	
// Create a TTree
TTree *DEtree = new TTree("Data"," Data Event");
//Define variables to make the Reco event
ALEvent *de=new ALEvent();
// Create a branch with event
DEtree->Branch("Calevent",&de); 	

int nentries=tree->GetEntries();
cout << "Number  of events: " << nentries << endl;  
	
for(int i=0; i<nentries; i++) 
{
	
	tree->GetEntry(i);
    if(i%100==0) cout << "Event " << i << endl;
	de=new ALEvent();
    de->Copy(e);  		
	
//Calibrate barometers
	float B1Pres = de->get_PressB1();
	float B2Pres = de->get_PressB2();
	float B1Cal = CalBar1(B1Pres);
	float B2Cal = CalBar2(B2Pres);
	de->set_PressB1(B1Cal);
	de->set_PressB2(B2Cal);

//Calibrate PHAs
	double T1 = de->get_EneT1().at(0);
	double T2 = de->get_EneT2().at(0);
	double T3 = de->get_EneT3().at(0);
	double T4 = de->get_EneT4().at(0);
	double G = de->get_Eneg().at(0);
	double T1Cal = CalT1(T1);
	double T2Cal = CalT2(T2);
	double T3Cal = CalT3(T3);
	double T4Cal = CalT4(T4);
	double GCal = CalG(G);
	if(T1>0) {
		de->get_EneT1().clear();
		de->add_EneT1(T1Cal);
	}
	if(T2>0) {
		de->get_EneT2().clear();
		de->add_EneT2(T2Cal);
	}
	if(T3>0) {
		de->get_EneT3().clear();
		de->add_EneT3(T3Cal);
	}
	if(T4>0) {
		de->get_EneT4().clear();
		de->add_EneT4(T4Cal);
	}
	if(G>0) {
		de->get_Eneg().clear();
		de->add_Eneg(GCal);
	}
      DEtree->Fill();
        //Free memory
        delete de;

} //end i
	
	
	//Write tree in output file
 	fileout->cd();
	DEtree->Fill();
 	DEtree->Write(0,TObject::kOverwrite);
 	//Close files
 	fileout->Close();
 	filein->Close();    
	} //end loop on flight files
	
}

void PressGraph(string RecoID) {

		
string ID = RecoID;	

//define vectors

vector <float> B2Day;
vector <float> B2Night;
vector <float> TimeDay;
vector <float> TimeNight;

double Launch = 22 + (7/60);
double Launchday = 15;

//ROOT filein
	
TChain*chainDT=new TChain("DayTime");
chainDT->Add(Form("/home/sarah/AESOPLITE/FlightData/BBR2/TimeCuts/*.root"));
ALEvent *dt = new ALEvent();
chainDT->SetBranchAddress("DayTime",&dt); 
int DTentries=chainDT->GetEntries();	
cout << "Number  of events: " << DTentries << endl;  

TChain*chainNT=new TChain("NightTime");
chainNT->Add(Form("/home/sarah/AESOPLITE/FlightData/BBR2/TimeCuts/*.root"));
ALEvent *nt = new ALEvent();
chainNT->SetBranchAddress("NightTime",&nt); 
int NTentries=chainNT->GetEntries();	
cout << "Number  of events: " << NTentries << endl;  
	

	
for(int i=0; i<NTentries; i++) 
{
	
	chainNT->GetEntry(i);
    if(i%100==0) cout << "Event " << i << endl;
		
	
//Get Barometer

	float B2PressNight = nt->get_PressB2();
	double yNCT1 = (double)nt->get_yCT1();
	double mNCT1 = (double)nt->get_mCT1();
	double dNCT1 = (double)nt->get_dCT1();
	double hNCT1= (double)nt->get_hCT1();
	double miNCT1 = (double)nt->get_miCT1();
	double sNCT1 = (double)nt->get_sCT1();
	float tmptimeNight = hNCT1 + (miNCT1/60) + (sNCT1/3600);
	tmptimeNight = tmptimeNight-Launch;
	tmptimeNight = tmptimeNight + 24*(dNCT1 - Launchday);
	cout << " NightTime t = " << tmptimeNight << "h" << endl;


	if (tmptimeNight > 3) {
		B2Night.push_back(B2PressNight);
		TimeNight.push_back(tmptimeNight);
	}
} //end i

for(int i=0; i<DTentries; i++) 
{
	
	chainDT->GetEntry(i);
    if(i%100==0) cout << "Event " << i << "h" << endl;
		
	
//Get Barometer

	float B2PressDay = dt->get_PressB2();
	double yCT1 = (double)dt->get_yCT1();
	double mCT1 = (double)dt->get_mCT1();
	double dCT1 = (double)dt->get_dCT1();
	double hCT1= (double)dt->get_hCT1();
	double miCT1 = (double)dt->get_miCT1();
	double sCT1 = (double)dt->get_sCT1();
	float tmptimeDay = hCT1 + (miCT1/60) + (sCT1/3600);
	tmptimeDay = tmptimeDay-Launch;
	tmptimeDay = tmptimeDay + 24*(dCT1 - Launchday);
	cout << " Daytime t = " << tmptimeDay << "h" << endl;


	if (tmptimeDay > 3) {
		B2Day.push_back(B2PressDay);
		TimeDay.push_back(tmptimeDay);
	}
} //end i
	
//Graph 
TCanvas *c1 = new TCanvas("c1","c1",4000,1000);
TMultiGraph *PressTime = new TMultiGraph();
int nsizeD = TimeDay.size();
TGraph *gB2Day = new TGraph(nsizeD,&TimeDay[0],&B2Day[0]);
gB2Day->SetMarkerStyle(kFullCircle);
gB2Day->SetMarkerSize(0.1);
gB2Day->SetMarkerColor(kRed);
int nsizeN = TimeNight.size();
TGraph *gB2Night = new TGraph(nsizeN,&TimeNight[0],&B2Night[0]);
gB2Night->SetMarkerStyle(kFullCircle);
gB2Night->SetMarkerSize(0.1);
gB2Night->SetMarkerColor(kGreen);
PressTime->Add(gB2Day);
PressTime->Add(gB2Night);
c1->cd();
PressTime->Draw("AP");
c1->SaveAs("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/TimeCutsTest.eps");


} 	//end function
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		


