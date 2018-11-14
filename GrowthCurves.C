////////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Sarah Mechbal, smechbal@ucsc.edu
///    Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, July 3rd, 2018
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "ALEvent.h"
#include "TMath.h"
#include "TGaxis.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)


double Tstart = 0;
double Tend = 134;	//total time of flight in hours)
double Launch = 22 + (7/60);
double Launchday = 15;
//Make time steps of 20mn hour to start:
double binsperhour = 3;
int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
const int NPBins = 10;
double PBins[NPBins] = {-300,-150,-80,-40,-20,20,40,80,150,300};
int NFiles = 23;
double TCut[7] = {19.8,30,42,56,65,82,90};

void GrowthCurves(string RecoID) {
cout << "There are " << NTimeBin << " time bins " << endl;
	
TFile *fileRaw;
TFile *fileCal;
string ID = RecoID;	
	
//Create histograms
TH1F**B2 = new TH1F*[NTimeBin];
TH1F***P0Bin = new TH1F**[NPBins];
for(int k=0; k<NTimeBin;k++) 
	{
		B2[k] = new TH1F(Form("B2 h%d",k),Form("B2 h%d",k), 10000, 0, 800);
	}

//Momentum histograms
for(int i=0; i<NPBins; i++) 
{
		P0Bin[i] = new TH1F*[NTimeBin];
	
		for(int j=0; j<NTimeBin; j++) 
		 {
			P0Bin[i][j] = new TH1F(Form("signed p0 h%d%d",i,j), Form("signed p0 h%d%d",i,j), 1,0,1);
		   }
	}
	
double Launch = 22 + (7/60);
double LaunchDay = 15 + (Launch*0.041667);

//ROOT fileout
 TFile *fileout = new TFile(Form("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/GrowthCurves_%s.root", ID.c_str()),"RECREATE");
for (int j=0; j<NFiles; j++) {
 if (j<10) {
fileRaw=new TFile(Form("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly/18A1_00%d.BPD.EVENT_%s.root",j,ID.c_str()),"READ");
fileCal=new TFile(Form("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/18A1_00%d.BPD.EVENT_%s_Calibrated.root",j,ID.c_str()),"READ");
 }
	 if (j>=10) {
fileRaw=new TFile(Form("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly/18A1_0%d.BPD.EVENT_%s.root",j,ID.c_str()),"READ");
fileCal=new TFile(Form("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/18A1_0%d.BPD.EVENT_%s_Calibrated.root",j,ID.c_str()),"READ");
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
	tmptime = tmptime + 24*(dCT1 - Launchday);
//PHA signals (calibrated)
    double T1PHA = de->get_EneT1().at(0);
	double T2PHA = de->get_EneT2().at(0);
	double T3PHA = de->get_EneT3().at(0);
	double T4PHA = de->get_EneT4().at(0);
	double GPHA = de->get_Eneg().at(0);

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
//Determine the correct time and energy bin
	int tbin = tmptime*binsperhour;
	int pbin;
    //	cout << "Assign time bin " << tbin <<endl;

	for(int k=0;k<NPBins;k++) 
		{
		if((PBins[k] <= P0sign) & (P0sign < PBins[k+1])) 
			{
			pbin = k;
	     //  cout << "NEXT" << endl;
		   // cout << PBins[k] << " =< " << P0sign << " < " << PBins[k+1] << endl;
		 }
		}
	if(tmptime>0) 
		{
	B2[tbin]->Fill(B2Pres);
//Apply cuts to only get clean events
    if(T1 & T2 & T3 & (T2PHA>0) & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=7) && (nhits <15))
	   {
    //	cout << "Assign bin number energy bin " << pbin << ", time bin " << tbin <<endl;
	//	cout << "Fill event P0sign = " << P0sign << "MeV, at time  =  " << tmptime << " hours" <<endl;
	//	cout << "Low edge of energy bin is " << P0Bin[pbin][tbin]->GetBinLowEdge(1) << endl;
		P0Bin[pbin][tbin]->Fill(P0sign);

		}	   
	}
	
		
	}     //end i, loop on events
	cout << "finished reading a file" << endl;
}		//end j, loop on files

//Make growth curves
int ngraphs = (NPBins/2)-1;
TGraph**CountRate = new TGraph*[NPBins];
TGraph**TimeSeries = new TGraph*[NPBins];
TGraph *PressB2vsTime = new TGraph();
PressB2vsTime->SetMarkerStyle(kOpenCircle);
//PressB2vsTime->SetMarkerSize(0.1);
PressB2vsTime->SetMarkerColor(kBlue);	

for (int j=0; j<NPBins;j++) 
	{
		CountRate[j] = new TGraph();
	    CountRate[j]->SetMarkerStyle(kOpenCircle);
		CountRate[j]->SetMarkerSize(0.8);
		TimeSeries[j] = new TGraph();
	    TimeSeries[j]->SetMarkerStyle(kDot);
		TimeSeries[j]->SetMarkerSize(0.8);
		if (j<4)  { 
		CountRate[j]->SetMarkerColor(kBlue);
		TimeSeries[j]->SetMarkerColor(kBlue);
		TimeSeries[j]->SetLineColor(kBlue);


			}
		if (j>4)  { 
		CountRate[j]->SetMarkerColor(kRed);
		TimeSeries[j]->SetMarkerColor(kRed);
		TimeSeries[j]->SetLineColor(kRed);

			}
		int npPos = 1;
		int npEle = 1;
		int np2 = 1;
		for (int k=0;k<NTimeBin;k++) 
		{
			double RatePos;
			double RateEle;
			double B2mean = B2[k]->GetMean();
			double time = Tstart + (k+1)/(binsperhour);
//Time since launch expressed in units of days
			double timeInDays = LaunchDay + time*(0.041667);
			if (time>0) {
			PressB2vsTime->SetPoint(np2, time, B2mean);
			np2++;
			}
			if (j<4) 
				{
				RateEle = P0Bin[j][k]->GetEntries(); //count rate for 30 min.
				CountRate[j]->SetPoint(npEle,B2mean,RateEle);
				//TimeSeries[j]->SetPoint(npEle,timeInDays,RateEle);
				TimeSeries[j]->SetPoint(npEle,time,RateEle);
				npEle++;
				}
		   if(j>4)
		   	{
			   RatePos = P0Bin[j][k]->GetEntries(); //count rate for 30 min.
			   CountRate[j]->SetPoint(npPos,B2mean,RatePos);
     		   TimeSeries[j]->SetPoint(npPos,time,RatePos);
			   npPos++;
			  
		   	}
		}
    }

	/*
cout << "number of points in one graph = " << CountRate[0]->GetN() << endl;
	
gStyle->SetTextFont(85);
TCanvas *c0 = new TCanvas("c0","c0",1100,2000);
//First subplot
TMultiGraph*GrowthCurves0 = new TMultiGraph();
TLegend *leg0 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves0->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves0->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves0->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurves0->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves0->Add(CountRate[4]);
GrowthCurves0->Add(CountRate[6]);
leg0->AddEntry(CountRate[4],"e- 10-30 MeV","p");
leg0->AddEntry(CountRate[6],"e+ 10-30 MeV","p");
c0->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurves0->SetMinimum(1);
GrowthCurves0->Draw("AP");	
leg0->Draw("same");
c0->Modified(); 
c0->Update();
GrowthCurves0->GetXaxis()->SetRangeUser(1,1000);
//Second subplot
TCanvas *c1 = new TCanvas("c1","c1",1100,2000);
TMultiGraph*GrowthCurves1 = new TMultiGraph();
TLegend *leg1 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves1->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves1->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurves1->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves1->Add(CountRate[3]);
GrowthCurves1->Add(CountRate[7]);
leg1->AddEntry(CountRate[3],"e- 30-50 MeV","p");
leg1->AddEntry(CountRate[7],"e+ 30-50 MeV","p");
c1->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurves1->SetMinimum(1);
GrowthCurves1->Draw("AP");		
leg1->Draw("same");
c1->Modified(); 
c1->Update();
GrowthCurves1->GetXaxis()->SetRangeUser(1,1000);
//Third subplot
TCanvas *c2 = new TCanvas("c2","c2",1100,2000);
TMultiGraph*GrowthCurves2 = new TMultiGraph();
TLegend *leg2 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves2->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves2->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves2->GetXaxis()->SetTitle("Atmospheric pressure (in g.cm^{-2})");
GrowthCurves2->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves2->Add(CountRate[2]);
GrowthCurves2->Add(CountRate[8]);
leg2->AddEntry(CountRate[2],"e- 50-80 MeV","p");
leg2->AddEntry(CountRate[8],"e+ 50-80 MeV","p");
c2->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1); 
gStyle->SetPadTickY(1);
GrowthCurves2->SetMinimum(1);
GrowthCurves2->Draw("AP");		
leg2->Draw("same");
c2->Modified(); 
c2->Update();
GrowthCurves2->GetXaxis()->SetRangeUser(1,1000);
//Fourth subplot
TCanvas *c3 = new TCanvas("c3","c3",1100,2000);
TMultiGraph*GrowthCurves3 = new TMultiGraph();
TLegend *leg3 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves3->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves3->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves3->GetXaxis()->SetTitle("Atmospheric pressure (in g.cm^{-2})");
GrowthCurves3->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves3->Add(CountRate[1]);
GrowthCurves3->Add(CountRate[9]);
leg3->AddEntry(CountRate[1],"e- 80-150 MeV","p");
leg3->AddEntry(CountRate[9],"e+ 80-150 MeV","p");
c3->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1); 
gStyle->SetPadTickY(1);
GrowthCurves3->SetMinimum(1);
GrowthCurves3->Draw("AP");		
leg3->Draw("same");
c3->Modified(); 
c3->Update();
GrowthCurves3->GetXaxis()->SetRangeUser(1,1000);
//Fifth and last subplot
TCanvas *c4 = new TCanvas("c4","c4",1100,2000);
TMultiGraph*GrowthCurves4 = new TMultiGraph();
TLegend *leg4 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves4->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves4->GetYaxis()->SetTitleOffset(0.8);
GrowthCurves4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurves4->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves4->Add(CountRate[0]);
GrowthCurves4->Add(CountRate[10]);
leg4->AddEntry(CountRate[0],"e- 150-300 MeV","p");
leg4->AddEntry(CountRate[10],"e+ 150-300 MeV","p");
c4->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1); 
gStyle->SetPadTickY(1);
GrowthCurves4->SetMinimum(1);
GrowthCurves4->Draw("AP");		
leg4->Draw("same");
c4->Modified(); 
c4->Update();
GrowthCurves4->GetXaxis()->SetRangeUser(1,1000);

c0->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf(");
c1->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf");
c2->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf");
c3->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf");
c4->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf)");
	*/
	
//Time series plot
gStyle->SetTextFont(85);
gStyle->SetPadTickX(0); 
gStyle->SetPadTickY(0);
/*
TCanvas *c00 = new TCanvas("c00","c00",4000,2000);
//First subplot
TMultiGraph*TimeSeries0 = new TMultiGraph();
TLegend *leg00 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries0->GetXaxis()->SetTitleOffset(1.2);
TimeSeries0->GetYaxis()->SetTitleOffset(0.7);
TimeSeries0->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries0->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries0->Add(TimeSeries[4]);
TimeSeries0->Add(TimeSeries[6]);
leg00->AddEntry(TimeSeries[4],"e- 10-30 MeV","l");
leg00->AddEntry(TimeSeries[6],"e+ 10-30 MeV","l");
c00->cd();
TimeSeries0->GetYaxis()->SetRangeUser(0,100);
TimeSeries0->GetXaxis()->SetLimits(0,133);
TimeSeries0->Draw("AP");		
leg00->Draw("same");
*/
	
//Second subplot
TCanvas *c11 = new TCanvas("c11","c11",4000,2000);
TMultiGraph*TimeSeries1 = new TMultiGraph();
TLegend *leg11 = new TLegend(0.6,0.85,0.9,0.9);
//Lines to show time cuts
TLine**LineCut = new TLine*[7];
TimeSeries1->GetXaxis()->SetTitleOffset(1.2);
TimeSeries1->GetYaxis()->SetTitleOffset(0.7);
TimeSeries1->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries1->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries1->Add(TimeSeries[3]);
TimeSeries1->Add(TimeSeries[5]);
leg11->AddEntry(TimeSeries[3],"e- 20-40 MeV","l");
leg11->AddEntry(TimeSeries[5],"e+ 20-40 MeV","l");
c11->cd();
TimeSeries1->GetYaxis()->SetRangeUser(0,100);
TimeSeries1->GetXaxis()->SetLimits(0,133);
TimeSeries1->Draw("ALP");		
leg11->Draw("same");
for (int i=0;i<7;i++) {
	LineCut[i] = new TLine(TCut[i],0.,1.001*TCut[i],100);
	LineCut[i]->Draw("same");
}
//Third subplot
TCanvas *c22 = new TCanvas("c22","c22",4000,2000);
TMultiGraph*TimeSeries2 = new TMultiGraph();
TLegend *leg22 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries2->GetXaxis()->SetTitleOffset(1.2);
TimeSeries2->GetYaxis()->SetTitleOffset(0.7);
TimeSeries2->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries2->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries2->Add(TimeSeries[2]);
TimeSeries2->Add(TimeSeries[6]);
leg22->AddEntry(TimeSeries[2],"e- 40-80 MeV","l");
leg22->AddEntry(TimeSeries[6],"e+ 40-80 MeV","l");
c22->cd();
TimeSeries2->GetYaxis()->SetRangeUser(0,100);
TimeSeries2->GetXaxis()->SetLimits(0,133);
TimeSeries2->Draw("ALP");		
leg22->Draw("same");
//Fourth subplot
TCanvas *c33 = new TCanvas("c33","c33",4000,2000);
TMultiGraph*TimeSeries3 = new TMultiGraph();
TLegend *leg33 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries3->GetXaxis()->SetTitleOffset(1.2);
TimeSeries3->GetYaxis()->SetTitleOffset(0.7);
TimeSeries3->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries3->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries3->Add(TimeSeries[1]);
TimeSeries3->Add(TimeSeries[7]);
leg33->AddEntry(TimeSeries[1],"e- 80-150 MeV","l");
leg33->AddEntry(TimeSeries[7],"e+ 80-150 MeV","l");
c33->cd();
TimeSeries3->GetYaxis()->SetRangeUser(0,100);
TimeSeries3->GetXaxis()->SetLimits(0,133);
TimeSeries3->Draw("ALP");		
leg33->Draw("same");
//Fifth subplot
TCanvas *c44 = new TCanvas("c44","c44",4000,2000);
TMultiGraph*TimeSeries4 = new TMultiGraph();
TLegend *leg44 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries4->GetXaxis()->SetTitleOffset(1.2);
TimeSeries4->GetYaxis()->SetTitleOffset(0.7);
TimeSeries4->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries4->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries4->Add(TimeSeries[0]);
TimeSeries4->Add(TimeSeries[8]);
leg44->AddEntry(TimeSeries[0],"e- 150-300 MeV","l");
leg44->AddEntry(TimeSeries[8],"e+ 150-300 MeV","l");
c44->cd();
TimeSeries4->GetYaxis()->SetRangeUser(0,100);
TimeSeries4->GetXaxis()->SetLimits(0,133);
TimeSeries4->Draw("ALP");		
leg44->Draw("same");

//c00->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeries.pdf(");
c11->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf(");
c22->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf");
c33->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf");
c44->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf)");

			
//Pressure vs time
TCanvas *CanPres = new TCanvas("CanPres","CanPress",4000,2000);
CanPres->cd();
PressB2vsTime->GetYaxis()->SetRangeUser(2.0,4.5);
PressB2vsTime->GetXaxis()->SetLimits(0,133);
PressB2vsTime->GetYaxis()->SetTitle("B2 calibrated pressure (g.cm^{-2})");
PressB2vsTime->GetXaxis()->SetTitle("Time since launch (in hours)");
PressB2vsTime->Draw("AP");
CanPres->SaveAs("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/PressVsTime.eps");
	
//Write histograms in ROOT file
	
fileout->cd();	
/*
c0->Write();
GrowthCurves0->Write();
c1->Write();
GrowthCurves1->Write();
c2->Write();
GrowthCurves2->Write();
c3->Write();
GrowthCurves3->Write();
c4->Write();
GrowthCurves4->Write();
*/

c11->Write();
TimeSeries1->Write();
c22->Write();
TimeSeries2->Write();
c33->Write();
TimeSeries3->Write();
c44->Write();
TimeSeries4->Write();

/*
for(int k=0;k<NPBins;k++) {
	for (int l=0; l<NTimeBin;l++) {
		P0Bin[k][l]->Write();
		}
	}
for (int l=0; l<NTimeBin;l++) 
	{
  B2[l]->Write();
	}
*/	
cout << "now closing all files" << endl;
fileout->Close();
fileRaw->Close();
fileCal->Close(); 	
	
} //end function	
	
void GrowthCurvesTimeCut(string RecoID) {
cout << "There are " << NTimeBin << " time bins " << endl;

string ID = RecoID;	
	
//Create histograms
TH1F**B2 = new TH1F*[NTimeBin];
TH1F***P0Bin = new TH1F**[NPBins];
for(int k=0; k<NTimeBin;k++) 
	{
		B2[k] = new TH1F(Form("B2 h%d",k),Form("B2 h%d",k), 10000, 0, 800);
	}

//Momentum histograms
for(int i=0; i<NPBins; i++) 
{
		P0Bin[i] = new TH1F*[NTimeBin];
	
		for(int j=0; j<NTimeBin; j++) 
		 {
			P0Bin[i][j] = new TH1F(Form("signed p0 h%d%d",i,j), Form("signed p0 h%d%d",i,j), 1,0,1);
		   }
	}
	
double Launch = 22 + (7/60);
double LaunchDay = 15 + (Launch*0.041667);

//ROOT fileout
 TFile *fileout = new TFile(Form("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/GrowthCurves_%s.root", ID.c_str()),"RECREATE");

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
	
/*	
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
	tmptime = tmptime + 24*(dCT1 - Launchday);
//PHA signals (calibrated)
    double T1PHA = de->get_EneT1().at(0);
	double T2PHA = de->get_EneT2().at(0);
	double T3PHA = de->get_EneT3().at(0);
	double T4PHA = de->get_EneT4().at(0);
	double GPHA = de->get_Eneg().at(0);

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
//Determine the correct time and energy bin
	int tbin = tmptime*binsperhour;
	int pbin;
    //	cout << "Assign time bin " << tbin <<endl;

	for(int k=0;k<NPBins;k++) 
		{
		if((PBins[k] <= P0sign) & (P0sign < PBins[k+1])) 
			{
			pbin = k;
	     //  cout << "NEXT" << endl;
		   // cout << PBins[k] << " =< " << P0sign << " < " << PBins[k+1] << endl;
		 }
		}
	if(tmptime>0) 
		{
	B2[tbin]->Fill(B2Pres);
//Apply cuts to only get clean events
    if(T1 & T2 & T3 & (T2PHA>0) & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=7) && (nhits <15))
	   {
    //	cout << "Assign bin number energy bin " << pbin << ", time bin " << tbin <<endl;
	//	cout << "Fill event P0sign = " << P0sign << "MeV, at time  =  " << tmptime << " hours" <<endl;
	//	cout << "Low edge of energy bin is " << P0Bin[pbin][tbin]->GetBinLowEdge(1) << endl;
		P0Bin[pbin][tbin]->Fill(P0sign);

		}	   
	}
	
		
	}     //end i, loop on events
	cout << "finished reading a file" << endl;
}		//end j, loop on files

//Make growth curves
int ngraphs = (NPBins/2)-1;
TGraph**CountRate = new TGraph*[NPBins];
TGraph**TimeSeries = new TGraph*[NPBins];
TGraph *PressB2vsTime = new TGraph();
PressB2vsTime->SetMarkerStyle(kOpenCircle);
//PressB2vsTime->SetMarkerSize(0.1);
PressB2vsTime->SetMarkerColor(kBlue);	

for (int j=0; j<NPBins;j++) 
	{
		CountRate[j] = new TGraph();
	    CountRate[j]->SetMarkerStyle(kOpenCircle);
		CountRate[j]->SetMarkerSize(0.8);
		TimeSeries[j] = new TGraph();
	    TimeSeries[j]->SetMarkerStyle(kDot);
		TimeSeries[j]->SetMarkerSize(0.8);
		if (j<4)  { 
		CountRate[j]->SetMarkerColor(kBlue);
		TimeSeries[j]->SetMarkerColor(kBlue);
		TimeSeries[j]->SetLineColor(kBlue);


			}
		if (j>4)  { 
		CountRate[j]->SetMarkerColor(kRed);
		TimeSeries[j]->SetMarkerColor(kRed);
		TimeSeries[j]->SetLineColor(kRed);

			}
		int npPos = 1;
		int npEle = 1;
		int np2 = 1;
		for (int k=0;k<NTimeBin;k++) 
		{
			double RatePos;
			double RateEle;
			double B2mean = B2[k]->GetMean();
			double time = Tstart + (k+1)/(binsperhour);
//Time since launch expressed in units of days
			double timeInDays = LaunchDay + time*(0.041667);
			if (time>0) {
			PressB2vsTime->SetPoint(np2, time, B2mean);
			np2++;
			}
			if (j<4) 
				{
				RateEle = P0Bin[j][k]->GetEntries(); //count rate for 30 min.
				CountRate[j]->SetPoint(npEle,B2mean,RateEle);
				//TimeSeries[j]->SetPoint(npEle,timeInDays,RateEle);
				TimeSeries[j]->SetPoint(npEle,time,RateEle);
				npEle++;
				}
		   if(j>4)
		   	{
			   RatePos = P0Bin[j][k]->GetEntries(); //count rate for 30 min.
			   CountRate[j]->SetPoint(npPos,B2mean,RatePos);
     		   TimeSeries[j]->SetPoint(npPos,time,RatePos);
			   npPos++;
			  
		   	}
		}
    }

	/*
cout << "number of points in one graph = " << CountRate[0]->GetN() << endl;
	
gStyle->SetTextFont(85);
TCanvas *c0 = new TCanvas("c0","c0",1100,2000);
//First subplot
TMultiGraph*GrowthCurves0 = new TMultiGraph();
TLegend *leg0 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves0->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves0->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves0->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurves0->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves0->Add(CountRate[4]);
GrowthCurves0->Add(CountRate[6]);
leg0->AddEntry(CountRate[4],"e- 10-30 MeV","p");
leg0->AddEntry(CountRate[6],"e+ 10-30 MeV","p");
c0->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurves0->SetMinimum(1);
GrowthCurves0->Draw("AP");	
leg0->Draw("same");
c0->Modified(); 
c0->Update();
GrowthCurves0->GetXaxis()->SetRangeUser(1,1000);
//Second subplot
TCanvas *c1 = new TCanvas("c1","c1",1100,2000);
TMultiGraph*GrowthCurves1 = new TMultiGraph();
TLegend *leg1 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves1->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves1->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurves1->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves1->Add(CountRate[3]);
GrowthCurves1->Add(CountRate[7]);
leg1->AddEntry(CountRate[3],"e- 30-50 MeV","p");
leg1->AddEntry(CountRate[7],"e+ 30-50 MeV","p");
c1->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurves1->SetMinimum(1);
GrowthCurves1->Draw("AP");		
leg1->Draw("same");
c1->Modified(); 
c1->Update();
GrowthCurves1->GetXaxis()->SetRangeUser(1,1000);
//Third subplot
TCanvas *c2 = new TCanvas("c2","c2",1100,2000);
TMultiGraph*GrowthCurves2 = new TMultiGraph();
TLegend *leg2 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves2->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves2->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves2->GetXaxis()->SetTitle("Atmospheric pressure (in g.cm^{-2})");
GrowthCurves2->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves2->Add(CountRate[2]);
GrowthCurves2->Add(CountRate[8]);
leg2->AddEntry(CountRate[2],"e- 50-80 MeV","p");
leg2->AddEntry(CountRate[8],"e+ 50-80 MeV","p");
c2->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1); 
gStyle->SetPadTickY(1);
GrowthCurves2->SetMinimum(1);
GrowthCurves2->Draw("AP");		
leg2->Draw("same");
c2->Modified(); 
c2->Update();
GrowthCurves2->GetXaxis()->SetRangeUser(1,1000);
//Fourth subplot
TCanvas *c3 = new TCanvas("c3","c3",1100,2000);
TMultiGraph*GrowthCurves3 = new TMultiGraph();
TLegend *leg3 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves3->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves3->GetYaxis()->SetTitleOffset(0.7);
GrowthCurves3->GetXaxis()->SetTitle("Atmospheric pressure (in g.cm^{-2})");
GrowthCurves3->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves3->Add(CountRate[1]);
GrowthCurves3->Add(CountRate[9]);
leg3->AddEntry(CountRate[1],"e- 80-150 MeV","p");
leg3->AddEntry(CountRate[9],"e+ 80-150 MeV","p");
c3->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1); 
gStyle->SetPadTickY(1);
GrowthCurves3->SetMinimum(1);
GrowthCurves3->Draw("AP");		
leg3->Draw("same");
c3->Modified(); 
c3->Update();
GrowthCurves3->GetXaxis()->SetRangeUser(1,1000);
//Fifth and last subplot
TCanvas *c4 = new TCanvas("c4","c4",1100,2000);
TMultiGraph*GrowthCurves4 = new TMultiGraph();
TLegend *leg4 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurves4->GetXaxis()->SetTitleOffset(1.2);
GrowthCurves4->GetYaxis()->SetTitleOffset(0.8);
GrowthCurves4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurves4->GetYaxis()->SetTitle("Count rate/20mn");
GrowthCurves4->Add(CountRate[0]);
GrowthCurves4->Add(CountRate[10]);
leg4->AddEntry(CountRate[0],"e- 150-300 MeV","p");
leg4->AddEntry(CountRate[10],"e+ 150-300 MeV","p");
c4->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1); 
gStyle->SetPadTickY(1);
GrowthCurves4->SetMinimum(1);
GrowthCurves4->Draw("AP");		
leg4->Draw("same");
c4->Modified(); 
c4->Update();
GrowthCurves4->GetXaxis()->SetRangeUser(1,1000);

c0->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf(");
c1->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf");
c2->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf");
c3->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf");
c4->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/GrowthCurves.pdf)");
	
	
//Time series plot
gStyle->SetTextFont(85);
gStyle->SetPadTickX(0); 
gStyle->SetPadTickY(0);

TCanvas *c00 = new TCanvas("c00","c00",4000,2000);
//First subplot
TMultiGraph*TimeSeries0 = new TMultiGraph();
TLegend *leg00 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries0->GetXaxis()->SetTitleOffset(1.2);
TimeSeries0->GetYaxis()->SetTitleOffset(0.7);
TimeSeries0->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries0->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries0->Add(TimeSeries[4]);
TimeSeries0->Add(TimeSeries[6]);
leg00->AddEntry(TimeSeries[4],"e- 10-30 MeV","l");
leg00->AddEntry(TimeSeries[6],"e+ 10-30 MeV","l");
c00->cd();
TimeSeries0->GetYaxis()->SetRangeUser(0,100);
TimeSeries0->GetXaxis()->SetLimits(0,133);
TimeSeries0->Draw("AP");		
leg00->Draw("same");

	
//Second subplot
TCanvas *c11 = new TCanvas("c11","c11",4000,2000);
TMultiGraph*TimeSeries1 = new TMultiGraph();
TLegend *leg11 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries1->GetXaxis()->SetTitleOffset(1.2);
TimeSeries1->GetYaxis()->SetTitleOffset(0.7);
TimeSeries1->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries1->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries1->Add(TimeSeries[3]);
TimeSeries1->Add(TimeSeries[5]);
leg11->AddEntry(TimeSeries[3],"e- 20-40 MeV","l");
leg11->AddEntry(TimeSeries[5],"e+ 20-40 MeV","l");
c11->cd();
TimeSeries1->GetYaxis()->SetRangeUser(0,100);
TimeSeries1->GetXaxis()->SetLimits(0,133);
TimeSeries1->Draw("ALP");		
leg11->Draw("same");
//Third subplot
TCanvas *c22 = new TCanvas("c22","c22",4000,2000);
TMultiGraph*TimeSeries2 = new TMultiGraph();
TLegend *leg22 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries2->GetXaxis()->SetTitleOffset(1.2);
TimeSeries2->GetYaxis()->SetTitleOffset(0.7);
TimeSeries2->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries2->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries2->Add(TimeSeries[2]);
TimeSeries2->Add(TimeSeries[6]);
leg22->AddEntry(TimeSeries[2],"e- 40-80 MeV","l");
leg22->AddEntry(TimeSeries[6],"e+ 40-80 MeV","l");
c22->cd();
TimeSeries2->GetYaxis()->SetRangeUser(0,100);
TimeSeries2->GetXaxis()->SetLimits(0,133);
TimeSeries2->Draw("ALP");		
leg22->Draw("same");
//Fourth subplot
TCanvas *c33 = new TCanvas("c33","c33",4000,2000);
TMultiGraph*TimeSeries3 = new TMultiGraph();
TLegend *leg33 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries3->GetXaxis()->SetTitleOffset(1.2);
TimeSeries3->GetYaxis()->SetTitleOffset(0.7);
TimeSeries3->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries3->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries3->Add(TimeSeries[1]);
TimeSeries3->Add(TimeSeries[7]);
leg33->AddEntry(TimeSeries[1],"e- 80-150 MeV","l");
leg33->AddEntry(TimeSeries[7],"e+ 80-150 MeV","l");
c33->cd();
TimeSeries3->GetYaxis()->SetRangeUser(0,100);
TimeSeries3->GetXaxis()->SetLimits(0,133);
TimeSeries3->Draw("ALP");		
leg33->Draw("same");
//Fifth subplot
TCanvas *c44 = new TCanvas("c44","c44",4000,2000);
TMultiGraph*TimeSeries4 = new TMultiGraph();
TLegend *leg44 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries4->GetXaxis()->SetTitleOffset(1.2);
TimeSeries4->GetYaxis()->SetTitleOffset(0.7);
TimeSeries4->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries4->GetYaxis()->SetTitle("Count rate/20mn");
TimeSeries4->Add(TimeSeries[0]);
TimeSeries4->Add(TimeSeries[8]);
leg44->AddEntry(TimeSeries[0],"e- 150-300 MeV","l");
leg44->AddEntry(TimeSeries[8],"e+ 150-300 MeV","l");
c44->cd();
TimeSeries4->GetYaxis()->SetRangeUser(0,100);
TimeSeries4->GetXaxis()->SetLimits(0,133);
TimeSeries4->Draw("ALP");		
leg44->Draw("same");

//c00->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeries.pdf(");
c11->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf(");
c22->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf");
c33->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf");
c44->Print("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/TimeSeriesHours.pdf)");

			
//Pressure vs time
TCanvas *CanPres = new TCanvas("CanPres","CanPress",4000,2000);
CanPres->cd();
PressB2vsTime->GetYaxis()->SetRangeUser(2.0,4.5);
PressB2vsTime->GetXaxis()->SetLimits(0,133);
PressB2vsTime->GetYaxis()->SetTitle("B2 calibrated pressure (g.cm^{-2})");
PressB2vsTime->GetXaxis()->SetTitle("Time since launch (in hours)");
PressB2vsTime->Draw("AP");
CanPres->SaveAs("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly_Calibrated/PressVsTime.eps");
	
//Write histograms in ROOT file
	
fileout->cd();	
/*
c0->Write();
GrowthCurves0->Write();
c1->Write();
GrowthCurves1->Write();
c2->Write();
GrowthCurves2->Write();
c3->Write();
GrowthCurves3->Write();
c4->Write();
GrowthCurves4->Write();

c11->Write();
TimeSeries1->Write();
c22->Write();
TimeSeries2->Write();
c33->Write();
TimeSeries3->Write();
c44->Write();
TimeSeries4->Write();

for(int k=0;k<NPBins;k++) {
	for (int l=0; l<NTimeBin;l++) {
		P0Bin[k][l]->Write();
		}
	}
for (int l=0; l<NTimeBin;l++) 
	{
  B2[l]->Write();
	}
	
cout << "now closing all files" << endl;
fileout->Close();
fileRaw->Close();
fileCal->Close(); 	
	*/
} //end function	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
