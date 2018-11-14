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
//Make time steps of 30 minutes
double binsperhour = 2;
const int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
const int NPBins = 10;
double PBins[NPBins] = {-300,-150,-80,-40,-20,20,40,80,150,300};
int NFiles = 1;


void TimeSeries(string RecoID) {
cout << "There are " << NTimeBin << " time bins " << endl;
	
TFile *fileRaw;
TFile *fileCal;
string ID = RecoID;	

//Create histograms
TH1F***P0BinDay = new TH1F**[NPBins];
TH1F***P0BinNight = new TH1F**[NPBins];
TH1F**B2 = new TH1F*[NTimeBin];

//Calibrated pressure from barometer 2
for(int k=0; k<NTimeBin;k++) 
	{
		B2[k] = new TH1F(Form("B2 h%d",k),Form("B2 h%d",k), 10000, 0, 800);
	}
//Momentum histograms
for(int i=0; i<NPBins; i++) 
{
		P0BinDay[i] = new TH1F*[NTimeBin];
		P0BinNight[i] = new TH1F*[NTimeBin];
	
		for(int j=0; j<NTimeBin; j++) 
		 {
			P0BinDay[i][j] = new TH1F(Form("signed p0 h%d%d Day",i,j), Form("signed p0 h%d%d Day",i,j), 1,0,1);
			P0BinNight[i][j] = new TH1F(Form("signed p0 h%d%d Night",i,j), Form("signed p0 h%d%d Night",i,j), 1,0,1);

		  
		}
	}
	
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
	
///////////////////////////////////////////////////////////
//////////////////DAYTIME DATA/////////////////////////////
///////////////////////////////////////////////////////////
	
for(int i=0; i<DTentries; i++) 
{
	chainDT->GetEntry(i);

    if(i%10000==0) cout << "Event " << i << endl;
//Trigger information and layers information
	uint8_t Ti=(uint8_t)dt->get_Ti();
    int  nhits = dt->get_Nhits();
	int NL= dt->get_NLayers();
    bool T1=dt->get_T1();
	bool T2=dt->get_T2();
	bool T3=dt->get_T3();
	bool T4=dt->get_T4();
	bool G=dt->get_guard();
//Time of event from CT1
	double yCT1 = (double)dt->get_yCT1();
	double mCT1 = (double)dt->get_mCT1();
	double dCT1 = (double)dt->get_dCT1();
	double hCT1= (double)dt->get_hCT1();
	double miCT1 = (double)dt->get_miCT1();
	double sCT1 = (double)dt->get_sCT1();
	float tmptimeDay = hCT1 + (miCT1/60) + (sCT1/3600);
	tmptimeDay = tmptimeDay-Launch;
	tmptimeDay = tmptimeDay + 24*(dCT1 - Launchday);
//PHA signals (calibrated)
    double T1PHA = dt->get_EneT1().at(0);
	double T2PHA = dt->get_EneT2().at(0);
	double T3PHA = dt->get_EneT3().at(0);
	double T4PHA = dt->get_EneT4().at(0);
	double GPHA = dt->get_Eneg().at(0);

//Barometer data (calibrated)
	float B2Pres = dt->get_PressB2();
//Convert pressure to g/cm^-2
	B2Pres = B2Pres* 1.3595;
//Reconstructed energy 
	float pPR=1000*fabs(dt->get_p0PR());   //in MeV
	double deflecPR = dt->get_deflecPR();
	float P0sign = pPR*TMath::Sign(1,deflecPR);
	double chiNBPR = dt->get_chi2NBPR();
	double chiBPR = dt->get_chi2BPR();	
//Determine the correct time and energy bin
	int tbin = tmptimeDay*binsperhour;
	//cout << " DayTime t = " << tmptimeDay << "h,  tbin = " << tbin << endl;

	int pbin;
    //	cout << "Assign time bin " << tbin <<endl;

	for(int k=0;k<NPBins;k++) 
		{
		if((PBins[k] <= P0sign) & (P0sign < PBins[k+1])) 
			{
			pbin = k;
		  //  cout << PBins[k] << " =< " << P0sign << " < " << PBins[k+1] << endl;
		 }
		}
	if(tmptimeDay>0) 
		{
	B2[tbin]->Fill(B2Pres);
    if((T1PHA>=65) & (T1PHA<=85) & (T2PHA>180) & (T3PHA>=65) & (T3PHA<=95) & (P0sign >= -300) & (P0sign <= 300) & (!G) & (NL==7) & (nhits >=7) && (nhits <15))
	   {
		P0BinDay[pbin][tbin]->Fill(P0sign);

		}	   
	}		
}     //end i, loop on events

	
///////////////////////////////////////////////////////////
//////////////////NIGHTTIME DATA/////////////////////////////
///////////////////////////////////////////////////////////
	
for(int i=0; i<NTentries; i++) 
{
	chainNT->GetEntry(i);

    if(i%10000==0) cout << "Event " << i << endl;
//Trigger information and layers information
	uint8_t Ti=(uint8_t)nt->get_Ti();
    int  nhits = nt->get_Nhits();
	int NL= nt->get_NLayers();
    bool T1=nt->get_T1();
	bool T2=nt->get_T2();
	bool T3=nt->get_T3();
	bool T4=nt->get_T4();
	bool G=nt->get_guard();
//Time of event from CT1
	double yNCT1 = (double)nt->get_yCT1();
	double mNCT1 = (double)nt->get_mCT1();
	double dNCT1 = (double)nt->get_dCT1();
	double hNCT1= (double)nt->get_hCT1();
	double miNCT1 = (double)nt->get_miCT1();
	double sNCT1 = (double)nt->get_sCT1();
	float tmptimeNight = hNCT1 + (miNCT1/60) + (sNCT1/3600);
	tmptimeNight = tmptimeNight-Launch;
	tmptimeNight = tmptimeNight + 24*(dNCT1 - Launchday);
//PHA signals (calibrated)
    double T1PHA = nt->get_EneT1().at(0);
	double T2PHA = nt->get_EneT2().at(0);
	double T3PHA = nt->get_EneT3().at(0);
	double T4PHA = nt->get_EneT4().at(0);
	double GPHA = nt->get_Eneg().at(0);

//Barometer data (calibrated)
	float B2Pres = nt->get_PressB2();
//Convert pressure to g/cm^-2
	B2Pres = B2Pres* 1.3595;
//Reconstructed energy 
	float pPR=1000*fabs(nt->get_p0PR());   //in MeV
	double deflecPR = nt->get_deflecPR();
	float P0sign = pPR*TMath::Sign(1,deflecPR);
	double chiNBPR = nt->get_chi2NBPR();
	double chiBPR = nt->get_chi2BPR();	
//Determine the correct time and energy bin
	int tbin = tmptimeNight*binsperhour;
	int pbin;
    //	cout << "Assign time bin " << tbin <<endl;

	for(int k=0;k<NPBins;k++) 
		{
		if((PBins[k] <= P0sign) & (P0sign < PBins[k+1])) 
			{
			pbin = k;
	//	    cout << PBins[k] << " =< " << P0sign << " < " << PBins[k+1] << endl;
		 }
		}
	if(tmptimeNight>0) 
		{
	B2[tbin]->Fill(B2Pres);
    if((T1PHA>=65) & (T1PHA<=85) & (T2PHA>180) & (T3PHA>=65) & (T3PHA<=95) & (P0sign >= -300) & (P0sign <= 300) & (!G) & (NL==7) & (nhits >=7) && (nhits <15))
	   {
		P0BinNight[pbin][tbin]->Fill(P0sign);
		}	   
	}		
}     //end i, loop on events

	

//Make growth curves
int ngraphs = (NPBins/2)-1;
TGraph**TimeSeriesDay = new TGraph*[NPBins];
TGraph**TimeSeriesNight = new TGraph*[NPBins];
TGraph**CountRateDay = new TGraph*[NPBins];
TGraph**CountRateNight = new TGraph*[NPBins];




for (int j=0; j<NPBins;j++) 
	{

		TimeSeriesDay[j] = new TGraph();
	    TimeSeriesDay[j]->SetMarkerStyle(kDot);
		TimeSeriesDay[j]->SetMarkerSize(0.8);
	
		TimeSeriesNight[j] = new TGraph();
	    TimeSeriesNight[j]->SetMarkerStyle(kDot);
		TimeSeriesNight[j]->SetMarkerSize(0.8);
		
	    CountRateDay[j] = new TGraph();
	    CountRateDay[j]->SetMarkerStyle(25);
	//	CountRateDay[j]->SetMarkerSize(0.8);
	
		CountRateNight[j] = new TGraph();
	    CountRateNight[j]->SetMarkerStyle(kOpenCircle);
	//	CountRateNight[j]->SetMarkerSize(0.8);

//electrons
		if (j<4)  { 
		TimeSeriesDay[j]->SetMarkerColor(kBlue);
		TimeSeriesDay[j]->SetLineColor(kBlue);
	    TimeSeriesNight[j]->SetMarkerColor(kGreen);
		TimeSeriesNight[j]->SetLineColor(kGreen);
		CountRateDay[j]->SetMarkerColor(kBlue);
		CountRateNight[j]->SetMarkerColor(kBlue);

//positrons
			}
		if (j>4)  { 
		TimeSeriesDay[j]->SetMarkerColor(kRed);
		TimeSeriesDay[j]->SetLineColor(kRed);
		TimeSeriesNight[j]->SetMarkerColor(kBlack);
		TimeSeriesNight[j]->SetLineColor(kBlack);
		CountRateDay[j]->SetMarkerColor(kRed);
		CountRateNight[j]->SetMarkerColor(kRed);
			}
		int npPosDay = 1;
		int npEleDay = 1;
		int npPosNight = 1;
		int npEleNight = 1;
		for (int k=0;k<NTimeBin;k++) 
		{
			double RatePos;
			double RateEle;
			double time = Tstart + (k+1)/(binsperhour);
			double B2mean = B2[k]->GetMean();

			//cout << "daytime t = " << time << "h" << endl;
			if (j<4) 
				{
				RateEle = P0BinDay[j][k]->GetEntries(); //count rate for 30 min.
				TimeSeriesDay[j]->SetPoint(npEleDay,time,RateEle);
				CountRateDay[j]->SetPoint(npEleDay,B2mean,RateEle);
				npEleDay++;
				
				}
		   if(j>4)
		   	{
			   RatePos = P0BinDay[j][k]->GetEntries(); //count rate for 30 min.
			   TimeSeriesDay[j]->SetPoint(npPosDay,time,RatePos);
			   CountRateDay[j]->SetPoint(npPosDay,B2mean,RatePos);
			   npPosDay++;
			   	
			   }
		   	}
		for (int k=0;k<NTimeBin;k++) 
		{
			double RatePos;
			double RateEle;
			double time = Tstart + (k+1)/(binsperhour);
			double B2mean = B2[k]->GetMean();
			//cout << "nightime t = " << time << "h" << endl;
			
			if (j<4) 
				{
				RateEle = P0BinNight[j][k]->GetEntries(); //count rate for 30 min.
				TimeSeriesNight[j]->SetPoint(npEleNight,time,RateEle);
			    CountRateNight[j]->SetPoint(npEleNight,B2mean,RateEle);			
				npEleNight++;
				
				}
		   if(j>4)
		   	{
			   RatePos = P0BinNight[j][k]->GetEntries(); //count rate for 30 min.
			   TimeSeriesNight[j]->SetPoint(npPosNight,time,RatePos);
			   CountRateNight[j]->SetPoint(npPosNight,B2mean,RatePos);			
			   npPosNight++;
			   
			   }
		   	}
		}
    
	
						 
//Time series plot
gStyle->SetTextFont(85);
gStyle->SetPadTickX(0); 
gStyle->SetPadTickY(0);


TCanvas *c11 = new TCanvas("c11","c11",4000,2000);
TMultiGraph*TimeSeries1 = new TMultiGraph();
TLegend *leg11 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries1->GetXaxis()->SetTitleOffset(1.2);
TimeSeries1->GetYaxis()->SetTitleOffset(0.7);
TimeSeries1->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries1->GetYaxis()->SetTitle("Count rate/30mn");
TimeSeries1->Add(TimeSeriesDay[3]);
TimeSeries1->Add(TimeSeriesDay[5]);
leg11->AddEntry(TimeSeriesDay[3],"e- 20-40 MeV Day","l");
leg11->AddEntry(TimeSeriesDay[5],"e+ 20-40 MeV Day","l");
TimeSeries1->Add(TimeSeriesNight[3]);
TimeSeries1->Add(TimeSeriesNight[5]);
leg11->AddEntry(TimeSeriesNight[3],"e- 20-40 MeV Night","l");
leg11->AddEntry(TimeSeriesNight[5],"e+ 20-40 MeV Night","l");
c11->cd();
TimeSeries1->GetYaxis()->SetRangeUser(0,100);
TimeSeries1->GetXaxis()->SetLimits(0,135);
TimeSeries1->Draw("ALP");		
leg11->Draw("same");
//
TCanvas *c22 = new TCanvas("c22","c22",4000,2000);
TMultiGraph*TimeSeries2 = new TMultiGraph();
TLegend *leg22 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries2->GetXaxis()->SetTitleOffset(1.2);
TimeSeries2->GetYaxis()->SetTitleOffset(0.7);
TimeSeries2->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries2->GetYaxis()->SetTitle("Count rate/30mn");
TimeSeries2->Add(TimeSeriesDay[2]);
TimeSeries2->Add(TimeSeriesDay[6]);
leg22->AddEntry(TimeSeriesDay[2],"e- 40-80 MeV Day","l");
leg22->AddEntry(TimeSeriesDay[6],"e+ 40-80 MeV Day","l");
TimeSeries2->Add(TimeSeriesNight[2]);
TimeSeries2->Add(TimeSeriesNight[6]);
leg22->AddEntry(TimeSeriesNight[2],"e- 40-80 MeV Night","l");
leg22->AddEntry(TimeSeriesNight[6],"e+ 40-80 MeV Night","l");
c22->cd();
TimeSeries2->GetYaxis()->SetRangeUser(0,100);
TimeSeries2->GetXaxis()->SetLimits(0,135);
TimeSeries2->Draw("ALP");		
leg22->Draw("same");
//
TCanvas *c33 = new TCanvas("c33","c33",4000,2000);
TMultiGraph*TimeSeries3 = new TMultiGraph();
TLegend *leg33 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries3->GetXaxis()->SetTitleOffset(1.2);
TimeSeries3->GetYaxis()->SetTitleOffset(0.7);
TimeSeries3->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries3->GetYaxis()->SetTitle("Count rate/30mn");
TimeSeries3->Add(TimeSeriesDay[1]);
TimeSeries3->Add(TimeSeriesDay[7]);
leg33->AddEntry(TimeSeriesDay[1],"e- 80-150 MeV Day","l");
leg33->AddEntry(TimeSeriesDay[7],"e+ 80-150 MeV Day","l");
TimeSeries3->Add(TimeSeriesNight[1]);
TimeSeries3->Add(TimeSeriesNight[7]);
leg33->AddEntry(TimeSeriesNight[1],"e- 80-150 MeV Night","l");
leg33->AddEntry(TimeSeriesNight[7],"e+ 80-150 MeV Night","l");
c33->cd();
TimeSeries3->GetYaxis()->SetRangeUser(0,100);
TimeSeries3->GetXaxis()->SetLimits(0,135);
TimeSeries3->Draw("ALP");		
leg33->Draw("same");
//
TCanvas *c44 = new TCanvas("c44","c44",4000,2000);
TMultiGraph*TimeSeries4 = new TMultiGraph();
TLegend *leg44 = new TLegend(0.6,0.85,0.9,0.9);
TimeSeries4->GetXaxis()->SetTitleOffset(1.2);
TimeSeries4->GetYaxis()->SetTitleOffset(0.7);
TimeSeries4->GetXaxis()->SetTitle("Time from launch (in hours)");
TimeSeries4->GetYaxis()->SetTitle("Count rate/30mn");
TimeSeries4->Add(TimeSeriesDay[0]);
TimeSeries4->Add(TimeSeriesDay[8]);
leg44->AddEntry(TimeSeriesDay[0],"e- 150-300 MeV Day","l");
leg44->AddEntry(TimeSeriesDay[8],"e+ 150-300 MeV Day","l");
TimeSeries4->Add(TimeSeriesNight[0]);
TimeSeries4->Add(TimeSeriesNight[8]);
leg44->AddEntry(TimeSeriesNight[0],"e- 150-300 MeV Night","l");
leg44->AddEntry(TimeSeriesNight[8],"e+ 150-300 MeV Night","l");
c44->cd();
TimeSeries4->GetYaxis()->SetRangeUser(0,100);
TimeSeries4->GetXaxis()->SetLimits(0,135);
TimeSeries4->Draw("ALP");		
leg44->Draw("same");


c11->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/TimeSeries.pdf(");
c22->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/TimeSeries.pdf");
c33->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/TimeSeries.pdf");
c44->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/TimeSeries.pdf)");
	
//Growth Curves
TCanvas *c1 = new TCanvas("c1","c1",1100,2000);
TMultiGraph*GrowthCurve1 = new TMultiGraph();
GrowthCurve1->SetTitle("DayTime Growth Curves");
TLegend *leg1 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurve1->GetXaxis()->SetTitleOffset(1.2);
GrowthCurve1->GetYaxis()->SetTitleOffset(0.7);
GrowthCurve1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurve1->GetYaxis()->SetTitle("Count rate/30mn");
GrowthCurve1->Add(CountRateDay[3]);
GrowthCurve1->Add(CountRateDay[5]);
leg1->AddEntry(CountRateDay[3],"e- 20-40 MeV Day","p");
leg1->AddEntry(CountRateDay[5],"e+ 20-40 MeV Day","p");
/*
GrowthCurve1->Add(CountRateNight[3]);
GrowthCurve1->Add(CountRateNight[5]);
leg1->AddEntry(CountRateNight[3],"e- 20-40 MeV Night","p");
leg1->AddEntry(CountRateNight[5],"e+ 20-40 MeV Night","p");
*/
c1->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurve1->SetMinimum(1);
GrowthCurve1->Draw("AP");		
leg1->Draw("same");
c1->Modified(); 
c1->Update();
GrowthCurve1->GetXaxis()->SetRangeUser(1,1000);
//
TCanvas *c2 = new TCanvas("c2","c2",1100,2000);
TMultiGraph*GrowthCurve2 = new TMultiGraph();
GrowthCurve2->SetTitle("DayTime Growth Curves");	
TLegend *leg2 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurve2->GetXaxis()->SetTitleOffset(1.2);
GrowthCurve2->GetYaxis()->SetTitleOffset(0.7);
GrowthCurve2->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurve2->GetYaxis()->SetTitle("Count rate/30mn");
GrowthCurve2->Add(CountRateDay[2]);
GrowthCurve2->Add(CountRateDay[6]);
leg2->AddEntry(CountRateDay[2],"e- 40-80 MeV Day","p");
leg2->AddEntry(CountRateDay[6],"e+ 40-80 MeV Day","p");
/*
GrowthCurve2->Add(CountRateNight[2]);
GrowthCurve2->Add(CountRateNight[6]);
leg2->AddEntry(CountRateNight[2],"e- 40-80 MeV Night","p");
leg2->AddEntry(CountRateNight[6],"e+ 40-80 MeV Night","p");
*/
c2->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurve2->SetMinimum(1);
GrowthCurve2->Draw("AP");		
leg2->Draw("same");
c1->Modified(); 
c1->Update();
GrowthCurve2->GetXaxis()->SetRangeUser(1,1000);
//
TCanvas *c3 = new TCanvas("c3","c3",1100,2000);
TMultiGraph*GrowthCurve3 = new TMultiGraph();
GrowthCurve3->SetTitle("DayTime Growth Curves");
TLegend *leg3 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurve3->GetXaxis()->SetTitleOffset(1.2);
GrowthCurve3->GetYaxis()->SetTitleOffset(0.7);
GrowthCurve3->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurve3->GetYaxis()->SetTitle("Count rate/30mn");
GrowthCurve3->Add(CountRateDay[1]);
GrowthCurve3->Add(CountRateDay[7]);
leg3->AddEntry(CountRateDay[1],"e- 80-150 MeV Day","p");
leg3->AddEntry(CountRateDay[7],"e+ 80-150 MeV Day","p");
/*
GrowthCurve3->Add(CountRateNight[1]);
GrowthCurve3->Add(CountRateNight[7]);
leg3->AddEntry(CountRateNight[1],"e- 80-150 MeV Night","p");
leg3->AddEntry(CountRateNight[7],"e+ 80-150 MeV Night","p");
*/
c3->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurve3->SetMinimum(1);
GrowthCurve3->Draw("AP");		
leg3->Draw("same");
c1->Modified(); 
c1->Update();
GrowthCurve3->GetXaxis()->SetRangeUser(1,1000);
//
TCanvas *c4 = new TCanvas("c4","c4",1100,2000);
TMultiGraph*GrowthCurve4 = new TMultiGraph();
TLegend *leg4 = new TLegend(0.6,0.85,0.9,0.9);
GrowthCurve4->SetTitle("DayTime Growth Curves");
GrowthCurve4->GetXaxis()->SetTitleOffset(1.2);
GrowthCurve4->GetYaxis()->SetTitleOffset(0.7);
GrowthCurve4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
GrowthCurve4->GetYaxis()->SetTitle("Count rate/30mn");
GrowthCurve4->Add(CountRateDay[0]);
GrowthCurve4->Add(CountRateDay[8]);
leg4->AddEntry(CountRateDay[0],"e- 150-300 MeV Day","p");
leg4->AddEntry(CountRateDay[8],"e+ 150-300 MeV Day","p");
/*
GrowthCurve4->Add(CountRateNight[0]);
GrowthCurve4->Add(CountRateNight[8]);
leg4->AddEntry(CountRateNight[0],"e- 150-300 MeV Night","p");
leg4->AddEntry(CountRateNight[8],"e+ 150-300 MeV Night","p");
*/
c4->cd();
gPad->SetLogy();
gPad->SetLogx();
gStyle->SetPadTickX(1);  
gStyle->SetPadTickY(1);
GrowthCurve4->SetMinimum(1);
GrowthCurve4->Draw("AP");		
leg4->Draw("same");
c1->Modified(); 
c1->Update();
GrowthCurve4->GetXaxis()->SetRangeUser(1,1000);


c1->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/GrowthCurveTimeCut.pdf(");
c2->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/GrowthCurveTimeCut.pdf");
c3->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/GrowthCurveTimeCut.pdf");
c4->Print("/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/GrowthCurveTimeCut.pdf)");
	
} //end function	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	