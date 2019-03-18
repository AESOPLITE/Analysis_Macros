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

void LoadDSparameters(string filename, float**DScoeff) {
// Read the file filename and load the 7 coefficients 
//of the Daniel and Stephens growth curve fit functions for 7 energy bins

//Number of lines in the file
int n=7; 
//Prefix chracters at beginning of line
string pre[7]={"20MEV","30MEV", "44.8MEV", "70MEV", "111MEV", "177.5MEV","279MEV"};
	
//Create stream from filename 
ifstream file;
file.open(filename, ios_base::in); // open file
int j=0;
for(string line; getline(file, line); )   //read stream line by line
   {
//    cout << line << endl;
    istringstream in(line);      //make a stream for the line itself
    string prefix;
    in >> prefix; //and read the first whitespace-separated token
    float coeff0;
    in >> coeff0;
    float coeff1;
    in >> coeff1;
    float coeff2;
    in >> coeff2;
    float coeff3;
    in >> coeff3;
    float coeff4;
    in >> coeff4;
    float coeff5;
    in >> coeff5;
    float coeff6;
    in >> coeff6;

	
    //check prefix to load the appropriate region variable	
	for(int i=0;i<n;i++) //what line we are on
      {
       if(prefix.compare(pre[i]) == 0)
        {
		   DScoeff[i][0]=coeff0;
		   DScoeff[i][1]=coeff1;
		   DScoeff[i][2]=coeff2;
		   DScoeff[i][3]=coeff3;
		   DScoeff[i][4]=coeff4;
		   DScoeff[i][5]=coeff5;
		   DScoeff[i][6]=coeff6;
		   j++;
	   }
	}
 }
	  if (j!=n) cout << "Error when loading the DC coeff parameters: Wrong number of parameters in the file " << filename << endl;
}
	
	
//From Paul's FORTRAN code	
double DanAndSteph(Double_t *x, Double_t *coeff) {
	double z8=0.0;	
	double x8=x[0];	
	for(int i=0; i<7;i++) 
	{
		int k = 5 + 1 - i;
		//cout << " k = " << k << ", coeff[k] = " << coeff[k] << endl;
		z8 = z8*TMath::Log(x8) + coeff[k];
	}
	//cout << "z8 = " << TMath::Exp(z8) << endl;
	return coeff[7]*TMath::Exp(z8);	
	
}


void GrowthCurvesTimeBin(string RecoID, double cutT2, double cutT3,bool ts) {

	double Tstart = 0;
	double Tend = 145;	//total time of flight in hours)
	double Launch = 22 + (7/60);
	double Launchday = 15;
	//Make time steps of 30mn
	double binsperhour = 2;
	int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
	const int NPBins = 12;
	double PBins[NPBins] = {-300,-150,-80,-60,-30,-10,10,30,60,80,150,300};
	double GeoFactor[5] = {2.68,10.43,11.55,12.02,14.64};
	//for the ascent all data is in file 000
	//for the final leg of the flight (no cut-off), data starts at file 15
	int StartFile = 0;
	int EndFile = 23;
	double TCut[7] = {19.8,30,42,56,65,82,90};
	int mcolor[7]={1,2,3,4,800+1,6,7};


	string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2/";
	string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves";
	string CutConfig = "AllFlightT2Above70T3Less110";
	string ID = RecoID;	

	//string CutConfig="testfit";
	//txt output file for positrons
	ofstream totxtPos;
	totxtPos.open(Form("%s/TimeSeriesPos_%s_%s.txt",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	totxtPos << " Time in hours " << Form("  Pos CountRate %2.0f/mn",60/binsperhour) << endl;
	//txt output file for electrons
	ofstream totxtEle;
	totxtEle.open(Form("%s/TimeSeriesEle_%s_%s.txt",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	totxtEle << " Time in hours " << Form("  Ele CountRate %2.0f/mn",60/binsperhour) << endl;

	//Axis limit for plots
	float ZMax=100000;		//growth curve y-axis
	float ZMin=100;			//growth curve y-axis
	float PMin=1;			//growth curve x-axis
	float PMax=1000;			//growth curve x-axis
	float TMin=0;			//time serie x-axis
	float TMax=135;			//time serie x-axis	
	float CountMax = 300;      //time serie y-axis
	
	//Initialize array
	string DSConfig[7] = {"20 MeV","30 MeV", "44.8 MeV", "70 MeV", "111 MeV", "177.5 MeV","279 MeV"};
	float**DScoeff=new float*[7];
	for(int i=0; i<7;i++) 
	 {
	 DScoeff[i] =  new float[7];
	 for(int j=0; j<7;j++)			//7 coefficients for each array 
	  {
		 DScoeff[i][j]=0.;
	 }
    }
	 string DSparamfile="./DScoefficients.dat"; 
	 LoadDSparameters(DSparamfile,DScoeff);	

	bool boolts = ts;		//boolean to draw, or not timeseries	
	double T2Cut = cutT2;	//value of cut on T2
	double T3Cut = cutT3;	//value of cut on T3
	TFile *fileRaw;
	TFile *fileCal;

	//Create histograms
	TH1F***B2 = new TH1F**[NPBins];
	TH1F***P0Bin = new TH1F**[NPBins];
    float***B2vec = new float**[NPBins];
    int**B2vecindex = new int*[NPBins];

	//Momentum histograms and pressure
	for(int i=0; i<NPBins; i++) 
	{
			P0Bin[i] = new TH1F*[NTimeBin];
			B2[i] = new TH1F*[NTimeBin];
		    B2vec[i]=new float*[NTimeBin];
	        B2vecindex[i]=new int[NTimeBin];
			for(int j=0; j<NTimeBin; j++) 
			 {
				P0Bin[i][j] = new TH1F(Form("signed p0 h%d%d",i,j), Form("signed p0 h%d%d",i,j), 1,0,1);
                B2[i][j] = new TH1F(Form("B2 h%d%d",i,j),Form("B2 h%d%d",i,j), 10000, 0, 800);
				B2vec[i][j]=new float[1000];
				B2vecindex[i][j]=0;
                for(int k =0;k<1000;k++) 
				{
					B2vec[i][j][k]=0.;
				}	
			   }
		}

	double LaunchDay = 15 + (Launch*0.041667);

	//ROOT fileout
	 TFile *fileout = new TFile(Form("%s/GrowthCurves_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
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
	//	if(i%1000==0) cout << "Event " << i << endl;
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
	//Time of event from PHA
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
		float tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
	 // cout << "tmptimeCT1 = " << tmptime << ", tmptimePHA = " << tmptimePHA << endl; 
		
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
		//cout << "B2Pres = " << B2Pres << "g.cm^-2 " << endl;
	//Reconstructed energy 
		float pPR=1000*fabs(de->get_p0PR());   //in MeV
		double deflecPR = de->get_deflecPR();
		float preco=de->get_p0reco();
	    float P0sign;
		if (strcmp(ID.c_str(),"PRonly")==0)  P0sign = pPR*TMath::Sign(1,deflecPR);
		else P0sign = preco*TMath::Sign(1,deflecPR);	
		double chiNBPR = de->get_chi2NBPR();
		double chiBPR = de->get_chi2BPR();	
	//Determine the correct time and energy bin
		int tbin = tmptimePHA*binsperhour;
		int pbin;

		for(int k=0;k<NPBins;k++) 
			{
			if((PBins[k] <= P0sign) & (P0sign < PBins[k+1])) 
				{
				pbin = k;
			 }
			}
		if(tmptimePHA>0) 
			{
	//Apply cuts to only get clean events
		if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
		   {
			//cout << "Assign bin number energy bin " << pbin << ", time bin " << tbin <<endl;
		// cout << "Fill event P0sign = " << P0sign << "MeV, at time  =  " << tmptimePHA << " hours" <<endl;
			P0Bin[pbin][tbin]->Fill(P0sign);
		    B2[pbin][tbin]->Fill(B2Pres);
			int index = B2vecindex[pbin][tbin];
			B2vec[pbin][tbin][index]=B2Pres;
			B2vecindex[pbin][tbin]++;
			}	   
		}
	}     //end i, loop on events
		cout << "finished reading a file" << endl;
	}		//end j, loop on files

	
		  for (int j=0; j<NPBins-1;j++) 
	     {
		   
		  for (int k=0;k<NTimeBin;k++) 
		 {
		   double mean =  B2[j][k]->GetMean();
		   double meanerr =  B2[j][k]->GetMeanError();
		   double nentries =  B2vecindex[j][k];
		   float B2min = 2000;
		   float B2max = 0;
 	       for(int l=0;l<B2vecindex[j][k];l++) 
		   {
			   double val = B2vec[j][k][l];
			   if (val > B2max) B2max = val;
			   if (val < B2min) B2min = val;
		   }
    B2[j][k]->SetMinimum(B2min);
    B2[j][k]->SetMaximum(B2max);
		  }
		 }
	
	//Make growth curves
	int ngraphs = (NPBins/2)-1;
	TGraphAsymmErrors**CountRate = new TGraphAsymmErrors*[NPBins];
    TGraphAsymmErrors**FluxRate = new TGraphAsymmErrors*[NPBins];
	TGraphAsymmErrors**TimeSeries = new TGraphAsymmErrors*[NPBins];
	TGraph *PressB2vsTime = new TGraph();
	PressB2vsTime->SetMarkerStyle(kOpenCircle);
	//PressB2vsTime->SetMarkerSize(0.1);
	PressB2vsTime->SetMarkerColor(kBlue);	

	for (int j=0; j<NPBins-1;j++) 
	{
		CountRate[j] = new TGraphAsymmErrors();
	    CountRate[j]->SetMarkerStyle(kOpenCircle);
		CountRate[j]->SetMarkerSize(0.8);
		FluxRate[j] = new TGraphAsymmErrors();
	    FluxRate[j]->SetMarkerStyle(kOpenCircle);
		FluxRate[j]->SetMarkerSize(0.8);
		TimeSeries[j] = new TGraphAsymmErrors();
	    TimeSeries[j]->SetMarkerStyle(kDot);
		TimeSeries[j]->SetMarkerSize(0.8);
		if (j<5)  { 
		CountRate[j]->SetMarkerColor(kBlue);
		FluxRate[j]->SetMarkerColor(kBlue);
		FluxRate[j]->SetLineColor(kBlue);
		TimeSeries[j]->SetMarkerColor(kBlue);
		TimeSeries[j]->SetLineColor(kBlue);
			}
		if (j>5)  { 
		CountRate[j]->SetMarkerColor(kRed);
		FluxRate[j]->SetMarkerColor(kRed);
		FluxRate[j]->SetLineColor(kRed);
		TimeSeries[j]->SetMarkerColor(kRed);
		TimeSeries[j]->SetLineColor(kRed);
			}
		int npPos = 1;
		int npEle = 1;
		int np2 = 1;


		for (int k=0;k<NTimeBin;k++) 
		{
			double RatePos,RateEle;
			double FluxPos,FluxEle;
			double RateEleErr,RatePosErr;
			double FluxPosErr,FluxEleErr;			
			double time = Tstart + (k+1)/(binsperhour);
			double timeWidth = (60/binsperhour)*60;				//width of time bin in seconds
			double MomWidth;
		    double B2mean = B2[j][k]->GetMean();
		    double B2max = B2[j][k]->GetMaximum();
		    double B2min = B2[j][k]->GetMinimum();		    
//Time since launch expressed in units of days
			double timeInDays = LaunchDay + time*(0.041667);
			if (j<5) 
				{
				MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
				RateEle = P0Bin[j][k]->GetEntries(); 
				RateEleErr = TMath::Sqrt(P0Bin[j][k]->GetEntries());
				FluxEle = RateEle/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidth);
				FluxEleErr = RateEleErr/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidth);
				if(B2mean!=0) {
		             cout << "j = " << j <<", RateEle = " << RateEle << endl;
					CountRate[j]->SetPoint(npEle,B2mean,RateEle);
					FluxRate[j]->SetPoint(npEle,B2mean,FluxEle);
				    FluxRate[j]->SetPointError(npEle,0,0,FluxEleErr,FluxEleErr);
					}
				TimeSeries[j]->SetPoint(npEle,time,RateEle);
				if(j==3) totxtEle << time << "   " << RateEle << endl;
				npEle++;
				}
		   if(j>5)
		   	{
               if (j!=NPBins-1) MomWidth = TMath::Abs(PBins[j] - PBins[j+1]);
			   RatePos = P0Bin[j][k]->GetEntries(); 
			   FluxPos = RatePos/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidth);
			   RatePosErr = TMath::Sqrt(P0Bin[j][k]->GetEntries());
			   FluxPosErr = RatePosErr/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidth);
			   if(B2mean!=0) {
		             cout << "j = " << j <<", RatePos = " << RatePos << endl;
				   CountRate[j]->SetPoint(npPos,B2mean,RatePos);
				   FluxRate[j]->SetPoint(npPos,B2mean,FluxPos);
				  FluxRate[j]->SetPointError(npPos,0,0,FluxPosErr,FluxPosErr);

			   }
     		   TimeSeries[j]->SetPoint(npPos,time,RatePos);
			   if(j==7) totxtPos << time << "   " << RatePos << endl;
			   npPos++;
			  
		   	}
		}
    }
	
	//Define Daniel & Stephens functions
	TF1**DS = new TF1*[7];
	for (int k=0;k<7;k++) 
	{
		DS[k] = new TF1(Form("ds%d",k),DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DS[k]->FixParameter(l,DScoeff[k][l]);
		DS[k]->SetParameter(7,1);
		DS[k]->SetLineColor(mcolor[k]);
		DS[k]->SetLineStyle(6);

	}
	
	//////////////////////////////////////////////////
	/////////////////DRAW CANVASES////////////////////
	//////////////////////////////////////////////////
	
	gStyle->SetTextFont(85);
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1);	
	gStyle->SetTitleSize(0.04,"XY");

	//Draw All Daniel & Stephens functions
	TCanvas *c = new TCanvas("c","c",1000,1000);
	TMultiGraph*Functions = new TMultiGraph();
	CountRate[4]->SetMarkerColor(kWhite);
	Functions->Add(CountRate[4]);
	TLegend *leg = new TLegend(0.20,0.70,0.3,0.85);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	c->cd();
	Functions->Draw("A");
	for (int k=0;k<7;k++) 
	{
		DS[k]->Draw("same");
	    leg->AddEntry(DS[k],Form("D&S %s", DSConfig[k].c_str()),"l");
	}
	leg->Draw("same");
	Functions->GetXaxis()->SetTitleOffset(1.2);
	Functions->GetYaxis()->SetTitleOffset(1.3);
	Functions->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	Functions->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	Functions->GetYaxis()->CenterTitle();	
	Functions->GetXaxis()->SetLimits(PMin,2000);
	gPad->SetLogy();
	gPad->SetLogx();
	Functions->SetMinimum(ZMin);
	Functions->SetMaximum(100000);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 
	
	//First growth curve subplot
	TCanvas *c0 = new TCanvas("c0","c0",1000,1000);
	TMultiGraph*GrowthCurves0 = new TMultiGraph();
	TLegend *leg0 = new TLegend(0.55,0.70,0.9,0.85);
	leg0->SetBorderSize(0);
	leg0->SetFillStyle(0);
	CountRate[4]->SetMarkerColor(kBlue);
	GrowthCurves0->Add(FluxRate[4]);
	GrowthCurves0->Add(FluxRate[6]);
	leg0->AddEntry(FluxRate[4],"e- 10-30 MeV","p");
    leg0->AddEntry(FluxRate[6],"e+ 10-30 MeV","p");	
	leg0->AddEntry(DS[0],"Daniel & Stephens 20 MeV","l");
	leg0->AddEntry(DS[1],"Daniel & Stephens 30 MeV","l");	
	c0->cd();
	GrowthCurves0->Draw("AP");
	FluxRate[4]->Fit("ds0","R");
	//DS[0]->Draw("same");
	//DS[1]->Draw("same");
	leg0->Draw("same");
	GrowthCurves0->GetXaxis()->SetTitleOffset(1.2);
	GrowthCurves0->GetYaxis()->SetTitleOffset(1.3);
	GrowthCurves0->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves0->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves0->GetYaxis()->CenterTitle();
	//GrowthCurves0->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	GrowthCurves0->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves0->SetMinimum(ZMin);
	GrowthCurves0->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Second subplot
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TMultiGraph*GrowthCurves1 = new TMultiGraph();
	TLegend *leg1 = new TLegend(0.55,0.70,0.9,0.85);
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	GrowthCurves1->Add(FluxRate[3]);
	GrowthCurves1->Add(FluxRate[7]);
	leg1->AddEntry(FluxRate[3],"e- 30-60 MeV","p");
	leg1->AddEntry(FluxRate[7],"e+ 30-60 MeV","p");
	leg1->AddEntry(DS[2],"Daniel & Stephens 50 MeV","l");
	c1->cd();
	GrowthCurves1->Draw("AP");		
	FluxRate[3]->Fit("ds2","","",80,300);	
	DS[2]->DrawF1(1,1000,"same");
	leg1->Draw("same");
	GrowthCurves1->GetXaxis()->SetTitleOffset(1.2);
	GrowthCurves1->GetYaxis()->SetTitleOffset(1.3);
	GrowthCurves1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves1->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves1->GetYaxis()->CenterTitle();
	GrowthCurves1->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves1->SetMinimum(ZMin);
	GrowthCurves1->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Third subplot
	TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
	TMultiGraph*GrowthCurves2 = new TMultiGraph();
	TLegend *leg2 = new TLegend(0.55,0.70,0.9,0.85);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);

	GrowthCurves2->Add(FluxRate[2]);
	GrowthCurves2->Add(FluxRate[8]);
	leg2->AddEntry(FluxRate[2],"e- 60-80 MeV","p");
	leg2->AddEntry(FluxRate[8],"e+ 60-80 MeV","p");
	leg2->AddEntry(DS[3],"Daniel & Stephens 70 MeV","l");
	c2->cd();
	GrowthCurves2->Draw("AP");	
	DS[3]->Draw("same");
	leg2->Draw("same");
	GrowthCurves2->GetXaxis()->SetTitleOffset(1.2);
	GrowthCurves2->GetYaxis()->SetTitleOffset(1.3);
	GrowthCurves2->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves2->GetYaxis()->CenterTitle();
	GrowthCurves2->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves2->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves2->SetMinimum(ZMin);
	GrowthCurves2->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fourth subplot
	TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
	TMultiGraph*GrowthCurves3 = new TMultiGraph();
	TLegend *leg3 = new TLegend(0.55,0.70,0.9,0.85);
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);	
	GrowthCurves3->Add(FluxRate[1]);
	GrowthCurves3->Add(FluxRate[9]);
	leg3->AddEntry(FluxRate[1],"e- 80-150 MeV","p");
	leg3->AddEntry(FluxRate[9],"e+ 80-150 MeV","p");
	leg3->AddEntry(DS[4],"Daniel & Stephens 111 MeV","l");
	c3->cd();
	GrowthCurves3->Draw("AP");	
	DS[4]->Draw("same");
	leg3->Draw("same");
	GrowthCurves3->GetXaxis()->SetTitleOffset(1.2);
	GrowthCurves3->GetYaxis()->SetTitleOffset(1.3);
	GrowthCurves3->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
//	GrowthCurves3->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	GrowthCurves3->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves3->GetYaxis()->CenterTitle();
	GrowthCurves3->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves3->SetMinimum(ZMin);
	GrowthCurves3->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fifth and last subplot
	TCanvas *c4 = new TCanvas("c4","c4",1000,1000);
	TMultiGraph*GrowthCurves4 = new TMultiGraph();
	TLegend *leg4 = new TLegend(0.55,0.70,0.9,0.85);
	leg4->SetBorderSize(0);
	leg4->SetFillStyle(0);		
	GrowthCurves4->Add(FluxRate[0]);
	GrowthCurves4->Add(FluxRate[10]);
	leg4->AddEntry(FluxRate[0],"e- 150-300 MeV","p");
	leg4->AddEntry(FluxRate[10],"e+ 150-300 MeV","p");
	leg4->AddEntry(DS[5],"Daniel & Stephens 177.5 MeV","l");
	leg4->AddEntry(DS[6],"Daniel & Stephens 279 MeV","l");
	c4->cd();
	GrowthCurves4->Draw("AP");	
	DS[5]->Draw("same");
	DS[6]->Draw("same");	
	leg4->Draw("same");
	GrowthCurves4->GetXaxis()->SetTitleOffset(1.2);
	GrowthCurves4->GetYaxis()->SetTitleOffset(1.3);
	GrowthCurves4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
//	GrowthCurves4->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	GrowthCurves4->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves4->GetXaxis()->SetLimits(PMin,PMax);
	GrowthCurves4->GetYaxis()->CenterTitle();
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves4->SetMinimum(ZMin);
	GrowthCurves4->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();
	
//	c->Print(Form("%s/GrowthCurves_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c0->Print(Form("%s/GrowthCurves_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c1->Print(Form("%s/GrowthCurves_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c2->Print(Form("%s/GrowthCurves_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c3->Print(Form("%s/GrowthCurves_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c4->Print(Form("%s/GrowthCurves_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	

if(boolts) {
	//Time series plot
	gStyle->SetTextFont(85);
	gStyle->SetPadTickX(0); 
	gStyle->SetPadTickY(0);

	TCanvas *c00 = new TCanvas("c00","c00",4000,2000);
	//First subplot
	TMultiGraph*TimeSeries0 = new TMultiGraph();
	TLegend *leg00 = new TLegend(0.75,0.85,0.9,0.9);
	TimeSeries0->GetXaxis()->SetTitleOffset(1.2);
	TimeSeries0->GetYaxis()->SetTitleOffset(1.3);
	TimeSeries0->GetXaxis()->SetTitle("Time from launch (in hours)");
	TimeSeries0->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	TimeSeries0->Add(TimeSeries[4]);
	TimeSeries0->Add(TimeSeries[6]);
	leg00->AddEntry(TimeSeries[4],"e- 10-30 MeV","l");
	leg00->AddEntry(TimeSeries[6],"e+ 10-30 MeV","l");
	c00->cd();
	TimeSeries0->GetYaxis()->SetRangeUser(0,CountMax);
	TimeSeries0->GetXaxis()->SetLimits(TMin,TMax);
	TimeSeries0->Draw("ALP");		
	leg00->Draw("same");
	

	//Second subplot
	TCanvas *c11 = new TCanvas("c11","c11",4000,2000);
	TMultiGraph*TimeSeries1 = new TMultiGraph();
	TLegend *leg11 = new TLegend(0.75,0.85,0.9,0.9);
	//Lines to show time cuts
	TLine**LineCut = new TLine*[7];
	TimeSeries1->GetXaxis()->SetTitleOffset(1.2);
	TimeSeries1->GetYaxis()->SetTitleOffset(1.3);
	TimeSeries1->GetXaxis()->SetTitle("Time from launch (in hours)");
	TimeSeries1->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	TimeSeries1->Add(TimeSeries[3]);
	TimeSeries1->Add(TimeSeries[7]);
	leg11->AddEntry(TimeSeries[3],"e- 30-60 MeV","l");
	leg11->AddEntry(TimeSeries[7],"e+ 30-60 MeV","l");
	c11->cd();
	TimeSeries1->GetYaxis()->SetRangeUser(0,CountMax);
	TimeSeries1->GetXaxis()->SetLimits(TMin,TMax);
	TimeSeries1->Draw("ALP");		
	leg11->Draw("same");
	for (int i=0;i<7;i++) {
		LineCut[i] = new TLine(TCut[i],0.,1.001*TCut[i],100);
		//LineCut[i]->Draw("same");
	}
	//Third subplot
	TCanvas *c22 = new TCanvas("c22","c22",4000,2000);
	TMultiGraph*TimeSeries2 = new TMultiGraph();
	TLegend *leg22 = new TLegend(0.75,0.85,0.9,0.9);
	TimeSeries2->GetXaxis()->SetTitleOffset(1.2);
	TimeSeries2->GetYaxis()->SetTitleOffset(1.3);
	TimeSeries2->GetXaxis()->SetTitle("Time from launch (in hours)");
	TimeSeries2->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	TimeSeries2->Add(TimeSeries[2]);
	TimeSeries2->Add(TimeSeries[8]);
	leg22->AddEntry(TimeSeries[2],"e- 60-80 MeV","l");
	leg22->AddEntry(TimeSeries[8],"e+ 60-80 MeV","l");
	c22->cd();
	TimeSeries2->GetYaxis()->SetRangeUser(0,CountMax);
	TimeSeries2->GetXaxis()->SetLimits(TMin,TMax);
	TimeSeries2->Draw("ALP");		
	leg22->Draw("same");
	//Fourth subplot
	TCanvas *c33 = new TCanvas("c33","c33",4000,2000);
	TMultiGraph*TimeSeries3 = new TMultiGraph();
	TLegend *leg33 = new TLegend(0.75,0.85,0.9,0.9);
	TimeSeries3->GetXaxis()->SetTitleOffset(1.2);
	TimeSeries3->GetYaxis()->SetTitleOffset(1.3);
	TimeSeries3->GetXaxis()->SetTitle("Time from launch (in hours)");
	TimeSeries3->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	TimeSeries3->Add(TimeSeries[2]);
	TimeSeries3->Add(TimeSeries[8]);
	leg33->AddEntry(TimeSeries[1],"e- 80-150 MeV","l");
	leg33->AddEntry(TimeSeries[9],"e+ 80-150 MeV","l");
	c33->cd();
	TimeSeries3->GetYaxis()->SetRangeUser(0,CountMax);
	TimeSeries3->GetXaxis()->SetLimits(TMin,TMax);
	TimeSeries3->Draw("ALP");		
	leg33->Draw("same");
	//Fifth subplot
	TCanvas *c44 = new TCanvas("c44","c44",4000,2000);
	TMultiGraph*TimeSeries4 = new TMultiGraph();
	TLegend *leg44 = new TLegend(0.75,0.85,0.9,0.9);
	TimeSeries4->GetXaxis()->SetTitleOffset(1.2);
	TimeSeries4->GetYaxis()->SetTitleOffset(1.3);
	TimeSeries4->GetXaxis()->SetTitle("Time from launch (in hours)");
	TimeSeries4->GetYaxis()->SetTitle(Form("Count rate/%2.0fmn",60/binsperhour));
	TimeSeries4->Add(TimeSeries[0]);
	TimeSeries4->Add(TimeSeries[8]);
	leg44->AddEntry(TimeSeries[0],"e- 150-300 MeV","l");
	leg44->AddEntry(TimeSeries[10],"e+ 150-300 MeV","l");
	c44->cd();
	TimeSeries4->GetYaxis()->SetRangeUser(0,CountMax);
	TimeSeries4->GetXaxis()->SetLimits(TMin,TMax);
	TimeSeries4->Draw("ALP");		
	leg44->Draw("same");

	c00->Print(Form("%s/TimeSeriesHours_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c11->Print(Form("%s/TimeSeriesHours_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c22->Print(Form("%s/TimeSeriesHours_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c33->Print(Form("%s/TimeSeriesHours_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c44->Print(Form("%s/TimeSeriesHours_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));


	fileout->cd();
	c11->Write();
	TimeSeries1->Write();
	c22->Write();
	TimeSeries2->Write();
	c33->Write();
	TimeSeries3->Write();
	c44->Write();
	TimeSeries4->Write();
 } //if timesrie boolean	


	//Write histograms in ROOT file	
	fileout->cd();	
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

	cout << "now closing all files" << endl;
	fileout->Close();
	fileRaw->Close();
	fileCal->Close(); 	
	totxtEle.close();
	totxtPos.close();

		
} //end function	
	

void GrowthCurvesAllElectronsPressureBins(string RecoID, double cutT2, double cutT3) {

	double Tstart = 0;
	double Tend = 134;	//total time of flight in hours)
	double Launch = 22 + (7./60.);
	double Launchday = 15;
	double LaunchDay = 15 + (Launch*0.041667);
	double binsperhour = 2;
	int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
	const int NPressureBins = 19;
	double 	PressureBins[NPressureBins] = {1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.0,3.25,3.50,5.0,10,30,70,100,300,600,1000};
	const int NPBins = 6;
	double PBins[NPBins] = {10,30,60,80,150,300};
	double GeoFactor[5] = {2.68,10.43,11.55,12.02,14.64};
	//for the ascent all data is in file 000
	//for the final leg of the flight (no cut-off), data starts at file 15
	int StartFile = 0;
	int EndFile = 23;
	double TCut[7] = {19.8,30,42,56,65,82,90};
	int mcolor[7]={1,2,3,4,800+1,6,7};
    double ELoss = 4;				//4 MeV energy loss in shell + T1T2T3

	string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2/";
	string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves";
	string CutConfig = "T2Above120T3Less70";
	//string CutConfig="testfit";

	//Axis limit for plots
	float ZMax=100000;		//growth curve y-axis
	float ZMin=10;			//growth curve y-axis
	float PMin=1;			//growth curve x-axis
	float PMax=1000;			//growth curve x-axis
	float TMin=0;			//time serie x-axis
	float TMax=133;			//time serie x-axis	
	
	
	//Initialize array
	string DSConfig[7] = {"20 MeV","30 MeV", "44.8 MeV", "70 MeV", "111 MeV", "177.5 MeV","279 MeV"};
	float**DScoeff=new float*[7];
	for(int i=0; i<7;i++) 
	 {
	 DScoeff[i] =  new float[7];
	 for(int j=0; j<7;j++)			//7 coefficients for each array 
	  {
		 DScoeff[i][j]=0.;
	 }
    }
	 string DSparamfile="./DScoefficients.dat"; 
	 LoadDSparameters(DSparamfile,DScoeff);	

	double T2Cut = cutT2;	//value of cut on T2
	double T3Cut = cutT3;	//value of cut on T3
	TFile *fileRaw;
	TFile *fileCal;
	string ID = RecoID;	

	//Create histograms
	TH1F***B2Day = new TH1F**[NPBins];
	TH1F***B2Night = new TH1F**[NPBins];
	TH1F***P0BinDay = new TH1F**[NPBins];
	TH1F***P0BinNight = new TH1F**[NPBins];	


	//Momentum histograms and pressure
	for(int i=0; i<NPBins; i++) 
	{
			P0BinDay[i] = new TH1F*[NPressureBins];
			P0BinNight[i] = new TH1F*[NPressureBins];
			B2Day[i] = new TH1F*[NPressureBins];
			B2Night[i] = new TH1F*[NPressureBins];
			for(int j=0; j<NPressureBins; j++) 
			 {
				P0BinDay[i][j] = new TH1F(Form("unsigned p0 h%d%d",i,j), Form("unsigned p0 h%d%d",i,j), 1,0,1);
				P0BinNight[i][j] = new TH1F(Form("unsigned p0 night h%d%d",i,j), Form("unsigned p0 night h%d%d",i,j), 1,0,1);
                B2Day[i][j] = new TH1F(Form("B2Day h%d%d",i,j),Form("B2Day h%d%d",i,j), 10000, 0, 1000);
                B2Night[i][j] = new TH1F(Form("B2Night h%d%d",i,j),Form("B2Night h%d%d",i,j), 10000, 0, 1000);
			}
		}

	
   //Determine the duration payload spends in each pressure bin
	
	double durationDay[NPressureBins];
	double durationNight[NPressureBins];
	for(int k=0;k<NPressureBins;k++) 
	{
		durationDay[k]=0;
		durationNight[k]=0;
	}
	float startperiod;			//initialize timefirst event
	int prevbinlayer=17;			//pressure bin of first even
	float prevtimelayer=0.0058328;		//first positive time measurement made
    float tmptimePHA;	
	//ROOT fileout
	 //TFile *fileout = new TFile(Form("%s/GrowthCurvesAllElectrons_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
     for (int j=StartFile; j<EndFile; j++)
	 {
	 	if (j<10) fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 	if (j>=10) fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 
	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	
	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
     {
		treeCal->GetEntry(i);
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
	    tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
		double yCT1 = (double)de->get_yCT1();
		double mCT1 = (double)de->get_mCT1();
		double dCT1 = (double)de->get_dCT1();
		double hCT1= (double)de->get_hCT1();
		double miCT1 = (double)de->get_miCT1();
		double sCT1 = (double)de->get_sCT1();
		float tmptimeCT1 = hCT1 + (miCT1/60) + (sCT1/3600);
		tmptimeCT1 = tmptimeCT1-Launch;
		tmptimeCT1 = tmptimeCT1 + 24*(dCT1 - Launchday);
		if (tmptimePHA<0) continue;
	//Barometer data (calibrated)
		float B2Pres = de->get_PressB2();
	//Convert pressure to g/cm^-2
		B2Pres = B2Pres* 1.3595;
		//cout << "t = " << tmptimePHA <<  ", B2Pres = " << B2Pres << " g.cm^-2 " << endl;
	//Determine the correct time, energy bin and pressure bin
		int pressurebin;
		for(int k=0;k<NPressureBins;k++) {
			if((PressureBins[k] <= B2Pres) & (B2Pres < PressureBins[k+1]))
			{
				pressurebin = k;
				//cout <<  "t = " << tmptimePHA << "h, B2 = " << B2Pres << " bin = " << pressurebin << " ["<<PressureBins[k] <<","<<PressureBins[k+1]<<"]"<<endl;
				break;
			}	
		}
		
		  //  if(pressurebin == prevbinlayer) cout << "same pressure bin " << pressurebin << endl;
			if(pressurebin != prevbinlayer) {
			//	cout << "different pressure bin !" << endl;
					        if (tmptimePHA>0 & tmptimePHA<5) {
								durationDay[prevbinlayer]+= tmptimePHA - startperiod;
								//cout << "pressbin = " << prevbinlayer << " durationDay = " << durationDay[prevbinlayer] << " h " << endl;
						   }
				            else if (tmptimePHA>90 ) {
								durationNight[prevbinlayer]+= tmptimePHA - startperiod;
							//   cout << "pressbin = " << prevbinlayer << " durationNight = " << durationNight[prevbinlayer] << " h " << endl;	
							}
				startperiod = tmptimePHA;
				prevbinlayer = pressurebin;
			}

	} // end i
		 //Fill in last bin
         if (tmptimePHA>0 & tmptimePHA<5) durationDay[prevbinlayer]+= tmptimePHA - startperiod;
		 else if (tmptimePHA>90 )  durationNight[prevbinlayer]+= tmptimePHA - startperiod;
		 startperiod = tmptimePHA;
		// cout << "end of file " << j << ", t = " << tmptimePHA << endl; 
		 double sumDay=0;
		 double sumNight=0;
		 for(int k=0;k<NPressureBins;k++) {
			cout << " k = " << k << " durationDay[i] = " << durationDay[k] << endl;
			 cout << " k = " << k << " durationNight[i] =  " << durationNight[k] << endl;
		     sumDay+=durationDay[k];
			 sumNight+=durationNight[k];
		 }
		 cout << "sumDay = " << sumDay << endl;
		 cout << "sumNight= " << sumNight <<endl;
 } //end j

	
	//ROOT fileout
	TFile *fileout = new TFile(Form("%s/GrowthCurvesAllElectronsPressureBins_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
     for (int j=StartFile; j<EndFile; j++)
	 {
	 	if (j<10) fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 	if (j>=10) fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 
	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	
	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
     {
		
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
	//Time of event from PHA
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
		float tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
	 // cout << "tmptimeCT1 = " << tmptime << ", tmptimePHA = " << tmptimePHA << endl; 
		
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
		//cout << "B2Pres = " << B2Pres << "g.cm^-2 " << endl;
	//Reconstructed energy 
		float pPR=1000*fabs(de->get_p0PR());   //in MeV
		double deflecPR = de->get_deflecPR();
		float preco=de->get_p0reco();
	    float P0sign;
		float ePR = 1000*(de->get_EkPR());		//total energy in MeV
		if (strcmp(ID.c_str(),"PRonly")==0)  P0sign = pPR*TMath::Sign(1,deflecPR);
		else P0sign = preco*TMath::Sign(1,deflecPR);	
		double chiNBPR = de->get_chi2NBPR();
		double chiBPR = de->get_chi2BPR();	
	    pPR = pPR + ELoss;					//correct for energy loss in shell + T1 + T2 + T3
	    ePR = ePR + ELoss;					//correct for energy loss in shell + T1 + T2 + T3
	//Determine the correct time, energy bin and pressure bin
		int tbin = tmptimePHA*binsperhour;
		int pbin;
		int pressurebin;
		for(int k=0;k<NPBins;k++) {if((PBins[k] <= pPR) & (pPR < PBins[k+1])) pbin = k; }
		for(int k=0;k<NPressureBins;k++) {if((PressureBins[k] <= B2Pres) & (B2Pres < PressureBins[k+1])) pressurebin = k;}			
		
		
	//Get DAYTIME values
		if(tmptimePHA>0 & tmptimePHA<5) 
			{
	//Apply cuts to only get clean events
		if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
		   {
	      //  cout << "t = " << tmptimePHA << ", ePR = " << ePR << " MeV, at B2 = " << B2Pres << ", tbin = " << tbin <<", pressurebin = " << pressurebin << endl;
			P0BinDay[pbin][pressurebin]->Fill(ePR);
		    B2Day[pbin][pressurebin]->Fill(B2Pres);
			}	   
		}

		
	//Get NIGHTIME values
	if(tmptimePHA>90) 
		{
//Apply cuts to only get clean events
	if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
	   {
	    //cout << "t = " << tmptimePHA << ", pPR = " << pPR << "MeV, at ePR  =  " << ePR << " B2 = " << B2Pres << ", tbin = " << tbin <<", pressurebin = " << pressurebin << endl;
			P0BinNight[pbin][pressurebin]->Fill(ePR);
		    B2Night[pbin][pressurebin]->Fill(B2Pres);

		}	   
	  }
	}     //end i, loop on events
}		//end j, loop on files

	
	//Make growth curves
	int ngraphs = (NPBins/2)-1;
	TGraphAsymmErrors**CountRate = new TGraphAsymmErrors*[NPBins];
    TGraphAsymmErrors**FluxRate = new TGraphAsymmErrors*[NPBins];
	TGraphAsymmErrors**TimeSeries = new TGraphAsymmErrors*[NPBins];
    TGraphAsymmErrors**FluxRateNight = new TGraphAsymmErrors*[NPBins];


	for (int j=0; j<NPBins;j++) 
	{
		CountRate[j] = new TGraphAsymmErrors();
	    CountRate[j]->SetMarkerStyle(kOpenCircle);
		CountRate[j]->SetMarkerSize(0.8);
		FluxRate[j] = new TGraphAsymmErrors();
	    FluxRate[j]->SetMarkerStyle(kOpenCircle);
		FluxRate[j]->SetMarkerSize(0.8);
		TimeSeries[j] = new TGraphAsymmErrors();
	    TimeSeries[j]->SetMarkerStyle(kDot);
		TimeSeries[j]->SetMarkerSize(0.8);
		FluxRateNight[j] = new TGraphAsymmErrors();
	    FluxRateNight[j]->SetMarkerStyle(kOpenCircle);
		FluxRateNight[j]->SetMarkerSize(0.8); 
		CountRate[j]->SetMarkerColor(kBlue);
		FluxRate[j]->SetMarkerColor(kBlue);
		FluxRate[j]->SetLineColor(kBlue);
		TimeSeries[j]->SetMarkerColor(kBlue);
		TimeSeries[j]->SetLineColor(kBlue);
		FluxRateNight[j]->SetMarkerColor(kRed);
		FluxRateNight[j]->SetLineColor(kRed);

		int npNight = 1;
		int npDay = 1;

		for (int k=0;k<NPressureBins;k++) 
		{
			double RateNight,RateDay;
			double FluxNight,FluxDay;
			double RateDayErr,RateNightErr;
			double FluxNightErr,FluxDayErr;			
			double time = Tstart + (k+1)/(binsperhour);
			double timeWidthDay = durationDay[k]*60*60;		       	//width of time bin in seconds
			double timeWidthNight = durationNight[k]*60*60;		       	//width of time bin in seconds
			double MomWidth;
		    double B2meanDay = B2Day[j][k]->GetMean();
		    double B2meanNight = B2Night[j][k]->GetMean();
			cout << "k = " << k << " B2Day mean pressure " << B2meanDay << "gcm-2 " << " , duration = " << timeWidthDay << endl;
			cout << "k = " << k << " B2Night mean pressure " << B2meanNight << "gcm-2 " << " , duration = " << timeWidthNight << endl;
		    double timeInDays = LaunchDay + time*(0.041667);
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			//RateDay
			RateDay = P0BinDay[j][k]->GetEntries(); 
			RateDayErr = TMath::Sqrt(P0BinDay[j][k]->GetEntries());
			FluxDay = RateDay/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthDay);
			FluxDayErr = RateDayErr/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthDay);
			if(B2meanDay!=0 && FluxDay!=0) {
				cout << "j = " << j <<", RateDay = " << RateDay << endl;
				CountRate[j]->SetPoint(npDay,B2meanDay,RateDay);
				FluxRate[j]->SetPoint(npDay,B2meanDay,FluxDay);
				FluxRate[j]->SetPointError(npDay,0,0,FluxDayErr,FluxDayErr);
				npDay++;
				}
			//Rate Night
			RateNight = P0BinNight[j][k]->GetEntries(); 
			RateNightErr = TMath::Sqrt(P0BinNight[j][k]->GetEntries());
			FluxNight = RateNight/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthNight);
			FluxNightErr = RateNightErr/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthNight);
			if(B2meanNight!=0 && FluxNight!=0) {
				cout << "j = " << j <<", RateNight = " << RateNight << endl;
				FluxRateNight[j]->SetPoint(npNight,B2meanNight,FluxNight);
				FluxRateNight[j]->SetPointError(npNight,0,0,FluxNightErr,FluxNightErr);
				npNight++;
				}
				
			    
		}
    }
	
	//Define Daniel & Stephens functions
	TF1**DS = new TF1*[7];
	for (int k=0;k<7;k++) 
	{
		DS[k] = new TF1(Form("ds%d",k),DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DS[k]->FixParameter(l,DScoeff[k][l]);
		DS[k]->SetParameter(7,1);
		DS[k]->SetLineColor(kBlack);
		DS[k]->SetLineStyle(6);

	}
	
	//////////////////////////////////////////////////
	/////////////////DRAW CANVASES////////////////////
	//////////////////////////////////////////////////
	
	gStyle->SetTextFont(85);
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1);	
	
	//First growth curve subplot
	TCanvas *c0 = new TCanvas("c0","c0",1000,1000);
	TMultiGraph*GrowthCurves0 = new TMultiGraph();
	TLegend *leg0 = new TLegend(0.55,0.70,0.9,0.85);
	TPaveText *pt0 = new TPaveText(0.75,0.65,0.9,0.75);
	leg0->SetBorderSize(0);
	leg0->SetFillStyle(0);
	pt0->SetBorderSize(0);
	pt0->SetFillStyle(0);
	CountRate[4]->SetMarkerColor(kBlue);
	GrowthCurves0->Add(FluxRate[0]);
	GrowthCurves0->Add(FluxRateNight[0]);
	leg0->AddEntry(FluxRate[0],"e^{-}+e^{+} 10-30 MeV, ascent","p");
	leg0->AddEntry(FluxRateNight[0],"e^{-}+e^{+} 10-30 MeV, float","p");
	leg0->AddEntry(DS[0],"Daniel & Stephens 20 MeV","l");
	c0->cd();
	GrowthCurves0->Draw("AP");
	FluxRate[0]->Fit("ds0","","",50,1000);
	pt0->AddText(Form("Normalization factor %2.2f",DS[0]->GetParameter(7)));
	DS[0]->DrawF1(1,1000,"same");
	leg0->Draw("same");
    pt0->Draw();
	GrowthCurves0->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves0->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves0->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves0->GetYaxis()->SetTitleSize(0.05);	
	GrowthCurves0->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves0->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves0->GetYaxis()->CenterTitle();
	GrowthCurves0->GetXaxis()->CenterTitle();	
	GrowthCurves0->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves0->SetMinimum(ZMin);
	GrowthCurves0->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Second subplot
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TMultiGraph*GrowthCurves1 = new TMultiGraph();
	TLegend *leg1 = new TLegend(0.55,0.70,0.9,0.85);
	TPaveText *pt1 = new TPaveText(0.75,0.65,0.9,0.75,"NDC");
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	pt1->SetBorderSize(0);
	pt1->SetFillStyle(0);
	GrowthCurves1->Add(FluxRate[1]);
	leg1->AddEntry(FluxRate[1],"e^{-}+e^{+} 30-60 MeV, ascent","p");
	GrowthCurves1->Add(FluxRateNight[1]);
	leg1->AddEntry(FluxRateNight[1],"e^{-}+e^{+} 30-60 MeV, float","p");
	leg1->AddEntry(DS[2],"Daniel & Stephens 45 MeV","l");
	c1->cd();
	GrowthCurves1->Draw("AP");		
	FluxRate[1]->Fit("ds2","","",50,1000);
	DS[2]->DrawF1(1,1000,"same");
	pt1->AddText(Form("Normalization factor %2.2f",DS[2]->GetParameter(7)));
	leg1->Draw("same");
    pt1->Draw();	
	GrowthCurves1->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves1->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves1->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves1->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves1->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves1->GetYaxis()->CenterTitle();
	GrowthCurves1->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves1->SetMinimum(ZMin);
	GrowthCurves1->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Third subplot
	TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
	TMultiGraph*GrowthCurves2 = new TMultiGraph();
	TLegend *leg2 = new TLegend(0.55,0.70,0.9,0.85);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	TPaveText *pt2 = new TPaveText(0.75,0.65,0.9,0.75);
	pt2->SetBorderSize(0);
	pt2->SetFillStyle(0);
	GrowthCurves2->Add(FluxRate[2]);
	leg2->AddEntry(FluxRate[2],"e^{-}+e^{+} 60-80 MeV, ascent","p");
	GrowthCurves2->Add(FluxRateNight[2]);
	leg2->AddEntry(FluxRateNight[2],"e^{-}+e^{+} 60-80 MeV, float","p");
	leg2->AddEntry(DS[3],"Daniel & Stephens 70 MeV","l");
	c2->cd();
	GrowthCurves2->Draw("AP");	
	FluxRate[2]->Fit("ds3","","",50,1000);
	DS[3]->DrawF1(1,1000,"same");
	pt2->AddText(Form("Normalization factor %2.2f",DS[3]->GetParameter(7)));
	leg2->Draw("same");
    pt2->Draw();
	GrowthCurves2->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves2->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves2->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves2->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves2->GetXaxis()->CenterTitle();
	GrowthCurves2->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves2->GetYaxis()->CenterTitle();
	GrowthCurves2->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves2->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves2->SetMinimum(ZMin);
	GrowthCurves2->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fourth subplot
	TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
	TMultiGraph*GrowthCurves3 = new TMultiGraph();
	TLegend *leg3 = new TLegend(0.55,0.70,0.9,0.85);
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);	
	TPaveText *pt3 = new TPaveText(0.75,0.65,0.9,0.75,"NDC");
	pt3->SetBorderSize(0);
	pt3->SetFillStyle(0);
	GrowthCurves3->Add(FluxRate[3]);
	leg3->AddEntry(FluxRate[3],"e^{-}+e^{+} 80-150 MeV, ascent","p");
	GrowthCurves3->Add(FluxRateNight[3]);
	leg3->AddEntry(FluxRateNight[3],"e^{-}+e^{+} 80-150 MeV, float","p");
	leg3->AddEntry(DS[4],"Daniel & Stephens 111 MeV","l");
	c3->cd();
	GrowthCurves3->Draw("AP");	
	FluxRate[3]->Fit("ds4","","",50,1000);
	DS[4]->DrawF1(1,1000,"same");
	pt3->AddText(Form("Normalization factor %2.2f",DS[4]->GetParameter(7)));
	leg3->Draw("same");
    pt3->Draw("same");
	GrowthCurves3->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves3->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves3->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves3->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves3->GetXaxis()->CenterTitle();	
	GrowthCurves3->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves3->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves3->GetYaxis()->CenterTitle();
	GrowthCurves3->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves3->SetMinimum(ZMin);
	GrowthCurves3->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fifth and last subplot
	TCanvas *c4 = new TCanvas("c4","c4",1000,1000);
	TMultiGraph*GrowthCurves4 = new TMultiGraph();
	TLegend *leg4 = new TLegend(0.55,0.70,0.9,0.85);
	leg4->SetBorderSize(0);
	leg4->SetFillStyle(0);
	TPaveText *pt4 = new TPaveText(2,20000,30,30000);
	pt4->SetBorderSize(0);
	pt4->SetFillStyle(0);
	GrowthCurves4->Add(FluxRate[4]);
	leg4->AddEntry(FluxRate[4],"e^{-}+e^{+} 150-300 MeV, ascent","p");
	GrowthCurves4->Add(FluxRateNight[4]);
	leg4->AddEntry(FluxRateNight[4],"e^{-}+e^{+} 150-300 MeV, float","p");
	leg4->AddEntry(DS[5],"Daniel & Stephens 177.5 MeV","l");
	c4->cd();
	GrowthCurves4->Draw("AP");	
	FluxRate[4]->Fit("ds5","","",30,500);
	DS[5]->DrawF1(1,1000,"same");
	pt4->AddText(Form("Normalization factor %2.2f",DS[5]->GetParameter(7)));
	cout << "p7 = " << DS[5]->GetParameter(7) << endl;
	leg4->Draw("same");
    pt4->Draw();
	GrowthCurves4->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves4->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves4->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves4->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves4->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves4->GetXaxis()->SetLimits(PMin,PMax);
	GrowthCurves4->GetYaxis()->CenterTitle();
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves4->SetMinimum(ZMin);
	GrowthCurves4->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();
	
	//////////////////////////////////////////////////
	/////////////////DRAW CANVASES////////////////////
	//////////////////////////////////////////////////
	
	gStyle->SetTextFont(85);
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1);	

	//Draw All Daniel & Stephens functions
	TCanvas *c = new TCanvas("c","c",1000,1000);
	TMultiGraph*Functions = new TMultiGraph();
	CountRate[4]->SetMarkerColor(kWhite);
	Functions->Add(CountRate[4]);
	TLegend *leg = new TLegend(0.20,0.70,0.3,0.85);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	c->cd();
	Functions->Draw("AP");
	for (int k=0;k<7;k++) 
	{
		DS[k]->Draw("same");
	    leg->AddEntry(DS[k],Form("D&S %s", DSConfig[k].c_str()),"l");
	}
	leg->Draw("same");
	Functions->GetXaxis()->SetTitleOffset(1.2);
	Functions->GetYaxis()->SetTitleOffset(1.3);
	Functions->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	Functions->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	Functions->GetYaxis()->CenterTitle();	
	Functions->GetXaxis()->SetLimits(PMin,2000);
	gPad->SetLogy();
	gPad->SetLogx();
	Functions->SetMinimum(ZMin);
	Functions->SetMaximum(100000);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 
	
	c0->Print(Form("%s/GrowthCurvesAllElectrons_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c1->Print(Form("%s/GrowthCurvesAllElectrons_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c2->Print(Form("%s/GrowthCurvesAllElectrons_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c3->Print(Form("%s/GrowthCurvesAllElectrons_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c4->Print(Form("%s/GrowthCurvesAllElectrons_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
    c->Print(Form("%s/GrowthCurvesAllElectrons_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));




	//Write histograms in ROOT file	
	fileout->cd();	
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

	cout << "now closing all files" << endl;
	//fileout->Close();
	//fileCal->Close(); 	
			
} //end function		
	
void GrowthCurvesSignedPressureBins(string RecoID, double cutT2, double cutT3) {

	double Tstart = 0;
	double Tend = 134;	//total time of flight in hours)
	double Launch = 22 + (7./60.);
	double Launchday = 15;
	double LaunchDay = 15 + (Launch*0.041667);
	double binsperhour = 2;
	int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
	const int NPressureBins = 18;
	double 	PressureBins[NPressureBins] = {1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.0,3.25,3.50,5.0,10,30,70,100,300,1000};
	const int NPBins = 12;
	double PBins[NPBins] = {-300,-150,-80,-60,-30,-10,10,30,60,80,150,300};
	const int NPBinsAll = 6;
	double PBinsAll[NPBinsAll] = {10,30,60,80,150,300};	
	double GeoFactor[5] = {2.68,10.43,11.55,12.02,14.64};
	//for the ascent all data is in file 000
	//for the final leg of the flight (no cut-off), data starts at file 15
	int StartFile = 0;
	int EndFile = 23;
	double TCut[7] = {19.8,30,42,56,65,82,90};
	int mcolor[7]={1,2,3,4,800+1,6,7};
    double ELoss = 4;				//4 MeV energy loss in shell + T1T2T3

	string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2/";
	string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves";
	string OutpathText = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/PythonScripts";
	string CutConfig = "T2Above120T3Less110";
	//string CutConfig="testfit";

	//txt output file for positrons
	ofstream totxtPos;
	totxtPos.open(Form("%s/GrowthCurvePos_60to80MeV.txt",OutpathText.c_str()));
	totxtPos << " Pressure in g.cm-2 " << " Flux " << " Std Dev " << endl;
	//txt output file for electrons
	ofstream totxtEle;
	totxtEle.open(Form("%s/GrowthCurveEle_60to80MeV.txt",OutpathText.c_str()));
	totxtEle << " Pressure in g.cm-2 " << " Flux " << " Std Dev " << endl;
	ofstream totxtAll;
	totxtAll.open(Form("%s/GrowthCurveAll_60to80MeV.txt",OutpathText.c_str()));
	totxtAll << " Pressure in g.cm-2 " << " Flux " << " Std Dev " << endl;
	
	//Axis limit for plots
	float ZMax=100000;		//growth curve y-axis
	float ZMin=1;			//growth curve y-axis
	float PMin=1;			//growth curve x-axis
	float PMax=1000;			//growth curve x-axis
	float TMin=0;			//time serie x-axis
	float TMax=133;			//time serie x-axis	
	
	
	//Initialize array
	string DSConfig[7] = {"20 MeV","30 MeV", "44.8 MeV", "70 MeV", "111 MeV", "177.5 MeV","279 MeV"};
	float**DScoeff=new float*[7];
	for(int i=0; i<7;i++) 
	 {
	 DScoeff[i] =  new float[7];
	 for(int j=0; j<7;j++)			//7 coefficients for each array 
	  {
		 DScoeff[i][j]=0.;
	 }
    }
	 string DSparamfile="./DScoefficients.dat"; 
	 LoadDSparameters(DSparamfile,DScoeff);	

	double T2Cut = cutT2;	//value of cut on T2
	double T3Cut = cutT3;	//value of cut on T3
	TFile *fileRaw;
	TFile *fileCal;
	string ID = RecoID;	

	//Create histograms
	TH1F***B2Day = new TH1F**[NPBins];
	TH1F***B2Night = new TH1F**[NPBins];
	TH1F***P0BinDay = new TH1F**[NPBins];
	TH1F***P0BinNight = new TH1F**[NPBins];	
	TH1F***P0BinAll = new TH1F**[NPBinsAll];
	TH1F***B2All = new TH1F**[NPBinsAll];

	
	//Momentum histograms and pressure
	for(int i=0; i<NPBinsAll; i++) 
	{
			P0BinAll[i] = new TH1F*[NPressureBins];
			B2All[i] = new TH1F*[NPressureBins];
			for(int j=0; j<NPressureBins; j++) 
			 {
				P0BinAll[i][j] = new TH1F(Form("all p0 h%d%d",i,j), Form("all p0 h%d%d",i,j), 1,0,1);
                B2All[i][j] = new TH1F(Form("B2All h%d%d",i,j),Form("B2All h%d%d",i,j), 10000, 0, 1000);
			}
		}


	//Momentum histograms and pressure
	for(int i=0; i<NPBins; i++) 
	{
			P0BinDay[i] = new TH1F*[NPressureBins];
			P0BinNight[i] = new TH1F*[NPressureBins];
			B2Day[i] = new TH1F*[NPressureBins];
			B2Night[i] = new TH1F*[NPressureBins];
			for(int j=0; j<NPressureBins; j++) 
			 {
				P0BinDay[i][j] = new TH1F(Form("unsigned p0 h%d%d",i,j), Form("unsigned p0 h%d%d",i,j), 1,0,1);
				P0BinNight[i][j] = new TH1F(Form("unsigned p0 night h%d%d",i,j), Form("unsigned p0 night h%d%d",i,j), 1,0,1);
                B2Day[i][j] = new TH1F(Form("B2Day h%d%d",i,j),Form("B2Day h%d%d",i,j), 10000, 0, 1000);
                B2Night[i][j] = new TH1F(Form("B2Night h%d%d",i,j),Form("B2Night h%d%d",i,j), 10000, 0, 1000);
			}
		}

	
   //Determine the duration payload spends in each pressure bin
	
	double durationDay[NPressureBins];
	double durationNight[NPressureBins];
	for(int k=0;k<NPressureBins;k++) 
	{
		durationDay[k]=0;
		durationNight[k]=0;
	}
	float startperiod;			//initialize timefirst event
	int prevbinlayer=17;			//pressure bin of first even
	float prevtimelayer=0.0058328;		//first positive time measurement made
    float tmptimePHA;	
	//ROOT fileout
	 //TFile *fileout = new TFile(Form("%s/GrowthCurvesAllElectrons_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
     for (int j=StartFile; j<EndFile; j++)
	 {
	 	if (j<10) fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 	if (j>=10) fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 
	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	
	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
     {
		treeCal->GetEntry(i);
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
	    tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
		if (tmptimePHA<0) continue;
	//Barometer data (calibrated)
		float B2Pres = de->get_PressB2();
	//Convert pressure to g/cm^-2
		B2Pres = B2Pres* 1.3595;
		//cout << "t = " << tmptimePHA <<  ", B2Pres = " << B2Pres << " g.cm^-2 " << endl;
	//Determine the correct time, energy bin and pressure bin
		int pressurebin;
		for(int k=0;k<NPressureBins;k++) {
			if((PressureBins[k] <= B2Pres) & (B2Pres < PressureBins[k+1]))
			{
				pressurebin = k;
				//cout <<  "t = " << tmptimePHA << "h, B2 = " << B2Pres << " bin = " << pressurebin << " ["<<PressureBins[k] <<","<<PressureBins[k+1]<<"]"<<endl;
				break;
			}	
		}
		
		  //  if(pressurebin == prevbinlayer) cout << "same pressure bin " << pressurebin << endl;
			if(pressurebin != prevbinlayer) {
			//	cout << "different pressure bin !" << endl;
					        if (tmptimePHA>0 & tmptimePHA<2.5) {
								durationDay[prevbinlayer]+= tmptimePHA - startperiod;
								//cout << "pressbin = " << prevbinlayer << " durationDay = " << durationDay[prevbinlayer] << " h " << endl;
						   }
				            else if (tmptimePHA>90)  {
								durationNight[prevbinlayer]+= tmptimePHA - startperiod;
							//   cout << "pressbin = " << prevbinlayer << " durationNight = " << durationNight[prevbinlayer] << " h " << endl;	
							}
				startperiod = tmptimePHA;
				prevbinlayer = pressurebin;
			}

	} // end i
		 //Fill in last bin
         if (tmptimePHA>0 & tmptimePHA<5) durationDay[prevbinlayer]+= tmptimePHA - startperiod;
		 else if (tmptimePHA>90 )  durationNight[prevbinlayer]+= tmptimePHA - startperiod;
		 startperiod = tmptimePHA;
		// cout << "end of file " << j << ", t = " << tmptimePHA << endl; 
		 double sumDay=0;
		 double sumNight=0;
		 for(int k=0;k<NPressureBins;k++) {
			//cout << " k = " << k << " durationDay[i] = " << durationDay[k] << endl;
			 //cout << " k = " << k << " durationNight[i] =  " << durationNight[k] << endl;
		     sumDay+=durationDay[k];
			 sumNight+=durationNight[k];
		 }
 } //end j
	double durationAll[NPressureBins];
	for(int k=0;k<NPressureBins;k++) 
	{
		durationAll[k]=durationDay[k]+durationNight[k];
	}
	
	//ROOT fileout
	TFile *fileout = new TFile(Form("%s/GrowthCurvesSignedPressureBins_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
     for (int j=StartFile; j<EndFile; j++)
	 {
	 	if (j<10) fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 	if (j>=10) fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 
	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	
	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
     {
		
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
	//Time of event from PHA
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
		float tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
	 // cout << "tmptimeCT1 = " << tmptime << ", tmptimePHA = " << tmptimePHA << endl; 
		
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
		//cout << "B2Pres = " << B2Pres << "g.cm^-2 " << endl;
	//Reconstructed energy 
		float pPR=1000*fabs(de->get_p0PR());   //in MeV
		double deflecPR = de->get_deflecPR();
		float preco=de->get_p0reco();
	    float P0sign;
		float ePR = 1000*(de->get_EkPR());		//total energy in MeV
		if (strcmp(ID.c_str(),"PRonly")==0)  P0sign = pPR*TMath::Sign(1,deflecPR);
		else P0sign = preco*TMath::Sign(1,deflecPR);	
		double chiNBPR = de->get_chi2NBPR();
		double chiBPR = de->get_chi2BPR();	
	    pPR = pPR + ELoss;					//correct for energy loss in shell + T1 + T2 + T3
	    ePR = ePR + ELoss;					//correct for energy loss in shell + T1 + T2 + T3
		float ePRsign = ePR*TMath::Sign(1,deflecPR);
	//Determine the correct time, energy bin and pressure bin
		int tbin = tmptimePHA*binsperhour;
		int pbin;
		int pressurebin;
		int pbinall;
		for(int k=0;k<NPBins;k++) {if((PBins[k] <= ePRsign) & (ePRsign < PBins[k+1])) pbin = k;}
		for(int k=0;k<NPressureBins;k++) {if((PressureBins[k] <= B2Pres) & (B2Pres < PressureBins[k+1])) pressurebin = k;}			
		for(int k=0;k<NPBinsAll;k++) {if((PBinsAll[k] <= ePR) & (ePR < PBinsAll[k+1])) pbinall = k;}
		
		
	//Get DAYTIME values
		if(tmptimePHA>0 & tmptimePHA<2.5) 
			{
	//Apply cuts to only get clean events
		if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
		   {
	      //  cout << "t = " << tmptimePHA << ", ePR = " << ePR << " MeV, at B2 = " << B2Pres << ", tbin = " << tbin <<", pressurebin = " << pressurebin << endl;
			P0BinDay[pbin][pressurebin]->Fill(ePRsign);
		    B2Day[pbin][pressurebin]->Fill(B2Pres);
			P0BinAll[pbinall][pressurebin]->Fill(ePR);
			B2All[pbinall][pressurebin]->Fill(B2Pres);

			}	   
		}
	
	//Get NIGHTIME values
	if(tmptimePHA>90) 
		{
//Apply cuts to only get clean events
	if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
	   {
	    //cout << "t = " << tmptimePHA << ", pPR = " << pPR << "MeV, at ePR  =  " << ePR << " B2 = " << B2Pres << ", tbin = " << tbin <<", pressurebin = " << pressurebin << endl;
			P0BinNight[pbin][pressurebin]->Fill(ePRsign);
		    B2Night[pbin][pressurebin]->Fill(B2Pres);
			P0BinAll[pbinall][pressurebin]->Fill(ePR);
			B2All[pbinall][pressurebin]->Fill(B2Pres);
		}	   
	  }
	}     //end i, loop on events
}		//end j, loop on files

	
	//Make growth curves
	int ngraphs = (NPBins/2)-1;
    TGraphAsymmErrors**FluxRateDay = new TGraphAsymmErrors*[NPBins];
    TGraphAsymmErrors**FluxRateNight = new TGraphAsymmErrors*[NPBins];
 
	for (int j=0; j<NPBins-1;j++) 
	{
		FluxRateDay[j] = new TGraphAsymmErrors();
	    FluxRateDay[j]->SetMarkerStyle(kOpenCircle);
		FluxRateDay[j]->SetMarkerSize(0.8);	
		FluxRateNight[j] = new TGraphAsymmErrors();
	    FluxRateNight[j]->SetMarkerStyle(kOpenCircle);
		FluxRateNight[j]->SetMarkerSize(0.8); 

//Electrons
		if (j<5)  { 
		FluxRateDay[j]->SetMarkerColor(kBlue);
		FluxRateDay[j]->SetLineColor(kBlue);
		FluxRateNight[j]->SetMarkerColor(kBlue);
		FluxRateNight[j]->SetLineColor(kBlue);
		}
//Positrons	
		if (j>5)  { 
		FluxRateDay[j]->SetMarkerColor(kRed);
		FluxRateDay[j]->SetLineColor(kRed);
		FluxRateNight[j]->SetMarkerColor(kRed);
		FluxRateNight[j]->SetLineColor(kRed);
		}
		int npNightPos = 1;
		int npDayPos = 1;
		int npNightEle = 1;
		int npDayEle = 1;
		for (int k=0;k<NPressureBins;k++) 
		{
			double RateNightEle,RateDayEle;
			double RateDayPos,RateNightPos;
			double FluxNightEle,FluxDayEle;
			double FluxNightPos,FluxDayPos;
			double RateDayEleErr,RateNightEleErr;
			double RateDayPosErr,RateNightPosErr;
			double FluxDayEleErr,FluxNightEleErr;
			double FluxDayPosErr,FluxNightPosErr;			
			double time = Tstart + (k+1)/(binsperhour);
			double timeWidthDay = durationDay[k]*60*60;		       	//width of time bin in seconds
			double timeWidthNight = durationNight[k]*60*60;		       	//width of time bin in seconds
			double MomWidth;
		    double B2meanDay = B2Day[j][k]->GetMean();
		    double B2meanNight = B2Night[j][k]->GetMean();
			//cout << "j = " << k << ", k = " << k << " B2Day mean pressure " << B2meanDay << "gcm-2 " << " , duration = " << timeWidthDay << endl;
		//	cout << "j = " << j << ", k = " << k << " B2Night mean pressure " << B2meanNight << "gcm-2 " << " , duration = " << timeWidthNight << endl;
			
			
			if (j<5) 
				{			
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			//RateDay
			RateDayEle = P0BinDay[j][k]->GetEntries(); 
			RateDayEleErr = TMath::Sqrt(P0BinDay[j][k]->GetEntries());
			FluxDayEle = RateDayEle/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthDay);
			FluxDayEleErr = RateDayEleErr/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthDay);
          cout << "j = " << k << ", k = " << k << " RateDayEle " << RateDayEle << ", FluxDayEle = " << FluxDayEle << endl;

			if(B2meanDay!=0 && FluxDayEle!=0) {
				FluxRateDay[j]->SetPoint(npDayEle,B2meanDay,FluxDayEle);
				FluxRateDay[j]->SetPointError(npDayEle,0,0,FluxDayEleErr,FluxDayEleErr);
				}
			//Rate Night
			RateNightEle = P0BinNight[j][k]->GetEntries(); 
			RateNightEleErr = TMath::Sqrt(P0BinNight[j][k]->GetEntries());
			FluxNightEle = RateNightEle/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthNight);
			FluxNightEleErr = RateNightEleErr/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidthNight);
          cout << "j = " << k << ", k = " << k << " B2meanNight " << B2meanNight << ", FluxNightEle = " << FluxNightEle << endl;

			if(B2meanNight!=0 && FluxNightEle!=0) {
				FluxRateNight[j]->SetPoint(npNightEle,B2meanNight,FluxNightEle);
				FluxRateNight[j]->SetPointError(npNightEle,0,0,FluxNightEleErr,FluxNightEleErr);
				}	
				npDayEle++;
				npNightEle++;
			} // if j < 5
		
			if (j>5) 
				{			
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			//RateDay
			RateDayPos = P0BinDay[j][k]->GetEntries(); 
			RateDayPosErr = TMath::Sqrt(P0BinDay[j][k]->GetEntries());
			FluxDayPos = RateDayPos/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidthDay);
			FluxDayPosErr = RateDayPosErr/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidthDay);
	          cout << "j = " << k << ", k = " << k << " B2meanDay " << B2meanDay << ", FluxDayPos = " << FluxDayPos << endl;
			if(B2meanDay!=0 && FluxDayPos!=0) {
				FluxRateDay[j]->SetPoint(npDayPos,B2meanDay,FluxDayPos);
				FluxRateDay[j]->SetPointError(npDayPos,0,0,FluxDayPosErr,FluxDayPosErr);
				}
			//Rate Night
			RateNightPos = P0BinNight[j][k]->GetEntries(); 
			RateNightPosErr = TMath::Sqrt(P0BinNight[j][k]->GetEntries());
			FluxNightPos = RateNightPos/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidthNight);
			FluxNightPosErr = RateNightPosErr/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidthNight);
          cout << "j = " << k << ", k = " << k << " B2meanNight " << B2meanNight << ", FluxNightPos = " << FluxNightPos << endl;

			if(B2meanNight!=0 && FluxNightPos!=0) {
				FluxRateNight[j]->SetPoint(npNightPos,B2meanNight,FluxNightPos);
				FluxRateNight[j]->SetPointError(npNightPos,0,0,FluxNightPosErr,FluxNightPosErr);
				}	
				npDayPos++;
				npNightPos++;
			} // if j > 5		
		}
    }
	
	
	TGraphAsymmErrors**FluxRateAllEle = new TGraphAsymmErrors*[NPBins];
		for (int j=0; j<NPBinsAll-1;j++) 
	  {	
		FluxRateAllEle[j] = new TGraphAsymmErrors();
	    FluxRateAllEle[j]->SetMarkerStyle(kOpenCircle);
		FluxRateAllEle[j]->SetMarkerSize(0.8);
		FluxRateAllEle[j]->SetMarkerColor(kBlack);
		int npAll = 1;
		for (int k=0;k<NPressureBins;k++) 
		{
			double Rate;
			double Flux;
			double RateErr;
			double FluxErr;
			double time = Tstart + (k+1)/(binsperhour);
			double timeWidth = durationAll[k]*60*60;		       	//width of time bin in seconds
			double MomWidth;
		    double B2mean = B2All[j][k]->GetMean();
			MomWidth = TMath::Abs(PBinsAll[j] - PBinsAll[j+1]); //width of momentum bin in MeV
			//RateDay
			Rate = P0BinAll[j][k]->GetEntries(); 
			RateErr = TMath::Sqrt(P0BinAll[j][k]->GetEntries());
			Flux = Rate/(1e-4*GeoFactor[j]*1e-3*MomWidth*timeWidth);
			FluxErr = RateErr/(1e-4*GeoFactor[j]*1e-3*MomWidth*timeWidth);
			if(B2mean!=0 && Flux!=0) {
				FluxRateAllEle[j]->SetPoint(npAll,B2mean,Flux);
				FluxRateAllEle[j]->SetPointError(npAll,0,0,FluxErr,FluxErr);
				if (j==2) totxtAll << B2mean << " "  << Flux << " " << FluxErr << endl;
				npAll++;
				}
		}
	}
	
	//Save points for 150-300MeV bin in txt file
	int npointsNight = FluxRateNight[2]->GetN();
	int npointsDay = FluxRateDay[2]->GetN();	
	
	for(int i=0;i<npointsDay;i++) 
		{
		double B2Ele, B2Pos;
		double FluxEle, FluxPos;
		double FluxErrEle, FluxErrPos;
		FluxRateDay[2]->GetPoint(i,B2Ele,FluxEle);
		FluxErrEle = FluxRateDay[2]->GetErrorY(i);
		if(FluxErrEle!=0)
		{
		cout << "DAY i = " << i << " B2Ele = " << B2Ele << " FluxEle = " << FluxEle << " FluxErrEle = " << FluxErrEle << endl;
		totxtEle << B2Ele << " "  << FluxEle << " " << FluxErrEle << endl;
		}
		FluxRateDay[8]->GetPoint(i,B2Pos,FluxPos);
		FluxErrPos = FluxRateDay[8]->GetErrorY(i);
		if(FluxErrPos!=0) 
		{
		cout << "DAY i = " << i << " B2Pos = " << B2Pos << " FluxPos = " << FluxPos << " FluxErrPos = " << FluxErrPos << endl;
		totxtPos << B2Pos << " "  << FluxPos << " " << FluxErrPos << endl;
		}
		}
		for(int i=0;i<npointsNight;i++) 
		{
		double B2Ele, B2Pos;
		double FluxEle, FluxPos;
		double FluxErrEle, FluxErrPos;
		FluxRateNight[2]->GetPoint(i,B2Ele,FluxEle);
		FluxErrEle = FluxRateNight[2]->GetErrorY(i);
		if(FluxErrEle!=0)
		{
		cout << "NIGHT i = " << i << " B2Ele = " << B2Ele << " FluxEle = " << FluxEle << " FluxErrEle = " << FluxErrEle << endl;
		totxtEle << B2Ele << " "  << FluxEle << " " << FluxErrEle << endl;
		}
		FluxRateNight[8]->GetPoint(i,B2Pos,FluxPos);
		FluxErrPos = FluxRateNight[8]->GetErrorY(i);
		if(FluxErrPos!=0) 
		{
		cout << "NIGHT i = " << i << " B2Ele = " << B2Ele << " FluxEle = " << FluxEle << " FluxErrEle = " << FluxErrEle << endl;
		totxtPos << B2Pos << " "  << FluxPos << " " << FluxErrPos << endl;
		 }	
		}
	
 
	
	//Define Daniel & Stephens functions
	TF1**DS = new TF1*[7];
	for (int k=0;k<7;k++) 
	{
		DS[k] = new TF1(Form("ds%d",k),DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DS[k]->FixParameter(l,DScoeff[k][l]);
		DS[k]->SetParameter(7,1);
		DS[k]->SetLineColor(kBlack);
		DS[k]->SetLineStyle(6);

	}
	
	//////////////////////////////////////////////////
	/////////////////DRAW CANVASES////////////////////
	//////////////////////////////////////////////////
	
	gStyle->SetTextFont(85);
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1);	

	
	//First growth curve subplot
	TCanvas *c0 = new TCanvas("c0","c0",1000,1000);
	TMultiGraph*GrowthCurves0 = new TMultiGraph();
	TLegend *leg0 = new TLegend(0.55,0.70,0.9,0.85);
	TPaveText *pt0 = new TPaveText(0.75,0.65,0.9,0.75);
	leg0->SetBorderSize(0);
	leg0->SetFillStyle(0);
	pt0->SetBorderSize(0);
	pt0->SetFillStyle(0);
    GrowthCurves0->Add(FluxRateNight[4]);
    GrowthCurves0->Add(FluxRateNight[6]);
    GrowthCurves0->Add(FluxRateDay[4]);
    GrowthCurves0->Add(FluxRateDay[6]);
    GrowthCurves0->Add(FluxRateAllEle[0]);
	leg0->AddEntry(FluxRateAllEle[0],"e^{-} + e^{+}  10-30 MeV","p");	
	leg0->AddEntry(FluxRateDay[4],"e^{-} 10-30 MeV","p");
	leg0->AddEntry(FluxRateDay[6],"e^{+} 10-30 MeV","p");
	leg0->AddEntry(DS[0],"Daniel & Stephens 20 MeV","l");
	c0->cd();
	GrowthCurves0->Draw("AP");
	FluxRateAllEle[0]->Fit("ds0","","",50,1000);
	pt0->AddText(Form("Normalization factor %2.2f",DS[0]->GetParameter(7)));
	DS[0]->DrawF1(1,1000,"same");
	leg0->Draw("same");
    pt0->Draw();
	GrowthCurves0->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves0->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves0->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves0->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves0->GetXaxis()->CenterTitle();	
	GrowthCurves0->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves0->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves0->GetYaxis()->CenterTitle();
	GrowthCurves0->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves0->SetMinimum(ZMin);
	GrowthCurves0->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Second subplot
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TMultiGraph*GrowthCurves1 = new TMultiGraph();
	TLegend *leg1 = new TLegend(0.55,0.70,0.9,0.85);
	TPaveText *pt1 = new TPaveText(0.75,0.65,0.9,0.75,"NDC");
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	pt1->SetBorderSize(0);
	pt1->SetFillStyle(0);
    GrowthCurves1->Add(FluxRateNight[3]);
    GrowthCurves1->Add(FluxRateNight[7]);
    GrowthCurves1->Add(FluxRateDay[3]);
    GrowthCurves1->Add(FluxRateDay[7]);
    GrowthCurves1->Add(FluxRateAllEle[1]);
	leg1->AddEntry(FluxRateAllEle[1],"e^{-} + e^{+} 30-60 MeV","p");
	leg1->AddEntry(FluxRateDay[3],"e^{-} 30-60 MeV","p");
	leg1->AddEntry(FluxRateDay[7],"e^{+} 30-60 MeV","p");
	leg1->AddEntry(DS[2],"Daniel & Stephens 45 MeV","l");
	c1->cd();
	GrowthCurves1->Draw("AP");		
	FluxRateAllEle[1]->Fit("ds2","","",50,1000);
	DS[2]->DrawF1(1,1000,"same");
	pt1->AddText(Form("Normalization factor %2.2f",DS[2]->GetParameter(7)));
	leg1->Draw("same");
    pt1->Draw();	
	GrowthCurves1->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves1->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves1->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves1->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves1->GetXaxis()->CenterTitle();	
	GrowthCurves1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves1->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves1->GetYaxis()->CenterTitle();
	GrowthCurves1->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves1->SetMinimum(ZMin);
	GrowthCurves1->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Third subplot
	TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
	TMultiGraph*GrowthCurves2 = new TMultiGraph();
	TLegend *leg2 = new TLegend(0.55,0.70,0.9,0.85);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	TPaveText *pt2 = new TPaveText(0.75,0.65,0.9,0.75);
	pt2->SetBorderSize(0);
	pt2->SetFillStyle(0);
    GrowthCurves2->Add(FluxRateNight[2]);
    GrowthCurves2->Add(FluxRateNight[8]);
    GrowthCurves2->Add(FluxRateDay[2]);
    GrowthCurves2->Add(FluxRateDay[8]);
    GrowthCurves2->Add(FluxRateAllEle[2]);
	leg2->AddEntry(FluxRateAllEle[2],"e^{-} + e^{+} 60-80 MeV","p");
	leg2->AddEntry(FluxRateDay[2],"e^{-} 60-80 MeV","p");
	leg2->AddEntry(FluxRateDay[8],"e^{+} 60-80 MeV","p");
	leg2->AddEntry(DS[3],"Daniel & Stephens 70 MeV","l");
	c2->cd();
	GrowthCurves2->Draw("AP");	
	FluxRateAllEle[2]->Fit("ds3","","",50,1000);
	DS[3]->DrawF1(1,1000,"same");
	pt2->AddText(Form("Normalization factor %2.2f",DS[3]->GetParameter(7)));
	leg2->Draw("same");
    pt2->Draw();
	GrowthCurves2->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves2->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves2->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves2->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves2->GetXaxis()->CenterTitle();	
	GrowthCurves2->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves2->GetYaxis()->CenterTitle();
	GrowthCurves2->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves2->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves2->SetMinimum(ZMin);
	GrowthCurves2->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fourth subplot
	TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
	TMultiGraph*GrowthCurves3 = new TMultiGraph();
	TLegend *leg3 = new TLegend(0.55,0.70,0.9,0.85);
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);	
	TPaveText *pt3 = new TPaveText(0.75,0.65,0.9,0.75,"NDC");
	pt3->SetBorderSize(0);
	pt3->SetFillStyle(0);
    GrowthCurves3->Add(FluxRateNight[1]);
    GrowthCurves3->Add(FluxRateNight[9]);
    GrowthCurves3->Add(FluxRateDay[1]);
    GrowthCurves3->Add(FluxRateDay[9]);
    GrowthCurves3->Add(FluxRateAllEle[3]);
	leg3->AddEntry(FluxRateAllEle[3],"e^{-} + e^{+} 80-150 MeV","p");
	leg3->AddEntry(FluxRateNight[1],"e^{-} 80-150 MeV","p");
	leg3->AddEntry(FluxRateNight[9],"e^{+} 80-150 MeV","p");
	leg3->AddEntry(DS[4],"Daniel & Stephens 111 MeV","l");
	c3->cd();
	GrowthCurves3->Draw("AP");	
	FluxRateAllEle[3]->Fit("ds4","","",50,1000);
	DS[4]->DrawF1(1,1000,"same");
	pt3->AddText(Form("Normalization factor %2.2f",DS[4]->GetParameter(7)));
	leg3->Draw("same");
    pt3->Draw("same");
	GrowthCurves3->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves3->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves3->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves3->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves3->GetXaxis()->CenterTitle();	
	GrowthCurves3->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves3->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves3->GetYaxis()->CenterTitle();
	GrowthCurves3->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves3->SetMinimum(ZMin);
	GrowthCurves3->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fifth and last subplot
	TCanvas *c4 = new TCanvas("c4","c4",1000,1000);
	TMultiGraph*GrowthCurves4 = new TMultiGraph();
	TLegend *leg4 = new TLegend(0.55,0.70,0.9,0.85);
	leg4->SetBorderSize(0);
	leg4->SetFillStyle(0);
	TPaveText *pt4 = new TPaveText(2,20000,30,30000);
	pt4->SetBorderSize(0);
	pt4->SetFillStyle(0);
    GrowthCurves4->Add(FluxRateNight[0]);
    GrowthCurves4->Add(FluxRateNight[10]);
    GrowthCurves4->Add(FluxRateDay[0]);
    GrowthCurves4->Add(FluxRateDay[10]);
    GrowthCurves4->Add(FluxRateAllEle[4]);
	leg4->AddEntry(FluxRateAllEle[4],"e^{-} + e^{+} 150-300 MeV","p");
	leg4->AddEntry(FluxRateDay[0],"e^{-} 150-300 MeV","p");
	leg4->AddEntry(FluxRateDay[10],"e^{+} 150-300 MeV","p");
	leg4->AddEntry(DS[6],"Daniel & Stephens 279 MeV","l");
	c4->cd();
	GrowthCurves4->Draw("AP");	
	FluxRateAllEle[4]->Fit("ds6","","",10,700);
	DS[6]->DrawF1(1,1000,"same");
	pt4->AddText(Form("Normalization factor %2.2f",DS[5]->GetParameter(7)));
	cout << "p7 = " << DS[5]->GetParameter(7) << endl;
	leg4->Draw("same");
    pt4->Draw();
	GrowthCurves4->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves4->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves4->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves4->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves4->GetXaxis()->CenterTitle();	
	GrowthCurves4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves4->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves4->GetXaxis()->SetLimits(PMin,PMax);
	GrowthCurves4->GetYaxis()->CenterTitle();
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves4->SetMinimum(ZMin);
	GrowthCurves4->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();
	

	
	c0->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c1->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c2->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c3->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c4->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));




	//Write histograms in ROOT file	
	fileout->cd();	
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

	cout << "now closing all files" << endl;
	//fileout->Close();
	//fileCal->Close(); 	
			
} //end function		
		
void GrowthCurvesSignedPressureBinsAllFlight(string RecoID, double cutT2, double cutT3) {

	double Tstart = 0;
	double Tend = 134;	//total time of flight in hours)
	double Launch = 22 + (7./60.);
	double Launchday = 15;
	double LaunchDay = 15 + (Launch*0.041667);
	double binsperhour = 2;
	int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
	const int NPressureBins = 19;
	double 	PressureBins[NPressureBins] = {1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.0,3.25,3.50,5.0,10,30,70,100,300,600,1000};
	const int NPBins = 12;
	double PBins[NPBins] = {-300,-150,-80,-60,-30,-10,10,30,60,80,150,300};
	double GeoFactor[5] = {2.68,10.43,11.55,12.02,14.64};
	//for the ascent all data is in file 000
	//for the final leg of the flight (no cut-off), data starts at file 15
	int StartFile = 0;
	int EndFile = 23;
	double TCut[7] = {19.8,30,42,56,65,82,90};
	int mcolor[7]={1,2,3,4,800+1,6,7};
    double ELoss = 4;				//4 MeV energy loss in shell + T1T2T3

	string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2/";
	string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves";
	string CutConfig = "T2Above120T3Less110";
	//string CutConfig="testfit";

	//Axis limit for plots
	float ZMax=10000000;		//growth curve y-axis
	float ZMin=100;			//growth curve y-axis
	float PMin=1;			//growth curve x-axis
	float PMax=1000;			//growth curve x-axis
	float TMin=0;			//time serie x-axis
	float TMax=133;			//time serie x-axis	
	
	
	//Initialize array
	string DSConfig[7] = {"20 MeV","30 MeV", "44.8 MeV", "70 MeV", "111 MeV", "177.5 MeV","279 MeV"};
	float**DScoeff=new float*[7];
	for(int i=0; i<7;i++) 
	 {
	 DScoeff[i] =  new float[7];
	 for(int j=0; j<7;j++)			//7 coefficients for each array 
	  {
		 DScoeff[i][j]=0.;
	 }
    }
	 string DSparamfile="./DScoefficients.dat"; 
	 LoadDSparameters(DSparamfile,DScoeff);	

	double T2Cut = cutT2;	//value of cut on T2
	double T3Cut = cutT3;	//value of cut on T3
	TFile *fileRaw;
	TFile *fileCal;
	string ID = RecoID;	

	//Create histograms
	TH1F***B2 = new TH1F**[NPBins];
	TH1F***P0Bin = new TH1F**[NPBins];


	//Momentum histograms and pressure
	for(int i=0; i<NPBins; i++) 
	{
			P0Bin[i] = new TH1F*[NPressureBins];
			B2[i] = new TH1F*[NPressureBins];
			for(int j=0; j<NPressureBins; j++) 
			 {
				P0Bin[i][j] = new TH1F(Form("unsigned p0 h%d%d",i,j), Form("unsigned p0 h%d%d",i,j), 1,0,1);
                B2[i][j] = new TH1F(Form("B2 h%d%d",i,j),Form("B2 h%d%d",i,j), 10000, 0, 1000);
			}
		}

	
   //Determine the duration payload spends in each pressure bin
	
	double duration[NPressureBins];
	for(int k=0;k<NPressureBins;k++) 
	{
		duration[k]=0;
	}
	float startperiod;			//initialize timefirst event
	int prevbinlayer=17;			//pressure bin of first even
	float prevtimelayer=0.0058328;		//first positive time measurement made
    float tmptimePHA;	
	//ROOT fileout
	 //TFile *fileout = new TFile(Form("%s/GrowthCurvesAllElectrons_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
     for (int j=StartFile; j<EndFile; j++)
	 {
	 	if (j<10) fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 	if (j>=10) fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 
	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	
	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
     {
		treeCal->GetEntry(i);
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
	    tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
		if (tmptimePHA<0) continue;
	//Barometer data (calibrated)
		float B2Pres = de->get_PressB2();
	//Convert pressure to g/cm^-2
		B2Pres = B2Pres* 1.3595;
		//cout << "t = " << tmptimePHA <<  ", B2Pres = " << B2Pres << " g.cm^-2 " << endl;
	//Determine the correct time, energy bin and pressure bin
		int pressurebin;
		for(int k=0;k<NPressureBins;k++) {
			if((PressureBins[k] <= B2Pres) & (B2Pres < PressureBins[k+1]))
			{
				pressurebin = k;
				break;
			}	
		}
			if(pressurebin != prevbinlayer) {
								duration[prevbinlayer]+= tmptimePHA - startperiod;
							}
				startperiod = tmptimePHA;
				prevbinlayer = pressurebin;
			

	} // end i
		 //Fill in last bin
 		duration[prevbinlayer]+= tmptimePHA - startperiod;
		 startperiod = tmptimePHA;
 } //end j
		 double sum=0;
		 for(int k=0;k<NPressureBins;k++) {
			cout << " k = " << k << " duration[i] = " << duration[k] << endl;
		     sum+=duration[k];
		 }
		 cout << "sum = " << sum << endl;
	
	//ROOT fileout
	TFile *fileout = new TFile(Form("%s/GrowthCurvesSignedPressureBinsAllFlight_%s_%s.root", Outpath.c_str(),ID.c_str(),CutConfig.c_str()),"RECREATE");
     for (int j=StartFile; j<EndFile; j++)
	 {
	 	if (j<10) fileCal=new TFile(Form("%s/%s/18A1_00%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 	if (j>=10) fileCal=new TFile(Form("%s/%s/18A1_0%d.BPD.EVENT_%s_Calibrated.root",Inpath.c_str(),ID.c_str(),j,ID.c_str()),"READ");
	 
	//Get Tree from the calibrated file	
	TTree *treeCal = (TTree*)fileCal->Get("Data");
	ALEvent *de = new ALEvent();
	treeCal->SetBranchAddress("Calevent",&de);  
	int nentries=treeCal->GetEntries();	
	cout << "Number  of events: " << nentries << endl;  

	for(int i=0; i<nentries; i++) 
     {
		
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
	//Time of event from PHA
		double yPHA = (double)de->get_yPHA();
		double mPHA = (double)de->get_mPHA();
		double dPHA = (double)de->get_dPHA();
		double hPHA= (double)de->get_hPHA();
		double miPHA = (double)de->get_miPHA();
		double sPHA = (double)de->get_sPHA();
		float tmptimePHA = hPHA + (miPHA/60) + (sPHA/3600);
		tmptimePHA = tmptimePHA-Launch;
		tmptimePHA = tmptimePHA + 24*(dPHA - Launchday);
	 // cout << "tmptimeCT1 = " << tmptime << ", tmptimePHA = " << tmptimePHA << endl; 
		
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
		//cout << "B2Pres = " << B2Pres << "g.cm^-2 " << endl;
	//Reconstructed energy 
		float pPR=1000*fabs(de->get_p0PR());   //in MeV
		double deflecPR = de->get_deflecPR();
		float preco=de->get_p0reco();
	    float P0sign;
		float ePR = 1000*(de->get_EkPR());		//total energy in MeV
		if (strcmp(ID.c_str(),"PRonly")==0)  P0sign = pPR*TMath::Sign(1,deflecPR);
		else P0sign = preco*TMath::Sign(1,deflecPR);	
		double chiNBPR = de->get_chi2NBPR();
		double chiBPR = de->get_chi2BPR();	
	    pPR = pPR + ELoss;					//correct for energy loss in shell + T1 + T2 + T3
	    ePR = ePR + ELoss;					//correct for energy loss in shell + T1 + T2 + T3
		float ePRsign = ePR*TMath::Sign(1,deflecPR);
	//Determine the correct time, energy bin and pressure bin
		int tbin = tmptimePHA*binsperhour;
		int pbin;
		int pressurebin;
		for(int k=0;k<NPBins;k++) {if((PBins[k] <= ePRsign) & (ePRsign < PBins[k+1])) pbin = k;}
		for(int k=0;k<NPressureBins;k++) {if((PressureBins[k] <= B2Pres) & (B2Pres < PressureBins[k+1])) pressurebin = k;}			

	//Apply cuts to only get clean events
		if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -300) & (P0sign <= 300) & (!G) & (nhits >=5) && (nhits <15))
		   {
	      //  cout << "t = " << tmptimePHA << ", ePR = " << ePR << " MeV, at B2 = " << B2Pres << ", tbin = " << tbin <<", pressurebin = " << pressurebin << endl;
			P0Bin[pbin][pressurebin]->Fill(ePRsign);
		    B2[pbin][pressurebin]->Fill(B2Pres);
			}	   

	}     //end i, loop on events
}		//end j, loop on files

	
	//Make growth curves
	int ngraphs = (NPBins/2)-1;
    TGraphAsymmErrors**FluxRate = new TGraphAsymmErrors*[NPBins];


	for (int j=0; j<NPBins-1;j++) 
	{

		FluxRate[j] = new TGraphAsymmErrors();
	    FluxRate[j]->SetMarkerStyle(kOpenCircle);
		FluxRate[j]->SetMarkerSize(0.8);	

//Electrons
		if (j<5)  { 
		FluxRate[j]->SetMarkerColor(kBlue);
		FluxRate[j]->SetLineColor(kBlue);
		}
//Positrons	
		if (j>5)  { 
		FluxRate[j]->SetMarkerColor(kRed);
		FluxRate[j]->SetLineColor(kRed);
		}
		int npPos = 1;
		int npEle = 1;
		for (int k=0;k<NPressureBins;k++) 
		{
			double RateEle,RatePos;
			double FluxEle,FluxPos;
			double RateEleErr,RatePosErr;
			double FluxEleErr,FluxPosErr;
			double time = Tstart + (k+1)/(binsperhour);
			double timeWidth= duration[k]*60*60;		       	//width of time bin in seconds
			double MomWidth;
		    double B2mean = B2[j][k]->GetMean();			
			
			if (j<5) 
				{			
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			RateEle = P0Bin[j][k]->GetEntries(); 
			RateEleErr = TMath::Sqrt(P0Bin[j][k]->GetEntries());
			FluxEle = RateEle/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidth);
			FluxEleErr = RateEleErr/(1e-4*GeoFactor[4-j]*1e-3*MomWidth*timeWidth);
            cout << "j = " << k << ", k = " << k << " RateEle " << RateEle << ", FluxEle = " << FluxEle << endl;
			if(B2mean!=0 && FluxEle!=0) {
				FluxRate[j]->SetPoint(npEle,B2mean,FluxEle);
				FluxRate[j]->SetPointError(npEle,0,0,FluxEleErr,FluxEleErr);
				}
			npEle++;
			} // if j < 5
		
			if (j>5) 
				{			
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			RatePos = P0Bin[j][k]->GetEntries(); 
			RatePosErr = TMath::Sqrt(P0Bin[j][k]->GetEntries());
			FluxPos = RatePos/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidth);
			FluxPosErr = RatePosErr/(1e-4*GeoFactor[j-6]*1e-3*MomWidth*timeWidth);
	        cout << "j = " << k << ", k = " << k << " RatePos " << RatePos << ", FluxPos = " << FluxPos << endl;
			if(B2mean!=0 && FluxPos!=0) {
				FluxRate[j]->SetPoint(npPos,B2mean,FluxPos);
				FluxRate[j]->SetPointError(npPos,0,0,FluxPosErr,FluxPosErr);
				}
				npPos++;
			} // if j > 5		
		}
    }
	
	//Define Daniel & Stephens functions
	TF1**DS = new TF1*[7];
	for (int k=0;k<7;k++) 
	{
		DS[k] = new TF1(Form("ds%d",k),DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DS[k]->FixParameter(l,DScoeff[k][l]);
		DS[k]->SetParameter(7,1);
		DS[k]->SetLineColor(kBlack);
		DS[k]->SetLineStyle(6);

	}
	
	//////////////////////////////////////////////////
	/////////////////DRAW CANVASES////////////////////
	//////////////////////////////////////////////////
	
	gStyle->SetTextFont(85);
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1);	

	
	//First growth curve subplot
	TCanvas *c0 = new TCanvas("c0","c0",1000,1000);
	TMultiGraph*GrowthCurves0 = new TMultiGraph();
	TLegend *leg0 = new TLegend(0.55,0.70,0.9,0.85);
	TPaveText *pt0 = new TPaveText(0.75,0.65,0.9,0.75);
	leg0->SetBorderSize(0);
	leg0->SetFillStyle(0);
	pt0->SetBorderSize(0);
	pt0->SetFillStyle(0);
    GrowthCurves0->Add(FluxRate[4]);
    GrowthCurves0->Add(FluxRate[6]);
	leg0->AddEntry(FluxRate[4],"e^{-} 10-30 MeV","p");
	leg0->AddEntry(FluxRate[6],"e^{+} 10-30 MeV","p");
	leg0->AddEntry(DS[0],"Daniel & Stephens 20 MeV","l");
	c0->cd();
	GrowthCurves0->Draw("AP");
	FluxRate[4]->Fit("ds0","","",50,1000);
	pt0->AddText(Form("Normalization factor %2.2f",DS[0]->GetParameter(7)));
	DS[0]->DrawF1(1,1000,"same");
	leg0->Draw("same");
    pt0->Draw();
	GrowthCurves0->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves0->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves0->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves0->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves0->GetXaxis()->CenterTitle();
	GrowthCurves0->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves0->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves0->GetYaxis()->CenterTitle();
	GrowthCurves0->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves0->SetMinimum(ZMin);
	GrowthCurves0->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Second subplot
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TMultiGraph*GrowthCurves1 = new TMultiGraph();
	TLegend *leg1 = new TLegend(0.55,0.70,0.9,0.85);
	TPaveText *pt1 = new TPaveText(0.75,0.65,0.9,0.75,"NDC");
	leg1->SetBorderSize(0);
	leg1->SetFillStyle(0);
	pt1->SetBorderSize(0);
	pt1->SetFillStyle(0);
    GrowthCurves1->Add(FluxRate[3]);
    GrowthCurves1->Add(FluxRate[7]);
	leg1->AddEntry(FluxRate[3],"e^{-} 30-60 MeV","p");
	leg1->AddEntry(FluxRate[7],"e^{+} 30-60 MeV","p");
	leg1->AddEntry(DS[2],"Daniel & Stephens 45 MeV","l");
	c1->cd();
	GrowthCurves1->Draw("AP");		
	FluxRate[3]->Fit("ds2","","",50,1000);
	DS[2]->DrawF1(1,1000,"same");
	pt1->AddText(Form("Normalization factor %2.2f",DS[2]->GetParameter(7)));
	leg1->Draw("same");
    pt1->Draw();	
	GrowthCurves1->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves1->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves1->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves1->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves1->GetXaxis()->CenterTitle();
	GrowthCurves1->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves1->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves1->GetYaxis()->CenterTitle();
	GrowthCurves1->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves1->SetMinimum(ZMin);
	GrowthCurves1->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update(); 

	
	//Third subplot
	TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
	TMultiGraph*GrowthCurves2 = new TMultiGraph();
	TLegend *leg2 = new TLegend(0.55,0.70,0.9,0.85);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	TPaveText *pt2 = new TPaveText(0.75,0.65,0.9,0.75);
	pt2->SetBorderSize(0);
	pt2->SetFillStyle(0);
    GrowthCurves2->Add(FluxRate[2]);
    GrowthCurves2->Add(FluxRate[8]);
	leg2->AddEntry(FluxRate[2],"e^{-} 60-80 MeV","p");
	leg2->AddEntry(FluxRate[8],"e^{+} 60-80 MeV","p");
	leg2->AddEntry(DS[3],"Daniel & Stephens 70 MeV","l");
	c2->cd();
	GrowthCurves2->Draw("AP");	
	FluxRate[2]->Fit("ds3","","",50,1000);
	DS[3]->DrawF1(1,1000,"same");
	pt2->AddText(Form("Normalization factor %2.2f",DS[3]->GetParameter(7)));
	leg2->Draw("same");
    pt2->Draw();
	GrowthCurves2->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves2->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves2->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves2->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves2->GetXaxis()->CenterTitle();
	GrowthCurves2->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves2->GetYaxis()->CenterTitle();
	GrowthCurves2->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves2->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves2->SetMinimum(ZMin);
	GrowthCurves2->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fourth subplot
	TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
	TMultiGraph*GrowthCurves3 = new TMultiGraph();
	TLegend *leg3 = new TLegend(0.55,0.70,0.9,0.85);
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);	
	TPaveText *pt3 = new TPaveText(0.75,0.65,0.9,0.75,"NDC");
	pt3->SetBorderSize(0);
	pt3->SetFillStyle(0);
    GrowthCurves3->Add(FluxRate[1]);
    GrowthCurves3->Add(FluxRate[9]);
	leg3->AddEntry(FluxRate[1],"e^{-} 80-150 MeV","p");
	leg3->AddEntry(FluxRate[9],"e^{+} 80-150 MeV","p");
	leg3->AddEntry(DS[4],"Daniel & Stephens 111 MeV","l");
	c3->cd();
	GrowthCurves3->Draw("AP");	
	FluxRate[1]->Fit("ds4","","",50,1000);
	DS[4]->DrawF1(1,1000,"same");
	pt3->AddText(Form("Normalization factor %2.2f",DS[4]->GetParameter(7)));
	leg3->Draw("same");
    pt3->Draw("same");
	GrowthCurves3->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves3->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves3->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves3->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves3->GetXaxis()->CenterTitle();
	GrowthCurves3->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves3->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves3->GetYaxis()->CenterTitle();
	GrowthCurves3->GetXaxis()->SetLimits(PMin,PMax);
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves3->SetMinimum(ZMin);
	GrowthCurves3->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();	
	
	//Fifth and last subplot
	TCanvas *c4 = new TCanvas("c4","c4",1000,1000);
	TMultiGraph*GrowthCurves4 = new TMultiGraph();
	TLegend *leg4 = new TLegend(0.55,0.70,0.9,0.85);
	leg4->SetBorderSize(0);
	leg4->SetFillStyle(0);
	TPaveText *pt4 = new TPaveText(2,20000,30,30000);
	pt4->SetBorderSize(0);
	pt4->SetFillStyle(0);
    GrowthCurves4->Add(FluxRate[0]);
    GrowthCurves4->Add(FluxRate[10]);
	leg4->AddEntry(FluxRate[0],"e^{-} 150-300 MeV","p");
	leg4->AddEntry(FluxRate[10],"e^{+} 150-300 MeV","p");
	leg4->AddEntry(DS[5],"Daniel & Stephens 177.5 MeV","l");
	c4->cd();
	GrowthCurves4->Draw("AP");	
	FluxRate[0]->Fit("ds5","","",50,1000);
	DS[5]->DrawF1(1,1000,"same");
	pt4->AddText(Form("Normalization factor %2.2f",DS[5]->GetParameter(7)));
	cout << "p7 = " << DS[5]->GetParameter(7) << endl;
	leg4->Draw("same");
    pt4->Draw();
	GrowthCurves4->GetXaxis()->SetTitleOffset(0.87);
	GrowthCurves4->GetYaxis()->SetTitleOffset(0.87);
	GrowthCurves4->GetXaxis()->SetTitleSize(0.05);
	GrowthCurves4->GetYaxis()->SetTitleSize(0.05);
	GrowthCurves4->GetXaxis()->CenterTitle();
	GrowthCurves4->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
	GrowthCurves4->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
	GrowthCurves4->GetXaxis()->SetLimits(PMin,PMax);
	GrowthCurves4->GetYaxis()->CenterTitle();
	gPad->SetLogy();
	gPad->SetLogx();
	GrowthCurves4->SetMinimum(ZMin);
	GrowthCurves4->SetMaximum(ZMax);
	gPad->SetLeftMargin(0.15);
	gPad->Update();
	

	
	c0->Print(Form("%s/GrowthCurvesBinnedSignedAllFlight_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c1->Print(Form("%s/GrowthCurvesBinnedSignedAllFlight_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c2->Print(Form("%s/GrowthCurvesBinnedSignedAllFlight_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c3->Print(Form("%s/GrowthCurvesBinnedSignedAllFlight_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
	c4->Print(Form("%s/GrowthCurvesBinnedSignedAllFlight_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));




	//Write histograms in ROOT file	
	fileout->cd();	
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

	cout << "now closing all files" << endl;
	//fileout->Close();
	//fileCal->Close(); 	
			
} //end function	
	
	
	
	
	
	
	
	
	
	
