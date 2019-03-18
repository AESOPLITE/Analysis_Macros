////////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Sarah Mechbal, smechbal@ucsc.edu
///    Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, July 3rd, 2018
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "ALEvent.h"
#include "TMath.h"
#include "TGaxis.h"
#include "LoadMCGrowthCurves.h"
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
	
	return coeff[7]*TMath::Exp(z8);		
}

TF1* fitele;
TF1* fitpos;
TF1* DSpos = new TF1("dspos",DanAndSteph,1,1000,8);
TF1* DSele = new TF1("dsele",DanAndSteph,1,1000,8);

double sumMC(Double_t *x, Double_t *par) {
	return fitele->Eval(x[0]) + fitpos->Eval(x[0]);
	}
double sumData(Double_t *x, Double_t *par) {
	return DSpos->Eval(x[0]) + DSele->Eval(x[0]);
	}
void FluxCurves(string RecoID, double cutT2, double cutT3) {

	double Tstart = 0;
	double Tend = 134;	//total time of flight in hours)
	double Launch = 22 + (7./60.);
	double Launchday = 15;
	double LaunchDay = 15 + (Launch*0.041667);
	double binsperhour = 2;
	int NTimeBin = floor((Tend - Tstart)*(binsperhour)) + 1;
	const int NPressureBins = 18;
	double 	PressureBins[NPressureBins] = {1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.0,3.25,3.50,5.0,10,30,70,100,300,1000};
	
	
	//Define Energy Bins
	const int NPBinsAll = 13;
	double PBinsAll[NPBinsAll] = {0,10,20,30,40,50,60,80,100,200,300,400,500};	
	double PBinsWidth[NPBinsAll-1] = {10,10,10,10,20,20,50,50,50,50,100,100};	
	const int NPBins = NPBinsAll*2;
	double PBins[NPBins] = {-500,-400,-300,-200,-100,-80,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,80,100,200,300,400,500};	
	double*MeanEneEle= new double[NPBinsAll-1];
	double*MeanEnePos= new double[NPBinsAll-1];
	double*MeanEneAll= new double[NPBinsAll-1];
	//double GeoFactor[NPBinsAll-1] = {0.001,2.68,6.24,10.43,11.11,11.55,12.02,12.02,12.02,14.64,14.64,14.64};					//old values
	double GeoFactor[NPBinsAll-1] = {0.0001,1.07,3.53,6.60,7.36,7.93,8.72,8.67,10.64,11.18,11.18,11.57};					//new values 
	float FloatDepth[7] = {2.16,2.39,2.65,2.86,3.16,3.35,3.67};
	for(int i=0;i<NPBinsAll-1;i++) 
	{
		MeanEneEle[i]=0.;
		MeanEnePos[i]=0.;
		MeanEneAll[i]=0.;
	}
	
	//for the ascent all data is in file 000
	//for the final leg of the flight (no cut-off), data starts at file 15
	int StartFile = 0;
	int EndFile = 23;
	double TCut[7] = {19.8,30,42,56,65,82,90};
	int mcolor[7]={1,2,3,4,800+1,6,7};
    double ELoss = 4;				//4 MeV energy loss in shell + T1T2T3

	string Inpath = "/home/sarah/AESOPLITE/FlightData/BBR2/";
	string Outpath = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves";
	string OutpathText = "/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/MCdata";
	string CutConfig = "T2Above120T3Less110";
	//string CutConfig="testfit";

	//txt output file for positrons
	ofstream totxtPos;
	totxtPos.open(Form("%s/GrowthCurvePos.txt",OutpathText.c_str()));
	totxtPos << " Pressure in g.cm-2 ";
	for(int i=0;i<NPBinsAll-1;i++) 
	{
		totxtPos << Form(" Flux e+ %3.0f MeV ",(PBinsAll[i]+PBinsAll[i+1])/2) << Form(" StdDev e- %3.0f MeV ",(PBinsAll[i]+PBinsAll[i+1])/2);
		if(i==NPBinsAll-2) totxtPos << endl;
	}	
	//txt output file for electrons
	ofstream totxtEle;
	totxtEle.open(Form("%s/GrowthCurveEle.txt",OutpathText.c_str()));
	totxtEle << " Pressure in g.cm-2 ";
	for(int i=0;i<NPBinsAll-1;i++) 
	{
		totxtEle << Form(" Flux e- %3.0f MeV ",(PBinsAll[i]+PBinsAll[i+1])/2) << Form(" StdDev e- %3.0f MeV ",(PBinsAll[i]+PBinsAll[i+1])/2);
		if(i==NPBinsAll-2) totxtEle << endl;
	}	
	ofstream totxtAll;
	totxtAll.open(Form("%s/GrowthCurveAll.txt",OutpathText.c_str()));
	totxtAll << " Pressure in g.cm-2 ";
	for(int i=0;i<NPBinsAll-1;i++) 
	{
		totxtAll << Form(" Flux e+e- %3.0f MeV ",(PBinsAll[i]+PBinsAll[i+1])/2) <<  Form(" StdDev e- %3.0f MeV ",(PBinsAll[i]+PBinsAll[i+1])/2);
		if(i==NPBinsAll-2) totxtAll << endl;
	}
	ofstream TeXTables;
	TeXTables.open(Form("%s/FluxToP_TeXTablesNorm.txt",OutpathText.c_str()));
	ofstream TeXTableFluxes[7];
	for(int i=0; i<7;i++) TeXTableFluxes[i].open(Form("%s/TeXTableFluxes_%2.2f.txt",OutpathText.c_str(),FloatDepth[i]));	


	
	//Axis limit for plots
	float ZMax=100000;		//growth curve y-axis
	float ZMin=1;			//growth curve y-axis
	float PMin=1;			//growth curve x-axis
	float PMax=3000;			//growth curve x-axis
	float TMin=0;			//time serie x-axis
	float TMax=133;			//time serie x-axis	
	
	
	//Initialize arrays
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

	//Load results of MC atmospheric simulation
	const int ndepth=24;
	const int nMCbin = 16;
	float**eMC = new float*[nMCbin];
	float**pMC = new float*[nMCbin];
	float**allMC = new float*[nMCbin];
	float*Depth = new float[ndepth];
	for(int i=0; i<nMCbin;i++) 
	 {
	 eMC[i] =  new float[ndepth];
	 pMC[i] =  new float[ndepth];
	 allMC[i] =  new float[ndepth];
	 for(int j=0; j<ndepth;j++)			//25 depths measurement for each energy bin
	  {	
		 eMC[i][j]=0.;
		 pMC[i][j]=0.;
		 allMC[i][j]=0.;
		 Depth[j]=0.;
	 }
	}
	 string MCGrowthCurvefile="/home/sarah/AESOPLITE/Analysis/FlightAnalysis/GrowthCurves/MCdata/MC_GrowthCurves_60cycles_Pcut0.1_downonly_AllZenPP_cosZinf011.txt"; 
	 LoadMCGrowthCurves(MCGrowthCurvefile,Depth,eMC,pMC);	
	for(int i=0; i<nMCbin;i++) 
	 {
	 for(int j=0; j<ndepth;j++)			//25 depths measurement for each energy bin
	  {	
		 allMC[i][j]= eMC[i][j] + pMC[i][j];
		// cout << "Ebin " << i << " Depth bin " <<j << " Depth = " << Depth[j] << ", eMC = " <<eMC[i][j] << ", pMC = " << pMC[i][j] << endl;
	 }
	}	

	
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

	//Momentum histograms and pressure for all electrons
	for(int i=0; i<NPBinsAll; i++) 
	{
			P0BinAll[i] = new TH1F*[NPressureBins];
			B2All[i] = new TH1F*[NPressureBins];
			for(int j=0; j<NPressureBins; j++) 
			 {
				P0BinAll[i][j] = new TH1F(Form("all p0 h%d%d",i,j), Form("all p0 h%d%d",i,j), 1000,0,500);
                B2All[i][j] = new TH1F(Form("B2All h%d%d",i,j),Form("B2All h%d%d",i,j), 10000, 0, 1000);
			}
		}

	//Momentum histograms and pressure for ele and pos separately
	for(int i=0; i<NPBins; i++) 
	{
			P0BinDay[i] = new TH1F*[NPressureBins];
			P0BinNight[i] = new TH1F*[NPressureBins];
			B2Day[i] = new TH1F*[NPressureBins];
			B2Night[i] = new TH1F*[NPressureBins];
			for(int j=0; j<NPressureBins; j++) 
			 {
				P0BinDay[i][j] = new TH1F(Form("unsigned p0 h%d%d",i,j), Form("unsigned p0 h%d%d",i,j), 1000,-500,500);
				P0BinNight[i][j] = new TH1F(Form("unsigned p0 night h%d%d",i,j), Form("unsigned p0 night h%d%d",i,j), 1000,-500,500);
                B2Day[i][j] = new TH1F(Form("B2Day h%d%d",i,j),Form("B2Day h%d%d",i,j), 10000, 0, 1000);
                B2Night[i][j] = new TH1F(Form("B2Night h%d%d",i,j),Form("B2Night h%d%d",i,j), 10000, 0, 1000);
			}
		}
	
	/////////////////////////////////////////////////////////////
    //Determine the duration payload spends in each pressure bin//
   //////////////////////////////////////////////////////////////
	
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
	//cout << "Number  of events: " << nentries << endl;  
	//Process all entries from file
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
	
	/////////////////////////////////////////////////////////////
    ////////////Process all files, fill histograms///////////////
   //////////////////////////////////////////////////////////////
		
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
		//if(i%1000==0) cout << "Event " << i << endl;
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
		if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -500) & (P0sign <= 500) & (!G) & (nhits >=5) && (nhits <15))
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
	if(T1 & T2 & T3 & T3PHA < T3Cut & T2PHA>T2Cut & (P0sign >= -500) & (P0sign <= 500) & (!G) & (nhits >=5) && (nhits <15))
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

	
	///////////////////////////////////////////////////////
	////////////Make Growth Curves TGraphs for DATA////////
	///////////////////////////////////////////////////////
	
	
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
		if (j<(NPBins/2)-1)  { 
		FluxRateDay[j]->SetMarkerColor(kBlue);
		FluxRateDay[j]->SetLineColor(kBlue);
		FluxRateNight[j]->SetMarkerColor(kBlue);
		FluxRateNight[j]->SetLineColor(kBlue);
		}
//Positrons	
		if (j>(NPBins/2)-2)  { 
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
			double MeanEne = P0BinNight[j][k]->GetMean();


			//cout << "j = " << k << ", k = " << k << " B2Day mean pressure " << B2meanDay << "gcm-2 " << " , duration = " << timeWidthDay << endl;
		//	cout << "j = " << j << ", k = " << k << " B2Night mean pressure " << B2meanNight << "gcm-2 " << " , EneMean = " << MeanEne << endl;
			
			
			if (j<(NPBins/2)-1) 
				{			
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			//RateDay
			RateDayEle = P0BinDay[j][k]->GetEntries(); 
			RateDayEleErr = TMath::Sqrt(P0BinDay[j][k]->GetEntries());
			FluxDayEle = RateDayEle/(1e-4*GeoFactor[(NPBins/2)-2-j]*1e-3*MomWidth*timeWidthDay);
			FluxDayEleErr = RateDayEleErr/(1e-4*GeoFactor[(NPBins/2)-2-j]*1e-3*MomWidth*timeWidthDay);
          //  cout << "j = " << j << ", k = " << k << " B2meanDay " << B2meanDay << ", MeanEneEle = " << MeanEne <<", FluxDayEle = " << FluxDayEle << "GeoFac = " << GeoFactor[(NPBins/2)-2-j] << endl;

			if(B2meanDay>1.0 && FluxDayEle!=0) {
				FluxRateDay[j]->SetPoint(npDayEle,B2meanDay,FluxDayEle);
				FluxRateDay[j]->SetPointError(npDayEle,0,0,FluxDayEleErr,FluxDayEleErr);
				}
			//Rate Night
			RateNightEle = P0BinNight[j][k]->GetEntries(); 
			RateNightEleErr = TMath::Sqrt(P0BinNight[j][k]->GetEntries());
			FluxNightEle = RateNightEle/(1e-4*GeoFactor[(NPBins/2)-2-j]*1e-3*MomWidth*timeWidthNight);
			FluxNightEleErr = RateNightEleErr/(1e-4*GeoFactor[(NPBins/2)-2-j]*1e-3*MomWidth*timeWidthNight);
         //   cout << "j = " << j << ", k = " << k << " B2meanNight " << B2meanNight << ", MeanEneEle = " << MeanEne <<", FluxNightEle = " << FluxNightEle << "GeoFac = " << GeoFactor[(NPBins/2)-2-j] << endl;
			if(k==5) {
          //   cout << "index MeanEneEle = " << (NPBins/2)-2-j << endl;
			   MeanEneEle[(NPBins/2)-2-j] = MeanEne;
			}

			if(B2meanNight>1.0 && FluxNightEle!=0) {
				FluxRateNight[j]->SetPoint(npNightEle,B2meanNight,FluxNightEle);
				FluxRateNight[j]->SetPointError(npNightEle,0,0,FluxNightEleErr,FluxNightEleErr);
				}	
				npDayEle++;
				npNightEle++;
			} // if j < 5
		
			if (j>(NPBins/2)-2) 
				{			
			MomWidth = TMath::Abs(PBins[j] - PBins[j+1]); //width of momentum bin in MeV
			//RateDay
			RateDayPos = P0BinDay[j][k]->GetEntries(); 
			RateDayPosErr = TMath::Sqrt(P0BinDay[j][k]->GetEntries());
			FluxDayPos = RateDayPos/(1e-4*GeoFactor[j-((NPBins/2)-1)]*1e-3*MomWidth*timeWidthDay);
			FluxDayPosErr = RateDayPosErr/(1e-4*GeoFactor[j-((NPBins/2)-1)]*1e-3*MomWidth*timeWidthDay);
        //    cout << "j = " << j << ", k = " << k << " B2meanDay " << B2meanDay << ", MeanEnePos = " << MeanEne <<", FluxDayPos = " << FluxDayPos << "GeoFac = " << GeoFactor[j-((NPBins/2)-1)] << endl;
			if(B2meanDay>1.0 && FluxDayPos!=0) {
				FluxRateDay[j]->SetPoint(npDayPos,B2meanDay,FluxDayPos);
				FluxRateDay[j]->SetPointError(npDayPos,0,0,FluxDayPosErr,FluxDayPosErr);
				}
			//Rate Night
			RateNightPos = P0BinNight[j][k]->GetEntries(); 
			RateNightPosErr = TMath::Sqrt(P0BinNight[j][k]->GetEntries());
			FluxNightPos = RateNightPos/(1e-4*GeoFactor[j-((NPBins/2)-1)]*1e-3*MomWidth*timeWidthNight);
			FluxNightPosErr = RateNightPosErr/(1e-4*GeoFactor[j-((NPBins/2)-1)]*1e-3*MomWidth*timeWidthNight);
           // cout << "j = " << j << ", k = " << k << " B2meanNight " << B2meanNight << ", MeanEne = " << MeanEne << ", FluxNightPos = " << FluxNightPos << ", GeoFactor = " << GeoFactor[j-((NPBins/2)-1)] << endl;
			if(k==5) {
           //  cout << "index MeanEnePos = " << j-((NPBins/2) -1)<< endl;
			   MeanEnePos[j-((NPBins/2)-1)] = MeanEne;
			}
			if(B2meanNight>1.0 && FluxNightPos!=0) {
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
			double MeanEne = P0BinAll[j][k]->GetMean(); 

			MomWidth = TMath::Abs(PBinsAll[j] - PBinsAll[j+1]); //width of momentum bin in MeV
			//RateDay
			Rate = P0BinAll[j][k]->GetEntries(); 
			RateErr = TMath::Sqrt(P0BinAll[j][k]->GetEntries());
			Flux = Rate/(1e-4*GeoFactor[j]*1e-3*MomWidth*timeWidth);
			FluxErr = RateErr/(1e-4*GeoFactor[j]*1e-3*MomWidth*timeWidth);
			if(B2mean>1.0 && Flux!=0) {
				FluxRateAllEle[j]->SetPoint(npAll,B2mean,Flux);
				FluxRateAllEle[j]->SetPointError(npAll,0,0,FluxErr,FluxErr);
				if(k==5) {
			   // cout << "PressureBin " << k << " point in graph = " << npAll << "  B2all = " << B2mean << " FluxAll = " << Flux << ", MeanEne = " << MeanEne << endl;
				MeanEneAll[j]=MeanEne;
					}
				if(j==0) totxtAll << B2mean << " "; 
				totxtAll << Flux << " " << FluxErr << " ";
				if(j==NPBinsAll-2) totxtAll << endl;	
				npAll++;	
				}
		}
	}	
	/////////////////////////////////////////////////////
	////////////////Save points in txt file/////////////
	///////////////////////////////////////////////////	
	
	int npointsNight = FluxRateNight[1]->GetN();
	int npointsDay = FluxRateDay[1]->GetN();	
	
	for(int i=0;i<npointsDay;i++) 
		{
			for(int j=0; j<(NPBins/2)-1;j++)
			{	
			double B2Ele, B2Pos;
			double FluxEle, FluxPos;
			double FluxErrEle, FluxErrPos;		
			int indexEle = (NPBins/2)-2-j;
		    int indexPos = (NPBins/2)-1+j;
			FluxRateDay[indexEle]->GetPoint(i,B2Ele,FluxEle);
			FluxErrEle = FluxRateDay[indexEle]->GetErrorY(i);
			if(FluxErrEle!=0)
			{
			//cout << "DAY i = " << i << ", j = " << j << ", B2Ele = " << B2Ele << " FluxEle = " << FluxEle << " FluxErrEle = " << FluxErrEle << endl;
			if(j==0)totxtEle << B2Ele << " "; 
			totxtEle << FluxEle << " " << FluxErrEle << " ";
			if(j==(NPBins/2)-2) totxtEle << endl;
			}
			FluxRateDay[indexPos]->GetPoint(i,B2Pos,FluxPos);
			FluxErrPos = FluxRateDay[indexPos]->GetErrorY(i);
			if(FluxErrPos!=0) 
				{
			//cout << "DAY i = " << i << ", j = " << j << ", B2Pos = " << B2Pos << " FluxPos = " << FluxPos << " FluxErrPos = " << FluxErrPos << endl;
			if(j==0)totxtPos << B2Pos << " "; 
			totxtPos << FluxPos << " " << FluxErrPos << " ";
			if(j==(NPBins/2)-2) totxtPos << endl;			
				}
			   } //end j
			} //end i
			for(int i=0;i<npointsNight;i++) 
			{
				for(int j=0; j<(NPBins/2)-1;j++)
				{	
				double B2Ele, B2Pos;
				double FluxEle, FluxPos;
				double FluxErrEle, FluxErrPos;
				int indexEle = (NPBins/2)-2-j;
		        int indexPos = (NPBins/2)-1+j;
				FluxRateNight[indexEle]->GetPoint(i,B2Ele,FluxEle);
				FluxErrEle = FluxRateNight[indexEle]->GetErrorY(i);
				if(FluxErrEle!=0)
					{
			//	cout << "NIGHT i = " << i << ", j = " << j << ", B2Ele = " << B2Ele << " FluxEle = " << FluxEle << " FluxErrEle = " << FluxErrEle << endl;
				if(j==0) totxtEle << B2Ele << " "; 
				totxtEle << FluxEle << " " << FluxErrEle << " ";
				if(j==(NPBins/2)-2) totxtEle << endl;
					}
				FluxRateNight[indexPos]->GetPoint(i,B2Pos,FluxPos);
				FluxErrPos = FluxRateNight[indexPos]->GetErrorY(i);
				if(FluxErrPos!=0) 
					{
				//cout << "NIGHT i = " << i << ", j = " << j << ", B2Pos = " << B2Pos << " FluxPos = " << FluxPos << " FluxErrPos = " << FluxErrPos << endl;
				if(j==0) totxtPos << B2Pos << " "; 
				totxtPos << FluxPos << " " << FluxErrPos << " ";
				if(j==(NPBins/2)-2) totxtPos << endl;			
					}
				}   //end j
	 	  } 		//end i 
 
	/////////////////////////////////////////////////////
	////////////Make Growth Curves TGraphs for MC///////
	///////////////////////////////////////////////////
	TGraphAsymmErrors**FluxEleMC = new TGraphAsymmErrors*[nMCbin];
	TGraphAsymmErrors**FluxPosMC = new TGraphAsymmErrors*[nMCbin];
	TGraphAsymmErrors**FluxAllMC = new TGraphAsymmErrors*[nMCbin];

	for(int i=0;i<nMCbin;i++) 
		{
		FluxEleMC[i] = new TGraphAsymmErrors();
		FluxEleMC[i]->SetMarkerStyle(22);
		FluxEleMC[i]->SetMarkerColor(kBlue);
		FluxPosMC[i] = new TGraphAsymmErrors();
		FluxPosMC[i]->SetMarkerStyle(22);
		FluxPosMC[i]->SetMarkerColor(kRed);	
		FluxAllMC[i] = new TGraphAsymmErrors();
		FluxAllMC[i]->SetMarkerStyle(22);
		FluxAllMC[i]->SetMarkerColor(kBlack);		
		 for(int j=0;j<ndepth;j++)
		 {
			// cout << "Ebin " << i << " Depth bin " <<j << " Depth = " << Depth[j] << ", eMC = " <<eMC[i][j] << ", pMC = " << pMC[i][j] << ", allMC = " << allMC[i][j] << endl;
			 FluxEleMC[i]->SetPoint(j,Depth[j],eMC[i][j]);
			 FluxPosMC[i]->SetPoint(j,Depth[j],pMC[i][j]);
			 FluxAllMC[i]->SetPoint(j,Depth[j],allMC[i][j]);
		 } //end j
	} //end i
	
 
	/////////////////////////////////////////////////////////////////////////
	//////Define Daniel&Stephens-like function, with free parameters/////////
	/////////////////////////////////////////////////////////////////////////


	//Define Daniel & Stephens functions
	TF1**DS = new TF1*[7];
	for (int k=0;k<7;k++) 
	{
		DS[k] = new TF1(Form("ds%d",k),DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DS[k]->SetParameter(l,DScoeff[k][l]);
		DS[k]->SetParameter(7,1);
		DS[k]->SetLineColor(mcolor[k]);
		DS[k]->SetLineStyle(6);

	}
	

	/////////////////////////////////////////////////////////////////////
	////////////Normalize MC point to data using bin at ~50g.cm-2///////
	////////////////////////////////////////////////////////////////////
/*
	for(int i=0;i<NPBins/2 -1;i++)
	{
		double B2Ele, B2Pos;
		double FluxEle, FluxPos;
		double NormEle, NormPos,NormAll;
		double B2EleMC, B2PosMC;
		double EleMC, PosMC;	
		double B2All, FluxAll;
		double B2MC, allMC;
		int indexEle = (NPBins/2)-2-i;
		int indexPos = (NPBins/2)-1+i;
		FluxRateDay[indexEle]->GetPoint(14,B2Ele,FluxEle);								// ~50 g.cm^-2
		FluxRateDay[indexPos]->GetPoint(14,B2Pos,FluxPos);								// ~50 g.cm^-2
		FluxRateAllEle[i]->GetPoint(10,B2All,FluxAll);									// ~50 g.cm^-2
		FluxEleMC[i]->GetPoint(17,B2EleMC,EleMC);
		FluxPosMC[i]->GetPoint(17,B2PosMC,PosMC);
		FluxAllMC[i]->GetPoint(17,B2MC,allMC);
		cout << "Point 17 i = " << i <<  " indexEle = " << indexEle << ", indexPos = " << indexPos << endl;		
		cout << "B2Ele = " << B2Ele << "  B2Pos = " << B2Pos << "  B2All = " << B2All << " B2MC = " << B2MC << endl;
	    cout << "FluxEle = " << FluxEle << "  FluxPos = " << FluxPos << "  FluxAll = " << FluxAll << endl;
		NormEle = FluxEle/EleMC;
		NormPos = FluxPos/PosMC;	
		NormAll = FluxAll/allMC;	
     //   cout << "NormEle = " << NormEle << ", NormPos = " << NormPos << ", NormAll = " << NormAll << endl;
     //   cout << "DiffEle = " << 100*(FluxEle-EleMC)/FluxEle << ",% DiffPos = " << 100*(FluxPos-PosMC)/FluxPos << ",% DifAll = " << 100*(FluxAll-allMC)/FluxAll << "%"<< endl;
		TeXTables << Form("%2.0f-%2.0f MeV & %2.2f & %2.2f &&",PBinsAll[i],PBinsAll[i+1], NormEle, NormPos);
		FluxRateDay[indexEle]->GetPoint(15,B2Ele,FluxEle);								// ~90 g.cm^-2
		FluxRateDay[indexPos]->GetPoint(15,B2Pos,FluxPos);								// ~90 g.cm^-2
		FluxRateAllEle[i]->GetPoint(11,B2All,FluxAll);									// ~90 g.cm^-2
		FluxEleMC[i]->GetPoint(18,B2EleMC,EleMC);
		FluxPosMC[i]->GetPoint(18,B2PosMC,PosMC);
		FluxAllMC[i]->GetPoint(18,B2MC,allMC);
		cout << " Point 18 i = " << i <<  " indexEle = " << indexEle << ", indexPos = " << indexPos << endl;		
		cout << "B2Ele = " << B2Ele << "  B2Pos = " << B2Pos << "  B2All = " << B2All << " B2MC = " << B2MC << endl;
	    cout << "FluxEle = " << FluxEle << "  FluxPos = " << FluxPos << "  FluxAll = " << FluxAll << endl;	
		TeXTables << Form(" %2.2f & %2.2f &&",FluxEle/EleMC, FluxPos/PosMC);
		FluxRateDay[indexEle]->GetPoint(16,B2Ele,FluxEle);								// ~50 g.cm^-2
		FluxRateDay[indexPos]->GetPoint(16,B2Pos,FluxPos);								// ~50 g.cm^-2
		FluxRateAllEle[i]->GetPoint(12,B2All,FluxAll);									// ~50 g.cm^-2
		FluxEleMC[i]->GetPoint(19,B2EleMC,EleMC);
		FluxPosMC[i]->GetPoint(19,B2PosMC,PosMC);
		FluxAllMC[i]->GetPoint(19,B2MC,allMC);
		cout << "Point 19 i = " << i <<  " indexEle = " << indexEle << ", indexPos = " << indexPos << endl;		
		cout << "B2Ele = " << B2Ele << "  B2Pos = " << B2Pos << "  B2All = " << B2All << " B2MC = " << B2MC << endl;
	    cout << "FluxEle = " << FluxEle << "  FluxPos = " << FluxPos << "  FluxAll = " << FluxAll << endl;

		TeXTables << Form(" %2.2f & %2.2f \\\\", FluxEle/EleMC, FluxPos/PosMC) << endl;
		 for (int j=0;j<FluxEleMC[i]->GetN();j++) 
		 {
			 FluxEleMC[i]->GetY()[j] *= NormEle;
			 FluxPosMC[i]->GetY()[j] *= NormPos;
			 FluxAllMC[i]->GetY()[j] *= NormAll;
		 }
	}
     
	TeXTables << " \n \n \n ";
	
	/////////////////////////////////////////////////////////////////////
	////////////Subtract Data - MC for data point at 2.40g.cm-2/////////
	/////////// Fill points for differential flux plots//////////////////
	////////////////////////////////////////////////////////////////////
	TGraphAsymmErrors*DiffFluxEle = new TGraphAsymmErrors();
	TGraphAsymmErrors*DiffFluxPos = new TGraphAsymmErrors();
	TGraphAsymmErrors*DiffFluxAll = new TGraphAsymmErrors();
	DiffFluxEle[i]->SetMarkerColor(kBlue);
	DiffFluxEle[i]->SetLineColor(kBlue);
	DiffFluxEle[i]->SetMarkerStyle(21);
	DiffFluxEle[i]->SetMarkerSize(0.8);
	DiffFluxPos[i]->SetMarkerColor(kRed);
	DiffFluxPos[i]->SetLineColor(kRed);
	DiffFluxPos[i]->SetMarkerStyle(21);
	DiffFluxPos[i]->SetMarkerSize(0.8);
	DiffFluxAll[i]->SetMarkerColor(kBlack);
	DiffFluxAll[i]->SetLineColor(kBlack);
	DiffFluxAll[i]->SetMarkerStyle(21);
	DiffFluxAll[i]->SetMarkerSize(0.8);
	for(int i=0;i<(NPBins/2)-1;i++)
	{
		double B2Ele, B2Pos;
		double FluxEle, FluxPos;
		double FluxEleErr,FluxPosErr,FluxAllErr;
		double NormEle, NormPos,NormAll;
		double B2EleMC, B2PosMC;
		double EleMC, PosMC;	
		double B2All, FluxAll;
		double B2MC, allMC;
		double MeanEle, MeanPos, MeanAll;
		int indexEle = (NPBins/2)-2-i;
		int indexPos = (NPBins/2)-1+i;
		FluxRateNight[indexEle]->GetPoint(6,B2Ele,FluxEle);	
		FluxRateNight[indexPos]->GetPoint(6,B2Pos,FluxPos);	
		FluxRateAllEle[i]->GetPoint(2,B2All,FluxAll);
		FluxEleErr = FluxRateNight[indexEle]->GetErrorY(6);
		FluxPosErr = FluxRateNight[indexPos]->GetErrorY(6);
		FluxAllErr = FluxRateAllEle[i]->GetErrorY(2);
		FluxEleMC[i]->GetPoint(6,B2EleMC,EleMC);
		FluxPosMC[i]->GetPoint(6,B2PosMC,PosMC);
		FluxAllMC[i]->GetPoint(6,B2MC,allMC);
		MeanEle = TMath::Abs(MeanEneEle[i]);
		MeanPos = MeanEnePos[i];
		MeanAll = MeanEneAll[i];
		cout << " i = " << i <<  " indexEle = " << indexEle << ", indexPos = " << indexPos << endl;			   
		cout << "B2Ele = " << B2Ele << " FluxEle = " << FluxEle << ", MeanEneEle = " << MeanEle << endl;
		cout << "B2Pos = " << B2Pos << " FluxPos = " << FluxPos << ", MeanEnePos = " << MeanPos << endl;
	//    cout << "B2All = " << B2All << " FluxAll = " << FluxAll << ", MeanEneAll = " << MeanAll << endl;
		//cout << "B2EleMC = " << B2EleMC << " FluxEleMC = " << EleMC << endl;
		//cout << "B2PosMC = " << B2PosMC << " FluxPosMC = " << PosMC << endl;
		//cout << "B2MC = " << B2MC << " FluxAllMC = " << allMC << endl;
		TeXTables << Form("%3.2f & %4.2f $\\pm$ %4.2f & %4.2f $\\pm$ %4.2f && %4.2f & %4.2f \\\\",MeanAll,FluxEle,FluxEleErr,FluxPos,FluxPosErr,EleMC,PosMC) << endl;
		DiffFluxEle[i]->SetPoint(i,MeanEle,(FluxEle-EleMC));
		DiffFluxPos[i]->SetPoint(i,MeanPos,(FluxPos-PosMC));
		DiffFluxAll[i]->SetPoint(i,MeanAll,(FluxAll-allMC));
		DiffFluxEle[i]->SetPointError(i,0.,0.,FluxEleErr,FluxEleErr);
		DiffFluxPos[i]->SetPointError(i,0.,0.,FluxPosErr,FluxPosErr);
		DiffFluxAll[i]->SetPointError(i,0.,0.,FluxAllErr,FluxAllErr);
		//cout << "Setting Ele point " << i << " MeanEle = " << MeanEle << ", Flux = " << (FluxEle-EleMC) << " , FluxError " << DiffFluxEle[i]->GetErrorY(i) << endl;
		//cout << "Setting Pos point " << i << " MeanPos = " << MeanPos << ", Flux = " << (FluxPos-PosMC) << " , FluxError " << DiffFluxPos[i]->GetErrorY(i)  << endl;
		//cout << "Setting All point " << i << " MeanAll = " << MeanAll << ", Flux = " << (FluxAll-allMC) << " , FluxError " << DiffFluxAll[i]->GetErrorY(i)  << endl;
		

	}	
*/
	
	////////////////////////////////////////////////////////////////////
	/////////////////////Draw Growth Curves////////////////////////////
	///////////////////////////////////////////////////////////////////
	
	int font = 12;
	gStyle->SetTitleFont(font,"xyz");
	gStyle->SetLabelFont(font,"xyz");
	gStyle->SetLegendFont(font);
	gStyle->SetPadTickX(1); 
	gStyle->SetPadTickY(1);	
//	gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);
			//Draw all DS functions
	TCanvas *cDS = new TCanvas("DS","DS",1000,1000);
	TMultiGraph*Functions = new TMultiGraph();
	for (int i=0;i<7;i++)
	{
		cDS->cd();
		if(i==0) DS[i]->Draw();
		else {
			DS[i]->Draw("same");
		}
		Functions->GetXaxis()->SetLimits(PMin,PMax);
		gPad->SetLogy();
		gPad->SetLogx();
		Functions->SetMinimum(ZMin);
		Functions->SetMaximum(ZMax);
		gPad->SetLeftMargin(0.15);
		gPad->Update(); 
	}
		 cDS->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));


	
		//////////////////////////////////////////////////////
	//////Differential Flux at Top of the Payload/////////
	//////////////////////////////////////////////////////
	TCanvas **cdiff = new TCanvas*[7];
	TMultiGraph**DiffFluxToP = new TMultiGraph*[7]; 
	TLegend **legdiff = new TLegend*[7];
	TGraphAsymmErrors**DiffFluxEle = new TGraphAsymmErrors*[7];
	TGraphAsymmErrors**DiffFluxPos = new TGraphAsymmErrors*[7];
	TGraphAsymmErrors**DiffFluxAll = new TGraphAsymmErrors*[7];
	for(int i=0; i<7;i++) {
		cdiff[i] = new TCanvas("cdiff","cdiff",1000,1000);
		DiffFluxToP[i]= new TMultiGraph();
		DiffFluxEle[i] = new TGraphAsymmErrors();
		DiffFluxPos[i] = new TGraphAsymmErrors();	
		DiffFluxAll[i] = new TGraphAsymmErrors();	
		legdiff[i] = new TLegend(0.68,0.65,0.85,0.85);
		legdiff[i]->SetBorderSize(0);
		legdiff[i]->SetFillStyle(0);
		////////////////////////////////////
		////////////make flux plots //////////
		////////////////////////////////////

		DiffFluxEle[i]->SetMarkerColor(kBlue);
		DiffFluxEle[i]->SetLineColor(kBlue);
		DiffFluxEle[i]->SetMarkerStyle(21);
		DiffFluxEle[i]->SetMarkerSize(0.8);
		DiffFluxPos[i]->SetMarkerColor(kRed);
		DiffFluxPos[i]->SetLineColor(kRed);
		DiffFluxPos[i]->SetMarkerStyle(21);
		DiffFluxPos[i]->SetMarkerSize(0.8);
		DiffFluxAll[i]->SetMarkerColor(kBlack);
		DiffFluxAll[i]->SetLineColor(kBlack);
		DiffFluxAll[i]->SetMarkerStyle(21);
		DiffFluxAll[i]->SetMarkerSize(0.8);
					}
	TCanvas **cGC = new TCanvas*[(NPBins/2)-1];
	TMultiGraph**GrowthCurves = new TMultiGraph*[(NPBins/2)-1]; 
	TLegend **leg = new TLegend*[(NPBins/2)-1];
	fileout->cd();	
	for (int i=1;i<(NPBins/2)-1;i++)
	{
		cGC[i] = new TCanvas(Form("c%d",i),Form("c%d",i),1000,1000);
		GrowthCurves[i] = new TMultiGraph();
		if (i<8) {leg[i] = new TLegend(0.18,0.15,0.65,0.35);}
		else {leg[i] = new TLegend(0.18,0.65,0.65,0.85);}
		leg[i]->SetBorderSize(0);
		leg[i]->SetFillStyle(0);
		int indexEle = (NPBins/2)-2-i;
		int indexPos = (NPBins/2)-1+i;
		cout << Form("%2.0f-%2.0f MeV",PBinsAll[i],PBinsAll[i+1]) << endl;
		GrowthCurves[i]->Add(FluxRateNight[indexEle]);
		GrowthCurves[i]->Add(FluxRateNight[indexPos]);
		GrowthCurves[i]->Add(FluxRateDay[indexEle]);
		GrowthCurves[i]->Add(FluxRateDay[indexPos]);
		GrowthCurves[i]->Add(FluxRateAllEle[i]);
		GrowthCurves[i]->Add(FluxEleMC[i]);
		GrowthCurves[i]->Add(FluxPosMC[i]);
		GrowthCurves[i]->Add(FluxAllMC[i]);
		leg[i]->AddEntry(FluxRateAllEle[i],"e^{-} + e^{+} data","p");	
		leg[i]->AddEntry(FluxRateDay[indexEle],"e^{-}","p");
		leg[i]->AddEntry(FluxRateDay[indexPos],"e^{+} data","p");
		leg[i]->AddEntry(FluxEleMC[i],"e^{-} MC","p");
		leg[i]->AddEntry(FluxPosMC[i],"e^{+} MC","p");
		leg[i]->AddEntry(FluxAllMC[i],"e^{-} + e^{+}","p");
		cGC[i]->cd();
		GrowthCurves[i]->Draw("AP");
		
		//////////////////////////////
		/////for all electrons////////
		//////////////////////////////		
	
		//fit function to MC points
		cout << "all electrons" << endl;
		FluxAllMC[i]->Fit("ds0","Q0","",0.5,1000);
		FluxAllMC[i]->GetFunction("ds0")->SetLineColor(kBlack);
		//FluxAllMC[i]->GetFunction("ds0")->DrawF1(0.5,1000,"same");	
		//get fit function 
		TF1* fitfunc = FluxAllMC[i]->GetFunction("ds0");
		//create new function with fixed parameters
		TF1* DSdata = new TF1("dsdata",DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DSdata->FixParameter(l,fitfunc->GetParameter(l));
		DSdata->SetLineStyle(9);
		DSdata->SetLineColor(kBlack);
		//fit data from 50 to 300 g.cm-2
		FluxRateAllEle[i]->Fit("dsdata","Q0","",50,500);
	    DSdata->DrawF1(0.5,1000,"same");
	//	leg[i]->AddEntry(DS[0], "DS function fit to MC");
		leg[i]->AddEntry(DSdata, "DS fit, normalized to data");
		//calculate normalization factor and propagate error for A/B
		float normdata =  DSdata->GetParameter(7);
		float dataerr =  DSdata->GetParError(7);
		float normmc = fitfunc->GetParameter(7);
		float mcerr =  fitfunc->GetParError(7);
		float NormAll = normdata/normmc;		
		float errNormAll = TMath::Abs(NormAll)*TMath::Sqrt((dataerr/normdata)*(dataerr/normdata) + (mcerr/normmc)*(mcerr/normmc));
		//fill TeX table
		TeXTables << Form("%2.0f-%2.0f MeV & %2.3f $\\pm$ %2.3f &",PBinsAll[i],PBinsAll[i+1],NormAll,errNormAll);

		//////////////////////////////
		/////for electrons only///////
		//////////////////////////////
		cout << "electrons" << endl;
		FluxEleMC[i]->Fit("ds0","Q0","",0.5,1000);
		FluxEleMC[i]->GetFunction("ds0")->SetLineColor(kBlue);
		//FluxEleMC[i]->GetFunction("ds0")->DrawF1(0.5,1000,"same");
		
		//get fit function 
		fitele = FluxEleMC[i]->GetFunction("ds0");	
		//create new function with fixed parameters
		TF1* DSele = new TF1("dsele",DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DSele->FixParameter(l,fitele->GetParameter(l));
		DSele->SetLineStyle(9);
		DSele->SetLineColor(kBlue);
		//fit data from 50 to 300 g.cm-2
		FluxRateDay[indexEle]->Fit("dsele","Q","",50,500);
	    DSele->DrawF1(0.5,1000,"same");
		//calculate normalization factor and propagate error for A/B
		normdata =  DSele->GetParameter(7);
		dataerr =  DSele->GetParError(7);
		normmc = fitele->GetParameter(7);
		mcerr =  fitele->GetParError(7);
		NormAll = normdata/normmc;		
		errNormAll = TMath::Abs(NormAll)*TMath::Sqrt((dataerr/normdata)*(dataerr/normdata) + (mcerr/normmc)*(mcerr/normmc));
		TeXTables << Form(" %2.3f $\\pm$ %2.3f &", NormAll,errNormAll);
		
		//////////////////////////////
		/////for positrons only///////
		//////////////////////////////
		cout << "positrons" << endl;
		FluxPosMC[i]->Fit("ds0","0","",0.5,1000);
		FluxPosMC[i]->GetFunction("ds0")->SetLineColor(kRed);
	//	FluxPosMC[i]->GetFunction("ds0")->DrawF1(0.5,1000,"same");
		//get fit function 
		fitpos = FluxPosMC[i]->GetFunction("ds0");
		//create new function with fixed parameters
		TF1* DSpos = new TF1("dspos",DanAndSteph,1,1000,8);
		for(int l=0;l<7;l++) DSpos->FixParameter(l,fitpos->GetParameter(l));
		DSpos->SetLineStyle(9);
		DSpos->SetLineColor(kRed);
		//fit data from 50 to 300 g.cm-2
		FluxRateDay[indexPos]->Fit("dspos","Q0","",50,500);
	    DSpos->DrawF1(0.5,1000,"same");
		//calculate normalization factor and propagate error for A/B
		normdata =  DSpos->GetParameter(7);
		dataerr =  DSpos->GetParError(7);
		normmc = fitpos->GetParameter(7);
		mcerr =  fitpos->GetParError(7);
		NormAll = normdata/normmc;		
		errNormAll = TMath::Abs(NormAll)*TMath::Sqrt((dataerr/normdata)*(dataerr/normdata) + (mcerr/normmc)*(mcerr/normmc));
        TeXTables << Form(" %2.3f $\\pm$ %2.3f \\\\", NormAll,errNormAll) << endl;


		
		/////////////////////////////
		/////sum of two fits////////
		////////////////////////////
/*

		TF1* Datasum = new TF1("sumdata",sumData,0.5,1000,16);
		Datasum->SetLineStyle(9);
		Datasum->SetLineColor(kBlack);	
		Datasum->DrawF1(0.5,1000,"same");
	*/	
		
		////////////////////////////
		///////subtract points//////
		////////////////////////////
		
		int npoints = FluxRateNight[indexEle]->GetN();
		cout << " n points = " << npoints << endl;
		for(int l=1;l<npoints;l++) {
			double FluxEle, FluxPos, FluxAll;
			double FluxEleErr,FluxPosErr,FluxAllErr;
			double B2Ele, B2Pos, B2All;
			double MeanEle, MeanPos, MeanAll;
			FluxRateNight[indexEle]->GetPoint(l,B2Ele,FluxEle);	
			FluxEleErr = FluxRateNight[indexEle]->GetErrorY(l);
			FluxRateNight[indexPos]->GetPoint(l,B2Pos,FluxPos);	
			FluxPosErr = FluxRateNight[indexPos]->GetErrorY(l);
			FluxRateNight[indexEle]->GetPoint(l,B2Ele,FluxEle);	
			FluxEleErr = FluxRateNight[indexEle]->GetErrorY(l);
			FluxRateAllEle[i]->GetPoint(l,B2All,FluxAll);
			FluxAllErr = FluxRateAllEle[i]->GetErrorY(l);
			MeanEle = TMath::Abs(MeanEneEle[i]);
			MeanPos = MeanEnePos[i];
			MeanAll = MeanEneAll[i];
			double MCfluxAll = DSdata->Eval(B2All);

			if(B2All >1 && B2All < 4) {
			cout << "l = " << l << " MeanAll = " << MeanAll << ", B2All = " << B2All << ", Data-MC = " << (FluxAll-MCfluxAll) << ", yError = " << FluxAllErr << endl;
			DiffFluxAll[l-1]->SetPoint(i,MeanAll,(FluxAll-MCfluxAll));
			DiffFluxAll[l-1]->SetPointError(i,0.,0.,FluxAllErr,FluxAllErr);
			}

			if(B2Ele>1) {
				double MCfluxEle = DSele->Eval(B2Ele);
				double MCfluxPos = DSpos->Eval(B2Pos);
				DiffFluxEle[l-5]->SetPoint(i,MeanEle,(FluxEle-MCfluxEle));
			    DiffFluxPos[l-5]->SetPoint(i,MeanPos,(FluxPos-MCfluxPos));
				DiffFluxEle[l-5]->SetPointError(i,0.,0.,FluxEleErr,FluxEleErr);
		    	DiffFluxPos[l-5]->SetPointError(i,0.,0.,FluxPosErr,FluxPosErr);	
				cout << "l = " << l << " MeanEle = " << MeanEle << ", B2Ele = " << B2Ele << ", Data-MC = " << (FluxEle-MCfluxEle) << ", FluxEleErr = " << FluxEleErr << endl;
				TeXTableFluxes[l-5] << Form("%3.2f & %4.2f $\\pm$ %4.2f & %4.2f $\\pm$ %4.2f && %4.2f & %4.2f \\\\",MeanAll,FluxEle,FluxEleErr,FluxPos,FluxPosErr,MCfluxEle,MCfluxPos) << endl;
	


		   }
		}
			
		leg[i]->Draw("same");
		GrowthCurves[i]->GetXaxis()->SetTitleOffset(0.87);
		GrowthCurves[i]->GetYaxis()->SetTitleOffset(0.87);
		GrowthCurves[i]->GetXaxis()->SetTitleSize(0.05);
		GrowthCurves[i]->GetYaxis()->SetTitleSize(0.05);
		GrowthCurves[i]->GetXaxis()->CenterTitle();	
		GrowthCurves[i]->GetXaxis()->SetTitle("Atmospheric depth (in g.cm^{-2})");
		GrowthCurves[i]->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
		GrowthCurves[i]->GetYaxis()->CenterTitle();
		GrowthCurves[i]->GetHistogram()->SetTitleFont(font,"t");
		GrowthCurves[i]->GetHistogram()->SetTitle(Form("%2.0f-%2.0f MeV",PBinsAll[i],PBinsAll[i+1]));
		GrowthCurves[i]->GetXaxis()->SetLimits(PMin,PMax);
		gPad->SetLogy();
		gPad->SetLogx();
		GrowthCurves[i]->SetMinimum(ZMin);
		GrowthCurves[i]->SetMaximum(ZMax);
		gPad->SetLeftMargin(0.15);
		gPad->Update(); 
	
		//if(i==1) {cGC[i]->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf(",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));}
		 //if(i==(NPBins/2)-2) {cGC[i]->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));}
		 cGC[i]->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));

		//Write histograms in ROOT file	
		cGC[i]->Write();

 } //end loop on energy bins
	


	for(int i=0; i<7;i++) {
		
					
		DiffFluxToP[i]->Add(DiffFluxEle[i]);
		DiffFluxToP[i]->Add(DiffFluxPos[i]);
		DiffFluxToP[i]->Add(DiffFluxAll[i]);
		legdiff[i]->AddEntry(DiffFluxAll[i],"e^{-} + e^{+} ","p");	
		legdiff[i]->AddEntry(DiffFluxEle[i],"e^{-}","p");
		legdiff[i]->AddEntry(DiffFluxPos[i],"e^{+} ","p");
		cdiff[i]->cd();
		DiffFluxToP[i]->Draw("APW");
		legdiff[i]->Draw("same");
		DiffFluxToP[i]->GetXaxis()->SetTitleOffset(0.87);
		DiffFluxToP[i]->GetYaxis()->SetTitleOffset(0.87);
		DiffFluxToP[i]->GetXaxis()->SetTitleSize(0.05);
		DiffFluxToP[i]->GetYaxis()->SetTitleSize(0.05);
		DiffFluxToP[i]->GetXaxis()->CenterTitle();	
		DiffFluxToP[i]->GetXaxis()->SetTitle("Energy (in MeV)");
		DiffFluxToP[i]->GetYaxis()->SetTitle("Flux (m^{2} sr s GeV)^{-1}");
		DiffFluxToP[i]->GetYaxis()->CenterTitle();
		DiffFluxToP[i]->GetHistogram()->SetTitleFont(font,"t");
		DiffFluxToP[i]->GetHistogram()->SetTitle(Form("Flux at the %2.2f g.cm^{-2}",FloatDepth[i]));
		DiffFluxToP[i]->GetXaxis()->SetLimits(0,1000);
		gPad->SetLogy();
		gPad->SetLogx();
		gPad->SetGrid(1,1); 
		DiffFluxToP[i]->SetMinimum(ZMin);
		DiffFluxToP[i]->SetMaximum(10000);
		gPad->SetLeftMargin(0.15);
		gPad->Update(); 
		if(i==6) cdiff[i]->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf)",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
		else cdiff[i]->Print(Form("%s/GrowthCurvesBinnedSigned_%s_%s.pdf",Outpath.c_str(),ID.c_str(),CutConfig.c_str()));
		
	}


	cout << "now closing all files" << endl;
	fileout->Close();
	fileCal->Close(); 	
		
} //end function		