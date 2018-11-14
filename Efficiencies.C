#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void Efficiencies(string s, int geoconf)
{
	
  //Load configuration parameter
float* zL=new float[7];
float*OffsetLL=new float[7];
float*OffsetRL=new float[7];
float*TrigThresh=new float[5];
for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
for(int i=0;i<5;i++)TrigThresh[i]=0;
string paramfile=Form("../ANALYSIS_SOFT/src/ALSim/Dataparameters%d.dat",geoconf); 
LoadDataparameters(paramfile,zL,OffsetLL,OffsetRL,TrigThresh);	

string file=s;
	
//Number of reconstruction types to study 
const int NReco=1;
string RecoInd[NReco]={"KFone"};

//Number of configurations to study
const int NConfig=6;
string TConfig[NConfig]={"T1&T4", "T1&T4&reco", "T1&T3&T4","T1&T4&T3&reco", "T1&T2&T4&reco", "T1&noT2&T4&reco"};					   
const int nEfficiencies=13; 	 		//count and record number of entries for all events, T1&T4 no hits, T1&T4 1 layer with hit, T1&T4 2 layer with hits, etc
string sEff[nEfficiencies] = {"T1&T4","T1&T4&hits","T1&T4 & 1+ layers","T1&T4 & 2+ layers", " T1&T4 & 3+ layers", "T1&T4 & 4+ layers", "T1&T4 & 5+ layers ", "T1T4 & 5+ layers & track","T1&T4 & 6+ layers","T1&T4 & 6+ layers & track", "T1&T4 & 7+ layers", "T1&T4 7 layers & track","T1&T4&reco"};	
//Input file 
//Input file 
string Inppath="/home/sarah/AESOPLITE/MCProduction/Detector/ProcessedFiles";	
string Outpath="/home/sarah/AESOPLITE/Analysis/MCAnalysis/Efficiencies/V4/20cmSource";
string startfile="aesopliteNonUniB_V4"; 
string endfile="_fort.99";;


 int mcolor[4]={432-3,400-9,616-9,632}; 
//zz0
float zz0=0;	
//Number of entries per reconstruction type per energy 

int* nevents= new int[NReco];
	
//Rate Histograms
TH1F**nEff= new TH1F*[NReco];	
//Reco efficiency histograms
TH1F**nEffReco= new TH1F*[NReco];


 for(int i=0;i<NReco;i++) 
   {
	
			nevents[i]=0;
		        //Rate histogram
                        nEff[i] = new TH1F(Form("%s, Rate efficiencies",RecoInd[i].c_str()),Form("%s, Rate efficiencies",RecoInd[i].c_str()),NConfig,0,NConfig);
                        nEff[i]->SetStats(0);
                        nEff[i]->GetXaxis()->SetAlphanumeric();
                        for (int l=1;l<=NConfig;l++) nEff[i]->GetXaxis()->SetBinLabel(l,TConfig[l-1].c_str());
		         //Reconstruction Rate histogram
                        nEffReco[i] = new TH1F(Form("%s, Reconstruction Rate efficiencies",RecoInd[i].c_str()),Form("%s Reconstruction Rate efficiencies",RecoInd[i].c_str()),nEfficiencies,0,nEfficiencies);
                        nEffReco[i]->SetStats(0);
                        nEffReco[i]->GetXaxis()->SetAlphanumeric();
                        for (int l=1;l<=nEfficiencies;l++) nEffReco[i]->GetXaxis()->SetBinLabel(l,sEff[l-1].c_str());
 	} // end i
cout <<"done creating histograms"<< endl;
//output ROOT file  
TFile **fileout = new TFile*[NReco];
TChain **chain=new TChain*[NReco]; 

for(int i=0; i<NReco; i++) 
	{

	fileout[i]=new TFile(Form("/Efficiencies_%s_%s.root", file.c_str(),RecoInd[i].c_str()),"RECREATE");
	//Create TChain to merge the files for a given energy
	//chain[i] =new TChain("Data"); 
    //chain[i]->Add(Form("%s/NL2084.BPD.EVENT_%s.root",Inppath.c_str(),RecoInd[i].c_str()));


   //Define variables to read event
   ALEvent *e = new ALEvent();      
   //Set address to access event data
   chain[i]->SetBranchAddress("event",&e); 
 
   // Get number of event in Tree
   int nentries=chain[i]->GetEntries();

   cout << "Number  of events: " << nentries << endl;  
   for (int k=0;k<nentries;k++)
   {
	chain[i]->GetEntry(k); //Load the entry i in the variable e 
	if(k%10000==0)   cout << "Reco Type: " << RecoInd[i] << endl;
	uint8_t Ti=(uint8_t)e->get_Ti();
        int  nhits = e->get_Nhits();
	//Number of layers with hit(s)
	int NL= e->get_NLayers();
	//check trigger requirements
	bool T1=e->get_T1();
	bool T2=e->get_T2();
	bool T3=e->get_T3();
	bool T4=e->get_T4();
	//Number of layers with hit(s) in bending/non-bending plane
	int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	float pPR=1000*fabs(e->get_p0PR());   //in MeV
	float ePR=e->get_EkPR();
	double deflecPR = e->get_deflecPR();	  
	double chiNBPR = e->get_chi2NBPR();
	double chiBPR = e->get_chi2BPR();
	double Ekreco = e->get_Ekreco();			//in MeV
	double pReco = e->get_p0reco();
    double ndf = e->get_ndf();
	double chireco = e->get_chi2();
	double phi0 = e->get_phi0();
	double tanl = e->get_tanl();
	double cpa_reco = e->get_cpa();
	double sigma_cpa = TMath::Sqrt(e->get_cpaerr2());	//in GeV
    int sitesUsed=0;

//Fill counter histograms
         if(T1 && T4){
		nEff[i]->Fill(TConfig[0].c_str(),1);
         	nevents[i]++;
		   }
    if(T1 && T4 && pReco>0){nEff[i]->Fill(TConfig[1].c_str(),1);}    
	if(T1 && T4 && T3 ){nEff[i]->Fill(TConfig[2].c_str(),1);} 
	if(T1 && T4 && T3 && pReco>0){nEff[i]->Fill(TConfig[3].c_str(),1);} 
	if(T1 && T4 && T2 && pReco>0){nEff[i]->Fill(TConfig[4].c_str(),1);} 
	if(T1 && T4 && (!T2) && pReco>0){nEff[i]->Fill(TConfig[5].c_str(),1);} 

//Fill histograms

	if(T1 && T4 && pReco>0) {nEffReco[i]->Fill(sEff[0].c_str(),1);}
    if(T1 && T3 && T4 && pReco>0) {nEffReco[i]->Fill(sEff[1].c_str(),1);}
	if(T1 && T4 && T3 && pReco>0 && chireco>=0) {nEffReco[i]->Fill(sEff[2].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=0){nEffReco[i]->Fill(sEff[3].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=1){nEffReco[i]->Fill(sEff[4].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=2){nEffReco[i]->Fill(sEff[5].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=3){nEffReco[i]->Fill(sEff[6].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=4){nEffReco[i]->Fill(sEff[7].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=5){nEffReco[i]->Fill(sEff[8].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed>=6){nEffReco[i]->Fill(sEff[9].c_str(),1);}
    if(T1 && T4 && T3 && pReco>0 && sitesUsed==7){nEffReco[i]->Fill(sEff[10].c_str(),1);}
	    } // end k, number of entries 
    } // end i, reco type

//////////////////////////////
// Display histograms////////
/////////////////////////////   

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);	
gROOT->SetStyle("Modern");
//////////////////////////////
////Efficiencies/////////////
/////////////////////////////
cout << "start display canvas" << endl;
TCanvas**CanEff=new TCanvas*[NReco]; 
TLegend**LegEff=new TLegend*[NReco];
TText***txt=new TText**[NReco];
double barwidth=0.2;
double baroffset=0.1;
 for(int i=0;i<NReco;i++)
   {
       CanEff[i] = new TCanvas(Form("%s, efficiencies", RecoInd[i].c_str()),Form("%s, efficiencies", RecoInd[i].c_str()),200,10,1200,800);
	   LegEff[i] = new TLegend(0.1,0.9,0.2,1.0);
       CanEff[i]->cd();
       CanEff[i]->SetBottomMargin(0.32);
   	   gPad->SetGridy(1);
   	   nEff[i]->SetLineWidth(0);
   	   nEff[i]->SetMarkerStyle(20);
   	   nEff[i]->SetMarkerSize(0);
	   nEff[i]->SetFillColor(kGreen);
       nEff[i]->Scale(100./nevents[i]);
       nEff[i]->SetBarWidth(barwidth);
       nEff[i]->LabelsOption("h","X"); 
       nEff[i]->GetYaxis()->SetRangeUser(0,110);
       nEff[i]->GetYaxis()->SetNdivisions(511,0);
   	   nEff[i]->GetYaxis()->SetTitle("Normalized rates");
	   nEff[i]->Draw("bar");
	   nEff[i]->SetTitle(Form("%s Trigger Rates",RecoInd[i].c_str()));

	   txt[i]=new TText*[NConfig];
	  for(int l=0;l<NConfig;l++)
     		 {
       txt[i][l]=new TText(nEff[i]->GetXaxis()->GetBinCenter(l+1),nEff[i]->GetBinContent(l+1)+1,Form("%d %%",(int)nEff[i]->GetBinContent(l+1))); 
       txt[i][l]->SetTextSize(0.025);
       txt[i][l]->SetTextAngle(55);    
       txt[i][l]->Draw();    
      }  // end l

     fileout[i]->cd();
    // LegEff[i]->Draw();
     gPad->Update();
     CanEff[i]->Write();
     CanEff[i]->SaveAs(Form("GroundRuns/Efficiencies.eps"));
    if (i==0) CanEff[i]->Print(Form("Efficiencies.pdf"));
	} //end i

//////////////////////////
//Reconstruction Rates////
/////////////////////////

TCanvas**CanEffReco=new TCanvas*[NReco];
TLegend**LegEffReco=new TLegend*[NReco];
TText****txtReco=new TText***[NReco];
 for(int i=0;i<NReco;i++)
   {
           CanEffReco[i]= new TCanvas(Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),Form("%s,reconstruction efficiencies", RecoInd[i].c_str()),200,10,1200,800);
    	   txtReco[i]=new TText**[Nene];
           LegEffReco[i] = new TLegend(0.1,0.9,0.2,1.0);
           CanEffReco[i]->cd();
           CanEffReco[i]->SetBottomMargin(0.32);

    for (int j=0; j<Nene;j++)
        {
           gPad->SetGridy(1);
           nEffReco[i]->SetLineWidth(0);
           nEffReco[i]->SetMarkerStyle(20);
           nEffReco[i]->SetFillColor(mcolor[j]);
           nEffReco[i]->SetMarkerSize(0);
           nEffReco[i]->Scale(100./p0reco[i][0]->GetEntries());
           nEffReco[i]->SetBarWidth(barwidth);
           nEffReco[i]->SetBarOffset(baroffset+j*barwidth);
           nEffReco[i]->LabelsOption("v","X");
           nEffReco[i]->GetYaxis()->SetRangeUser(0,110);
           nEffReco[i]->GetYaxis()->SetNdivisions(511,0);
           nEffReco[i]->GetYaxis()->SetTitle("Normalized rates");
 if (j==0) {
                nEffReco[i]->Draw("bar");
                nEffReco[i]->SetTitle(Form("%s reconstruction efficiencies",RecoInd[i].c_str()));
                }
         else {
                nEffReco[i]->Draw("bar same");
                }
           LegEffReco[i]->AddEntry(nEffReco[i],Form("%d%s",Ene[j],UNIT.c_str()),"f");
           txtReco[i]=new TText*[NConfig];
          for(int l=0;l<NConfig;l++)
                 {
       if(nEffReco[i]->GetBinContent(l+1)>=100&&l==0) txtReco[i][l]=new TText(nEffReco[i]->GetXaxis()->GetBinCenter(l+1),nEffReco[i]->GetBinContent(l+1)+1,Form("%d",(int)nEffReco[i]->GetBinContent(l+1)));
       else txtReco[i][l]=new TText((nEffReco[i]->GetXaxis()->GetBinCenter(l+1)+(j-1)*barwidth)-0.15,nEffReco[i]->GetBinContent(l+1)+1,Form("%.1f",nEffReco[i]->GetBinContent(l+1)));
       txtReco[i][l]->SetTextSize(0.015);
       txtReco[i][l]->SetTextAngle(55);
       txtReco[i][l]->Draw();
      }  // end l
   } //end j
     fileout[i]->cd();
     LegEffReco[i]->Draw();
     gPad->Update();
     CanEffReco[i]->Write();
     CanEffReco[i]->SaveAs(Form("%s/%d/KF/%s_EfficienciesReco.eps", directory.c_str(),type,RecoInd[i].c_str()));
      if (i==NReco-1)CanEffReco[i]->Print(Form("%s/%d/KF/Efficiencies.pdf)", directory.c_str(),type));
     else CanEffReco[i]->Print(Form("%s/%d/KF/Efficiencies.pdf", directory.c_str(),type));
        } //end i




} //end function


