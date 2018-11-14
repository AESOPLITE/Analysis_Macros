#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "TChain.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TPaveStats.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void EfficienciesData(string s, int geoconf)
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
string RecoInd[NReco]={"PRonly"};

const int nEfficiencies=7;
//const int nEfficiencies=13; 	 		//count and record number of entries for all events, T1&T4 no hits, T1&T4 1 layer with hit, T1&T4 2 layer with hits, etc
string sEff[nEfficiencies] = {"TopShell","TopShell&T1","TopShell&T1&T2","TopShell&T1&T2&T3",
							 "TopShell&T1&T2&T3&T4","TopShell&T1&T2&T3&T4&PR", 
							  "TopShell&T1&T2&T3&T4&PR&reco"};	
//string sEff[nEfficiencies] = {"T1&T4", "T1&T4&reco", "T1&T4&T3&reco","T1&T4 & 1+ layers","T1&T4 & 2+ layers", " T1&T4 & 3+ layers", "T1&T4 & 4+ layers", "T1&T4 & 5+ layers ", "T1&T4 & 5+ layers & track","T1&T4 & 6+ layers","T1&T4 & 6+ layers & track", "T1&T4 & 7+ layers", "T1&T4 7 layers & track"};	
//const int NConfig=12;
//string TConfig[NConfig]={"T1&T4", "T1&T4&reco","T1&T4&T3&reco","T1&T4&T3&reco chi>0" ,"T1&T2&T3&reco & 0+ sites used","T1&T4&T3&reco & 1+ sites","T1&T4&T3&reco & 2+ sites","T1&T4&T3&reco & 3+ sites","T1&T4&T3&reco & 4+sites","T1&T4&T3&reco & 5+ sites","T1&T4&T3&reco & 6+ sites","T1&T4&T3&reco & 7 sites"};					   
const int NConfig=4;
string TConfig[NConfig]={"T1&T3", "T1&T3&5layers+",
						 "T1&T3&6layers+","T1&T3&7layers"};					   

	
//Input file 
//Input file 
string Inppath="/home/sarah/AESOPLITE/GroundRuns/Sweden2018";
string directory= "/home/sarah/AESOPLITE/Analysis/GroundRuns";

 int mcolor[4]={432-3,400-9,616-9,632}; 
//zz0
float zz0=0;	
//Number of entries per reconstruction type per energy 

int* nT1T4= new int[NReco];
int* nevents= new int[NReco];
int *nT1T3 = new int[NReco];
//Rate Histograms
TH1F**nEff= new TH1F*[NReco];
//Reco efficiency histograms
TH1F**nEffReco= new TH1F*[NReco];

 for(int i=0;i<NReco;i++)
   {
                        nevents[i]=0;
	 					nT1T4[i]=0;
	 					nT1T3[i]=0;
                        //Rate histogram
                        nEff[i] = new TH1F(Form("%s, Rate efficiencies",RecoInd[i].c_str()),Form("%s, Rate efficiencies",RecoInd[i].c_str()),nEfficiencies,0,nEfficiencies);
                        nEff[i]->SetStats(0);
                        nEff[i]->GetXaxis()->SetAlphanumeric();
                        for (int l=1;l<=NConfig;l++) nEff[i]->GetXaxis()->SetBinLabel(l,sEff[l-1].c_str());
                         //Reconstruction Rate histogram
                        nEffReco[i] = new TH1F(Form("%s, Reconstruction Rate efficiencies",RecoInd[i].c_str()),Form("%s Reconstruction Rate efficiencies",RecoInd[i].c_str()),NConfig,0,NConfig);
                        nEffReco[i]->SetStats(0);
                        nEffReco[i]->GetXaxis()->SetAlphanumeric();
                        for (int l=1;l<=nEfficiencies;l++) nEffReco[i]->GetXaxis()->SetBinLabel(l,TConfig[l-1].c_str());
        } // end i

cout <<"done creating histograms"<< endl;
//output ROOT file  
TFile **fileout = new TFile*[NReco];
TChain **chain=new TChain*[NReco]; 

for(int i=0; i<NReco; i++) 
	{
        fileout[i]=new TFile(Form("%s/Efficiencies_%s_%s.root", directory.c_str(),file.c_str(),RecoInd[i].c_str()),"RECREATE");
        //Create TChain to merge the files for a given energy
      chain[i] =new TChain("Data");
    chain[i]->Add(Form("%s/NLAll.BPD.EVENT_PRonly.root",Inppath.c_str()));


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
	double preco = e->get_p0reco();
    double ndf = e->get_ndf();
	double chireco = e->get_chi2();
	double phi0 = e->get_phi0();
	double tanl = e->get_tanl();
	double cpa_reco = e->get_cpa();
	double sigma_cpa = TMath::Sqrt(e->get_cpaerr2());	//in GeV
	int sitesUsed=0;
	nevents[i]++;
	 //Fill counter histograms
/*
         if(TShell){nEff[i]->Fill(sEff[0].c_str(),1);}
         if(TShell && T1){nEff[i]->Fill(sEff[1].c_str(),1);} 
         if(TShell && T1 && T2){nEff[i]->Fill(sEff[2].c_str(),1);} 
		 if(TShell && T1 && T2 && T3){nEff[i]->Fill(sEff[3].c_str(),1);}
         if(TShell && T1 && T2 && T3 && T4){nEff[i]->Fill(sEff[4].c_str(),1);}
         if(TShell && T1 && T2 && T3 && T4 && deflecPR>0){nEff[i]->Fill(sEff[5].c_str(),1);}
         if(TShell && T1 && T2 && T3 && T4 && deflecPR>0 & preco>0){nEff[i]->Fill(sEff[6].c_str(),1);}
*/
	   if(T1 && T3 ){
		   nEffReco[i]->Fill(TConfig[0].c_str(),1);
	   	   nT1T3[i]++;} 	   
         if(T1 && T3 && NL>=5){nEffReco[i]->Fill(TConfig[1].c_str(),1);} 
		 if(T1 && T3 && NL>=6){nEffReco[i]->Fill(TConfig[2].c_str(),1);}
         if(T1 && T3 && NL>=7){nEffReco[i]->Fill(TConfig[3].c_str(),1);}
	  /* 
//Fill counter histograms
         if(T1 && T4){
		nEff[i]->Fill(sEff[0].c_str(),1);
         	nT1T4[i]++;
		   }
         if(T1 && T4 && deflecPR>0){nEff[i]->Fill(sEff[1].c_str(),1);} 
         if(T1 && T4 && T3 && deflecPR>0){nEff[i]->Fill(sEff[2].c_str(),1);} 
		 if(T1 && T4 && nhits>0 && NL>=1){nEff[i]->Fill(sEff[3].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=2){nEff[i]->Fill(sEff[4].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=3){nEff[i]->Fill(sEff[5].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=4){nEff[i]->Fill(sEff[6].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=5){nEff[i]->Fill(sEff[7].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=5 && chiNBPR>=0 && chiBPR>=0){nEff[i]->Fill(sEff[8].c_str(),1);}	
         if(T1 && T4 && nhits>0 && NL>=6){nEff[i]->Fill(sEff[9].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL>=6 && chiNBPR>=0 && chiBPR>=0){nEff[i]->Fill(sEff[10].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL==7){nEff[i]->Fill(sEff[11].c_str(),1);}
         if(T1 && T4 && nhits>0 && NL==7 && chiNBPR>=0 && chiBPR>=0){nEff[i]->Fill(sEff[12].c_str(),1);}

//Fill histograms

	if(T1 && T4 ) 
	  {
	  nEffReco[i]->Fill(TConfig[0].c_str(),1);
	  p0PR[i][0]->Fill(1/pPR);
	  p0reco[i][0]->Fill(1/preco);
	  pMC[i][0]->Fill(1/p0MC);
	  chi2reco[i][0]->Fill(chireco);
	  chi2NBPR[i][0]->Fill(chiNBPR);
	  chi2BPR[i][0]->Fill(chiBPR);		
	  ResoReco[i]->Fill((cpa_reco-cpa_MC)/(sigma_cpa)); 
	 }
	if(T1 && T4 && deflecPR>0 ) 
	  {
	  nEffReco[i]->Fill(TConfig[1].c_str(),1);
	  p0PR[i][1]->Fill(1/pPR);
	  p0reco[i][1]->Fill(1/preco);
	  pMC[i][1]->Fill(1/p0MC);
	  chi2reco[i][1]->Fill(chireco);
	  chi2NBPR[i][1]->Fill(chiNBPR);
	  chi2BPR[i][1]->Fill(chiBPR);		
	  ResoReco[i]->Fill((cpa_reco-cpa_MC)/(sigma_cpa)); 
	 }
        if(T1 && T3 && T4 && deflecPR>0)
          {
          nEffReco[i]->Fill(TConfig[2].c_str(),1);
          p0PR[i][2]->Fill(1/pPR);
          p0reco[i][2]->Fill(1/preco);
          pMC[i][2]->Fill(1/p0MC);
          chi2reco[i][2]->Fill(chireco);
          chi2NBPR[i][2]->Fill(chiNBPR);
          chi2BPR[i][2]->Fill(chiBPR);
          }
	if(T1 && T4 && T3 && preco>0 && deflecPR>=0)
          {
          nEffReco[i]->Fill(TConfig[3].c_str(),1);
          p0PR[i][3]->Fill(1/pPR);
          p0reco[i][3]->Fill(1/preco);
          pMC[i][3]->Fill(1/p0MC);
          chi2reco[i][3]->Fill(chireco);
          chi2NBPR[i][3]->Fill(chiNBPR);
          chi2BPR[i][3]->Fill(chiBPR);
          }
       
	  if(T1 && T4 && T3 && deflecPR>0 && deflecPR>=0 )
          {
          nEffReco[i]->Fill(TConfig[4].c_str(),1);
          p0PR[i][4]->Fill(1/pPR);
          p0reco[i][4]->Fill(1/preco);
          pMC[i][4]->Fill(1/p0MC);
          chi2reco[i][4]->Fill(chireco);
          chi2NBPR[i][4]->Fill(chiNBPR);
          chi2BPR[i][4]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=1 )
          {
          nEffReco[i]->Fill(TConfig[5].c_str(),1);
          p0PR[i][5]->Fill(1/pPR);
          p0reco[i][5]->Fill(1/preco);
          pMC[i][5]->Fill(1/p0MC);
          chi2reco[i][5]->Fill(chireco);
          chi2NBPR[i][5]->Fill(chiNBPR);
          chi2BPR[i][5]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=2 )
          {
          nEffReco[i]->Fill(TConfig[6].c_str(),1);
          p0PR[i][6]->Fill(1/pPR);
          p0reco[i][6]->Fill(1/preco);
          pMC[i][6]->Fill(1/p0MC);
          chi2reco[i][6]->Fill(chireco);
          chi2NBPR[i][6]->Fill(chiNBPR);
          chi2BPR[i][6]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=3 )
          {
          nEffReco[i]->Fill(TConfig[7].c_str(),1);
          p0PR[i][7]->Fill(1/pPR);
          p0reco[i][7]->Fill(1/preco);
          pMC[i][7]->Fill(1/p0MC);
          chi2reco[i][7]->Fill(chireco);
          chi2NBPR[i][7]->Fill(chiNBPR);
          chi2BPR[i][7]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=4 )
          {
          nEffReco[i]->Fill(TConfig[8].c_str(),1);
          p0PR[i][8]->Fill(1/pPR);
          p0reco[i][8]->Fill(1/preco);
          pMC[i][8]->Fill(1/p0MC);
          chi2reco[i][8]->Fill(chireco);
          chi2NBPR[i][8]->Fill(chiNBPR);
          chi2BPR[i][8]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=5 )
          {
          nEffReco[i]->Fill(TConfig[9].c_str(),1);
          p0PR[i][9]->Fill(1/pPR);
          p0reco[i][9]->Fill(1/preco);
          pMC[i][9]->Fill(1/p0MC);
          chi2reco[i][9]->Fill(chireco);
          chi2NBPR[i][9]->Fill(chiNBPR);
          chi2BPR[i][9]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed>=6 )
          {
          nEffReco[i]->Fill(TConfig[10].c_str(),1);
          p0PR[i][10]->Fill(1/pPR);
          p0reco[i][10]->Fill(1/preco);
          pMC[i][10]->Fill(1/p0MC);
          chi2reco[i][10]->Fill(chireco);
          chi2NBPR[i][10]->Fill(chiNBPR);
          chi2BPR[i][10]->Fill(chiBPR);
          }
        if(T1 && T4 && T3 && preco>0 && sitesUsed==7 )
          {
          nEffReco[i]->Fill(TConfig[11].c_str(),1);
          p0PR[i][11]->Fill(1/pPR);
          p0reco[i][11]->Fill(1/preco);
          pMC[i][11]->Fill(1/p0MC);
          chi2reco[i][11]->Fill(chireco);
          chi2NBPR[i][11]->Fill(chiNBPR);
          chi2BPR[i][11]->Fill(chiBPR);
          }
*/
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
TCanvas**CanEffReco=new TCanvas*[NReco]; 
TLegend**LegEffReco=new TLegend*[NReco];
TText***txtReco=new TText**[NReco];
double barwidth=0.2;
double baroffset=0.3;
 for(int i=0;i<NReco;i++)
   {
		   CanEffReco[i] = new TCanvas(Form("%s, efficiencies", RecoInd[i].c_str()),Form("%s, efficiencies", RecoInd[i].c_str()),200,10,1200,800);
		   LegEffReco[i] = new TLegend(0.1,0.9,0.2,1.0);
		   CanEffReco[i]->cd();
		   CanEffReco[i]->SetBottomMargin(0.32);
           gPad->SetGridy(1);
           nEffReco[i]->SetLineWidth(0);
           nEffReco[i]->SetMarkerStyle(20);
           nEffReco[i]->SetMarkerSize(0);
           nEffReco[i]->SetFillColor(mcolor[i]);
		   nEffReco[i]->Scale(100./nT1T3[i]);
		   nEffReco[i]->SetBarWidth(barwidth);
	 	   nEffReco[i]->SetBarOffset(baroffset);
		   nEffReco[i]->LabelsOption("v","X");
		   nEffReco[i]->GetYaxis()->SetRangeUser(0,110);
		   nEffReco[i]->GetYaxis()->SetNdivisions(511,0);
           nEffReco[i]->GetYaxis()->SetTitle("Rates normalized to T1T3");
           nEffReco[i]->Draw("bar");
           nEffReco[i]->SetTitle(Form("%s %s Trigger Rates",file.c_str(),RecoInd[i].c_str()));
           txtReco[i]=new TText*[NConfig];
           for(int l=0;l<NConfig;l++)
                 {
		//	   if(nEff[i]->GetBinContent(l+1)>=100&&l==0) txt[i][l]=new TText((nEff[i]->GetXaxis()->GetBinCenter(l+1)-0.15),nEff[i]->GetBinContent(l+1)+1,Form("(%.2f)",100*(float(nevents[i])/float(nT1T4[i])))); 
			   txtReco[i][l]=new TText((nEffReco[i]->GetXaxis()->GetBinCenter(l+1))-0.15,nEffReco[i]->GetBinContent(l+1)+1,Form("%.1f",nEffReco[i]->GetBinContent(l+1))); 			 
			   txtReco[i][l]->SetTextSize(0.015);
			   txtReco[i][l]->SetTextAngle(55);
			   txtReco[i][l]->Draw();
      }  // end l
	
     fileout[i]->cd();
     gPad->Update();
     CanEffReco[i]->Write();
  //  if (i==0) CanEff[i]->Print(Form("%s/Efficiencies_%s.pdf(", directory.c_str(),file.c_str()));
   // if (i==NReco-1)CanEff[i]->Print(Form("%s/Efficiencies_%s.pdf)", directory.c_str(),file.c_str()));
//	else  CanEff[i]->Print(Form("%s/Efficiencies_%s.pdf", directory.c_str(),file.c_str()));

	} //end i

	
//////////////////////////
//Reconstruction Rates////
/////////////////////////



  //   CanEffReco[i]->SaveAs(Form("%s/%sEfficienciesReco.eps", directory.c_str(),RecoInd[i].c_str()));
   //   if (i==NReco-1)CanEffReco[i]->Print(Form("%s/Efficiencies_%d.pdf)", directory.c_str(),type));
    // else CanEffReco[i]->Print(Form("%s/Efficiencies_%d.pdf", directory.c_str(),type));




} //end function



