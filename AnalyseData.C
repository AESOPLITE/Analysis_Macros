////////////////////////////////////////////////////////////////////////////////////////// 
///    Authors: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 14 , 2017
///		&&
///            Sarah Mechbal, smechbal@ucsc.edu
///            Department of Physics, University of California, Santa Cruz, January 22, 2018 
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "TChain.h"
#include "TPaveStats.h"
#include "TGaxis.h"


ClassImp(ALTckhit)
ClassImp(ALEvent)

void AnalyseDataEvent(string s,string c,string RecoID, int geoconf, int ndf);
void DeflectionPR(string s,string c, int);
void DeflectionReco(string s, string c, int);
void PullPR(string s, string c, int geoconf);
void PullReco(string s, string c, int geoconf);
void EReco(string s, string c, int geoconf);

void AnalyseDataEvent(string s,string c, string RecoID, int geoconf, int ndf)
{
  //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("/home/sarah/AESOPLITE/ANALYSIS_SOFT/src/ALSim/Dataparameters%d.dat",geoconf); 

 LoadDataparameters(paramfile,zL,OffsetLL,OffsetRL,TrigThresh);

 for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", zL:" << zL[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 cout << "T1 threshold: " << TrigThresh[0] <<endl;
 cout << "T2 threshold: " << TrigThresh[1] <<endl;
 cout << "T3 threshold: " << TrigThresh[2] <<endl;
 cout << "T4 threshold: " << TrigThresh[3] <<endl;
 cout << "Guard threshold: " << TrigThresh[4] <<endl;

 int ndfcut = ndf;
 string file[1]={""};
 string coinc[1]={""};
 string ID = RecoID;
 file[0]=s;
 coinc[0]=c;   
 
 //////////////////////////// 
//Resolution histograms
//////////////////////////// 

 TH1F**XResoPR=new TH1F*[7];
 TH1F**YResoPR=new TH1F*[7];
 TH1F**ZResoPR=new TH1F*[7];

 TH1F**XResoReco=new TH1F*[7];
 TH1F**YResoReco=new TH1F*[7];
 TH1F**ZResoReco=new TH1F*[7];


 
 for(int j=0;j<7;j++)
   {
    XResoPR[j]=new TH1F(Form("Xreso L%d from PR",j),Form("Xreso L%d from PR",j),200,-1,1); 
    YResoPR[j]=new TH1F(Form("Yreso L%d from PR",j),Form("Yreso L%d from PR",j),200,-1,1);
    ZResoPR[j]=new TH1F(Form("Zreso L%d from PR",j),Form("Zreso L%d from PR",j),600,-3,3);

    XResoReco[j]=new TH1F(Form("Xreso L%d from Reco",j),Form("Xreso L%d from Reco",j),2000,-10,10); 
    YResoReco[j]=new TH1F(Form("Yreso L%d from Reco",j),Form("Yreso L%d from Reco",j),2000,-10,10);
    ZResoReco[j]=new TH1F(Form("Zreso L%d from Reco",j),Form("Zreso L%d from Reco",j),2000,-30,30);
   }
 
 ////////////////////////////////
 //Deflection histograms from PR
 ////////////////////////////////
 
 TH1F*DefT1T3T4noT2=new TH1F("Deflection  T1&T3&T4&noT2","Deflection  T1&T3&T4&noT2",80,-0.2,0.2); 
 DefT1T3T4noT2->GetXaxis()->SetTitle("#Delta#Theta=#Theta_{out}-#Theta_{in} in rad");
 DefT1T3T4noT2->GetYaxis()->SetTitle("entries");
 DefT1T3T4noT2->SetLineColor(kRed);
 DefT1T3T4noT2->SetLineWidth(1);
 
 TH1F*DefT1T3T4T2=new TH1F("Deflection  T1&T3&T4&T2","Deflection  T1&T3&T4&T2",80,-0.2,0.2); 
 DefT1T3T4T2->GetXaxis()->SetTitle("#Delta#Theta=#Theta_{out}-#Theta_{in} in rad");
 DefT1T3T4T2->GetYaxis()->SetTitle("entries");
 DefT1T3T4T2->SetLineColor(kBlue);
 DefT1T3T4T2->SetLineWidth(1);

 
 ////////////////////////////////
 //Deflection histograms from KF
 ////////////////////////////////
 
 TH1F*DefT1T3T4noT2KF=new TH1F("Deflection  T1&T3&T4&noT2 KF","Deflection  T1&T3&T4&noT2 KF",80,-0.2,0.2); 
 DefT1T3T4noT2KF->GetXaxis()->SetTitle("#Delta#Theta=#Theta_{out}-#Theta_{in} in rad");
 DefT1T3T4noT2KF->GetYaxis()->SetTitle("entries");
 DefT1T3T4noT2KF->SetLineColor(kRed);
 DefT1T3T4noT2KF->SetLineWidth(1);
 
 TH1F*DefT1T3T4T2KF=new TH1F("Deflection  T1&T3&T4&T2 KF","Deflection  T1&T3&T4&T2 KF",80,-0.2,0.2); 
 DefT1T3T4T2KF->GetXaxis()->SetTitle("#Delta#Theta=#Theta_{out}-#Theta_{in} in rad");
 DefT1T3T4T2KF->GetYaxis()->SetTitle("entries");
 DefT1T3T4T2KF->SetLineColor(kBlue);
 DefT1T3T4T2KF->SetLineWidth(1);
 
 
 ///////////////////////////////////////////////////////////
 //Reconstructed momentum (independent of mass) from PR & KF
 ///////////////////////////////////////////////////////////
 
 TH1F*pRecoT1TT3T4noT2=new TH1F("pReco T1&T3&T4&NoT2","pReco  T1&T3&T4&NoT2",500,0.,5000); 
 pRecoT1TT3T4noT2->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoT1TT3T4noT2->GetYaxis()->SetTitle("entries");
 pRecoT1TT3T4noT2->SetLineColor(kBlue);
 pRecoT1TT3T4noT2->SetLineWidth(1);
 
  TH1F*pPRT1TT3T4noT2=new TH1F("pPR T1&T3&T4&NoT2","pPR  T1&T3&T4&NoT2",500,0.,5000); 
 pPRT1TT3T4noT2->GetXaxis()->SetTitle("Reconstructed momentum from PR in MeV/c");
 pPRT1TT3T4noT2->GetYaxis()->SetTitle("entries");
 pPRT1TT3T4noT2->SetLineColor(kRed);
 pPRT1TT3T4noT2->SetLineWidth(1);
 
  
 TH1F*pRecoT1TT3T4noT2cut=new TH1F("pReco T1&T3&T4&NoT2 cut","pReco  T1&T3&T4&NoT2 cut",500,0.,5000); 
 pRecoT1TT3T4noT2cut->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoT1TT3T4noT2cut->GetYaxis()->SetTitle("entries");
 pRecoT1TT3T4noT2cut->SetLineColor(kBlue);
 pRecoT1TT3T4noT2cut->SetLineWidth(1);
 
 TH1F*pRecoT1TT3T4T2=new TH1F("pReco  T1&T3&T4&T2","pReco  T1&T3&T4&T2",500,0,5000); 
 pRecoT1TT3T4T2->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoT1TT3T4T2->GetYaxis()->SetTitle("entries");
 pRecoT1TT3T4T2->SetLineColor(kRed);
 pRecoT1TT3T4T2->SetLineWidth(1);
 
  TH1F*pPRT1TT3T4T2=new TH1F("pPR T1&T3&T4&T2","pPR  T1&T3&T4&T2",500,0.,5000); 
 pPRT1TT3T4T2->GetXaxis()->SetTitle("Reconstructed momentum from PR in MeV/c");
 pPRT1TT3T4T2->GetYaxis()->SetTitle("entries");
 pPRT1TT3T4T2->SetLineColor(kRed);
 pPRT1TT3T4T2->SetLineWidth(1);
 
 TH1F*pRecoT1TT3T4T2cut=new TH1F("pReco  T1&T3&T4&T2 cut","pReco  T1&T3&T4&T2 cut",500,0,5000); 
 pRecoT1TT3T4T2cut->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoT1TT3T4T2cut->GetYaxis()->SetTitle("entries");
 pRecoT1TT3T4T2cut->SetLineColor(kBlue);
 pRecoT1TT3T4T2cut->SetLineWidth(1);

 TH1F*pRecoElecT2=new TH1F("pReco  T1&T3&T4&T2 elec","pReco  T1&T3&T4&T2 elec",500,0,5000);
 pRecoElecT2->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoElecT2->GetYaxis()->SetTitle("entries");
 pRecoElecT2->SetLineColor(kBlue);
 pRecoElecT2->SetLineWidth(1);

 TH1F*pRecoPosT2=new TH1F("pReco  T1&T3&T4&T2 pos","pReco  T1&T3&T4&T2 pos",500,0,5000);
 pRecoPosT2->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoPosT2->GetYaxis()->SetTitle("entries");
 pRecoPosT2->SetLineColor(kRed);
 pRecoPosT2->SetLineWidth(1);

 TH1F*pPRElecT2=new TH1F("pPR T1&T3&T4&T2 elec","pPR  T1&T3&T4&T2 elec",500,0,5000);
 pPRElecT2->SetLineColor(kBlue);
 pPRElecT2->SetLineWidth(1);
 pPRElecT2->SetLineStyle(10);

 TH1F*pPRPosT2=new TH1F("pPR  T1&T3&T4&T2 pos","pPR  T1&T3&T4&T2 pos",500,0,5000);
 pPRPosT2->SetLineColor(kRed);
 pPRPosT2->SetLineWidth(1);
 pPRPosT2->SetLineStyle(10);

 TH1F*pRecoElecnoT2=new TH1F("pReco  T1&T3&T4&noT2 elec","pReco  T1&T3&T4&noT2 elec",500,0,5000);
 pRecoElecnoT2->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoElecnoT2->GetYaxis()->SetTitle("entries");
 pRecoElecnoT2->SetLineColor(kBlue);
 pRecoElecnoT2->SetLineWidth(1);

 TH1F*pRecoPosnoT2=new TH1F("pReco  T1&T3&T4&noT2 pos","pReco  T1&T3&T4&noT2 pos",500,0,5000);
 pRecoPosnoT2->GetXaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoPosnoT2->GetYaxis()->SetTitle("entries");
 pRecoPosnoT2->SetLineColor(kRed);
 pRecoPosnoT2->SetLineWidth(1);

 TH1F*pPRElecnoT2=new TH1F("pPR T1&T3&T4&noT2 elec","pPR  T1&T3&T4&noT2 elec",500,0,5000);
 pPRElecnoT2->SetLineColor(kBlue);
 pPRElecnoT2->SetLineWidth(1);
 pPRElecnoT2->SetLineStyle(10);

 TH1F*pPRPosnoT2=new TH1F("pPR  T1&T3&T4&noT2 pos","pPR  T1&T3&T4&noT2 pos",500,0,5000);
 pPRPosnoT2->SetLineColor(kRed);
 pPRPosnoT2->SetLineWidth(1);
 pPRPosnoT2->SetLineStyle(10);
	
 ////////////////////////////////////////////////////////////
 //TH2F Reconstructed momentum as a function of chi2 and ndf
 ////////////////////////////////////////////////////////////
 
 TH2F*pRecoVSchi2=new TH2F("pReco vs chi2", "pReco vs chi2",30,0,30,5000,0,5000);
 pRecoVSchi2->GetXaxis()->SetTitle("chi2 from KF");
 pRecoVSchi2->GetYaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoVSchi2->SetTitle("Reconstructed momentum as a function of chi2 T1&T3&T4&noT2");
 
  
 TH2F*pRecoVSndf=new TH2F("pReco vs ndf", "pReco vs ndf",20,-10,10,5000,0,5000);
 pRecoVSndf->GetXaxis()->SetTitle("ndf from KF");
 pRecoVSndf->GetYaxis()->SetTitle("Reconstructed momentum in MeV/c");
 pRecoVSndf->SetTitle("Reconstructed momentum as a function of ndf T1&T3&T4&noT2");

///////////////////////////////////
// chi2 of KF fit /////////////////
///////////////////////////////////
 TH1F*KFchi2=new TH1F("KF chi2","KF chi2",1000,0,20);
 KFchi2->GetXaxis()->SetTitle("chi2");
 KFchi2->GetYaxis()->SetTitle("entries");
 KFchi2->SetLineColor(kBlack);
 KFchi2->SetLineWidth(1);

 TH1F*KFchi2cut=new TH1F("KF chi2 cut","KF chi2 cut",1000,0,20);
 KFchi2cut->GetXaxis()->SetTitle("chi2");
 KFchi2cut->GetYaxis()->SetTitle("entries");
 KFchi2cut->SetLineColor(kBlack);
 KFchi2cut->SetLineWidth(1);



 
 // Make cut on hypothetical particle type
 // Compare reconstructed energy from PR and KF
 
//ROOT output file 
TFile*fileout=new TFile(Form("/home/sarah/AESOPLITE/GroundRuns/Sweden2018/%s_%s_RecoAnalysis.root", file[0].c_str(),ID.c_str()),"RECREATE"); 
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");


chain->Add("/home/sarah/AESOPLITE/GroundRuns/Sweden2018/NLAll_PRonly.root");

 
 //Define variables to read event
 ALEvent *e = new ALEvent();      
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 int nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;  
      
        double PosPartnoT2=0;
        double NegPartnoT2=0;
        double PosPartT2=0;
        double NegPartT2=0;


 for (int j=0;j<nentries;j++)
   {
    chain->GetEntry(j); //Load the entry i in the variable e 

    if(j%10000==0) cout << "Event: " << j <<endl;
    
    int nnhits = (int) e->get_Nhits();
   	
	//apply external trigger requirements
	
	bool T1=e->get_T1();
	bool T2=e->get_T2();
	bool T3=e->get_T3();
	bool T4=e->get_T4();
        if(e->get_deflecPR()!=0) {
	//extract PR and reco variables
	double preco=e->get_p0reco();
	double Ekreco=e->get_Ekreco();
	double deflecPR=e->get_deflecPR();
	double chi2reco=e->get_chi2();
	double ndf=e->get_ndf();
	double chi2NBPR=e->get_chi2NBPR();
	double chi2BPR=e->get_chi2BPR();
	double pPR=e->get_p0PR()*1000;
    double slopeRecoL1=0;
    double thetaRecoL1=0;
    double slopeRecoL5=0;
    double thetaRecoL5=0;
	double deflecReco=0;
			
	//fill pull histograms (in mm)
	for(int k=0;k<nnhits;k++)
         {
		int Lindex=(int)e->get_hits().at(k)->get_L();

//record directional cosines at L1 and L5
	if(Lindex==1)
		      {								
	   //Reco variables at L1
	   if((int)e->get_hits().at(k)->get_fUsed()==1){
	   double cxL1Reco = e->get_hits().at(k)->get_cxreco();
	   double cyL1Reco = e->get_hits().at(k)->get_cyreco();
	   double czL1Reco = e->get_hits().at(k)->get_czreco();
	   slopeRecoL1= (cyL1Reco/czL1Reco);
	   thetaRecoL1 = TMath::ATan(slopeRecoL1); 	
			    }
		      }	
		//LAYER L5, B PLANE
	if(Lindex==5)
		  {						
		   //Reco variables at L5
			if((int)e->get_hits().at(k)->get_fUsed()==1) {
		   double cxL5Reco = e->get_hits().at(k)->get_cxreco();
		   double cyL5Reco = e->get_hits().at(k)->get_cyreco();
		   double czL5Reco = e->get_hits().at(k)->get_czreco();
		   slopeRecoL5= (cyL5Reco/czL5Reco);
		   thetaRecoL5 = TMath::ATan(slopeRecoL5); 	
		   deflecReco = thetaRecoL5 - thetaRecoL1;
		   }
	}
	//Coordinate Pulls
	      if(e->get_hits().at(k)->get_flagPR() && e->get_hits().at(k)->get_fUsed())
           {
			float tmp=0;
			int Lindex=(int)e->get_hits().at(k)->get_L();
			tmp=(e->get_hits().at(k)->get_x())*10;
			XResoReco[Lindex]->Fill(tmp-e->get_hits().at(k)->get_xreco());
			tmp=(e->get_hits().at(k)->get_y())*10;
			YResoReco[Lindex]->Fill(tmp-e->get_hits().at(k)->get_yreco());
			tmp=(e->get_hits().at(k)->get_z())*10;
			ZResoReco[Lindex]->Fill(tmp-e->get_hits().at(k)->get_zreco());
                   } //if flagPR && if used by KF	 
		     if(e->get_hits().at(k)->get_flagPR())
           {
			float tmp=0;
			tmp=(e->get_hits().at(k)->get_x())*10;
			XResoPR[Lindex]->Fill(tmp-(e->get_hits().at(k)->get_xPR()*10));
			tmp=(e->get_hits().at(k)->get_y())*10;
			YResoPR[Lindex]->Fill(tmp-(e->get_hits().at(k)->get_yPR()*10));
			tmp=(e->get_hits().at(k)->get_z())*10;
			ZResoPR[Lindex]->Fill(tmp-(e->get_hits().at(k)->get_zPR()*10));
		   } //if flagPR
		}	//end k

/////////////////////////////////////////////
//FILL HISTOGRAMS///////////////////////////
///////////////////////////////////////////*
			
			
	//fill T1&T3&T4&NoT2 histograms
    
   if(T1 && T3 && T4 && (!T2)) 
   {
	   pRecoT1TT3T4noT2->Fill(preco);
	   DefT1T3T4noT2->Fill(deflecPR);
	   pPRT1TT3T4noT2->Fill(pPR);
	  	  DefT1T3T4noT2KF->Fill(deflecReco);
	  if (Ekreco!=-999) KFchi2->Fill(chi2reco);
           if(deflecPR>0) {
			   PosPartnoT2++;
			   pRecoPosnoT2->Fill(preco);
			   pPRPosnoT2->Fill(pPR);
		   }
          else {
			  pRecoElecnoT2->Fill(preco);
			  pPRElecnoT2->Fill(pPR);
			  NegPartnoT2++;
		  }
	   if (ndf>ndfcut) {
		 pRecoT1TT3T4noT2cut->Fill(preco);
         KFchi2cut->Fill(chi2reco);
	}       
   }
   
   	//fill T1&T3&T4&T2 histograms

	if(T1 && T2 && T3 && T4) 
   {
	   pRecoT1TT3T4T2->Fill(preco);
 	   DefT1T3T4T2->Fill(deflecPR);
	   pPRT1TT3T4T2->Fill(pPR);
	  pRecoVSchi2->Fill(chi2reco,preco);
	  pRecoVSndf->Fill(ndf,preco);
	  DefT1T3T4T2KF->Fill(deflecReco);
          if(Ekreco!=-999) KFchi2->Fill(chi2reco);
	  if(deflecPR>0) PosPartT2++;
	  else NegPartT2++;
	  if(deflecPR>0) {
		 pRecoPosT2->Fill(preco);
		 pPRPosT2->Fill(pPR);
		}
	  else if(deflecPR<0) {
		 pRecoElecT2->Fill(preco);
		 pPRElecT2->Fill(pPR);
		} 
	  if(ndf>ndfcut) {
		 pRecoT1TT3T4T2cut->Fill(preco);
	         KFchi2cut->Fill(chi2reco);
	}
   }
  
   } //if deflecPR!=0
  } //end loop on events
		 
/////////////////////
//	DISPLAY
/////////////////////

  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(11111111);
 
//////////////////////////////////////////////////////////////////
// Plot momentum for all events that have T1 & T3 & T4
  
	 TCanvas*CanpRecoT1TT3T4=new TCanvas("CanpRecoT1TT3T4","CanpRecoT1TT3T4",200,10,1200,800);
    CanpRecoT1TT3T4->cd(1);
    TLegend*legT1T3T4=new TLegend(0.1,0.9,0.2,1.0);
    legT1T3T4->SetBorderSize(0);
    legT1T3T4->SetFillStyle(0);
    legT1T3T4->AddEntry(pRecoT1TT3T4noT2,"No T2 signal","l");
    legT1T3T4->AddEntry(pRecoT1TT3T4T2,"With T2 signal","l");
    float MaxT1T3T4=0;
    float tmp1T1T3T4=pRecoT1TT3T4noT2->GetBinContent(pRecoT1TT3T4noT2->GetMaximumBin());
    float tmp2T1T3T4=pRecoT1TT3T4T2->GetBinContent(pRecoT1TT3T4T2->GetMaximumBin());
    if(tmp1T1T3T4>tmp2T1T3T4) MaxT1T3T4=tmp1T1T3T4;
    else MaxT1T3T4=tmp2T1T3T4;
    pRecoT1TT3T4noT2->SetMaximum(1.05*MaxT1T3T4);
    pRecoT1TT3T4noT2->Draw("hist");
    pPRT1TT3T4T2->Draw("sames");
    legT1T3T4->Draw();
    pRecoT1TT3T4noT2->SetTitle("T1&T3&T4 all events");
    CanpRecoT1TT3T4->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf(", file[0].c_str(),ID.c_str()));
    fileout->cd();
    pRecoT1TT3T4noT2->Write();
    CanpRecoT1TT3T4->Write();
 /*
    TCanvas*CanpRecoT1TT3T4noT2ndfcut=new TCanvas("CanpRecoT1TT3T4noT2 ndfcut","CanpRecoT1TT3T4noT2 ndfcut",200,10,1200,800);
    CanpRecoT1TT3T4noT2ndfcut->cd(1);
    pRecoT1TT3T4noT2cut->Draw("hist");
    pRecoT1TT3T4noT2cut->SetTitle(Form("T1&T3&T4&NoT2: Reconstructed momentum ndf > %d", ndfcut));
    fileout->cd();
    pRecoT1TT3T4noT2cut->Write();
    CanpRecoT1TT3T4noT2ndfcut->Write();
    CanpRecoT1TT3T4noT2ndfcut->Print(Form("Data/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
///////////////////////////////////////////////////////////////////
    TCanvas*CanpRecoT1TT3T4T2=new TCanvas("CanpRecoT1TT3T4T2","CanpRecoT1TT3T4T2",200,10,1200,800);
    CanpRecoT1TT3T4T2->cd(1);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    TPaveText*ptT2=new TPaveText(0.7,0.8,0.8,0.9,"NDC");
    ptT2->AddText(Form("#frac{positive particles}{negative particles}: = %f ", (PosPartT2/NegPartT2)));
    ptT2->SetBorderSize(0);
     TLegend*leg1=new TLegend(0.1,0.9,0.2,1.0);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(pRecoT1TT3T4T2,"From KF reco","l");
    leg1->AddEntry(pPRT1TT3T4T2,"From PR reco","l");
    float Max1=0;
    float tmp11=pRecoT1TT3T4T2->GetBinContent(pRecoT1TT3T4T2->GetMaximumBin());
    float tmp12=pPRT1TT3T4T2->GetBinContent(pPRT1TT3T4T2->GetMaximumBin());
    if(tmp11>tmp12) Max1=tmp11;
    else Max1=tmp12;
    pRecoT1TT3T4T2->SetMaximum(1.05*Max1);
    pRecoT1TT3T4T2->Draw("hist");
    pPRT1TT3T4T2->Draw("same");
    leg1->Draw();
    ptT2->Draw();
    pRecoT1TT3T4T2->SetTitle("T1&T1&T3&T4&T2: Reconstructed momentum");
    fileout->cd();
    pRecoT1TT3T4T2->Write();
    CanpRecoT1TT3T4T2->Write();
    CanpRecoT1TT3T4T2->Print(Form("Data/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
*/ 
///////////////////////////////////////////////////////////////////////////////////////////////////
//Log-Log T1T2T3T4, separates histograms for positive/negative tracks 
////////////////////////////////////////////////////////////////////////////////////////////////// 
    TCanvas*CanpRecoT1TT3T4noT2=new TCanvas("CanpRecoT1TT3T4noT2","CanpRecoT1TT3T4noT2",200,10,1200,800);
    //gPad->SetLogy(1);
    //gPad->SetLogx(1);
    CanpRecoT1TT3T4noT2->cd(1);
    TLegend*leg0=new TLegend(0.1,0.9,0.2,1.0);
    TPaveText*ptnoT2=new TPaveText(0.7,0.8,0.8,0.9,"NDC");
    ptnoT2->AddText(Form("#frac{positive particles}{negative particles}: = %f ", (pRecoPosnoT2->GetEntries())/(pRecoElecnoT2->GetEntries())));
    ptnoT2->SetBorderSize(0);
	ptnoT2->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetFillStyle(0);
    leg0->AddEntry(pRecoElecnoT2,"From KF reco, negative tracks","l");
    leg0->AddEntry(pRecoPosnoT2,"From KF reco, positive tracks","l");
    float Max0=0;
    float tmp01=pRecoElecnoT2->GetBinContent(pRecoElecnoT2->GetMaximumBin());
    float tmp02=pRecoPosnoT2->GetBinContent(pRecoPosnoT2->GetMaximumBin());
    if(tmp01>tmp02) Max0=tmp01;
    else Max0=tmp02;
    pRecoElecnoT2->SetMaximum(1.05*Max0);
    pRecoElecnoT2->Draw("hist");
    pRecoPosnoT2->Draw("same");
    leg0->Draw();
    ptnoT2->Draw();
    pRecoElecnoT2->SetTitle("T1&T3&T4&NoT2: Reconstructed momentum");
    CanpRecoT1TT3T4noT2->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
    fileout->cd();
    pRecoElecnoT2->Write();
    CanpRecoT1TT3T4noT2->Write();
	
//with T2
    TCanvas*CanpRecoT1TT3T4T2=new TCanvas("CanpRecoT1TT3T4T2","CanpRecoT1TT3T4T2",200,10,1200,800);
   // gPad->SetLogy(1);
   // gPad->SetLogx(1);
    CanpRecoT1TT3T4T2->cd(1);
    TLegend*leg1=new TLegend(0.1,0.9,0.2,1.0);
    TPaveText*ptT2=new TPaveText(0.7,0.8,0.8,0.9,"NDC");
    ptT2->AddText(Form("#frac{positive particles}{negative particles}: = %f ", (pRecoPosT2->GetEntries())/(pRecoElecT2->GetEntries())));
    ptT2->SetBorderSize(0);
	ptT2->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(pRecoElecT2,"From KF reco, negative tracks","l");
    leg1->AddEntry(pRecoPosT2,"From KF reco, positive tracks","l");
    float Max1=0;
    float tmp11=pRecoElecT2->GetBinContent(pRecoElecT2->GetMaximumBin());
    float tmp12=pRecoPosT2->GetBinContent(pRecoPosT2->GetMaximumBin());
    if(tmp11>tmp12) Max1=tmp11;
    else Max1=tmp12;
    pRecoElecT2->SetMaximum(1.05*Max1);
    pRecoElecT2->Draw("hist");
    pRecoPosT2->Draw("same");
    leg1->Draw();
    ptT2->Draw();
    pRecoElecT2->SetTitle("T1&T3&T4&T2: Reconstructed momentum");
    CanpRecoT1TT3T4T2->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
    fileout->cd();
    pRecoElecT2->Write();
    CanpRecoT1TT3T4T2->Write();
///////////////////////////////////////////////////////////////////////////////////////////////////	 
	 
    TCanvas*CanpRecoT1TT3T4noT2ndfcut=new TCanvas("CanpRecoT1TT3T4noT2 ndfcut","CanpRecoT1TT3T4noT2 ndfcut",200,10,1200,800);
    CanpRecoT1TT3T4noT2ndfcut->cd(1);
    pRecoT1TT3T4noT2cut->Draw("hist");
    pRecoT1TT3T4noT2cut->SetTitle(Form("T1&T3&T4&NoT2: Reconstructed momentum ndf > %d", ndfcut));
    fileout->cd();
    pRecoT1TT3T4noT2cut->Write();
    CanpRecoT1TT3T4noT2ndfcut->Write();
    CanpRecoT1TT3T4noT2ndfcut->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
	
		ofstream toTeX;
	toTeX.open("TestToTeX.txt");
	toTeX << "Hello World " << endl;
	toTeX.close();
	 
////////////////////////////////////////////////////////////////////////////////////////////////////////////:
    
    TCanvas*CanpRecoT1TT3T4T2ndfcut=new TCanvas("CanpRecoT1TT3T4T2 ndfcut","CanpRecoT1TT3T4T2 ndfcut",200,10,1200,800);
    CanpRecoT1TT3T4T2ndfcut->cd(1);
    pRecoT1TT3T4T2cut->Draw("hist");
    pRecoT1TT3T4T2cut->SetTitle(Form("T1&T3&T4&T2: Reconstructed momentum ndf > %d", ndfcut));
    fileout->cd();
    pRecoT1TT3T4T2cut->Write();
    CanpRecoT1TT3T4T2ndfcut->Write();
    CanpRecoT1TT3T4T2ndfcut->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
 
    TCanvas*CanpRecoT1TT3T4T2cut=new TCanvas("CanpRecoT1TT3T4T2cut","CanpRecoT1TT3T4T2cut",200,10,1200,800);
    CanpRecoT1TT3T4T2cut->cd(1);
    pRecoT1TT3T4T2->GetXaxis()->SetRangeUser(0,200);
    pRecoT1TT3T4T2->Draw("hist");
    pRecoT1TT3T4T2->SetTitle("T1&T3&T4&T2: Reconstructed momentum");
    fileout->cd();
    pRecoT1TT3T4T2cut->Write();
    CanpRecoT1TT3T4T2cut->Write();
    CanpRecoT1TT3T4T2cut->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
/*
    TCanvas*CanpRecoElecT2=new TCanvas("CanpRecoT1TT3T4T2Elec","CanpRecoT1TT3T4T2Elec",200,10,1200,800);
    CanpRecoElecT2->cd(1);
    TLegend*legElec=new TLegend(0.1,0.9,0.2,1.0);
    legElec->SetBorderSize(0);
    legElec->SetFillStyle(0);
    legElec->AddEntry(pRecoElecT2,"From KF reco, e^{-}","l");
    legElec->AddEntry(pRecoPosT2,"From KF reco, e^{+}","l");
    legElec->AddEntry(pPRElecT2,"From PR reco, e^{-}","l");
    legElec->AddEntry(pPRPosT2,"From PR reco, e^{+}","l");
    float MaxElec=0;
    float tmpElec=pRecoElecT2->GetBinContent(pRecoElecT2->GetMaximumBin());
    float tmpPos=pRecoPosT2->GetBinContent(pRecoPosT2->GetMaximumBin());
    if(tmpElec>tmpPos) MaxElec=tmpElec;
    else MaxElec=tmpPos;
    pRecoT1TT3T4T2->SetMaximum(1.05*MaxElec);
    pRecoElecT2->Draw("hist");
    pRecoPosT2->Draw("sames");
    pPRElecT2->Draw("sames");
    pPRPosT2->Draw("sames");
    legElec->Draw();
    gPad->Update();
    TPaveStats *sElec = (TPaveStats*)pRecoElecT2->FindObject("stats");
    sElec->SetTextColor(kBlue);
    double y1 = sElec->GetY1NDC();
    double y2 = sElec->GetY2NDC();
    gPad->Update();
    TPaveStats *sPos = (TPaveStats*)pRecoPosT2->FindObject("stats");
    gPad->Update();
    sPos->SetTextColor(kRed);
    sPos->SetY2NDC(y1);
    sPos->SetY1NDC(y1-fabs(y2-y1));
    gPad->Update();

    pRecoElecT2->SetTitle("T1&T3&T4&T2: Reconstructed momentum e^{+} & e^{-}");
    fileout->cd();
    pRecoElecT2->Write(); 
    pRecoPosT2->Write();
    CanpRecoElecT2->Write();
    CanpRecoElecT2->Print(Form("Data/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
  */  
    
    TCanvas*CanpRecoVSchi2=new TCanvas("CanpRecoVSchi2","CanpRecoVSchi2",200,10,1200,800);
    CanpRecoVSchi2->cd(1);
    pRecoVSchi2->Draw("colz");
    fileout->cd();
    pRecoVSchi2->Write();
    CanpRecoVSchi2->Write();
    CanpRecoVSchi2->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
    
    TCanvas*CanpRecoVSndf=new TCanvas("CanpRecoVSndf","CanpRecoVSndf",200,10,1200,800);
    CanpRecoVSndf->cd(1);
    pRecoVSndf->Draw("colz");
    fileout->cd();
    pRecoVSndf->Write();
    CanpRecoVSndf->Write();
    CanpRecoVSndf->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));

    TCanvas*CanKFchi2=new TCanvas("CanKFchi2","CanKFchi2",200,10,1200,800);
    CanKFchi2->cd(1);
    KFchi2->Draw("hist");
    KFchi2->SetTitle("chi2 from KF fit, all reconstructed event");
    fileout->cd();
    KFchi2->Write();
    CanKFchi2->Write();
    CanKFchi2->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));

    TCanvas*CanKFchi2cut=new TCanvas("CanKFchi2cut","CanKFchi2cut",200,10,1200,800);
    CanKFchi2cut->cd(1);
    KFchi2cut->Draw("hist");
    KFchi2cut->SetTitle(Form("chi2 from KF fit, reconstructed event with ndf > %d",ndfcut));
    fileout->cd();
    KFchi2cut->Write();
    CanKFchi2cut->Write();
    CanKFchi2cut->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));

	
////////////////////////////////////////////////////////////////////	
	TCanvas*CanDeflecT1T3T4noT2=new TCanvas("CanDeflecT1T3T4noT2","CanDeflecT1T3T4noT2",200,10,1200,800);
	TLegend*leg=new TLegend(0.0,0.7,0.45,1.0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(DefT1T3T4noT2,"T1&T3&T4&noT2","l");
	leg->AddEntry(DefT1T3T4T2,"T1&T3&T4&T2","l");
        CanDeflecT1T3T4noT2->cd(1);
	float Max=0;
	float tmp1=DefT1T3T4noT2->GetBinContent(DefT1T3T4noT2->GetMaximumBin());
	float tmp2=DefT1T3T4T2->GetBinContent(DefT1T3T4T2->GetMaximumBin());
	if(tmp1>tmp2) Max=tmp1;
	else Max=tmp2;
	DefT1T3T4noT2->SetMaximum(1.05*Max);
	DefT1T3T4noT2->Draw("hist");
	DefT1T3T4T2->Draw("same");
	leg->Draw();
        CanDeflecT1T3T4noT2->SetTitle("T1&T1&T3&T4&NoT2: Deflection from PR");
        fileout->cd();
        DefT1T3T4noT2->Write();
	DefT1T3T4T2->Write();
        CanDeflecT1T3T4noT2->Write();
	CanDeflecT1T3T4noT2->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
/////////////////////////////////////////////////////
////////////////////Deflection from KF //////////////
////////////////////////////////////////////////////
	 
	TCanvas*CanDeflecT1T3T4T2KF=new TCanvas("CanDeflecT1T3T4T2 KF","CanDeflecT1T3T4T2 KF",200,10,1200,800);
	CanDeflecT1T3T4T2KF->cd(1);
        CanDeflecT1T3T4T2KF->SetTitle("T1&T1&T3&T4&T2: Deflection from KF");
	CanDeflecT1T3T4T2KF->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
	TCanvas*CanDeflecT1T3T4noT2KF=new TCanvas("CanDeflecT1T3T4noT2 KF","CanDeflecT1T3T4noT2 KF",200,10,1200,800);
	TLegend*legKF=new TLegend(0.0,0.7,0.45,1.0);
	legKF->SetBorderSize(0);
	legKF->SetFillStyle(0);
	legKF->AddEntry(DefT1T3T4noT2,"T1&T3&T4&noT2","l");
	legKF->AddEntry(DefT1T3T4T2,"T1&T3&T4&T2","l");
        CanDeflecT1T3T4noT2KF->cd(1);
	float MaxKF=0;
	float tmp1KF=DefT1T3T4noT2->GetBinContent(DefT1T3T4noT2KF->GetMaximumBin());
	float tmp2KF=DefT1T3T4T2->GetBinContent(DefT1T3T4T2KF->GetMaximumBin());
	if(tmp1KF>tmp2KF) MaxKF=tmp1KF;
	else MaxKF=tmp2KF;
	DefT1T3T4noT2->SetMaximum(1.05*MaxKF);
	DefT1T3T4noT2KF->Draw("hist");
	DefT1T3T4T2KF->Draw("same");
	legKF->Draw();
	 
        CanDeflecT1T3T4noT2KF->SetTitle("T1&T1&T3&T4&NoT2: Deflection from KF");
        fileout->cd();
        DefT1T3T4noT2->Write();
	DefT1T3T4T2->Write();
        CanDeflecT1T3T4noT2KF->Write();
	CanDeflecT1T3T4noT2KF->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
	
	TCanvas*CanDeflecT1T3T4T2=new TCanvas("CanDeflecT1T3T4T2","CanDeflecT1T3T4T2",200,10,1200,800);
	CanDeflecT1T3T4T2->cd(1);
	DefT1T3T4T2->Draw("hist");
        CanDeflecT1T3T4T2->SetTitle("T1&T1&T3&T4&T2: Deflection from PR");
	CanDeflecT1T3T4T2->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));	
	///XYZ resolution pulls
 double nsigma =4;
 TCanvas**CanXYZReso=new TCanvas*[7];
 for(int i=0;i<7;i++) 
 {
     CanXYZReso[i]=new TCanvas(Form("XYZReso Layer %d", i),Form("XYZReso Layer %d",i),200,10,1200,800);     
     CanXYZReso[i]->Divide(2,2);
     
   if(i==0 || i==4 || i==6)
	   {
	CanXYZReso[i]->cd(1);
	XResoReco[i]->SetLineColor(kRed);
	XResoReco[i]->SetTitle(Form("Pull X from KF, Tracker layer %d (non-bending)",i));
	XResoReco[i]->GetXaxis()->SetRangeUser(XResoReco[i]->GetMean()-nsigma*XResoReco[i]->GetRMS(),XResoReco[i]->GetMean()+nsigma*XResoReco[i]->GetRMS());
	XResoReco[i]->GetXaxis()->SetTitle("data-reco (mm)");
	XResoReco[i]->Draw("hist");
	
	CanXYZReso[i]->cd(3);
	XResoPR[i]->SetLineColor(kRed);
	XResoPR[i]->SetTitle(Form("Pull X from PR, Tracker layer (non-bending) %d",i));
	XResoPR[i]->GetXaxis()->SetRangeUser(XResoPR[i]->GetMean()-nsigma*XResoPR[i]->GetRMS(),XResoPR[i]->GetMean()+nsigma*XResoPR[i]->GetRMS());
	XResoPR[i]->GetXaxis()->SetTitle("data-PR (mm)");
	XResoPR[i]->Draw("hist");
	   }
	else {
	CanXYZReso[i]->cd(1);
	YResoReco[i]->SetLineColor(kRed);
	YResoReco[i]->SetTitle(Form("Pull Y from KF, Tracker layer %d (bending)",i));
	YResoReco[i]->GetXaxis()->SetRangeUser(YResoReco[i]->GetMean()-nsigma*YResoReco[i]->GetRMS(),YResoReco[i]->GetMean()+nsigma*YResoReco[i]->GetRMS());
	YResoReco[i]->GetXaxis()->SetTitle("data-reco (mm)");
	YResoReco[i]->Draw("hist");
	
	CanXYZReso[i]->cd(3);
	YResoPR[i]->SetLineColor(kRed);
	YResoPR[i]->SetTitle(Form("Pull X from PR, Tracker layer %d (bending)",i));
	YResoPR[i]->GetXaxis()->SetRangeUser(YResoPR[i]->GetMean()-nsigma*YResoPR[i]->GetRMS(),YResoPR[i]->GetMean()+nsigma*YResoPR[i]->GetRMS());
	YResoPR[i]->GetXaxis()->SetTitle("data-PR (mm)");
	YResoPR[i]->Draw("hist");
	}
	CanXYZReso[i]->cd(2);
	ZResoReco[i]->SetLineColor(kBlue);
	ZResoReco[i]->SetTitle(Form("Pull Z from KF, Tracker layer %d",i));
	ZResoReco[i]->GetXaxis()->SetRangeUser(ZResoReco[i]->GetMean()-nsigma*ZResoReco[i]->GetRMS(),ZResoReco[i]->GetMean()+nsigma*ZResoReco[i]->GetRMS());
	ZResoReco[i]->GetXaxis()->SetTitle("data-reco (mm)");
	ZResoReco[i]->Draw("hist");
	
	CanXYZReso[i]->cd(4);
	ZResoPR[i]->SetLineColor(kBlue);
	ZResoPR[i]->SetTitle(Form("Pull Z from PR, Tracker layer %d",i));
	ZResoPR[i]->GetXaxis()->SetRangeUser(ZResoPR[i]->GetMean()-nsigma*ZResoPR[i]->GetRMS(),ZResoPR[i]->GetMean()+nsigma*ZResoPR[i]->GetRMS());
	ZResoPR[i]->GetXaxis()->SetTitle("data-PR (mm)");
	ZResoPR[i]->Draw("hist");
   
        if(i<6)CanXYZReso[i]->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf", file[0].c_str(),ID.c_str()));
	if(i==6)  CanXYZReso[i]->Print(Form("GroundRuns/Reconstruction/%s_%s_RecoAnalysis.pdf)", file[0].c_str(),ID.c_str()));
	} //i
	
} //end function AnalyseDataEvent
