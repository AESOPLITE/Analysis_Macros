
////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, May , 2017
////   Modififed by Sarah Mechbal, smechbal@ucsc.edu
////   Department of Physics, University of California Santa Cruz, December 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadMCparameters.h"
#include "TChain.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void AnalyseTestOneEnergy(int t, int Ene, int ndf, int thetmin, int thetmax, int cycles, string s);

void AnalyseTestOneEnergy(int t, int Ene, int ndf, int thetmin, int thetmax, int cycles, string s)
{

  //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;

 //string MCparamfile="./MCparametersV1.dat"; 
string MCparamfile="./MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh);

 int NTConfig=1;//Number of configuration of triggers
 string TConfig[1]={"T1&T4"};
 //double TrigThresh[4]={0.3, 0.29, 0.57,0.1};		//trigger thresholds in MeV determined for T1/T3/T4/Guard
 double RestMass   = 0.51099;				//mass electron in MeV

 //Input file  
 string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V3/Disc"; 
 //filename structure
 string startfile="aesopliteNonUniB_V3";
 string endfile="_fort.99";
 string directory="Disc/V3/NonUniB/toshow";

 string ID = s;
	
 //type of particle
 int type=t;
 int Nene=1;
 //Energies
 int EneMuon = 1;
 int ene = Ene;
 double Enek;
 //Number of cycles per energy
 int Ncycles=cycles;  
//cut values of renconstructed degrees of freedom 
int ndf_cut = ndf;
float eMC_cut = 0.9;
//make cuts as a function of incident angle
int thetamin=thetmin;
int thetamax=thetmax; 	
//Define histogram outputs

 TH1F *EkReso=new TH1F;
 TH1F *EkPR=new TH1F;
 TH1F *EkResoCut = new TH1F;
 TH1F *EkrecoCutL0 = new TH1F;
 TH1F *eMCL0 = new TH1F;
 TH1F *inversep0reco=new TH1F;
 TH1F *inversep0cut=new TH1F;
 TH1F *inversep0PR=new TH1F;
 TH1F *inversep0CutL0=new TH1F;
 TH1F *chi2ndf=new TH1F;
 TH1F *Phi0Pull=new TH1F;
 TH1F *TanLPull=new TH1F;
 TH1F *KappaPull=new TH1F;
 TH1F ***LayerPull=new TH1F**;
 TH1F ***LayerPullPR=new TH1F**;
 TH1F **Occupancy = new TH1F*;
 

//Define scatter plot output
 TH2F *ErecovsEL0 = new TH2F;
 TH2F *chi2vsEkreco = new TH2F;
 TH2F *ndfvsEkreco = new TH2F;
 TH2F  *EkRecovsDeltaPhi0 = new TH2F;
 TH2F  *EkRecovsDeltaCpa = new TH2F; 
 TH2F  *EkRecovsDeltaTanL = new TH2F;
 TH2F  *CosZ0vsndf = new TH2F;
 TH2F  *CosZ0vsEkreco = new TH2F;
 //Create TChain to merge the files for a given energy
 TChain *chain=new TChain; 

//ROOT output file 
string UNIT;
float Enemin=0;
float Enemax=500;
int Nbins=1;

 if(type==3 || type==4) UNIT="MeV";
 if(type==10 || type ==11) {
 
UNIT="GeV";
Enemin=10;
Enemax=5000;
Nbins=5;
RestMass= 105.6584;				//mass muon in MeV
}

//Nbins=TMath::Nint(Nbins*(Ncycles/100));
TFile*fileout=new TFile(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.root", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax),"RECREATE");

 
   
    cout << "Energy: " << Ene << " " << UNIT <<endl;
    chain=new TChain("MC");
    
   for(int j=0;j<Ncycles;j++)//Number of cycles
      
      { 
		  
     // 
  //chain->Add(Form("%s/%d/RecoEvent_%s_%d_%dGeV%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,EneMuon,j+1,endfile.c_str(),ID.c_str()));
     chain->Add(Form("%s/%d/RecoEvent_%s_%d_%d%s%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,UNIT.c_str(),j+1,endfile.c_str(),ID.c_str()));

	  } 
    
   
    //Define histograms 
        EkReso=new TH1F(Form("Ekreco %d %s",Ene,UNIT.c_str()),Form("Ekreco %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),Enemin,Enemax);
	EkPR=new TH1F(Form("EkPR %d %s",Ene,UNIT.c_str()),Form("EkPR %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),Enemin,Enemax);
	EkrecoCutL0 = new TH1F(Form("Ekreco ndf&eMC cut %d %s",Ene,UNIT.c_str()),Form("Ekreco ndf&eMC cut %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),Enemin,Enemax);
	EkResoCut=new TH1F(Form("Ekreco chi2 cut %d %s",Ene,UNIT.c_str()),Form("Ekreco chi2 cut %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),Enemin,Enemax);
        eMCL0 = new TH1F("eMC at L0", "eMC at L0", Nbins*(Enemax-Enemin),Enemin,Enemax);
	inversep0reco=new TH1F(Form("1/p0reco %d %s",Ene,UNIT.c_str()),Form("1/p0reco %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),1./Enemin, 1./Enemax);
	inversep0PR=new TH1F(Form("1/p0PR %d %s",Ene,UNIT.c_str()),Form("1/p0PR %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),1./Enemin, 1./Enemax);
	inversep0CutL0 = new TH1F(Form("1/p0 ndf&eMC cut %d %s",Ene,UNIT.c_str()),Form("1/p0 ndf&eMC %d %s",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),1./Enemin, 1.0/Enemax);
	inversep0cut=new TH1F(Form("1/p0reco %d %s cut",Ene,UNIT.c_str()),Form("1/p0reco %d %s cut",Ene,UNIT.c_str()),Nbins*(Enemax-Enemin),1./Enemin, 1./Enemax);
	Phi0Pull = new TH1F(Form("#Helix parameter #phi_0 error"),Form("Helix parameter #phi_{0}"),500,-5,5);
        TanLPull = new TH1F(Form("#Helix parameter Tan#lambda error"),Form("Helix parameter Tan#lambda"),500,-5,5);
        KappaPull = new TH1F(Form("#Helix parameter #kappa error"),Form("Helix parameter #kappa"),500,-5,5);
        chi2ndf=new TH1F(Form("chi2ndf %d %s",Ene,UNIT.c_str()),Form("chi2ndf %d %s", Ene,UNIT.c_str()),100,0,1);
        chi2vsEkreco = new TH2F(Form("Ekreco vs chi2/ndf %d %s", Ene,UNIT.c_str()),Form("Ekreco vs chi2/ndf %d %s", Ene,UNIT.c_str()), 100, -15, 15, Nbins*(Enemax-Enemin), Enemin, 5*Enemax); 
        ErecovsEL0  = new TH2F("Ekreco vs eMC at L0", "Ekreco vs eMC at L0", Nbins*(Enemax-Enemin), Enemin, Enemax, Nbins*(Enemax-Enemin), Enemin, Enemax);
	ndfvsEkreco =  new TH2F(Form("Ekreco vs ndf %d %s", Ene,UNIT.c_str()),Form("Ekreco vs ndf %d %s",Ene, UNIT.c_str()), 10, -10, 10, Nbins*(Enemax-Enemin), Enemin, 5*Enemax); 
	CosZ0vsndf = new TH2F("CosZ0 vs ndf","CosZ0 vs ndf", 10,-10,10,30,0,30);
	CosZ0vsEkreco = new TH2F("CosZ0 vs Ekreco","CosZ0 vs Ekreco", (Enemax-Enemin),Enemin,Enemax,30,0,30);
    int npoints = 0;
 
    for(int j=0;j<7;j++)
      {
       LayerPull[j]= new TH1F*[3];
       LayerPullPR[j]=new TH1F*[3];
       Occupancy[j] = new TH1F(Form("Occupancy L%d,%d %s",j,Ene,UNIT.c_str()),Form("Occupancy L%d,%d %s",j,Ene,UNIT.c_str()),500,-10,10); 
		  for(int k=0;k<3;k++) 
		 {
		 LayerPull[j][k] = new TH1F(Form("Pulls L%d,%d %s",j,Ene,UNIT.c_str()),Form("Pulls L%d,%d %s",j,Ene,UNIT.c_str()),10000,-10,10); 
		 LayerPullPR[j][k] = new TH1F(Form("Pulls PR L%d,%d %s",j,Ene,UNIT.c_str()),Form("Pulls PR L%d,%d %s",j,Ene,UNIT.c_str()),10000,-10,10); 
     		 }
	  }
    
    //Define variables to read event
    ALEvent *e = new ALEvent();      
    //Set address to access event data
   // chain->SetBranchAddress("RecoEvent",&e); 
	chain->SetBranchAddress("Revent",&e); 
  
    // Get number of event in Tree
    int nentries=chain->GetEntries();
    cout << "Number  of events: " << nentries << endl;  
      
    for (int j=0;j<nentries;j++)
      {
       chain->GetEntry(j); //Load the entry i in the variable e 
       if(j==0)Enek=1000*e->get_EkMC();
       bool* t=new bool[NTConfig];
       for(int k=0;k<NTConfig;k++) t[k]=true;
       if(j%1000000==0) cout << "Event: " << j <<endl;
    
       //No threshold is applied on the signal in the scintillators
       //Should be done here in the future  
       if(!e->get_T1())t[1]=t[4]=t[5]=t[6]=t[7]=t[8]=false;
       if(!e->get_T3())t[2]=t[4]=t[6]=t[7]=t[8]=false;
       if(!e->get_T4())t[3]=t[5]=t[6]=t[7]=t[8]=false;
	   //Check on the resonctruction Ek between  0 and 2 times the truth Ek
       if(e->get_Ekreco()<0){
      // cout << "Ek < 0 " << endl;
       continue;
       } 
       if(e->get_Ekreco()>5*1000*e->get_EkMC()){
       //cout << "Ekreco 5 times higher than EkMC" << endl;
       continue;
       }
     //   cout << "we're in!" << endl;
	 
      //Fill total energy deposited in the scintillators
	  //Internal trigger requirement: at least 5 hits, event chosen by PR

		//apply external trigger requirements
	bool T1=e->get_T1();
	bool T3=e->get_T3();
	bool T4=e->get_T4();

	if(T1 && T4) {
		//cout <<"internal trig "  << endl;
	 if(e->get_Ekreco()!=-999) {				//check that event was reconstructed	
	
           EkReso->Fill(e->get_Ekreco());
	   inversep0reco->Fill(1/(e->get_p0reco())); 
		 EkPR->Fill(e->get_EkPR()*1000);
	   chi2vsEkreco->Fill((e->get_chi2())/(e->get_ndf()),e->get_Ekreco());
           ndfvsEkreco->Fill(e->get_ndf(),e->get_Ekreco());
          
           //helix parameters pull distribution
           TMatrixD CovLast = e->get_Cov_last();
           double phi0init = e->get_phi0_init();
           double tanLinit = e->get_tanl_init();
           double cpainit = e->get_cpa_init();
           double phi0last = e->get_phi0();
           double tanLlast = e->get_tanl();
           double cpalast = e->get_cpa();
           double phierr2 = e->get_phi0err2();
           double tanlerr2 = e->get_tanlerr2();
           double cpaerr2 = e->get_cpaerr2();
           int ndf = e->get_ndf();
           double chi2=e->get_chi2();
           double chi2_ndf=chi2/ndf;
           double eMC, pMC;
	       double CosZ0 = e->get_CZ0MC();
	   double Theta0 = 180-((180*TMath::ACos(CosZ0))/(TMath::Pi()));
	  // cout << "CosZ0 = " << CosZ0 << endl;
	  // cout << "Theta0 = "<<Theta0 << endl;
	   	   CosZ0vsndf->Fill(ndf,Theta0);
	       CosZ0vsEkreco->Fill(e->get_Ekreco(),Theta0);
	     //make cuts as function of incident angle
           if((Theta0 > thetamin) && (Theta0 < thetamax)) {
	         Phi0Pull->Fill((phi0last-phi0init)/TMath::Sqrt(phierr2));
      	     TanLPull->Fill((tanLlast-tanLinit)/TMath::Sqrt(tanlerr2));
             KappaPull->Fill((cpalast-cpainit)/TMath::Sqrt(cpaerr2));
	    chi2ndf->Fill(e->get_cl());
	       	for(int k=0;k<e->get_Nhits();k++)
        	 {
	    	int Lindex=(int)e->get_hits().at(k)->get_L();
            if(e->get_hits().at(k)->get_flagPR() && e->get_hits().at(k)->get_fUsed() &&(!e->get_hits().at(k)->get_fGhost()))
             {
				float tmp=0;
				tmp=(e->get_hits().at(k)->get_x())*10;
				LayerPull[Lindex][0]->Fill(tmp-e->get_hits().at(k)->get_xreco());
				tmp=(e->get_hits().at(k)->get_y())*10;				
				LayerPull[Lindex][1]->Fill(tmp-e->get_hits().at(k)->get_yreco());
				tmp=(e->get_hits().at(k)->get_z())*10;
				LayerPull[Lindex][2]->Fill(tmp-e->get_hits().at(k)->get_zreco());

			  //Pulls truth - PR
				tmp=(e->get_hits().at(k)->get_xPR())*10;
				LayerPullPR[Lindex][0]->Fill((e->get_hits().at(k)->get_x()*10) -tmp);
				tmp=(e->get_hits().at(k)->get_yPR())*10;				
				LayerPullPR[Lindex][1]->Fill((e->get_hits().at(k)->get_y()*10) - tmp);
				tmp=(e->get_hits().at(k)->get_zPR())*10;
				LayerPullPR[Lindex][2]->Fill((e->get_hits().at(k)->get_z()*10) -tmp);
				if(Lindex==0 || Lindex==4 || Lindex==6) Occupancy[Lindex]->Fill(e->get_hits().at(k)->get_x());
				else Occupancy[Lindex]->Fill(e->get_hits().at(k)->get_y());
				if(Lindex==0) {
					eMC = 1000*(e->get_hits().at(k)->get_eMC());
					eMCL0->Fill(eMC);
					pMC=TMath::Sqrt((eMC*eMC)-(RestMass*RestMass));
					} //end if L0			   	
				} //if hit used
			}	//k
	
		   EkResoCut->Fill(e->get_Ekreco());
		   ErecovsEL0->Fill(eMC, e->get_Ekreco());
		  if(eMC>eMC_cut*Ene) {
		  EkrecoCutL0->Fill(e->get_Ekreco());
		  inversep0CutL0->Fill(1/pMC);
		  inversep0cut->Fill(1/(e->get_p0reco()));
		  inversep0PR->Fill(1/(e->get_p0PR()*1000));
	     
		  	} //end if condition on eMC
		  } //end condition on incident angle
		} //check event reco
		}  //end external/internal trigger requirement 
      }//j
  
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
   
  double sigma=3;
  
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1111);
  TCanvas*CanEKReso=new TCanvas("CanEKReso","CanEKReso",200,10,1200,800);
  CanEKReso->cd(1);
  double max =0;
  double maxEkReso = EkReso->GetMaximum();
  double  maxEkPR = EkPR->GetMaximum();
	if(maxEkReso>maxEkPR) max=maxEkReso;
	else max=maxEkPR;
  EkReso->SetMaximum(1.10*max);
  EkReso->SetTitle(Form("%s: %d %s",TConfig[0].c_str(),Ene,UNIT.c_str()));
  EkReso->SetLineColor(kBlue);
  EkPR->SetLineColor(kRed);
  EkReso->GetXaxis()->SetTitle("Reconstruction E_{k} in MeV");
  EkReso->GetYaxis()->SetTitle("#entries");
  EkReso->GetXaxis()->SetRangeUser(EkReso->GetMean()-sigma*EkReso->GetRMS(),EkReso->GetMean()+sigma*EkReso->GetRMS());
// EkReso->Fit("landau");	
  EkReso->Draw("hist");
  EkPR->Draw("same");
  fileout->cd();
  EkReso->Write();
  CanEKReso->Write();
  CanEKReso->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf(", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

 
  //////////////////////////////////////////////////////////////////
	  TCanvas*CanInversep0reco=new TCanvas("CanInversep0reco","CanInversep0reco",200,10,1200,800);
     TLine *p0MCline=new TLine;

     CanInversep0reco->cd(1);
     inversep0reco->Draw("hist");
     inversep0reco->SetTitle(Form("%s, all events: Inverse momentum %d MeV",TConfig[0].c_str(),Ene));
     inversep0reco->SetLineColor(kBlue);
	 inversep0PR->SetLineColor(kRed);
     inversep0reco->GetXaxis()->SetTitle("Reconstruction 1/p_{0} in MeV^{-1}");
     inversep0reco->GetXaxis()->SetRangeUser(inversep0reco->GetMean()-sigma*inversep0reco->GetRMS(),inversep0reco->GetMean()+sigma*inversep0reco->GetRMS());
     inversep0reco->GetYaxis()->SetTitle("#entries");
     inversep0reco->Draw("same");
     fileout->cd();
     inversep0reco->Write();
     CanInversep0reco->Write();
     CanInversep0reco->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

  
 //////////////////////////////////////////////////////////////////
     TCanvas*CanInversep0recocut=new TCanvas("CanInversep0recocut","CanInversep0recocut",200,10,1200,800);
     TLegend*leg1=new TLegend(0.1,0.8,0.3,0.7,"NDC");
     TPaveText *PTp0cut = new TPaveText(0.15,0.8,0.25,0.9,"NDC");
     PTp0cut->SetBorderSize(0);
     PTp0cut->SetFillStyle(0);
     leg1->SetBorderSize(0);leg1->SetFillStyle(0);
     leg1->AddEntry(inversep0cut,"1/preco for #theta_{0} & eMC cut","l"); 
     leg1->AddEntry(inversep0PR,"1/p0PR for #theta_{0} & eMC cut","l");
     CanInversep0recocut->cd(1);
	  double maxP0 =0;
   double maxP0Reco = inversep0cut->GetMaximum();
  double  maxp0PR = inversep0PR->GetMaximum();
   if(maxP0Reco>maxp0PR) maxP0=maxP0Reco;
     else maxP0=maxp0PR;
     inversep0cut->SetMaximum(1.10*maxP0);
     inversep0cut->SetTitle(Form("%s: Inverse momentum %d MeV for %d#circ< #theta_{0} < %d#circ",TConfig[0].c_str(), Ene,thetamin,thetamax)); 
     inversep0cut->GetXaxis()->SetTitle("Reconstruction 1/p_{0} in MeV^{-1}");
     inversep0cut->GetYaxis()->SetTitle("#entries");
     inversep0cut->SetStats(1);
     inversep0cut->SetLineColor(kBlue);
     inversep0PR->SetLineColor(kRed);
     double meanp0 = inversep0cut->GetMean();
     double  sigmap0 = inversep0cut->GetRMS();
     TF1 *gaussianp0 = new TF1("gaussianp0", "gaus", meanp0-sigma*sigmap0, meanp0+sigma*sigmap0) ;
     inversep0cut->Fit(gaussianp0,"R0");	
     inversep0cut->Draw("hist");
     double meanPR = inversep0PR->GetMean();
     double  sigmaPR = inversep0PR->GetRMS();
     TF1 *gaussianPR = new TF1("gaussianPR", "gaus", meanPR-sigma*sigmaPR, meanPR+sigma*sigmaPR) ;
     inversep0PR->Fit(gaussianPR,"R0"); 
     inversep0PR->Draw("sames"); 
     inversep0cut->GetXaxis()->SetRangeUser(inversep0cut->GetMean()-sigma*inversep0cut->GetRMS(),inversep0cut->GetMean()+sigma*inversep0cut->GetRMS());
     leg1->Draw("same");
     PTp0cut->AddText(Form(" KF fit #sigma = %4.1f%% ", (100*((gaussianp0->GetParameter(2))/(gaussianp0->GetParameter(1))))));
     PTp0cut->AddText(Form(" PR fit #sigma = %4.1f%% ", (100*((gaussianPR->GetParameter(2))/(gaussianPR->GetParameter(1))))));
     gPad->Update();
     TPaveStats *sKF = (TPaveStats*)inversep0cut->FindObject("stats");
     double y1 = sKF->GetY1NDC();
     double y2 = sKF->GetY2NDC();
     gPad->Update();
     TPaveStats *sPR = (TPaveStats*)inversep0PR->FindObject("stats");
     gPad->Update();
     sKF->SetTextColor(kBlue);
     sPR->SetTextColor(kRed);
     sPR->SetY2NDC(y1);
     sPR->SetY1NDC(y1-fabs(y2-y1));
     gPad->Update();

     PTp0cut->Draw();
     gPad->Update();
     fileout->cd();
    inversep0cut->Write();
    inversep0CutL0->Write();
     CanInversep0recocut->Write();
    CanInversep0recocut->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

  
  
  ////////////////////////////////////////////////////////////////////
  int nsig=3;
  TCanvas *CanHelixPull = new TCanvas("CanHelixPull", "CanHelixPull",20,10,1200,800);
  CanHelixPull->Divide(3,1);
  CanHelixPull->cd(1);
  Phi0Pull->SetTitle(Form("Helix parameters phi0 errors"));
  Phi0Pull->SetLineColor(kBlack);
  Phi0Pull->GetXaxis()->SetTitle("(#phi_{0fit}-#phi_{0PR})/#sigma_{#phi0}");
  Phi0Pull->GetYaxis()->SetTitle("#entries");
  double meanphi0 = Phi0Pull->GetMean();
  double  sigmaphi0 = Phi0Pull->GetRMS();
  Phi0Pull->GetXaxis()->SetRangeUser(meanphi0-nsig*sigmaphi0,meanphi0+nsig*sigmaphi0);
   Phi0Pull->SetMaximum(Phi0Pull->GetMaximum()*1.10);
  TF1 *gaussianphi0 = new TF1("gaussianphi0", "gaus", meanphi0-sigmaphi0, meanphi0+sigmaphi0) ;
   // Phi0Pull->Fit(gaussianphi0,"R");	
  Phi0Pull->Draw("same");
  
  CanHelixPull->cd(2);
  TanLPull->SetTitle(Form("Helix parameters tanL errors"));
  TanLPull->SetLineColor(kBlack);
  TanLPull->GetXaxis()->SetTitle("(TanL_{0fit}-TanL_{0PR})/#sigma_{TanL}");
  TanLPull->GetYaxis()->SetTitle("#entries");
  double meanTanL = TanLPull->GetMean();
  double  sigmaTanL = TanLPull->GetRMS();
  TanLPull->GetXaxis()->SetRangeUser(meanTanL-nsig*sigmaTanL, meanTanL+nsig*sigmaTanL);
  TanLPull->SetMaximum(TanLPull->GetMaximum()*1.10);
  TF1 *gaussianTanL = new TF1("gaussianTanL", "gaus", meanTanL-sigmaTanL, meanTanL+sigmaTanL);
 // gaussianTanL->SetParameter(5,6.6);
 // gaussianTanL->SetParameter(2,0.01);
 //   TanLPull->Fit(gaussianTanL,"R");	
  TanLPull->Draw("same");
  CanHelixPull->cd(3);
  KappaPull->SetTitle(Form("Helix parameters kappa errors"));
  KappaPull->SetLineColor(kBlack);
  KappaPull->GetXaxis()->SetTitle("(#kappa_{fit}-#kappa_{PR})/#sigma_{#kappa}");
  KappaPull->GetYaxis()->SetTitle("#entries");
  double meankappa = KappaPull->GetMean();
  double  sigmakappa = KappaPull->GetRMS();
  KappaPull->GetXaxis()->SetRangeUser(meankappa-nsig*sigmakappa,meankappa+nsig+sigmakappa);
  KappaPull->SetMaximum(KappaPull->GetMaximum()*1.10); 
  TF1 *gaussiankappa = new TF1("gaussiankappa", "gaus", meankappa-nsig*sigmakappa, meankappa+nsig*sigmakappa) ;
    KappaPull->Fit(gaussiankappa,"R0");	
  KappaPull->Draw("same");
  fileout->cd();
  Phi0Pull->Write();
  TanLPull->Write();
  KappaPull->Write();

   CanHelixPull->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

////////////////////////////////////////////////////////////////////////////////////

TCanvas *CanOccupancy = new TCanvas("Occupancy per layer", "Occupancy per layer", 200,10,1200,800);

CanOccupancy->Divide(4,2);

for(int j=0;j<7;j++) {
CanOccupancy->cd(j+1);
Occupancy[j]->Draw("hist");
if(j==0||j==4||j==6) {
 Occupancy[j]->SetLineColor(kBlue);
 Occupancy[j]->GetXaxis()->SetTitle("X_{used} in cm");
 }
else {
 Occupancy[j]->SetLineColor(kRed);
 Occupancy[j]->GetXaxis()->SetTitle("Y_{used} in cm");
 }
Occupancy[j]->SetTitle(Form("Layer %d, coordinate used",j));
fileout->cd();
Occupancy[j]->Write();
}
     CanOccupancy->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));


TCanvas *Canchi2vsEkreco = new TCanvas("cabchi2vsekreco","canchi2vsekreco",200,10,1200,800);
Canchi2vsEkreco->cd(1);
chi2vsEkreco->Draw("colz");
float tmp2=0;
tmp2=chi2vsEkreco->GetBinContent(chi2vsEkreco->GetMaximumBin());
chi2vsEkreco->SetMaximum(1.10*tmp2);
chi2vsEkreco->SetTitle("Ekreco vs chi2 all events");
chi2vsEkreco->GetYaxis()->SetTitle("Ekreco in MeV");
chi2vsEkreco->GetXaxis()->SetTitle("chi2/ndf");
if(type==3||type==4) chi2vsEkreco->GetYaxis()->SetRangeUser(0,2.0*Ene);
else chi2vsEkreco->GetYaxis()->SetRangeUser(0,2000*Ene);
fileout->cd();
chi2vsEkreco->Write();
Canchi2vsEkreco->Write();
Canchi2vsEkreco->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

  //////////////////////////////////////////////////////////////////
 
 TCanvas *CanCosZ0vsEkreco=new TCanvas("CanCosZ0vsEkreco","CanCosZ0vsEkreco",200,10,1200,800);

CanCosZ0vsEkreco->cd(1);
CosZ0vsEkreco->Draw("colz");
CosZ0vsEkreco->SetTitle("Incident angle vs Ekreco");
CosZ0vsEkreco->GetYaxis()->SetTitle("Incident normal angle (deg)");
CosZ0vsEkreco->GetXaxis()->SetTitle("Ekreco in MeV");
if(type==3||type==4) CosZ0vsEkreco->GetXaxis()->SetRangeUser(0,3.0*Ene);
else CosZ0vsEkreco->GetXaxis()->SetRangeUser(0,3000*Ene);
fileout->cd();
CosZ0vsEkreco->Write();
CanCosZ0vsEkreco->Write();
CanCosZ0vsEkreco->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

  //////////////////////////////////////////////////////////////////
 
 TCanvas *CanCosZ0vsndf=new TCanvas("CanCosZ0vsEkndf","CanCosZ0vsEkndf",200,10,1200,800);

CanCosZ0vsndf->cd(1);
CosZ0vsndf->Draw("colz");
CosZ0vsndf->SetTitle("Incident angle vs ndf");
CosZ0vsndf->GetYaxis()->SetTitle("Incident normal angle (deg)");
CosZ0vsndf->GetXaxis()->SetTitle("Reconstructed ndf");
fileout->cd();
CosZ0vsndf->Write();
CanCosZ0vsndf->Write();
CanCosZ0vsndf->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

///////////////////////////////////////////////////////////////////////////////////////////

TCanvas *CanEkrecovsEL02D=new TCanvas("CanEkrecovsEL0","CanEkrecovsEL0",200,10,800,800);

TF1 *straightline = new TF1("line","x",0,2*Ene);
CanEkrecovsEL02D->cd(1);
ErecovsEL0->SetStats(0);
ErecovsEL0->Draw("colz");
tmp2=ErecovsEL0->GetBinContent(ErecovsEL0->GetMaximumBin());
//ErecovsEL0->SetMaximum(1.10*tmp2);
ErecovsEL0->SetTitle(Form("Ekreco vs eMC at L0, %d %s at injection point, incident angle cut",Ene,UNIT.c_str()));
ErecovsEL0->GetYaxis()->SetTitle("Ekreco in MeV");
ErecovsEL0->GetXaxis()->SetTitle("eMC at L0 in MeV");
if(type==3||type==4) {
	ErecovsEL0->GetYaxis()->SetRangeUser(0,2*Ene);
	ErecovsEL0->GetXaxis()->SetRangeUser(0,2*Ene);
	}
else {
	ErecovsEL0->GetYaxis()->SetRangeUser(0,2000*Ene);
	ErecovsEL0->GetXaxis()->SetRangeUser(0,2000*Ene);
	}
straightline->Draw("same");
fileout->cd();
ErecovsEL0->Write();
CanEkrecovsEL02D->Write();
CanEkrecovsEL02D->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

///////////////////////////////////////////////////////////////////////////////////////////


TCanvas *CanndfvsEkreco=new TCanvas("CanndfvsEkreco","CanndfvsEkreco",200,10,1200,800);

CanndfvsEkreco->cd(1);
gStyle->SetOptStat(0);
ndfvsEkreco->Draw("colz");
tmp2=ndfvsEkreco->GetBinContent(ndfvsEkreco->GetMaximumBin());
ndfvsEkreco->SetMaximum(1.10*tmp2);
ndfvsEkreco->SetTitle(Form(" %s, all events: Ekreco %d MeV vs ndf",TConfig[0].c_str(), Ene));
ndfvsEkreco->GetYaxis()->SetTitle("Ekreco in MeV");
ndfvsEkreco->GetXaxis()->SetTitle("ndf");
if(type==3||type==4) ndfvsEkreco->GetYaxis()->SetRangeUser(0,2*Ene);
else ndfvsEkreco->GetYaxis()->SetRangeUser(0,2000*Ene);
fileout->cd();
ndfvsEkreco->Write();
CanndfvsEkreco->Write();
CanndfvsEkreco->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

  //////////////////////////////////////////////////////////////////
  TCanvas*Canchi2ndf=new TCanvas("Canchi2ndf","Canchi2ndf",200,10,1200,800);
	
 
 
     Canchi2ndf->cd(1);
     chi2ndf->Draw("hist");
     gStyle->SetOptStat(1111);
     chi2ndf->SetTitle(Form("%s: %d %s",TConfig[0].c_str(), Ene, UNIT.c_str()));
     chi2ndf->SetLineColor(kBlack);
     chi2ndf->GetXaxis()->SetTitle(Form("Confidence level for %d<#theta_{0}<%d",thetamin,thetamax));
     chi2ndf->GetYaxis()->SetTitle("#entries");
     fileout->cd();
     chi2ndf->Write();
     Canchi2ndf->Write();
     Canchi2ndf->Print(Form("%s/%d/RecoTest_%s_%d_%d%s_%s_%dto%d.pdf)", directory.c_str(),type,startfile.c_str(), type, Ene,UNIT.c_str(),ID.c_str(),thetamin,thetamax));

  cout << "Last canvas printed" << endl;


//close ROOT file
     fileout->Close();
     cout << "Close fileout " << endl;
 
}


