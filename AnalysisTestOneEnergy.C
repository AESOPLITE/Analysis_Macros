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

void AnalyseTestOneEnergyOld();

void AnalyseTestOneEnergyOld()
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
  double B = 0.3;		//average magnetic field in T
  double c = TMath::C();
 string MCparamfile="./MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh);
 	


 int NTConfig=1;//Number of configuration of triggers
 string TConfig[1]={"T1&T3&T4&NoGuard"};
// double TrigThresh[4]={0.3, 0.29, 0.57,0.1};		//trigger thresholds in MeV determined for T1/T3/T4/Guard
 double RestMass   = 0.51099;				//mass electron in MeV

 //Input file  
 string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V2/Disc"; 
 //filename structure
 string startfile="aesopliteNonUniB_V2";
 string endfile="_fort.99";
 //string RecoInd="7hits";
 string ID = "OneIterForwardPRinitNonUniB"; 
 string directory = "OldMacro";
	
 //type of particle
 int type=3;
 //Number of energies
 //int Nene=12;
 int Nene=1;
 //Energies
 int EneMuon = 1;
 int Ene = 40;
 double Enek;
 //Number of cycles per energy
 int Ncycles=10;  
//cut values of renconstructed degrees of freedom 
int ndf_cut = 4;
	
//Define histogram outputs

 TH1F *EkReso=new TH1F;
 TH1F *EkResoCut = new TH1F;
 TH1F *EkrecoCutL0 = new TH1F;
 TH1F *eMCL0 = new TH1F;
 TH1F *inversep0reco=new TH1F;
 TH1F *inversep0cut=new TH1F;
 TH1F *inversep0CutL0=new TH1F;
 TH1F *chi2ndf=new TH1F;
 TH1F ***LayerPull=new TH1F**;

//Define scatter plot output
 TH2F *ErecovsEL0 = new TH2F;
 TH2F *chi2vsEkreco = new TH2F;
 TH2F *ndfvsEkreco = new TH2F;
 TH2F  *EkRecovsDeltaPhi0 = new TH2F;
 TH2F  *EkRecovsDeltaCpa = new TH2F; 
 TH2F  *EkRecovsDeltaTanL = new TH2F;
 //Create TChain to merge the files for a given energy
 TChain *chain=new TChain; 

//ROOT output file 
TFile*fileout=new TFile(Form("%s/RecoTest_%s_%d_%dMeV_%s.root", directory.c_str(),startfile.c_str(), type, Ene,ID.c_str()),"RECREATE");

 
   
    cout << "Energy: " << Ene << "MeV" <<endl;
    chain=new TChain("MC");
    
   for(int j=0;j<Ncycles;j++)//Number of cycles
      
      { 
		  
     // chain->Add(Form("%s/%d/RecoEvent_%s_%d_%dMeV%03d%s_%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j+1,endfile.c_str(),RecoInd.c_str(),ID.c_str()));
  //chain->Add(Form("%s/%d/RecoEvent_%s_%d_%dGeV%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,EneMuon,j+1,endfile.c_str(),ID.c_str()));
     chain->Add(Form("%s/%d/RecoEvent_%s_%d_%dMeV%03d%s_%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j+1,endfile.c_str(),ID.c_str()));

	  } 
    
   
    //Define histograms 
    EkReso=new TH1F(Form("Ekreco %d MeV",Ene),Form("Ekreco %d MeV",Ene),2*Ene,0,5*Ene);
	EkrecoCutL0 = new TH1F(Form("Ekreco ndf&eMC cut %d MeV",Ene),Form("Ekreco ndf&eMC cut %d MeV",Ene),2*Ene,0,5*Ene);
	EkResoCut=new TH1F(Form("Ekreco chi2 cut %d MeV",Ene),Form("Ekreco chi2 cut %d MeV",Ene),2*Ene,0,5*Ene);
    eMCL0 = new TH1F("eMC at L0", "eMC at L0", 2*Ene, 0, Ene);
	inversep0reco=new TH1F(Form("1/p0reco %d MeV",Ene),Form("1/p0reco %d MeV",Ene),2*Ene,0, 2.0/Ene);
	inversep0CutL0 = new TH1F(Form("1/p0 ndf&eMC cut %d MeV",Ene),Form("1/p0 ndf&eMC %d MeV",Ene),2*Ene,0, 2.0/Ene);
	inversep0cut=new TH1F(Form("1/p0reco %d MeV cut",Ene),Form("1/p0reco %d MeV cut",Ene),2*Ene,0, 2.0/Ene);
    chi2ndf=new TH1F(Form("chi2ndf %d MeV",Ene),Form("chi2ndf %d MeV", Ene),20,-10,10);
    chi2vsEkreco = new TH2F(Form("Ekreco vs chi2/ndf %d MeV", Ene),Form("Ekreco vs chi2/ndf %d MeV", Ene), 100, -15, 15, 2*Ene, 0, 5*Ene); 
    ErecovsEL0  = new TH2F("Ekreco vs eMC at L0", "Ekreco vs eMC at L0", 2*Ene, 0, Ene, Ene, 0, 2*Ene);
	EkRecovsDeltaPhi0 = new TH2F("Ekreco vs #Delta_{#phi_{0}}", "Ekreco vs #Delta_{#phi_{0}}", 100, -1, 1, 2*Ene, 0, 5*Ene);
	EkRecovsDeltaCpa = new TH2F("Ekreco vs #Delta_{#kappa}}", "Ekreco vs #Delta_{#kappa}", 1000, -5, 5, 2*Ene, 0, 5*Ene);
	EkRecovsDeltaTanL = new TH2F("Ekreco vs #Delta_{TanL}}", "Ekreco vs #Delta_{TanL}", 100, -1, 1, 2*Ene, 0, 5*Ene);
	ndfvsEkreco =  new TH2F(Form("Ekreco vs ndf %d MeV", Ene),Form("Ekreco vs ndf %d MeV", Ene), 60, -10, 50, 2*Ene, 0, 5*Ene); 
    int npoints = 0;
 
    for(int j=0;j<7;j++)
      {
       LayerPull[j]= new TH1F*[3];
		  for(int k=0;k<3;k++) 
		 {
		 LayerPull[j][k] = new TH1F(Form("Pulls L%d,%d MeV",j+1,Ene),Form("Pulls L%d,%d MeV",j+1,Ene),200,-2,2); 
     		 }
	  }
    
    //Define variables to read event
    ALEvent *e = new ALEvent();      
    //Set address to access event data
    chain->SetBranchAddress("Revent",&e); 
	//chain->SetBranchAddress("Revent",&e); 
  
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
       
//apply external trigger requirements
	   int nT1=e->get_EneT1().size();
	   double tmpT1=0;
	   for(int l=0;l<nT1;l++) tmpT1+=1000*e->get_EneT1().at(l);
	   int nT3=e->get_EneT3().size();
	   double tmpT3=0;
	   for(int l=0;l<nT3;l++) tmpT3+=1000*e->get_EneT3().at(l);
		int nT4=e->get_EneT4().size();
	   double tmpT4=0;
	   for(int l=0;l<nT4;l++) tmpT4+=1000*e->get_EneT4().at(l);
	   int nTG=e->get_Eneg().size();
	   double tmpG=0;
	   for(int l=0;l<nTG;l++) tmpG+=1000*e->get_Eneg().at(l);
		 int NL= e->get_NLayers();  
	//apply T1&T3&T4&NoGuard&7hits requirement
       if( tmpT1>TrigThresh[0] && tmpT3>TrigThresh[1] && tmpT4>TrigThresh[2]&& tmpG<TrigThresh[3] && NL==7) 
	   {
		if(e->get_Ekreco()!=999) {				//check that event was reconstructed	
       EkReso->Fill(e->get_Ekreco());
	   inversep0reco->Fill(1/(e->get_p0reco())); 
	   chi2vsEkreco->Fill((e->get_chi2())/(e->get_ndf()),e->get_Ekreco());
       ndfvsEkreco->Fill(e->get_ndf(),e->get_Ekreco());
       EkRecovsDeltaPhi0->Fill(e->get_phi0(),e->get_Ekreco());
	  EkRecovsDeltaCpa->Fill(e->get_cpa(),e->get_Ekreco());
       EkRecovsDeltaTanL->Fill(e->get_tanl(),e->get_Ekreco());
	
       if(e->get_ndf()!=0) {
		   chi2ndf->Fill(e->get_chi2()/e->get_ndf());
	   }
	
double eMC, pMC;
	         for(int k=0;k<e->get_Nhits();k++)
         {
          if(e->get_hits().at(k)->get_xreco()!=0 && e->get_hits().at(k)->get_yreco()!=0)//Check that the hit was used for reconstruction: Will be change to -999
           {
	    float tmp=0;
	    int Lindex=(int)e->get_hits().at(k)->get_mregMC()%11;
	    tmp=(e->get_hits().at(k)->get_xin()+e->get_hits().at(k)->get_xout())*10/2;
	    LayerPull[Lindex][0]->Fill(tmp-e->get_hits().at(k)->get_xreco());
	    tmp=(e->get_hits().at(k)->get_yin()+e->get_hits().at(k)->get_yout())*10/2;
	    LayerPull[Lindex][1]->Fill(tmp-e->get_hits().at(k)->get_yreco());
	    tmp=(e->get_hits().at(k)->get_zin()+e->get_hits().at(k)->get_zout())*10/2;
	    LayerPull[Lindex][2]->Fill(tmp-e->get_hits().at(k)->get_zreco());
			   
	  		if(Lindex==0) {
				eMC = 1000*(e->get_hits().at(k)->get_eMC());
				eMCL0->Fill(eMC);
				} //end if L0	
		   	}//end if 	
	 	}	//k
	  //make cuts for given ndf
	   if(e->get_ndf()>=ndf_cut) {
		   EkResoCut->Fill(e->get_Ekreco());
		   inversep0cut->Fill(1/(e->get_p0reco()));
		   ErecovsEL0->Fill(eMC, e->get_Ekreco());

		  if(eMC>0.9*Ene) {
		  EkrecoCutL0->Fill(e->get_Ekreco());
		  inversep0CutL0->Fill(1/(e->get_p0reco()));
		  	} //end if condition on eMC
		  } //end ndf condition
		} //check event reco
		}  //end external/internal trigger requirement 
      }//j
  
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
   
  
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(11111111);
  TCanvas*CanEKReso=new TCanvas("CanEKReso","CanEKReso",200,10,1200,800);
  TLine *EkMCline=new TLine;
  //TF1 *gaussian = new TF1("gaussian", "gaus",0, 200) ;

     CanEKReso->cd(1);
     EkReso->Draw("hist");
     EkReso->SetTitle(Form("%s: %d MeV",TConfig[0].c_str(),Ene));
     EkReso->SetLineColor(kBlack);
     EkReso->GetXaxis()->SetTitle("Reconstruction E_{k} in MeV");
     EkReso->GetYaxis()->SetTitle("#entries");
    // EkReso->Fit("landau");	
	 EkReso->Draw("same");
     EkMCline=new TLine(Enek,0,Enek,1.01*EkReso->GetBinContent(EkReso->GetMaximumBin()));
     EkMCline->SetLineColor(kBlue);
     EkMCline->Draw("same");
     fileout->cd();
     EkReso->Write();
     CanEKReso->Write();
  CanEKReso->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf(",directory.c_str(), startfile.c_str(), type, Ene,ID.c_str()));

  
  
  //////////////////////////////////////////////////////////////////
	
  TCanvas*CanEKResoCut=new TCanvas("CanEKResoCut","CanEKResoCut",200,10,1200,800);
  TLegend*leg=new TLegend(0.45,0.6,0.85,0.9);
  leg->SetBorderSize(0);leg->SetFillStyle(0);
  leg->AddEntry(EkResoCut,Form("Ekreco for ndf > %d", ndf_cut),"l");
  leg->AddEntry(EkrecoCutL0,Form("Ekreco for ndf > %d && eMC at L0>0.9*EkMC", ndf_cut),"l");

     CanEKResoCut->cd(1);
     EkResoCut->Draw("hist");
	// gPad->Update();
	//TPaveStats *stReco = (TPaveStats*)EkResoCut->FindObject("stats");
	//double y1=stReco->GetY1NDC(); //new x start position
	//double y2 =stReco->GetY2NDC(); //new x end position	
	
     EkResoCut->SetStats(0);
     EkResoCut->SetTitle(Form("%s: %d MeV, ndf >= %d",TConfig[0].c_str(),Ene,ndf_cut));
     EkResoCut->SetLineColor(kBlack);
     EkResoCut->GetXaxis()->SetTitle("Reconstruction E_{k} in MeV");
     EkResoCut->GetYaxis()->SetTitle("#entries");

	 EkrecoCutL0->SetLineColor(kBlue);
	 EkrecoCutL0->Draw("sames");
	 eMCL0->SetLineColor(kRed);
	 eMCL0->Draw("sames");
	 eMCL0->SetStats(0);
	 //gPad->Update();
	 //TPaveStats *stcut = (TPaveStats*)EkrecoCutL0->FindObject("stats");
	 //stcut->SetY1NDC(y2-y1); 
	 //stcut->SetY2NDC(y1); 
	leg->Draw("same");
   fileout->cd();
   EkResoCut->Write();
   EkrecoCutL0->Write();
   CanEKResoCut->Write();
   
  CanEKResoCut->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf",directory.c_str(), startfile.c_str(), type, Ene,ID.c_str()));

  
  
  //////////////////////////////////////////////////////////////////
	  TCanvas*CanInversep0reco=new TCanvas("CanInversep0reco","CanInversep0reco",200,10,1200,800);
     TLine *p0MCline=new TLine;

     CanInversep0reco->cd(1);
     inversep0reco->Draw("hist");
     inversep0reco->SetTitle(Form("%s: Inverse momentum %d MeV",TConfig[0].c_str(),Ene));
     inversep0reco->SetLineColor(kBlack);
     inversep0reco->GetXaxis()->SetTitle("Reconstruction 1/p_{0} in MeV^{-1}");
     inversep0reco->GetYaxis()->SetTitle("#entries");
     inversep0reco->Draw("same");
     fileout->cd();
     inversep0reco->Write();
     CanInversep0reco->Write();
  CanInversep0reco->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf",directory.c_str(), startfile.c_str(), type, Ene,ID.c_str()));

  
 //////////////////////////////////////////////////////////////////
	  TCanvas*CanInversep0recocut=new TCanvas("CanInversep0recocut","CanInversep0recocut",200,10,1200,800);
  TLine *p0MClinecut=new TLine;
    TLegend*leg1=new TLegend(0.05,0.7,0.45,1.0);

  //leg1->SetBorderSize(0);leg->SetFillStyle(0);
  leg1->AddEntry(inversep0cut,Form("1/p0 for ndf > %d", ndf_cut),"l");
  leg1->AddEntry(inversep0CutL0,Form("1/p0 for ndf > %d && eMC at L0>0.9*EkMC", ndf_cut),"l");
     CanInversep0recocut->cd(1);
     inversep0cut->SetTitle(Form("%s: Inverse momentum %d MeV at injection point for ndf >= %d",TConfig[0].c_str(), Ene,ndf_cut));
     inversep0cut->SetLineColor(kBlack);
     inversep0cut->GetXaxis()->SetTitle("Reconstruction 1/p_{0} in MeV^{-1}");
     inversep0cut->GetYaxis()->SetTitle("#entries");

     inversep0cut->SetStats(0);
     inversep0CutL0->SetLineColor(kBlue);

     double meanp0 = inversep0CutL0->GetMean();
     double  sigmap0 = inversep0CutL0->GetRMS();
     TF1 *gaussianp0 = new TF1("gaussianp0", "gaus", meanp0-sigmap0, meanp0+sigmap0) ;
     inversep0CutL0->Fit(gaussianp0,"R");	
     inversep0cut->Draw("same");
     inversep0CutL0->Draw("sames");
     leg1->Draw("same");

     p0MClinecut=new TLine((1/Enek),0,(1/Enek),1.01*EkReso->GetBinContent(EkReso->GetMaximumBin()));
     p0MClinecut->SetLineColor(kGreen);
     p0MClinecut->Draw("same");
     fileout->cd();
    inversep0cut->Write();
    inversep0CutL0->Write();
     CanInversep0recocut->Write();
  CanInversep0recocut->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf",directory.c_str(), startfile.c_str(), type, Ene,ID.c_str()));

  
  //////////////////////////////////////////////////////////////////
  TCanvas **CanXYZReso=new TCanvas*[7];
	/*
  TLegend*leg=new TLegend(0.15,0.6,0.45,0.9);
  leg->SetBorderSize(0);leg->SetFillStyle(0);
  leg->AddEntry(XReso[0],"#Delta X","l");
  leg->AddEntry(YReso[0],"#Delta Y","l");
  leg->AddEntry(ZReso[0],"#Delta Z","l");
 */
 for(int i=0;i<7;i++) 
 {
    
     CanXYZReso[i]=new TCanvas(Form("XYZReso Layer %d, %dMeV", i, Ene),Form("XYZReso Layer %d, %dMeV",i,Ene),200,10,1200,800);     
     CanXYZReso[i]->Divide(3,1);
	 for(int j=0;j<3;j++) 
	 {
      CanXYZReso[i]->cd(j+1);
   if(j==0) {
   LayerPull[i][j]->SetLineColor(kRed);
	LayerPull[i][j]->SetTitle(Form("Pull X ndf >= %d, %d MeV, Tracker layer %d",ndf_cut,Ene,i+1));
	LayerPull[i][j]->GetXaxis()->SetTitle("truth-reco (mm)");
   }
		    if(j==1) {
   LayerPull[i][j]->SetLineColor(kBlue);
	LayerPull[i][j]->SetTitle(Form("Pull Y ndf >= %d, %d MeV, Tracker layer %d",ndf_cut,Ene,i+1));
	LayerPull[i][j]->GetXaxis()->SetTitle("truth-reco (mm)");
   }
		    if(j==2) {
   LayerPull[i][j]->SetLineColor(kBlack);
	LayerPull[i][j]->SetTitle(Form("Pull Z ndf >= %d, %d MeV, Tracker layer %d",ndf_cut,Ene,i+1));
	LayerPull[i][j]->GetXaxis()->SetTitle("truth-reco (mm)");
   }
        
	float tmp2=0;
	tmp2=LayerPull[i][j]->GetBinContent(LayerPull[i][j]->GetMaximumBin());
	LayerPull[i][j]->SetMaximum(1.10*tmp2);
      LayerPull[i][j]->Draw("hist");
	LayerPull[i][j]->GetXaxis()->SetRangeUser(LayerPull[i][j]->GetMean()-3*LayerPull[i][j]->GetStdDev(), LayerPull[i][j]->GetMean()+3*LayerPull[i][j]->GetStdDev());
	fileout->cd();
	LayerPull[i][j]->Write();
	//if(j==0)leg->Draw("same");
       }//j   
 
      CanXYZReso[i]->Write();
      CanXYZReso[i]->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf",directory.c_str(), startfile.c_str(), type, Ene,ID.c_str()));
} //i
	
///////////////////////////////////////////////////////////////////////////////////
TCanvas *Canchi2vsEkreco=new TCanvas("Canchi2vsEkreco","Canchi2vsEkreco",200,10,1200,800);

Canchi2vsEkreco->cd(1);
chi2vsEkreco->Draw("colz");
chi2vsEkreco->SetTitle("Ekreco MeV vs chi2");
chi2vsEkreco->GetYaxis()->SetTitle("Ekreco in MeV");
chi2vsEkreco->GetXaxis()->SetTitle("chi2/ndf");
fileout->cd();
chi2vsEkreco->Write();
Canchi2vsEkreco->Write();
Canchi2vsEkreco->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf",directory.c_str(),startfile.c_str(), type, Ene,ID.c_str()));
	
  //////////////////////////////////////////////////////////////////
	
TCanvas *CanEkrecovsEL03D=new TCanvas("CanEkrecovsEL0","CanEkrecovsEL0",200,10,1200,800);

CanEkrecovsEL03D->cd(1);
ErecovsEL0->SetStats(0);
ErecovsEL0->Draw("lego2z");
ErecovsEL0->SetTitle(Form("Ekreco vs eMC at L0, %d MeV at injection point, ndf > %d",Ene, ndf_cut));
ErecovsEL0->GetYaxis()->SetTitle("Ekreco in MeV");
ErecovsEL0->GetXaxis()->SetTitle("eMC at L0 in MeV");
fileout->cd();
ErecovsEL0->Write();
CanEkrecovsEL03D->Write();
CanEkrecovsEL03D->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf", directory.c_str(),startfile.c_str(), type, Ene,ID.c_str()));

  //////////////////////////////////////////////////////////////////
	
TCanvas *CanEkrecovsEL02D=new TCanvas("CanEkrecovsEL0","CanEkrecovsEL0",200,10,1200,800);

CanEkrecovsEL02D->cd(1);
ErecovsEL0->SetStats(0);
ErecovsEL0->Draw("colz");
ErecovsEL0->SetTitle(Form("Ekreco vs eMC at L0, %d MeV at injection point, ndf > %d",Ene, ndf_cut));
ErecovsEL0->GetYaxis()->SetTitle("Ekreco in MeV");
ErecovsEL0->GetXaxis()->SetTitle("eMC at L0 in MeV");
fileout->cd();
ErecovsEL0->Write();
CanEkrecovsEL02D->Write();
CanEkrecovsEL02D->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf", directory.c_str(),startfile.c_str(), type, Ene,ID.c_str()));

	
  //////////////////////////////////////////////////////////////////
TCanvas *CanndfvsEkreco=new TCanvas("CanndfvsEkreco","CanndfvsEkreco",200,10,1200,800);

CanndfvsEkreco->cd(1);
ndfvsEkreco->Draw("colz");
ndfvsEkreco->SetTitle(Form(" %s: Ekreco %d MeV vs ndf",TConfig[0].c_str(), Ene));
ndfvsEkreco->GetYaxis()->SetTitle("Ekreco in MeV");
ndfvsEkreco->GetXaxis()->SetTitle("ndf");
fileout->cd();
ndfvsEkreco->Write();
CanndfvsEkreco->Write();
CanndfvsEkreco->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf", directory.c_str(),startfile.c_str(), type, Ene,ID.c_str()));
  //////////////////////////////////////////////////////////////////
  TCanvas*Canchi2ndf=new TCanvas("Canchi2ndf","Canchi2ndf",200,10,1200,800);
	
 
 
     Canchi2ndf->cd(1);
     chi2ndf->Draw("hist");
 gStyle->SetOptStat(1100);
     chi2ndf->SetTitle(Form("%s: %d MeV",TConfig[0].c_str(), Ene));
     chi2ndf->SetLineColor(kBlack);
     chi2ndf->GetXaxis()->SetTitle("Chi2/ndf");
     chi2ndf->GetYaxis()->SetTitle("#entries");
     fileout->cd();
     chi2ndf->Write();
     Canchi2ndf->Write();
  Canchi2ndf->Print(Form("%s/RecoTest_%s_%d_%dMeV_%s.pdf)", directory.c_str(),startfile.c_str(), type, Ene,ID.c_str()));

//close ROOT file
     fileout->Close();
 
}


