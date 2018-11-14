////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 5 , 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "LoadMCparameters.h"
#include "TChain.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TF2.h"
#include <algorithm>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TArc.h"
#include "Fit/Fitter.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)



void DisplayEventsGIF(int t, int ene, string s)
{
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<2;i++)ShellReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
  double B = 0.3;		//average magnetic field in T
  double c = TMath::C();
 string MCparamfile="./MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);

//Input file 
string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4";
string startfile="aesopliteNonUniB_V4";
string endfile="_fort.99";
string source = s;
string directory= "/home/smechbal/Documents/AESOPLITE/Analysis/MCAnalysis/EventDisplay";
string RecoInd="KFone";


//type of particle
int type=t;
string stype;	
if(t==3)	stype="e^{#minus}";
if(t==4)	stype="e^{#plus}";
if(t==11)	stype="#mu^{#minus}";
if(t==10)	stype="#mu^{#plus}";
if(t==1)	stype="protons";
if(t==6)	stype="#alpha";

//Number of energies
int Ene= ene;
 string UNIT="MeV";	
	
float mass=0.000511;//electon mass in GeV
if(t==11 || t==10)     mass=0.10566;//muon mass in GeV
if(t==1)     mass=0.93827;//proton mass in GeV
if(t==6)     mass=3.7273;//alpha-particle mass in GeV

int mcolor[4]={4,2,3,7}; 
	
//TChain to merge the files for a given energy
 TChain*chain=new TChain("MC");
 int nentries=0;
 chain->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s0001%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene,UNIT.c_str(),endfile.c_str(),RecoInd.c_str()));
 TFile *fileout=new TFile(Form("%s/DisplayEventMC_%s_%d_%s.root", directory.c_str(),startfile.c_str(), type,RecoInd.c_str()),"RECREATE");
 TCanvas *can=new TCanvas();
 
 int firstpage=0;
 ALEvent *e = new ALEvent();      
 //Define variables to read event
 //Set address to access event data
 chain->SetBranchAddress("Revent",&e); 
 
 // Get number of event in Tree
 nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;
 
 float ZMax=68;
 float ZMin=-40;
	int ncount=0;

 
 //Loop over events
 for(int i=0;i<nentries;i++)
   {
    chain->GetEntry(i);
    if(i%100000==0) cout << "Event: " << i <<endl;
    int nnhits=e->get_Nhits();
    uint8_t Ti=(uint8_t)e->get_Ti();
    //cout << "Extract number of hits done: "<<nnhits  <<endl;
    double deflection=e->get_deflecPR();

    ///////////////////////////////    
    //Load MC variables
    ///////////////////////////////
//Plot events that go through T1/T3, but not through all 7 layers
    int  nhits = e->get_Nhits();
	//check trigger requirements
	bool TShell;
	bool T1=e->get_T1();
	bool T2=e->get_T2();
	bool T3=e->get_T3();
	bool T4=e->get_T4();
	bool guard=e->get_guard();
	bool Lhit[7] = {false,false,false,false,false,false,false};
	double EFoam=0;
	double EShell=0;
	double tFoam=0;
	double tShell=0;
	double E1=0;
	double tT1=0;
	double E2=0;
	double E3=0;
	double E4=0;
	double EG=0;
	//Insulation Foam
	int nTFoam=e->get_EneIsofoam().size();
	int nTShell=e->get_EneShell().size();
	if(nTFoam>0 && nTShell>0)TShell=true;
	else TShell=false;	   
	//Number of layers with hit(s)
	int NL= e->get_NLayers();
	int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
	int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	//Reconstruction variables
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
	double e0MC;
	double p0MC;
	double cpa_MC;
	int typeMC = e->get_typeMC();
	int QMC;
	if(typeMC==3 || typeMC==11) QMC = -1;
	if(typeMC==4 || typeMC==10) QMC = 1;    

	///////////////////////////////
	//Energy Deposits
	///////////////////////////////
	//Insulation Foam
   if(nTFoam==0) continue;
	for(int j=0;j<1;j++) {
			EFoam+=1000.*e->get_EneIsofoam().at(j);
			tFoam=e->get_timeIsofoam().at(j);
	}

	//Shell
   if(nTShell==0) continue;
	for(int j=0;j<1;j++) {
			EShell+=1000.*e->get_EneShell().at(j);
	}
	//T1
	int nT1=e->get_EneT1().size();
	for(int j=0;j<nT1;j++) {
			E1+=1000.*e->get_EneT1().at(j);
	}
	//T2
	int nT2=e->get_EneT2().size();
	for(int j=0;j<nT2;j++) {E2+=1000.*e->get_EneT2().at(j);}
	//T3
	int nT3=e->get_EneT3().size();
	for(int j=0;j<nT3;j++) {E3+=1000.*e->get_EneT3().at(j);}
	//T4
	int nT4=e->get_EneT4().size();
	for(int j=0;j<nT4;j++) {E4+=1000.*e->get_EneT4().at(j);}
	//Guard
	int nTG=e->get_Eneg().size();
	for(int j=0;j<nTG;j++) {EG+=1000.*e->get_Eneg().at(j);}
    for(int j=0;j<nnhits;j++)
      {
       int L= e->get_hits().at(j)->get_L();
	   Lhit[L] = true;
	}
	 

    ///////////////////////////////    
    //Events selection
    ///////////////////////////////
	 
	//select events that FAIL to go through instrument
  //if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2] && NL < 5)
	 if(TShell && T1 && T3)

	{	 
	cout << "Event " << i << ", NL = " << NL << endl;	 
	ncount++;
    if(ncount>100) {
		cout << "100 event reached " << endl;
		break;
	}
    //////////////////////////////   
    //Extract hits information
    //////////////////////////////

    can=new TCanvas(Form("Event %d",e->get_eventnumber()),Form("Event %d",e->get_eventnumber()),10,10,1000,1000);
	can->SetTitle(Form("%s %d%s, %s", stype.c_str(), Ene, UNIT.c_str(), source.c_str()));
	 can->Divide(2,1);	
//array of particle age for each branch graph
	  
//Count number of traces
    float oldage=0;
	float tstart;
	float tend;
	int nbranch=1;
	int nZpos = e->get_posZ().size();	  
	for(int j=0;j<nZpos-1;j++)
	  {
	 float age = e->get_posAge().at(j);
	 if (age<oldage) {
		 nbranch++;
		 tstart=age;
	 }	 
	 oldage=age;
	  }
	  tend = oldage;
	  //cout << "number of traces = " << nbranch << endl;

//Save the particle's age, for each trace separately
	float**TimeBranch = new float*[nbranch];
	int*PointBranch = new int[nbranch];
	  for(int j=0;j<nbranch;j++) {
		  TimeBranch[j] =  new float[nZpos];
		  PointBranch[j]=0;
		  for(int k=0;k<nZpos;k++) {
			  TimeBranch[j][k]=0;
		  }	  
	  }
	float oldage2=0;
	int nbranch3=0;
	//Loop over unordered MC hits
	for(int j=0;j<nZpos-1;j++)
	 {	 
		float xpos = e->get_posX().at(j);
		float ypos = e->get_posY().at(j);		 
		float zpos = e->get_posZ().at(j);
		float typepos = e->get_posType().at(j);
		float age = e->get_posAge().at(j);
		float mompos = e->get_posP().at(j);
		if (age<oldage2) {nbranch3++;}
		PointBranch[nbranch3]= PointBranch[nbranch3]+1;	 
		TimeBranch[nbranch3][PointBranch[nbranch3]] = age; 
		oldage2=age;
		} // loop over all unordered MC hits
	  
//Sort vectors by chronological order
     vector<float> SortedAge = e->get_posAge();
	 int nSort = SortedAge.size();
     vector<float> SortedPosX(nSort); 
     vector<float> SortedPosY(nSort);
     vector<float> SortedPosZ(nSort);
     vector<float> SortedPosType(nSort);
	 vector<float> SortedPosMom(nSort);  
     sort(SortedAge.begin(),SortedAge.end());
     for(int j=0;j<nSort-1;j++)
	 	{
		 for(int k=0; k<nSort-1; k++)
		 	{
			 float sortedTime = SortedAge.at(j);
			 float unsortedTime = e->get_posAge().at(k);
			 float xpos = e->get_posX().at(k);
			 float ypos = e->get_posY().at(k);		 
			 float zpos = e->get_posZ().at(k);
			 float typepos = e->get_posType().at(k);
			 float age = e->get_posAge().at(k);
			 float mompos = e->get_posP().at(k);
			 if(sortedTime==unsortedTime) {
				 SortedPosX.at(j)=xpos;
				 SortedPosY.at(j)=ypos;
				 SortedPosZ.at(j)=zpos;
				 SortedPosType.at(j)=typepos;
				 SortedPosMom.at(j)=mompos;			 
			 }
		 }
	 }	  
			 
    TGraph **grALLNB=new TGraph*[nbranch];
    TGraph **grALLB=new TGraph*[nbranch];
	TPaveText**Pmom= new TPaveText*[nbranch];
    for(int j=0; j<nbranch;j++) 
	{
		grALLNB[j] = new TGraph();
		grALLNB[j]->SetMarkerStyle(2);
		grALLNB[j]->SetMarkerSize(0.3);
		grALLNB[j]->SetMarkerColor(mcolor[j]);
		grALLB[j] = new TGraph();
		grALLB[j]->SetMarkerStyle(2);
		grALLB[j]->SetMarkerSize(0.3);
		grALLB[j]->SetMarkerColor(mcolor[j]);
		Pmom[j]=new TPaveText(25,60-4*j,39,64-4*j);
	//mom[j]=new TPaveText(25,60,39,64);
		Pmom[j]->SetFillStyle(0);
		Pmom[j]->SetBorderSize(0);
		Pmom[j]->SetTextColor(mcolor[j]);
	}
	  
//ORDERED LOOP OVER ALL MC HITS 
	 int nbranch2=0;
	 int nstep=0;
	 //	  for(int j=0; j<100;j++)
	  for(int j=0;j<nZpos-1;j++)  	 		
	  {	 
	 float xpos = SortedPosX.at(j);
	 float ypos = SortedPosY.at(j);		 
	 float zpos = SortedPosZ.at(j);
	 float typepos =SortedPosType.at(j);
	 float age = SortedAge.at(j);
	 float mompos = SortedPosMom.at(j);
   // cout << "j = " << j << "t = " << age << "s, mompos = " << mompos*1000 << "MeV" << endl;
	  for(int k=0; k<nbranch;k++) {
		  for(int l=0;l< PointBranch[k];l++) {
		  if (age==TimeBranch[k][l] && age!=0) {
    //   cout << " branch " << k << ", point " << l << ", mom = " << mompos*1000 << "MeV, age = " << TimeBranch[k][l] << endl;
	 	 grALLNB[k]->SetPoint(grALLNB[k]->GetN(),xpos,zpos);
		 grALLB[k]->SetPoint(grALLB[k]->GetN(),ypos,zpos);	
		 if(typepos==3){
			 grALLNB[k]->SetMarkerStyle(2);
			 grALLB[k]->SetMarkerStyle(2);
		 }
		 else if(typepos==4){
			 grALLNB[k]->SetMarkerStyle(25);
			 grALLB[k]->SetMarkerStyle(25);
		 }
		 else if(typepos==7){
			 grALLNB[k]->SetMarkerStyle(4);
			 grALLB[k]->SetMarkerStyle(4);
		 }
		 Pmom[k]->Clear();		 
		 Pmom[k]->AddText(Form("P = %4.2f MeV ",mompos*1000));
		 if(k>nstep){
			 nbranch2++;
			 nstep++;
		  }
		 }
	 else continue;
	}
	  
    ///////////////////////////////////////////// 
    //Extract Reconstructed track information
    ////////////////////////////////////////////    
    
    TF1* fNB=new TF1("fNB","pol1",-20,20);
    fNB->SetLineColor(kBlue);
    fNB->SetLineWidth(2);
    fNB->SetLineStyle(1);    
    float zz0=0;
    //TF1* fB=new TF1("fB","pol2",-20,20);
    TF1* fB= new TF1("fB","[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,40);
    fB->FixParameter(3,zz0);
    fB->SetLineColor(kBlue);
    fB->SetLineWidth(2);
    fB->SetLineStyle(1);    
  // cout << "Extraction of non-bending plane fit parameters done" <<endl;
//Non bending plane
    float p0=e->get_interPR();
    float p1=e->get_slopePR();
	// cout << "p0 = "  << p0 << ", p1 = " << p1 <<endl;
    //if(p1==0) continue;
    fNB->FixParameter(0,-p0/p1);
    fNB->FixParameter(1,1./p1);
//Bending plane
    float a=e->get_aPR();
    float b=e->get_bPR();
    float c=e->get_cPR();
     //ax^2+bx+c
    fB->FixParameter(0,c);
    fB->FixParameter(1,b);
    fB->FixParameter(2,a); 	
	//cout << "a = "  << a << ", b = " << b << ", c = " << c << endl;
//   cout << "Extraction of bending plane fit parameters done" <<endl;

		  
    //Incoming Straight particle  
    float lim=TckZPos[1];//z position of 2nd layer
    TF1*incomingB=new TF1("incomingB","pol1",-25,25);
    incomingB->SetLineStyle(2);
    incomingB->SetLineWidth(2);
    incomingB->SetLineColor(kBlue);
    float aa=fB->Eval(lim);
    float diff=2*a*lim+2*a*zz0+b;
    incomingB->FixParameter(0,-(aa-diff*lim)/diff);
    incomingB->FixParameter(1,1./diff);
    
    TF1*incomingBout=new TF1("incomingBout","pol1",-25,25);
    incomingBout->SetLineStyle(2);
    incomingBout->SetLineWidth(2);
    incomingBout->SetLineColor(kGray);
    incomingBout->FixParameter(0,-(aa-diff*lim)/diff);
    incomingBout->FixParameter(1,1./diff);    
    
    //Outcoming Straight particle    
    float limo=TckZPos[5];//z position of 6th layer
    TF1*outcomingB=new TF1("outcomingB","pol1",-25,25);
    float aaout=fB->Eval(limo);
    float diffout=2*a*limo+2*a*zz0+b;
    outcomingB->FixParameter(0,-(aaout-diffout*limo)/diffout);
    outcomingB->FixParameter(1,1./diffout);
    outcomingB->SetLineStyle(2);
    outcomingB->SetLineWidth(2);
    outcomingB->SetLineColor(kBlue);
         
    //Signed curvature at the 3 points in and around magnets
    float zzz[3]={TckZPos[3],TckZPos[2],TckZPos[1]};
    float curv[3]={3,3,3};
    TF1* fcurv=new TF1("fcurv","2*[0]/TMath::Power(1+TMath::Power(2*[0]*x+[1],2),3./2.)",-20,20);
    fcurv->SetParameter(0,a);
    fcurv->SetParameter(1,2*a*zz0+b);
    for(int ij=0;ij<3;ij++)
      {
       curv[ij]= fcurv->Eval(zzz[ij]);   
      }  

    //Detector layout
    //Layers
    TLine**Line=new TLine*[7];

    for(int ijk=0;ijk<7;ijk++)
      {
       Line[ijk]=new TLine(-9,TckZPos[ijk],9,TckZPos[ijk]);
       Line[ijk]->SetLineColor(kBlack);  
       Line[ijk]->SetLineWidth(1);  
      }       
    //Magnets 
    TBox*boxM1=new TBox(-15,-9.61682,-6.703,-9.61682+4.59994);
    TBox*boxM2=new TBox(-15,-16.2513,-6.703,-16.2513+4.59994);
    TBox*boxM3=new TBox(6.703,-9.61682,15,-9.61682+4.59994);
    TBox*boxM4=new TBox(6.703,-16.2513,15,-16.2513+4.59994);
    boxM1->SetFillColor(kGray);
    boxM2->SetFillColor(kGray);
    boxM3->SetFillColor(kGray);
    boxM4->SetFillColor(kGray);
    //T1
    TBox*boxT1=new TBox(-13.,33.49968,13.,33.49968+0.5);
    boxT1->SetFillColor(kRed);
    if(T1)boxT1->SetFillColor(kGreen);
    TPaveText*PHT1=new TPaveText(15,32,25,36);
    PHT1->AddText(Form("T1=%4.2f MeV",E1));
    PHT1->SetFillStyle(0);PHT1->SetBorderSize(0);
    //T2
    TGraph* grT2 = new TGraph();
    grT2->SetPoint(0,-6.5,2.024);//Bottom
    grT2->SetPoint(1,6.5,2.024);//Botton
    grT2->SetPoint(2,13.5,29.964);//Top
    grT2->SetPoint(3,-13.5,29.964);//Top
    grT2->SetPoint(4,-6.5,2.024);//Bottom again
    grT2->SetFillStyle(0);
    grT2->SetLineWidth(3);
    grT2->SetLineColor(kRed);
    if(T2)grT2->SetLineColor(kGreen); 
    TPaveText*PHT2=new TPaveText(15,11,25,15);
    PHT2->AddText(Form("T2=%4.2f MeV",E2));
    PHT2->SetFillStyle(0);PHT2->SetBorderSize(0);     
	 //T3
    TBox*boxT3=new TBox(-3.5,0.,3.5,0.5);
    boxT3->SetFillColor(kRed);
    if(T3)boxT3->SetFillColor(kGreen);
    TPaveText*PHT3=new TPaveText(15,-1,25,3);
    PHT3->AddText(Form("T3=%4.2f MeV",E3));
    PHT3->SetFillStyle(0);PHT3->SetBorderSize(0);
    //Guard
    TBox*boxG1=new TBox(-13.5,-0.5588,-3.5,-0.5588+0.5);
    boxG1->SetFillColor(kRed);
    if(guard)boxG1->SetFillColor(kGreen);
    TBox*boxG2=new TBox(+3.5,-0.5588,13.5,-0.5588+0.5);
    boxG2->SetFillColor(kRed);
    if(guard)boxG2->SetFillColor(kGreen);     
    TPaveText*PHG=new TPaveText(15,-3,25,1);
    PHT1->AddText(Form("G=%4.2f MeV", EG));
    PHG->SetFillStyle(0);PHG->SetBorderSize(0);   
    //T4
    TBox*boxT4=new TBox(-18.,-25.59012,18.,-25.59012+1);
    boxT4->SetFillColor(kRed);
    if(T4)boxT4->SetFillColor(kGreen);
    TPaveText*PHT4=new TPaveText(15,-25,25,-21);
    PHT4->AddText(Form("T4=%4.2f MeV",E4));
    PHT4->SetFillStyle(0);PHT4->SetBorderSize(0);
    
	/////////////////////////////////////
	/////////////Non-bending plot////////
	/////////////////////////////////////
	
 	can->cd(1);	 
    TMultiGraph* multi=new TMultiGraph();
    multi->SetTitle(Form("Event: %d, Non bending plane",i));	 
   // multi->Add(grNB,"p");
    multi->Add(grT2,"");
    multi->Draw("a");
	for (int k=0;k<nbranch2+1;k++) {
		multi->Add(grALLNB[k],"p");
		}
	
	multi->GetXaxis()->SetTitle("X (cm)");
	multi->GetYaxis()->SetTitle("Z (cm)");
	multi->GetYaxis()->SetTitleOffset(1.1);
	multi->GetXaxis()->SetLimits(-40,40);
	multi->SetMaximum(ZMax);
	multi->SetMinimum(ZMin);
    gPad->Update();  		
    		 
	TLine * Line0=new TLine(0,ZMin,0,ZMax);
    Line0->SetLineColor(kGray);
    Line0->SetLineStyle(3);
    Line0->SetLineWidth(1);
   // Line0->Draw("same"); 
    //T1
    boxT1->Draw("l same");
    PHT1->Draw();
    //T2
    PHT2->Draw();
    //T3
    boxT3->SetLineWidth(1);
    boxT3->SetLineColor(1);
    boxT3->Draw("l same");
    PHT3->Draw();
    //Guard
    boxG1->SetLineWidth(1);
    boxG2->SetLineWidth(1);
    boxG1->SetLineColor(1);
    boxG2->SetLineColor(1);
    boxG1->Draw("l same");
    boxG2->Draw("l same");
    PHG->Draw();
    //T4
    boxT4->Draw("l same");
    PHT4->Draw();	 
    //Layers
    for(int ijk=0;ijk<7;ijk++) Line[ijk]->Draw("same");  
    //Non beding stright line
  //  if(fNB->GetParameter(1)>0) fNB->DrawF1(fNB->GetX(-30),fNB->GetX(40),"same");
  //  else fNB->DrawF1(fNB->GetX(40),fNB->GetX(-30),"same");
    //Magnet
    boxM1->Draw("l same");
    boxM2->Draw("l same");
    boxM3->Draw("l same");
    boxM4->Draw("l same");       

    TPaveText*ThetaNB;
    int sh=15;
    if(fNB->GetX(40)+sh<25)ThetaNB=new TPaveText(fNB->GetX(40),40,fNB->GetX(40)+sh,44);
    else ThetaNB=new TPaveText(fNB->GetX(40)-sh,40,fNB->GetX(40),44);
    ThetaNB->SetFillStyle(0);ThetaNB->SetBorderSize(0);
    ThetaNB->AddText(Form("#theta_{NB}= %6.4f",TMath::ATan(p1)));
   // ThetaNB->Draw();

	TPaveText*Shell=new TPaveText(15,42,35,55);
    Shell->AddText(Form("Eloss in top shell: %4.2f MeV",EShell+EFoam));
    Shell->SetFillStyle(0);Shell->SetBorderSize(0);   
    Shell->Draw();   

	TPaveText*PNL=new TPaveText(16,-10,36,-20);
	PNL->AddText(Form("%d layers with hit",NL));    
    PNL->SetFillStyle(0);PNL->SetBorderSize(0);   
    PNL->Draw();  		  

    //Legend for particle type markers
	TMarker *melec = new TMarker();
	melec->SetMarkerColor(kBlack);
	melec->SetMarkerStyle(2);
	melec->SetMarkerSize(0.3);
	TMarker *mpos = new TMarker();
	mpos->SetMarkerColor(kBlack);
	mpos->SetMarkerStyle(25);
	mpos->SetMarkerSize(0.3);
	TMarker *mphoton = new TMarker();
	mphoton->SetMarkerColor(kBlack);
	mphoton->SetMarkerStyle(4);
	mphoton->SetMarkerSize(0.3);
	 
	TLegend*Leg=new TLegend(-39,50,-15,66,"","");
	Leg->AddEntry(melec,"electron","p");
	Leg->AddEntry(mpos,"positron" ,"p");
	Leg->AddEntry(mphoton,"photon" ,"p");
    Leg->SetFillStyle(0);
	Leg->SetBorderSize(0);   
    Leg->Draw("same"); 
		  
	//Draw momentum reading box
	for (int k=0;k<nbranch2+1;k++) {
		     Pmom[k]->Draw();
			}


    /////////////////////////////////////
	/////////////Bending plot////////////
	/////////////////////////////////////
		 
	//	cout << "bending side plot" << endl;

    can->cd(2);
	TMultiGraph *multiB=new TMultiGraph(); 	 
   	multiB->SetTitle("Bending plane");
    multiB->Add(grT2,"");
    multiB->Draw("a");
    for (int k=0;k<nbranch2+1;k++) {
	    multiB->Add(grALLB[k],"p");
		}	
    multiB->GetXaxis()->SetTitle("Y (cm)");
    multiB->GetYaxis()->SetTitle("Z (cm)");
    multiB->GetYaxis()->SetTitleOffset(1.1);
    multiB->GetXaxis()->SetLimits(-40,40);
    multiB->SetMaximum(ZMax);
    multiB->SetMinimum(ZMin);
	gPad->Update();
		 	 
    //Triclk to the Inverse bending track for display
    int Nn=500;
    Double_t*x=new Double_t[Nn]; 
    Double_t*y=new Double_t[Nn];
    Double_t dx=(lim-limo)/Nn;
    Double_t x1=limo;
    x[0] = x1;
    y[0] = fB->Eval(x[0]);
    for (int ij=1; ij<Nn; ij++) 
      {
       x1   = x1+dx;
       x[ij] = x1;
       y[ij] = fB->Eval(x[ij]);
      }

    TGraph *gr = new TGraph(Nn,y,x);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);   
   // Line0->Draw("same");

    //Magnet
    boxM1->Draw("l same");
    boxM2->Draw("l same");
    boxM3->Draw("l same");
    boxM4->Draw("l same");       
    //T1
    boxT1->Draw("l same");
    //T2
    //T3
    boxT3->SetLineWidth(1);
    boxT3->SetLineColor(1);
    boxT3->Draw("l same");
    //Guard
    boxG1->SetLineWidth(1);
    boxG1->SetLineColor(1);
    boxG2->SetLineWidth(1);
    boxG2->SetLineColor(1);
    boxG1->Draw("l same");
    boxG2->Draw("l same");
    //T4
    boxT4->Draw("l same");
    //Layers
    for(int ijk=0;ijk<7;ijk++) Line[ijk]->Draw("same");      

    TPaveText*ThetaIn;
    if(incomingB->GetX(40)+sh<25)ThetaIn=new TPaveText(incomingB->GetX(40),40,incomingB->GetX(40)+sh,44);
    else ThetaIn=new TPaveText(incomingB->GetX(40)-sh,40,incomingB->GetX(40),44);
    ThetaIn->SetFillStyle(0);ThetaIn->SetBorderSize(0);
    ThetaIn->AddText(Form("#theta_{Bin}= %6.3f",TMath::ATan(diff)));
    
    TPaveText*ThetaOut;
     if(outcomingB->GetX(-30)+sh<25) ThetaOut=new TPaveText(outcomingB->GetX(-30),-34,outcomingB->GetX(-30)+sh,-30);
     else  ThetaOut=new TPaveText(outcomingB->GetX(-30)-sh,-34,outcomingB->GetX(-30),-30);
    ThetaOut->SetFillStyle(0);ThetaOut->SetBorderSize(0);
    ThetaOut->AddText(Form("#theta_{Bout}= %6.3f",TMath::ATan(diffout)));
  //  ThetaIn->Draw();
    //ThetaOut->Draw();
    
    TPaveText*Def=new TPaveText(-23,-38,10,-33);
    //Def->AddText(Form("Deflection: %6.4f",(float)deflection));
    float Def2=0;
    Def2=TMath::ATan(diffout)-TMath::ATan(diff);
    deflection=Def2;
    Def->AddText(Form("Deflection: #theta_{Bout}-#theta_{Bin}= %6.3f",(float)deflection));
    Def->SetFillStyle(0);Def->SetBorderSize(0);
  //  Def->Draw();
    TPaveText*Cur=new TPaveText(10,-23,25,-16);
    Cur->AddText(Form("Curvature at L1=%5.4f",curv[2]));
    Cur->AddText(Form("Curvature at L2=%5.4f",curv[1]));
    Cur->AddText(Form("Curvature at L3=%5.4f",curv[0]));
    Cur->SetFillStyle(0);Cur->SetBorderSize(0);   
    //Cur->Draw();
    //if(incomingB->GetX(33.25)<-14 ||incomingB->GetX(33.25)>14)continue;
    
	 } //end k
		  
	////////////////////////////////////
	////Print and save canvas///////////
	////////////////////////////////////
		  
	can->Update();		
   // can->Print(Form("%s/EventDisplayMC_%d_%s_Event%d_%s.gif+10", directory.c_str(),type,source.c_str(),i,RecoInd.c_str()));
   // can->Print(Form("%s/EventDisplayMC_%d_%s_Event%d_%d_%s.png", directory.c_str(),type,source.c_str(),i,ncount,RecoInd.c_str())); 	
	  } //END LOOP OVER ALL EVENT HITS
	  
  // cout << "done with loop over all hits" << endl;
  // can->Print(Form("%s/EventDisplayMC_%d_%s_Event%d_%s.gif++", directory.c_str(),type,source.c_str(),i,RecoInd.c_str()));
   if(firstpage==0) can->Print(Form("%s/EventDisplayMC_%d_%s_%d%s_%s.pdf(", directory.c_str(),type,source.c_str(),Ene, UNIT.c_str(),RecoInd.c_str()),Form("Title:Event %d",i));
   else can->Print(Form("%s/EventDisplayMC_%d_%s_%d%s_%s.pdf", directory.c_str(),type,source.c_str(),Ene, UNIT.c_str(),RecoInd.c_str()),Form("Title:Event %d",i));
   fileout->cd();
   can->Write();
   firstpage=-1;
   delete can;  

	} // END OF EVENT SELECTION CONDITION
 }//i
   
 fileout->Close();  
 can=new TCanvas();
 can->Print(Form("%s/EventDisplayMC_%d_%s_%d%s_%s.pdf)", directory.c_str(),type,source.c_str(),Ene, UNIT.c_str(),RecoInd.c_str()));

}//end function


