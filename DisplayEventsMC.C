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

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TArc.h"
#include "Fit/Fitter.h"

ClassImp(ALTckhit)
ClassImp(ALEvent)

void DisplayEventsDefault(string file,string  coinc,int);
void DisplayEvents(string file,string tag,string  coinc,int);
void DisplayOneEvent(string file,string tag,string  coinc,int,int);

void TrackinTrigger(string file,string  coinc,int);


int doy(int y,int m,int d);


int doy(int y,int m,int d)
{
 int iuruu=0;
 int day=0;
 int monday[2][12]={{31,28,31,30,31,30,31,31,30,31,30,31},
                    {31,29,31,30,31,30,31,31,30,31,30,31}};
 if((y % 4) ==0 )iuruu=1 ;
 if(m!=1)
  {
   for(int i=0;i<m-1;i++)
    {
       day+=monday[iuruu][i]; 
    }  //m
  }//if
 day+=d; 
 return day;
}


void DisplayEventsMC(int t, int ene)
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
string Inppath="/home/sarah/AESOPLITE/MCProduction/Detector/ProcessedFiles";	
string startfile="aesopliteNonUniB_V4";
string endfile="_fort.99";
string source = "test";
string directory= "/home/sarah/AESOPLITE/Analysis/MCAnalysis/Efficiencies/V4";
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

	
//TChain to merge the files for a given energy
 TChain*chain=new TChain("MC");
 int Nevents=0;
 chain->Add(Form("%s/%d/%s/RecoEvent_%s_%d_%d%s%s_%s.root", Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene,UNIT.c_str(),endfile.c_str(),RecoInd.c_str()));
 TFile *fileout=new TFile(Form("%s/%s/DisplayEventMC_%s_%d_%s.root", directory.c_str(),source.c_str(),startfile.c_str(), type,RecoInd.c_str()),"RECREATE");
 TCanvas *can=new TCanvas();
 
 int firstpage=0;
 ALEvent *e = new ALEvent();      
 //Define variables to read event
 //Set address to access event data
 chain->SetBranchAddress("Revent",&e); 
 
 // Get number of event in Tree
 Nevents=chain->GetEntries();
 cout << "Number  of events: " << Nevents << endl;
 
 float ZMax=50;
 float ZMin=-40;
 
 //Loop over events
 for(int i=0;i<Nevents;i++)
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
	 
	//All positions vectors
	 

	 ///////////////////////////////    
    //Events selection
    ///////////////////////////////
  if(TShell && T1 && T3 && Lhit[0] && Lhit[1] && Lhit[2]  && Lhit[3]  && Lhit[4]  && (!Lhit[5] || Lhit[6]))
	//if(TShell && T1 && T3 && NL > 5)


	 {	 
	//cout << "Event " << i << ", NL = " << NL << endl;	 


    ///////////////////////////////    
    //Extract hits information
    ///////////////////////////////
   
    TGraph *grALLNB=new TGraph();
    grALLNB->SetMarkerStyle(5);
    grALLNB->SetMarkerColor(kBlack);
    TGraph *grALLB=new TGraph();
    grALLB->SetMarkerStyle(5);
    grALLB->SetMarkerColor(kBlack);
    TGraph *grNB=new TGraph();
    grNB->SetMarkerStyle(kCircle);
    grNB->SetMarkerColor(kBlue);
    TGraph *grB=new TGraph();
    grB->SetMarkerStyle(kCircle);
    grB->SetMarkerColor(kBlue);

    for(int j=0;j<nnhits;j++)
      {
       int L= e->get_hits().at(j)->get_L();
       float x=e->get_hits().at(j)->get_x();
       float y=e->get_hits().at(j)->get_y();
       float z=e->get_hits().at(j)->get_z();
	   bool flagPR=e->get_hits().at(j)->get_flagPR();
       if(L==0 ||L==4||L==6)//Non bending plane
        {
         if(flagPR)grNB->SetPoint(grNB->GetN(),x,z);
        }
       else if(L==1 ||L==2||L==3 ||L==5) //Bending plane
        {
         if(flagPR)grB->SetPoint(grB->GetN(),y,z);
        }
      } 	 
	
//LOOP OVER ALL MC HITS 		 

	int nZpos = e->get_posZ().size();
	for(int j=0;j<nZpos;j++)
	  {
	 float xpos = e->get_posX().at(j);
	 float ypos = e->get_posY().at(j);		 
	 float zpos = e->get_posZ().at(j);
	 cout << "xpos  = " << xpos << ", ypos = " << ypos << ", zpos = " << zpos << endl;
	 grALLNB->SetPoint(grALLNB->GetN(),xpos,zpos);
	 grALLB->SetPoint(grALLB->GetN(),ypos,zpos);
	}
	// cout << "done with hit extraction" << endl;
  
    ///////////////////////////////    
    //Extract Reconstructed track information
    ///////////////////////////////    
    
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
    
    ///////////////////////////////    
    //Event Display 
    ///////////////////////////////
    can=new TCanvas(Form("Event %d",e->get_eventnumber()),Form("Event %d",e->get_eventnumber()),200,10,1600,1440);
    can->Divide(2,1);
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
    // cout << "could it be here?" << endl;
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
    
    //Non Bending plot   
    can->cd(1);
    
    TMultiGraph* multi=new TMultiGraph();multi->SetTitle(Form("Event: %d, Non bending plane",i));
    multi->Add(grALLNB,"p");
    multi->Add(grNB,"p");
    multi->Add(grT2,"");
    multi->Draw("a");
    multi->GetXaxis()->SetTitle("X (cm)");
    multi->GetYaxis()->SetTitle("Z (cm)");
    multi->GetYaxis()->SetTitleOffset(1.1);
    multi->GetXaxis()->SetLimits(-25,25);
    multi->SetMaximum(ZMax);
    multi->SetMinimum(ZMin);
    gPad->Update();
    
    TLine * Line0=new TLine(0,ZMin,0,ZMax);
    Line0->SetLineColor(kGray);
    Line0->SetLineStyle(3);
    Line0->SetLineWidth(1);
    Line0->Draw("same");
    
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
    if(fNB->GetParameter(1)>0) fNB->DrawF1(fNB->GetX(-30),fNB->GetX(40),"same");
    else fNB->DrawF1(fNB->GetX(40),fNB->GetX(-30),"same");
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
    ThetaNB->Draw();
    
	TPaveText*Shell=new TPaveText(-23,45,10,50);
    Shell->AddText(Form("Energy loss in top shell: %4.2f MeV",EShell+EFoam));
    Shell->AddText(Form("Number of layers with hit NL = %d",NL));    
    Shell->SetFillStyle(0);Shell->SetBorderSize(0);
    Shell->Draw();   
	 
	//	cout << "bending side plot" << endl;
    //Bending plot   
    can->cd(2);
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

    TMultiGraph *multiB=new TMultiGraph(); multiB->SetTitle("Bending plane");
    multiB->Add(grALLB,"p");
    multiB->Add(grB,"p");
    multiB->Add(grT2,"");
    multiB->Add(gr,"l");
    multiB->Draw("a");
    multiB->GetXaxis()->SetTitle("Y (cm)");
    multiB->GetYaxis()->SetTitle("Z (cm)");
    multiB->GetYaxis()->SetTitleOffset(1.1);
    multiB->GetXaxis()->SetLimits(-25,25);
    multiB->SetMaximum(ZMax);
    multiB->SetMinimum(ZMin);
    gPad->Update();
    
    Line0->Draw("same");

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

    if(incomingB->GetParameter(1)>0) incomingB->DrawF1(incomingB->GetX(lim),incomingB->GetX(40),"same");
    else incomingB->DrawF1(incomingB->GetX(40),incomingB->GetX(lim),"same");
    
    if(incomingBout->GetParameter(1)>0) incomingBout->DrawF1(incomingBout->GetX(lim),incomingBout->GetX(-30),"same");
    else incomingBout->DrawF1(incomingBout->GetX(-30),incomingBout->GetX(lim),"same");
       
    if(outcomingB->GetParameter(1)<0) outcomingB->DrawF1(outcomingB->GetX(limo),outcomingB->GetX(-30),"same");
    else outcomingB->DrawF1(outcomingB->GetX(-30),outcomingB->GetX(limo),"same");

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
    ThetaIn->Draw();
    ThetaOut->Draw();
    
    TPaveText*Def=new TPaveText(-23,-38,10,-33);
    //Def->AddText(Form("Deflection: %6.4f",(float)deflection));
    float Def2=0;
    Def2=TMath::ATan(diffout)-TMath::ATan(diff);
    deflection=Def2;
    Def->AddText(Form("Deflection: #theta_{Bout}-#theta_{Bin}= %6.3f",(float)deflection));
    Def->SetFillStyle(0);Def->SetBorderSize(0);
    Def->Draw();
    TPaveText*Cur=new TPaveText(10,-23,25,-16);
    Cur->AddText(Form("Curvature at L1=%5.4f",curv[2]));
    Cur->AddText(Form("Curvature at L2=%5.4f",curv[1]));
    Cur->AddText(Form("Curvature at L3=%5.4f",curv[0]));
    Cur->SetFillStyle(0);Cur->SetBorderSize(0);   
    //Cur->Draw();
   
    //if(incomingB->GetX(33.25)<-14 ||incomingB->GetX(33.25)>14)continue;
    fileout->cd();
    can->Write();
    
    if(firstpage==0) can->Print(Form("%s/EventDisplayMC_%d_%s_%s.pdf(", directory.c_str(),type,source.c_str(),RecoInd.c_str()),Form("Title:Event %d",i));
    else can->Print(Form("%s/EventDisplayMC_%d_%s_%s.pdf", directory.c_str(),type,source.c_str(),RecoInd.c_str()),Form("Title:Event %d",i));

	firstpage=-1;
    delete can;
	} //end of event selection condition
   }//i
   
 fileout->Close();
 can=new TCanvas();
 can->Print(Form("%s/EventDisplayMC_%d_%s_%s.pdf)", directory.c_str(),type,source.c_str(),RecoInd.c_str()));

}//end function


void DisplayEvent(string file,string tag,string coinc,int geoconf,int EVTNum)
{
  //Load configuration parameter
 float* TckZPos=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)TckZPos[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

 LoadDataparameters(paramfile,TckZPos,OffsetLL,OffsetRL,TrigThresh);

 for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", TckZPos:" << TckZPos[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 cout << "T1 threshold: " << TrigThresh[0] <<endl;
 cout << "T2 threshold: " << TrigThresh[1] <<endl;
 cout << "T3 threshold: " << TrigThresh[2] <<endl;
 cout << "T4 threshold: " << TrigThresh[3] <<endl;
 cout << "Guard threshold: " << TrigThresh[4] <<endl;
    
    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");
 
 int Nevents=0;
 
 cout << Form("../Data/%s.EVENT_%s.root",file.c_str(),tag.c_str()) <<endl;
 chain->Add(Form("../Data/%s.EVENT_%s.root",file.c_str(),tag.c_str()));
 
 TFile*fileout=new TFile(Form("%s_%s.EVTNum%d.root",file.c_str(),tag.c_str(),EVTNum),"RECREATE");
 TCanvas *can=new TCanvas();
 
 int firstpage=0;
 ALEvent *e = new ALEvent();      
 //Define variables to read event
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 Nevents=chain->GetEntries();
 cout << "Number  of events: " << Nevents << endl;
 
 float ZMax=50;
 float ZMin=-40;
 int* Nh=new int[7];

 //Loop over events
 for(int i=0;i<Nevents;i++)
   {
    chain->GetEntry(i);
    if(i%10000==0) cout << "Event: " << i <<endl;
   // if(i!=757-1)continue;
    
    int nnhits=e->get_Nhits();
    uint8_t Ti=(uint8_t)e->get_Ti();
    //Number of layers wih hit(s)
    int NL=0;
    for(int i=0;i<7;i++) NL+=(int)((Ti >>i) & 0x01);
     
    if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
    int* Lay=new int[7];
    for(int i=0;i<7;i++) Lay[i]=(int)((Ti >>i) & 0x01);    //cout << "Extract number of hits done: "<<nnhits  <<endl;
    double deflection=e->get_deflecPR();
    ///////////////////////////////    
    //Selection of events (Cuts)
    ///////////////////////////////
    //if(nnhits<10)continue;
  //  if(NL<5) continue;
 //   if(e->get_chi2NBPR()<0) continue;
//    if(e->get_chi2BPR()<0) continue;
//    if(e->get_chi2()>0) continue;
    //cout << "Extract internal tracker trigger information done" <<endl;
    //if(abs(deflection)>2||abs(deflection)<1.5) continue;
    //if(e->get_EneT2().at(0)>0) continue;
   // if(e->get_EneT3().at(0)<0) continue;
   // if(e->get_Eneg().at(0)>0) continue;
    //if(e->get_chi2NB()<=0||e->get_chi2NB()<10) continue;
   // if(e->get_chi2B()<=0||e->get_chi2B()<10) continue;
    
    ///////////////////////////////    
    //Extract hits information
    ///////////////////////////////
   
    TGraph *grALLNB=new TGraph();
    grALLNB->SetMarkerStyle(5);
    grALLNB->SetMarkerColor(kBlack);
    TGraph *grALLB=new TGraph();
    grALLB->SetMarkerStyle(5);
    grALLB->SetMarkerColor(kBlack);
    TGraph *grNB=new TGraph();
    grNB->SetMarkerStyle(kCircle);
    grNB->SetMarkerColor(kBlue);
    TGraph *grB=new TGraph();
    grB->SetMarkerStyle(kCircle);
    grB->SetMarkerColor(kBlue);
    //cout << "Graphs are defined" <<endl;
    for(int j=0;j<7;j++)Nh[j]=0;

    for(int j=0;j<nnhits;j++)
      {
       int L= e->get_hits().at(j)->get_L();
       float x=e->get_hits().at(j)->get_x();
       float y=e->get_hits().at(j)->get_y();
       float z=e->get_hits().at(j)->get_z();
       bool flagPR=e->get_hits().at(j)->get_flagPR();
       Nh[L]++;
       if(L==0 ||L==4||L==6)//Non bending plane
        {
         grALLNB->SetPoint(grALLNB->GetN(),x,z);
         if(flagPR)grNB->SetPoint(grNB->GetN(),x,z);
        }
       else if(L==1 ||L==2||L==3 ||L==5) //Bending plane
        {
         grALLB->SetPoint(grALLB->GetN(),y,z);
         if(flagPR)grB->SetPoint(grB->GetN(),y,z);
        }
      }  

    ///////////////////////////////    
    //Fit Circle in the bending plane
    ///////////////////////////////     
      
    auto chi2Function = [&](const Double_t *par)
      {
       //minimisation function computing the sum of squares of residuals
       // looping at the graph points
       Int_t np = grB->GetN();
       Double_t f = 0;
       Double_t *x = grB->GetX();
       Double_t *y = grB->GetY();
       for (Int_t j=0;j<np;j++) {
          Double_t u = x[j] - par[0];
          Double_t v = y[j] - par[1];
          Double_t dr = par[2] - std::sqrt(u*u+v*v);
          f += dr*dr;
       }
       return f;
      };     
    // wrap chi2 function in a function object for the fit
    // 3 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function,3);
    ROOT::Fit::Fitter  fitter;
    double pStart[3] = {0,0,1};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("x0");
    fitter.Config().ParSettings(1).SetName("y0");
    fitter.Config().ParSettings(2).SetName("R");
    // do the fit 
    bool ok = fitter.FitFCN();
    if (!ok)
     {
      Error("line3Dfit","Line3D Fit failed");
     }   
    const ROOT::Fit::FitResult & result = fitter.Result();
    result.Print(std::cout);      
    
    //Draw the circle on top of the points
    TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
    arc->SetLineColor(kRed);
    arc->SetLineWidth(1);
    arc->SetFillStyle(0);
    //arc->DrawArc();
      
    ///////////////////////////////    
    //Extract Reconstructed track information
    ///////////////////////////////    
    
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

    //Non bending plane
    float p0=e->get_interPR();
    float p1=e->get_slopePR();
    if(p1==0) continue;
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

    //cout << "Extraction of bending plane fit parameters done" <<endl;

    
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
    
    ///////////////////////////////    
    //Event Display 
    ///////////////////////////////
    can=new TCanvas(Form("Event %d",e->get_eventnumber()),Form("Event %d",e->get_eventnumber()),200,10,1600,1440);
    can->Divide(2,1);
    
    TPaveText**TNh=new TPaveText*[7];
    
    for(int ijk=0;ijk<7;ijk++)
      {
       TNh[ijk]=new TPaveText(-20,TckZPos[ijk]-2,-18,TckZPos[ijk]+2);
       TNh[ijk]->SetFillStyle(0);TNh[ijk]->SetBorderSize(0);
       TNh[ijk]->SetTextColor(kBlack);
       TNh[ijk]->AddText(Form("%d",Nh[ijk]));
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
    if(e->get_EneT1().at(0)>0)boxT1->SetFillColor(kGreen);
    TPaveText*PHT1=new TPaveText(15,32,25,36);
    PHT1->AddText(Form("PHT1=%d",(int)e->get_EneT1().at(0)));
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
    if(e->get_EneT2().at(0)>0)grT2->SetLineColor(kGreen); 
    TPaveText*PHT2=new TPaveText(15,11,25,15);
    PHT2->AddText(Form("PHT2=%d",(int)e->get_EneT2().at(0)));
    PHT2->SetFillStyle(0);PHT2->SetBorderSize(0);     
    //T3
    TBox*boxT3=new TBox(-3.5,0.,3.5,0.5);
    boxT3->SetFillColor(kRed);
    if(e->get_EneT3().at(0)>0)boxT3->SetFillColor(kGreen);
    TPaveText*PHT3=new TPaveText(15,-1,25,3);
    PHT3->AddText(Form("PHT3=%d",(int)e->get_EneT3().at(0)));
    PHT3->SetFillStyle(0);PHT3->SetBorderSize(0);
    //Guard
    TBox*boxG1=new TBox(-13.5,-0.5588,-3.5,-0.5588+0.5);
    boxG1->SetFillColor(kRed);
    if(e->get_Eneg().at(0)>0)boxG1->SetFillColor(kGreen);
    TBox*boxG2=new TBox(+3.5,-0.5588,13.5,-0.5588+0.5);
    boxG2->SetFillColor(kRed);
    if(e->get_Eneg().at(0)>0)boxG2->SetFillColor(kGreen);     
    TPaveText*PHG=new TPaveText(15,-3,25,1);
    PHG->AddText(Form("PHG=%d",(int)e->get_Eneg().at(0)));
    PHG->SetFillStyle(0);PHG->SetBorderSize(0);
    
    //T4
    TBox*boxT4=new TBox(-18.,-25.59012,18.,-25.59012+1);
    boxT4->SetFillColor(kRed);
    if(e->get_EneT4().at(0)>0)boxT4->SetFillColor(kGreen);
    TPaveText*PHT4=new TPaveText(15,-25,25,-21);
    PHT4->AddText(Form("PHT4=%d",(int)e->get_EneT4().at(0)));
    PHT4->SetFillStyle(0);PHT4->SetBorderSize(0);
    
    //Non Bending plot   
    can->cd(1);
    
    TMultiGraph* multi=new TMultiGraph();multi->SetTitle(Form("Event: %d, Non bending plane",(int)e->get_eventnumber()));
    multi->Add(grALLNB,"p");
    multi->Add(grNB,"p");
    multi->Add(grT2,"");
    multi->Draw("a");
    multi->GetXaxis()->SetTitle("X (cm)");
    multi->GetYaxis()->SetTitle("Z (cm)");
    multi->GetXaxis()->SetLimits(-25,25);
    multi->SetMaximum(ZMax);
    multi->SetMinimum(ZMin);
    gPad->Update();
    
    TLine * Line0=new TLine(0,ZMin,0,ZMax);
    Line0->SetLineColor(kGray);
    Line0->SetLineStyle(3);
    Line0->SetLineWidth(1);
    Line0->Draw("same");
    
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
    if(fNB->GetParameter(1)>0) fNB->DrawF1(fNB->GetX(-30),fNB->GetX(40),"same");
    else fNB->DrawF1(fNB->GetX(40),fNB->GetX(-30),"same");
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
    ThetaNB->Draw();

    TPaveText*Tchi2NB;
    Tchi2NB=new TPaveText(15,44,25,50);
    Tchi2NB->SetFillStyle(0);Tchi2NB->SetBorderSize(0);
    Tchi2NB->SetTextColor(kBlue);
    Tchi2NB->AddText(Form("#chi^{2}= %6.4f",e->get_chi2NBPR()));
    Tchi2NB->Draw();
    
    TNh[0]->Draw();
    TNh[4]->Draw();
    TNh[6]->Draw();
    
    
    
    //Bending plot   
    can->cd(2);
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

    TMultiGraph *multiB=new TMultiGraph(); multiB->SetTitle("Bending plane");
    multiB->Add(grALLB,"p");
    multiB->Add(grB,"p");
    multiB->Add(grT2,"");
    multiB->Add(gr,"l");
    multiB->Draw("a");
    multiB->GetXaxis()->SetTitle("Y (cm)");
    multiB->GetYaxis()->SetTitle("Z (cm)");
    multiB->GetXaxis()->SetLimits(-25,25);
    multiB->SetMaximum(ZMax);
    multiB->SetMinimum(ZMin);
    gPad->Update();
    
    Line0->Draw("same");

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

    if(incomingB->GetParameter(1)>0) incomingB->DrawF1(incomingB->GetX(lim),incomingB->GetX(40),"same");
    else incomingB->DrawF1(incomingB->GetX(40),incomingB->GetX(lim),"same");
    
    if(incomingBout->GetParameter(1)>0) incomingBout->DrawF1(incomingBout->GetX(lim),incomingBout->GetX(-30),"same");
    else incomingBout->DrawF1(incomingBout->GetX(-30),incomingBout->GetX(lim),"same");
       
    if(outcomingB->GetParameter(1)<0) outcomingB->DrawF1(outcomingB->GetX(limo),outcomingB->GetX(-30),"same");
    else outcomingB->DrawF1(outcomingB->GetX(-30),outcomingB->GetX(limo),"same");

    TPaveText*ThetaIn;
    if(incomingB->GetX(40)+sh<25)ThetaIn=new TPaveText(incomingB->GetX(40),40,incomingB->GetX(40)+sh,44);
    else ThetaIn=new TPaveText(incomingB->GetX(40)-sh,40,incomingB->GetX(40),44);
    ThetaIn->SetFillStyle(0);ThetaIn->SetBorderSize(0);
    ThetaIn->AddText(Form("#theta_{Bin}= %6.4f",TMath::ATan(diff)));
    
    TPaveText*ThetaOut;
     if(outcomingB->GetX(-30)+sh<25) ThetaOut=new TPaveText(outcomingB->GetX(-30),-34,outcomingB->GetX(-30)+sh,-30);
     else  ThetaOut=new TPaveText(outcomingB->GetX(-30)-sh,-34,outcomingB->GetX(-30),-30);
    ThetaOut->SetFillStyle(0);ThetaOut->SetBorderSize(0);
    ThetaOut->AddText(Form("#theta_{Bout}= %6.4f",TMath::ATan(diffout)));
    ThetaIn->Draw();
    ThetaOut->Draw();
    
    TPaveText*Def=new TPaveText(-23,-38,10,-33);
    //Def->AddText(Form("Deflection: %6.4f",(float)deflection));
    float Def2=0;
    Def2=TMath::ATan(diffout)-TMath::ATan(diff);
    deflection=Def2;
    Def->AddText(Form("Deflection: #theta_{Bout}-#theta_{Bin}= %6.4f",(float)deflection));
    Def->SetFillStyle(0);Def->SetBorderSize(0);
    Def->Draw();
    TPaveText*Cur=new TPaveText(10,-23,25,-16);
    Cur->AddText(Form("Curvature at L1=%5.4f",curv[2]));
    Cur->AddText(Form("Curvature at L2=%5.4f",curv[1]));
    Cur->AddText(Form("Curvature at L3=%5.4f",curv[0]));
    Cur->SetFillStyle(0);Cur->SetBorderSize(0);
    //Cur->Draw();
    
    //arc->Draw();
    //TPaveText*RArc=new TPaveText(-20,45,-5,49);
    //RArc->AddText(Form("R_{circ}= %6.4f cm",result.Parameter(2)));
    //RArc->SetFillStyle(0);RArc->SetBorderSize(0);RArc->SetTextColor(kRed);
    //RArc->Draw();   
    
    
    TPaveText*Tchi2B;
    Tchi2B=new TPaveText(15,44,25,50);
    Tchi2B->SetFillStyle(0);Tchi2B->SetBorderSize(0);
    Tchi2B->SetTextColor(kBlue);
    Tchi2B->AddText(Form("#chi^{2}= %6.4f",e->get_chi2BPR()));
    Tchi2B->Draw();

    TNh[1]->Draw();
    TNh[2]->Draw();
    TNh[3]->Draw();
    TNh[5]->Draw();
   
    
    //if(incomingB->GetX(33.25)<-14 ||incomingB->GetX(33.25)>14)continue;
    fileout->cd();
    can->Write();
    //if(firstpage==0)
    // can->Print(Form("%s.EVENTDISPLAYCirc.pdf(",file.c_str()),Form("Title:Event %d",e->get_eventnumber()));
   // else can->Print(Form("%s.EVENTDISPLAYCirc.pdf",file.c_str()),Form("Title:Event %d",e->get_eventnumber()));
    //firstpage=-1;
    delete can;

   }//i
   
 fileout->Close();  
 //can=new TCanvas();
 //can->Print(Form("%s.EVENTDISPLAYCirc.pdf)",file.c_str()));

}//end function

void DisplayOneEvent(string file,string tag,string coinc,int geoconf,int EVTNum)
{
  //Load configuration parameter
 float* TckZPos=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)TckZPos[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

 LoadDataparameters(paramfile,TckZPos,OffsetLL,OffsetRL,TrigThresh);

 for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", TckZPos:" << TckZPos[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 cout << "T1 threshold: " << TrigThresh[0] <<endl;
 cout << "T2 threshold: " << TrigThresh[1] <<endl;
 cout << "T3 threshold: " << TrigThresh[2] <<endl;
 cout << "T4 threshold: " << TrigThresh[3] <<endl;
 cout << "Guard threshold: " << TrigThresh[4] <<endl;
    
    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");
 
 int Nevents=0;
 
 cout << Form("../Data/%s.EVENT_%s.root",file.c_str(),tag.c_str()) <<endl;
 chain->Add(Form("../Data/%s.EVENT_%s.root",file.c_str(),tag.c_str()));
 
 TFile*fileout=new TFile(Form("%s_%s.EVTNum%d.root",file.c_str(),tag.c_str(),EVTNum),"RECREATE");
 TCanvas *can=new TCanvas();
 
 int firstpage=0;
 ALEvent *e = new ALEvent();      
 //Define variables to read event
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 Nevents=chain->GetEntries();
 cout << "Number  of events: " << Nevents << endl;
 
 float ZMax=50;
 float ZMin=-40;
 int* Nh=new int[7];

 //Loop over events
 for(int i=0;i<Nevents;i++)
   {
    chain->GetEntry(i);
    if(i%10000==0) cout << "Event: " << i <<endl;
   // if(i!=757-1)continue;
    if((int)e->get_eventnumber()!=EVTNum) continue;
    int nnhits=e->get_Nhits();
    uint8_t Ti=(uint8_t)e->get_Ti();
    //Number of layers wih hit(s)
    int NL=0;
    for(int i=0;i<7;i++) NL+=(int)((Ti >>i) & 0x01);
     
    if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
    int* Lay=new int[7];
    for(int i=0;i<7;i++) Lay[i]=(int)((Ti >>i) & 0x01);    //cout << "Extract number of hits done: "<<nnhits  <<endl;
    double deflection=e->get_deflecPR();

    ///////////////////////////////    
    //Selection of events (Cuts)
    ///////////////////////////////
//    if(nnhits<10)continue;
  //   if(NL<5) continue;
    //cout << "Extract internal tracker trigger information done" <<endl;
    //if(abs(deflection)>2||abs(deflection)<1.5) continue;
    //if(e->get_EneT2().at(0)>0) continue;
   // if(e->get_EneT3().at(0)<0) continue;
   // if(e->get_Eneg().at(0)>0) continue;
    //if(e->get_chi2NB()<=0||e->get_chi2NB()<10) continue;
   // if(e->get_chi2B()<=0||e->get_chi2B()<10) continue;
    
    ///////////////////////////////    
    //Extract hits information
    ///////////////////////////////
   
    TGraph *grALLNB=new TGraph();
    grALLNB->SetMarkerStyle(5);
    grALLNB->SetMarkerColor(kBlack);
    TGraph *grALLB=new TGraph();
    grALLB->SetMarkerStyle(5);
    grALLB->SetMarkerColor(kBlack);
    TGraph *grNB=new TGraph();
    grNB->SetMarkerStyle(kCircle);
    grNB->SetMarkerColor(kBlue);
    TGraph *grB=new TGraph();
    grB->SetMarkerStyle(kCircle);
    grB->SetMarkerColor(kBlue);
    //cout << "Graphs are defined" <<endl;
    for(int j=0;j<7;j++)Nh[j]=0;

    for(int j=0;j<nnhits;j++)
      {
       int L= e->get_hits().at(j)->get_L();
       float x=e->get_hits().at(j)->get_x();
       float y=e->get_hits().at(j)->get_y();
       float z=e->get_hits().at(j)->get_z();
       bool flagPR=e->get_hits().at(j)->get_flagPR();
       Nh[L]++;
       if(L==0 ||L==4||L==6)//Non bending plane
        {
         grALLNB->SetPoint(grALLNB->GetN(),x,z);
         if(flagPR)grNB->SetPoint(grNB->GetN(),x,z);
        }
       else if(L==1 ||L==2||L==3 ||L==5) //Bending plane
        {
         grALLB->SetPoint(grALLB->GetN(),y,z);
         if(flagPR)grB->SetPoint(grB->GetN(),y,z);
        }
      }  

    ///////////////////////////////    
    //Fit Circle in the bending plane
    ///////////////////////////////     
      
    auto chi2Function = [&](const Double_t *par)
      {
       //minimisation function computing the sum of squares of residuals
       // looping at the graph points
       Int_t np = grB->GetN();
       Double_t f = 0;
       Double_t *x = grB->GetX();
       Double_t *y = grB->GetY();
       for (Int_t j=0;j<np;j++) {
          Double_t u = x[j] - par[0];
          Double_t v = y[j] - par[1];
          Double_t dr = par[2] - std::sqrt(u*u+v*v);
          f += dr*dr;
       }
       return f;
      };     
    // wrap chi2 function in a function object for the fit
    // 3 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function,3);
    ROOT::Fit::Fitter  fitter;
    double pStart[3] = {0,0,1};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("x0");
    fitter.Config().ParSettings(1).SetName("y0");
    fitter.Config().ParSettings(2).SetName("R");
    // do the fit 
    bool ok = fitter.FitFCN();
    if (!ok)
     {
      Error("line3Dfit","Line3D Fit failed");
     }   
    const ROOT::Fit::FitResult & result = fitter.Result();
    result.Print(std::cout);      
    
    //Draw the circle on top of the points
    TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
    arc->SetLineColor(kRed);
    arc->SetLineWidth(1);
    arc->SetFillStyle(0);
    //arc->DrawArc();
      
    ///////////////////////////////    
    //Extract Reconstructed track information
    ///////////////////////////////    
    
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

    //Non bending plane
    float p0=e->get_interPR();
    float p1=e->get_slopePR();
    if(p1==0) continue;
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

    //cout << "Extraction of bending plane fit parameters done" <<endl;

    
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
    
    ///////////////////////////////    
    //Event Display 
    ///////////////////////////////
    can=new TCanvas(Form("Event %d",e->get_eventnumber()),Form("Event %d",e->get_eventnumber()),200,10,1600,1440);
    can->Divide(2,1);
    
    TPaveText**TNh=new TPaveText*[7];
    
    for(int ijk=0;ijk<7;ijk++)
      {
       TNh[ijk]=new TPaveText(-20,TckZPos[ijk]-2,-18,TckZPos[ijk]+2);
       TNh[ijk]->SetFillStyle(0);TNh[ijk]->SetBorderSize(0);
       TNh[ijk]->SetTextColor(kBlack);
       TNh[ijk]->AddText(Form("%d",Nh[ijk]));
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
    if(e->get_EneT1().at(0)>0)boxT1->SetFillColor(kGreen);
    TPaveText*PHT1=new TPaveText(15,32,25,36);
    PHT1->AddText(Form("PHT1=%d",(int)e->get_EneT1().at(0)));
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
    if(e->get_EneT2().at(0)>0)grT2->SetLineColor(kGreen); 
    TPaveText*PHT2=new TPaveText(15,11,25,15);
    PHT2->AddText(Form("PHT2=%d",(int)e->get_EneT2().at(0)));
    PHT2->SetFillStyle(0);PHT2->SetBorderSize(0);     
    //T3
    TBox*boxT3=new TBox(-3.5,0.,3.5,0.5);
    boxT3->SetFillColor(kRed);
    if(e->get_EneT3().at(0)>0)boxT3->SetFillColor(kGreen);
    TPaveText*PHT3=new TPaveText(15,-1,25,3);
    PHT3->AddText(Form("PHT3=%d",(int)e->get_EneT3().at(0)));
    PHT3->SetFillStyle(0);PHT3->SetBorderSize(0);
    //Guard
    TBox*boxG1=new TBox(-13.5,-0.5588,-3.5,-0.5588+0.5);
    boxG1->SetFillColor(kRed);
    if(e->get_Eneg().at(0)>0)boxG1->SetFillColor(kGreen);
    TBox*boxG2=new TBox(+3.5,-0.5588,13.5,-0.5588+0.5);
    boxG2->SetFillColor(kRed);
    if(e->get_Eneg().at(0)>0)boxG2->SetFillColor(kGreen);     
    TPaveText*PHG=new TPaveText(15,-3,25,1);
    PHG->AddText(Form("PHG=%d",(int)e->get_Eneg().at(0)));
    PHG->SetFillStyle(0);PHG->SetBorderSize(0);
    
    //T4
    TBox*boxT4=new TBox(-18.,-25.59012,18.,-25.59012+1);
    boxT4->SetFillColor(kRed);
    if(e->get_EneT4().at(0)>0)boxT4->SetFillColor(kGreen);
    TPaveText*PHT4=new TPaveText(15,-25,25,-21);
    PHT4->AddText(Form("PHT4=%d",(int)e->get_EneT4().at(0)));
    PHT4->SetFillStyle(0);PHT4->SetBorderSize(0);
    
    //Non Bending plot   
    can->cd(1);
    
    TMultiGraph* multi=new TMultiGraph();multi->SetTitle(Form("Event: %d, Non bending plane",(int)e->get_eventnumber()));
    multi->Add(grALLNB,"p");
    multi->Add(grNB,"p");
    multi->Add(grT2,"");
    multi->Draw("a");
    multi->GetXaxis()->SetTitle("X (cm)");
    multi->GetYaxis()->SetTitle("Z (cm)");
    multi->GetXaxis()->SetLimits(-25,25);
    multi->SetMaximum(ZMax);
    multi->SetMinimum(ZMin);
    gPad->Update();
    
    TLine * Line0=new TLine(0,ZMin,0,ZMax);
    Line0->SetLineColor(kGray);
    Line0->SetLineStyle(3);
    Line0->SetLineWidth(1);
    Line0->Draw("same");
    
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
    if(fNB->GetParameter(1)>0) fNB->DrawF1(fNB->GetX(-30),fNB->GetX(40),"same");
    else fNB->DrawF1(fNB->GetX(40),fNB->GetX(-30),"same");
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
    ThetaNB->Draw();

    TPaveText*Tchi2NB;
    Tchi2NB=new TPaveText(15,44,25,50);
    Tchi2NB->SetFillStyle(0);Tchi2NB->SetBorderSize(0);
    Tchi2NB->SetTextColor(kBlue);
    Tchi2NB->AddText(Form("#chi^{2}= %6.4f",e->get_chi2NBPR()));
    Tchi2NB->Draw();
    
    TNh[0]->Draw();
    TNh[4]->Draw();
    TNh[6]->Draw();
    
    
    
    //Bending plot   
    can->cd(2);
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

    TMultiGraph *multiB=new TMultiGraph(); multiB->SetTitle("Bending plane");
    multiB->Add(grALLB,"p");
    multiB->Add(grB,"p");
    multiB->Add(grT2,"");
    multiB->Add(gr,"l");
    multiB->Draw("a");
    multiB->GetXaxis()->SetTitle("Y (cm)");
    multiB->GetYaxis()->SetTitle("Z (cm)");
    multiB->GetXaxis()->SetLimits(-25,25);
    multiB->SetMaximum(ZMax);
    multiB->SetMinimum(ZMin);
    gPad->Update();
    
    Line0->Draw("same");

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

    if(incomingB->GetParameter(1)>0) incomingB->DrawF1(incomingB->GetX(lim),incomingB->GetX(40),"same");
    else incomingB->DrawF1(incomingB->GetX(40),incomingB->GetX(lim),"same");
    
    if(incomingBout->GetParameter(1)>0) incomingBout->DrawF1(incomingBout->GetX(lim),incomingBout->GetX(-30),"same");
    else incomingBout->DrawF1(incomingBout->GetX(-30),incomingBout->GetX(lim),"same");
       
    if(outcomingB->GetParameter(1)<0) outcomingB->DrawF1(outcomingB->GetX(limo),outcomingB->GetX(-30),"same");
    else outcomingB->DrawF1(outcomingB->GetX(-30),outcomingB->GetX(limo),"same");

    TPaveText*ThetaIn;
    if(incomingB->GetX(40)+sh<25)ThetaIn=new TPaveText(incomingB->GetX(40),40,incomingB->GetX(40)+sh,44);
    else ThetaIn=new TPaveText(incomingB->GetX(40)-sh,40,incomingB->GetX(40),44);
    ThetaIn->SetFillStyle(0);ThetaIn->SetBorderSize(0);
    ThetaIn->AddText(Form("#theta_{Bin}= %6.4f",TMath::ATan(diff)));
    
    TPaveText*ThetaOut;
     if(outcomingB->GetX(-30)+sh<25) ThetaOut=new TPaveText(outcomingB->GetX(-30),-34,outcomingB->GetX(-30)+sh,-30);
     else  ThetaOut=new TPaveText(outcomingB->GetX(-30)-sh,-34,outcomingB->GetX(-30),-30);
    ThetaOut->SetFillStyle(0);ThetaOut->SetBorderSize(0);
    ThetaOut->AddText(Form("#theta_{Bout}= %6.4f",TMath::ATan(diffout)));
    ThetaIn->Draw();
    ThetaOut->Draw();
    
    TPaveText*Def=new TPaveText(-23,-38,10,-33);
    //Def->AddText(Form("Deflection: %6.4f",(float)deflection));
    float Def2=0;
    Def2=TMath::ATan(diffout)-TMath::ATan(diff);
    deflection=Def2;
    Def->AddText(Form("Deflection: #theta_{Bout}-#theta_{Bin}= %6.4f",(float)deflection));
    Def->SetFillStyle(0);Def->SetBorderSize(0);
    Def->Draw();
    TPaveText*Cur=new TPaveText(10,-23,25,-16);
    Cur->AddText(Form("Curvature at L1=%5.4f",curv[2]));
    Cur->AddText(Form("Curvature at L2=%5.4f",curv[1]));
    Cur->AddText(Form("Curvature at L3=%5.4f",curv[0]));
    Cur->SetFillStyle(0);Cur->SetBorderSize(0);
    //Cur->Draw();
    
    arc->Draw();
    TPaveText*RArc=new TPaveText(-20,45,-5,49);
    RArc->AddText(Form("R_{circ}= %6.4f cm",result.Parameter(2)));
    RArc->SetFillStyle(0);RArc->SetBorderSize(0);RArc->SetTextColor(kRed);
    RArc->Draw();   
    
    
    TPaveText*Tchi2B;
    Tchi2B=new TPaveText(15,44,25,50);
    Tchi2B->SetFillStyle(0);Tchi2B->SetBorderSize(0);
    Tchi2B->SetTextColor(kBlue);
    Tchi2B->AddText(Form("#chi^{2}= %6.4f",e->get_chi2BPR()));
    Tchi2B->Draw();

    TNh[1]->Draw();
    TNh[2]->Draw();
    TNh[3]->Draw();
    TNh[5]->Draw();
   
    
    //if(incomingB->GetX(33.25)<-14 ||incomingB->GetX(33.25)>14)continue;
    fileout->cd();
    can->Write();
   // if(firstpage==0)
   //  can->Print(Form("%s.EVENTDISPLAYCirc.pdf(",file.c_str()),Form("Title:Event %d",e->get_eventnumber()));
   // else can->Print(Form("%s.EVENTDISPLAYCirc.pdf",file.c_str()),Form("Title:Event %d",e->get_eventnumber()));
   // firstpage=-1;
    delete can;

   }//i
   
 fileout->Close();  
 //can=new TCanvas();
 //can->Print(Form("%s.EVENTDISPLAYCirc.pdf)",file.c_str()));

}//end function


void TrackinTrigger(string file,string coinc,int geoconf)
{
 
    
  gStyle->SetPalette(kBird);
   
    
 //Load configuration parameter
 float* TckZPos=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)TckZPos[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

 LoadDataparameters(paramfile,TckZPos,OffsetLL,OffsetRL,TrigThresh);

 for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", TckZPos:" << TckZPos[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 cout << "T1 threshold: " << TrigThresh[0] <<endl;
 cout << "T2 threshold: " << TrigThresh[1] <<endl;
 cout << "T3 threshold: " << TrigThresh[2] <<endl;
 cout << "T4 threshold: " << TrigThresh[3] <<endl;
 cout << "Guard threshold: " << TrigThresh[4] <<endl;
 

 
 
 //TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");
 
 int Nevents=0;
 
 //cout << Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/18A1_000to022.EVENT_PRonly.root") <<endl;
 //chain->Add(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/18A1_000to022.EVENT_PRonly.root"));
 chain->Add(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/FromSarah_10072018Calib/18A1_000to022.EVENT_FromSarah_10072018_PRonly_Calibrated.root"));

 //TFile*fileout=new TFile(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/FromSarah_10072018Calib/18A1_000to022.EVENT_FromSarah_10072018_PRonly_Calibrated_EVENTDISPLAY.root"),"RECREATE");
 //TCanvas *can=new TCanvas();
 
 int firstpage=0;
 ALEvent *e = new ALEvent();      
 //Define variables to read event
 //Set address to access event data
 chain->SetBranchAddress("Calevent",&e); 
 
 // Get number of event in Tree
 Nevents=chain->GetEntries();
 cout << "Number  of events: " << Nevents << endl;
 
 float ZMax=50;
 float ZMin=-40;
 
 //Timing variables
 double Overflowtime=TMath::Power(2,32);	
 double Overflowgo=TMath::Power(2,16);		
 double OverflowtimeTagPHA=0;	
 double OverflowgoTagPHA=0;		
 double OverflowtimeTagEVT=0;	
 double OverflowgoTagEVT=0;	
 double timeEVT=0;
 double prevtimeEVT=-100;
 double timePHA=0;
 double prevtimePHA=-100;
 double goEVT=0;
 double prevgoEVT=-100;
 double goPHA=0;
 double prevgoPHA=-100;	
 double firsttPHA=0; 	   
 double firstgoPHA=0; 	   
 double firsttEVT=0; 	   
 double firstgoEVT=0; 	   
	 
 for(int j=0;j<Nevents;j++)
   {
    chain->GetEntry(j);   
    if(e->get_GoEVT()<0)continue;
    firsttPHA=e->get_tPHA(); 	   
    firstgoPHA=e->get_GoPHA(); 	   
    firsttEVT=e->get_tEVT(); 	   
    firstgoEVT=e->get_GoEVT(); 	   
	j=Nevents;  
   }   

   
 int nPHA=6;
 TH1F** TrigX=new TH1F*[nPHA];
 TH1F** TrigY=new TH1F*[nPHA];
 TH2F** TrigXY=new TH2F*[nPHA];
 float* Ztrig=new float[nPHA]; 
 TH1F** PHAs=new TH1F*[nPHA];
 string TrigCros[6]={"T1","T2 top","T2 bottom","T3","Guard","T4"};
 Ztrig[0]=33.49968 ;
 Ztrig[1]=27.94 ;
 Ztrig[2]=2.024 ;
 Ztrig[3]=0. ;
 Ztrig[4]=-0.5588 ; 
 Ztrig[5]=-25.59012 ; 
 

 TEllipse** Cir= new TEllipse*[nPHA]; 
 float* Ra=new float[nPHA]; 
 Ra[0]=13.;
 Ra[1]=13.;
 Ra[2]=6.5;
 Ra[3]=3.5;
 Ra[4]=13.5;
 Ra[5]=18.;
 
   for(int i=0;i<nPHA;i++)
     {
      TrigX[i]=new TH1F(Form("X, %s",TrigCros[i].c_str()),Form("X, %s",TrigCros[i].c_str()),180,-30,30);
      TrigY[i]=new TH1F(Form("Y, %s",TrigCros[i].c_str()),Form("Y, %s",TrigCros[i].c_str()),180,-30,30);
      TrigXY[i]=new TH2F(Form("XY, %s",TrigCros[i].c_str()),Form("XY, %s",TrigCros[i].c_str()),180,-30,30,180,-30,30);
        
      PHAs[i]=new TH1F(Form("PHA, %s",TrigCros[i].c_str()),Form("PHA, %s",TrigCros[i].c_str()),5001,-1,5000);

      
      TrigX[i]->GetXaxis()->SetTitle("X(cm)");
      TrigY[i]->GetXaxis()->SetTitle("Y(cm)");
      TrigX[i]->GetYaxis()->SetTitle("Entries");
      TrigY[i]->GetYaxis()->SetTitle("Entries");
      TrigXY[i]->GetXaxis()->SetTitle("X(cm)");
      TrigXY[i]->GetYaxis()->SetTitle("Y(cm)") ;       
      PHAs[i]->GetXaxis()->SetTitle("PHA");
      PHAs[i]->GetYaxis()->SetTitle("Entries");
      
      
      Cir[i]=new TEllipse(0.,0., Ra[i],Ra[i]);   
      Cir[i]->SetLineColor(kBlack);
      Cir[i]->SetFillStyle(0);
      Cir[i]->SetLineWidth(2);
    }
  TH1F* p0s=new TH1F(Form("p0PR"),Form("p0PR"),4096,0,4096);
  p0s->GetXaxis()->SetTitle("Momentum p0PR");
  p0s->GetYaxis()->SetTitle("Entries");
  
    
 //Loop over events
 for(int i=0;i<Nevents;i++)
   {
    chain->GetEntry(i);
    if(i%100000==0) cout << "Event: " << i <<endl;
    //if(i<500000)continue;
   // if(i>505000)continue;
    
    int nnhits=e->get_Nhits();
    uint8_t Ti=(uint8_t)e->get_Ti();
    //cout << "Extract number of hits done: "<<nnhits  <<endl;
    double deflection=e->get_deflecPR();

    ///////////////////////////////    
    //Selection of events (Cuts)
    ///////////////////////////////
    //if(nnhits!=7)continue;
    //if(Ti!=127) continue;
    //cout << "Extract internal tracker trigger information done" <<endl;
    //if(abs(deflection)>2||abs(deflection)<0.5) continue;
    //if(abs(deflection)<0.7) continue;
    //if(e->get_EneT1().at(0)>200) continue;
    if(e->get_EneT1().at(0)!=-1) continue;
    //if(e->get_EneT2().at(0)<0) continue;
    //if(e->get_EneT3().at(0)<0) continue;
    //if(e->get_EneT3().at(0)>200) continue;
    //if(e->get_EneT4().at(0)<0) continue;
    //if(e->get_EneT4().at(0)>200) continue;
   //  if(e->get_Eneg().at(0)>0) continue;
  //  if(e->get_p0PR()<0) continue;
   // if(e->get_dPHA()>=19) continue;
    //if(e->get_PressB2()>1.6) continue;
  //  if(e->get_chi2NBPR()<=0||e->get_chi2NBPR()>2) continue;
 //   if(e->get_chi2BPR()<=0||e->get_chi2BPR()>2) continue;

    //Global timing	
	if(e->get_GoEVT()<prevgoEVT&& e->get_GoEVT()<0.01*Overflowgo&& prevgoEVT>0.9 *Overflowgo)//overflow go
	 {
	  OverflowgoTagEVT++;
	 }  
	prevgoEVT=	e->get_GoEVT();
	goEVT	=e->get_GoEVT()+OverflowgoTagEVT*Overflowgo;
		 
	if(e->get_GoPHA()<prevgoPHA&& e->get_GoPHA()<0.01*Overflowtime&& prevgoPHA>0.9 *Overflowtime)//overflow go PH A same as time 2^32
	 {
	  OverflowgoTagPHA++;
	  //cout << "Overflow Go PHA: event " <<  e->get_eventnumber()	 << endl; 
	 }  
	prevgoPHA=	e->get_GoPHA();
	goPHA	=e->get_GoPHA()+OverflowgoTagPHA*Overflowtime;
		 
	if(e->get_tEVT()<prevtimeEVT&& e->get_tEVT()<0.01*Overflowtime&& prevtimeEVT>0.9 *Overflowtime)//overflow time
	 {
	  OverflowtimeTagEVT++;
	  //cout << "Overflow time EVT: event " <<  e->get_eventnumber()	 << endl; 
	 }  
	prevtimeEVT=	e->get_tEVT();
	timeEVT	=e->get_tEVT()+OverflowtimeTagEVT*Overflowtime;
		 
	if(e->get_tPHA()<prevtimePHA && e->get_tPHA()<0.01*Overflowtime&& prevtimePHA>0.9 *Overflowtime)//overflow time
	 {
	  OverflowtimeTagPHA++;
	  //cout << "Overflow time PHA: event " <<  e->get_eventnumber()	 << endl; 
	 }  
	prevtimePHA=	e->get_tPHA();
	timePHA	=e->get_tPHA()+OverflowtimeTagPHA*Overflowtime;   
	double difftEVTtPHA=timeEVT-firsttEVT-timePHA+firsttPHA;
	//if(abs(difftEVTtPHA)>1000)continue;   
	   
    ///////////////////////////////    
    //Extract hits information
    ///////////////////////////////
   
    TGraph *grALLNB=new TGraph();
    grALLNB->SetMarkerStyle(5);
    grALLNB->SetMarkerColor(kBlack);
    TGraph *grALLB=new TGraph();
    grALLB->SetMarkerStyle(5);
    grALLB->SetMarkerColor(kBlack);
    TGraph *grNB=new TGraph();
    grNB->SetMarkerStyle(kCircle);
    grNB->SetMarkerColor(kBlue);
    TGraph *grB=new TGraph();
    grB->SetMarkerStyle(kCircle);
    grB->SetMarkerColor(kBlue);
    //cout << "Graphs are defined" <<endl;

    for(int j=0;j<nnhits;j++)
      {
       int L= e->get_hits().at(j)->get_L();
       float x=e->get_hits().at(j)->get_x();
       float y=e->get_hits().at(j)->get_y();
       float z=e->get_hits().at(j)->get_z();
       bool flagPR=e->get_hits().at(j)->get_flagPR();
       
       if(L==0 ||L==4||L==6)//Non bending plane
        {
         grALLNB->SetPoint(grALLNB->GetN(),x,z);
         if(flagPR)grNB->SetPoint(grNB->GetN(),x,z);
        }
       else if(L==1 ||L==2||L==3 ||L==5) //Bending plane
        {
         grALLB->SetPoint(grALLB->GetN(),y,z);
         if(flagPR)grB->SetPoint(grB->GetN(),y,z);
        }
      }  
  
    ///////////////////////////////    
    //Extract Reconstructed track information
    ///////////////////////////////    
    
    TF1* fNB=new TF1("fNB","pol1",-400,400);
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

    //Non bending plane
    float p0=e->get_interPR();
    float p1=e->get_slopePR();
    if(p1==0) continue;
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

    //cout << "Extraction of bending plane fit parameters done" <<endl;

    //Incoming Straight particle  
    float lim=TckZPos[1];//z position of 2nd layer
    TF1*incomingB=new TF1("incomingB","pol1",-400,400);
    incomingB->SetLineStyle(2);
    incomingB->SetLineWidth(2);
    incomingB->SetLineColor(kBlue);
    float aa=fB->Eval(lim);
    float diff=2*a*lim+2*a*zz0+b;
    incomingB->FixParameter(0,-(aa-diff*lim)/diff);
    incomingB->FixParameter(1,1./diff);
    
    TF1*incomingBout=new TF1("incomingBout","pol1",-400,400);
    incomingBout->SetLineStyle(2);
    incomingBout->SetLineWidth(2);
    incomingBout->SetLineColor(kGray);
    incomingBout->FixParameter(0,-(aa-diff*lim)/diff);
    incomingBout->FixParameter(1,1./diff);    
    
    //Outcoming Straight particle    
    float limo=TckZPos[5];//z position of 6th layer
    TF1*outcomingB=new TF1("outcomingB","pol1",-400,400);
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
    
    
    ///////////////////////////////    
    //Fill 2d histograms of crossing point of track in Triggers 
    ///////////////////////////////   
    
    for(int i=0;i<nPHA;i++)
      {
       TrigX[i]->Fill(fNB->GetX(Ztrig[i]));
       
       if(Ztrig[i]>TckZPos[6])//Need to use incoming line
         {
          TrigY[i]->Fill(incomingB->GetX(Ztrig[i]));
          TrigXY[i]->Fill(fNB->GetX(Ztrig[i]),incomingB->GetX(Ztrig[i]));
         } 
       else //Need to use outcoming line
        {
         TrigY[i]->Fill(outcomingB->GetX(Ztrig[i]));
         TrigXY[i]->Fill(fNB->GetX(Ztrig[i]),outcomingB->GetX(Ztrig[i]));
        }  
    }
    
//     if(fNB->GetX(Ztrig[0]) <13  &&fNB->GetX(Ztrig[0]) >8 && incomingB->GetX(Ztrig[0])>-20&& incomingB->GetX(Ztrig[0])<-15)
//      {
//        PHAs[0]->Fill(e->get_EneT1().at(0)) ;    
//        PHAs[1]->Fill(e->get_EneT2().at(0)) ;    
//        PHAs[2]->Fill(e->get_EneT2().at(0)) ;    
//        PHAs[3]->Fill(e->get_EneT3().at(0)) ;    
//        PHAs[4]->Fill(e->get_Eneg().at(0)) ;    
//        PHAs[5]->Fill(e->get_EneT4().at(0)) ;    
//        p0s->Fill(e->get_p0PR());  
//      }
    PHAs[0]->Fill(e->get_EneT1().at(0)) ;    
    PHAs[1]->Fill(e->get_EneT2().at(0)) ;    
    PHAs[2]->Fill(e->get_EneT2().at(0)) ;    
    PHAs[3]->Fill(e->get_EneT3().at(0)) ;    
    PHAs[4]->Fill(e->get_Eneg().at(0)) ;    
    PHAs[5]->Fill(e->get_EneT4().at(0)) ;    
    p0s->Fill(e->get_p0PR());     
    ///////////////////////////////    
    //Event Display 
    ///////////////////////////////
//     can=new TCanvas(Form("Event %d",e->get_eventnumber()),Form("Event %d",e->get_eventnumber()),200,10,1600,1440);
//     can->Divide(2,1);
//     //Detector layout
//     //Layers
//     TLine**Line=new TLine*[7];
// 
//     for(int ijk=0;ijk<7;ijk++)
//       {
//        Line[ijk]=new TLine(-9,TckZPos[ijk],9,TckZPos[ijk]);
//        Line[ijk]->SetLineColor(kBlack);  
//        Line[ijk]->SetLineWidth(1);  
//       }
//        
//     //Magnets 
//     TBox*boxM1=new TBox(-15,-9.61682,-6.703,-9.61682+4.59994);
//     TBox*boxM2=new TBox(-15,-16.2513,-6.703,-16.2513+4.59994);
//     TBox*boxM3=new TBox(6.703,-9.61682,15,-9.61682+4.59994);
//     TBox*boxM4=new TBox(6.703,-16.2513,15,-16.2513+4.59994);
//     boxM1->SetFillColor(kGray);
//     boxM2->SetFillColor(kGray);
//     boxM3->SetFillColor(kGray);
//     boxM4->SetFillColor(kGray);
//      
//     //T1
//     TBox*boxT1=new TBox(-13.,33.49968,13.,33.49968+0.5);
//     boxT1->SetFillColor(kRed);
//     if(e->get_EneT1().at(0)>0)boxT1->SetFillColor(kGreen);
//     TPaveText*PHT1=new TPaveText(15,32,25,36);
//     PHT1->AddText(Form("PHT1=%d",(int)e->get_EneT1().at(0)));
//     PHT1->SetFillStyle(0);PHT1->SetBorderSize(0);
//     //T2
//     TGraph* grT2 = new TGraph();
//     grT2->SetPoint(0,-6.5,2.024);//Bottom
//     grT2->SetPoint(1,6.5,2.024);//Botton
//     grT2->SetPoint(2,13.5,29.964);//Top
//     grT2->SetPoint(3,-13.5,29.964);//Top
//     grT2->SetPoint(4,-6.5,2.024);//Bottom again
//     grT2->SetFillStyle(0);
//     grT2->SetLineWidth(3);
//     grT2->SetLineColor(kRed);
//     if(e->get_EneT2().at(0)>0)grT2->SetLineColor(kGreen); 
//     TPaveText*PHT2=new TPaveText(15,11,25,15);
//     PHT2->AddText(Form("PHT2=%d",(int)e->get_EneT2().at(0)));
//     PHT2->SetFillStyle(0);PHT2->SetBorderSize(0);     
//     //T3
//     TBox*boxT3=new TBox(-3.5,0.,3.5,0.5);
//     boxT3->SetFillColor(kRed);
//     if(e->get_EneT3().at(0)>0)boxT3->SetFillColor(kGreen);
//     TPaveText*PHT3=new TPaveText(15,-1,25,3);
//     PHT3->AddText(Form("PHT3=%d",(int)e->get_EneT3().at(0)));
//     PHT3->SetFillStyle(0);PHT3->SetBorderSize(0);
//     //Guard
//     TBox*boxG1=new TBox(-13.5,-0.5588,-3.5,-0.5588+0.5);
//     boxG1->SetFillColor(kRed);
//     if(e->get_Eneg().at(0)>0)boxG1->SetFillColor(kGreen);
//     TBox*boxG2=new TBox(+3.5,-0.5588,13.5,-0.5588+0.5);
//     boxG2->SetFillColor(kRed);
//     if(e->get_Eneg().at(0)>0)boxG2->SetFillColor(kGreen);     
//     TPaveText*PHG=new TPaveText(15,-3,25,1);
//     PHG->AddText(Form("PHG=%d",(int)e->get_Eneg().at(0)));
//     PHG->SetFillStyle(0);PHG->SetBorderSize(0);
//     
//     //T4
//     TBox*boxT4=new TBox(-18.,-25.59012,18.,-25.59012+1);
//     boxT4->SetFillColor(kRed);
//     if(e->get_EneT4().at(0)>0)boxT4->SetFillColor(kGreen);
//     TPaveText*PHT4=new TPaveText(15,-25,25,-21);
//     PHT4->AddText(Form("PHT4=%d",(int)e->get_EneT4().at(0)));
//     PHT4->SetFillStyle(0);PHT4->SetBorderSize(0);
//     
//     //Non Bending plot   
//     can->cd(1);
//     
//     TMultiGraph* multi=new TMultiGraph();multi->SetTitle(Form("Event: %d, Non bending plane",i));
//     multi->Add(grALLNB,"p");
//     multi->Add(grNB,"p");
//     multi->Add(grT2,"");
//     multi->Draw("a");
//     multi->GetXaxis()->SetTitle("X (cm)");
//     multi->GetYaxis()->SetTitle("Z (cm)");
//     multi->GetYaxis()->SetTitleOffset(1.1);
//     multi->GetXaxis()->SetLimits(-25,25);
//     multi->SetMaximum(ZMax);
//     multi->SetMinimum(ZMin);
//     gPad->Update();
//     
//     TLine * Line0=new TLine(0,ZMin,0,ZMax);
//     Line0->SetLineColor(kGray);
//     Line0->SetLineStyle(3);
//     Line0->SetLineWidth(1);
//     Line0->Draw("same");
//     
//     //T1
//     boxT1->Draw("l same");
//     PHT1->Draw();
//     //T2
//     PHT2->Draw();
//     //T3
//     boxT3->SetLineWidth(1);
//     boxT3->SetLineColor(1);
//     boxT3->Draw("l same");
//     PHT3->Draw();
//     //Guard
//     boxG1->SetLineWidth(1);
//     boxG2->SetLineWidth(1);
//     boxG1->SetLineColor(1);
//     boxG2->SetLineColor(1);
// 
//     boxG1->Draw("l same");
//     boxG2->Draw("l same");
//     PHG->Draw();
//     //T4
//     boxT4->Draw("l same");
//     PHT4->Draw();
//     //Layers
//     for(int ijk=0;ijk<7;ijk++) Line[ijk]->Draw("same");  
//     //Non beding stright line
//     if(fNB->GetParameter(1)>0) fNB->DrawF1(fNB->GetX(-30),fNB->GetX(40),"same");
//     else fNB->DrawF1(fNB->GetX(40),fNB->GetX(-30),"same");
//     //Magnet
//     boxM1->Draw("l same");
//     boxM2->Draw("l same");
//     boxM3->Draw("l same");
//     boxM4->Draw("l same");       
// 
//     TPaveText*ThetaNB;
//     int sh=15;
//     if(fNB->GetX(40)+sh<25)ThetaNB=new TPaveText(fNB->GetX(40),40,fNB->GetX(40)+sh,44);
//     else ThetaNB=new TPaveText(fNB->GetX(40)-sh,40,fNB->GetX(40),44);
//     ThetaNB->SetFillStyle(0);ThetaNB->SetBorderSize(0);
//     ThetaNB->AddText(Form("#theta_{NB}= %6.4f",TMath::ATan(p1)));
//     ThetaNB->Draw();
//     
//     //Bending plot   
//     can->cd(2);
//     //Triclk to the Inverse bending track for display
//     int Nn=500;
//     Double_t*x=new Double_t[Nn]; 
//     Double_t*y=new Double_t[Nn];
//     Double_t dx=(lim-limo)/Nn;
//     Double_t x1=limo;
//     x[0] = x1;
//     y[0] = fB->Eval(x[0]);
//     for (int ij=1; ij<Nn; ij++) 
//       {
//        x1   = x1+dx;
//        x[ij] = x1;
//        y[ij] = fB->Eval(x[ij]);
//       }
// 
//     TGraph *gr = new TGraph(Nn,y,x);
//     gr->SetLineColor(kBlue);
//     gr->SetLineWidth(2);
// 
//     TMultiGraph *multiB=new TMultiGraph(); multiB->SetTitle("Bending plane");
//     multiB->Add(grALLB,"p");
//     multiB->Add(grB,"p");
//     multiB->Add(grT2,"");
//     multiB->Add(gr,"l");
//     multiB->Draw("a");
//     multiB->GetXaxis()->SetTitle("Y (cm)");
//     multiB->GetYaxis()->SetTitle("Z (cm)");
//     multiB->GetYaxis()->SetTitleOffset(1.1);
//     multiB->GetXaxis()->SetLimits(-25,25);
//     multiB->SetMaximum(ZMax);
//     multiB->SetMinimum(ZMin);
//     gPad->Update();
//     
//     Line0->Draw("same");
// 
//     //Magnet
//     boxM1->Draw("l same");
//     boxM2->Draw("l same");
//     boxM3->Draw("l same");
//     boxM4->Draw("l same");       
// 
//     //T1
//     boxT1->Draw("l same");
//     //T2
//     //T3
//     boxT3->SetLineWidth(1);
//     boxT3->SetLineColor(1);
//     boxT3->Draw("l same");
//     //Guard
//     boxG1->SetLineWidth(1);
//     boxG1->SetLineColor(1);
//     boxG2->SetLineWidth(1);
//     boxG2->SetLineColor(1);
//     boxG1->Draw("l same");
//     boxG2->Draw("l same");
//     //T4
//     boxT4->Draw("l same");
//     //Layers
//     for(int ijk=0;ijk<7;ijk++) Line[ijk]->Draw("same");      
// 
//     if(incomingB->GetParameter(1)>0) incomingB->DrawF1(incomingB->GetX(lim),incomingB->GetX(40),"same");
//     else incomingB->DrawF1(incomingB->GetX(40),incomingB->GetX(lim),"same");
//     
//     if(incomingBout->GetParameter(1)>0) incomingBout->DrawF1(incomingBout->GetX(lim),incomingBout->GetX(-30),"same");
//     else incomingBout->DrawF1(incomingBout->GetX(-30),incomingBout->GetX(lim),"same");
//        
//     if(outcomingB->GetParameter(1)<0) outcomingB->DrawF1(outcomingB->GetX(limo),outcomingB->GetX(-30),"same");
//     else outcomingB->DrawF1(outcomingB->GetX(-30),outcomingB->GetX(limo),"same");
// 
//     TPaveText*ThetaIn;
//     if(incomingB->GetX(40)+sh<25)ThetaIn=new TPaveText(incomingB->GetX(40),40,incomingB->GetX(40)+sh,44);
//     else ThetaIn=new TPaveText(incomingB->GetX(40)-sh,40,incomingB->GetX(40),44);
//     ThetaIn->SetFillStyle(0);ThetaIn->SetBorderSize(0);
//     ThetaIn->AddText(Form("#theta_{Bin}= %6.3f",TMath::ATan(diff)));
//     
//     TPaveText*ThetaOut;
//      if(outcomingB->GetX(-30)+sh<25) ThetaOut=new TPaveText(outcomingB->GetX(-30),-34,outcomingB->GetX(-30)+sh,-30);
//      else  ThetaOut=new TPaveText(outcomingB->GetX(-30)-sh,-34,outcomingB->GetX(-30),-30);
//     ThetaOut->SetFillStyle(0);ThetaOut->SetBorderSize(0);
//     ThetaOut->AddText(Form("#theta_{Bout}= %6.3f",TMath::ATan(diffout)));
//     ThetaIn->Draw();
//     ThetaOut->Draw();
//     
//     TPaveText*Def=new TPaveText(-23,-38,10,-33);
//     //Def->AddText(Form("Deflection: %6.4f",(float)deflection));
//     float Def2=0;
//     Def2=TMath::ATan(diffout)-TMath::ATan(diff);
//     deflection=Def2;
//     Def->AddText(Form("Deflection: #theta_{Bout}-#theta_{Bin}= %6.3f",(float)deflection));
//     Def->SetFillStyle(0);Def->SetBorderSize(0);
//     Def->Draw();
//     TPaveText*Cur=new TPaveText(10,-23,25,-16);
//     Cur->AddText(Form("Curvature at L1=%5.4f",curv[2]));
//     Cur->AddText(Form("Curvature at L2=%5.4f",curv[1]));
//     Cur->AddText(Form("Curvature at L3=%5.4f",curv[0]));
//     Cur->SetFillStyle(0);Cur->SetBorderSize(0);
//     
//     
//     TPaveText*Alt=new TPaveText(-23,45,10,50);
//     Alt->AddText(Form("Altitude: %4.2f mmHg",e->get_PressB2()));
//     Alt->SetFillStyle(0);Alt->SetBorderSize(0);
//     Alt->Draw();
//     
//     //Cur->Draw();
//     
//     //if(incomingB->GetX(33.25)<-14 ||incomingB->GetX(33.25)>14)continue;
//    // fileout->cd();
//    // can->Write();
//     
//   //  if(firstpage==0)
//   //   can->Print(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/FromSarah_10072018Calib/18A1_000to022.EVENT_PRonly_Calibrated.root.FromSarah_10072018.EVENTDISPLAY.pdf("),Form("Title:Event %d",i));
//   //  else can->Print(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/FromSarah_10072018Calib/18A1_000to022.EVENT_PRonly_Calibrated.root.FromSarah_10072018.EVENTDISPLAY.pdf"),Form("Title:Event %d",i));
//     
// 
// 	firstpage=-1;
//     delete can;

   }//i
   
 //fileout->Close();  
// can=new TCanvas();
// can->Print(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/FromSarah_10072018Calib/18A1_000to022.EVENT_PRonly_Calibrated.root.FromSarah_10072018.EVENTDISPLAY.pdf)"));

TCanvas *canX=new TCanvas("canX","canX",200,10,1500,1000);
canX->Divide(3,2);

for(int i=0;i<nPHA;i++)
  {
   canX->cd(i+1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gPad->SetTopMargin(0.15);
   gPad->SetBottomMargin(0.15);
   TrigX[i]->Draw();
  }

TCanvas *canY=new TCanvas("canY","canY",200,10,1500,1000);
canY->Divide(3,2);

for(int i=0;i<nPHA;i++)
  {
   canY->cd(i+1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gPad->SetTopMargin(0.15);
   gPad->SetBottomMargin(0.15);
   TrigY[i]->Draw();
  }

TCanvas *canXY=new TCanvas("canXY","canXY",200,10,1500,1000);
canXY->Divide(3,2);




for(int i=0;i<nPHA;i++)
  {
   canXY->cd(i+1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gPad->SetTopMargin(0.15);
   gPad->SetBottomMargin(0.15);
   TrigXY[i]->Draw("colz1");
   Cir[i]->Draw("same");
   if(i==4) Cir[i-1]->Draw("same");
  }
  
  
TCanvas *canP=new TCanvas("canP","canP",200,10,1500,1000);
canP->Divide(3,2);

for(int i=0;i<nPHA;i++)
  {
   canP->cd(i+1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gPad->SetTopMargin(0.15);
   gPad->SetBottomMargin(0.15);
   PHAs[i]->Draw();
  }

  
TCanvas *canPC=new TCanvas("canPC","canPC",200,10,1500,1000);
canPC->Divide(3,2);

for(int i=0;i<nPHA;i++)
  {
   canPC->cd(i+1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gPad->SetTopMargin(0.15);
   gPad->SetBottomMargin(0.15);
   PHAs[i]->GetCumulative()->Draw();
  }  
  
  
}//end function





void fitCircle(Int_t ) ;
void fitCircle(Int_t n=10000) {
   //generates n points around a circle and fit them
   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->SetGrid();
   TGraph* gr = new TGraph(n);
   if (n> 999) gr->SetMarkerStyle(1);
   else        gr->SetMarkerStyle(3);
   TRandom3 r;
   Double_t x,y;
   for (Int_t i=0;i<n;i++) {
      r.Circle(x,y,r.Gaus(4,0.3));
      gr->SetPoint(i,x,y);
   }
   c1->DrawFrame(-5,-5,5,5);
   gr->Draw("p");
   auto chi2Function = [&](const Double_t *par) {
      //minimisation function computing the sum of squares of residuals
      // looping at the graph points
      Int_t np = gr->GetN();
      Double_t f = 0;
      Double_t *x = gr->GetX();
      Double_t *y = gr->GetY();
      for (Int_t i=0;i<np;i++) {
         Double_t u = x[i] - par[0];
         Double_t v = y[i] - par[1];
         Double_t dr = par[2] - std::sqrt(u*u+v*v);
         f += dr*dr;
      }
      return f;
   };
   // wrap chi2 funciton in a function object for the fit
   // 3 is the number of fit parameters (size of array par)
   ROOT::Math::Functor fcn(chi2Function,3);
   ROOT::Fit::Fitter  fitter;
   double pStart[3] = {0,0,1};
   fitter.SetFCN(fcn, pStart);
   fitter.Config().ParSettings(0).SetName("x0");
   fitter.Config().ParSettings(1).SetName("y0");
   fitter.Config().ParSettings(2).SetName("R");
   // do the fit 
   bool ok = fitter.FitFCN();
   if (!ok) {
      Error("line3Dfit","Line3D Fit failed");
   }   
   const ROOT::Fit::FitResult & result = fitter.Result();
   result.Print(std::cout);
   //Draw the circle on top of the points
   TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
   arc->SetLineColor(kRed);
   arc->SetLineWidth(4);
   arc->SetFillStyle(0);
   arc->Draw();
}

