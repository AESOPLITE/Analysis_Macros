////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 14 , 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "TChain.h"
#include "TPaveStats.h"
#include "TGaxis.h"


ClassImp(ALTckhit)
ClassImp(ALEvent)

void AnalyseDataEvent(string s,string c,int geoconf);
void PullDistribution(string s,string c,int geoconf);
void Occupancy(string s,string c);
void Pulses(string s,string c);
void RealRates();
void RecoEff();
void RecoEffOne(string s,string c);
void Deflection(string s,string c, int);
void Timing(string s,string c);
void NoiseInvestigation(string s,string c,int geoconf);
void PalestineChip(string s,string c,int geoconf);
void PalestinePulsePHADisc();
void PalestineLOGICDisc();
void DeflectionT1T2T3(string s,string c,int geoconf);

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


void AnalyseDataEvent(string s,string c,int geoconf)
{
  //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

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

 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;   

//////////////////////////// 
//Occupancy histograms
//////////////////////////// 
    
 TH1F* LOcc=new TH1F("Layer occupancy","Layer occupancy",10,0,10);
 LOcc->GetXaxis()->SetTitle("Layer (0 to 6)");
 LOcc->GetYaxis()->SetTitle("event number");
 LOcc->SetLineColor(kRed);
 LOcc->SetLineWidth(1);  
 LOcc->SetMinimum(0);  

 TH1F* ClusOcc=new TH1F("Cluster occupancy","Cluster occupancy",10,0,10);
 ClusOcc->GetXaxis()->SetTitle("Layer (0 to 6)");
 ClusOcc->GetYaxis()->SetTitle("cluster number");
 ClusOcc->SetLineColor(kRed);
 ClusOcc->SetLineWidth(1);  
 ClusOcc->SetMinimum(0);  

 
 TH1F** ChipOcc=new TH1F*[7];
 TH1F** fStripOcc=new TH1F*[7];
 TH1F** CoordOcc=new TH1F*[7];
 
 for(int i=0;i<7;i++)
   {
    ChipOcc[i]=new TH1F(Form("Chip ID L%d",i),Form("Chip ID L%d",i),15,0,15);
    ChipOcc[i]->GetXaxis()->SetTitle("Chip ID");
    ChipOcc[i]->GetYaxis()->SetTitle("Entries");
    ChipOcc[i]->SetLineColor(kRed);
    ChipOcc[i]->SetLineWidth(1);
    ChipOcc[i]->GetXaxis()->SetRangeUser(0,15);  
    ChipOcc[i]->SetMinimum(0);  
    fStripOcc[i]=new TH1F(Form("First strip ID L%d",i),Form("First strip ID L%d",i),769,0,769);
    fStripOcc[i]->GetXaxis()->SetTitle("First strip ID");
    fStripOcc[i]->GetYaxis()->SetTitle("Entries");
    fStripOcc[i]->SetLineColor(kRed);
    fStripOcc[i]->SetLineWidth(1);
    fStripOcc[i]->GetXaxis()->SetRangeUser(0,769);  
    fStripOcc[i]->SetMinimum(0);  
    CoordOcc[i]=new TH1F(Form("Cluster coord. L%d",i),Form("Cluster coord. L%d",i),200,-100,100);
    if(i==0||i==4||i==6)CoordOcc[i]->GetXaxis()->SetTitle("X (mm)");
    else CoordOcc[i]->GetXaxis()->SetTitle("Y (mm)");
    CoordOcc[i]->GetYaxis()->SetTitle("Entries");
    CoordOcc[i]->SetLineColor(kRed);
    CoordOcc[i]->SetLineWidth(1);
    CoordOcc[i]->GetXaxis()->SetRangeUser(-100,100);    
    CoordOcc[i]->SetMinimum(0);      
   }//i   
    
    
//////////////////////////// 
//Trigger related histograms
//////////////////////////// 
 
 TH1F**AllPH=new TH1F*[5];
 string sPH[5]={"T1","T2","T3","T4","Guard"};
 string XtitlePH="Pulse height";
 string YtitlePH="Entries";
 
 for(int i=0;i<5;i++)
   {
    AllPH[i]=new TH1F(Form("Pulse heights %s",sPH[i].c_str()),Form("Pulse heights %s",sPH[i].c_str()),250,-1,500);
    AllPH[i]->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
    AllPH[i]->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
    AllPH[i]->SetLineColor(kRed);
    AllPH[i]->SetLineWidth(3);
    AllPH[i]->GetXaxis()->SetRangeUser(-1,500);
   } //i  

 TH1F*T3woG=new TH1F("Pulse heights T3, without Guard","Pulse heights T3, no Guard",250,-1,500);
 TH1F*T3wG=new TH1F("Pulse heights T3, with Guard","Pulse heights T3, with Guard",250,-1,500); 
 TH1F*GwoT3=new TH1F("Pulse heights Guard, without T3","Pulse heights Guard, without T3",250,-1,500);
 TH1F*GwT3=new TH1F("Pulse heights Guard, with T3","Pulse heights Guard, with T3",250,-1,500); 

 T3woG->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
 T3woG->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
 T3woG->SetLineColor(kRed);T3woG->SetLineWidth(3);
 T3woG->GetXaxis()->SetRangeUser(-1,500);
 T3wG->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
 T3wG->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
 T3wG->SetLineColor(kRed);T3wG->SetLineWidth(3);
 T3wG->GetXaxis()->SetRangeUser(-1,500);
 GwoT3->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
 GwoT3->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
 GwoT3->SetLineColor(kRed);GwoT3->SetLineWidth(3);
 GwoT3->GetXaxis()->SetRangeUser(-1,500);
 GwT3->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
 GwT3->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
 GwT3->SetLineColor(kRed);GwT3->SetLineWidth(3);
 GwT3->GetXaxis()->SetRangeUser(-1,500);

//////////////////////////// 
//Deflection histograms
////////////////////////////   
 TH1F*hDef=new TH1F("Deflection  T1&T3&T4&NoT2","Deflection  T1&T3&T4&NoT2",80,-0.2,0.2); 
 hDef->GetXaxis()->SetTitle("#DeltaSlope=S_{out}-S_{in}");
 hDef->GetYaxis()->SetTitle("entries");
 hDef->SetLineColor(kRed);hDef->SetLineWidth(2);
  
 
 TH2F*DefvsT1T3T4=new TH2F("Deflection vs  T1&T3&T4","Deflection vs T1&T3&T4",100,0,1000,200,-0.5,0.5); 
 DefvsT1T3T4->GetYaxis()->SetTitle("#DeltaSlope=S_{out}-S_{in}");
 DefvsT1T3T4->GetXaxis()->SetTitle("Sum of pulse heights: T1+T3+T4");

 TH1F*hNBslope=new TH1F("NB slope T1&T3&T4","NB slope T1&T3&T4",100,-0.5,0.5); 
 hNBslope->GetXaxis()->SetTitle("Non bending plane: slope of x=slope#timesz+intersection");
 hNBslope->GetYaxis()->SetTitle("entries");
 hNBslope->SetLineColor(kRed);hNBslope->SetLineWidth(2);
  
 TH1F*hNBinter=new TH1F("NB intersection T1&T3&T4","NB intersection T1&T3&T4",100,-10,10); 
 hNBinter->GetXaxis()->SetTitle("Non bending plane: intersection of x=slope#timesz+intersection");
 hNBinter->GetYaxis()->SetTitle("entries");
 hNBinter->SetLineColor(kRed);hNBinter->SetLineWidth(2);  
  
 TH1F*hNBXT1=new TH1F("X extrapolated in T1, T1&T3&T4","X extrapolated in T1, T1&T3&T4",120,-15,15); 
 hNBXT1->GetXaxis()->SetTitle("Non bending plane: X (cm) extrapolated in T1");
 hNBXT1->GetYaxis()->SetTitle("entries");
 hNBXT1->SetLineColor(kRed);hNBXT1->SetLineWidth(2);   

 
 TH1F*hNBXLEET4=new TH1F("X extrapolated in old T4, T1&T4","X extrapolated in old T4, T1&T4",120,-15,15); 
 hNBXLEET4->GetXaxis()->SetTitle("Non bending plane: X (cm) extrapolated in LEE T4");
 hNBXLEET4->GetYaxis()->SetTitle("entries");
 hNBXLEET4->SetLineColor(kRed);hNBXLEET4->SetLineWidth(2);   
 
 TH1F*hNBXT4=new TH1F("X extrapolated in T4, T1&T4","X extrapolated in T4, T1&T4",120,-15,15); 
 hNBXT4->GetXaxis()->SetTitle("Non bending plane: X (cm) extrapolated in  T4");
 hNBXT4->GetYaxis()->SetTitle("entries");
 hNBXT4->SetLineColor(kRed);hNBXT4->SetLineWidth(2);   

 TH1F*hNBYT4=new TH1F("Y extrapolated in T4, T1&T4","Y extrapolated in T4, T1&T4",120,-15,15); 
 hNBYT4->GetXaxis()->SetTitle("Non bending plane: Y (cm) extrapolated in old T4");
 hNBYT4->GetYaxis()->SetTitle("entries");
 hNBYT4->SetLineColor(kRed);hNBYT4->SetLineWidth(2);   

 TH1F*hNBYLEET4=new TH1F("Y extrapolated in old T4, T1&T4","Y extrapolated in old T4, T1&T4",120,-15,15); 
 hNBYLEET4->GetXaxis()->SetTitle("Non bending plane: Y (cm) extrapolated in LEE T4");
 hNBYLEET4->GetYaxis()->SetTitle("entries");
 hNBYLEET4->SetLineColor(kRed);hNBYLEET4->SetLineWidth(2); 
//////////////////////////// 
//Resolution histograms
//////////////////////////// 

 TH1F**XReso=new TH1F*[7];
 TH1F**YReso=new TH1F*[7];
 TH1F**ZReso=new TH1F*[7];
 TH2F**XResoX=new TH2F*[7];
 TH2F**YResoY=new TH2F*[7];
 TH2F**ZResoXY=new TH2F*[7];
 
 for(int j=0;j<7;j++)
   {
    XReso[j]=new TH1F(Form("Xreso L%d",j),Form("Xreso L%d",j),200,-1,1); 
    YReso[j]=new TH1F(Form("Yreso L%d",j),Form("Yreso L%d",j),200,-1,1);
    ZReso[j]=new TH1F(Form("Zreso L%d",j),Form("Zreso L%d",j),200,-3,3);
    XResoX[j]=new TH2F(Form("Xreso vs X L%d",j),Form("Xreso vs X L%d",j),1000,-100,100,400,-1,1); 
    YResoY[j]=new TH2F(Form("Yreso vs Y L%d",j),Form("Yreso vs Y L%d",j),1000,-100,100,400,-1,1);
    if(j==0||j==4||j==6)ZResoXY[j]=new TH2F(Form("Zreso vs X L%d",j),Form("Zreso vs X L%d",j),1000,-100,100,400,-2,2);
    else ZResoXY[j]=new TH2F(Form("Zreso vs Y L%d",j),Form("Zreso vs Y L%d",j),1000,-100,100,400,-2,2);
   }
    
    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");

 chain->Add(Form("../Data/%s.EVENT.root",file[0].c_str()));
    
 //Define variables to read event
 ALEvent *e = new ALEvent();      
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 int nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;  
      
 for (int j=0;j<nentries;j++)
   {
    chain->GetEntry(j); //Load the entry i in the variable e 

    if(j%10000==0) cout << "Event: " << j <<endl;
    
    int nnhits = (int) e->get_Nhits();

    ////////////////////////////////////
    //TRIGGER INFO
    ////////////////////////////////////  
    
    AllPH[0]->Fill(e->get_EneT1().at(0));
    AllPH[1]->Fill(e->get_EneT2().at(0));
    AllPH[2]->Fill(e->get_EneT3().at(0));
    AllPH[3]->Fill(e->get_EneT4().at(0));
    AllPH[4]->Fill(e->get_Eneg().at(0));
    
    if(e->get_EneT3().at(0)>0)GwT3->Fill(e->get_Eneg().at(0));
    else GwoT3->Fill(e->get_Eneg().at(0));
    if(e->get_Eneg().at(0)>0)T3wG->Fill(e->get_EneT3().at(0));
    else T3woG->Fill(e->get_EneT3().at(0));
    
    if(e->get_EneT1().at(0)>0 && e->get_EneT3().at(0)>0 && e->get_EneT2().at(0)<0 &&e->get_EneT4().at(0)>0&&e->get_deflecPR()!=0) hDef->Fill(e->get_deflecPR());
    if(e->get_EneT1().at(0)>0 && e->get_EneT3().at(0)>0 && e->get_EneT2().at(0)<0 &&e->get_EneT4().at(0)>0&&e->get_deflecPR()!=0) 
     {  
      DefvsT1T3T4->Fill(e->get_EneT1().at(0)+e->get_EneT3().at(0)+e->get_EneT4().at(0),e->get_deflecPR());
     }
    
     if(e->get_EneT1().at(0)>0 && e->get_EneT3().at(0)>0 &&e->get_EneT4().at(0)>0&&e->get_slopePR()!=0) 
     {  
      hNBslope->Fill(e->get_slopePR());
      hNBinter->Fill(e->get_interPR());
      hNBXT1->Fill(e->get_interPR()+e->get_slopePR()*33.5);
     }
     if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&e->get_PHA6().at(0)>0&&e->get_slopePR()!=0&e->get_chi2NBPR()<10) 
     {  
      hNBXLEET4->Fill(e->get_interPR()+e->get_slopePR()*(-25.59012-2.54*1.5));//Middle of old LEE T4 
     }
     if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&e->get_slopePR()!=0&e->get_chi2NBPR()<10) 
     {  
      hNBXT4->Fill(e->get_interPR()-e->get_slopePR()*(25.09012));//Middle of T4 
     }     
    float zz0=0;
    float a=e->get_aPR();
    float b=e->get_bPR();
    float c=e->get_cPR();
         //TF1* fB=new TF1("fB","pol2",-20,20);
    TF1* fB= new TF1("fB","[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,40);
    fB->FixParameter(3,zz0);
    fB->FixParameter(0,c);
    fB->FixParameter(1,b);
    fB->FixParameter(2,a);

    float limo=zL[5];//z position of 6th layer
    TF1*outcomingB=new TF1("outcomingB","pol1",-25,25);
    float aaout=fB->Eval(limo);
    float diffout=2*a*limo+2*a*zz0+b;
    outcomingB->FixParameter(0,-(aaout-diffout*limo)/diffout);
    outcomingB->FixParameter(1,1./diffout);

    if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&e->get_slopePR()!=0&e->get_chi2NBPR()<10) 
     {  
      hNBYT4->Fill(outcomingB->GetX(-25.09012));//Middle of T4 
     }    
    if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&e->get_PHA6().at(0)>0&&e->get_slopePR()!=0&e->get_chi2NBPR()<10) 
     {  
      hNBYLEET4->Fill(outcomingB->GetX(-25.59012-2.54*1.5));//Middle of old LEE T4 
     }    
        
    ////////////////////////////////////
    //OCCUPANCY
    //////////////////////////////////// 
    //if(e->get_Ti()!=127)      continue;

    //if(e->get_EneT1().at(0)<0 ||e->get_EneT4().at(0)<0)
    // {   
    //  continue;
    // }
    
    int* Ltmp=new int[7];
    for(int i=0;i<7;i++)Ltmp[i]=0;
    for(int i=0;i<nnhits;i++)
      { 
       if(e->get_hits().at(i)->get_flagPR()!=1)continue;
       int L=e->get_hits().at(i)->get_L();
       int fStripID=e->get_hits().at(i)->get_fstripID();
       fStripOcc[L]->Fill(fStripID);
       int chip=e->get_hits().at(i)->get_chip();
       ChipOcc[L]->Fill(chip);
       float x=e->get_hits().at(i)->get_x();
       float y=e->get_hits().at(i)->get_y();
       float coord=y;//B plan per default
       if(L==0||L==4||L==6)  coord=x;//NB plane
       CoordOcc[L]->Fill(10*coord);
       //Cluster occupancy
       ClusOcc->Fill(L);
       //Layer occupancy
       Ltmp[L]++;
       if(Ltmp[L]==1) LOcc->Fill(L);
       
      }    
     
     
    ////////////////////////////////////
    //TRIGGER CUT
    ////////////////////////////////////  
      
    //cout << "Internal trigger: " << de->get_Ti() <<endl;
    if(e->get_Ti()!=127)      continue;
    
    //T1&T3&T4
    if(e->get_EneT1().at(0)<0 ||e->get_EneT4().at(0)<0||e->get_EneT3().at(0)<0)
     {   
      continue;
     }
    
    ////////////////////////////////////
    //PULL 
    ////////////////////////////////////      

    for(int i=0;i<nnhits;i++)
      { 
       //Check that the hit is used    
       if(e->get_hits().at(i)->get_flagPR()!=1)continue;
       int L=e->get_hits().at(i)->get_L();//0 to 6
       float x=e->get_hits().at(i)->get_x();
       float y=e->get_hits().at(i)->get_y();
       float z=e->get_hits().at(i)->get_z();
       float xPR=e->get_hits().at(i)->get_xPR();
       float yPR=e->get_hits().at(i)->get_yPR();
       float zPR=e->get_hits().at(i)->get_zPR();  
       XReso[L]->Fill(10*(xPR-x));
       YReso[L]->Fill(10*(yPR-y));
       ZReso[L]->Fill(10*(zPR-z));   
       XResoX[L]->Fill(10*x,10*(xPR-x));
       YResoY[L]->Fill(10*y,10*(yPR-y));
       if(L==0||L==4||L==6)ZResoXY[L]->Fill(10*x,10*(zPR-z)); 
       else ZResoXY[L]->Fill(10*y,10*(zPR-z)); 
      } 
   }//j
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 string PDFFile=Form("Data_%s_%s.pdf",file[0].c_str(),coinc[0].c_str());

 //string PDFFile="Data_NL0089.T1T3T4_FIRSTGEO.pdf";
// string PDFFile="Data_NL0086.T1T4.pdf";
// string PDFFile="Data_17_09_18to20.T1T4.pdf";
 //string PDFFile="Data_17_09_20.T1T4.pdf";
 
 TCanvas*CanPH=new TCanvas("Pulse heights","Pulse heights",200,10,1200,800);  
 CanPH->Divide(3,2);
 
 for(int j=0;j<5;j++)
   {
    CanPH->cd(j+1);  
    gPad->SetLogy();
    AllPH[j]->Draw();
       
   }
  CanPH->Print(Form("%s(",PDFFile.c_str()),Form("Title:Pulse heights"));    


//////////////////////////////    
 TCanvas*CanPHT3G=new TCanvas("Pulse heights T3 vs Guard","Pulse heights T3 vs Guard",200,10,1200,800);  
 CanPHT3G->Divide(2,2);
 CanPHT3G->cd(1);     
 gPad->SetLogy();
 T3woG->Draw();
 CanPHT3G->cd(2);     
 gPad->SetLogy();
 T3wG->Draw();
 CanPHT3G->cd(3);     
 gPad->SetLogy();
 GwoT3->Draw();
 CanPHT3G->cd(4);     
 gPad->SetLogy();
 GwT3->Draw();
 CanPHT3G->Print(Form("%s",PDFFile.c_str()),Form("Title:T3 and Guard"));    


//////////////////////////////    
 TCanvas*CanLOcc=new TCanvas("Layer occupancy","Layer occupancy",200,10,1400,600);  
 CanLOcc->Divide(2,1);
 CanLOcc->cd(1);
 LOcc->GetYaxis()->SetTitleOffset(1.1);
 LOcc->Draw();
 CanLOcc->cd(2);
 ClusOcc->GetYaxis()->SetTitleOffset(1.1);
 ClusOcc->Draw();
 CanLOcc->Print(Form("%s",PDFFile.c_str()),Form("Title:Occupancy"));    
 
 TCanvas**CanHOcc=new TCanvas*[7];
 
 for(int i=0;i<7;i++)
   {
    CanHOcc[i]=new TCanvas(Form("Hit occupancy, L%d",i),Form("Hit occupancy, L%d",i),200,10,1400,1200);  
    CanHOcc[i]->Divide(2,2);
    CanHOcc[i]->cd(1);
    ChipOcc[i]->GetYaxis()->SetTitleOffset(1.1);
    ChipOcc[i]->Draw();
    CanHOcc[i]->cd(2);
    fStripOcc[i]->GetYaxis()->SetTitleOffset(1.1);
    fStripOcc[i]->Draw();
    CanHOcc[i]->cd(3);
    CoordOcc[i]->GetYaxis()->SetTitleOffset(1.1);
    CoordOcc[i]->Draw();
    CanHOcc[i]->Print(Form("%s",PDFFile.c_str()),Form("Title:Occupancy L%d",i));    
   }
 

//////////////////////////////    
 TCanvas*CanDef=new TCanvas("Deflections","Deflections",200,10,1400,600);  
 CanDef->Divide(2,1);
 CanDef->cd(1);
 hDef->GetYaxis()->SetTitleOffset(1.1);
 hDef->Draw();
 CanDef->cd(2);
 DefvsT1T3T4->GetYaxis()->SetTitleOffset(1.1);
 DefvsT1T3T4->Draw("colz");
 gPad->Update();
 TPaveStats* ps = (TPaveStats *)DefvsT1T3T4->GetListOfFunctions()->FindObject("stats");
 ps->SetX1NDC(0.68);
 ps->SetY1NDC(0.60);
 ps->SetX2NDC(0.88);
 ps->SetY2NDC(0.88);
 gPad->Update();
 gPad->Modified();
 CanDef->Print(Form("%s",PDFFile.c_str()),Form("Title:Deflection"));    

//////////////////////////////    

 TCanvas*CanNB=new TCanvas("Non Bending plane","Non Bending plane",200,10,1800,1200);  
 CanNB->Divide(3,2);
 CanNB->cd(1);
 hNBslope->GetYaxis()->SetTitleOffset(1.2);
 hNBslope->Draw();
 CanNB->cd(2);
 hNBinter->GetYaxis()->SetTitleOffset(1.2);
 hNBinter->Draw();
 CanNB->cd(3);
 hNBXT1->GetYaxis()->SetTitleOffset(1.2);
 hNBXT1->Draw();
 CanNB->cd(4);
 hNBXT4->GetYaxis()->SetTitleOffset(1.2);
 hNBXT4->Draw();
 hNBXLEET4->SetLineColor(kBlack);
 hNBXLEET4->Draw("same");
 CanNB->cd(5);
 hNBYT4->GetYaxis()->SetTitleOffset(1.2);
 hNBYT4->Draw();
 hNBYLEET4->SetLineColor(kBlack);
 hNBYLEET4->Draw("same");

 CanNB->Print(Form("%s",PDFFile.c_str()),Form("Title:Non Bending plane"));    

 
//////////////////////////////    
 TCanvas**CanXYZReso=new TCanvas*[7];  
 int Nsig=5;

 for(int j=0;j<7;j++)
   {
    CanXYZReso[j]=new TCanvas(Form("CanXYZReso L%d",j),Form("CanXYZReso L%d",j),200,10,1200,1200);     
    CanXYZReso[j]->Divide(2,2);
    CanXYZReso[j]->cd(1);
    if(j==1||j==2||j==3||j==5)
     {
      YReso[j]->SetLineColor(kRed);
      YReso[j]->SetTitle(Form("Tracker B layer %d",j));
      YReso[j]->GetXaxis()->SetTitle("#Delta Y (mm)"); 
      YReso[j]->GetXaxis()->SetRangeUser(YReso[j]->GetMean()-Nsig*YReso[j]->GetRMS(),YReso[j]->GetMean()+Nsig*YReso[j]->GetRMS());
      YReso[j]->Draw();
     }  
    else
     {
      XReso[j]->SetLineColor(kRed);
      XReso[j]->SetTitle(Form("Tracker NB layer %d",j));
      XReso[j]->GetXaxis()->SetTitle("#Delta X (mm)");             
      XReso[j]->GetXaxis()->SetRangeUser(XReso[j]->GetMean()-Nsig*XReso[j]->GetRMS(),XReso[j]->GetMean()+Nsig*XReso[j]->GetRMS());
      XReso[j]->Draw();
     }
    CanXYZReso[j]->cd(2);
    ZReso[j]->SetLineColor(kRed);
    ZReso[j]->SetTitle(Form("Tracker layer %d",j));
    ZReso[j]->GetXaxis()->SetTitle("#Delta Z (mm)");             
    ZReso[j]->GetXaxis()->SetRangeUser(ZReso[j]->GetMean()-Nsig*ZReso[j]->GetRMS(),ZReso[j]->GetMean()+Nsig*ZReso[j]->GetRMS());
    ZReso[j]->Draw();
    CanXYZReso[j]->cd(3);
    if(j==1||j==2||j==3||j==5)
     {
      YResoY[j]->SetTitle(Form("Tracker B layer %d",j));
      YResoY[j]->GetYaxis()->SetTitle("#Delta Y (mm)");      
      YResoY[j]->GetXaxis()->SetTitle("Y (mm)");      
      YResoY[j]->Draw("colz");
     }  
    else
     {
      XResoX[j]->SetTitle(Form("Tracker NB layer %d",j));
      XResoX[j]->GetYaxis()->SetTitle("#Delta X (mm)");             
      XResoX[j]->GetXaxis()->SetTitle("X (mm)");             
      XResoX[j]->Draw("colz");
     }
    CanXYZReso[j]->cd(4);
    if(j==1||j==2||j==3||j==5)
     {
      ZResoXY[j]->SetTitle(Form("Tracker B layer %d",j));
      ZResoXY[j]->GetYaxis()->SetTitle("#Delta Z (mm)");             
      ZResoXY[j]->GetXaxis()->SetTitle("Y (mm)");             
      ZResoXY[j]->Draw("colz");
     }
    else
     {
      ZResoXY[j]->SetTitle(Form("Tracker NB layer %d",j));
      ZResoXY[j]->GetYaxis()->SetTitle("#Delta Z (mm)");             
      ZResoXY[j]->GetXaxis()->SetTitle("X (mm)");             
      ZResoXY[j]->Draw("colz");
         
     } 
    if(j==6)CanXYZReso[j]->Print(Form("%s)",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
    else CanXYZReso[j]->Print(Form("%s",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
     
   }//j   
    
  
    
}

void PullDistribution(string s,string c,int geoconf)
{
  //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

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

 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;     

//////////////////////////// 
//Resolution histograms
//////////////////////////// 

 TH1F**XReso=new TH1F*[7];
 TH1F**YReso=new TH1F*[7];
 TH1F**ZReso=new TH1F*[7];
 TH2F**XResoX=new TH2F*[7];
 TH2F**YResoY=new TH2F*[7];
 TH2F**ZResoXY=new TH2F*[7];
 
 for(int j=0;j<7;j++)
   {
    XReso[j]=new TH1F(Form("Xreso L%d",j),Form("Xreso L%d",j),200,-1,1); 
    YReso[j]=new TH1F(Form("Yreso L%d",j),Form("Yreso L%d",j),200,-1,1);
    ZReso[j]=new TH1F(Form("Zreso L%d",j),Form("Zreso L%d",j),200,-3,3);
    XResoX[j]=new TH2F(Form("Xreso vs X L%d",j),Form("Xreso vs X L%d",j),1000,-100,100,400,-1,1); 
    YResoY[j]=new TH2F(Form("Yreso vs Y L%d",j),Form("Yreso vs Y L%d",j),1000,-100,100,400,-1,1);
    if(j==0||j==4||j==6)ZResoXY[j]=new TH2F(Form("Zreso vs X L%d",j),Form("Zreso vs X L%d",j),1000,-100,100,400,-2,2);
    else ZResoXY[j]=new TH2F(Form("Zreso vs Y L%d",j),Form("Zreso vs Y L%d",j),1000,-100,100,400,-2,2);
   }
    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");
 chain->Add(Form("../Data/NL0248.BPD.EVENT.root"));
 chain->Add(Form("../Data/%s.EVENT.root",file[0].c_str()));
 chain->Add(Form("../Data/NL0253.BPD.EVENT.root"));
 chain->Add(Form("../Data/NL0254.BPD.EVENT.root"));
    
 //Define variables to read event
 ALEvent *e = new ALEvent();      
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 int nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;  
      
 for (int j=0;j<nentries;j++)
   {
    chain->GetEntry(j); //Load the entry i in the variable e 

    if(j%10000==0) cout << "Event: " << j <<endl;
    
    int nnhits = (int) e->get_Nhits();
    // if(nnhits<7) continue;
    uint8_t Ti=(uint8_t)e->get_Ti();
    //Number of layers wih hit(s)
    int NL=0;
    for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
      
    if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
    int* Lay=new int[7];
    for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);
      
    
    if(NL!=7)continue;
    
    ////////////////////////////////////
    //TRIGGER CUT
    ////////////////////////////////////  
      
    //cout << "Internal trigger: " << de->get_Ti() <<endl;
   // if(e->get_Ti()!=127)      continue;
    
    //T1&T3&T4
    if(e->get_EneT1().at(0)<0 ||e->get_EneT4().at(0)<0)
     {   
      continue;
     }
    
    ////////////////////////////////////
    //PULL 
    ////////////////////////////////////      

    for(int i=0;i<nnhits;i++)
      { 
       //Check that the hit is used    
       if(e->get_hits().at(i)->get_flagPR()!=1)continue;
       int L=e->get_hits().at(i)->get_L();//0 to 6
       float x=e->get_hits().at(i)->get_x();
       float y=e->get_hits().at(i)->get_y();
       float z=e->get_hits().at(i)->get_z();
       float xPR=e->get_hits().at(i)->get_xPR();
       float yPR=e->get_hits().at(i)->get_yPR();
       float zPR=e->get_hits().at(i)->get_zPR();  
       XReso[L]->Fill(10*(xPR-x));
       YReso[L]->Fill(10*(yPR-y));
       ZReso[L]->Fill(10*(zPR-z));   
       XResoX[L]->Fill(10*x,10*(xPR-x));
       YResoY[L]->Fill(10*y,10*(yPR-y));
       if(L==0||L==4||L==6)ZResoXY[L]->Fill(10*x,10*(zPR-z)); 
       else ZResoXY[L]->Fill(10*y,10*(zPR-z)); 
      } 
   }//j
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 string PDFFile=Form("PullDistribution%s_%s.pdf",file[0].c_str(),coinc[0].c_str());

//////////////////////////////    
 TCanvas**CanXYZReso=new TCanvas*[7];  
 int Nsig=5;

 for(int j=0;j<7;j++)
   {
    CanXYZReso[j]=new TCanvas(Form("CanXYZReso L%d",j),Form("CanXYZReso L%d",j),200,10,1200,1200);     
    CanXYZReso[j]->Divide(2,2);
    CanXYZReso[j]->cd(1);
    if(j==1||j==2||j==3||j==5)
     {
      cout << "L="<<j << "PullY mean= " << YReso[j]->GetMean() << " mm"<<endl;
      YReso[j]->SetLineColor(kRed);
      YReso[j]->SetTitle(Form("Tracker B layer %d",j));
      YReso[j]->GetXaxis()->SetTitle("#Delta Y (mm)"); 
      //YReso[j]->GetXaxis()->SetRangeUser(YReso[j]->GetMean()-Nsig*YReso[j]->GetRMS(),YReso[j]->GetMean()+Nsig*YReso[j]->GetRMS());
      YReso[j]->Draw();
     }  
    else
     {
      cout << "L="<<j << "PullX mean= " << XReso[j]->GetMean() << " mm"<<endl;
      XReso[j]->SetLineColor(kRed);
      XReso[j]->SetTitle(Form("Tracker NB layer %d",j));
      XReso[j]->GetXaxis()->SetTitle("#Delta X (mm)");             
      //XReso[j]->GetXaxis()->SetRangeUser(XReso[j]->GetMean()-Nsig*XReso[j]->GetRMS(),XReso[j]->GetMean()+Nsig*XReso[j]->GetRMS());
      XReso[j]->Draw();
     }
    CanXYZReso[j]->cd(2);
    cout << "L="<<j << "PullZ mean= " << ZReso[j]->GetMean() << " mm"<<endl;
    ZReso[j]->SetLineColor(kRed);
    ZReso[j]->SetTitle(Form("Tracker layer %d",j));
    ZReso[j]->GetXaxis()->SetTitle("#Delta Z (mm)");             
    //ZReso[j]->GetXaxis()->SetRangeUser(ZReso[j]->GetMean()-Nsig*ZReso[j]->GetRMS(),ZReso[j]->GetMean()+Nsig*ZReso[j]->GetRMS());
    ZReso[j]->Draw();
    CanXYZReso[j]->cd(3);
    if(j==1||j==2||j==3||j==5)
     {
      YResoY[j]->SetTitle(Form("Tracker B layer %d",j));
      YResoY[j]->GetYaxis()->SetTitle("#Delta Y (mm)");      
      YResoY[j]->GetXaxis()->SetTitle("Y (mm)");      
      YResoY[j]->Draw("colz");
     }  
    else
     {
      XResoX[j]->SetTitle(Form("Tracker NB layer %d",j));
      XResoX[j]->GetYaxis()->SetTitle("#Delta X (mm)");             
      XResoX[j]->GetXaxis()->SetTitle("X (mm)");             
      XResoX[j]->Draw("colz");
     }
    CanXYZReso[j]->cd(4);
    if(j==1||j==2||j==3||j==5)
     {
      ZResoXY[j]->SetTitle(Form("Tracker B layer %d",j));
      ZResoXY[j]->GetYaxis()->SetTitle("#Delta Z (mm)");             
      ZResoXY[j]->GetXaxis()->SetTitle("Y (mm)");             
      ZResoXY[j]->Draw("colz");
     }
    else
     {
      ZResoXY[j]->SetTitle(Form("Tracker NB layer %d",j));
      ZResoXY[j]->GetYaxis()->SetTitle("#Delta Z (mm)");             
      ZResoXY[j]->GetXaxis()->SetTitle("X (mm)");             
      ZResoXY[j]->Draw("colz");
         
     } 
    if(j==6)CanXYZReso[j]->Print(Form("%s)",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
    else if(j==0)CanXYZReso[j]->Print(Form("%s(",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
    else CanXYZReso[j]->Print(Form("%s",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
     
   }//j   
    
  
    
}



void Occupancy(string s,string c)
{
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;   

 int Ncut=4;
 
 string scut[4]={"All (1+ hit)","T1&T4& 1+ hit","T1&T4& 7 layers","T1&T4& 7 layers & used hits"};
// string scut[4]={"All","T1&T4&T3","T1&T4&T3&7 layers","T1&T4&T3&7 layers & used hits"};
 
 //string sBoard[7]={"A","B","C","D","E","H","G"};	
 //string sBoard[7]={"H","B","C","D","E","I","A"};	
 //string sBoard[7]={"G","B","H","D","E","I","A"};	
 string sBoard[7]={"C","B","A","D","E","I","G"};	
	

//////////////////////////// 
//Occupancy histograms
//////////////////////////// 
    
 TH1F** LOcc=new TH1F*[Ncut];
 TH1F** ClusOcc=new TH1F*[Ncut];
 TH1F*** ChipOcc=new TH1F**[Ncut];
 TH1F*** fStripOcc=new TH1F**[Ncut];
 TH1F*** StripOcc=new TH1F**[Ncut];
 TH1F*** CoordOcc=new TH1F**[Ncut];
 

 double*** LvsChipperL=new double**[7];
 for(int i=0;i<7;i++) 
   {
    LvsChipperL[i]=new double*[7];
    for(int j=0;j<7;j++) 
      { 
       LvsChipperL[i][j]=new double[12]; 
       for(int k=0;k<12;k++) 
         { 
          LvsChipperL[i][j][k]=0;   
         }       
      }
   } 
   
 TH2F**ChipperLEff=new TH2F*[7];
 for(int i=0;i<7;i++) 
   {
    ChipperLEff[i]=new TH2F(Form("Efficiency for events with layer %i hitted",i),Form("Efficiency for events with layer %i hitted",i),12,0,12,7,0,7);
    ChipperLEff[i]->GetXaxis()->SetTitle("Chip (0 to 11)");
    ChipperLEff[i]->GetYaxis()->SetTitle("Layer (0 to 6)");
    ChipperLEff[i]->SetTitle(Form("Efficiency for events with hit(s) on L%i (Board %s)",i,sBoard[i].c_str()));   
   } 
    
    
 TH2F*Chip2Doccupancy=new TH2F(Form("%% of Occupancy"),Form("%% of Occupancy"),12,0,12,7,0,7);
 Chip2Doccupancy->GetXaxis()->SetTitle("Chip (0 to 11)");
 Chip2Doccupancy->GetYaxis()->SetTitle("Layer (0 to 6)");
 Chip2Doccupancy->SetTitle(Form("%% of Occupancy per chip"));   

 
 for(int i=0;i<Ncut;i++)
   {  
    LOcc[i]=new TH1F(Form("Layer occupancy, %s",scut[i].c_str()),Form("Layer occupancy, %s",scut[i].c_str()),10,0,10);
    LOcc[i]->GetXaxis()->SetTitle("Layer (0 to 6)");
    LOcc[i]->GetYaxis()->SetTitle("% of total number of events");
    LOcc[i]->SetLineColor(1+i);
    LOcc[i]->SetLineWidth(1);  
    LOcc[i]->SetMinimum(0);  

    ClusOcc[i]=new TH1F(Form("Cluster occupancy, %s",scut[i].c_str()),Form("Cluster occupancy, %s",scut[i].c_str()),10,0,10);
    ClusOcc[i]->GetXaxis()->SetTitle("Layer (0 to 6)");
    ClusOcc[i]->GetYaxis()->SetTitle("Cluster number");
    ClusOcc[i]->SetLineColor(1+i);
    ClusOcc[i]->SetLineWidth(1);  
    ClusOcc[i]->SetMinimum(0);  
    
    ChipOcc[i]=new TH1F*[7];
    fStripOcc[i]=new TH1F*[7];
    StripOcc[i]=new TH1F*[7];
    CoordOcc[i]=new TH1F*[7];
    for(int j=0;j<7;j++)
      {    
       ChipOcc[i][j]=new TH1F(Form("Chip ID L%d, %s",j,scut[i].c_str()),Form("Chip ID L%d, %s",j,scut[i].c_str()),15,0,15);
       ChipOcc[i][j]->GetXaxis()->SetTitle("Chip ID");
       ChipOcc[i][j]->GetYaxis()->SetTitle("Entries");
       ChipOcc[i][j]->SetLineColor(1+i);
       ChipOcc[i][j]->SetLineWidth(1);
       ChipOcc[i][j]->GetXaxis()->SetRangeUser(0,15);  
       ChipOcc[i][j]->SetMinimum(0);  
       fStripOcc[i][j]=new TH1F(Form("First strip ID L%d, %s",j,scut[i].c_str()),Form("First strip ID L%d, %s",j,scut[i].c_str()),769,0,769);
       fStripOcc[i][j]->GetXaxis()->SetTitle("First strip ID");
       fStripOcc[i][j]->GetYaxis()->SetTitle("Entries");
       fStripOcc[i][j]->SetLineColor(1+i);
       fStripOcc[i][j]->SetLineWidth(1);
       fStripOcc[i][j]->GetXaxis()->SetRangeUser(0,769);  
       fStripOcc[i][j]->SetMinimum(0);  
       StripOcc[i][j]=new TH1F(Form("Strip L%d, %s",j,scut[i].c_str()),Form("Strip L%d, %s",j,scut[i].c_str()),769,0,769);
       StripOcc[i][j]->GetXaxis()->SetTitle("Strip");
       StripOcc[i][j]->GetYaxis()->SetTitle("Entries");
       StripOcc[i][j]->SetLineColor(1+i);
       StripOcc[i][j]->SetLineWidth(1);
       StripOcc[i][j]->GetXaxis()->SetRangeUser(0,769);  
       StripOcc[i][j]->SetMinimum(0);  
       CoordOcc[i][j]=new TH1F(Form("Cluster coord. L%d, %s",j,scut[i].c_str()),Form("Cluster coord. L%d, %s",j,scut[i].c_str()),200,-100,100);
       if(j==0||j==4||j==6)CoordOcc[i][j]->GetXaxis()->SetTitle("X (mm)");
       else CoordOcc[i][j]->GetXaxis()->SetTitle("Y (mm)");
       CoordOcc[i][j]->GetYaxis()->SetTitle("Entries");
       CoordOcc[i][j]->SetLineColor(1+i);
       CoordOcc[i][j]->SetLineWidth(1);
       CoordOcc[i][j]->GetXaxis()->SetRangeUser(-100,100);    
       CoordOcc[i][j]->SetMinimum(0);
      }//j
   }//i

    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");

 // chain->Add(Form("../Data/%s.EVENT.root",file[0].c_str()));
 //Esrange
 chain->Add(Form("../Data/NL2084.BPD.EVENT.root"));
 chain->Add(Form("../Data/NL2085.BPD.EVENT.root"));
 chain->Add(Form("../Data/NL2088.BPD.EVENT.root"));
    
 //Define variables to read event
 ALEvent *e = new ALEvent();      
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 int nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;  
      
 for (int j=0;j<nentries;j++)
   {
    chain->GetEntry(j); //Load the entry i in the variable e 

    if(j%50000==0) cout << "Event: " << j <<endl;
    
    int nnhits = (int) e->get_Nhits();
     
    ////////////////////////////////////
    //OCCUPANCY
    //////////////////////////////////// 
    //if(e->get_Ti()!=127)      continue;

    //if(e->get_EneT1().at(0)<0 ||e->get_EneT4().at(0)<0)
    // {   
    //  continue;
    // }
    
    int* Ltmp=new int[7];
    for(int i=0;i<7;i++)Ltmp[i]=0;
    int** LCtmp=new int*[7];
    for(int i=0;i<7;i++)
      {
       LCtmp[i]=new int[12]; 
       for(int jk=0;jk<12;jk++)
         {
          LCtmp[i][jk]=0;
         }
      }   
    unsigned int  Ti=(unsigned int) e->get_Ti(); 
    //internal trigger
    int* Lay=new int[7];
    for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);   
    
    for(int i=0;i<nnhits;i++)
      { 
       //if(e->get_hits().at(i)->get_flagPR()!=1)continue;
       int L=e->get_hits().at(i)->get_L();
       int fStripID=e->get_hits().at(i)->get_fstripID();
       fStripOcc[0][L]->Fill(fStripID);
       StripOcc[0][L]->Fill(fStripID);
       int nstrips=e->get_hits().at(i)->get_nstrips();
       for(int ij=1;ij<=nstrips;ij++) StripOcc[0][L]->Fill(fStripID+ij);
      
       int chip=e->get_hits().at(i)->get_chip();
       ChipOcc[0][L]->Fill(chip);
       float x=e->get_hits().at(i)->get_x();
       float y=e->get_hits().at(i)->get_y();
       float coord=y;//B plan per default
       if(L==0||L==4||L==6)  coord=x;//NB plane
       CoordOcc[0][L]->Fill(10*coord);
       //Cluster occupancy
       ClusOcc[0]->Fill(L);
       //Layer occupancy
       Ltmp[L]++;
       LCtmp[L][chip]++;
       if(Ltmp[L]==1)
        {
         LOcc[0]->Fill(L);
        } 
       if(LCtmp[L][chip]==1)
        {
         Chip2Doccupancy->Fill(chip,L,1);
         for(int ijk=0;ijk<7;ijk++)
           {
            if(Lay[ijk]==1)   LvsChipperL[ijk][L][chip]++;
           } 
        }
       
       //Second cut
       //T1 and T4 and 1+ hit
       if(e->get_EneT1().at(0)>0 && e->get_EneT4().at(0)>0)
        {   
         fStripOcc[1][L]->Fill(fStripID);
         StripOcc[1][L]->Fill(fStripID);
         for(int ij=1;ij<=nstrips;ij++) StripOcc[1][L]->Fill(fStripID+ij);
         ChipOcc[1][L]->Fill(chip);
         CoordOcc[1][L]->Fill(10*coord);
         ClusOcc[1]->Fill(L);
         if(Ltmp[L]==1) LOcc[1]->Fill(L);		
		}	
	      		
        //Third cut
        //7 layers with hits
        if(e->get_EneT1().at(0)>0 && e->get_EneT4().at(0)>0 && e->get_Ti()==127)
         {  
          fStripOcc[2][L]->Fill(fStripID);
          StripOcc[2][L]->Fill(fStripID);
          for(int ij=1;ij<=nstrips;ij++) StripOcc[2][L]->Fill(fStripID+ij);
          ChipOcc[2][L]->Fill(chip);
          CoordOcc[2][L]->Fill(10*coord);
          ClusOcc[2]->Fill(L);
          if(Ltmp[L]==1) LOcc[2]->Fill(L);
		 }	  
        //Fourth cut
        //only hits used in reconstructionBoard[i].c_str()
        if(e->get_EneT1().at(0)>0 && e->get_EneT4().at(0)>0 && e->get_Ti()==127&&e->get_hits().at(i)->get_flagPR()==1)
         {   
          fStripOcc[3][L]->Fill(fStripID);
          StripOcc[3][L]->Fill(fStripID);
          for(int ij=1;ij<=nstrips;ij++) StripOcc[3][L]->Fill(fStripID+ij);
          ChipOcc[3][L]->Fill(chip);
          CoordOcc[3][L]->Fill(10*coord);
          ClusOcc[3]->Fill(L);
          if(Ltmp[L]==1) LOcc[3]->Fill(L);
         }//fourth cut 
      }//i 
   }//j
   

//Normalize to get hardware trigger efficiency   
 //float scaler=100./nentries;
 float scaler=1;
 for(int i=0;i<Ncut;i++)
   {    
    LOcc[i]->Scale(scaler);
  //  LOcc[i]->SetMaximum(scaler);
    //ClusOcc[i]->Scale(scaler);
    // ClusOcc[i]->SetMaximum(scaler);
   
    for(int j=0;j<7;j++)
       {
        //ChipOcc[i][j]->Scale(scaler);
       // ChipOcc[i][j]->SetMaximum(scaler);
        //fStripOcc[i][j]->Scale(scaler);
       // fStripOcc[i][j]->SetMaximum(scaler);
       // CoordOcc[i][j]->Scale(scaler);
       // CoordOcc[i][j]->SetMaximum(scaler);
        //StripOcc[i][j]->Scale(scaler);
        //StripOcc[i][j]->SetMaximum(scaler);
      }    
   }
  
  
  for(int i=0;i<7;i++) 
   {
    for(int j=0;j<7;j++) 
      { 
       double tmpChip=0;
       for(int k=0;k<12;k++) 
         { 
          ChipperLEff[i]->Fill(k,j,100.*LvsChipperL[i][j][k]/LOcc[0]->GetBinContent(1));
          tmpChip+=LvsChipperL[i][j][k];
         }     
       //ChipperLEff[i]->Fill(12,j,100.*tmpChip/LOcc[0]->GetBinContent(1));
         
      }
   } 
 
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
 gStyle->SetPalette( kDarkBodyRadiator);
 gStyle->SetPalette(kTemperatureMap);
 TGaxis::SetMaxDigits(4);

//////////////////////////////    

 string PDFFile=Form("Occupancies_%s_%s.pdf",file[0].c_str(),coinc[0].c_str());
 //string PDFFile="Occupancies_17_09_18to20_T1T4.pdf";
 //string PDFFile="Occupancies_17_09_20_T1T4.pdf";

//////////////////////////////    
 TCanvas*CanLOcc=new TCanvas("Layer occupancy","Layer occupancy",200,10,1400,600);  
 TLegend* leg=new TLegend(0.12,0.12,0.6,0.4);
 leg->SetBorderSize(0);leg->SetFillStyle(0);
 for(int i=0;i<Ncut;i++) leg->AddEntry(LOcc[i],Form("%s",scut[i].c_str()),"l");
 CanLOcc->Divide(2,1);
 CanLOcc->cd(1);
 LOcc[0]->GetYaxis()->SetTitleOffset(1.5);
 LOcc[0]->Draw("hist");
 for(int i=1;i<Ncut;i++)LOcc[i]->Draw("hist same");
 leg->Draw("same");
 CanLOcc->cd(2);
 ClusOcc[0]->GetYaxis()->SetTitleOffset(1.5);
 ClusOcc[0]->Draw("hist");
 for(int i=1;i<Ncut;i++)ClusOcc[i]->Draw("hist same");
 CanLOcc->Print(Form("%s(",PDFFile.c_str()),Form("Title:Occupancy"));    
 
 TCanvas**CanHOcc=new TCanvas*[7];
 
 for(int i=0;i<7;i++)
   {
    CanHOcc[i]=new TCanvas(Form("Hit occupancy, L%d",i),Form("Hit occupancy, L%d",i),200,10,1400,1200);  
    CanHOcc[i]->Divide(2,2);
    CanHOcc[i]->cd(1);
    ChipOcc[0][i]->GetYaxis()->SetTitleOffset(1.5);
    ChipOcc[0][i]->Draw("hist");
    for(int j=1;j<Ncut;j++)ChipOcc[j][i]->Draw("hist same");
    leg->Draw("same");
    CanHOcc[i]->cd(2);
    fStripOcc[0][i]->GetYaxis()->SetTitleOffset(1.5);
    fStripOcc[0][i]->SetMaximum(3*fStripOcc[0][i]->GetEntries()/fStripOcc[0][i]->GetNbinsX());
    fStripOcc[0][i]->Draw("hist");
    for(int j=1;j<Ncut;j++)fStripOcc[j][i]->Draw("hist same");
    CanHOcc[i]->cd(3);
    CoordOcc[0][i]->GetYaxis()->SetTitleOffset(1.5);
    CoordOcc[0][i]->SetMaximum(3*CoordOcc[0][i]->GetEntries()/CoordOcc[0][i]->GetNbinsX());
    CoordOcc[0][i]->Draw("hist");
    for(int j=1;j<Ncut;j++)CoordOcc[j][i]->Draw("hist same");
    CanHOcc[i]->cd(4);
    StripOcc[0][i]->GetYaxis()->SetTitleOffset(1.5);
    StripOcc[0][i]->SetMaximum(3*StripOcc[0][i]->GetEntries()/StripOcc[0][i]->GetNbinsX());
    StripOcc[0][i]->Draw("hist");
    for(int j=1;j<Ncut;j++)StripOcc[j][i]->Draw("hist same");
    //if(i==6)CanHOcc[i]->Print(Form("%s)",PDFFile.c_str()),Form("Title:Occupancy L%d (B%s)",i,sBoard[i].c_str()));    
     CanHOcc[i]->Print(Form("%s",PDFFile.c_str()),Form("Title:Occupancy L%d, (B%s)",i,sBoard[i].c_str()));    
   }
 
  gStyle->SetPaintTextFormat("4.1f");
  TCanvas**CanLvsCperL=new TCanvas*[7];
  
  for(int i=0;i<7;i++)
   {
    CanLvsCperL[i]=new TCanvas(Form("CanLvsCperL %d",i),Form("CanLvsCperL %d",i),200,10,1500,900);
    CanLvsCperL[i]->cd();
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    ChipperLEff[i]->GetXaxis()->SetNdivisions(12,0,0);
    ChipperLEff[i]->GetYaxis()->SetNdivisions(7,0,0);
    ChipperLEff[i]->Draw("colz text");
    ChipperLEff[i]->SetContour(60);
    ChipperLEff[i]->SetMarkerSize(2);
    ChipperLEff[i]->SetMaximum(15);
    ChipperLEff[i]->SetMinimum(0);
    CanLvsCperL[i]->Print(Form("%s",PDFFile.c_str()),Form("Title:Chip occupancy, fired L%d, (B%s)",i,sBoard[i].c_str()));    
   }    
 
 
 Chip2Doccupancy->Scale(100./nentries);
 TCanvas*CanChip2Doccupancy=new TCanvas("2D Chip occupancy","2D Chip occupancy",200,10,1500,900);
 CanChip2Doccupancy->cd();
 gPad->SetGridx(1);
 gPad->SetGridy(1);
 Chip2Doccupancy->GetXaxis()->SetNdivisions(12,0,0);
 Chip2Doccupancy->GetYaxis()->SetNdivisions(7,0,0);
 Chip2Doccupancy->Draw("colz text");
 Chip2Doccupancy->SetContour(60);
 Chip2Doccupancy->SetMarkerSize(2);
 Chip2Doccupancy->SetMaximum(10);
 Chip2Doccupancy->SetMinimum(0);
 CanChip2Doccupancy->Print(Form("%s)",PDFFile.c_str()),Form("Title:2D Chip occupancy"));    
    

 
 
}


void Pulses(string s,string c)
{

 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;   
    
    
//////////////////////////// 
//Trigger related histograms
//////////////////////////// 
 
 TH1F**AllPH=new TH1F*[5];
 string sPH[5]={"T1","T2","T3","T4","Guard"};
 string XtitlePH="Pulse height";
 string YtitlePH="Entries";
	
 const Int_t nx=20;
 string T[nx]={" All","T1","T2","T3","T4","G","T1T4","T1T4 & hits","T1T2","T1T2 & hits","T1T3","T1T3 & hits","T1T2T3","T1T2T3 & hits","T1T3T4","T1T3T4 & hits","T1T2T3T4","T1T2T3T4 & hits","T1T2T3T4 & NoG","T1T2T3T4 & NoG & hits"};
 TH1F*TEff=new TH1F("Rates","Rates",20,0,20);
 TEff->SetStats(0);
 TEff->SetFillColor(38);

 TEff->GetXaxis()->SetAlphanumeric();
 for(int i=1;i<=nx;i++)TEff->GetXaxis()->SetBinLabel(i,T[i-1].c_str());
 
 for(int i=0;i<5;i++)
   {
    AllPH[i]=new TH1F(Form("%s",sPH[i].c_str()),Form("%s",sPH[i].c_str()),1000,-1,4096);
    AllPH[i]->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
    AllPH[i]->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
    AllPH[i]->SetLineColor(kRed);
    AllPH[i]->SetLineWidth(3);
    AllPH[i]->GetXaxis()->SetRangeUser(-1,500);
   } //i  
    
 TH2F*GvsT3=new TH2F("Guard vs T3","Guard vs T3",551,-51,500,551,-51,500);	
	
	
	
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");

 chain->Add(Form("../Data/%s.EVENT.root",file[0].c_str()));

 //Define variables to read event
 ALEvent *e = new ALEvent();      
 //Set address to access event data
 chain->SetBranchAddress("event",&e); 
 
 // Get number of event in Tree
 int nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;  
      
 for (int j=0;j<nentries;j++)
   {
    chain->GetEntry(j); //Load the entry i in the variable e 

    if(j%10000==0) cout << "Event: " << j <<endl;

    ////////////////////////////////////
    //TRIGGER HEIGHTS 
    ////////////////////////////////////  
    
    AllPH[0]->Fill(e->get_EneT1().at(0));
    //if(e->get_EneT3().at(0)>=0)
    AllPH[1]->Fill(e->get_EneT2().at(0));
    AllPH[2]->Fill(e->get_EneT3().at(0));
    AllPH[3]->Fill(e->get_EneT4().at(0));
    AllPH[4]->Fill(e->get_Eneg().at(0));
    if(e->get_EneT3().at(0)>=0||e->get_Eneg().at(0)>=0) GvsT3->Fill(e->get_EneT3().at(0),e->get_Eneg().at(0));

    int nnhits=e->get_Nhits();
     
    ////////////////////////////////////
    //TRIGGER CUT
    ////////////////////////////////////  
     
    //No Cut
    TEff->Fill(T[0].c_str(),1);
    //T1     
    if(e->get_EneT1().at(0)>0)TEff->Fill(T[1].c_str(),1);
    //T2
    if(e->get_EneT2().at(0)>0)TEff->Fill(T[2].c_str(),1);
    //T3
    if(e->get_EneT3().at(0)>0)TEff->Fill(T[3].c_str(),1);
    //T4
    if(e->get_EneT4().at(0)>0)TEff->Fill(T[4].c_str(),1);
    //G
    if(e->get_Eneg().at(0)>0)TEff->Fill(T[5].c_str(),1);
    
    //T1&T4     
    if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0)TEff->Fill(T[6].c_str(),1);
    //T1&T4&hits     
    if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0)TEff->Fill(T[7].c_str(),1);

    //T1&T2     
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0)TEff->Fill(T[8].c_str(),1);
    //T1&T2&hits     
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&nnhits>0)TEff->Fill(T[9].c_str(),1);

    //T1&T3     
    if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0)TEff->Fill(T[10].c_str(),1);
    //T1&T3&hits     
    if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&nnhits>0)TEff->Fill(T[11].c_str(),1);
 
    //"T1&T2&T3"
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0)TEff->Fill(T[12].c_str(),1);
    //"T1&T2&T3&hits
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0)TEff->Fill(T[13].c_str(),1);

    //"T1&T3&T4"
    if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0)TEff->Fill(T[14].c_str(),1);
    //"T1&T3&T4&hits
    if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0)TEff->Fill(T[15].c_str(),1);

    //"T1&T2&T3&T4"
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0)TEff->Fill(T[16].c_str(),1);
    //"T1&T2&T3&T4&hits
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0)TEff->Fill(T[17].c_str(),1);

    //"T1&T2&T3&T4&NoG"
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0)TEff->Fill(T[18].c_str(),1);
    //"T1&T2&T3&T4&NoG&hits
    if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0&&nnhits>0)TEff->Fill(T[19].c_str(),1);
     
   }//j
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 TCanvas*CanGvsT3=new TCanvas("Pulse heights GvsT3","Pulse heights GvsT3",200,10,1000,1000);  
 CanGvsT3->cd();
 gPad->SetLeftMargin(0.15);
 gPad->SetTopMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.15);
 GvsT3->GetXaxis()->SetTitle("Pulse height in T3");
 GvsT3->GetYaxis()->SetTitle("Pulse height in guard");
 GvsT3->GetYaxis()->SetTitleOffset(1.9);
 GvsT3->GetXaxis()->SetRangeUser(-50,500);
 GvsT3->GetYaxis()->SetRangeUser(-50,500);
 GvsT3->Draw("colz");
	
 
 TCanvas*CanPH=new TCanvas("Pulse heights","Pulse heights",200,10,1200,800);  
 CanPH->Divide(3,2);
 
 for(int j=0;j<5;j++)
   {
    CanPH->cd(j+1);  
    gPad->SetLogy();
    AllPH[j]->Draw();
   }
 CanPH->Print(Form("Pulses_%s_%s.pdf(",file[0].c_str(),coinc[0].c_str()),Form("Title:Pulse heights All"));    
 for(int j=0;j<5;j++)
   {
    CanPH->cd(j+1);  
    gPad->SetLogy(0);
	gPad->SetLeftMargin(0.15);
	AllPH[j]->SetLineColor(kBlack);
	AllPH[j]->SetLineWidth(2);
    //AllPH[j]->GetXaxis()->SetRangeUser(1.,500);
    AllPH[j]->GetYaxis()->SetTitleOffset(1.9);
    gStyle->SetOptStat(0);
    AllPH[j]->Draw();
   }
 CanPH->Print(Form("Pulses_%s_%s.pdf",file[0].c_str(),coinc[0].c_str()),Form("Title:Pulse heights"));    
  

 TCanvas**CanPHi=new TCanvas*[5];

 for(int j=0;j<5;j++)
   {
  	CanPHi[j]=new TCanvas(Form("Pulse heights %s",sPH[j].c_str()),Form("Pulse heights %s",sPH[j].c_str()),200,10,1200,800);  
    CanPHi[j]->cd();
	gPad->SetLeftMargin(0.15);
    AllPH[j]->Draw();
    AllPH[j]->GetXaxis()->SetRangeUser(-1.,4096);
    CanPHi[j]->Print(Form("Pulses_%s_%s.pdf",file[0].c_str(),coinc[0].c_str()),Form("Title:Pulse heights %s",sPH[j].c_str()));    
    AllPH[j]->GetXaxis()->SetRangeUser(1.,500);
    CanPHi[j]->Print(Form("Pulses_%s_%s.pdf",file[0].c_str(),coinc[0].c_str()),Form("Title:Pulse heights %s Zoom",sPH[j].c_str()));    
    gPad->SetLogy(1);
    CanPHi[j]->Print(Form("Pulses_%s_%s.pdf",file[0].c_str(),coinc[0].c_str()),Form("Title:Pulse heights %s ZoomLog",sPH[j].c_str()));    
   }	   
	
	
	
  
 TCanvas*CanEff=new TCanvas("CanEff","CanEff",200,10,1300,800);  
 CanEff->cd();
 CanEff->SetBottomMargin(0.27);
 gPad->SetGridy(1);

 TEff->SetLineWidth(0);
 TEff->SetMarkerStyle(20);
 TEff->SetMarkerSize(0);
 TEff->Scale(100./nentries);
 TEff->SetBarWidth(0.8);
 TEff->SetBarOffset(0.1);
 TEff->GetYaxis()->SetRangeUser(0,110);
 TEff->GetYaxis()->SetNdivisions(511,0);
 TEff->GetYaxis()->SetTitle("Normalized rates");
 TEff->LabelsOption("v","X");
 TEff->Draw("bar");
 
 TText** txt=new TText*[nx];
 
 for(int i=0;i<nx;i++)
   {
    if(TEff->GetBinContent(i+1)>=100&&i==0) txt[i]=new TText(i+0.2,TEff->GetBinContent(i+1)+1,Form("%d",(int)TEff->GetBinContent(i+1))); 
    else if(TEff->GetBinContent(i+1)>=100) txt[i]=new TText(i+0.1,TEff->GetBinContent(i+1)+1,Form("%3.2f",TEff->GetBinContent(i+1))); 
    else if(TEff->GetBinContent(i+1)>=10) txt[i]=new TText(i+0.15,TEff->GetBinContent(i+1)+1,Form("%3.2f",TEff->GetBinContent(i+1))); 
    else txt[i]=new TText(i+0.25,TEff->GetBinContent(i+1)+1,Form("%3.2f",TEff->GetBinContent(i+1))); 
    txt[i]->SetTextSize(0.02);    
    txt[i]->Draw();    
   } 
 CanEff->Print(Form("Pulses_%s_%s.pdf)",file[0].c_str(),coinc[0].c_str()),Form("Title:Rates"));    

 
 
 
}



void RealRates()
{
 int Nf=4;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[4]={"NL0086","NL0087","NL0088","NL0089"};
 string coinc[4]={"T1&T4","T1&T4","T1&T4","T1&T3&T4"};

 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=20;
 string T[nx]={" All","T1","T2","T3","T4","G","T1T4","T1T4 & hits","T1T2","T1T2 & hits","T1T3","T1T3 & hits","T1T2T3","T1T2T3 & hits","T1T3T4","T1T3T4 & hits","T1T2T3T4","T1T2T3T4 & hits","T1T2T3T4 & NoG","T1T2T3T4 & NoG & hits"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];
 
 for(int i=0;i<Nf;i++)
   {
    cout << Form("/data/psmangeard/AESOPLite/Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("../Data/%s.BPD.EVENT.root",file[i].c_str()));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
     
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);
       //T1     
       if(e->get_EneT1().at(0)>0)TEff[i]->Fill(T[1].c_str(),1);
       //T2
       if(e->get_EneT2().at(0)>0)TEff[i]->Fill(T[2].c_str(),1);
       //T3
       if(e->get_EneT3().at(0)>0)TEff[i]->Fill(T[3].c_str(),1);
       //T4
       if(e->get_EneT4().at(0)>0)TEff[i]->Fill(T[4].c_str(),1);
       //G
       if(e->get_Eneg().at(0)>0)TEff[i]->Fill(T[5].c_str(),1);
       
       //T1&T4     
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0)TEff[i]->Fill(T[6].c_str(),1);
       //T1&T4&hits     
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0)TEff[i]->Fill(T[7].c_str(),1);
       
       //T1&T2     
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0)TEff[i]->Fill(T[8].c_str(),1);
       //T1&T2&hits     
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&nnhits>0)TEff[i]->Fill(T[9].c_str(),1);
       
        //T1&T3     
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0)TEff[i]->Fill(T[10].c_str(),1);
       //T1&T3&hits     
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&nnhits>0)TEff[i]->Fill(T[11].c_str(),1);

       //"T1&T2&T3"
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0)TEff[i]->Fill(T[12].c_str(),1);
       //"T1&T2&T3&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0)TEff[i]->Fill(T[13].c_str(),1);
       
       //"T1&T3&T4"
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0)TEff[i]->Fill(T[14].c_str(),1);
       //"T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0)TEff[i]->Fill(T[15].c_str(),1);

       //"T1&T2&T3&T4"
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0)TEff[i]->Fill(T[16].c_str(),1);
       //"T1&T2&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0)TEff[i]->Fill(T[17].c_str(),1);

       //"T1&T2&T3&T4&NoG"
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0)TEff[i]->Fill(T[18].c_str(),1);
       //"T1&T2&T3&T4&NoG&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0&&nnhits>0)TEff[i]->Fill(T[19].c_str(),1);
 
    }   //j
  }//i

   
   
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 string PDFFile="GroundRuns_Rates.pdf";

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 
 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.27);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      } 
      
    if(i==0) Can[i]->Print(Form("%s(",PDFFile.c_str()),Form("Title:Rates"));    
    else if(i==Nf-1) Can[i]->Print(Form("%s)",PDFFile.c_str()),Form("Title:Rates"));    
    else Can[i]->Print(Form("%s",PDFFile.c_str()),Form("Title:Rates"));    
   }
  
  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
   } //i

      
 
 
}



void RecoEff()
{

 int Nf=4;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[4]={"NL0086","NL0087","NL0088","NL0089"};
 string coinc[4]={"T1&T4","T1&T4","T1&T4","T1&T3&T4"};

 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=11;
 string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 6+ layers","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];


//////////////////////////// 
//Occupancy histograms
//////////////////////////// 
    
 TH1F*** LOcc=new TH1F**[Nf];
 TH1F*** ClusOcc=new TH1F**[Nf];
 TH1F**** ChipOcc=new TH1F***[Nf];
 TH1F**** fStripOcc=new TH1F***[Nf];
 TH1F**** StripOcc=new TH1F***[Nf];
 TH1F**** CoordOcc=new TH1F***[Nf];
 
 for(int h=0;h<Nf;h++)
   { 
     
    LOcc[h]=new TH1F*[nx];
    ClusOcc[h]=new TH1F*[nx];
    ChipOcc[h]=new TH1F**[nx];
    fStripOcc[h]=new TH1F**[nx];
    StripOcc[h]=new TH1F**[nx];
    CoordOcc[h]=new TH1F**[nx];
 
    for(int i=0;i<nx;i++)
      {  
       LOcc[h][i]=new TH1F(Form("Run: %s, Layer occupancy, %s",file[h].c_str(),T[i].c_str()),Form("Run: %s, Layer occupancy, %s",file[h].c_str(),T[i].c_str()),9,0,9);
       LOcc[h][i]->GetXaxis()->SetTitle("Layer (0 to 6)");
       LOcc[h][i]->GetYaxis()->SetTitle("Event number");
       LOcc[h][i]->SetLineColor(i+1);
       LOcc[h][i]->SetLineWidth(1);  
       LOcc[h][i]->SetMinimum(0);  

       ClusOcc[h][i]=new TH1F(Form("Run: %s, Cluster occupancy, %s",file[h].c_str(),T[i].c_str()),Form("Run: %s, Cluster occupancy, %s",file[h].c_str(),T[i].c_str()),9,0,9);
       ClusOcc[h][i]->GetXaxis()->SetTitle("Layer (0 to 6)");
       ClusOcc[h][i]->GetYaxis()->SetTitle("Cluster number");
       ClusOcc[h][i]->SetLineColor(i+1);
       ClusOcc[h][i]->SetLineWidth(1);  
       ClusOcc[h][i]->SetMinimum(0);  
    
       ChipOcc[h][i]=new TH1F*[7];
       fStripOcc[h][i]=new TH1F*[7];
       StripOcc[h][i]=new TH1F*[7];
       CoordOcc[h][i]=new TH1F*[7];
       for(int j=0;j<7;j++)
         {    
          ChipOcc[h][i][j]=new TH1F(Form("Run: %s, Chip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, Chip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),15,0,15);
          ChipOcc[h][i][j]->GetXaxis()->SetTitle("Chip ID");
          ChipOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          ChipOcc[h][i][j]->SetLineColor(i+1);
          ChipOcc[h][i][j]->SetLineWidth(1);
          ChipOcc[h][i][j]->GetXaxis()->SetRangeUser(0,15);  
          ChipOcc[h][i][j]->SetMinimum(0);  
          fStripOcc[h][i][j]=new TH1F(Form("Run: %s, First strip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, First strip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),769,0,769);
          fStripOcc[h][i][j]->GetXaxis()->SetTitle("First strip ID");
          fStripOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          fStripOcc[h][i][j]->SetLineColor(i+1);
          fStripOcc[h][i][j]->SetLineWidth(1);
          fStripOcc[h][i][j]->GetXaxis()->SetRangeUser(0,769);  
          fStripOcc[h][i][j]->SetMinimum(0);  
          StripOcc[h][i][j]=new TH1F(Form("Run: %s, Strip L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, Strip L%d, %s",file[h].c_str(),j,T[i].c_str()),769,0,769);
          StripOcc[h][i][j]->GetXaxis()->SetTitle("Strip");
          StripOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          StripOcc[h][i][j]->SetLineColor(i+1);
          StripOcc[h][i][j]->SetLineWidth(1);
          StripOcc[h][i][j]->GetXaxis()->SetRangeUser(0,769);  
          StripOcc[h][i][j]->SetMinimum(0);  
          CoordOcc[h][i][j]=new TH1F(Form("Run: %s, Cluster coord. L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, Cluster coord. L%d, %s",file[h].c_str(),j,T[i].c_str()),200,-100,100);
          if(j==0||j==4||j==6)CoordOcc[h][i][j]->GetXaxis()->SetTitle("X (mm)");
          else CoordOcc[h][i][j]->GetXaxis()->SetTitle("Y (mm)");
          CoordOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          CoordOcc[h][i][j]->SetLineColor(i+1);
          CoordOcc[h][i][j]->SetLineWidth(1);
          CoordOcc[h][i][j]->GetXaxis()->SetRangeUser(-100,100);    
          CoordOcc[h][i][j]->SetMinimum(0);
         }//j
      }//i
   }//h


 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/%s.BPD.EVENT.root",file[i].c_str()));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
       
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
       //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}
       
       int*Ltmp=new int[7];
       for(int k=0;k<7;k++) Ltmp[k]=0;
       
       for(int k=0;k<nnhits;k++)
         {
          int L=e->get_hits().at(k)->get_L();
          int fStripID=e->get_hits().at(k)->get_fstripID();
          int nstrips=e->get_hits().at(k)->get_nstrips();
          int chip=e->get_hits().at(k)->get_chip();
          float x=e->get_hits().at(k)->get_x();
          float y=e->get_hits().at(k)->get_y();
          float coord=y;//B plan per default
          if(L==0||L==4||L==6)  coord=x;//NB plane
          Ltmp[L]++;
          for(int ij=1;ij<nx;ij++)
           {
            if(Ncase[ij]==0) continue;
            if(ij==nx-1&&e->get_hits().at(k)->get_flagPR()==0) continue;
            fStripOcc[i][ij][L]->Fill(fStripID);
            StripOcc[i][ij][L]->Fill(fStripID);
            for(int ijk=1;ijk<=nstrips;ijk++) StripOcc[i][ij][L]->Fill(fStripID+ijk);
            ChipOcc[i][ij][L]->Fill(chip);
            CoordOcc[i][ij][L]->Fill(10*coord);
            ClusOcc[i][ij]->Fill(L);
            if(Ltmp[L]==1) LOcc[i][ij]->Fill(L);
           }               
         }  //k
      } //j
   }//i

 
  
  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("/data/psmangeard/AESOPLite/Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 //gStyle->SetOptStat("emruo");
 gStyle->SetOptStat(0);
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 string PDFFile="GroundRuns_RecoEff.pdf";

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 
 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.29);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      } 
      
    //if(i==0) Can[i]->Print(Form("%s(",PDFFile.c_str()),Form("Title:Rates"));    
   // else if(i==Nf-1) Can[i]->Print(Form("%s)",PDFFile.c_str()),Form("Title:Rates"));    
    //else Can[i]->Print(Form("%s",PDFFile.c_str()),Form("Title:Rates"));    
    if (i==Nf-1)Can[i]->Print(Form("Run%s_Tckrates.pdf",file[i].c_str()),Form("Title:Rates"));    
   }

//////////////////////////////    
//Strip Occupancy
//////////////////////////////    
 TCanvas**CanStrip=new TCanvas*[Nf]; 
  
 for(int i=0;i<Nf;i++)
   {
    CanStrip[i]=new TCanvas(Form("Strip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Strip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanStrip[i]->Divide(4,2);
    for(int j=0;j<7;j++)
      {  
       CanStrip[i]->cd(j+1);
       gPad->SetBottomMargin(0.1);
       StripOcc[i][2][j]->SetLineColor(kBlack);
       StripOcc[i][9][j]->SetLineColor(kBlue);
       StripOcc[i][10][j]->SetLineColor(kRed);
       StripOcc[i][2][j]->SetTitle(Form("Strip Occupancy, Run %s, Layer %d ",file[i].c_str(),j));
       StripOcc[i][2][j]->GetYaxis()->SetTitleOffset(1.5);
       StripOcc[i][2][j]->Draw("hist");
       
       for(int k=2;k<nx;k++)
         {  
          if(k!=9 && k!=10) continue;
          StripOcc[i][k][j]->Draw("same");
         }//k     
      }//j 
   } //  i

//////////////////////////////    
//Chip Occupancy
//////////////////////////////  
 TCanvas**CanChip=new TCanvas*[Nf]; 
  
 for(int i=0;i<Nf;i++)
   {
    CanChip[i]=new TCanvas(Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    CanChip[i]->Divide(4,2);
    for(int j=0;j<7;j++)
      {  
       CanChip[i]->cd(j+1);
       gPad->SetBottomMargin(0.1);
       ChipOcc[i][2][j]->SetLineColor(kBlack);
       ChipOcc[i][9][j]->SetLineColor(kBlue);
       ChipOcc[i][10][j]->SetLineColor(kRed);
       ChipOcc[i][2][j]->GetYaxis()->SetTitleOffset(1.5);
       ChipOcc[i][2][j]->SetTitle(Form("Chip Occupancy, Run %s, Layer %d ",file[i].c_str(),j));
       ChipOcc[i][2][j]->Draw("hist");
       for(int k=2;k<nx;k++)
         {  
          if(k!=9 && k!=10) continue;
          ChipOcc[i][k][j]->Draw("same");
         }//k     
      }//j 
   } //  i
   
   
}//end function

void RecoEffOne(string s,string c)
{

 int Nf=1;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;
 
 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=15;
 string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 6+ layers","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track","T1T3T4 & L0L4L6","T1T3T4 & L1L3L5","T1T3T4 & (L0L4L6&L1L3L5)","T1T3T4 & (L0L4L6 | L1L3L5)"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];
 TH1F** hChi2B=new TH1F*[Nf];
 TH1F** hChi2NB=new TH1F*[Nf];

//////////////////////////// 
//Occupancy histograms
//////////////////////////// 
    
 TH1F*** LOcc=new TH1F**[Nf];
 TH1F*** ClusOcc=new TH1F**[Nf];
 TH1F**** ChipOcc=new TH1F***[Nf];
 TH1F**** fStripOcc=new TH1F***[Nf];
 TH1F**** StripOcc=new TH1F***[Nf];
 TH1F**** CoordOcc=new TH1F***[Nf];

 for(int h=0;h<Nf;h++)
   { 
     
    LOcc[h]=new TH1F*[nx];
    ClusOcc[h]=new TH1F*[nx];
    ChipOcc[h]=new TH1F**[nx];
    fStripOcc[h]=new TH1F**[nx];
    StripOcc[h]=new TH1F**[nx];
    CoordOcc[h]=new TH1F**[nx];

    for(int i=0;i<nx;i++)
      {  
       LOcc[h][i]=new TH1F(Form("Run: %s, Layer occupancy, %s",file[h].c_str(),T[i].c_str()),Form("Run: %s, Layer occupancy, %s",file[h].c_str(),T[i].c_str()),9,0,9);
       LOcc[h][i]->GetXaxis()->SetTitle("Layer (0 to 6)");
       LOcc[h][i]->GetYaxis()->SetTitle("Event number");
       LOcc[h][i]->SetLineColor(i+1);
       LOcc[h][i]->SetLineWidth(1);  
       LOcc[h][i]->SetMinimum(0);  

       ClusOcc[h][i]=new TH1F(Form("Run: %s, Cluster occupancy, %s",file[h].c_str(),T[i].c_str()),Form("Run: %s, Cluster occupancy, %s",file[h].c_str(),T[i].c_str()),9,0,9);
       ClusOcc[h][i]->GetXaxis()->SetTitle("Layer (0 to 6)");
       ClusOcc[h][i]->GetYaxis()->SetTitle("Cluster number");
       ClusOcc[h][i]->SetLineColor(i+1);
       ClusOcc[h][i]->SetLineWidth(1);  
       ClusOcc[h][i]->SetMinimum(0);  
    
       ChipOcc[h][i]=new TH1F*[7];
       fStripOcc[h][i]=new TH1F*[7];
       StripOcc[h][i]=new TH1F*[7];
       CoordOcc[h][i]=new TH1F*[7];
       for(int j=0;j<7;j++)
         {    
          ChipOcc[h][i][j]=new TH1F(Form("Run: %s, Chip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, Chip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),15,0,15);
          ChipOcc[h][i][j]->GetXaxis()->SetTitle("Chip ID");
          ChipOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          ChipOcc[h][i][j]->SetLineColor(i+1);
          ChipOcc[h][i][j]->SetLineWidth(1);
          ChipOcc[h][i][j]->GetXaxis()->SetRangeUser(0,15);  
          ChipOcc[h][i][j]->SetMinimum(0);  
          fStripOcc[h][i][j]=new TH1F(Form("Run: %s, First strip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, First strip ID L%d, %s",file[h].c_str(),j,T[i].c_str()),769,0,769);
          fStripOcc[h][i][j]->GetXaxis()->SetTitle("First strip ID");
          fStripOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          fStripOcc[h][i][j]->SetLineColor(i+1);
          fStripOcc[h][i][j]->SetLineWidth(1);
          fStripOcc[h][i][j]->GetXaxis()->SetRangeUser(0,769);  
          fStripOcc[h][i][j]->SetMinimum(0);  
          StripOcc[h][i][j]=new TH1F(Form("Run: %s, Strip L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, Strip L%d, %s",file[h].c_str(),j,T[i].c_str()),769,0,769);
          StripOcc[h][i][j]->GetXaxis()->SetTitle("Strip");
          StripOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          StripOcc[h][i][j]->SetLineColor(i+1);
          StripOcc[h][i][j]->SetLineWidth(1);
          StripOcc[h][i][j]->GetXaxis()->SetRangeUser(0,769);  
          StripOcc[h][i][j]->SetMinimum(0);  
          CoordOcc[h][i][j]=new TH1F(Form("Run: %s, Cluster coord. L%d, %s",file[h].c_str(),j,T[i].c_str()),Form("Run: %s, Cluster coord. L%d, %s",file[h].c_str(),j,T[i].c_str()),200,-100,100);
          if(j==0||j==4||j==6)CoordOcc[h][i][j]->GetXaxis()->SetTitle("X (mm)");
          else CoordOcc[h][i][j]->GetXaxis()->SetTitle("Y (mm)");
          CoordOcc[h][i][j]->GetYaxis()->SetTitle("Entries");
          CoordOcc[h][i][j]->SetLineColor(i+1);
          CoordOcc[h][i][j]->SetLineWidth(1);
          CoordOcc[h][i][j]->GetXaxis()->SetRangeUser(-100,100);    
          CoordOcc[h][i][j]->SetMinimum(0);
         }//j
      }//i
   }//h


 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
     //Chi2 histograms
   
    hChi2NB[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),20,0,.01);
    hChi2NB[i]->SetStats(1);
    
     hChi2B[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),20,0,.01);
     hChi2B[i]->SetStats(1);
   
    
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/%s.EVENT.root",file[i].c_str()));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int* Lay=new int[7];
       for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);

       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
       
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
       //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>0&&e->get_chi2BPR()>0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}
        //T1&T3&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T3&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
       //T1&T3&T4& L0L1L3||L4L5L6
       if(Ncase[11]==1&&Ncase[12]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T3&T4& L0L1L3L4L5L6
       if(Ncase[11]==1||Ncase[12]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
      
       int**Ltmp=new int*[7];
       for(int k=0;k<7;k++) 
         {
          Ltmp[k]=new int[nx];
          for(int kl=0;kl<nx;kl++) Ltmp[k][kl]=0;
         } 
       for(int k=0;k<nnhits;k++)
         {
          int L=e->get_hits().at(k)->get_L();
          int fStripID=e->get_hits().at(k)->get_fstripID();
          int nstrips=e->get_hits().at(k)->get_nstrips();
          int chip=e->get_hits().at(k)->get_chip();
          float x=e->get_hits().at(k)->get_x();
          float y=e->get_hits().at(k)->get_y();
          float coord=y;//B plan per default
          if(L==0||L==4||L==6)  coord=x;//NB plane
          
          for(int ij=1;ij<nx;ij++)
           {
            if(Ncase[ij]==0) continue;
            if(ij==10&&e->get_hits().at(k)->get_flagPR()<=0) continue;
            fStripOcc[i][ij][L]->Fill(fStripID);
            StripOcc[i][ij][L]->Fill(fStripID);
            for(int ijk=1;ijk<=nstrips;ijk++) StripOcc[i][ij][L]->Fill(fStripID+ijk);
            ChipOcc[i][ij][L]->Fill(chip);
            CoordOcc[i][ij][L]->Fill(10*coord);
            ClusOcc[i][ij]->Fill(L);
            Ltmp[L][ij]++;
            if(Ltmp[L][ij]==1) LOcc[i][ij]->Fill(L);
           }               
         }  //k
       if(Ncase[10]==1)
        {
         hChi2B[i]->Fill(e->get_chi2BPR());   
         hChi2NB[i]->Fill(e->get_chi2NBPR());  
        }   
      } //j
   }//i

 
  
  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 //gStyle->SetOptStat("emruo");
 gStyle->SetOptStat(0);
 gStyle->SetPalette( kDarkBodyRadiator);
 TGaxis::SetMaxDigits(4);

//////////////////////////////    
 

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 
 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.32);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      }    
 
    if (i==Nf-1)Can[i]->Print(Form("Efficiency_%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }

//////////////////////////////    
//Layer Occupancy
//////////////////////////////    
  TCanvas**CanLayer=new TCanvas*[Nf]; 
  TLegend*LegLayOcc=new TLegend(0.,0.4,1.,0.6); 
  LegLayOcc->SetBorderSize(0);
  LegLayOcc->SetFillColor(0);
  for(int i=0;i<Nf;i++)
    {
     CanLayer[i]=new TCanvas(Form("Layer Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Layer Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
     CanLayer[i]->Divide(2,1);
     CanLayer[i]->cd(1);
     gPad->SetBottomMargin(0.1);    
     
     LOcc[i][2]->SetLineColor(kBlack);
     LOcc[i][9]->SetLineColor(kBlue);
     LOcc[i][10]->SetLineColor(kRed);
     LegLayOcc->AddEntry(LOcc[i][2],Form("%s",T[2].c_str()),"l");
     LegLayOcc->AddEntry(LOcc[i][9],Form("%s",T[9].c_str()),"l");
     LegLayOcc->AddEntry(LOcc[i][10],Form("%s (used hits only)",T[10].c_str()),"l");
   
     LOcc[i][2]->SetMaximum(1.2*LOcc[i][2]->GetEntries()/7);
     LOcc[i][2]->SetTitle(Form("Layer Occupancy, Run %s",file[i].c_str()));
     LOcc[i][2]->GetYaxis()->SetTitleOffset(1.5);
     LOcc[i][2]->Draw("hist");
        
     for(int k=2;k<nx;k++)
       {  
        if(k!=9 && k!=10) continue;
        LOcc[i][k]->Draw("same");
       }//k     
     
     CanLayer[i]->cd(2);
     LegLayOcc->Draw();
     if (i==Nf-1)CanLayer[i]->Print(Form("Efficiency_%s_%s.pdf",file[i].c_str(),coinc[i].c_str()),Form("Title:Layer Occupancy"));    
    } //  i  
  
//////////////////////////////    
//Strip Occupancy
//////////////////////////////    
 TCanvas**CanStrip=new TCanvas*[Nf]; 
 TLegend*LegStripOcc=new TLegend(0.,0.4,1.,0.6); 
 LegStripOcc->SetBorderSize(0);
 LegStripOcc->SetFillColor(0);
 for(int i=0;i<Nf;i++)
   {
    CanStrip[i]=new TCanvas(Form("Strip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Strip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanStrip[i]->Divide(4,2);
    for(int j=0;j<7;j++)
      {  
       CanStrip[i]->cd(j+1);
       gPad->SetBottomMargin(0.1);
       StripOcc[i][2][j]->SetLineColor(kBlack);
       StripOcc[i][9][j]->SetLineColor(kBlue);
       StripOcc[i][10][j]->SetLineColor(kRed);
       if(j==0)
        {
         LegStripOcc->AddEntry(StripOcc[i][2][j],Form("%s",T[2].c_str()),"l");
         LegStripOcc->AddEntry(StripOcc[i][9][j],Form("%s",T[9].c_str()),"l");
         LegStripOcc->AddEntry(StripOcc[i][10][j],Form("%s (used hits only)",T[10].c_str()),"l");
        }
       StripOcc[i][2][j]->SetMaximum(3.4*StripOcc[i][2][j]->GetEntries()/768.);
       if(j>2)StripOcc[i][2][j]->SetMaximum(2.*StripOcc[i][2][j]->GetEntries()/768.);
       StripOcc[i][2][j]->SetTitle(Form("Strip Occupancy, Run %s, Layer %d ",file[i].c_str(),j));
       StripOcc[i][2][j]->GetYaxis()->SetTitleOffset(1.5);
       StripOcc[i][2][j]->Draw("hist");
       
       for(int k=2;k<nx;k++)
         {  
          if(k!=9 && k!=10) continue;
          StripOcc[i][k][j]->Draw("same");
         }//k     
      }//j 
     CanStrip[i]->cd(8);
     LegStripOcc->Draw();
    if (i==Nf-1)CanStrip[i]->Print(Form("Efficiency_%s_%s.pdf",file[i].c_str(),coinc[i].c_str()),Form("Title:Strip Occupancy"));    
  
   } //  i

//////////////////////////////    
//Chip Occupancy
//////////////////////////////  
 TCanvas**CanChip=new TCanvas*[Nf]; 
 TLegend*LegChipOcc=new TLegend(0.,0.4,1.,0.6); 
 LegChipOcc->SetBorderSize(0);
 LegChipOcc->SetFillColor(0);
  
 for(int i=0;i<Nf;i++)
   {
    CanChip[i]=new TCanvas(Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanChip[i]->Divide(4,2);
    for(int j=0;j<7;j++)
      {  
       CanChip[i]->cd(j+1);
       gPad->SetBottomMargin(0.1);
       
       ChipOcc[i][2][j]->SetLineColor(kBlack);
       ChipOcc[i][9][j]->SetLineColor(kBlue);
       ChipOcc[i][10][j]->SetLineColor(kRed);
       if(j==0)
        {
         LegChipOcc->AddEntry(StripOcc[i][2][j],Form("%s",T[2].c_str()),"l");
         LegChipOcc->AddEntry(StripOcc[i][9][j],Form("%s",T[9].c_str()),"l");
         LegChipOcc->AddEntry(StripOcc[i][10][j],Form("%s (used hits only)",T[10].c_str()),"l");
        }

       ChipOcc[i][2][j]->GetYaxis()->SetTitleOffset(1.5);
       ChipOcc[i][2][j]->SetMaximum(3.*ChipOcc[i][2][j]->GetEntries()/12.);
       if(j>2)ChipOcc[i][2][j]->SetMaximum(2.*ChipOcc[i][2][j]->GetEntries()/12.);
       ChipOcc[i][2][j]->SetTitle(Form("Chip Occupancy, Run %s, Layer %d ",file[i].c_str(),j));
       ChipOcc[i][2][j]->Draw("hist");
       for(int k=2;k<nx;k++)
         {  
          if(k!=9 && k!=10) continue;
          ChipOcc[i][k][j]->Draw("same");
         }//k     
      }//j 
     CanChip[i]->cd(8);
     LegChipOcc->Draw();
    CanChip[i]->Print(Form("Efficiency_%s_%s.pdf",file[i].c_str(),coinc[i].c_str()),Form("Title:Chip Occupancy"));    
   } //  i

//////////////////////////////    
//Chi2 distribution
////////////////////////////// 
 gStyle->SetOptStat(1);

 TCanvas**CanChi2=new TCanvas*[Nf]; 

  
 for(int i=0;i<Nf;i++)
   {
    CanChi2[i]=new TCanvas(Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanChi2[i]->Divide(2,1);
    CanChi2[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    hChi2NB[i]->GetXaxis()->SetTitle("Chi2 in Non-Bending plane");
    //hChi2NB[i]->GetYaxis()->SetTitle("arbitrary unit");
    hChi2NB[i]->Draw("hist");       
    CanChi2[i]->cd(2);
    gPad->SetBottomMargin(0.1);
    hChi2B[i]->Draw("hist");       
    hChi2B[i]->GetXaxis()->SetTitle("Chi2 in Bending plane");
    //hChi2B[i]->GetYaxis()->SetTitle("arbitrary unit");
    if (i==Nf-1)CanChi2[i]->Print(Form("Efficiency_%s_%s.pdf)",file[i].c_str(),coinc[i].c_str()),Form("Title:Chi2 distribution"));    
   } //  i   
   
   
}//end function


void Deflection(string s,string c,int geoconf)
{
 
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

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


 int Nf=1;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;
 
 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=21;
 //string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 5+ layers & track","T1T3T4 & 6+ layers","T1T3T4 & 6+ layers & track","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track","T1T3T4 & L0L4L6","T1T3T4 & L1L3L5","T1T3T4 & (L0L4L6&L1L3L5)","T1T3T4 & (L0L4L6 | L1L3L5)"};
 
 string T[nx]={"All","T1T4","T1T4 & hits","T1T4 & 1+ layers","T1T4 & 2+ layers","T1T4 & 3+ layers","T1T4 & 4+ layers","T1T4 & 5+ layers","T1T4 & 5+ layers & track","T1T4 & 6+ layers","T1T4 & 6+ layers & track","T1T4 & 7 layers","T1T4 & 7 layers & track","T1T4 & L0L4L6 (NB)","T1T4 & L1L2L3 (B)","T1T4 & (L0L4L6&L1L2L3)","T1T4 & (L0L4L6 | L1L2L3)","T1T4&Pattern 0","T1T4&Pattern 1 (B only)","T1T4&Pattern 2 (NB only)","T1T4&Pattern 3 (B&NB)"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];
 TH1F** hChi2B=new TH1F*[Nf];
 TH1F** hChi2NB=new TH1F*[Nf];

 //Deflection histograms
 TH1F**HDef=new TH1F*[Nf];
 TH1F**HDef2=new TH1F*[Nf];

 TH1F**HDefT2=new TH1F*[Nf];
 TH1F**HDefnoT2=new TH1F*[Nf];
 
 TH2F**HDefnoT2vsT4=new TH2F*[Nf];

 TH2F***HDdefvsstripID=new TH2F**[Nf];
 
 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
     //Chi2 histograms
   
    hChi2NB[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),20,0,20);
    hChi2NB[i]->SetStats(1);
    
    hChi2B[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),20,0,20);
    hChi2B[i]->SetStats(1);
   
    //Deflection histrograms
    HDef[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection",file[i].c_str(),coinc[i].c_str()),1000,-1,1);
    HDef[i]->SetStats(1);

    HDef2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection2",file[i].c_str(),coinc[i].c_str()),1800,TMath::Pi(),TMath::Pi());
    HDef2[i]->SetStats(1);
 
	   
    HDefT2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection with T2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection with T2",file[i].c_str(),coinc[i].c_str()),1800,TMath::Pi(),TMath::Pi());
    HDefT2[i]->SetStats(1);  
    HDefnoT2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection without T2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection without T2",file[i].c_str(),coinc[i].c_str()),1800,TMath::Pi(),TMath::Pi());
    HDefnoT2[i]->SetStats(1);  
	
    HDefnoT2vsT4[i]= new TH2F(Form("Coincidence %s, Deflection without T2 vs T4",coinc[i].c_str()),Form("Coincidence %s, Deflection without T2 vs T4",coinc[i].c_str()),2000,0,2000,900,TMath::Pi(),TMath::Pi());
    HDefnoT2vsT4[i]->SetStats(1);  
	     
    HDdefvsstripID[i]=new TH2F*[7];
    for(int j=0;j<7;j++)
      {
       HDdefvsstripID[i][j]=new TH2F(Form("Run %s, Coincidence %s, Deflection vs Chip (L%d)",file[i].c_str(),coinc[i].c_str(),j),Form("Run %s, Coincidence %s, Deflection vs Chip (L%d)",file[i].c_str(),coinc[i].c_str(),j),12,0,12,1800,-TMath::Pi(),TMath::Pi());
       HDdefvsstripID[i][j]->SetStats(1);  
      } 
    
    chain[i]=new TChain("Data");
   // chain[i]->Add(Form("../Data/%s.EVENT.root",file[i].c_str()));
    //Esrange
    chain[i]->Add(Form("../Data/NL2084.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2085.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2088.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2089.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2090.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2091.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2092.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2101.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2102.BPD.EVENT.root"));
    chain[i]->Add(Form("../Data/NL2103.BPD.EVENT.root"));
    //Palestine
    //chain[i]->Add(Form("../Data/NL0239.BPD.EVENT.root"));
   // chain[i]->Add(Form("../Data/NL0244.BPD.EVENT.root"));
    //chain[i]->Add(Form("../Data/NL0247.BPD.EVENT.root"));
   // chain[i]->Add(Form("../Data/NL0248.BPD.EVENT.root"));
   // chain[i]->Add(Form("../Data/NL0250.BPD.EVENT.root"));
    //chain[i]->Add(Form("../Data/NL0254.BPD.EVENT.root"));
  //  chain[i]->Add(Form("../Data/NL2004.BPD.EVENT.root"));
//chain[i]->Add(Form("../Data/NL2026.BPD.EVENT.root"));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
      // if(nnhits<7) continue;
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int* Lay=new int[7];
       for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);

       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
       
       
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
 /*
       //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
        //T1&T3&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T3&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T3&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T3&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
      */
      
       /*//T1&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
       //T1&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
     */
       
       //T1&T4 New internal trigger
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
       //T1&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T4& L1L2L3
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[2]==1&&Lay[3]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T4& L0L2L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T4& L0L1L2L3L4L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
       //Fill up Pattern of Online tracker internal trigger
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0 && e->get_PatternEVT()==0){TEff[i]->Fill(T[17].c_str(),1);Ncase[17]=1;}
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0 && e->get_PatternEVT()==1){TEff[i]->Fill(T[18].c_str(),1);Ncase[18]=1;}
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0 && e->get_PatternEVT()==2){TEff[i]->Fill(T[19].c_str(),1);Ncase[19]=1;}
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0 && e->get_PatternEVT()==3){TEff[i]->Fill(T[20].c_str(),1);Ncase[20]=1;}
       
       int**Ltmp=new int*[7];
       for(int k=0;k<7;k++) 
         {
          Ltmp[k]=new int[nx];
          for(int kl=0;kl<nx;kl++) Ltmp[k][kl]=0;
         }

       //Deflection
       ///////////////////////////////    
       //Extract PR track information
       ///////////////////////////////    
       if(nnhits<5)continue;
       if(NL<5)continue;
       
       //Non bending plane
       float p0=e->get_interPR();
       float p1=e->get_slopePR();
       TF1* fNB=new TF1("fNB","pol1",-20,20);
       if(p1==0) continue;
       fNB->FixParameter(0,-p0/p1);
       fNB->FixParameter(1,1./p1);

       //Bending plane //ax^2+bx+c
       float zz0=0;
       float a=e->get_aPR();
       float b=e->get_bPR();
       float c=e->get_cPR();
       TF1* fB= new TF1("fB","[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,40);
       fB->FixParameter(3,zz0);
       fB->FixParameter(0,c);
       fB->FixParameter(1,b);
       fB->FixParameter(2,a);
       float lim=zL[1];//z position of 2nd layer
       float aa=fB->Eval(lim);
       float diff=2*a*lim+2*a*zz0+b;
       float limo=zL[5];//z position of 6th layer
       float aaout=fB->Eval(limo);
       float diffout=2*a*limo+2*a*zz0+b;     
 
       //Deflection diff of slopes  
       float deflection=e->get_deflecPR();
       
       //Deflection2 diff of angles  
       float deflection2=TMath::ATan(diffout)-TMath::ATan(diff);
       if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)<0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0 )
         { 	
          for(int ijk=0;ijk<nnhits;ijk++)
            {
             if(e->get_hits().at(ijk)->get_flagPR()==true) HDdefvsstripID[i][(int)e->get_hits().at(ijk)->get_L()]->Fill(e->get_hits().at(ijk)->get_chip(),deflection2);
            
            } 
         }  
       if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0 )
        {  
         HDef[i]->Fill(deflection);
         HDef2[i]->Fill(deflection2);
        }
       
   	  if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0 )
        {  
         HDefT2[i]->Fill(deflection2);
			
        }		  
		  
	   if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)<=0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0&&e->get_chi2BPR()<=15&&e->get_chi2NBPR()<=15)
        {  
         HDefnoT2[i]->Fill(deflection2);
         HDefnoT2vsT4[i]->Fill(e->get_EneT4().at(0),deflection2);
        }	 

       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0) 
        {
         hChi2B[i]->Fill(e->get_chi2BPR());   
         hChi2NB[i]->Fill(e->get_chi2NBPR());  
        }   
      } //j
   }//i

  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
	
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 //gStyle->SetOptStat("emruo");
 gStyle->SetOptStat(1);
 //gStyle->SetPalette( kDarkBodyRadiator);
 gStyle->SetPalette( kTemperatureMap);
 TGaxis::SetMaxDigits(4);

//////////////////////////////    
 

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.32);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      }    
 
    if (i==Nf-1)Can[i]->Print(Form("Deflectionbis%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }


 
//////////////////////////////    
//Deflection
//////////////////////////////   
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(1);

 TF1**f=new TF1*[Nf];	
	
 TCanvas**CanDef=new TCanvas*[Nf]; 
 for(int i=0;i<Nf;i++)
   {
    CanDef[i]=new TCanvas(Form("Deflection, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Deflection, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanDef[i]->Divide(2,1);
    CanDef[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    HDefT2[i]->GetXaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
    HDefT2[i]->GetYaxis()->SetTitle("entries");
    HDefT2[i]->SetTitle("Coincidence T1.T2.T3.T4.noG");
    HDefT2[i]->GetYaxis()->SetTitleOffset(1.9);
    HDefT2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
    HDefT2[i]->Draw("hist");       
    f[i]= new TF1(Form("f%d",i),"gaus",-0.2,0.2);
    f[i]->SetLineColor(kBlack);
    HDefT2[i]->Fit(f[i],"R","",-0.02,0.02);  
    f[i]->DrawF1(-0.02,0.02,"same");
//HDefnoT2[i]->Draw("hist same");       
    CanDef[i]->cd(2);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    HDefnoT2[i]->Draw("hist");       
    HDefnoT2[i]->GetXaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
    HDefnoT2[i]->GetYaxis()->SetTitle("entries");
    HDefnoT2[i]->GetYaxis()->SetTitleOffset(1.9);
    HDefnoT2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
    HDefnoT2[i]->SetTitle("Coincidence T1.noT2.T3.T4.noG");
    CanDef[i]->Print(Form("Deflectionbis%s_%s.pdf",file[i].c_str(),coinc[i].c_str()),Form("Title:Deflection"));   
//   CanDef[i]->SaveAs("Deflection_newCKpressurewithfit.pdf");
   } //  i  
 
 TCanvas**canproton=new TCanvas*[Nf];  
 for(int i=0;i<Nf;i++)
   {
	canproton[i]=new TCanvas(Form("Deflection, Run %s, Coincidence %s  proton",file[i].c_str(),coinc[i].c_str()),Form("Deflection, Run %s, Coincidence %s  proton",file[i].c_str(),coinc[i].c_str()),200,10,1000,1000);  
    canproton[i]->cd();
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.15);  
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);  
	HDefnoT2vsT4[i]->GetXaxis()->SetTitle("Pulse height in T4");
	HDefnoT2vsT4[i]->GetYaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
	HDefnoT2vsT4[i]->GetYaxis()->SetTitleOffset(1.9);
	HDefnoT2vsT4[i]->GetYaxis()->SetRangeUser(-0.2,0.2);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
	HDefnoT2vsT4[i]->Draw("colz");
	   
   }
	
	
//////////////////////////////    
//Chi2 distribution
////////////////////////////// 

 TCanvas**CanChi2=new TCanvas*[Nf]; 

 TF1* chi2k1=new TF1("chi2k1","TMath::Power(x,-1./2.)*TMath::Exp(-x/2.)/(TMath::Power(2,1./2.)*TMath::Gamma(1./2.))",0,20); 
 for(int i=0;i<Nf;i++)
   {
    CanChi2[i]=new TCanvas(Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanChi2[i]->Divide(2,1);
    CanChi2[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    hChi2NB[i]->GetXaxis()->SetTitle("Chi2 in Non-Bending plane");
    //hChi2NB[i]->GetYaxis()->SetTitle("arbitrary unit");
    hChi2NB[i]->DrawNormalized("hist");
    chi2k1->Draw("same");
    CanChi2[i]->cd(2);
    gPad->SetBottomMargin(0.1);
    hChi2B[i]->DrawNormalized("hist");       
    hChi2B[i]->GetXaxis()->SetTitle("Chi2 in Bending plane");
    chi2k1->Draw("same");

    //hChi2B[i]->GetYaxis()->SetTitle("arbitrary unit");
    if (i==Nf-1)CanChi2[i]->Print(Form("Deflectionbis%s_%s.pdf)",file[i].c_str(),coinc[i].c_str()),Form("Title:Chi2 distribution"));    
   } //  i   
   
   
//////////////////////////////    
//StripID vs deflection distribution
//////////////////////////////  

 TCanvas***Canstrip=new TCanvas**[Nf]; 

 for(int i=0;i<Nf;i++)
   {
    Canstrip[i]=new TCanvas*[7];
    for(int j=0;j<7;j++)
      {
       Canstrip[i][j]=new TCanvas(Form("Deflection vs Chip, L%d, Run %s, Coincidence %s ",j,file[i].c_str(),coinc[i].c_str()),Form("Deflection vs Chips, L%d, Run %s, Coincidence %s ",j,file[i].c_str(),coinc[i].c_str()),200,10,1600,800);   
       Canstrip[i][j]->cd();
       HDdefvsstripID[i][j]->GetYaxis()->SetRangeUser(-0.2,0.2);
       HDdefvsstripID[i][j]->GetXaxis()->SetRangeUser(0,12);
       HDdefvsstripID[i][j]->GetXaxis()->SetTitle("Chip (0 to 11)");
       HDdefvsstripID[i][j]->GetYaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
       gPad->SetGridx(1);
       gPad->SetGridy(1);      
       HDdefvsstripID[i][j]->Draw("colz");
      }  
   }
   
}//end function


void DeflectionT1T2T3(string s,string c,int geoconf)
{
 
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

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


 int Nf=1;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;
 
 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=21;
 //string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 5+ layers & track","T1T3T4 & 6+ layers","T1T3T4 & 6+ layers & track","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track","T1T3T4 & L0L4L6","T1T3T4 & L1L3L5","T1T3T4 & (L0L4L6&L1L3L5)","T1T3T4 & (L0L4L6 | L1L3L5)"};
 
 string T[nx]={"All","T1T2T3","T1T2T3 & hits","T1T2T3 & 1+ layers","T1T2T3 & 2+ layers","T1T2T3 & 3+ layers","T1T2T3 & 4+ layers","T1T2T3 & 5+ layers","T1T2T3 & 5+ layers & track","T1T2T3 & 6+ layers","T1T2T3 & 6+ layers & track","T1T2T3 & 7 layers","T1T2T3 & 7 layers & track","T1T2T3 & L0L4L6 (NB)","T1T2T3 & L1L2L3 (B)","T1T2T3 & (L0L4L6&L1L2L3)","T1T2T3 & (L0L4L6 | L1L2L3)","T1T2T3&Pattern 0","T1T2T3&Pattern 1 (B only)","T1T2T3&Pattern 2 (NB only)","T1T2T3&Pattern 3 (B&NB)"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];
 TH1F** hChi2B=new TH1F*[Nf];
 TH1F** hChi2NB=new TH1F*[Nf];

 //Deflection histograms
 TH1F**HDef=new TH1F*[Nf];
 TH1F**HDef2=new TH1F*[Nf];

 TH1F**HDefT2=new TH1F*[Nf];
 TH1F**HDefnoT2=new TH1F*[Nf];
 
 TH2F**HDefnoT2vsT4=new TH2F*[Nf];

 TH2F***HDdefvsstripID=new TH2F**[Nf];
 
 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
     //Chi2 histograms
   
    hChi2NB[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),20,0,20);
    hChi2NB[i]->SetStats(1);
    
    hChi2B[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),20,0,20);
    hChi2B[i]->SetStats(1);
   
    //Deflection histrograms
    HDef[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection",file[i].c_str(),coinc[i].c_str()),1000,-1,1);
    HDef[i]->SetStats(1);

    HDef2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection2",file[i].c_str(),coinc[i].c_str()),1800,TMath::Pi(),TMath::Pi());
    HDef2[i]->SetStats(1);
 
	   
    HDefT2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection with T2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection with T2",file[i].c_str(),coinc[i].c_str()),14400,TMath::Pi(),TMath::Pi());
    HDefT2[i]->SetStats(1);  
    HDefnoT2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection without T2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection without T2",file[i].c_str(),coinc[i].c_str()),7200,TMath::Pi(),TMath::Pi());
    HDefnoT2[i]->SetStats(1);  
	
    HDefnoT2vsT4[i]= new TH2F(Form("Coincidence %s, Deflection without T2 vs T4",coinc[i].c_str()),Form("Coincidence %s, Deflection without T2 vs T4",coinc[i].c_str()),2000,0,2000,900,TMath::Pi(),TMath::Pi());
    HDefnoT2vsT4[i]->SetStats(1);  
	     
    HDdefvsstripID[i]=new TH2F*[7];
    for(int j=0;j<7;j++)
      {
       HDdefvsstripID[i][j]=new TH2F(Form("Run %s, Coincidence %s, Deflection vs Chip (L%d)",file[i].c_str(),coinc[i].c_str(),j),Form("Run %s, Coincidence %s, Deflection vs Chip (L%d)",file[i].c_str(),coinc[i].c_str(),j),12,0,12,1800,-TMath::Pi(),TMath::Pi());
       HDdefvsstripID[i][j]->SetStats(1);  
      } 
    
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("../Data/%s.EVENT.root",file[i].c_str()));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0086.BPD.EVENT.root"));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0088.BPD.EVENT.root"));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0089.BPD.EVENT.root"));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
      // if(nnhits<7) continue;
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int* Lay=new int[7];
       for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);

       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
       
       
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
 /*
       //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
        //T1&T3&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T3&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T3&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T3&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
      */
      
       /*//T1&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
       //T1&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
     */
       
       //T1&T4 New internal trigger
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
       //T1&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T4& L1L2L3
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[2]==1&&Lay[3]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T4& L0L2L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T4& L0L1L2L3L4L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
       //Fill up Pattern of Online tracker internal trigger
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&& e->get_PatternEVT()==0){TEff[i]->Fill(T[17].c_str(),1);Ncase[17]=1;}
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&& e->get_PatternEVT()==1){TEff[i]->Fill(T[18].c_str(),1);Ncase[18]=1;}
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&& e->get_PatternEVT()==2){TEff[i]->Fill(T[19].c_str(),1);Ncase[19]=1;}
       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&& e->get_PatternEVT()==3){TEff[i]->Fill(T[20].c_str(),1);Ncase[20]=1;}
       
       int**Ltmp=new int*[7];
       for(int k=0;k<7;k++) 
         {
          Ltmp[k]=new int[nx];
          for(int kl=0;kl<nx;kl++) Ltmp[k][kl]=0;
         }

       //Deflection
       ///////////////////////////////    
       //Extract PR track information
       ///////////////////////////////    
       if(nnhits<5)continue;
       
       //Non bending plane
       float p0=e->get_interPR();
       float p1=e->get_slopePR();
       TF1* fNB=new TF1("fNB","pol1",-20,20);
       if(p1==0) continue;
       fNB->FixParameter(0,-p0/p1);
       fNB->FixParameter(1,1./p1);

       //Bending plane //ax^2+bx+c
       float zz0=0;
       float a=e->get_aPR();
       float b=e->get_bPR();
       float c=e->get_cPR();
       TF1* fB= new TF1("fB","[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,40);
       fB->FixParameter(3,zz0);
       fB->FixParameter(0,c);
       fB->FixParameter(1,b);
       fB->FixParameter(2,a);
       float lim=zL[1];//z position of 2nd layer
       float aa=fB->Eval(lim);
       float diff=2*a*lim+2*a*zz0+b;
       float limo=zL[5];//z position of 6th layer
       float aaout=fB->Eval(limo);
       float diffout=2*a*limo+2*a*zz0+b;     
 
       //Deflection diff of slopes  
       float deflection=e->get_deflecPR();
       
       //Deflection2 diff of angles  
       float deflection2=TMath::ATan(diffout)-TMath::ATan(diff);
       if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_Eneg().at(0)<=0 )
         { 	
          for(int ijk=0;ijk<nnhits;ijk++)
            {
             if(e->get_hits().at(ijk)->get_flagPR()==true) HDdefvsstripID[i][(int)e->get_hits().at(ijk)->get_L()]->Fill(e->get_hits().at(ijk)->get_chip(),deflection2);
            
            } 
         }  
       if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_Eneg().at(0)<=0 )
        {  
         HDef[i]->Fill(deflection);
         HDef2[i]->Fill(deflection2);
        }
       
   	  if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_Eneg().at(0)<=0 )
        {  
         HDefT2[i]->Fill(deflection2);
			
        }		  
		  
	 if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_Eneg().at(0)<=0&&e->get_chi2BPR()<=5&&e->get_chi2NBPR()<=5)
        {  
         HDefnoT2[i]->Fill(deflection2);
         HDefnoT2vsT4[i]->Fill(e->get_EneT4().at(0),deflection2);
        }	 

       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_Eneg().at(0)<=0) 
        {
         hChi2B[i]->Fill(e->get_chi2BPR());   
         hChi2NB[i]->Fill(e->get_chi2NBPR());  
        }   
      } //j
   }//i

  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
	
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 gStyle->SetOptStat("emruo");
 //gStyle->SetOptStat(1);
 //gStyle->SetPalette( kDarkBodyRadiator);
 gStyle->SetPalette( kTemperatureMap);
 TGaxis::SetMaxDigits(4);

//////////////////////////////    
 

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.32);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      }    
 
    if (i==Nf-1)Can[i]->Print(Form("Deflectionbis%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }


 
//////////////////////////////    
//Deflection
//////////////////////////////   
 //gStyle->SetOptStat(0);
 gStyle->SetOptFit(1);

 TF1**f=new TF1*[Nf];	
	
 TCanvas**CanDef=new TCanvas*[Nf]; 
 for(int i=0;i<Nf;i++)
   {
    CanDef[i]=new TCanvas(Form("Deflection, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Deflection, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
   // CanDef[i]->Divide(2,1);
    CanDef[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    HDefT2[i]->GetXaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
    HDefT2[i]->GetYaxis()->SetTitle("entries");
    HDefT2[i]->SetTitle("Coincidence T1.T2.T3.noG");
    HDefT2[i]->GetYaxis()->SetTitleOffset(1.9);
    //HDefT2[i]->GetXaxis()->SetRangeUser(-0.05,0.05);
    HDefT2[i]->Draw("hist");       
    f[i]= new TF1(Form("f%d",i),"gaus",-0.2,0.2);
    f[i]->SetLineColor(kBlack);
    HDefT2[i]->Fit(f[i],"R","",-0.02,0.02);  
    f[i]->DrawF1(-0.02,0.02,"same");
    // HDefnoT2[i]->Draw("hist same");       
  //  CanDef[i]->cd(2);
   // gPad->SetBottomMargin(0.1);
    //gPad->SetLeftMargin(0.15);
    //HDefnoT2[i]->Draw("hist");       
    //HDefnoT2[i]->GetXaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
    //HDefnoT2[i]->GetYaxis()->SetTitle("entries");
    //HDefnoT2[i]->GetYaxis()->SetTitleOffset(1.9);
    //HDefnoT2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
    //HDefnoT2[i]->SetTitle("Coincidence T1.noT2.T3.T4.noG");
    CanDef[i]->Print(Form("Deflectionbis%s_%s.pdf",file[i].c_str(),coinc[i].c_str()),Form("Title:Deflection"));    
   } //  i  
 
 TCanvas**canproton=new TCanvas*[Nf];  
 for(int i=0;i<Nf;i++)
   {
	canproton[i]=new TCanvas(Form("Deflection, Run %s, Coincidence %s  proton",file[i].c_str(),coinc[i].c_str()),Form("Deflection, Run %s, Coincidence %s  proton",file[i].c_str(),coinc[i].c_str()),200,10,1000,1000);  
    canproton[i]->cd();
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.15);  
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);  
	HDefnoT2vsT4[i]->GetXaxis()->SetTitle("Pulse height in T4");
	HDefnoT2vsT4[i]->GetYaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
	HDefnoT2vsT4[i]->GetYaxis()->SetTitleOffset(1.9);
	HDefnoT2vsT4[i]->GetYaxis()->SetRangeUser(-0.2,0.2);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
	HDefnoT2vsT4[i]->Draw("colz");
	   
   }
	

//////////////////////////////    
//Chi2 distribution
////////////////////////////// 

 TCanvas**CanChi2=new TCanvas*[Nf]; 

 TF1* chi2k1=new TF1("chi2k1","TMath::Power(x,-1./2.)*TMath::Exp(-x/2.)/(TMath::Power(2,1./2.)*TMath::Gamma(1./2.))",0,20); 
 for(int i=0;i<Nf;i++)
   {
    CanChi2[i]=new TCanvas(Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanChi2[i]->Divide(2,1);
    CanChi2[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    hChi2NB[i]->GetXaxis()->SetTitle("Chi2 in Non-Bending plane");
    //hChi2NB[i]->GetYaxis()->SetTitle("arbitrary unit");
    hChi2NB[i]->DrawNormalized("hist");
    chi2k1->Draw("same");
    CanChi2[i]->cd(2);
    gPad->SetBottomMargin(0.1);
    hChi2B[i]->DrawNormalized("hist");       
    hChi2B[i]->GetXaxis()->SetTitle("Chi2 in Bending plane");
    chi2k1->Draw("same");

    //hChi2B[i]->GetYaxis()->SetTitle("arbitrary unit");
    if (i==Nf-1)CanChi2[i]->Print(Form("Deflectionbis%s_%s.pdf)",file[i].c_str(),coinc[i].c_str()),Form("Title:Chi2 distribution"));    
   } //  i   
   
   
//////////////////////////////    
//StripID vs deflection distribution
//////////////////////////////  

 TCanvas***Canstrip=new TCanvas**[Nf]; 

 for(int i=0;i<Nf;i++)
   {
    Canstrip[i]=new TCanvas*[7];
    for(int j=0;j<7;j++)
      {
       Canstrip[i][j]=new TCanvas(Form("Deflection vs Chip, L%d, Run %s, Coincidence %s ",j,file[i].c_str(),coinc[i].c_str()),Form("Deflection vs Chips, L%d, Run %s, Coincidence %s ",j,file[i].c_str(),coinc[i].c_str()),200,10,1600,800);   
       Canstrip[i][j]->cd();
       HDdefvsstripID[i][j]->GetYaxis()->SetRangeUser(-0.2,0.2);
       HDdefvsstripID[i][j]->GetXaxis()->SetRangeUser(0,12);
       HDdefvsstripID[i][j]->GetXaxis()->SetTitle("Chip (0 to 11)");
       HDdefvsstripID[i][j]->GetYaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
        gPad->SetGridx(1);
        gPad->SetGridy(1);      
       HDdefvsstripID[i][j]->Draw("colz");
      }  
   }
   
}//end function

void PalestineChip(string s,string c,int geoconf)
{
 
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

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


 int Nf=1;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;
 
 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=17;
 //string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 5+ layers & track","T1T3T4 & 6+ layers","T1T3T4 & 6+ layers & track","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track","T1T3T4 & L0L4L6","T1T3T4 & L1L3L5","T1T3T4 & (L0L4L6&L1L3L5)","T1T3T4 & (L0L4L6 | L1L3L5)"};
 
  string T[nx]={"All","T1T4","T1T4 & hits","T1T4 & 1+ layers","T1T4 & 2+ layers","T1T4 & 3+ layers","T1T4 & 4+ layers","T1T4 & 5+ layers","T1T4 & 5+ layers & track","T1T4 & 6+ layers","T1T4 & 6+ layers & track","T1T4 & 7 layers","T1T4 & 7 layers & track","T1T4 & L0L4L6","T1T4 & L1L3L5","T1T4 & (L0L4L6&L1L3L5)","T1T4 & (L0L4L6 | L1L3L5)"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];
 TH1F** hChi2B=new TH1F*[Nf];
 TH1F** hChi2NB=new TH1F*[Nf];

 //Deflection histograms
 TH1F**HDef=new TH1F*[Nf];
 TH1F**HDef2=new TH1F*[Nf];

 TH1F**HDefT2=new TH1F*[Nf];
 TH1F**HDefnoT2=new TH1F*[Nf];
 
 TH2F**HDefnoT2vsT4=new TH2F*[Nf];


 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
     //Chi2 histograms
   
    hChi2NB[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2NB",file[i].c_str(),coinc[i].c_str()),20,0,20);
    hChi2NB[i]->SetStats(1);
    
    hChi2B[i]= new TH1F(Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Chi2B",file[i].c_str(),coinc[i].c_str()),20,0,20);
    hChi2B[i]->SetStats(1);
   
    //Deflection histrograms
    HDef[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection",file[i].c_str(),coinc[i].c_str()),1000,-1,1);
    HDef[i]->SetStats(1);

    HDef2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection2",file[i].c_str(),coinc[i].c_str()),900,TMath::Pi(),TMath::Pi());
    HDef2[i]->SetStats(1);
 
	   
	HDefT2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection with T2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection with T2",file[i].c_str(),coinc[i].c_str()),900,TMath::Pi(),TMath::Pi());
    HDefT2[i]->SetStats(1);  
	HDefnoT2[i]= new TH1F(Form("Run %s, Coincidence %s, Deflection without T2",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Deflection without T2",file[i].c_str(),coinc[i].c_str()),900,TMath::Pi(),TMath::Pi());
    HDefnoT2[i]->SetStats(1);  
	
	HDefnoT2vsT4[i]= new TH2F(Form("Coincidence %s, Deflection without T2 vs T4",coinc[i].c_str()),Form("Coincidence %s, Deflection without T2 vs T4",coinc[i].c_str()),2000,0,2000,900,TMath::Pi(),TMath::Pi());
    HDefnoT2vsT4[i]->SetStats(1);  
	      
	   
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("../Data/%s.EVENT.root",file[i].c_str()));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0086.BPD.EVENT.root"));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0088.BPD.EVENT.root"));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0089.BPD.EVENT.root"));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int* Lay=new int[7];
       for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);

       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
       
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
 /*
       //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
        //T1&T3&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T3&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T3&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T3&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
      */
      
       //T1&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
       //T1&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
     
      
      
       int**Ltmp=new int*[7];
       for(int k=0;k<7;k++) 
         {
          Ltmp[k]=new int[nx];
          for(int kl=0;kl<nx;kl++) Ltmp[k][kl]=0;
         }

       //Deflection
       ///////////////////////////////    
       //Extract PR track information
       ///////////////////////////////    
       //Non bending plane
       float p0=e->get_interPR();
       float p1=e->get_slopePR();
       TF1* fNB=new TF1("fNB","pol1",-20,20);
       if(p1==0) continue;
       fNB->FixParameter(0,-p0/p1);
       fNB->FixParameter(1,1./p1);

       //Bending plane //ax^2+bx+c
       float zz0=0;
       float a=e->get_aPR();
       float b=e->get_bPR();
       float c=e->get_cPR();
       TF1* fB= new TF1("fB","[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,40);
       fB->FixParameter(3,zz0);
       fB->FixParameter(0,c);
       fB->FixParameter(1,b);
       fB->FixParameter(2,a);
       float lim=zL[1];//z position of 2nd layer
       float aa=fB->Eval(lim);
       float diff=2*a*lim+2*a*zz0+b;
       float limo=zL[5];//z position of 6th layer
       float aaout=fB->Eval(limo);
       float diffout=2*a*limo+2*a*zz0+b;     
 
       //Deflection diff of slopes  
       float deflection=e->get_deflecPR();
       
       //Deflection2 diff of angles  
       float deflection2=TMath::ATan(diffout)-TMath::ATan(diff);

	   if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0 &&nnhits==7&&Ti==127&& e->get_chi2BPR()<=2&& e->get_chi2NBPR()<=2)
        {  
         HDef[i]->Fill(deflection);
         HDef2[i]->Fill(deflection2);
			
        }
       
   	  if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0 &&nnhits==7&&Ti==127&& e->get_chi2BPR()<=10)
        {  
         HDefT2[i]->Fill(deflection2);
			
        }		  
		  
	   if(e->get_EneT1().at(0)>0&&e->get_EneT2().at(0)<=0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&e->get_Eneg().at(0)<=0 &&nnhits==7&&Ti==127&& e->get_chi2BPR()<=10)
        {  
         HDefnoT2[i]->Fill(deflection2);
		 HDefnoT2vsT4[i]->Fill(e->get_EneT4().at(0),deflection2);
        }			  

       if(e->get_EneT1().at(0)>0 &&e->get_EneT2().at(0)>0&&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0&&NL==7&&e->get_chi2NBPR()>0&&e->get_chi2BPR()>0) 
        {
         hChi2B[i]->Fill(e->get_chi2BPR());   
         hChi2NB[i]->Fill(e->get_chi2NBPR());  
        }   
      } //j
   }//i

  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
	
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 //gStyle->SetOptStat("emruo");
 gStyle->SetOptStat(0);
 gStyle->SetPalette( kDarkBodyRadiator);
 TGaxis::SetMaxDigits(4);

//////////////////////////////    
 

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.32);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      }    
 
    if (i==Nf-1)Can[i]->Print(Form("Deflectionbis%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }


 
//////////////////////////////    
//Deflection
//////////////////////////////   
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(1);

 TF1**f=new TF1*[Nf];	
	
 TCanvas**CanDef=new TCanvas*[Nf]; 
 for(int i=0;i<Nf;i++)
   {
    CanDef[i]=new TCanvas(Form("Deflection, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Deflection, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanDef[i]->Divide(2,1);
    CanDef[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    HDefT2[i]->GetXaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
    HDefT2[i]->GetYaxis()->SetTitle("entries");
    HDefT2[i]->SetTitle("Coincidence T1.T2.T3.T4.noG");
    HDefT2[i]->GetYaxis()->SetTitleOffset(1.9);
    HDefT2[i]->GetXaxis()->SetRangeUser(-0.05,0.05);
    HDefT2[i]->Draw("hist");       
    f[i]= new TF1(Form("f%d",i),"gaus",-0.2,0.2);
	f[i]->SetLineColor(kBlack);
	HDefT2[i]->Fit(f[i],"R","",-0.02,0.02);  
	f[i]->DrawF1(-0.02,0.02,"same");
	// HDefnoT2[i]->Draw("hist same");       
    CanDef[i]->cd(2);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    HDefnoT2[i]->Draw("hist");       
    HDefnoT2[i]->GetXaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
    HDefnoT2[i]->GetYaxis()->SetTitle("entries");
    HDefnoT2[i]->GetYaxis()->SetTitleOffset(1.9);
    HDefnoT2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
    HDefnoT2[i]->SetTitle("Coincidence T1.noT2.T3.T4.noG");
    CanDef[i]->Print(Form("Deflectionbis%s_%s.pdf",file[i].c_str(),coinc[i].c_str()),Form("Title:Deflection"));    
   } //  i  
 
 TCanvas**canproton=new TCanvas*[Nf];  
 for(int i=0;i<Nf;i++)
   {
	canproton[i]=new TCanvas(Form("Deflection, Run %s, Coincidence %s  proton",file[i].c_str(),coinc[i].c_str()),Form("Deflection, Run %s, Coincidence %s  proton",file[i].c_str(),coinc[i].c_str()),200,10,1000,1000);  
    canproton[i]->cd();
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.15);  
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);  
	HDefnoT2vsT4[i]->GetXaxis()->SetTitle("Pulse height in T4");
	HDefnoT2vsT4[i]->GetYaxis()->SetTitle("#theta_{Bout} #minus #theta_{Bin} (rad)");
	HDefnoT2vsT4[i]->GetYaxis()->SetTitleOffset(1.9);
	HDefnoT2vsT4[i]->GetYaxis()->SetRangeUser(-0.2,0.2);
	HDefnoT2vsT4[i]->Draw("colz");
	   
   }
	
	
//////////////////////////////    
//Chi2 distribution
////////////////////////////// 

 TCanvas**CanChi2=new TCanvas*[Nf]; 

 TF1* chi2k1=new TF1("chi2k1","TMath::Power(x,-1./2.)*TMath::Exp(-x/2.)/(TMath::Power(2,1./2.)*TMath::Gamma(1./2.))",0,20); 
 for(int i=0;i<Nf;i++)
   {
    CanChi2[i]=new TCanvas(Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Chip Occupancy, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1600,800);  
    CanChi2[i]->Divide(2,1);
    CanChi2[i]->cd(1);
    gPad->SetBottomMargin(0.1);
    hChi2NB[i]->GetXaxis()->SetTitle("Chi2 in Non-Bending plane");
    //hChi2NB[i]->GetYaxis()->SetTitle("arbitrary unit");
    hChi2NB[i]->DrawNormalized("hist");
    chi2k1->Draw("same");
    CanChi2[i]->cd(2);
    gPad->SetBottomMargin(0.1);
    hChi2B[i]->DrawNormalized("hist");       
    hChi2B[i]->GetXaxis()->SetTitle("Chi2 in Bending plane");
    chi2k1->Draw("same");

    //hChi2B[i]->GetYaxis()->SetTitle("arbitrary unit");
    if (i==Nf-1)CanChi2[i]->Print(Form("Deflectionbis%s_%s.pdf)",file[i].c_str(),coinc[i].c_str()),Form("Title:Chi2 distribution"));    
   } //  i   
   
   
}//end function





void NoiseInvestigation(string s,string c,int geoconf)
{
 
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("Dataparameters%d.dat",geoconf); 

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

  
 int Nf=1;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;
 
 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 
 const Int_t nx=17;
 string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 5+ layers & track","T1T3T4 & 6+ layers","T1T3T4 & 6+ layers & track","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track","T1T3T4 & L0L4L6","T1T3T4 & L1L3L5","T1T3T4 & (L0L4L6&L1L3L5)","T1T3T4 & (L0L4L6 | L1L3L5)"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];
 TH1F**MissingL=new TH1F*[Nf];

 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.EVENT.root",file[i].c_str()) <<endl;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
	
	//Missing Layer histograms   
    MissingL[i]= new TH1F(Form("Run %s, Coincidence %s, Missing Layer",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s, Missing Layer",file[i].c_str(),coinc[i].c_str()),7,0,7);
    MissingL[i]->SetFillColor(38);
	   
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("../Data/%s.EVENT.root",file[i].c_str()));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0086.BPD.EVENT.root"));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0088.BPD.EVENT.root"));
   // chain[i]->Add(Form("/data/psmangeard/AESOPLite/Data/NL0089.BPD.EVENT.root"));
    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%10000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int* Lay=new int[7];
       for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);

       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
       
	   if(NL==6)
	    {
		 for(int ij=0;ij<7;ij++)
		   {
		    if(Lay[ij]==0) MissingL[i]->Fill(ij);
		   } 
		}  
		  
		  
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
       //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
        //T1&T3&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T3&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT4().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T3&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T3&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
      
       int**Ltmp=new int*[7];
       for(int k=0;k<7;k++) 
         {
          Ltmp[k]=new int[nx];
          for(int kl=0;kl<nx;kl++) Ltmp[k][kl]=0;
         }
      } //j
   }//i

  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
	
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 //gStyle->SetOptStat("emruo");
 gStyle->SetOptStat(0);
 gStyle->SetPalette( kDarkBodyRadiator);
 TGaxis::SetMaxDigits(4);

//////////////////////////////    
 

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.32);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      }    
 
   // if (i==Nf-1)Can[i]->Print(Form("NoiseInvestigation%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }
	
 TCanvas**Can2=new TCanvas*[Nf]; 

 for(int i=0;i<Nf;i++)
   {
 
    Can2[i]=new TCanvas(Form("Can, Run %s, Coincidence %s, Missing Layer ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s, Missing Layer ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can2[i]->cd();
    gPad->SetGridy(1);
    MissingL[i]->Draw("");
    MissingL[i]->GetYaxis()->SetTitle("Entries");
    MissingL[i]->GetXaxis()->SetTitle("Layer (0 to 6)");
    MissingL[i]->SetMinimum(0);
   // if (i==Nf-1)Can[i]->Print(Form("NoiseInvestigation%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }
		


 
 
   
   
}//end function



void Timing(string s,string c)
{
 
 int Nf=1;
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[Nf];
 
 string file[1]={""};
 string coinc[1]={""};

 file[0]=s;
 coinc[0]=c;
 
 int* Nevents= new int[Nf];
 int*Duration= new int[Nf];//in second
 
 

	
 const Int_t nx=17;
 //string T[nx]={"All","T1T3T4","T1T3T4 & hits","T1T3T4 & 1+ layers","T1T3T4 & 2+ layers","T1T3T4 & 3+ layers","T1T3T4 & 4+ layers","T1T3T4 & 5+ layers","T1T3T4 & 5+ layers & track","T1T3T4 & 6+ layers","T1T3T4 & 6+ layers & track","T1T3T4 & 7 layers","T1T3T4 & 7 layers & track","T1T3T4 & L0L4L6","T1T3T4 & L1L3L5","T1T3T4 & (L0L4L6&L1L3L5)","T1T3T4 & (L0L4L6 | L1L3L5)"};
 string T[nx]={"All","T1T3T2","T1T3T2 & hits","T1T3T2 & 1+ layers","T1T3T2 & 2+ layers","T1T3T2 & 3+ layers","T1T3T2 & 4+ layers","T1T3T2 & 5+ layers","T1T3T2 & 5+ layers & track","T1T3T2 & 6+ layers","T1T3T2 & 6+ layers & track","T1T3T2 & 7 layers","T1T3T2 & 7 layers & track","T1T3T2 & L0L4L6","T1T3T2 & L1L3L5","T1T3T2 & (L0L4L6&L1L3L5)","T1T3T2 & (L0L4L6 | L1L3L5)"};
 
 //Rate histograms
 TH1F**TEff=new TH1F*[Nf];

 //Timing histograms
 
 
 TH2D**HtimeEVTvsPHA=new TH2D*[Nf];
 TH1D**DiffEVTvsPHA=new TH1D*[Nf];
 TH2D**DiffEVTvsPHAvsPHA=new TH2D*[Nf];
 TH2D**HtEVTvsPHA=new TH2D*[Nf];
 TH1D**DifftEVTvsPHA=new TH1D*[Nf];
 TH2D**DifftEVTvsPHAvsPHA=new TH2D*[Nf];
 TH2D**tEVTvsGoEVT=new TH2D*[Nf];
 TH2D**tPHAvsGoPHA=new TH2D*[Nf];
 TH1D**GlobaltEVTvsPHA=new TH1D*[Nf];
 TH1D**GlobalgoEVTvsPHA=new TH1D*[Nf];
 TH2D**GlobaltEVTvsEVTNum=new TH2D*[Nf];
 TH2D**GlobaltPHAvsEVTNum=new TH2D*[Nf];
 TH2D**GlobalDiffvsEVTNum=new TH2D*[Nf];
 TH2D**GlobalDiffvsQ1=new TH2D*[Nf];

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
 for(int i=0;i<Nf;i++)
   {
    cout << Form("../Data/%s.EVENT.root",file[i].c_str()) <<endl;
  
    chain[i]=new TChain("Data");
    chain[i]->Add(Form("../Data/%s.EVENT.root",file[i].c_str()));

    ALEvent *e = new ALEvent();      
    //Define variables to read event
    //Set address to access event data
    chain[i]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    Nevents[i]=chain[i]->GetEntries();
    cout << "Number  of events: " << Nevents[i] << endl;
    int timefirstevent=0; //in second from 2017
    chain[i]->GetEntry(0);
    int y=e->get_yPHA();
    int d=e->get_dPHA();
    int m=e->get_mPHA();
    int h=e->get_hPHA();
    int mi=e->get_miPHA();
    int s=e->get_sPHA();  
    timefirstevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
	
    double firsttPHA=0; 	   
    double firstgoPHA=0; 	   
    double firsttEVT=0; 	   
    double firstgoEVT=0; 	   
	 
	for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);   
	   if(e->get_GoEVT()<0)continue;
	   firsttPHA=e->get_tPHA(); 	   
       firstgoPHA=e->get_GoPHA(); 	   
       firsttEVT=e->get_tEVT(); 	   
       firstgoEVT=e->get_GoEVT(); 	   
	   j=Nevents[i];  
	  }   
	cout<< "First time PHA "  <<    firsttPHA <<endl;
	cout<< "First go PHA "  <<    firstgoPHA <<endl;
	cout<< "First time EVT "  <<    firsttEVT <<endl;
	cout<< "First go EVT "  <<    firstgoEVT <<endl;
	   
    int timelastevent=0; //in second from 2017
    chain[i]->GetEntry(Nevents[i]-1);
    y=e->get_yPHA();
    d=e->get_dPHA();
    m=e->get_mPHA();
    h=e->get_hPHA();
    mi=e->get_miPHA();
    s=e->get_sPHA();  
    timelastevent=doy(y,m,d)*24*3600+3600*h+60*mi+s;
    Duration[i]=timelastevent-timefirstevent;
    
    //Rate histograms
    TEff[i]= new TH1F(Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),nx,0,nx);
    TEff[i]->SetStats(0);
    TEff[i]->SetFillColor(38);
    TEff[i]->GetXaxis()->SetAlphanumeric();
    for(int j=1;j<=nx;j++)TEff[i]->GetXaxis()->SetBinLabel(j,T[j-1].c_str());
    
    //timing histograms
   	
	HtimeEVTvsPHA[i]= new TH2D(Form("Coincidence %s, Time Go EVT vs PHA",coinc[i].c_str()),Form("Coincidence %s, Time Go EVT vs PHA",coinc[i].c_str()),1000,firstgoPHA,1.5*Nevents[i]+firstgoPHA,1000,0,66000.);
    HtimeEVTvsPHA[i]->SetStats(1);  
	HtEVTvsPHA[i]= new TH2D(Form("Coincidence %s, Time EVT vs PHA",coinc[i].c_str()),Form("Coincidence %s, Time EVT vs PHA",coinc[i].c_str()),1000,0,4500000000,1000,0,4500000000);
    HtimeEVTvsPHA[i]->SetStats(1);  
	DiffEVTvsPHA[i]= new TH1D(Form("Coincidence %s, Time  Go EVT-PHA",coinc[i].c_str()),Form("Coincidence %s, Time Go EVT-PHA",coinc[i].c_str()),100,-500,50000);
	DiffEVTvsPHAvsPHA[i]= new TH2D(Form("Coincidence %s, Time Go EVT-PHA vs PHA",coinc[i].c_str()),Form("Coincidence %s, Time Go EVT-PHA vs PHA",coinc[i].c_str()),1000,firstgoPHA,1.5*Nevents[i]+firstgoPHA,100,-50,50);
    DiffEVTvsPHAvsPHA[i]->SetStats(1);
	   
	DifftEVTvsPHA[i]= new TH1D(Form("Coincidence %s, Time   EVT-PHA",coinc[i].c_str()),Form("Coincidence %s, Time EVT-PHA",coinc[i].c_str()),2000,-60000,+30000);
	DifftEVTvsPHAvsPHA[i]= new TH2D(Form("Coincidence %s, Time  EVT-PHA vs PHA",coinc[i].c_str()),Form("Coincidence %s, Time EVT-PHA vs PHA",coinc[i].c_str()),1000,0,4500000000,2000,-60000,30000);
    DifftEVTvsPHAvsPHA[i]->SetStats(1); 

	tEVTvsGoEVT[i]= new TH2D(Form("Coincidence %s, Time  EVT vs Go EVT",coinc[i].c_str()),Form("Coincidence %s, Time  EVT vs Go EVT",coinc[i].c_str()),1000,0,66000,1000,0,4500000000);
    tEVTvsGoEVT[i]->SetStats(1); 
	tPHAvsGoPHA[i]= new TH2D(Form("Coincidence %s, Time  PHA vs Go PHA",coinc[i].c_str()),Form("Coincidence %s, Time  PHA vs Go PHA",coinc[i].c_str()),1000,firstgoPHA,1.5*Nevents[i]+firstgoPHA,1000,0,4500000000);
    tPHAvsGoPHA[i]->SetStats(1); 
		     
	GlobaltEVTvsPHA[i]= new TH1D(Form("Coincidence %s, Global Time   EVT-PHA",coinc[i].c_str()),Form("Coincidence %s, Global Time EVT-PHA",coinc[i].c_str()),20000,-100000,+3000000);
	GlobalgoEVTvsPHA[i]= new TH1D(Form("Coincidence %s, Global Go   EVT-PHA",coinc[i].c_str()),Form("Coincidence %s, Global Go EVT-PHA",coinc[i].c_str()),400,-1000,100);
   
	GlobaltEVTvsEVTNum[i]= new TH2D(Form("Coincidence %s, Time  EVT vs Event number",coinc[i].c_str()),Form("Coincidence %s, Time  EVT vs Event number",coinc[i].c_str()),10000,0,Nevents[i],1000,0,450000000000);
    GlobaltEVTvsEVTNum[i]->SetStats(1); 
	GlobaltPHAvsEVTNum[i]= new TH2D(Form("Coincidence %s, Time  PHA vs Event number",coinc[i].c_str()),Form("Coincidence %s, Time  PHA vs Event number",coinc[i].c_str()),10000,0,Nevents[i],1000,0,450000000000);
    GlobaltPHAvsEVTNum[i]->SetStats(1); 
	GlobalDiffvsEVTNum[i]= new TH2D(Form("Coincidence %s, Difference  EVT-PHA vs Event number",coinc[i].c_str()),Form("Coincidence %s, Difference  EVT-PHA vs Event number",coinc[i].c_str()),10000,0,Nevents[i],1000,-500,500);
    GlobalDiffvsEVTNum[i]->SetStats(1);      
	GlobalDiffvsQ1[i]= new TH2D(Form("Coincidence %s, Difference  EVT-PHA vs Q1",coinc[i].c_str()),Form("Coincidence %s, Difference  EVT-PHA vs Q1",coinc[i].c_str()),10,0,10,20000,-300000,+3000000);

	   
	   
	   
    //Loop over events
    for(int j=0;j<Nevents[i];j++)
      {
       chain[i]->GetEntry(j);
       if(j%50000==0) cout << "Event: " << j <<endl;

       int nnhits=e->get_Nhits();
       uint8_t Ti=(uint8_t)e->get_Ti();
       //Number of layers wih hit(s)
       int NL=0;
       for(int ij=0;ij<7;ij++) NL+=(int)((Ti >>ij) & 0x01);
       
       if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
       int* Lay=new int[7];
       for(int ij=0;ij<7;ij++) Lay[ij]=(int)((Ti >>ij) & 0x01);

       int*Ncase=new int[nx];
       for(int k=0;k<nx;k++)Ncase[k]=0;
 
	   	   
	   
 	/*   if((e->get_eventnumber()>=0&&e->get_eventnumber()<=200 )||(e->get_eventnumber()>=562&&e->get_eventnumber()<=598))	
 	//   if((e->get_eventnumber()>=72930&&e->get_eventnumber()<=72990 ))	
		{
		 cout << "Event Number: "<< e->get_eventnumber() << " ";//<<endl;
		 cout << e->get_yPHA() << "/";	
		 cout << e->get_mPHA() << "/";	
		 cout << e->get_dPHA() << " ";	
		 cout << e->get_hPHA() << ":";	
		 cout << e->get_miPHA() << ":";	
		 cout << e->get_sPHA() << " ";	
		// cout << "time PHA: " <<  e->get_tPHA()	 << endl; 
		// cout << "go PHA: " <<  e->get_GoPHA()	 << endl; 
		// cout << "time EVT: " <<  e->get_tEVT()	 << endl; 
		// cout << "go EVT: " <<  e->get_GoEVT()	 << endl; 
  		 cout <<  Form("%12.0lf",(double)e->get_GoPHA())<< " " ;
  		 cout <<  Form("%12.0lf",(double)e->get_tPHA()) << " " ;
  		 cout <<  Form("%12.0lf",(double)e->get_GoEVT() )<< " " ;
  		 cout <<  Form("%12.0lf",(double)e->get_tEVT()) << " " ;
		 cout <<endl;
		}	  	  
	*/	  
       ////////////////////////////////////
       //Timing
       ////////////////////////////////////  
 	   if(e->get_GoEVT()>-1&&e->get_GoPHA()>-1&&e->get_tEVT()>-1&&e->get_tPHA()>-1)
	    {			
		 HtimeEVTvsPHA[i]->Fill(e->get_GoPHA(),e->get_GoEVT());
		 DiffEVTvsPHA[i]->Fill(e->get_GoEVT()-e->get_GoPHA());
		 DiffEVTvsPHAvsPHA[i]->Fill(e->get_GoPHA(),e->get_GoEVT()-e->get_GoPHA());
		 
		 HtEVTvsPHA[i]->Fill(e->get_tPHA(),e->get_tEVT());
		 DifftEVTvsPHA[i]->Fill(e->get_tEVT()-e->get_tPHA());
		 DifftEVTvsPHAvsPHA[i]->Fill(e->get_tPHA(),e->get_tEVT()-e->get_tPHA());	
						
		 tEVTvsGoEVT[i]->Fill(e->get_GoEVT(),e->get_tEVT());
		 tPHAvsGoPHA[i]->Fill(e->get_GoPHA(),e->get_tPHA());

		 //Global timing	
		 if(e->get_GoEVT()<prevgoEVT&& e->get_GoEVT()<0.01*Overflowgo&& prevgoEVT>0.9 *Overflowgo)//overflow go
		  {
		   OverflowgoTagEVT++;
		  }  
		 prevgoEVT=	e->get_GoEVT();
		 goEVT	=e->get_GoEVT()+OverflowgoTagEVT*Overflowgo;
		 
		 if(e->get_GoPHA()<prevgoPHA&& e->get_GoPHA()<0.01*Overflowtime&& prevgoPHA>0.9 *Overflowtime)//overflow go PHA same as time 2^32
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
		
 	    if(e->get_eventnumber()>=0&&e->get_eventnumber()<=40 )	
// 	   if((e->get_eventnumber()>=72930&&e->get_eventnumber()<=72990 ))	
	  	  {
		   cout << "Diff Evt Num: "<< e->get_eventnumber() << " ";//<<endl;
		   // cout << "time PHA: " <<  e->get_tPHA()	 << endl; 
		   // cout << "go PHA: " <<  e->get_GoPHA()	 << endl; 
		   // cout << "time EVT: " <<  e->get_tEVT()	 << endl; 
		  // cout << "go EVT: " <<  e->get_GoEVT()	 << endl; 
  		    cout << e->get_yPHA() << "/";	
		   cout << e->get_mPHA() << "/";	
		   cout << e->get_dPHA() << " ";	
		   cout << e->get_hPHA() << ":";	
		   cout << e->get_miPHA() << ":";	
		   cout << e->get_sPHA() << " ";
		   cout <<  Form("%12.0lf",(double)e->get_GoPHA())<< " " ;
  		   cout <<  Form("%12.0lf",(double)e->get_tPHA()) << " " ;
  		   cout <<  Form("%12.0lf",(double)e->get_GoEVT() )<< " " ;
  		   cout <<  Form("%12.0lf",(double)e->get_tEVT()) << " " ;  
		   cout << Form("%12.0lf",(double)(goPHA-firstgoPHA)) << " " ;
  		   cout << Form("%12.0lf",(double)(goEVT-firstgoEVT)) << " " ;
  		   cout << Form("%12.0lf",(double)(timePHA-firsttPHA))<< " " ; 	
  		   cout << Form("%12.0lf",(double)(timeEVT-firsttEVT))<< " " ; 	
  		   cout << Form("%12.0lf",(double)(goEVT-firstgoEVT-goPHA+firstgoPHA)) << " " ;
  		   cout << Form("%12.0lf",(double)(timeEVT-firsttEVT-timePHA+firsttPHA))<< " " ; 	
		   cout <<endl;
	  	  }
		
			
	     GlobalgoEVTvsPHA[i]->Fill(goEVT-firstgoEVT-goPHA+firstgoPHA);
	     GlobaltEVTvsPHA[i]->Fill(timeEVT-firsttEVT-timePHA+firsttPHA);
		 GlobalDiffvsQ1[i]->Fill(e->get_Q1EVT(),timeEVT-firsttEVT-timePHA+firsttPHA);
		 GlobaltEVTvsEVTNum[i]->Fill(e->get_eventnumber(),timeEVT);
		 GlobaltPHAvsEVTNum[i]->Fill(e->get_eventnumber(),timePHA);
	     GlobalDiffvsEVTNum[i]->Fill(e->get_eventnumber(),timeEVT-firsttEVT-timePHA+firsttPHA);
			
		} 
		  
       ////////////////////////////////////
       //TRIGGER CUT
       ////////////////////////////////////  
       //No Cut
       TEff[i]->Fill(T[0].c_str(),1);  Ncase[0]=1;   
	  // if(fabs(timeEVT-firsttEVT-timePHA+firsttPHA+26530)>10)continue;  
	   //T1&T3&T4
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0){TEff[i]->Fill(T[1].c_str(),1);Ncase[1]=1;}
       //T1&T3&T4&hits
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0){TEff[i]->Fill(T[2].c_str(),1);Ncase[2]=1;}
       if(e->get_Q1EVT()!=5)continue;
       if(fabs((int)(timeEVT-firsttEVT-timePHA+firsttPHA)%65536)>5)continue;
	   //T1&T3&T4&1+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=1){TEff[i]->Fill(T[3].c_str(),1);Ncase[3]=1;}
       //T1&T3&T4&2+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=2){TEff[i]->Fill(T[4].c_str(),1);Ncase[4]=1;}
       //T1&T3&T4&3+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=3){TEff[i]->Fill(T[5].c_str(),1);Ncase[5]=1;}
       //T1&T3&T4&4+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=4){TEff[i]->Fill(T[6].c_str(),1);Ncase[6]=1;}
       //T1&T3&T4&5+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=5){TEff[i]->Fill(T[7].c_str(),1);Ncase[7]=1;}
       //T1&T3&T4&5+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=5&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[8].c_str(),1);Ncase[8]=1;}       
       //T1&T3&T4&6+ layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=6){TEff[i]->Fill(T[9].c_str(),1);Ncase[9]=1;}
       //T1&T3&T4&6+ layers &track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL>=6&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[10].c_str(),1);Ncase[10]=1;}       
       //T1&T3&T4&7 layers
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL==7){TEff[i]->Fill(T[11].c_str(),1);Ncase[11]=1;}
       //T1&T3&T4&7 layers & track
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && NL==7&&e->get_chi2NBPR()>=0&&e->get_chi2BPR()>=0){TEff[i]->Fill(T[12].c_str(),1);Ncase[12]=1;}
        //T1&T3&T4& L0L4L6
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && Lay[0]==1&&Lay[4]==1&&Lay[6]==1){TEff[i]->Fill(T[13].c_str(),1);Ncase[13]=1;}
       //T1&T3&T4& L1L3L5
       if(e->get_EneT1().at(0)>0 &&e->get_EneT3().at(0)>0&&e->get_EneT2().at(0)>0&&nnhits>0 && Lay[1]==1&&Lay[3]==1&&Lay[5]==1){TEff[i]->Fill(T[14].c_str(),1);Ncase[14]=1;}
       //T1&T3&T4& L0L1L3||L4L5L6
       if(Ncase[13]==1&&Ncase[14]==1){TEff[i]->Fill(T[15].c_str(),1);Ncase[15]=1;}
       //T1&T3&T4& L0L1L3L4L5L6
       if(Ncase[13]==1||Ncase[14]==1){TEff[i]->Fill(T[16].c_str(),1);Ncase[16]=1;}
      
       int**Ltmp=new int*[7];
       for(int k=0;k<7;k++) 
         {
          Ltmp[k]=new int[nx];
          for(int kl=0;kl<nx;kl++) Ltmp[k][kl]=0;
         } 
      } //j
   }//i

  for(int i=0;i<Nf;i++)
   {
    cout << "File: " << Form("../Data/%s.BPD.EVENT.root",file[i].c_str()) <<endl;
    cout << "Number  of events: " << Nevents[i] << endl;
    cout << "Duration: " << Duration[i] <<endl;
    cout << "Trigger coincidence rate: " <<  (float) Nevents[i]/(float) Duration[i] << "Hz" << endl;
    for(int j=0;j<nx;j++)
      {
       cout << T[j].c_str() << ": " << TEff[i]->GetBinContent(j+1)<< endl;
      }
   } //i
	
//////////////////////////////////   
//   DISPLAY
//////////////////////////////////   
 //gStyle->SetOptStat("emruo");
 gStyle->SetOptStat(1);
 gStyle->SetPalette( kDarkBodyRadiator);
 //TGaxis::SetMaxDigits(4);

//////////////////////////////    
 

 TCanvas**Can=new TCanvas*[Nf]; 
 TText*** txt=new TText**[Nf];

 for(int i=0;i<Nf;i++)
   {
 
    Can[i]=new TCanvas(Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Can, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,1300,800);  
    Can[i]->cd();
    Can[i]->SetBottomMargin(0.32);
    gPad->SetGridy(1);
    if(Duration[i]>0)TEff[i]->SetTitle(Form("Run %s, Coincidence %s: %3.3f Hz",file[i].c_str(),coinc[i].c_str(),(float) Nevents[i]/(float) Duration[i]));
    TEff[i]->SetLineWidth(0);
    TEff[i]->SetMarkerStyle(20);
    TEff[i]->SetMarkerSize(0);
    TEff[i]->Scale(100./Nevents[i]);
    TEff[i]->SetBarWidth(0.8);
    TEff[i]->SetBarOffset(0.1);
    TEff[i]->GetYaxis()->SetRangeUser(0,110);
    TEff[i]->GetYaxis()->SetNdivisions(511,0);
    TEff[i]->GetYaxis()->SetTitle("Normalized rates");
    TEff[i]->LabelsOption("v","X");
    TEff[i]->Draw("bar");
 
    txt[i]=new TText*[nx]; 
    for(int j=0;j<nx;j++)
      {
       if(TEff[i]->GetBinContent(j+1)>=100&&j==0) txt[i][j]=new TText(j+0.2,TEff[i]->GetBinContent(j+1)+1,Form("%d",(int)TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=100) txt[i][j]=new TText(j+0.1,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else if(TEff[i]->GetBinContent(j+1)>=10) txt[i][j]=new TText(j+0.15,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       else txt[i][j]=new TText(j+0.25,TEff[i]->GetBinContent(j+1)+1,Form("%3.2f",TEff[i]->GetBinContent(j+1))); 
       txt[i][j]->SetTextSize(0.02);    
       txt[i][j]->Draw();    
      }    
 
    //if (i==Nf-1)Can[i]->Print(Form("timing%s_%s.pdf(",file[i].c_str(),coinc[i].c_str()),Form("Title:Rates"));    
   }

//////////////////////////////    

 TCanvas**CanT=new TCanvas*[Nf]; 
 for(int i=0;i<Nf;i++)
   {
    CanT[i]=new TCanvas(Form("Time, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Time, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,2100,1500);  
    CanT[i]->Divide(4,2);

	CanT[i]->cd(1);
    gPad->SetLeftMargin(0.15);
    tEVTvsGoEVT[i]->GetXaxis()->SetTitle("EVT Go");
    tEVTvsGoEVT[i]->GetYaxis()->SetTitleOffset(1.9);
    tEVTvsGoEVT[i]->GetYaxis()->SetTitle("EVT Time");
	tEVTvsGoEVT[i]->Draw("colz");     
	   
	CanT[i]->cd(2);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    HtimeEVTvsPHA[i]->GetXaxis()->SetTitle("PHA Go");
    HtimeEVTvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    HtimeEVTvsPHA[i]->GetYaxis()->SetTitle("EVT Go");
	HtimeEVTvsPHA[i]->Draw("colz");   
	   
    CanT[i]->cd(3);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    DiffEVTvsPHA[i]->GetYaxis()->SetTitle("Entries");
    DiffEVTvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    DiffEVTvsPHA[i]->GetXaxis()->SetTitle("EVT Go - PHA Go");
	DiffEVTvsPHA[i]->Draw("hist"); 	   

	CanT[i]->cd(4);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    DiffEVTvsPHAvsPHA[i]->GetXaxis()->SetTitle("PHA Go");
    DiffEVTvsPHAvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    DiffEVTvsPHAvsPHA[i]->GetYaxis()->SetTitle("EVT Go - PHA Go");
	DiffEVTvsPHAvsPHA[i]->Draw("colz");  

	CanT[i]->cd(5);
    gPad->SetLeftMargin(0.15);
    tPHAvsGoPHA[i]->GetXaxis()->SetTitle("PHA Go");
    tPHAvsGoPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    tPHAvsGoPHA[i]->GetYaxis()->SetTitle("PHA Time");
	tPHAvsGoPHA[i]->Draw("colz");     
	      
	CanT[i]->cd(6);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    HtEVTvsPHA[i]->GetXaxis()->SetTitle("PHA time");
    HtEVTvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    HtEVTvsPHA[i]->GetYaxis()->SetTitle("EVT time");
	HtEVTvsPHA[i]->Draw("colz");   

    CanT[i]->cd(7);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    gPad->SetLogy(1);	
    DifftEVTvsPHA[i]->GetYaxis()->SetTitle("Entries");
    DifftEVTvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    DifftEVTvsPHA[i]->GetXaxis()->SetTitle("EVT time - PHA time");
	DifftEVTvsPHA[i]->Draw("hist"); 	   

   	CanT[i]->cd(8);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);
	gPad->SetLogz(1);	
   
    DifftEVTvsPHAvsPHA[i]->GetXaxis()->SetTitle("PHA time");
    DifftEVTvsPHAvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    DifftEVTvsPHAvsPHA[i]->GetYaxis()->SetTitle("EVT time - PHA time");
	DifftEVTvsPHAvsPHA[i]->Draw("colz");  
   }   
//////////////////////////////    

 TCanvas**CanG=new TCanvas*[Nf]; 
 for(int i=0;i<Nf;i++)
   {
    CanG[i]=new TCanvas(Form("Global Time, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),Form("Global Time, Run %s, Coincidence %s ",file[i].c_str(),coinc[i].c_str()),200,10,2100,1500);  
    CanG[i]->Divide(3,2);

    CanG[i]->cd(1);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    GlobalgoEVTvsPHA[i]->GetYaxis()->SetTitle("Entries");
    GlobalgoEVTvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    GlobalgoEVTvsPHA[i]->GetXaxis()->SetTitle("EVT Go - PHA Go");
	GlobalgoEVTvsPHA[i]->Draw("hist"); 	   

    CanG[i]->cd(2);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    GlobaltEVTvsPHA[i]->GetYaxis()->SetTitle("Entries");
    GlobaltEVTvsPHA[i]->GetYaxis()->SetTitleOffset(1.9);
    GlobaltEVTvsPHA[i]->GetXaxis()->SetTitle("EVT Time - PHA Time");
	GlobaltEVTvsPHA[i]->Draw("hist"); 	   
	
    CanG[i]->cd(3);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    GlobaltEVTvsEVTNum[i]->GetYaxis()->SetTitle("EVT Global time");
    GlobaltEVTvsEVTNum[i]->GetYaxis()->SetTitleOffset(1.9);
    GlobaltEVTvsEVTNum[i]->GetXaxis()->SetTitle("Event number");
	GlobaltEVTvsEVTNum[i]->Draw("colz"); 	   
	
    CanG[i]->cd(4);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    GlobaltPHAvsEVTNum[i]->GetYaxis()->SetTitle("PHA Global time");
    GlobaltPHAvsEVTNum[i]->GetYaxis()->SetTitleOffset(1.9);
    GlobaltPHAvsEVTNum[i]->GetXaxis()->SetTitle("Event number");
    GlobaltPHAvsEVTNum[i]->Draw("colz"); 	   

    CanG[i]->cd(5);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    GlobalDiffvsQ1[i]->GetYaxis()->SetTitle("Global time difference EVT - PHA");
    GlobalDiffvsQ1[i]->GetYaxis()->SetTitleOffset(1.9);
    GlobalDiffvsQ1[i]->GetXaxis()->SetTitle("Q1");
    GlobalDiffvsQ1[i]->Draw("colz"); 	  
   
    CanG[i]->cd(6);
    gPad->SetLeftMargin(0.15);
    //gPad->SetGridy(1);	
    GlobalDiffvsEVTNum[i]->GetYaxis()->SetTitle("Global time difference EVT - PHA");
    GlobalDiffvsEVTNum[i]->GetYaxis()->SetTitleOffset(1.9);
    GlobalDiffvsEVTNum[i]->GetXaxis()->SetTitle("Event number");
    GlobalDiffvsEVTNum[i]->Draw("colz"); 	   

   }   


}//end function



void PalestinePulsePHADisc()
{

int NDisc=5;
int filename[9]={2007,2006,2008,2009,2010};
int Disc[9]={0,7,10,20,50};
    
//////////////////////////// 
//Trigger related histograms
//////////////////////////// 
 
 TH1F**AllPH=new TH1F*[NDisc];
 string XtitlePH="Pulse height";
 string YtitlePH="Entries";


 for(int i=0;i<NDisc;i++)
   {
    AllPH[i]=new TH1F(Form("T1 PHA Disc %d",Disc[i]),Form("T1 PHA Disc %d",Disc[i]),512,0,2048);
    AllPH[i]->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
    AllPH[i]->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
    AllPH[i]->SetLineColor(kBlack);
    AllPH[i]->SetLineWidth(1);
    AllPH[i]->GetXaxis()->SetRangeUser(0,500);
   } //i  
    
 //Create TChain to merge the files for a given energy
 TChain*chain;
 
 
 for(int i=0;i<NDisc;i++)
   {  
     
    chain=new TChain("Data");

    chain->Add(Form("../Data/NL%4d.BPD.EVENT.root",filename[i]));

    //Define variables to read event
    ALEvent *e = new ALEvent();      
    //Set address to access event data
    chain->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    int nentries=chain->GetEntries();
    cout << "Number  of events: " << nentries << endl;  
      
    for (int j=0;j<nentries;j++)
      {
       chain->GetEntry(j); //Load the entry i in the variable e 

       if(j%1000==0) cout << "Event: " << j <<endl;

       ////////////////////////////////////
       //TRIGGER HEIGHTS 
       ////////////////////////////////////  
    
       AllPH[i]->Fill(e->get_EneT1().at(0));
    
      }//j
   
   } 
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 
 TCanvas*CanPH=new TCanvas("Pulse heights","Pulse heights",200,10,1200,800);  
 CanPH->Divide(3,2);
 
 for(int j=0;j<NDisc;j++)
   {
    CanPH->cd(j+1);  
    AllPH[j]->GetXaxis()->SetRangeUser(1.,500);
    AllPH[j]->Draw();
   }
   
/* CanPH->Print(Form("Palestine18_DiscPHAT1.pdf("));    
 for(int j=0;j<NDisc;j++)
   {
    CanPH->cd(j+1);  
    gPad->SetLogy(0);
    gPad->SetLeftMargin(0.15);
    AllPH[j]->SetLineColor(kBlack);
    AllPH[j]->SetLineWidth(2);
    //AllPH[j]->GetXaxis()->SetRangeUser(1.,500);
    AllPH[j]->GetYaxis()->SetTitleOffset(1.9);
    gStyle->SetOptStat(0);
    AllPH[j]->Draw();
   }
 CanPH->Print(Form("Palestine18_DiscPHAT1.pdf"),Form("Title:Pulse heights"));    
  

 TCanvas**CanPHi=new TCanvas*[NDisc];

 for(int j=0;j<NDisc;j++)
   {
    CanPHi[j]=new TCanvas(Form("Pulse heights %s",sPH[j].c_str()),Form("Pulse heights %s",sPH[j].c_str()),200,10,1200,800);  
    CanPHi[j]->cd();
    gPad->SetLeftMargin(0.15);
    AllPH[j]->Draw();
    AllPH[j]->GetXaxis()->SetRangeUser(-1.,4096);
    CanPHi[j]->Print(Form("Palestine18_DiscPHA_T1.pdf"),Form("Title:PHA Disc %d",Disc[j]));    
    AllPH[j]->GetXaxis()->SetRangeUser(1.,500);
    CanPHi[j]->Print(Form("Palestine18_DiscPHA_T1.pdf"),Form("Title:PHA Disc %d Zoom",Disc[j])));    
    gPad->SetLogy(1);
    if(j!=NDisc-1)CanPHi[j]->Print(Form("Palestine18_DiscPHA_T1.pdf"),Form("Title:PHA Disc %s ZoomLog",Disc[j]));   
    else CanPHi[j]->Print(Form("Palestine18_DiscPHA_T1.pdf)"),Form("Title:PHA Disc %d ZoomLog",Disc[j]));   
    
   }	   
 
*/
 
 
 
}


void UnbiasPulses()
{

int N=5;
string Trig[5]={"T1","T2","T3","T4","G"};
string Coin[5]={"T3&T4","T1&T4","T1&T4","T1&T3","T1&T4"};
    



//////////////////////////// 
//Trigger related histograms
//////////////////////////// 
 
 TH1F**AllPH=new TH1F*[N];
 string XtitlePH="Pulse height";
 string YtitlePH="Entries";


 for(int i=0;i<N;i++)
   {
    if(i==0||i==4)AllPH[i]=new TH1F(Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),400,0,2048);
    else AllPH[i]=new TH1F(Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),512,0,2048);
    AllPH[i]->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
    AllPH[i]->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
    AllPH[i]->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
    AllPH[i]->SetLineColor(kBlack);
    AllPH[i]->SetLineWidth(1);
    AllPH[i]->GetXaxis()->SetRangeUser(0,500);
   } //i  
    
 //Create TChain to merge the files for a given energy
 TChain**chain=new TChain*[3];
 
 for(int i=0;i<3;i++)
   {
    chain [i]=new TChain("Data");  
    if(i==0 )chain[i]->Add(Form("../Data/NL2013.BPD.EVENT.root"));
    else if(i==1)
     {
      chain[i]->Add(Form("../Data/NL0212.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0214.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0218.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0221.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0250.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0252.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0253.BPD.EVENT.root"));
      chain[i]->Add(Form("../Data/NL0254.BPD.EVENT.root"));
    //  chain->Add(Form("../Data/NL0255.BPD.EVENT.root"));
     }
    else chain[i]->Add(Form("../Data/NL2014.BPD.EVENT.root"));
       
   }

 for(int i=0;i<N;i++)
   {  
    //Define variables to read event
    ALEvent *e = new ALEvent();    
    int kk=0;
    if(i==1||i==2||i==4) kk=1;
    if(i==3) kk=2;
    //Set address to access event data
    chain[kk]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    int nentries=chain[kk]->GetEntries();
    cout << "Number  of events: " << nentries << endl;  
      
    for (int j=0;j<nentries;j++)
      {
       chain[kk]->GetEntry(j); //Load the entry i in the variable e 

       if(j%1000==0) cout << "Event: " << j <<endl;

       ////////////////////////////////////
       //TRIGGER HEIGHTS 
       ////////////////////////////////////  
    
       if(i==0)AllPH[i]->Fill(e->get_EneT1().at(0));
       if(i==1)AllPH[i]->Fill(e->get_EneT2().at(0));
       if(i==2)AllPH[i]->Fill(e->get_EneT3().at(0));
       if(i==3)AllPH[i]->Fill(e->get_EneT4().at(0));
       if(i==4)AllPH[i]->Fill(e->get_Eneg().at(0));
    
      }//j
   
   } 
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 
 TCanvas*CanPH=new TCanvas("Pulse heights","Pulse heights",200,10,1200,800);  
 CanPH->Divide(3,2);
 
 for(int j=0;j<N;j++)
   {
    CanPH->cd(j+1);  
    gPad->SetLeftMargin(0.15);
    AllPH[j]->GetXaxis()->SetRangeUser(1.,500);
    AllPH[j]->Draw();
    int binmax = AllPH[j]->GetMaximumBin(); 
    double x = AllPH[j]->GetXaxis()->GetBinCenter(binmax);
    cout << Trig[j]<<", Maximum PHA (Bin center): "<< x << endl;    
    
   }
   

 
 
 
}



void PalestineLOGICDisc()
{

 int N=5;
 string Trig[5]={"T1","T2","T3","T4","G"};   
 int NDiscmax=0;   

 int NDisc[5]={8,9,7,10,9};
 for(int i=0;i<N;i++)
  {
    if(NDisc[i]>NDiscmax)  NDiscmax=NDisc[i];
  } 
 cout << "NDiscmax: " << NDiscmax   <<endl;
 int** filename=new int*[N];
 float** Disc=new float*[N];
 int** PHA=new int*[N];
  
 int color[10]={1,2,4,6,41,34,46,12,42,38};    
  cout << "Before set up filename etc.... "  <<endl;

 for(int i=0;i<N;i++)
   {
    filename[i]=new int[NDiscmax];
    Disc[i]=new float[NDiscmax];
    PHA[i]=new int[NDiscmax];
    for(int j=0;j<NDiscmax;j++)
      {
       filename[i][j]=0;PHA[i][j]=0;Disc[i][j]=0; 
       if(i==0 && j<NDisc[i])//T1
        {
         filename[i][0]=2021;PHA[i][0]=4;Disc[i][0]=4;
         filename[i][1]=2022;PHA[i][1]=4;Disc[i][1]=5;
         filename[i][2]=2023;PHA[i][2]=4;Disc[i][2]=6;
         filename[i][3]=2020;PHA[i][3]=4;Disc[i][3]=7;
         filename[i][4]=2006;PHA[i][4]=7;Disc[i][4]=7;
         filename[i][5]=2017;PHA[i][5]=7;Disc[i][5]=8;
         filename[i][6]=2019;PHA[i][6]=7;Disc[i][6]=9;
         filename[i][7]=2012;PHA[i][7]=7;Disc[i][7]=10;
        }                                      
      
      if(i==1 && j<NDisc[i])//T2
        {
         filename[i][0]=2048;PHA[i][0]=3;Disc[i][0]=3;
         filename[i][1]=2049;PHA[i][1]=3;Disc[i][1]=4;
         filename[i][2]=2050;PHA[i][2]=3;Disc[i][2]=5;
         filename[i][3]=2051;PHA[i][3]=3;Disc[i][3]=6;
         filename[i][4]=2052;PHA[i][4]=3;Disc[i][4]=7;
         filename[i][5]=2053;PHA[i][5]=7;Disc[i][5]=7;
         filename[i][6]=2057;PHA[i][6]=7;Disc[i][6]=8;
         filename[i][7]=2058;PHA[i][7]=7;Disc[i][7]=9;
         filename[i][8]=2059;PHA[i][8]=7;Disc[i][8]=10;
        }    
      if(i==2 && j<NDisc[i])//T3
        {
         filename[i][0]=2031;PHA[i][0]=5;Disc[i][0]=5;
         filename[i][1]=2032;PHA[i][1]=5;Disc[i][1]=6;
         filename[i][2]=2033;PHA[i][2]=5;Disc[i][2]=7;
         filename[i][3]=2034;PHA[i][3]=7;Disc[i][3]=7;
         filename[i][4]=2035;PHA[i][4]=7;Disc[i][4]=8;
         filename[i][5]=2036;PHA[i][5]=7;Disc[i][5]=9;
         filename[i][6]=2037;PHA[i][6]=7;Disc[i][6]=10;
        }                                        
      if(i==3 && j<NDisc[i])//T4                 
        {
         filename[i][0]=2038;PHA[i][0]=4;Disc[i][0]=4;
         filename[i][1]=2039;PHA[i][1]=4;Disc[i][1]=5;
         filename[i][2]=2040;PHA[i][2]=4;Disc[i][2]=6;
         filename[i][3]=2041;PHA[i][3]=4;Disc[i][3]=7;
         filename[i][4]=2042;PHA[i][4]=7;Disc[i][4]=7;
         filename[i][5]=2043;PHA[i][5]=7;Disc[i][5]=8;
         filename[i][6]=2044;PHA[i][6]=7;Disc[i][6]=9;
         filename[i][7]=2045;PHA[i][7]=7;Disc[i][7]=10;
         filename[i][8]=2046;PHA[i][8]=7;Disc[i][8]=11;
         filename[i][9]=2047;PHA[i][9]=7;Disc[i][9]=12;
        }    
       if(i==4 && j<NDisc[i])//guard                 
        {
         filename[i][0]=2060;PHA[i][0]=3;Disc[i][0]=3;
         filename[i][1]=2061;PHA[i][1]=3;Disc[i][1]=4;
         filename[i][2]=2062;PHA[i][2]=3;Disc[i][2]=5;
         filename[i][3]=2063;PHA[i][3]=3;Disc[i][3]=6;
         filename[i][4]=2064;PHA[i][4]=3;Disc[i][4]=7;
         filename[i][5]=2065;PHA[i][5]=7;Disc[i][5]=7;
         filename[i][6]=2066;PHA[i][6]=7;Disc[i][6]=8;
         filename[i][7]=2067;PHA[i][7]=7;Disc[i][7]=9;
         filename[i][8]=2068;PHA[i][8]=7;Disc[i][8]=10;
        }    
          
          
    }//j
  }//i
  
   cout << "After set up filenam etc.... "  <<endl;

//////////////////////////// 
//Trigger related histograms
//////////////////////////// 
 
 string XtitlePH="Pulse height";
 string YtitlePH="Entries";
 TH1F***AllPH=new TH1F**[N];

 cout << "Before  definition of histograms "  <<endl;

 for(int i=0;i<N;i++)
   {  
    AllPH[i]=new TH1F*[NDiscmax];
    for(int j=0;j<NDisc[i];j++)
      {
       AllPH[i][j]=new TH1F(Form("%s (PHA,LOGIC)=(%d,%d)",Trig[i].c_str(),PHA[i][j],(int)Disc[i][j]),Form("%s (PHA,LOGIC)=(%d,%d)",Trig[i].c_str(),PHA[i][j],(int)Disc[i][j]),512,0,512);
       AllPH[i][j]->GetXaxis()->SetTitle(Form("%s",XtitlePH.c_str()));
       AllPH[i][j]->GetYaxis()->SetTitle(Form("%s",YtitlePH.c_str()));
       AllPH[i][j]->SetLineColor(kBlack);
       AllPH[i][j]->SetLineWidth(1);
       AllPH[i][j]->GetYaxis()->SetTitleSize(0.05);
       AllPH[i][j]->GetXaxis()->SetTitleSize(0.05);
       AllPH[i][j]->GetYaxis()->SetLabelSize(0.05);
       AllPH[i][j]->GetXaxis()->SetLabelSize(0.05);  
       AllPH[i][j]->GetXaxis()->SetRangeUser(0,500);
      }
   } //i  
    
    
  cout << "After  definition of histograms "  <<endl;
   
    
 //Create TChain to merge the files for a given energy
 TChain*chain;

 cout << "Before loop over chains"  <<endl;
 
 for(int i=0;i<N;i++)
   {  
    for(int j=0;j<NDisc[i];j++)
      {
      
       chain=new TChain("Data");  
       chain->Add(Form("../Data/NL%4d.BPD.EVENT.root",filename[i][j]));
       //Define variables to read event
       ALEvent *e = new ALEvent();      
       //Set address to access event data
       chain->SetBranchAddress("event",&e); 
       // Get number of event in Tree
       int nentries=chain->GetEntries();
       cout << "Number  of events: " << nentries << endl;  
       
       for (int k=0;k<nentries;k++)
         {
          chain->GetEntry(k); //Load the entry i in the variable e 

          if(k%1000==0) cout << "Event: " << k <<endl;

          ////////////////////////////////////
          //TRIGGER HEIGHTS 
          ////////////////////////////////////  
    
          if(i==0)AllPH[i][j]->Fill(e->get_EneT1().at(0));
          if(i==1)AllPH[i][j]->Fill(e->get_EneT2().at(0));
          if(i==2)AllPH[i][j]->Fill(e->get_EneT3().at(0));
          if(i==3)AllPH[i][j]->Fill(e->get_EneT4().at(0));
          if(i==4)AllPH[i][j]->Fill(e->get_Eneg().at(0));
    
         }//k
      }//j
   } //i
   
   cout << "After loop over chains"  <<endl;

    
//////////////////////////// 
//Trigger related histograms
//////////////////////////// 
 
 string Coin[5]={"T3&T4","T1&T4","T1&T4","T1&T3","T1&T4"};
 
 TH1F**UAllPH=new TH1F*[N];
 string UXtitlePH="Pulse height";
 string UYtitlePH="Entries";


 for(int i=0;i<N;i++)
   {
    if(i==0||i==4)UAllPH[i]=new TH1F(Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),512,0,512);
    else UAllPH[i]=new TH1F(Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),Form("Unbias pulses %s, Coincidence %s",Trig[i].c_str(),Coin[i].c_str()),512,0,512);
    UAllPH[i]->GetXaxis()->SetTitle(Form("%s",UXtitlePH.c_str()));
    UAllPH[i]->GetXaxis()->SetTitle(Form("%s",UXtitlePH.c_str()));
    UAllPH[i]->GetYaxis()->SetTitle(Form("%s",UYtitlePH.c_str()));
    UAllPH[i]->GetYaxis()->SetTitleSize(0.05);
    UAllPH[i]->GetXaxis()->SetTitleSize(0.05);
    UAllPH[i]->GetYaxis()->SetLabelSize(0.05);
    UAllPH[i]->GetXaxis()->SetLabelSize(0.05);  
    UAllPH[i]->SetLineColor(kBlack);
    UAllPH[i]->SetLineWidth(1);
    UAllPH[i]->GetXaxis()->SetRangeUser(0,500);
   } //i  
    
 //Create TChain to merge the files for a given energy
 TChain**chainU=new TChain*[3];
 
 for(int i=0;i<N;i++)
   {
    chainU [i]=new TChain("Data");  
    if(i==0 )
     {
      chainU[i]->Add(Form("../Data/NL2013.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL2028.BPD.EVENT.root"));
     }
    else if(i==1)
     {
      chainU[i]->Add(Form("../Data/NL0212.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0214.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0218.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0221.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0250.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0252.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0253.BPD.EVENT.root"));
      chainU[i]->Add(Form("../Data/NL0254.BPD.EVENT.root"));
     }
    else chainU[i]->Add(Form("../Data/NL2014.BPD.EVENT.root"));
       
   }

 int* binmax=new int[N];
 double* x=new double[N];
 for(int i=0;i<N;i++)
   {  
    //Define variables to read event
    ALEvent *e = new ALEvent();    
    int kk=0;
    if(i==1||i==2||i==4) kk=1;
    if(i==3) kk=2;
    //Set address to access event data
    chainU[kk]->SetBranchAddress("event",&e); 
 
    // Get number of event in Tree
    int nentries=chainU[kk]->GetEntries();
    cout << "Number  of events: " << nentries << endl;  
      
    for (int j=0;j<nentries;j++)
      {
       chainU[kk]->GetEntry(j); //Load the entry i in the variable e 

       if(j%1000==0) cout << "Event: " << j <<endl;

       ////////////////////////////////////
       //TRIGGER HEIGHTS 
       ////////////////////////////////////  
    
       if(i==0)UAllPH[i]->Fill(e->get_EneT1().at(0));
       if(i==1)UAllPH[i]->Fill(e->get_EneT2().at(0));
       if(i==2)UAllPH[i]->Fill(e->get_EneT3().at(0));
       if(i==3)UAllPH[i]->Fill(e->get_EneT4().at(0));
       if(i==4)UAllPH[i]->Fill(e->get_Eneg().at(0));
    
      }//j
    binmax[i] = UAllPH[i]->GetMaximumBin(); 
    x[i] = UAllPH[i]->GetXaxis()->GetBinCenter(binmax[i]);
    cout << Trig[i]<<", Maximum PHA (Bin center): "<< x[i] << endl;    

   } 
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 //gStyle->SetPalette( kDarkBodyRadiator);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
//////////////////////////////    
  
 TCanvas** CanPHNorm=new TCanvas*[N];
 TCanvas** CanPH=new TCanvas*[N];
 TLegend** leg= new TLegend*[N];
 TGraph** tg=new TGraph*[N];
 TLine**Line=new TLine*[N];
 
 float** cutoff=new float*[N];

 TF1** fitT=new TF1*[N];
 TPad** subpad = new TPad*[N];
 TLegend** legfit =new TLegend*[N];

 TH1 **frame = new TH1*[N];

 float* y=new float[N];
 
 y[0]=0.04;
 y[1]=0.12;
 y[2]=0.1;
 y[3]=0.085;
 y[4]=0.1;
 
 for(int i=0;i<N;i++)
   {
    //CanPHNorm[i]= new TCanvas(Form("Pulse heights %s",Trig[i].c_str()),Form("Pulse heights %s",Trig[i].c_str()),200,10,1200,800); 
    CanPH[i]= new TCanvas(Form("Pulse heights bis %s",Trig[i].c_str()),Form("Pulse heights bis%s",Trig[i].c_str()),200,10,1700,1000); 
    //{8,9,7,10,9}
    if (i==2) CanPH[i]->Divide(3,3);    
    else CanPH[i]->Divide(4,3);    
    leg[i]=new TLegend(0.3,0.65,0.9,0.9);
    leg[i]->SetBorderSize(0);leg[i]->SetFillStyle(0);
    leg[i]->SetNColumns(3);
       
    //CanPHNorm[i]->cd();   
    
   // frame[i] = new TH1F(Form("frame %d",i),Form("frame %d",i),1000,1,200);
   // frame[i]->SetMaximum(y[i]);
   // frame[i]->Draw();
    CanPH[i]->cd(1);
    gPad->SetLeftMargin(0.15);

    //UAllPH[i]->GetXaxis()->SetRangeUser(1.,200);
    UAllPH[i]->Draw();
    //UAllPH[i]->GetYaxis()->SetRangeUser(0.,y[i]);

    for(int j=0;j<NDisc[i];j++)
      {
       CanPH[i]->cd(j+2);  
       gPad->SetLeftMargin(0.15);
      // gPad->SetBottomMargin(0.15);
       // leg[i]->AddEntry(AllPH[i][j],Form("%s (PHA,LOGIC)=(%d,%d)",Trig[i].c_str(),PHA[i][j],(int)Disc[i][j]),"l");
       AllPH[i][j]->GetXaxis()->SetRangeUser(1.,200);
       //AllPH[i][j]->GetYaxis()->SetRangeUser(0.,1000);
       AllPH[i][j]->SetLineColor(color[j]);
       //AllPH[i][j]->SetLineWidth(2);
       //AllPH[i][j]->DrawNormalized("same");
       AllPH[i][j]->Draw();
       gPad->Update();
      }//j   
    //leg[i]->AddEntry(UAllPH[i],Form("%s (7,7), Coin %s",Trig[i].c_str(),Coin[i].c_str()),"l");  
   // leg[i]->Draw("same");  
    Line[i]=new TLine(x[i],0.,x[i],0.6*y[i]);
    Line[i]->SetLineColor(kGreen+3);
    Line[i]->SetLineWidth(2);
   // Line[i]->Draw("same");
    
    cutoff[i]=new float[NDisc[i]];
    for(int j=0;j<NDisc[i];j++)
      {
       cutoff[i][j]=0;   
       float tmpt=0; 
       tmpt=AllPH[i][j]->GetBinContent(1);
       for(int k=2;k<AllPH[i][j]->GetNbinsX()-1;k++)
         {
          if( AllPH[i][j]->GetBinContent(k) - tmpt >10) {cutoff[i][j]=AllPH[i][j]->GetXaxis()->GetBinCenter(k);k=AllPH[i][j]->GetNbinsX();}
          else tmpt=AllPH[i][j]->GetBinContent(k);
         }//k
      }//j
    tg[i]=new TGraph(NDisc[i],&Disc[i][0],&cutoff[i][0]);
    tg[i]->SetTitle("");
    tg[i]->SetMarkerSize(1);
    tg[i]->SetMarkerStyle(20);
    tg[i]->SetMarkerColor(kRed);
    tg[i]->GetXaxis()->SetTitle("LOGIC Disc");
    tg[i]->GetYaxis()->SetTitle("PHA Cutoff");
    tg[i]->GetYaxis()->SetTitleSize(0.05);
    tg[i]->GetXaxis()->SetTitleSize(0.05);
    tg[i]->GetYaxis()->SetLabelSize(0.05);
    tg[i]->GetXaxis()->SetLabelSize(0.05);  
          
    fitT[i]=new TF1(Form("fit%s",Trig[i].c_str()),"pol1",0,12);
    tg[i]->Fit(fitT[i]);
    fitT[i]->SetLineColor(kRed);
    
   // subpad[i] = new TPad(Form("subpad %s",Trig[i].c_str()),Form("subpad %s",Trig[i].c_str()),0.53,0.25,0.93,0.65); 
   // subpad[i]->SetFillStyle(0);
   // subpad[i]->Draw(); 
   // subpad[i]->cd();
    CanPH[i]->cd(NDisc[i]+2);  
    gPad->SetLeftMargin(0.15);
//    gPad->SetBottomMargin(0.15);
    //subpad[i]->SetFillStyle(0);
    tg[i]->Draw("ap");
    
    legfit[i]=new TLegend(0.25,0.7,0.75,0.9);
    legfit[i]->SetBorderSize(0);legfit[i]->SetFillStyle(0);
    legfit[i]->AddEntry(fitT[i],Form("s=%3.2f, i=%3.2f",fitT[i]->GetParameter(1),fitT[i]->GetParameter(0)),"l");
    legfit[i]->Draw("same");
    
    
    
    
     if(i==0)CanPH[i]->Print(Form("PulseHeight_Disc.pdf("),Form("Title:%s",Trig[i].c_str()));    
     else if(i==N-1)CanPH[i]->Print(Form("PulseHeight_Disc.pdf)"),Form("Title:%s",Trig[i].c_str()));    
     else CanPH[i]->Print(Form("PulseHeight_Disc.pdf"),Form("Title:%s",Trig[i].c_str()));    

    
    
    
   }//i
     
}
