////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, Februay 9 , 2018
////////////////////////////////////////////////////////////////////////////////////////// 

#include "headers.h"
#include "ALEvent.h"
#include "LoadDataparameters.h"
#include "TChain.h"
#include "TPaveStats.h"
#include "TGaxis.h"


ClassImp(ALTckhit)
ClassImp(ALEvent)


void AnalyseFlight2018(string s,string c,int geoconf);

void PullDistribution(string s,string c,int geoconf);

void AnalyseFlight2018(string s,string c,int geoconf)
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

 //output file
 ofstream outputTXT;
 
 outputTXT.open(Form("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly/18A1_000to022.EVENT_PRonly.txt"));
    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");

 chain->Add(Form("/home/sarah/AESOPLITE/FlightData/BBR2/PRonly/18A1_000to022.EVENT_PRonly_Calibrated.root"));
    
 //Define variables to read event
 ALEvent *e = new ALEvent();      
 //Set address to access event data
 chain->SetBranchAddress("Calevent",&e); 
 
 // Get number of event in Tree
 int nentries=chain->GetEntries();
 cout << "Number  of events: " << nentries << endl;  
      
 outputTXT << "Event"  << " " ;
 outputTXT << "tPHA"  << " " ;
 outputTXT << "TempB1" << " " ;
 outputTXT << "TempB2" << " " ;
 outputTXT << "PressB1" << " " ;
 outputTXT << "PressB2" << " " ;
 outputTXT << "T1"  << " " ;
 outputTXT << "T2"  << " " ;
 outputTXT << "T3"  << " " ;
 outputTXT << "T4"  << " " ;
 outputTXT << "g"  << " " ;
 outputTXT << "GOPHA"  << " " ;
 outputTXT << "GOEVT"  << " " ;
 outputTXT << "Nhits"  << " " ;
 outputTXT << "Ti"  << " " ;
 outputTXT << "deflecPR"  << " " ;
 outputTXT << "p0PR"  << " " ;
 outputTXT << "chi2NB"  << " " ;
 outputTXT << "chi2B"  << " " ;
 outputTXT << "slopePR"  << " " ;
 outputTXT << "interPR"  << " " ;
 
 outputTXT <<endl;

 
 
 for (int j=0;j<nentries;j++)
   {
    chain->GetEntry(j); //Load the entry i in the variable e 

    if(j%10000==0) cout << "Event: " << j <<endl;
    
    ////////////////////////////////////
    //TRIGGER CUT
    ////////////////////////////////////  
    outputTXT << Form("%5d",e->get_eventnumber())  << " "; 
    outputTXT << Form("%8.6f",e->get_dPHA()+e->get_hPHA()/24.+e->get_miPHA()/(24.*60.)+e->get_sPHA()/(24.*3600.))  << " "; 
    outputTXT << Form("%5.2f",e->get_TempB1())  << " "; 
    outputTXT << Form("%5.2f",e->get_TempB2())  << " "; 
    outputTXT << Form("%5.2f",e->get_PressB1())  << " "; 
    outputTXT << Form("%5.2f",e->get_PressB2())  << " "; 
    outputTXT << Form("%4d",(int)e->get_EneT1().at(0))  << " "; 
    outputTXT << Form("%4d",(int)e->get_EneT2().at(0))  << " "; 
    outputTXT << Form("%4d",(int)e->get_EneT3().at(0))  << " "; 
    outputTXT << Form("%4d",(int)e->get_EneT4().at(0))  << " "; 
    outputTXT << Form("%4d",(int)e->get_Eneg().at(0))  << " "; 
    outputTXT << Form("%5d",(int)e->get_GoPHA())  << " "; 
    outputTXT << Form("%5d",(int)e->get_GoEVT())  << " "; 
    outputTXT << Form("%5d",(int)e->get_Nhits())  << " "; 
    outputTXT << Form("%3d",(int)e->get_Ti())  << " "; 
    outputTXT << Form("%6.1f",1000*e->get_deflecPR())  << " "; 
    outputTXT << Form("%6.1f",1000*e->get_p0PR())  << " "; 
    outputTXT << Form("%5.2f",e->get_chi2NBPR())  << " "; 
    outputTXT << Form("%5.2f",e->get_chi2BPR())  << " "; 
    outputTXT << Form("%5.2f",e->get_slopePR())  << " "; 
    outputTXT << Form("%5.2f",e->get_interPR())  << " "; 
    outputTXT <<endl;
   }//j
   
 outputTXT.close();
//////////////////////////////   
// Display  
//////////////////////////////   
 //gStyle->SetOptStat("emruo");
 //gStyle->SetPalette( kDarkBodyRadiator);
 
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
 TH1F*Deflection=new TH1F("Deflection","Deflection",100,-0.1,0.1);
 
 for(int j=0;j<7;j++)
   {
    XReso[j]=new TH1F(Form("Xreso L%d",j),Form("Xreso L%d",j),100,-1,1); 
    YReso[j]=new TH1F(Form("Yreso L%d",j),Form("Yreso L%d",j),100,-1,1);
    ZReso[j]=new TH1F(Form("Zreso L%d",j),Form("Zreso L%d",j),100,-3,3);
    XResoX[j]=new TH2F(Form("Xreso vs X L%d",j),Form("Xreso vs X L%d",j),1000,-100,100,400,-1,1); 
    YResoY[j]=new TH2F(Form("Yreso vs Y L%d",j),Form("Yreso vs Y L%d",j),1000,-100,100,400,-1,1);
    if(j==0||j==4||j==6)ZResoXY[j]=new TH2F(Form("Zreso vs X L%d",j),Form("Zreso vs X L%d",j),1000,-100,100,400,-2,2);
    else ZResoXY[j]=new TH2F(Form("Zreso vs Y L%d",j),Form("Zreso vs Y L%d",j),1000,-100,100,400,-2,2);
   }
    
 //Create TChain to merge the files for a given energy
 TChain*chain=new TChain("Data");

  chain->Add(Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/18A1_000to022.EVENT_PRonly.root"));  
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
    if(e->get_Ti()!=127)      continue;
    
    //T1&T3&T4
    if(e->get_EneT1().at(0)<300 || e->get_EneT1().at(0)>500)      continue;
    if(e->get_EneT2().at(0)<180)        continue;
    if(e->get_EneT3().at(0)<300 || e->get_EneT3().at(0)>500)        continue;
    if(e->get_EneT4().at(0)<300 || e->get_EneT4().at(0)>500)        continue;
    if(e->get_Eneg().at(0)>0 )        continue;

      
    
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
    ////////////////////////////////////
    //DEFLECTION 
    ////////////////////////////////////      
      
     Deflection->Fill(e->get_deflecPR());

     
      
      
   }//j
   
   
//////////////////////////////   
// Display  
//////////////////////////////   
 gStyle->SetOptStat("emruo");
 gStyle->SetPalette( kDarkBodyRadiator);

//////////////////////////////    
 
 string PDFFile=Form("/home/psm/Documents/UDEL/AESOPLITE/Esrange2018/FlightData/BBR2/18A1_SplitBPD/PullDistribution.pdf");

//////////////////////////////    
 TCanvas**CanXYZReso=new TCanvas*[7];  
 int Nsig=5;

 for(int j=0;j<7;j++)
   {
    CanXYZReso[j]=new TCanvas(Form("CanXYZReso L%d",j),Form("CanXYZReso L%d",j),200,10,1200,1200);     
    CanXYZReso[j]->Divide(2,1);
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
    
//if(j==6)CanXYZReso[j]->Print(Form("%s)",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
    if(j==0)CanXYZReso[j]->Print(Form("%s(",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
    else CanXYZReso[j]->Print(Form("%s",PDFFile.c_str()),Form("Title:Resolution L%d",j));    
     
   }//j   
    
   TCanvas*canDeflec=new TCanvas("canDeflec","canDeflec",200,10,800,800);     
   canDeflec->cd();
   Deflection->GetYaxis()->SetTitle("entries");             
   Deflection->GetXaxis()->SetTitle("deflection (rad)");             
   Deflection->Draw();
   canDeflec->Print(Form("%s)",PDFFile.c_str()),Form("Title:Deflection"));    
     
    
}


