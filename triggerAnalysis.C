

#include "headers.h"
#include "ALEvent.h"
#include <cmath>

ClassImp(ALTckhit)
ClassImp(ALEvent)

class LineFit {
		double slope; 
		double intercept;
	  double chi2;
	public:
		LineFit(int nPnts, double x[], double y[]) {
      if (nPnts < 2) {
        chi2 = -999.;
        return;
      }
      double res = 0.0228/sqrt(12.0);
			double xbar=0, ybar=0, x2bar=0, xybar=0;
			for (int i=0; i<nPnts; ++i) {
				xbar += x[i];
				ybar += y[i];
				x2bar += x[i]*x[i];
				xybar += x[i]*y[i];
			}
			xbar = xbar/nPnts;
			ybar = ybar/nPnts;
			x2bar = x2bar/nPnts;
			xybar = xybar/nPnts;
      double det = x2bar - xbar*xbar;
      if (det == 0.) {
        cout << "LineFit: zero determinate" << endl;
        chi2 = -999.;
        return;
      }
			slope = (xybar - xbar*ybar)/det;
			intercept = ybar - slope*xbar;
			chi2 = 0.;
			for (int i=0; i<nPnts; ++i) {
				chi2 += pow((intercept + slope*x[i] - y[i])/res,2);
			}
		}
		double eval(double x) {return (intercept + slope*x);}
		double getSlope() {return slope;}
		double getIntercept() {return intercept;}
		double getChi2() {return chi2;}
};

class ParabolaFit {
		double a, b, c;
		double chi2;
	public:
		ParabolaFit(int nPnts, double x[], double y[]) {

	    double sum = 0.0;
	    double sumx = 0.0;
	    double sumxy = 0.0;
	    double sumx2y = 0.0;
	    double sumy = 0.0;
	    double sumx2 = 0.0;
	    double sumx3 = 0.0;
	    double sumx4 = 0.0;
      double res = 0.0228/sqrt(12.0);
      
	    for (int ii = 0; ii < nPnts; ii++) {
		    sum = sum + 1.0;
		    sumx = sumx + x[ii];
		    sumx2 = sumx2 + x[ii] * x[ii];
		    sumx3 = sumx3 + x[ii] * x[ii] * x[ii];
		    sumx4 = sumx4 + x[ii] * x[ii] * x[ii] * x[ii];
		    sumy = sumy + y[ii];
		    sumxy = sumxy + x[ii] * y[ii];
		    sumx2y = sumx2y + x[ii] * x[ii] * y[ii];
	    }

		  double det = sum*sumx2*sumx4 + sumx*sumx3*sumx2 + sumx2*sumx*sumx3 - (sumx2*sumx2*sumx2 + sumx3*sumx3*sum + sumx4*sumx*sumx);
		  double a1 = sumy*sumx2*sumx4 + sumx*sumx3*sumx2y + sumx2*sumxy*sumx3 - (sumx2y*sumx2*sumx2 + sumx3*sumx3*sumy + sumx4*sumxy*sumx);
		  double a2 = sum*sumxy*sumx4 + sumy*sumx3*sumx2 + sumx2*sumx*sumx2y - (sumx2*sumxy*sumx2 + sumx2y*sumx3*sum + sumx4*sumx*sumy);
		  double a3 = sum*sumx2*sumx2y + sumx*sumxy*sumx2 + sumy*sumx*sumx3 - (sumx2*sumx2*sumy + sumx3*sumxy*sum + sumx2y*sumx*sumx);

		  if (det == 0.) {
			  cout << "ParabolaFit: zero determinate.  Abort!" << endl;
        chi2 = -999.;
			  return;
		  }

      a = a1 / det;
		  b = a2 / det;
		  c = a3 / det;
			chi2 = 0.;
			for (int i=0; i<nPnts; ++i) {
				chi2 += pow((a + (c*x[i] + b)*x[i] - y[i])/res,2);
			}
		}
		double eval(double x) {return (a + (c*x + b)*x);}
		double getA() {return a;}
		double getB() {return b;}
    double getC() {return c;}
		double getChi2() {return chi2;}
};

void triggerAnalysis() {
	
//Load all ROOT files onto a TChain

TChain*chain=new TChain("Data");
//chain->Add(Form("/home/robert/aesoplite/FlightData/BBR2/PRonly/18A1_007.BPD.EVENT_PRonly.root"));
//chain->Add(Form("/home/robert/aesoplite/FlightData/BBR2/PRonly/*.root"));
chain->Add(Form("/home/robert/aesoplite/GroundData/NL2268.BPD.EVENT_PRonly.root"));
ALEvent *e = new ALEvent();
chain->SetBranchAddress("event",&e); 
int entries=chain->GetEntries();	
cout << "Number  of events: " << entries << endl;  

//Histogram definitions
TH1F* hDay=new TH1F("Day","Day of the month",40,0,40);
TH1F* hHour=new TH1F("Hour","Hour of the day",26,0,26);
TH1F* hTrgPat=new TH1F("TrgPat","Trigger Pattern",10,0,10); 
TH1F* hTrgPatB=new TH1F("TrgPatB","Trigger pattern for expected trigger in bending view",10,0,10); 
TH1F* hTrgPatNB=new TH1F("TrgPatNB","Trigger pattern for expected trigger in non-bending view",10,0,10); 
TH1F* hTrgPatB0=new TH1F("TrgPatB0","Trigger pattern for trigger expected only in bending view",10,0,10); 
TH1F* hTrgPatNB0=new TH1F("TrgPatNB0","Trigger pattern for trigger expected only in non-bending view",10,0,10); 
TH1F* hL = new TH1F("L","Hit Layer Index",10,0,10);
TH1F* hB = new TH1F("B","Number of bending view hits",10,0,10);
TH1F* hNB = new TH1F("NB","Number of non-bending view hits",10,0,10);
TH1F* hBL = new TH1F("BL","Number of bending view trigger layers hit for pattern = 2 or 3",10,0,10);
TH1F* hNBL = new TH1F("NBL","Number of non-bending view trigger layers hit for pattern = 1 or 3",10,0,10);
TH1F* hBL1 = new TH1F("BL1","Number of bending view trigger layers hit for pattern = 1",10,0,10);
TH1F* hNBL1 = new TH1F("NBL1","Number of non-bending view trigger layers hit for pattern = 2",10,0,10);
TH1F* hBL2 = new TH1F("BL2","Number of bending vieneww trigger layers hit for pattern = 2",10,0,10);
TH1F* hNBL2 = new TH1F("NBL2","Number of non-bending view trigger layers hit for pattern = 1",10,0,10);
TH1F* hBL3 = new TH1F("BL3","Number of bending view trigger layers hit for pattern = 0",10,0,10);
TH1F* hNBL3 = new TH1F("NBL3","Number of non-bending view trigger layers hit for pattern = 0",10,0,10);
TH1F* hYbend = new TH1F("Ybend","Y coordinate in the bending view",100,-15.,15.);
TH1F* hZbend = new TH1F("Zbend","Z coordinate in the bending view",100,-25.,5.);
TH1F* hXnonBend = new TH1F("XnonBend","X coordinate in the non-bending view",100,-15.,15.);
TH1F* hZnonBend = new TH1F("ZnonBend","Z coordinate in the non-bending view",100,-25.,5.);
TH1F* hChi2Bend = new TH1F("Chi2Bend","chi^2 of parabola fit in the bending view",100,0.,1000.);
TH1F* hChi2nBend = new TH1F("Chi2nBend","chi^2 of line fit in the non-bending view",100,0.,1000.);

//loop through all entries
for(int i=0; i<entries; i++) {
  chain->GetEntry(i);

  if(i%10000==0) cout << "Event " << i << endl;
	
	//Trigger information and layers information
 	//Trigger pattern info
  int Day = e->get_dPHA();
  int Hour = e->get_hPHA();
  hDay->Fill(Day);
  hHour->Fill(Hour);
  //if (Day > 16) break;
  //if (Day == 16 && Hour > 21) break;
	int PatternEVT = e->get_PatternEVT();
  if (PatternEVT < 0) continue;  //Not an event
  hTrgPat->Fill(PatternEVT);
  int  nhits = e->get_Nhits();    //total number of hits
  int NL= e->get_NLayers();       //number of layers with hits
  //Number of layers with hit(s) in bending/non-bending plane
  int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
  int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	
	//Reconstructed energy 
	float pPR=1000*fabs(e->get_p0PR());   //in MeV
	double deflecPR = e->get_deflecPR();
	float P0sign = pPR*TMath::Sign(1,deflecPR);
	double chiNBPR = e->get_chi2NBPR();
	double chiBPR = e->get_chi2BPR();	
	
	//PMT trigger info
	bool T1=e->get_T1();
	bool T2=e->get_T2();
	bool T3=e->get_T3();
	bool T4=e->get_T4();
	bool G=e->get_guard();
	
  int nBend = 0;
  int nNonBend = 0;
  double xf[7], yf[7], x2f[7], y2f[7];
  bool lHit[7] = {false,false,false,false,false,false,false};
	for(int k=0;k<nhits;k++) {
		int Lindex=(int)e->get_hits().at(k)->get_L();  //from layer L0 on top to L6
    if (Lindex < 0 || Lindex > 6) continue;
    hL->Fill(Lindex);

    //hit coordinates
    double zreco = e->get_hits().at(k)->get_z();
    if (zreco < -990.) continue;
    double xreco = e->get_hits().at(k)->get_x();
	  double yreco = e->get_hits().at(k)->get_y();

    //cout << k << " " << Lindex << " " << xreco << " " << yreco << " " << zreco << endl;     
    if (Lindex == 0 || Lindex == 4 || Lindex == 6) {
      if (nNonBend < 6) nNonBend++;
      hXnonBend->Fill(xreco);
      hZnonBend->Fill(zreco);
      x2f[nNonBend] = zreco;
      y2f[nNonBend] = xreco;
    } else {
      if (nBend < 6) nBend++;
      xf[nBend] = zreco;
      yf[nBend] = yreco;
      hYbend->Fill(yreco);
      hZbend->Fill(zreco);
    }
    lHit[Lindex] = true;
  }	//end loop over hits
  hB->Fill(nBend);
  hNB->Fill(nNonBend);
 
  int nBendL = 0;
  int nBendLT = 0;
  int nNonBendL = 0;
  for (int lyr=0; lyr<7; lyr++) {
    if (lyr == 0 || lyr == 4 || lyr == 6) {
      if (lHit[lyr]) nNonBendL++;
    } else {
      if (lHit[lyr]) {
        nBendL++;
        if (lyr != 5) nBendLT++;
      }
    }
  }
  if (nBend == 4 && nBendL == nBend) {
    ParabolaFit *pFit = new ParabolaFit(nBend,xf,yf);
    hChi2Bend->Fill(pFit->getChi2());
    delete pFit;
  }
  if (nNonBend == 3 && nNonBendL == nNonBend) {
    LineFit *lFit = new LineFit(nNonBend,x2f,y2f);
    hChi2nBend->Fill(lFit->getChi2());
    delete lFit;
  }
    
  bool tBend = lHit[1] && lHit[2] && lHit[3];
  bool tNonBend = lHit[0] && lHit[4] && lHit[6];
  if (tNonBend) hTrgPatNB->Fill(PatternEVT);
  if (tNonBend && !tBend) hTrgPatNB0->Fill(PatternEVT);
  if (tBend) hTrgPatB->Fill(PatternEVT);
  if (tBend && !tNonBend) hTrgPatB0->Fill(PatternEVT);
  if (PatternEVT == 2 || PatternEVT == 3) hBL->Fill(nBendLT);
  if (PatternEVT == 1 || PatternEVT == 3) hNBL->Fill(nNonBendL);
  if (PatternEVT == 2) hNBL1->Fill(nNonBendL);
  if (PatternEVT == 1) hBL1->Fill(nBendLT);
  if (PatternEVT == 2) hBL2->Fill(nBendLT);
  if (PatternEVT == 1) hNBL2->Fill(nNonBendL);
  if (PatternEVT == 0) hBL3->Fill(nBendLT);
  if (PatternEVT == 0) hNBL3->Fill(nNonBendL);
} //end loop over entries

  TCanvas* c0 = new TCanvas("c0","c0");
  c0->Divide(2,2);
  c0->cd(1);
  hDay->Draw();
  c0->cd(2);
  hHour->Draw();
  c0->cd(3);
  hChi2Bend->Draw();
  c0->cd(4);
  hChi2nBend->Draw();

	TCanvas* c1 = new TCanvas("c1","c1");
	c1->Divide(2,2);
	c1->cd(1);
  hTrgPat->Draw();
  c1->cd(2);
  hL->Draw();
  c1->cd(3);
  hB->Draw();
  c1->cd(4);
  hNB->Draw();
  
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);
  c2->cd(1);
  hTrgPatNB->Draw();
  c2->cd(2);
  hTrgPatB->Draw();
  c2->cd(3);
  hTrgPatNB0->Draw();
  c2->cd(4);
  hTrgPatB0->Draw();
  
  TCanvas* c3 = new TCanvas("c3","c3");
  c3->Divide(2,2);
  c3->cd(1);
  hBL->Draw();
  c3->cd(2);
  hNBL->Draw();
  c3->cd(3);
  hBL1->Draw();
  c3->cd(4);
  hNBL1->Draw();
  	
  TCanvas* c4 = new TCanvas("c4","c4");
  c4->Divide(2,2);
  c4->cd(1);
  hBL2->Draw();
  c4->cd(2);
  hNBL2->Draw();
  c4->cd(3);
  hBL3->Draw();
  c4->cd(4);
  
  TCanvas* c5 = new TCanvas("c5","c5");
  c5->Divide(2,2);
  c5->cd(1);
  hXnonBend->Draw();
  c5->cd(2);
  hZnonBend->Draw();
  c5->cd(3);
  hYbend->Draw();
  c5->cd(4);
  hZbend->Draw();
};
