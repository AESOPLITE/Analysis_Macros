

#include "headers.h"
#include "ALEvent.h"
#include <cmath>

ClassImp(ALTckhit)
ClassImp(ALEvent)

void Example() {
	
//Load all ROOT files onto a TChain

TChain*chain=new TChain("Data");
chain->Add(Form("/home/robert/aesoplite/FlightData/BBR2/PRonly/*.root"));
ALEvent *e = new ALEvent();
chain->SetBranchAddress("event",&e); 
int entries=chain->GetEntries();	
cout << "Number  of events: " << entries << endl;  

//loop through all entries
for(int i=0; i<entries; i++) 
{
	chain->GetEntry(i);    //load up current event

    if(i%10000==0) cout << "Event " << i << endl;
	
	//functions to get tracker informations
	
	//Trigger information and layers information
    int  nhits = e->get_Nhits();    //total number of hits
	int NL= e->get_NLayers();       //number of layers with hits
	 //Number of layers with hit(s) in bending/non-bending plane
    int NLB = e->get_Layer(1) + e->get_Layer(2) + e->get_Layer(3) + e->get_Layer(5);
    int NLNB =  e->get_Layer(0) + e->get_Layer(4) + e->get_Layer(6);
	
	//Trigger pattern info
	int PatternEVT = e->get_PatternEVT();
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
	
	//To look at hits, you must run a loop through all ALTckhit 
		for(int k=0;k<nnhits;k++)
         {
		int Lindex=(int)e->get_hits().at(k)->get_L();  //from layer L0 on top to L6
      //hit coordinates
       double xreco = e->get_hits().at(k)->get_xreco();
	   double yreco = e->get_hits().at(k)->get_yreco();
	   double zreco = e->get_hits().at(k)->get_zreco();
			
		//directional cosines
	   double cxreco = e->get_hits().at(k)->get_cxreco();
	   double cyreco = e->get_hits().at(k)->get_cyreco();
	   double czreco = e->get_hits().at(k)->get_czreco();
			    }
					
	} //end loop on function
	
}		// end function
				