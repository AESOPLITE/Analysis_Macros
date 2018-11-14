void LoadMacros()

{

//Load shared libraries and compile them if necessary
//You can have "++" to force to compile the code every time

gROOT->ProcessLine(".L ALTckhit.cxx+");
gROOT->ProcessLine(".L ALEvent.cxx+");
gROOT->ProcessLine(".L LoadDataparameters.cxx+");
gROOT->ProcessLine(".L LoadMCparameters.cxx+");

//Load the main analysis files, compile if needed
//You can have "++" to force to compile the code every time

//gROOT->ProcessLine(".L AnalyseData.C+");
//gROOT->ProcessLine(".L Analysis_Reconstruction_MC.C+");
gROOT->ProcessLine(".L Analysis_PR_MC_OneEnergy.C+");
//gROOT->ProcessLine(".L AnalysisTestOneEnergy.C+");
//gROOT->ProcessLine(".L Analysis_Pulls_MC.C+");
gROOT->ProcessLine(".L EfficienciesMC.C+");
gROOT->ProcessLine(".L EfficienciesData.C+");
gROOT->ProcessLine(".L Calibrate_Flight.C+");
gROOT->ProcessLine(".L GrowthCurves.C+");
//gROOT->ProcessLine(".L TimeCuts.C+");
//gROOT->ProcessLine(".L TimeSeries.C+");
//gROOT->ProcessLine(".L DisplayEventsFlightRuns2018.C+");
gROOT->ProcessLine(".L EnergyDeposit.C+");
gROOT->ProcessLine(".L EneDep2D.C+");
gROOT->ProcessLine(".L Synthesis_MC.C+");
gROOT->ProcessLine(".L DisplayEventsMC.C+");
gROOT->ProcessLine(".L DisplayEventsGIF.C+");


}
