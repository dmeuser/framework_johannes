#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRandom3.h>

Config const &cfg=Config::get();
extern "C"

void run(){
	std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
	
	//~ for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
		//~ std::cout<<dss.name<<std::endl;
	//~ }

	hist::Histograms<TH2F> hs2d(vsDatasubsets);

	hs2d.addHist("test", ";bla;bla", {0,10,20,30},{10,10,10}, {0,100,200,300},{100,100,100});
	
	hs2d.setCurrentSample("SinglePhoton_03Feb2017");

	TRandom3 *rand = new TRandom3();

	for (int i=0; i<200; i++){
		hs2d.fill("test",rand->Gaus(15.0,20.0),rand->Gaus(150.0,200.0));
	}
	
	TH2F *test = hs2d.getHistogram("test","SinglePhoton_03Feb2017");
	hs2d.mergeOverflow();	
	TH2F *test_merged = hs2d.getHistogram("test","SinglePhoton_03Feb2017");
	
	//~ for (int i=1; i<=3; i++){
		//~ if (i==1){	
			//~ std::cout<<test->GetBinContent(i,3)+test->GetBinContent(i,4)+test->GetBinContent(i-1,3)<<std::endl;
			//~ std::cout<<test_merged->GetBinContent(i,3)<<"   "<<test_merged->GetBinContent(i,4)<<std::endl;
		//~ }
	//~ }
	std::cout<<test->GetBinContent(0,0)+test->GetBinContent(1,0)+test->GetBinContent(0,1)+test->GetBinContent(1,1)<<std::endl;
	std::cout<<test_merged->GetBinContent(1,1)<<std::endl;
	
}


