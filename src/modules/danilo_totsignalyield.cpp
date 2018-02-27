//Creating Histograms for GGM scan total signal yield and total stat. uncertainty
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph2D.h>
#include <cmath>


extern "C"
void run()
{     
	TFile file("../output/SignalTrees/SignalTree_GGM_M1_M2_inclusiv.root","read");
	TH2F* bin0 = (TH2F*)file.Get("SignalRate_bin0");
	TH2F* bin1 = (TH2F*)file.Get("SignalRate_bin1");
	TH2F* bin2 = (TH2F*)file.Get("SignalRate_bin2");
	TH2F* bin3 = (TH2F*)file.Get("SignalRate_bin3");
	TH2F* e_bin0 = (TH2F*)file.Get("SignalStat_bin0");
	TH2F* e_bin1 = (TH2F*)file.Get("SignalStat_bin1");
	TH2F* e_bin2 = (TH2F*)file.Get("SignalStat_bin2");
	TH2F* e_bin3 = (TH2F*)file.Get("SignalStat_bin3");
	
	TH2F* sum = new TH2F(*bin0);
	sum->Add(bin1);
	sum->Add(bin2);
	sum->Add(bin3);
	
	
	TH2F* e_sum = new TH2F(*bin0);
	float e_temp = 0;
	for (int i=1; i<=e_bin0->GetXaxis()->GetNbins()+1; i++) {
		for (int j=1; j<=e_bin0->GetYaxis()->GetNbins()+1; j++) {
			e_temp = sqrt(pow(bin0->GetBinContent(i,j)*(e_bin0->GetBinContent(i,j)-1),2.0)+pow(bin1->GetBinContent(i,j)*(e_bin1->GetBinContent(i,j)-1),2.0)+pow(bin2->GetBinContent(i,j)*(e_bin2->GetBinContent(i,j)-1),2.0)+pow(bin3->GetBinContent(i,j)*(e_bin3->GetBinContent(i,j)-1),2.0));
			e_temp = e_temp/sum->GetBinContent(i,j);
			e_sum->SetBinContent(i,j,e_temp);
			std::cout<<e_temp<<"  "<<e_sum->GetBinContent(i,j)<<std::endl;
			
		}
	}
	
	TFile out("../output/stuff/GGM_signalYieldTot.root","update");
	sum->Write("TotalYield",TObject::kOverwrite);
	e_sum->Write("TotalStatError",TObject::kOverwrite);
	out.Close();
	file.Close();
}
