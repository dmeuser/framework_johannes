#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

/*
 * Plotting multiple limts in one canvas
 */
 
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLine.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>

using namespace std;

static Config const &cfg=Config::get();

int plot(){
	
	io::RootFileSaver saver("danilo/plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/T5Wg_limits_mock/limits_T5Wg_johannes.root","read");
	TFile file_2("../input/T5Wg_limits_mock/limits_T5Wg_knut_2.root","read");
	TFile file_3("../input/T5Wg_limits_mock/limits_T5Wg.root","read");
	
	TGraph *johannes_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *knut_exp = (TGraph*) file_2.Get("exp");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combined_obs = (TGraph*) file_3.Get("gr_obsC_sm");
	
	johannes_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	johannes_exp->SetLineColor(kGreen);
	johannes_exp->SetLineWidth(3);
	johannes_exp->SetLineStyle(2);
	johannes_exp->Draw();
	
	knut_exp->SetLineColor(kBlue);
	knut_exp->SetLineWidth(3);
	knut_exp->SetLineStyle(2);
	knut_exp->Draw("same");
	
	combined_exp->SetLineColor(kRed);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	
	combined_obs->SetLineColor(kBlack);
	combined_obs->SetLineWidth(3);
	combined_obs->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*johannes_exp,"Photon+ST (exp)","l");
	legE.append(*knut_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	legE.append(*combined_obs,"Combined (obs)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_compare_mock",true,false);
	
	return 0;
}

extern "C"

void run(){
	plot();
}