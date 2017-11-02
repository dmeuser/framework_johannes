#include <string>
#include <iostream>
#include <math.h>
#include <algorithm>

#include <TH1.h>
#include <TStyle.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>


using namespace std;

int trigger_13TeV(){ 
	 
   TFile *file = new TFile("output/histograms_v19.root");
 
   TH1F *H_numerator= (TH1F*)file->Get("distributions100.0/effi/trigger_effi_n_match/MET");
   TH1F *H_denominator= (TH1F*)file->Get("distributions100.0/effi/trigger_effi_d_match/MET");

   TH1F *H_numerator_cut= (TH1F*)file->Get("distributions100.0/effi/trigger_effi_n_cut/MET");
   TH1F *H_denominator_cut= (TH1F*)file->Get("distributions100.0/effi/trigger_effi_d_cut/MET");

   const double effi[] = {0,180,1000};
   
 //  H_numerator = (TH1F*)H_numerator->Rebin(2);
 //  H_denominator = (TH1F*)H_denominator->Rebin(2);
      
 //  H_numerator = (TH1F*)H_numerator->Rebin(2,"H_numerator",effi);
 //  H_denominator = (TH1F*)H_denominator->Rebin(2,"H_denominator",effi);

   H_numerator_cut = (TH1F*)H_numerator->Rebin(2,"H_numerator",effi);
   H_denominator_cut = (TH1F*)H_denominator->Rebin(2,"H_denominator",effi);

   double efficiency = H_numerator->GetBinContent(2)/ H_denominator->GetBinContent(2);
      
	TEfficiency *T_effi = new TEfficiency(*H_numerator, *H_denominator);

   //needs Histogram only filled in range want to be used for efficiency calculation
   TGraphAsymmErrors *TAS_effi = new TGraphAsymmErrors();
   TAS_effi->Divide(H_numerator_cut,H_denominator_cut,"v cp");
   double *TAS_Y = TAS_effi->GetY();
   double TAS_Y_effi = TAS_Y[1];
   double TAS_error_high = TAS_effi->GetErrorYhigh(1);
   double TAS_error_low = TAS_effi->GetErrorYlow(1);

   cout << "Efficiency:  " << TAS_Y_effi << "  +  " << TAS_error_high << "  -  " << TAS_error_low << endl;

   

	 gStyle->SetMarkerSize(1);
	 gStyle->SetMarkerStyle(8);

	 
	 H_denominator->SetTitle("");
	 H_denominator->GetYaxis()->SetTitle("#varepsilon_{#gamma-req.}");
	 H_denominator->GetYaxis()->SetTitleOffset(0.5);	
	 H_denominator->GetXaxis()->SetTitleOffset(0.8);
	 H_denominator->GetYaxis()->SetTitleSize(0.6);	
	 H_denominator->GetXaxis()->SetTitleSize(0.15);	 		
	 H_denominator->GetXaxis()->SetTitle("1st gamma_{tight} p_{T} [GeV]");	 	

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetMarkerStyle(8);

	TCanvas *c1 = new TCanvas("c1","Plots",790,600);
	
	TPad *canvasDefault_1 = new TPad("canvasDefault_1", "newpad",0.0,0.0,1.0,0.41);//0.325);		

  canvasDefault_1->SetTopMargin(0.01);
  canvasDefault_1->SetBottomMargin(0.3);
  canvasDefault_1->SetRightMargin(0.05);
  canvasDefault_1->SetLeftMargin(0.11);	
  canvasDefault_1->SetFillStyle(0);
  canvasDefault_1->cd();
  T_effi->SetMarkerSize(1);
  T_effi->SetLineWidth(3);
  T_effi->SetMarkerStyle(8);
 
  T_effi->Draw("");
		
   c1->cd();
	canvasDefault_1->Draw();
	canvasDefault_1->cd();	
	gPad->Update();
	T_effi->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.05);		
	T_effi->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.1);	
   T_effi->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.13);	
   T_effi->GetPaintedGraph()->GetXaxis()->SetLabelSize(0.1);	
   T_effi->GetPaintedGraph()->GetYaxis()->SetLabelSize(0.1);		
   T_effi->GetPaintedGraph()->GetYaxis()->SetTitleOffset(0.38);
	T_effi->SetTitle("");//"#slash{E}_{T}-req. efficiency");	
	T_effi->GetPaintedGraph()->GetYaxis()->SetTitle("#varepsilon_{#gamma trigger}");
	T_effi->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.1);		
	T_effi->GetPaintedGraph()->GetXaxis()->SetTitle("#gamma p_{T} (GeV)");

   std::cout << "Trigger efficiency:   " << efficiency << std::endl;
   gPad->SetGrid();
   gStyle->SetGridStyle(3);
   gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	TH1F *ratio = (TH1F*)H_denominator->Clone("ratio");

	for ( int i = 1; i < H_denominator->GetNbinsX(); i++){
		ratio->SetBinContent(i,0);			
   }
   
  ratio->Draw("AXIS");
  ratio->Draw("same AXIG");	
	ratio->GetYaxis()->SetRangeUser(0,1.05);
	ratio->GetYaxis()->SetNdivisions(10);   
	ratio->GetXaxis()->SetTitleSize(0.125);	
  ratio->GetYaxis()->SetTitleSize(0.14);	
  ratio->GetXaxis()->SetLabelSize(0.1);	
  ratio->GetYaxis()->SetLabelSize(0.1);		
  ratio->GetYaxis()->SetTitleOffset(0.35);
	ratio->SetTitle("");//"#slash{E}_{T}-req. efficiency");	
	ratio->GetYaxis()->SetTitle("#varepsilon_{#gamma trigger}");
	ratio->GetXaxis()->SetTitleOffset(1);		
	ratio->GetXaxis()->SetTitle("#gamma p_{T} (GeV)"); 	
	ratio->GetYaxis()->SetRangeUser(0,1.05);
   T_effi->SetMarkerSize(1);
	T_effi->Draw("same");
	
//	TLatex *EffiPt = new TLatex(0.54,0.45,"#scale[1.32]{#varepsilon}_{#gamma trigger} = 98.4^{+0.7}_{-1.1}(stat.) %");
   TLatex *EffiPt = new TLatex(0.54,0.45,"#scale[1.32]{#varepsilon}_{#gamma trigger} = 97.4 %");
   EffiPt->SetNDC(true);
  EffiPt->SetTextSize(0.1);
 //  EffiPt->Draw();
	TLine *l_cut_pt = new TLine(180,0,180,1.05);
	l_cut_pt->SetLineWidth(4);
	l_cut_pt->Draw("same");
	TLatex *Arrow2 = new TLatex(180,0.6,"#rightarrow");
	Arrow2->SetTextSize(0.13);	
	Arrow2->Draw("same");	

	gPad->SetGrid();
	gStyle->SetGridStyle(3);				
	gPad->Update();		
	
	 c1->cd();	
			
	 TPad *canvasDefault_2 = new TPad("canvasDefault_2", "newpad",0.0,0.4,1,1); //0.32
	 canvasDefault_2->Draw(); 
 	 canvasDefault_2->cd();
 	 canvasDefault_2->SetTopMargin(0.1);
 	 canvasDefault_2->SetBottomMargin(0.01);
 	 canvasDefault_2->SetRightMargin(0.05);
 	 canvasDefault_2->SetLeftMargin(0.11);	 
 	 canvasDefault_2->SetFillStyle(0);

	 H_denominator->SetLineColor(kRed); 
	 H_denominator->SetLineStyle(1);
	 H_denominator->SetLineWidth(4);	 	 	
	 H_denominator->SetMinimum(5);
	 H_numerator->SetLineColor(kBlue);	
	 H_numerator->SetFillStyle(3003);	
	 H_numerator->SetFillColor(kBlue);	  	 
	 H_numerator->SetLineWidth(1);

	 H_denominator->Draw("same hist");
	 H_numerator->Draw("same f hist");		 
	 H_denominator->SetMaximum(25000); 
	 H_denominator->GetYaxis()->SetTitle("Events / bin");
	 H_denominator->GetYaxis()->SetTitleOffset(0.6);	
	 H_denominator->GetXaxis()->SetTitleOffset(0.5);		
	 H_denominator->GetXaxis()->SetTitle("#gamma p_{T} (GeV)"); 	 
	 H_denominator->GetYaxis()->SetTitleSize(0.08); 
	 H_denominator->GetYaxis()->SetLabelSize(0.071);
	 H_denominator->GetXaxis()->SetTitleSize(0.0); 
	 gPad->Update(); 	 
	 gPad->SetGrid();
	gStyle->SetGridStyle(3); 	 	
	 gStyle->SetOptStat(00000);
	 gPad->Update();  	 			 	
    canvasDefault_2->SetLogy(1);

   TLatex *CMS_text = new TLatex(0.11,28000,"#scale[1.3]{#scale[1.2]{CMS}#scale[0.5]{ }#scale[0.8]{#bf{#it{private work}}}#bf{, 35.86 fb^{-1} (13 TeV)}}");
   CMS_text->Draw();	
	TLegend *infoBox = new TLegend(0.23, 0.65, 0.95, 0.9,"");//0.75, 0.83, 0.99, 0.99, "");
  infoBox->AddEntry(H_denominator,"Events accepted by E_{T}^{miss} baseline trigger","l");
  infoBox->AddEntry(H_numerator,"Events accepted by signal and E_{T}^{miss} baseline trigger","f");	
  infoBox->SetShadowColor(0);  // 0 = transparent
  infoBox->SetFillColor(kWhite); 
	infoBox->SetTextSize(0.055);
  infoBox->Draw();	

   return 0;
}


