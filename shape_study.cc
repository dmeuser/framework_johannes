#include <string>
#include <iostream>
#include <math.h>
#include <algorithm>

#include <TH1.h>
#include <TGraph.h>
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

int shape_study(){

   TCanvas can1, can2, can3, can4, can5, can6;
   can1.cd();
   TGraph *g_pdf_Vg = new TGraph("/user/jschulz/2016/photonmet/output/pdf_SR_Vg_graph.tex");
   gStyle->SetPadBottomMargin(0.3);
   gStyle->SetPadLeftMargin(0.22);
   g_pdf_Vg->Draw("AP*");
   g_pdf_Vg->SetTitle("PDF uncertainty in shape for V#gamma");
   g_pdf_Vg->GetXaxis()->SetTitle("S_{T}^{#gamma} (GeV)");
   g_pdf_Vg->GetYaxis()->SetTitle("prediction variation in %");
   g_pdf_Vg->GetXaxis()->SetLabelSize(1.5*g_pdf_Vg->GetXaxis()->GetLabelSize());
   g_pdf_Vg->GetXaxis()->SetTitleOffset(0.85);
   
   g_pdf_Vg->GetYaxis()->SetLabelSize(1.5*g_pdf_Vg->GetYaxis()->GetLabelSize());   
   g_pdf_Vg->GetXaxis()->SetTitleSize(1.5*g_pdf_Vg->GetXaxis()->GetTitleSize());
   g_pdf_Vg->GetYaxis()->SetTitleSize(1.5*g_pdf_Vg->GetYaxis()->GetTitleSize());
   TGraph *g_pdf_Vg_mean = new TGraph();
   g_pdf_Vg_mean->SetPoint(0,700,0.0159707);
   g_pdf_Vg_mean->SetPoint(1,900,0.0200078);
   g_pdf_Vg_mean->SetPoint(2,1150,0.0244478);
   g_pdf_Vg_mean->SetPoint(3,1450,0.0382614);
   g_pdf_Vg_mean->SetPoint(4,700,-0.0159707);
   g_pdf_Vg_mean->SetPoint(5,900,-0.0200078);
   g_pdf_Vg_mean->SetPoint(6,1150,-0.0244478);
   g_pdf_Vg_mean->SetPoint(7,1450,-0.0382614);
   g_pdf_Vg_mean->SetFillStyle(3002);
   g_pdf_Vg_mean->SetMarkerStyle(25);  
   g_pdf_Vg_mean->SetMarkerColor(kRed);
   g_pdf_Vg_mean->SetFillColor(kRed);
   g_pdf_Vg_mean->SetMarkerSize(2);
   g_pdf_Vg_mean->Draw("B same");
   can1.SaveAs("Vg_pdf_shape.pdf");

   can2.cd();
   TGraph *g_pdf_gJ = new TGraph("/user/jschulz/2016/photonmet/output/pdf_SR_gJ_graph.tex");
   gStyle->SetPadBottomMargin(0.3);
   gStyle->SetPadLeftMargin(0.22);
   g_pdf_gJ->Draw("AP*");
   g_pdf_gJ->SetTitle("PDF uncertainty in shape for #gamma+jets");
   g_pdf_gJ->GetXaxis()->SetTitle("S_{T}^{#gamma} (GeV)");
   g_pdf_gJ->GetYaxis()->SetTitle("prediction variation in %");
   g_pdf_gJ->GetXaxis()->SetLabelSize(1.5*g_pdf_gJ->GetXaxis()->GetLabelSize());
   g_pdf_gJ->GetXaxis()->SetTitleOffset(0.85);
   
   g_pdf_gJ->GetYaxis()->SetLabelSize(1.5*g_pdf_gJ->GetYaxis()->GetLabelSize());   
   g_pdf_gJ->GetXaxis()->SetTitleSize(1.5*g_pdf_gJ->GetXaxis()->GetTitleSize());
   g_pdf_gJ->GetYaxis()->SetTitleSize(1.5*g_pdf_gJ->GetYaxis()->GetTitleSize());
   TGraph *g_pdf_gJ_mean = new TGraph();
   g_pdf_gJ_mean->SetPoint(0,700,0.0815288);
   g_pdf_gJ_mean->SetPoint(1,900,0.0255214);
   g_pdf_gJ_mean->SetPoint(2,1150,0.0192825);
   g_pdf_gJ_mean->SetPoint(3,1450,0.0512977);
   g_pdf_gJ_mean->SetPoint(4,700,-0.0815288);
   g_pdf_gJ_mean->SetPoint(5,900,-0.0255214);
   g_pdf_gJ_mean->SetPoint(6,1150,-0.0192825);
   g_pdf_gJ_mean->SetPoint(7,1450,-0.0512977);
   g_pdf_gJ_mean->SetFillStyle(3002);
   g_pdf_gJ_mean->SetMarkerStyle(25);  
   g_pdf_gJ_mean->SetMarkerColor(kRed);
   g_pdf_gJ_mean->SetFillColor(kRed);
   g_pdf_gJ_mean->SetMarkerSize(2);
   g_pdf_gJ_mean->Draw("B same");
   can2.SaveAs("gJ_pdf_shape.pdf");

   can3.cd();
   TGraph *g_mu_Vg = new TGraph("/user/jschulz/2016/photonmet/output/mu_SR_Vg_graph.tex");
   gStyle->SetPadBottomMargin(0.3);
   gStyle->SetPadLeftMargin(0.22);
   g_mu_Vg->Draw("AP*");
   g_mu_Vg->SetTitle("#mu_{R} and #mu_{F} uncertainty in shape for V#gamma");
   g_mu_Vg->GetXaxis()->SetTitle("S_{T}^{#gamma} (GeV)");
   g_mu_Vg->GetYaxis()->SetTitle("prediction variation in %");
   g_mu_Vg->GetXaxis()->SetLabelSize(1.5*g_mu_Vg->GetXaxis()->GetLabelSize());
   g_mu_Vg->GetXaxis()->SetTitleOffset(0.85);
   
   g_mu_Vg->GetYaxis()->SetLabelSize(1.5*g_mu_Vg->GetYaxis()->GetLabelSize());   
   g_mu_Vg->GetXaxis()->SetTitleSize(1.5*g_mu_Vg->GetXaxis()->GetTitleSize());
   g_mu_Vg->GetYaxis()->SetTitleSize(1.5*g_mu_Vg->GetYaxis()->GetTitleSize());
   can3.SaveAs("Vg_mu_shape.pdf");

   can4.cd();
   TGraph *g_mu_gJ = new TGraph("/user/jschulz/2016/photonmet/output/mu_SR_gJ_graph.tex");
   g_mu_gJ->Draw("AP*");
   g_mu_gJ->SetTitle("#mu_{R} and #mu_{F} uncertainty in shape for #gamma+jets");
   g_mu_gJ->GetXaxis()->SetTitle("S_{T}^{#gamma} (GeV)");
   g_mu_gJ->GetYaxis()->SetTitle("prediction variation in %");
   g_mu_gJ->GetXaxis()->SetLabelSize(1.5*g_mu_gJ->GetXaxis()->GetLabelSize());
   g_mu_gJ->GetXaxis()->SetTitleOffset(0.85);
   
   g_mu_gJ->GetYaxis()->SetLabelSize(1.5*g_mu_gJ->GetYaxis()->GetLabelSize());   
   g_mu_gJ->GetXaxis()->SetTitleSize(1.5*g_mu_gJ->GetXaxis()->GetTitleSize());
   g_mu_gJ->GetYaxis()->SetTitleSize(1.5*g_mu_gJ->GetYaxis()->GetTitleSize());
   can4.SaveAs("gJ_mu_shape.pdf");

   can5.cd();
   TGraph *g_JES_Vg = new TGraph("/user/jschulz/2016/photonmet/output/JES_SR_Vg_graph.tex");
   g_JES_Vg->SetTitle("Vgamma JES study in photon pT in CR");
   g_JES_Vg->Draw("A*");
   can5.SaveAs("Vg_JES_shape.pdf");

   can6.cd();
   TGraph *g_JES_gJ = new TGraph("/user/jschulz/2016/photonmet/output/JES_SR_gJ_graph.tex");
   g_JES_gJ->SetTitle("gamma+jets JES study in photon pT in CR");
   g_JES_gJ->Draw("A*");
   can6.SaveAs("gJ_JES_shape.pdf");
   
      
   return 0;
}


