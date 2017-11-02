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
   g_pdf_Vg->Draw("A*");
   TGraph *g_pdf_Vg_mean = new TGraph();
   g_pdf_Vg_mean->SetPoint(0,700,0.0165178);
   g_pdf_Vg_mean->SetPoint(1,900,0.0292389);
   g_pdf_Vg_mean->SetPoint(2,1150,0.0339644);
   g_pdf_Vg_mean->SetPoint(3,1450,0.0530947);
   g_pdf_Vg_mean->SetPoint(4,700,-0.0165178);
   g_pdf_Vg_mean->SetPoint(5,900,-0.0292389);
   g_pdf_Vg_mean->SetPoint(6,1150,-0.0339644);
   g_pdf_Vg_mean->SetPoint(7,1450,-0.0530947);
   g_pdf_Vg_mean->SetFillStyle(3002);
   g_pdf_Vg_mean->SetMarkerStyle(25);  
   g_pdf_Vg_mean->SetMarkerColor(kRed);
   g_pdf_Vg_mean->SetFillColor(kRed);
   g_pdf_Vg_mean->SetMarkerSize(2);
   g_pdf_Vg_mean->Draw("B same");
   can1.SaveAs("Vg_pdf_shape.pdf");

   can2.cd();
   TGraph *g_pdf_gJ = new TGraph("/user/jschulz/2016/photonmet/output/pdf_SR_gJ_graph.tex");
   g_pdf_gJ->Draw("A*");
   TGraph *g_pdf_gJ_mean = new TGraph();
   g_pdf_gJ_mean->SetPoint(0,700,0.0385493);
   g_pdf_gJ_mean->SetPoint(1,900,0.0175165);
   g_pdf_gJ_mean->SetPoint(2,1150,0.0707174);
   g_pdf_gJ_mean->SetPoint(3,1450,0.0648986);
   g_pdf_gJ_mean->SetPoint(4,700,-0.0385493);
   g_pdf_gJ_mean->SetPoint(5,900,-0.0175165);
   g_pdf_gJ_mean->SetPoint(6,1150,-0.0707174);
   g_pdf_gJ_mean->SetPoint(7,1450,-0.0648986);
   g_pdf_gJ_mean->SetFillStyle(3002);
   g_pdf_gJ_mean->SetMarkerStyle(25);  
   g_pdf_gJ_mean->SetMarkerColor(kRed);
   g_pdf_gJ_mean->SetFillColor(kRed);
   g_pdf_gJ_mean->SetMarkerSize(2);
   g_pdf_gJ_mean->Draw("B same");
   can2.SaveAs("gJ_pdf_shape.pdf");

   can3.cd();
   TGraph *g_mu_Vg = new TGraph("/user/jschulz/2016/photonmet/output/mu_SR_Vg_graph.tex");
   g_mu_Vg->Draw("A*");
   can3.SaveAs("Vg_mu_shape.pdf");

   can4.cd();
   TGraph *g_mu_gJ = new TGraph("/user/jschulz/2016/photonmet/output/mu_SR_gJ_graph.tex");
   g_mu_gJ->Draw("A*");
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


