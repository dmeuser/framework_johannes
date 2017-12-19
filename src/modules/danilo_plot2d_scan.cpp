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


extern "C"
void run()
{
   // First Scan
   TFile file("../input/GGM_scan/hist2d.root","read");
   io::RootFileSaver saver("plots.root","danilo_plot2d_scan");
   TCanvas can;
   
   for (TString svar :{"massNLSP","massChargino","NLO_ms_prospino"}) {
      
      TH2F *hist = (TH2F*) file.Get(svar);
      can.cd();
      if (svar == "NLO_ms_prospino") {
         can.SetLogz();
      }
      else can.SetLogz(false);
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.10);
      hist->GetXaxis()->SetTitle("M_{1} (GeV)");
      hist->GetYaxis()->SetTitle("M_{2} (GeV)");
      if (svar == "massNLSP") {
         hist->GetZaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
      }
      else if (svar == "massChargino") {
         hist->GetZaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
      }
      else if (svar == "NLO_ms_prospino") {
         hist->GetZaxis()->SetTitle("Cross section (fb)");
      }
      hist->GetYaxis()->SetTitleOffset(1.3);
      hist->GetXaxis()->SetTitleOffset(0.9);
      hist->GetZaxis()->SetTitleOffset(1.3);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->GetXaxis()->SetTitleSize(0.05);
      hist->GetZaxis()->SetTitleSize(0.05);
      hist->GetYaxis()->SetLabelSize(0.04);
      hist->GetXaxis()->SetLabelSize(0.04);
      hist->GetZaxis()->SetLabelSize(0.04);
      hist->SetAxisRange(250,1500,"X");
      hist->SetAxisRange(250,1500,"Y");
      hist->SetStats(false);
      hist->Draw("colz");
      can.RedrawAxis();
      saver.save(can,"M1_M2/"+svar,true,true);
      can.Clear();
      
   }
   
   
   //Second Scan
   TFile file2("../input/GGM_scan/hist2d_scan2.root","read");
   
   for (TString svar :{"massNLSP","massChargino","NLO_ms_prospino"}) {
      
      TH2F *hist = (TH2F*) file2.Get(svar);
      can.cd();
      if (svar == "NLO_ms_prospino") {
         can.SetLogz();
      }
      else can.SetLogz(false);
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.10);
      hist->GetXaxis()->SetTitle("M_{1} (GeV)");
      hist->GetYaxis()->SetTitle("M_{3} (GeV)");
      if (svar == "massNLSP") {
         hist->GetZaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
      }
      else if (svar == "massChargino") {
         hist->GetZaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
      }
      else if (svar == "NLO_ms_prospino") {
         hist->GetZaxis()->SetTitle("Cross section (fb)");
      }
      hist->GetYaxis()->SetTitleOffset(1.3);
      hist->GetXaxis()->SetTitleOffset(0.9);
      hist->GetZaxis()->SetTitleOffset(1.3);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->GetXaxis()->SetTitleSize(0.05);
      hist->GetZaxis()->SetTitleSize(0.05);
      hist->GetYaxis()->SetLabelSize(0.04);
      hist->GetXaxis()->SetLabelSize(0.04);
      hist->GetZaxis()->SetLabelSize(0.04);
      //~ hist->SetAxisRange(250,1500,"X");
      //~ hist->SetAxisRange(250,1500,"Y");
      hist->SetStats(false);
      hist->Draw("colz");
      can.RedrawAxis();
      saver.save(can,"M1_M3/"+svar,true,true);
      can.Clear();
      
   }
}
