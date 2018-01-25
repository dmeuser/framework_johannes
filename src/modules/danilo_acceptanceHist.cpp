//Transforming the Tgraph2d into a TH2D for the acceptance
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


extern "C"
void run()
{     
      for (TString scan :{"GGM_M1_M2","GGM_M1_M3","T5Wg"}) {
            
            for (TString selection :{"exclusiv","inclusiv"}) {
                  
                  TFile file("../output/signal_scan_"+selection+"_v03D.root","read");
                  
                  io::RootFileSaver saver("plots.root","danilo_acceptanceHist");
                  TCanvas can;
                  
                  TString path = scan+"/pre_ph165/c_MET300/MT300/STg/"+scan+"_acceptance";
                  TGraph2D *graph = (TGraph2D*) file.Get(path);
                  can.cd();
                  TH2D *hist = graph->GetHistogram();
                  gPad->SetRightMargin(0.2);
                  gPad->SetLeftMargin(0.13);
                  gPad->SetBottomMargin(0.10);
                  hist->SetBit(TH1::kNoTitle);
                  hist->GetXaxis()->SetTitle("M_{1} (GeV)");
                  if (scan=="GGM_M1_M2") hist->GetYaxis()->SetTitle("M_{2} (GeV)");
                  else if (scan=="GGM_M1_M3") hist->GetYaxis()->SetTitle("M_{3} (GeV)");
                  else if (scan=="T5Wg") {
                        hist->GetXaxis()->SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
                        hist->GetYaxis()->SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
                  }
                  hist->GetZaxis()->SetTitle("acceptance [%]");
                  hist->GetYaxis()->SetTitleOffset(1.3);
                  hist->GetXaxis()->SetTitleOffset(0.9);
                  hist->GetZaxis()->SetTitleOffset(1.3);
                  hist->GetYaxis()->SetTitleSize(0.05);
                  hist->GetXaxis()->SetTitleSize(0.05);
                  hist->GetZaxis()->SetTitleSize(0.05);
                  hist->GetYaxis()->SetLabelSize(0.04);
                  hist->GetXaxis()->SetLabelSize(0.04);
                  hist->GetZaxis()->SetLabelSize(0.04);
                  hist->GetZaxis()->SetLabelOffset(0.02);
                  hist->SetStats(false);
                  hist->SetMaximum(0.4);
                  hist->Draw("colz");
                  can.RedrawAxis();
                  saver.save(can,scan+"_"+selection+"_acceptance",true,true);
                  can.Clear();
            }
      }
}
