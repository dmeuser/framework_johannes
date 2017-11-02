/* check statistical precision of JES shift in tail*/

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("JES_studies%.1f",cfg.processFraction*100));

   TH1F hnorm(*histReader.read<TH1F>("pre_ph165/c_MET300/MT300/STg600/STg/GJets_DR"));
   hnorm=hist::rebinned(hnorm,{600,800,1000,1300,1600},{200,200,300,300});
   TH1F hUp(*histReader.read<TH1F>("pre_ph165/c_MET300/MT300/STg600/STg_JESu/GJets_DR"));
   hUp=hist::rebinned(hUp,{600,800,1000,1300,1600},{200,200,300,300});
   TH1F hDown(*histReader.read<TH1F>("pre_ph165/c_MET300/MT300/STg600/STg_JESd/GJets_DR"));
   hDown=hist::rebinned(hDown,{600,800,1000,1300,1600},{200,200,300,300});

   hnorm.Scale(cfg.trigger_eff_Ph*1.714741);
   hUp.Scale(cfg.trigger_eff_Ph*1.723989);
   hDown.Scale(cfg.trigger_eff_Ph*1.794585);

   debug << hUp.GetBinContent(4) << hDown.GetBinContent(4);

   io::RootFileSaver saver("plots.root","MET_shift");
   gfx::SplitCan spcan;
   spcan.cdUp();
   gPad->SetLogy();
   hnorm.SetMarkerSize(0);
   hnorm.SetLineColor(kBlack);
   hUp.SetMarkerSize(0);
   hUp.SetLineColor(kRed);
   hDown.SetLineColor(kBlue);
   hDown.SetMarkerSize(0);  
   hnorm.Draw("hist e");
   hUp.Draw("hist e same");
   hDown.Draw("hist e same");
   
   gfx::LegendEntries le;
   le.append(hnorm,"GJets norm","l");
   le.append(hUp,"GJets MET up","l");
   le.append(hDown,"GJets MET down","l");
   le.buildLegend(.6,.75).DrawClone();

   spcan.cdLow();
   TH1F histRatio=hist::getRatio(hUp,hnorm,"shift/norm",hist::COMB);
   TH1F histRatioDown=hist::getRatio(hDown,hnorm,"shift/norm",hist::COMB);
   histRatio.SetMarkerStyle(kOpenCircle);
   histRatioDown.SetMarkerStyle(kOpenCircle);
   histRatio.SetMarkerSize(1);
   histRatio.SetLineColor(kRed);
   histRatioDown.SetLineColor(kBlue);  
   histRatio.Draw("pe");
   histRatioDown.Draw("p e same");
   histRatio.SetMinimum(0.1);
   histRatio.SetMaximum(2.9);
   saver.save(spcan,"sr_STg");
}
