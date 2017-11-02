/* check ratio of W,Z */

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   TH1F hW(*histReader.read<TH1F>("pre_ph165/c_S80/MT300/STg/WGToLNuG"));
   hW=hist::rebinned(hW,{260,600,800,1000,1200},{340,200,200,200});
   TH1F hZ(*histReader.read<TH1F>("pre_ph165/c_S80/MT300/STg/ZNuNuGJets"));
   hZ=hist::rebinned(hZ,{260,600,800,1000,1200},{340,200,200,200});

   hW.Scale(cfg.trigger_eff_Ph*cfg.sf.Vg);
   hZ.Scale(cfg.trigger_eff_Ph*cfg.sf.Vg);

   io::RootFileSaver saver("plots.root","wz_ratio");
   gfx::SplitCan spcan;
   spcan.cdUp();
   gPad->SetLogy();
   hZ.SetMarkerSize(0);
   hW.SetMarkerSize(0);
   hZ.Draw("hist e");
   hW.Draw("hist e same");

   gfx::LegendEntries le;
   le.append(hW,cfg.datasets.getLabel("WGToLNuG"),"l");
   le.append(hZ,cfg.datasets.getLabel("ZNuNuGJets"),"l");
   le.buildLegend(.6,.75).DrawClone();

   spcan.cdLow();
   TH1F histRatio=hist::getRatio(hW,hZ,"W/Z",hist::COMB);
   histRatio.SetMarkerStyle(kOpenCircle);
   histRatio.SetMarkerSize(1);
   histRatio.Draw("pe");
   histRatio.SetMinimum(0.75);
   histRatio.SetMaximum(1.25);
   saver.save(spcan,"sr_STg");
}
