#ifndef TESTS_HPP__
#define TESTS_HPP__

#include "Config.hpp"
#include "tools/io.hpp"
#include "tools/physics.hpp"
#include "tools/gfx.hpp"
#include "tools/hist.hpp"

#include <TCanvas.h>
#include <TChain.h>
#include <TH1F.h>
#include <TTreeReader.h>

extern "C"
void run()
{
   Config const &cfg=Config::get();

   io::log<<"== Logger tests:";
   io::log*"Start"/"Space"*"Nospace"<<"Endnopace";
   io::log*"Start"/"Space"*"Nospace">>"Endspace";

   io::log<<""<<"== Datasets tests:";
   std::vector<TString> datasets;
   for (Datasubset const &dss:cfg.datasets.getDatasubsets()){
      datasets.push_back(dss.filename);
      io::log << TString::Format("  %s %i %g",dss.name.c_str(),dss.Ngen,dss.xsec);
   }
   io::log<<"Signal:";
   for (auto dss: cfg.datasets.getDataset("T5gg").subsets){
      io::log << TString::Format("  %s %i %g",dss.name.c_str(),dss.Ngen,dss.xsec);
   }

   TChain chain(Config::get().treeName);
   for(TString dataset: datasets){
      // std::cout << dataset << std::endl;
      chain.AddFile(Config::get().dataBasePath+dataset);
   }
   io::log*"\ntotal MC events "<<chain.GetEntries();

   io::log<<"Data:";
   for (auto dss: cfg.datasets.getDataset("SinglePhoton").subsets){
      io::log*"  "*dss.filename;
   }

   io::log<<""<<""<<"== Config tests:";
   io::log*"Lumi"/cfg.lumi*", process fraction"/cfg.processFraction
      *", release mode">>cfg.releaseMode;
   io::log*"Modules to run:";
   for (std::string m: cfg.modules) io::log/m;
   io::log<<"";

   io::log<<""<<"== Plotting/Decoration tests:";
   io::RootFileSaver saver("plots.root","tests");
   TCanvas can;
   TH1F h("",";x-title PT;EntriesBIN",100,-4,4);
   h.FillRandom("gaus");
   h.Draw();
   saver.save(can,"titleTestUnit");
   h.GetXaxis()->SetTitle("x-title");
   h.GetYaxis()->SetTitle("EntriesBIN");
   TLatex l1=gfx::cornerLabel("pos1",1);
   TLatex l2=gfx::cornerLabel("pos2",2);
   TLatex l3=gfx::cornerLabel("pos3",3);
   TLatex l4=gfx::cornerLabel("pos4",4);
   l1.DrawClone();
   l2.DrawClone();
   l3.DrawClone();
   l4.DrawClone();

   saver.save(can,"titleTestNoUnit");

   io::RootFileSaver saverNoSim("plots.root","tests",false);
   TCanvas canNoSim;
   h.Draw();
   saverNoSim.save(canNoSim,"noSim");

   can.cd();
   can.SetLogy();
   TH1F hVary=hist::fromWidths("hVary",";;Entries / bin",
                               {-5,-4,-2,0  ,3,10},
                               { 1,.5,.1,1.5,7   });
   hVary.FillRandom("gaus",1e6);
   hVary.Draw();
   saver.save(can,"variableBins");
}

#endif /* TESTS_HPP__ */
