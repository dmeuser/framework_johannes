/* check pileup reweighting */

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>

Config const &cfg=Config::get();


extern "C"
void run()
{
   std::vector<TString> datasetsToUse={"GJets","ZNuNuGJets","WLNuJets","SinglePhoton"};
   hist::Histograms<TH1F> hs(cfg.datasets.getDatasubsetNames(datasetsToUse));

   hs.addHist("Nvtx"   ,";N_{vtx};normalized events / bin",30,0-.5,30-.5);
   hs.addHist("Nvtx_w" ,";N_{vtx};normalized events / bin",30,0-.5,30-.5);

   for (auto const &dss: cfg.datasets.getDatasubsets(datasetsToUse)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      hs.setCurrentSample(dss.name);

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<Int_t> Nvtx(reader, "nGoodVertices");

      int iEv=0;
      while (reader.Next()) {
         iEv++;
         if (iEv>cfg.processFraction*dss.entries) break;

         float fEventWeight=*w_mc;
         hs.setFillWeight(fEventWeight);
         hs.fill("Nvtx",*Nvtx);
         fEventWeight*=*w_pu;
         hs.setFillWeight(fEventWeight);
         hs.fill("Nvtx_w",*Nvtx);
      } // evt loop

      // hs.scaleLumi();
      hs.mergeOverflow();
      file.Close();
   } // datasets

   hs.combineFromSubsamples({"GJets","ZNuNuGJets","WLNuJets","SinglePhoton"});
   io::RootFileSaver saver("plots.root",TString::Format("pu%.1f",cfg.processFraction*100));
   TCanvas can;
   can.SetLogy();
   for (TString sVar : {"GJets","ZNuNuGJets","WLNuJets"}) {
      gfx::LegendEntries le;
      TH1F h_data(*hs.getHistogram("Nvtx","SinglePhoton"));
      TH1F h_mc  (*hs.getHistogram("Nvtx",sVar));
      TH1F h_mc_w(*hs.getHistogram("Nvtx_w",sVar));

      h_mc.SetLineColor(kGray);
      h_mc.SetFillColor(kGray);
      h_mc.SetFillStyle(1001);
      h_mc_w.SetLineColor(kGray+2);
      h_mc_w.SetFillColor(kGray+2);

      h_data.Scale(1./h_data.Integral());
      h_mc.Scale(1./h_mc.Integral());
      h_mc_w.Scale(1./h_mc_w.Integral());

      h_mc.Draw("hist");
      h_mc_w.Draw("hist same");
      h_data.Draw("pe same");

      le.append(h_data,"data","pe");
      le.append(cfg.datasets.getLabel(sVar));
      le.append(h_mc,"unweighted","l");
      le.append(h_mc_w,"weighted","l");
      TLegend leg=le.buildLegend(.6,.7);
      leg.Draw();

      saver.save(can,sVar);
   }
}
