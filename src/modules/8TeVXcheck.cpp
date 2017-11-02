#include "Config.hpp"
#include "tools/io.hpp"
#include "tools/gfx.hpp"

#include <TFile.h>
#include <TTree.h>

extern "C"
void run(){
   Config cfg=Config::get();

   TFile fNew(cfg.dataBasePath+"../MINIAODSIM/ggm/M2_640_M1_630_MINIAODSIM8TeV.root","read");
   TFile fOld("/user/lange/data/run2/JS/JS_8TeV_WB_640_630_nTuple.root");
   if (fOld.IsZombie()||fNew.IsZombie()) {
      return;
   }
   TTree *tNew=(TTree*)fNew.Get("Events");
   TTree *tOld=(TTree*)fOld.Get("susyTree");

   TCanvas can;
   // can.SetLogy();
   io::RootFileSaver saver("plots.root","8TeVXcheck",true);
   for (int pdgId:{22,1000039}){
      TString label=pdgId==22?"#gamma":"gravitino";
      TH1F hNew("hNew",";gen "+label+" PT;normalized entries",50,0,1000);
      TH1F hOld("hOld",";gen "+label+" PT;normalized entries",50,0,1000);
      hNew.Sumw2(true);
      hOld.Sumw2(true);
      TString pdgIdCut=TString::Format("recoGenParticles_prunedGenParticles__PAT.obj.pdgId()==%d",pdgId);
      if (pdgId==22) pdgIdCut+=" && recoGenParticles_prunedGenParticles__PAT.obj.pt()>80";
      tNew->Draw("recoGenParticles_prunedGenParticles__PAT.obj.pt()>>hNew",
                 pdgIdCut,"goff");
      pdgIdCut=TString::Format("genParticles.pdgId==%d",pdgId);
      if (pdgId==22) pdgIdCut+=" && genParticles.momentum.Pt()>80";
      tOld->Draw("genParticles.momentum.Pt()>>hOld",pdgIdCut,"goff");

      hNew.Scale(1./hNew.Integral());
      hOld.Scale(1./hOld.Integral());
      hNew.Draw("e");
      hOld.Draw("e same");
      hNew.SetMarkerSize(0);
      hOld.SetMarkerSize(0);
      hNew.SetLineColor(kBlack);
      hOld.SetLineColor(kGray+1);
      TLegend leg=gfx::mkLegend(.6,.7);
      leg.SetHeader("8 TeV GGM");
      leg.AddEntry(&hOld, "JS (Pythia6)", "le");
      leg.AddEntry(&hNew, "JL (Pythia8)", "le");
      leg.Draw();

      float chi2=hOld.Chi2Test(&hNew,"WW CHI2/NDF P");
      hOld.KolmogorovTest(&hNew,"X D");
      float ks=hOld.KolmogorovTest(&hNew);
      TPaveText lab=gfx::paveText(.6,.35,.9,.5);
      lab.AddText(TString::Format("#chi^{2}/N_{df}=%.1f",chi2));
      lab.AddText(TString::Format("KS prob=%.2f",ks));
      lab.Draw();

      label=pdgId==22?"gamma":"gtino";
      saver.save(can, "GGM_"+label+"_pt",true,false);
   }
   fNew.Close();
   fOld.Close();
}
