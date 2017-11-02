#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"

#include <TChain.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TTreeReader.h>

static Config const &cfg=Config::get();

static void collectionMatch(TString const &datasetName,bool sim)
{
   Dataset const &dataset=cfg.datasets.getDataset(datasetName);
   TChain chain(cfg.treeName);
   for (auto const &dss: dataset.subsets){
      chain.Add(dss.getPath());
   }

   TTreeReader reader(&chain);
   TTreeReaderValue<std::vector<tree::Photon>>   photons  (reader, "photons");
   TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
   TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
   TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets"   );

   TH2F hPhJ      ("",";|p_{T}^{#gamma} #minus p_{T}^{jet}| / p_{T}^{#gamma};#DeltaR(#gamma, jet);Entries / bin",100,0,10,100,0,8);
   TH2F hPhJ_zoom ("",";|p_{T}^{#gamma} #minus p_{T}^{jet}| / p_{T}^{#gamma};#DeltaR(#gamma, jet);Entries / bin",100,0,2,100,0,.6);
   TH2F hElJ      ("",";|p_{T}^{e} #minus p_{T}^{jet}| / p_{T}^{e};#DeltaR(e, jet);Entries / bin",100,0,10,100,0,8);
   TH2F hElJ_zoom ("",";|p_{T}^{e} #minus p_{T}^{jet}| / p_{T}^{e};#DeltaR(e, jet);Entries / bin",100,0,2,100,0,.6);
   TH2F hMuJ      ("",";|p_{T}^{#mu} #minus p_{T}^{jet}| / p_{T}^{#mu};#DeltaR(#mu, jet);Entries / bin",100,0,10,100,0,8);
   TH2F hMuJ_zoom ("",";|p_{T}^{#mu} #minus p_{T}^{jet}| / p_{T}^{#mu};#DeltaR(#mu, jet);Entries / bin",100,0,10,100,0,.6);
   TH2F hPhEl     ("",";|p_{T}^{#gamma} #minus p_{T}^{e}| / p_{T}^{#gamma};#DeltaR(#gamma, e);Entries / bin",100,0,10,100,0,8);
   TH2F hPhEl_zoom("",";|p_{T}^{#gamma} #minus p_{T}^{e}| / p_{T}^{#gamma};#DeltaR(#gamma, e);Entries / bin",100,0,2,100,0,.6);

   float dR,dPt;
   while (reader.Next()){
      for (auto const &ph:*photons){
         // if (!ph.isLoose) continue; // only loose anyway
         for (auto const &jet:*jets){
            if (!jet.isLoose) continue;
            dR=ph.p.DeltaR(jet.p);
            dPt=fabs(ph.p.Pt()-jet.p.Pt())/ph.p.Pt();
            hPhJ.Fill(dPt, dR);
            hPhJ_zoom.Fill(dPt, dR);
         }
         for (auto const &el:*electrons){
            if (!el.isLoose) continue;
            dR=ph.p.DeltaR(el.p);
            dPt=fabs(ph.p.Pt()-el.p.Pt())/ph.p.Pt();
            hPhEl.Fill(dPt, dR);
            hPhEl_zoom.Fill(dPt, dR);
         }
      }
      for (auto const &jet:*jets){
         if (!jet.isLoose) continue;
         for (auto const &el:*electrons){
            if (!el.isLoose) continue;
            dR=el.p.DeltaR(jet.p);
            dPt=fabs(el.p.Pt()-jet.p.Pt())/el.p.Pt();
            hElJ.Fill(dPt, dR);
            hElJ_zoom.Fill(dPt, dR);
         }
         for (auto const &mu:*muons){
            dR=mu.p.DeltaR(jet.p);
            dPt=fabs(mu.p.Pt()-jet.p.Pt())/mu.p.Pt();
            hMuJ.Fill(dPt, dR);
            hMuJ_zoom.Fill(dPt, dR);
         }
      }
   } // event loop
   io::RootFileSaver saver("plots.root","matching");
   TCanvas can;
   can.SetRightMargin (.15);
   can.SetBottomMargin(.22);
   can.SetLogz();

   hPhJ.Draw("colz");
   gfx::setupZaxis(hPhJ);
   saver.save(can,"photon-jet_"+datasetName,true,sim);
   hPhJ_zoom.Draw("colz");
   gfx::setupZaxis(hPhJ_zoom);
   saver.save(can,"photon-jet(zoom)_"+datasetName,true,sim);

   hPhEl.Draw("colz");
   gfx::setupZaxis(hPhEl);
   saver.save(can,"photon-electron_"+datasetName,true,sim);
   hPhEl_zoom.Draw("colz");
   gfx::setupZaxis(hPhEl_zoom);
   saver.save(can,"photon-electron(zoom)_"+datasetName,true,sim);

   hElJ.Draw("colz");
   gfx::setupZaxis(hElJ);
   saver.save(can,"electron-jet_"+datasetName,true,sim);
   hElJ_zoom.Draw("colz");
   gfx::setupZaxis(hElJ_zoom);
   saver.save(can,"electron-jet(zoom)_"+datasetName,true,sim);

   hMuJ.Draw("colz");
   gfx::setupZaxis(hMuJ);
   saver.save(can,"muon-jet_"+datasetName,true,sim);
   hMuJ_zoom.Draw("colz");
   gfx::setupZaxis(hMuJ_zoom);
   saver.save(can,"muon-jet(zoom)_"+datasetName,true,sim);
}

static void genMatch(TString const &datasetName,bool sim)
{
   Dataset const &dataset=cfg.datasets.getDataset(datasetName);
   TChain chain(cfg.treeName);
   for (auto const &dss: dataset.subsets){
      chain.Add(dss.getPath());
   }

   TTreeReader reader(&chain);
   TTreeReaderValue<std::vector<tree::Photon>>     photons  (reader, "photons");
   TTreeReaderValue<std::vector<tree::GenParticle>>genPart  (reader, "genParticles");

   TH2F hPh      ("",";|p_{T}^{reco} #minus p_{T}^{gen}| / p_{T}^{reco};#DeltaR(reco, gen);Entries / bin",100,0,10,100,0,8);
   TH2F hPh_zoom ("",";|p_{T}^{reco} #minus p_{T}^{gen}| / p_{T}^{reco};#DeltaR(reco, gen);Entries / bin",100,0,2,100,0,.6);

   float dR,dPt;
   while (reader.Next()){
      for (auto const &ph:*photons){
         // if (!ph.isLoose) continue; // only loose anyway
         for (auto const &genp:*genPart){
            if (genp.pdgId!=22) continue; // no photon
            dR=ph.p.DeltaR(genp.p);
            dPt=fabs(ph.p.Pt()-genp.p.Pt())/ph.p.Pt();
            hPh.Fill(dPt, dR);
            hPh_zoom.Fill(dPt, dR);
         }
      }

   } // event loop
   io::RootFileSaver saver("plots.root","matching");
   TCanvas can;
   can.SetRightMargin (.15);
   can.SetBottomMargin(.22);
   can.SetLogz();

   hPh.Draw("colz");
   gfx::setupZaxis(hPh);
   saver.save(can,"photon-gen_"+datasetName,true,sim);
   hPh_zoom.Draw("colz");
   gfx::setupZaxis(hPh_zoom);
   saver.save(can,"photon-gen(zoom)_"+datasetName,true,sim);
}


extern "C"
void run()
{
   collectionMatch("SinglePhoton",false);
   collectionMatch("GJets",true);
   genMatch("WGToLNuG",true);
   genMatch("GJets",true);
}
