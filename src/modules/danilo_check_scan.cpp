#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

Config const &cfg=Config::get();

extern "C"
void run()
{
   TFile file("/user/dmeuser/master/root-files/testGGMCombi.root","read");

   TTreeReader reader(cfg.treeName, &file);
   TTreeReaderValue<float> w_pu(reader, "pu_weight");
   TTreeReaderValue<UInt_t> runNo(reader, "runNo");
   TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
   TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
   TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
   TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
   TTreeReaderValue<std::vector<tree::Photon>>   photons  (reader, "photons");
   TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
   TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
   TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
   TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
   TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
   TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
   TTreeReaderValue<tree::MET> MET(reader, "met");
   TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
   TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
   TTreeReaderValue<float> HTgen(reader, "genHt");
   TTreeReaderValue<bool> trigger_Ph   (reader, "HLT_Photon165_HE10_v");
   TTreeReaderValue<bool> baseMETTr(reader, "HLT_PFMET170_HBHECleaned_v");
   TTreeReaderValue<bool> trigger_PhMET(reader, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v");

   int iEv=0;
   //~ int events =0;
   while (reader.Next()){
      iEv++;
      for (tree::IntermediateGenParticle const &intgen: *intermediateGenParticles) {
         if (intgen.pdgId != 1000022) continue;
         std::cout<<"---------------------------"<<std::endl;
         std::cout<<intgen.pdgId<<std::endl;
         for (tree::GenParticle const &gen: intgen.daughters) {
            std::cout<<gen.pdgId<<std::endl;
         }
      }
   }
}
