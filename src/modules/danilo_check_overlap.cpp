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
   //~ TFile file("/user/dmeuser/master/data/v01D/SinglePhoton_03Feb2017.root","read");
   TFile file("/user/dmeuser/master/root-files/overlap_lepton_1.root","read");

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
   
   std::vector<long int> wrong = {1289901438,89000916,321243519};
   //~ std::vector<long int> wrong = {844732197,321243519};

   while (reader.Next()){
      if (std::find(wrong.begin(), wrong.end(), *evtNo) != wrong.end()){
         std::cout << *runNo << ":" << *lumNo << ":" << *evtNo <<std::endl;
         //~ std::cout << (*MET).p.Pt() << std::endl;
         //~ 
         //~ std::cout << phys::invmass((*electrons)[0],(*photons)[0]) << std::endl;
         //~ std::cout << phys::M_T((*electrons)[0],*MET) << std::endl;
         
         
         //~ for (tree::Muon const &mu: *muons) {
            //~ std::cout << "Muon" <<std::endl;
            //~ std::cout << mu.p.DeltaR((*photons)[0].p) << std::endl;
         //~ }
         //~ 
         //~ for (tree::Electron const &ele: *electrons) {
            //~ std::cout << "Electron" <<std::endl;
            //~ std::cout << ele.p.DeltaR((*photons)[0].p) << std::endl;
         //~ }
         
         //~ for (tree::Photon const &pho: *photons){
            //~ for (tree::Photon const &pho_2: *photons){
               //~ std::cout << fabs(pho.p.Eta()-pho_2.p.Eta()) << "," << fabs(pho.p.Phi()-pho_2.p.Phi()) << std::endl;
            //~ }
         //~ }
         std::cout << phys::invmass((*photons)[0],(*photons)[1]) << std::endl;
         //~ std::cout << phys::invmass((*photons)[0],(*photons)[1]) << std::endl;
         //~ std::cout << (*photons)[0].p.DeltaR((*photons)[1].p) << std::endl;
         
         //~ bool leptoVeto = false;
         //~ for (tree::Muon const &mu: *muons) {
            //~ if (mu.p.Pt() < 25 || mu.isMedium == 0) { continue; }
            //~ if (fabs(mu.p.Eta()) < 2.4 && mu.PFminiIso < 0.2 &&
                  //~ fabs(mu.d0) < 0.05 && fabs(mu.dZ) < 0.1) {
                     //~ leptoVeto = true;
            //~ }
         //~ }
         
         //~ for (tree::Electron const &ele: *electrons) {
            //~ if (ele.p.Pt() < 25 || ele.isMedium == 0) { continue; }
            //~ if (fabs(ele.p.Eta()) < 1.4442) {
               //~ if (ele.r9 > 0.5 && ele.SigmaIEtaIEtaFull5x5 < 0.00998 &&
                     //~ fabs(ele.dEtaAtVtx) < 0.00311 && fabs(ele.dPhiAtVtx) < 0.103 &&
                     //~ ele.HoverE < 0.253 && ele.MissHits <= 1 && fabs(ele.EoverPInv) < 0.134 &&
                     //~ ele.ConvVeto == 1 && ele.PFminiIso < 0.1) {
                        //~ leptoVeto = true;
               //~ }
            //~ }
            //~ if (fabs(ele.p.Eta()) > 1.56 && fabs(ele.p.Eta()) < 2.5) {
               //~ if (ele.r9 > 0.8 && ele.SigmaIEtaIEtaFull5x5 < 0.0298 &&
                     //~ fabs(ele.dEtaAtVtx) < 0.00609 && fabs(ele.dPhiAtVtx) < 0.045 &&
                     //~ ele.HoverE < 0.0878 && ele.MissHits <= 1 && fabs(ele.EoverPInv) < 0.13 &&
                     //~ ele.ConvVeto == 1 && ele.PFminiIso < 0.1) {
                        //~ leptoVeto = true;
               //~ }
            //~ }
         //~ }
         //~ std::cout << leptoVeto << std::endl;
      }
   }
}
