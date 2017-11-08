#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>

Config const &cfg=Config::get();

enum Scan_t
{
   T5gg,
   T5Wg,
   T6gg,
   T6Wg,
   TChiWG,
   TChiNg,   
   GGM,
};


void runScan(Scan_t scan)
{

   TString dsName;
   if (scan==GGM) dsName="GGM";
   if (scan==T5gg) dsName="T5gg";
   if (scan==T5Wg) dsName="T5Wg";
   if (scan==T6gg) dsName="T6gg";
   if (scan==T6Wg) dsName="T6Wg";
   if (scan==TChiWG) dsName="TChiWG";
   if (scan==TChiNg) dsName="TChiNg";
   

   float fPt=0;
   float fDR=0;
   float fPr=0;
   float fMET=0;
   float fMT=0;
   float fSR=0;
   float fSR1=0;
   float fSR2=0;
   float fSR3=0;
   float fSR4=0;  
   // int i=0;
   for (Datasubset const &dss: cfg.datasets.getDataset(dsName).subsets){
      // i++;
      // if (i>3) break;
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.name+"'... ");

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<tree::Photon>>   photons  (reader, "photons");
      TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
      TTreeReaderValue<std::vector<tree::Particle>> genJets  (reader, "genJets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles  (reader, "genParticles");
      TTreeReaderValue<tree::MET> MET(reader, "met");

      while (reader.Next()){
         float fEventWeight=*w_pu * *w_mc;
         std::vector<tree::Photon const *> lPho,mPho,tPho,lPixPho;
         for (tree::Photon const &ph: *photons){
            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
            if (fabs(ph.p.Eta())>1.4442) continue;
            if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
            if (ph.hasPixelSeed){
               lPixPho.push_back(&ph);
            } else {
               lPho.push_back(&ph);
               if (ph.isMedium) mPho.push_back(&ph);
               if (ph.isTight)  tPho.push_back(&ph);
            }
         }
         if (lPho.empty()) continue;


         std::vector<tree::Photon const*> const &pho = lPho;
         float const phoPt=pho[0]->p.Pt(); // set *before* wCalc->get() !
         float const MT=phys::M_T(*pho[0],*MET);
         float STg=MET->p.Pt();
         for (auto const &ph: pho){
            STg+=ph->p.Pt();
         }

         // jet related
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);

         bool clean_MET = true;

         for (auto const &jet: cjets) {
            if (jet.p.Pt() < 100) continue;            
            if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
         }
         if (clean_MET == false) continue;

         if (scan!=GGM) {
            // fast-sim strange jet -> strange met veto
            // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendationsICHEP16
            bool vetoEvent=false;
            for (tree::Jet const &j: cjets) {
               if (std::fabs(j.p.Eta())>2.5) continue;
               if (j.chf > 0.1) continue;
               bool matched=false;
               for (tree::Particle const &gj: *genJets) {
                  if (j.p.DeltaR(gj.p) < 0.3) matched=true;
               }
               if (!matched) {
                  vetoEvent=true;
                  break;
               }
            }
            if (vetoEvent) {
               continue;
            }
         }

         float minDR=std::numeric_limits<float>::max();
         for (tree::Jet const &jet: *jets){
            if (jet.isLoose && jet.p.Pt()>30){
               float const dr=jet.p.DeltaR(pho[0]->p);
               float const dpt=fabs(jet.p.Pt()-phoPt)/phoPt;
               if (dr<.1 && dpt<.5) continue;
               minDR=std::min(minDR,dr);
            }
         }
         if (phoPt<180) continue;
         fPt+=fEventWeight;
         if (minDR < .5) continue;
         fDR+=fEventWeight;
         fPr+=fEventWeight;

         if (MET->p.Pt() > 300) {
            fMET+=fEventWeight;
            if (MT>300) {
               fMT+=fEventWeight;
               if (STg>600){
                  fSR+=fEventWeight;
               }
               if (STg>1300)     fSR4+=fEventWeight;
               else if (STg>1000) fSR3+=fEventWeight;              
               else if (STg>800) fSR2+=fEventWeight;
               else if (STg>600) fSR1+=fEventWeight;
            }
         }
      } // evt loop
      TH1F *hcf=(TH1F*)file.Get("TreeWriter/hCutFlow");
      float fI=hcf->GetBinContent(3);
      io::log<<dsName;
      io::log<<dss.xsec;
      float const Ngen=fI;
      const float w=dss.xsec/Ngen*cfg.lumi*cfg.trigger_eff_Ph;
      // io::log<<TString::Format("pre %.2f\nSR  %.2f\nST  %.2f",fPr/fI,fSR/fI,fST/fI);
      io::log<<TString::Format("ini %.2f & (%.0f$\\%%$)\n"
                               "pre %.2f & (%.0f$\\%%$)\n"
                               "MET  %.2f & (%.0f$\\%%$)\n"
                               "MT  %.2f & (%.0f$\\%%$)\n"
                               "SR %.2f & (%.0f$\\%%$)\n"
                               "SR1 %.2f & (%.0f$\\%%$)\n"
                               "SR2 %.2f & (%.0f$\\%%$)\n"
                               "SR3 %.2f & (%.0f$\\%%$)\n"
                               "SR4 %.2f & (%.0f$\\%%$)",
                               fI*w,fI/fI*100,
                               fPr*w,fPr/fI*100,
                               fMET*w,fMET/fI*100,
                               fMT*w,fMT/fI*100,
                               fSR*w,fSR/fI*100,
                               fSR1*w,fSR1/fI*100,
                               fSR2*w,fSR2/fI*100,
                               fSR3*w,fSR3/fI*100,
                               fSR4*w,fSR4/fI*100
         );
      file.Close();
   } // datasets
}

extern "C"
void run()
{
//   runScan(GGM);
 //  runScan(T5gg);
   runScan(T5Wg);
   runScan(TChiWG);   
}
