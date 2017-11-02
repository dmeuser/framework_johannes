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

Config const &cfg=Config::get();

void saveHistograms(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH1F> &hs,hist::Histograms<TH1F> &hs_pix,bool saveData)
{
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR","T5Wg","TChiWG"}){
            saver_hist.save(*hs.getHistogram(sVar,sSample),sVar+"/"+sSample);
         }
         TH1F hEFake(*hs_pix.getHistogram(sVar,"SinglePhoton"));
         hEFake.Scale(cfg.efake.f/(1.-cfg.efake.f));
         hEFake.SetLineColor(cfg.efake.color);
         // hEFake.SetFillStyle(1001);
         saver_hist.save(hEFake,sVar+"/efake");
         if (saveData) saver_hist.save(*hs.getHistogram(sVar,"SinglePhoton"),sVar+"/SinglePhoton");
         if (saveData) saver_hist.save(*hs.getHistogram(sVar,"MET"),sVar+"/MET");        
      }
   }
}

void saveHistograms(TString &sVar, io::RootFileSaver const &saver_hist,hist::Histograms<TH1F> &hs,hist::Histograms<TH1F> &hs_pix,bool saveData)
{
   for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR","T5Wg","TChiWG"}){
      saver_hist.save(*hs.getHistogram(sVar,sSample),sVar+"/"+sSample);
   }
   TH1F hEFake(*hs_pix.getHistogram(sVar,"SinglePhoton"));
   hEFake.Scale(cfg.efake.f/(1.-cfg.efake.f));
   hEFake.SetLineColor(cfg.efake.color);
   saver_hist.save(hEFake,sVar+"/efake");
   if (saveData) saver_hist.save(*hs.getHistogram(sVar,"SinglePhoton"),sVar+"/SinglePhoton");
   if (saveData) saver_hist.save(*hs.getHistogram(sVar,"MET"),sVar+"/MET");        

}
// Macro to add histograms to both, the normal Histograms instance
// and the "pixel" Histograms
#define ADD_HIST(args...) do { \
      hs_notPix.addHist(args); \
      hs_pix.addHist(args);    \
   } while (0)

enum pass_t {pass_normal,pass_pixel};

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   hist::Histograms<TH1F> hs_notPix(vsDatasubsets);
   hist::Histograms<TH1F> hs_pix(vsDatasubsets);
   int nWeights = 109;

   for (int i = 0; i < nWeights; i++){
      ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_"+std::to_string(i),";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
      ADD_HIST("pre_ph165/c_MET300/MT300/STg600/STg_"+std::to_string(i),";STg;EventsBIN"           ,3000,0,3000);
   }


   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      hs_notPix.setCurrentSample(dss.name);
      hs_pix.setCurrentSample(dss.name);
      hs_pix.setCurrentSample(dss.name);

      bool const isData=dss.isData;
      bool const isSignal=dss.isSignal;

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
      TTreeReaderValue<std::vector<tree::Photon>>   photons  (reader, "photons");
      TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
      TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
      TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<float> HTgen(reader, "genHt");
      TTreeReaderValue<bool> trigger_Ph   (reader, "HLT_Photon165_HE10_v");

      float phoPt;

      std::shared_ptr<WeightCalculator> wCalc;
      if (dss.datasetName=="ZNuNuGJets" || dss.datasetName=="ZGTo2LG"){
         wCalc=std::make_shared<WeightCalculatorZgamma>(phoPt);
      }
      else if (dss.datasetName=="WGToLNuG"){
         wCalc=std::make_shared<WeightCalculatorWgamma>();
      }

      int iEv=0;
      int processEvents=cfg.processFraction*dss.entries;
      while (reader.Next()){
         iEv++;
         if (iEv>processEvents) break;
         if (iEv%(std::max(processEvents/10,1))==0){
            io::log*".";
            io::log.flush();
         }

         float fEventWeight=*w_pu * *w_mc;
         hs_notPix.setFillWeight(fEventWeight);
         hs_pix.setFillWeight(fEventWeight);

         std::vector<tree::Photon const *> lPho,mPho,tPho,lPixPho;
         for (tree::Photon const &ph: *photons){
            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
            if (fabs(ph.p.Eta())>1.4442) continue;
            if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
            if (ph.hasPixelSeed){
               lPixPho.push_back(&ph);
            } else {
               //Attention check lossePhoton 15 or new
               if (ph.isLoose) lPho.push_back(&ph);
               if (ph.isMedium) mPho.push_back(&ph);
               if (ph.isTight)  tPho.push_back(&ph);                       
            }
         }

         // independent of "normal"/"pixel" run
         // jet related
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         int const Njets=cjets.size();

         
         bool clean_MET = true;
         
         for (auto const &jet: cjets) {
            if (jet.p.Pt() < 100) continue;            
            if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
         }
   
         if (!clean_MET) continue;

         float dPhiMETnearJet=4;
         float iJet=0;
         for (auto const &j: cjets){
            iJet++;
            const float dPhi=MET->p.DeltaPhi(j.p);
            if (std::abs(dPhi) < std::abs(dPhiMETnearJet))
               dPhiMETnearJet=dPhi;
         }

         float dPhiMETnearJet_u=4;
         for (auto const &j: cjets){
            const float dPhi=MET_JESu->p.DeltaPhi(j.p);
            if (std::abs(dPhi) < std::abs(dPhiMETnearJet_u))
               dPhiMETnearJet_u=dPhi;
         }


         float dPhiMETnearJet_d=4;
         for (auto const &j: cjets){
            const float dPhi=MET_JESd->p.DeltaPhi(j.p);
            if (std::abs(dPhi) < std::abs(dPhiMETnearJet_d))
               dPhiMETnearJet_d=dPhi;
         }
         

         // one pass for "normal" reconstruction, one for "pixel-photons"
         for (pass_t pass:{pass_normal,pass_pixel}){
            if (pass==pass_normal){
               if (lPho.empty()) continue;
            } else if (pass==pass_pixel) {
               if (!isData) continue; // no pixel-pass for MC
               if (lPixPho.empty()) continue;
               // pixel-object needs to be leading EM object
               if (!lPho.empty() && lPho[0]->p.Pt()>lPixPho[0]->p.Pt()) continue;
            }
            std::vector<tree::Photon const*> const &pho = (pass==pass_pixel) ? lPixPho : lPho;
            hist::Histograms<TH1F> &hs = (pass==pass_pixel) ? hs_pix : hs_notPix;

            // reject e->photon fakes in bkg MC (in data genParticles is empty, so nothing is rejected)
            if (pass==pass_normal && !isSignal && phys::matchesGen(*lPho[0],*genParticles,11,0.1,0.5)) continue;

            phoPt=pho[0]->p.Pt(); // set *before* wCalc->get() !
            float const met=MET->p.Pt();


            //setting weights
            fEventWeight=*w_pu * *w_mc;
            if (wCalc) fEventWeight*=wCalc->get();
            // [1-8] are muF, muR weights
            // [4] is x2, x2; [8] is x0.5, x0.5
           // if (w_pdf->size()>=9) fEventWeight*=(*w_pdf)[1];
     //       if (w_pdf->size()>=1) debug << (*w_pdf)[5];
            
            hs_notPix.setFillWeight(fEventWeight);
            hs_pix.setFillWeight(fEventWeight);

            float minDR=std::numeric_limits<float>::max();
            for (tree::Jet const &jet: *jets){
               if (jet.isLoose && jet.p.Pt()>30){
                  float const dr=jet.p.DeltaR(pho[0]->p);
                  float const dpt=fabs(jet.p.Pt()-phoPt)/phoPt;
                  if (dr<.1 && dpt<.5) continue;
                  minDR=std::min(minDR,dr);
               }
            }
            //jet und photon separated
            if (minDR < .5) continue;
            
            float const MT=phys::M_T(*pho[0],*MET);
            float const dPhiMETph=MET->p.DeltaPhi(pho[0]->p);
            float const dPhiJetPh=Njets>0 ? pho[0]->p.DeltaPhi(cjets[0].p) : 4;
            float const METdotPh=MET->p.Dot(pho[0]->p);

            // nearest jet or photon
            float dPhiMETnearJetPh=dPhiMETnearJet; // nearest jet or photon
            for (auto const &ph: pho){
               const float dPhi=MET->p.DeltaPhi(ph->p);
               if (std::abs(dPhi) < std::abs(dPhiMETnearJetPh))
                  dPhiMETnearJetPh=dPhi;
            }

            float dPhiMETnearJetPh_JESu=dPhiMETnearJet_u; // nearest jet or photon
            for (auto const &ph: pho){
               const float dPhi=MET_JESu->p.DeltaPhi(ph->p);
               if (std::abs(dPhi) < std::abs(dPhiMETnearJetPh_JESu))
                  dPhiMETnearJetPh_JESu=dPhi;
            }

            float dPhiMETnearJetPh_JESd=dPhiMETnearJet_d; // nearest jet or photon
            for (auto const &ph: pho){
               const float dPhi=MET_JESd->p.DeltaPhi(ph->p);
               if (std::abs(dPhi) < std::abs(dPhiMETnearJetPh_JESd))
                  dPhiMETnearJetPh_JESd=dPhi;
            }

            float STg=met;
            for (auto const &ph: pho){
               STg+=ph->p.Pt();
            }

            bool triggerMatch = true;

            // different preselections (for different trigger options)
            if (phoPt>180 && (!isData || (*trigger_Ph && triggerMatch))){
               if (MT > 100 && met > 100) {
                  if (met < 300 || MT < 300){
                     for (int i = 0; i < nWeights; i++){
                        if (w_pdf->size()>=9) {
                           if (i > 0) fEventWeight/=(*w_pdf)[i-1];                    
                           fEventWeight*=(*w_pdf)[i];                       
                           hs_notPix.setFillWeight(fEventWeight);
                           hs_pix.setFillWeight(fEventWeight);
                        }                     
                        hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_"+std::to_string(i),std::abs(dPhiMETnearJetPh));
                     }
                  }
                  if (met > 300 && MT > 300 && STg > 600){
                     for (int i = 0; i < nWeights; i++){
                        if (w_pdf->size()>=9) {
                           if (i > 0) fEventWeight/=(*w_pdf)[i-1];                    
                           fEventWeight*=(*w_pdf)[i];                       
                           hs_notPix.setFillWeight(fEventWeight);
                           hs_pix.setFillWeight(fEventWeight);
                        }                     
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/STg_"+std::to_string(i),STg);
                     } 
                  }
               }                                  
            }
         } // normal/pixel pass
      } // evt loop
      io::log<<"";
      hs_notPix.scaleLumi();
      hs_notPix.mergeOverflow();
      hs_pix.scaleLumi();
      hs_pix.mergeOverflow();
      file.Close();
   } // datasets

   // calling the "normal" histogram "hs" from here, since it's the most used
   hist::Histograms<TH1F> &hs = hs_notPix;
   std::vector<TString> samplesToCombine={"QCD","ZNuNuGJets","ZGTo2LG","WGToLNuG","ZNuNuJets","WLNuJets","diboson","GJets_DR","T5Wg","TChiWG","SinglePhoton","MET"};
   hs    .combineFromSubsamples(samplesToCombine);
   hs_pix.combineFromSubsamples(samplesToCombine);
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("weight_studies%.1f",cfg.processFraction*100),false);
   TCanvas can;
   can.SetLogy();
   TString sVar;
   for (int i = 0; i < nWeights;i++){
       sVar= "pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_"+std::to_string(i);
       saveHistograms(sVar,saver_hist,hs,hs_pix,true);
       sVar= "pre_ph165/c_MET300/MT300/STg600/STg_"+std::to_string(i);
       saveHistograms(sVar,saver_hist,hs,hs_pix,true); 
   }
}
