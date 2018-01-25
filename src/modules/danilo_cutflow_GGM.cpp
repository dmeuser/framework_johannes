/* module to extract signal yields */

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

#include <regex>

Config const &cfg=Config::get();
io::Logger fit_result_weights(("PU_uncert.tex"));


void runScan_80X()
{
   //initialize cutflow histogram
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hs.addHist("CutFlow_M1_1000_M2_1000",";;Ereignisse / bin",12,0,12);
   
   hs.setCurrentSample("GGM_M11000_M21000");
   
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","all");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","N_{#gamma}>0");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","p_{T}^{#gamma}>180");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","FastSim-Veto");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","clean-MET-veto");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","dR(Photon,Jet)");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","MET-cut");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","MT-cut");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","STg-cut");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","HTG-Veto");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","Lepton-Veto");
   hs.fillbinFake("CutFlow_M1_1000_M2_1000","Diphoton-Veto");
   
   TString fname=cfg.dataBasePath;
   fname+="GGM_M11000_M21000.root";
   std::map<std::string,TH1F> hSR;
   std::map<std::string,TH1F> hCR;
   std::map<std::string,TH1F> hPresel;
   std::map<std::string,TH1F> hISRWeight;

   TFile file(fname,"read");
   if (file.IsZombie()) {
      return;
   }
   io::log * ("Processing '"+fname+"' ");

   float w_pu;
   Char_t w_mc;
   std::vector<float> *w_pdf=0;
   std::vector<tree::Photon> *photons=0;
   std::vector<tree::Electron> *electrons=0;
   std::vector<tree::Muon> *muons=0;
   std::vector<tree::Jet> *jets=0;
   std::vector<tree::Particle> *genJets=0;
   std::vector<tree::GenParticle> *genParticles=0;
   UShort_t signal_m1 = 0;
   UShort_t signal_m2 = 0;
   //~ UShort_t signal_nBinos = 0;
   Int_t nGoodVertices = 0; 
   tree::MET *MET=0;
   tree::MET *genMET=0;
   TTree *tree=(TTree*)file.Get(cfg.treeName);
   tree->SetBranchAddress("pu_weight", &w_pu);
   tree->SetBranchAddress("mc_weight", &w_mc);
   tree->SetBranchAddress("pdf_weights", &w_pdf);
   tree->SetBranchAddress("photons", &photons);
   tree->SetBranchAddress("electrons", &electrons);
   tree->SetBranchAddress("muons", &muons);
   tree->SetBranchAddress("jets", &jets);
   tree->SetBranchAddress("genJets", &genJets);
   tree->SetBranchAddress("genParticles", &genParticles);
 //  tree->SetBranchAddress("modelName", &modelName);
   tree->SetBranchAddress("signal_m1", &signal_m1);
   tree->SetBranchAddress("signal_m2", &signal_m2);
 //  tree->SetBranchAddress("signal_nBinos", &signal_nBinos);  
   tree->SetBranchAddress("met", &MET);
   tree->SetBranchAddress("met_gen", &genMET);
   tree->SetBranchAddress("nGoodVertices", &nGoodVertices);
   

   std::map<std::string,int> miAcc,PVlowAll,PVhighAll,PVlowSR,PVhighSR;
   int iFastSimVeto=0;
   int iBeforeVeto=0;

   Long64_t iEvents = tree->GetEntries();
   int processEvents=cfg.processFraction*iEvents;
   for (int iEvent=0; iEvent<iEvents; iEvent++){
      if (iEvent>processEvents) break;
      if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
      tree->GetEvent(iEvent);
      
      float fEventWeight=w_pu * w_mc;
      hs.setFillWeight(fEventWeight);
      hs.fillbin("CutFlow_M1_1000_M2_1000","all");
      
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

      if (lPho.empty()) {
         continue;
      }
      hs.fillbin("CutFlow_M1_1000_M2_1000","N_{#gamma}>0");
      
      std::vector<tree::Photon const*> const &pho = lPho;
      float const phoPt=pho[0]->p.Pt(); // set *before* wCalc->get() !
      int const Nph=pho.size();
      float const dPhiMETph=MET->p.DeltaPhi(pho[0]->p);
      if (phoPt<180) {
         continue;
      }  //photon pT cut
      hs.fillbin("CutFlow_M1_1000_M2_1000","p_{T}^{#gamma}>180");
      float const MT=phys::M_T(*pho[0],*MET);
      float STg=MET->p.Pt();
      float genSTg=genMET->p.Pt();
      for (auto const &ph: pho){
         STg+=ph->p.Pt();
         genSTg+=ph->p.Pt();
      }
      
      //Vetos for GGM combination
      ///////////////
      //Lepton Veto//
      ///////////////
      bool leptoVeto = false;
      std::vector<tree::Electron const *> el_comb;
      std::vector<tree::Muon const *> mu_comb;
      
      for (tree::Electron const &ele: *electrons) {
         if (ele.p.Pt() < 25 || ele.isMedium == 0) { continue; }
         if (fabs(ele.p.Eta()) < 1.4442) {
            if (ele.r9 > 0.5 && ele.SigmaIEtaIEtaFull5x5 < 0.00998 &&
                  fabs(ele.dEtaAtVtx) < 0.00311 && fabs(ele.dPhiAtVtx) < 0.103 &&
                  ele.HoverE < 0.253 && ele.MissHits <= 1 && fabs(ele.EoverPInv) < 0.134 &&
                  ele.ConvVeto == 1 && ele.PFminiIso < 0.1) {
                     leptoVeto = true;
                     el_comb.push_back(&ele);
            }
         }
         if (fabs(ele.p.Eta()) > 1.56 && fabs(ele.p.Eta()) < 2.5) {
            if (ele.r9 > 0.8 && ele.SigmaIEtaIEtaFull5x5 < 0.0298 &&
                  fabs(ele.dEtaAtVtx) < 0.00609 && fabs(ele.dPhiAtVtx) < 0.045 &&
                  ele.HoverE < 0.0878 && ele.MissHits <= 1 && fabs(ele.EoverPInv) < 0.13 &&
                  ele.ConvVeto == 1 && ele.PFminiIso < 0.1) {
                     leptoVeto = true;
                     el_comb.push_back(&ele);
            }
         }
      }
      for (tree::Muon const &mu: *muons) {
         if (mu.p.Pt() < 25 || mu.isMedium == 0) { continue; }
         if (fabs(mu.p.Eta()) < 2.4 && mu.PFminiIso < 0.2 &&
               fabs(mu.d0) < 0.05 && fabs(mu.dZ) < 0.1) {
                  leptoVeto = true;
                  mu_comb.push_back(&mu);
         }
      }
      //leading lepton
      tree::Particle const *leadLep;
      bool leadLep_ele = false;
      bool MuAndEle = false;
      if (leptoVeto == true) {
         if (el_comb.size() > 0) {
            if (mu_comb.size() > 0) {
               if (el_comb[0]->p.Pt() > mu_comb[0]->p.Pt()) {
                  leadLep = el_comb[0];
                  leadLep_ele = true;
               }
               else leadLep = mu_comb[0];
               MuAndEle = true;
            }
            else {
               leadLep = el_comb[0];
               leadLep_ele = true;
            }
         }
         else leadLep = mu_comb[0];
      }
      //additional cuts from lepton analysis
      float invmassPhoLep = 0;
      float MT_LepMet = 0;
      if (leptoVeto == true) {
         if (leadLep->p.DeltaR(pho[0]->p) < 0.8) leptoVeto = false;
         else {
            for (tree::Electron const &ele: *electrons) {
               if (ele.p.DeltaR(pho[0]->p) < 0.3 && ele.pUncorrected.Pt() > 2.0) leptoVeto = false;
            }
            for (tree::Muon const &mu: *muons) {
               if (mu.p.DeltaR(pho[0]->p) < 0.3 && mu.p.Pt() > 3.0) leptoVeto = false;
            }
         }
         if (leptoVeto == true) {
            invmassPhoLep = phys::invmass(*pho[0],*leadLep);
            if (fabs(invmassPhoLep-91.1876) < 10 && leadLep_ele == true) leptoVeto = false;
         }
         if (leptoVeto == true) {
            if (MuAndEle == true) {
               if (phys::M_T(*el_comb[0],*MET) < 100 && phys::M_T(*mu_comb[0],*MET) < 100) leptoVeto = false;
            }
            else {
               MT_LepMet = phys::M_T(*leadLep,*MET);
               if (MT_LepMet < 100) leptoVeto = false;
            }
         }
         if (leptoVeto == true) {
            if (mPho.size() > 1) {
               if (mPho[0]->p.Pt() > 40 && mPho[1]->p.Pt() > 40){
                  leptoVeto = false;
               }
            }
         }
      }
      
      /////////////////
      //Diphoton Veto//
      /////////////////
      bool diphotonVeto = false;
      if (pho[0]->isMedium && pho[0]->p.Pt() > 40 && Nph > 1 && MET->p.Pt() > 100) {
         if(pho[1]->isMedium && pho[1]->p.Pt() > 40 && phys::invmass(*pho[0],*pho[1]) > 105 && pho[0]->p.DeltaR(pho[1]->p) > 0.3) {
            diphotonVeto = true;
         }
      }
      
      /////////////
      //EMHT-Veto//
      /////////////
      bool emhtVeto = false;
      float emht = pho[0]->p.Pt();
      for (auto const &jet: *jets) {
         if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
            if (jet.p.DeltaR(pho[0]->p) > 0.3) {
               emht += jet.p.Pt();
            }
         }
      }
      if (emht > 700 && MET->p.Pt() > 350 && fabs(dPhiMETph) > 0.3 && fabs(fabs(dPhiMETph)-TMath::Pi()) > 0.3 && phoPt > 100)  {
         emhtVeto = true;
      }

      
      // jet related
      std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);

      // fast-sim strange jet -> strange met veto
      // https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendationsICHEP16
      //so far this veto not included in xs calculation -> smaller acceptance
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
      iBeforeVeto++;
      if (vetoEvent) {
         iFastSimVeto++;
         continue;
      }
      hs.fillbin("CutFlow_M1_1000_M2_1000","FastSim-Veto");
      

      bool clean_MET = true;
      
      for (auto const &jet: cjets) {
         if (jet.p.Pt() < 100) continue;            
         if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
      }

      if (!clean_MET) {
         continue;
      }
      hs.fillbin("CutFlow_M1_1000_M2_1000","clean-MET-veto");

      float minDR=std::numeric_limits<float>::max();
      for (tree::Jet const &jet: *jets){
         if (jet.isLoose && jet.p.Pt()>30){
            float const dr=jet.p.DeltaR(pho[0]->p);
            float const dpt=fabs(jet.p.Pt()-phoPt)/phoPt;
            if (dr<.1 && dpt<.5) continue;
            minDR=std::min(minDR,dr);
         }
      }
      if (minDR < .5){
         continue;
      }
      hs.fillbin("CutFlow_M1_1000_M2_1000","dR(Photon,Jet)");

      float dPhiMETnearJet=4;
      float iJet=0;
      for (auto const &j: cjets){
         iJet++;
         const float dPhi=MET->p.DeltaPhi(j.p);
         if (std::abs(dPhi) < std::abs(dPhiMETnearJet))
            dPhiMETnearJet=dPhi;
      }
      
      // nearest jet or photon
      float dPhiMETnearJetPh=dPhiMETnearJet; // nearest jet or photon
      for (auto const &ph: pho){
         const float dPhi=MET->p.DeltaPhi(ph->p);
         if (std::abs(dPhi) < std::abs(dPhiMETnearJetPh))
            dPhiMETnearJetPh=dPhi;
      }
      
      //Signal Region selection
      if (MET->p.Pt()>300) {
         hs.fillbin("CutFlow_M1_1000_M2_1000","MET-cut");
         if (MT>300) {
            hs.fillbin("CutFlow_M1_1000_M2_1000","MT-cut");
            if (STg>600) {
               hs.fillbin("CutFlow_M1_1000_M2_1000","STg-cut");
               if (emhtVeto == false) {
                  hs.fillbin("CutFlow_M1_1000_M2_1000","HTG-Veto");
                  if (leptoVeto == false) {
                     hs.fillbin("CutFlow_M1_1000_M2_1000","Lepton-Veto");
                     if (diphotonVeto == false) {
                        hs.fillbin("CutFlow_M1_1000_M2_1000","Diphoton-Veto");
                     }
                  }
               }
            }
         }
      }
      
   } // evt loop
   
   std::cout<<hs.getVariableNames()[0]<<std::endl;
   io::RootFileSaver saver("plots.root",TString::Format("danilo_cutflow_GGM%.1f",cfg.processFraction*100));
   TCanvas can;
   TH1F *cutflow;
   cutflow = hs.getHistogram("CutFlow_M1_1000_M2_1000","GGM_M11000_M21000");
   std::cout<<cutflow->GetBinContent(1)<<std::endl;
   cutflow->Draw("E1");
   saver.save(can,"CutFlow_M1_1000_M2_1000",false);
}

extern "C"
void run()
{
   runScan_80X();
}
