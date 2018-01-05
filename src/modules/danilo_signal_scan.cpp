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

enum Scan_t
{
   GGM,
   TChiWg,
   TChiNg,   
   T5gg,
   T5Wg,   
   T6Wg,
   T6gg,  
};


std::map<int,float> getXsecs(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==GGM)  path += "/xsec_wino-bino.csv";
   else if (scan==TChiWg) path += "/xsec_N2C1_wino.csv";
   else if (scan==TChiNg) path += "/xsec_comb_wino.csv";   
   else if (scan==T5Wg) path += "/xsec_gluglu.csv";
   else if (scan==T5gg) path += "/xsec_gluglu.csv";
   else if (scan==T6Wg) path += "/xsec_sqsq.csv";   //squarks! update
   else if (scan==T6gg) path += "/xsec_sqsq.csv";   //squarks! update
   unsigned const nCol=scan==GGM ? 4 : 3;
   // T5Wg also contains gg, WW -> actually only 1/2 of generated events
   // instead, multiply xsec by 2
   int const multiplier= 1;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==GGM) {
            // for key, combine M2 and M1 with a 0 in between
            // debug<<((int)values[0]*100000+(int)values[1]);
            m_M_XS[(int)values[0]*100000+(int)values[1]]=values[2]/1000.; // convert to pb
         } else if (scan==TChiWg || scan==TChiNg) {
            m_M_XS[(int)values[0]]=values[1]/1000.; // convert to pb
         } else {
            m_M_XS[(int)values[0]]=values[1]*multiplier;
         }
      }
   }
   return m_M_XS;
}

std::string getModelName(Scan_t scan, UShort_t signal_m1, UShort_t signal_m2 = 0, UShort_t signal_nBinos = 0)
{
   std::string modelName = "";
   if (scan==TChiWg)  modelName = "TChiWG_";
   else if (scan==TChiNg)  modelName = "TChiNG_";
   else if (scan==T5Wg)  modelName = "T5Wg_";
   else if (scan==T5gg)  modelName = "T5gg_";
   else if (scan==T6Wg)  modelName = "T6Wg_";
   else if (scan==T6gg)  modelName = "T6gg_";
   else if (scan==GGM)  modelName = "GGM_";   
   if (signal_m1) modelName += std::to_string(signal_m1);
   if (signal_m2) modelName += "_" + std::to_string(signal_m2);
  
   return modelName;
}


std::pair<int,int> getMasses(std::string fileName,Scan_t scan)
{
   std::smatch m;
   std::regex e;
  // if (scan==GGM)  e=".*_M2_(.*)_M1_(.*)\\.root";
   if (scan==GGM)  e=".*_(.*)_(.*)";
   else if (scan==TChiWg)  e="TChiWG_(.*)";
   else if (scan==TChiNg)  e="TChiNG_(.*)";
   else if (scan==T5Wg) e="T5.*_(.*)_(.*)";
   else if (scan==T5gg) e="T5.*_(.*)_(.*)";
   else if (scan==T6Wg) e="T6.*_(.*)_(.*)";
   else if (scan==T6gg) e="T6.*_(.*)_(.*)";

   std::regex_search (fileName,m,e);
   if (scan==TChiWg || scan ==TChiNg) {
      assert(m.size()==2);
      // io::log*m[1]>>m[2];
      return std::make_pair(std::stoi(m[1]),0);
   }

   assert(m.size()==3);
   // io::log*m[1]>>m[2];
   return std::make_pair(std::stoi(m[1]),std::stoi(m[2]));
}

void runScan_80X(Scan_t scan)
{
   std::map<int,float> mXsecs=getXsecs(scan);
   TString fname=cfg.dataBasePath;
   if (scan==TChiWg) fname+="SMS-TChiWG.root";
   else if (scan==TChiNg) fname+="SMS-TChiNG_BF50N50G.root";  
   else if (scan==T5Wg) fname+="SMS-T5Wg.root";
   else if (scan==T5gg) fname+="SMS-T5gg.root";
   else if (scan==T6Wg) fname+="SMS-T6Wg.root";
   else if (scan==T6gg) fname+="SMS-T6gg.root";
   else if (scan==GGM) fname+="GGM.root";
   else debug<<"unsupported scan!";
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

   std::string model= "";
   std::map<std::string,int> miAcc,PVlowAll,PVhighAll,PVlowSR,PVhighSR;
   int iFastSimVeto=0;
   int iBeforeVeto=0;

   Long64_t iEvents = tree->GetEntries();
   for (int iEvent=0; iEvent<iEvents; iEvent++){
      if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
      tree->GetEvent(iEvent);
      model = getModelName(scan, signal_m1, signal_m2);
      
 //              debug << ".......................hier";
 /*
      if (scan==T5Xg) {
         // find gg and Wg events, veto WW
         int Ng=0;
         int id;
         for (tree::GenParticle gp: *genParticles) {
            if (!gp.fromHardProcess) continue;
            id=std::abs(gp.pdgId);
            if (id==22) Ng++;
         }
         if (Ng==0) continue; // WW
         else if (Ng==2) { // gg
            // rename model
            model="T5gg"+model.substr(model.find('_'));
         } // else Wg: already called "T5Wg"
      }
      */
      if (hSR.count(model)<1) {
         
         hSR[model]=hist::fromWidths((model+"SR").c_str(),";STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         hSR[model+"SRErrISR"]=hist::fromWidths((model+"SRErrISR").c_str(),";STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});         
         hSR[model+"_gen"]=hist::fromWidths((model+"genSR").c_str(),";gen STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         hCR[model]=hist::fromWidths((model+"CR").c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});

         hPresel[model]=hist::fromWidths((model).c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});
         hPresel[model+"_mu2"]=hist::fromWidths((model+"_mu2").c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});
         hPresel[model+"_mu05"]=hist::fromWidths((model+"_mu05").c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});
         hSR[model+"_mu2"]=hist::fromWidths((model+"SR"+"_mu2").c_str(),";STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         hSR[model+"_mu05"]=hist::fromWidths((model+"SR"+"_mu05").c_str(),";STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         hISRWeight[model+"_before"]=hist::fromWidths((model+"_before").c_str(),";phoPt;EventsBIN",{0.,700.},{50.});
         hISRWeight[model+"_after"]=hist::fromWidths((model+"_after").c_str(),";phoPt;EventsBIN",{0.,700.},{50.});
         hISRWeight[model+"_afterErr"]=hist::fromWidths((model+"_afterErr").c_str(),";phoPt;EventsBIN",{0.,700.},{50.});
         

         miAcc[model]=0;
         PVlowAll[model]=0;
         PVhighAll[model]=0;
         PVlowSR[model]=0;              
         PVhighSR[model]=0;
      }
      float fEventWeight=w_pu * w_mc;
      float fEventWeightError = fEventWeight;

      //ISR weighting
      TVector3 EWKinoPair;
      double EWKinoPairPt = 0;
      if (scan == TChiWg || scan == TChiNg){
         for (tree::GenParticle &genP: *genParticles){
            if (fabs(genP.pdgId) > 1000022){
            EWKinoPair += genP.p;
            }
         }
         EWKinoPairPt = EWKinoPair.Pt();
         hISRWeight[model+"_before"].Fill(EWKinoPairPt,fEventWeight);
         if (EWKinoPairPt < 50) {
            fEventWeight*=1;
         }
         else if (EWKinoPairPt < 100) {
            fEventWeight*=1.052;
         }
         else if (EWKinoPairPt < 150) {
            fEventWeight*=1.179;
         }
         else if (EWKinoPairPt < 200) {
            fEventWeight*=1.150;
         }
         else if (EWKinoPairPt < 300) {
            fEventWeight*=1.057;
         }
         else if (EWKinoPairPt < 400) {
            fEventWeight*=1.;
         }
         else if (EWKinoPairPt < 600) {
            fEventWeight*=0.912;
         }
         else{
            fEventWeight*=0.783;
         }                 
      }
      else{
         hISRWeight[model+"_before"].Fill(EWKinoPairPt,fEventWeight);
      }

      if (nGoodVertices < 17){
         PVlowAll[model]++;
      }
      else {
         PVhighAll[model]++;
      }
      
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
         hISRWeight[model+"_after"].Fill(-1,fEventWeight);
         hISRWeight[model+"_afterErr"].Fill(-1,fEventWeightError);
         continue;
      }

      std::vector<tree::Photon const*> const &pho = lPho;
      float const phoPt=pho[0]->p.Pt(); // set *before* wCalc->get() !
      int const Nph=pho.size();
      float const dPhiMETph=MET->p.DeltaPhi(pho[0]->p);
      if (phoPt<180) {
         hISRWeight[model+"_after"].Fill(-1,fEventWeight);
         hISRWeight[model+"_afterErr"].Fill(-1,fEventWeightError);
         continue;
      }  //photon pT cut
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
      //Apply Vetos for GGM combination
      if (leptoVeto == true || diphotonVeto == true || emhtVeto == true) continue;


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
         hISRWeight[model+"_after"].Fill(-1,fEventWeight);
         hISRWeight[model+"_afterErr"].Fill(-1,fEventWeightError);
         continue;
      }

      bool clean_MET = true;
      
      for (auto const &jet: cjets) {
         if (jet.p.Pt() < 100) continue;            
         if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
      }

      if (!clean_MET) {
         hISRWeight[model+"_after"].Fill(-1,fEventWeight);
         hISRWeight[model+"_afterErr"].Fill(-1,fEventWeightError);
         continue;
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
      if (minDR < .5){
         hISRWeight[model+"_after"].Fill(-1,fEventWeight);
         hISRWeight[model+"_afterErr"].Fill(-1,fEventWeightError);
         continue;
      }

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
      hPresel[model].Fill(dPhiMETnearJetPh,fEventWeight);
      if (scan == GGM){
         hPresel[model+"_mu2"].Fill(dPhiMETnearJetPh,fEventWeight);
         hPresel[model+"_mu05"].Fill(dPhiMETnearJetPh,fEventWeight);
      } else {
         hPresel[model+"_mu2"].Fill(dPhiMETnearJetPh,fEventWeight*w_pdf->at(4));
         hPresel[model+"_mu05"].Fill(dPhiMETnearJetPh,fEventWeight*w_pdf->at(8));
      }
      
      if (MET->p.Pt()>100 && MT>100) {
         if (MET->p.Pt()<300 || MT<300) {
            hCR[model].Fill(dPhiMETnearJetPh,fEventWeight);
         }
      }
      
      //adjust genMET study for MET uncertainty -> selection with gen met
      if (MET->p.Pt()>300 && MT>300 && STg>600) {
         miAcc[model]++;
         if (nGoodVertices < 17){
            PVlowSR[model]++;
         }
         else {
            PVhighSR[model]++;
         }

         hISRWeight[model+"_after"].Fill(EWKinoPairPt,fEventWeight);
         hISRWeight[model+"_afterErr"].Fill(EWKinoPairPt,fEventWeightError);
         hSR[model].Fill(STg,fEventWeight);
         hSR[model+"SRErrISR"].Fill(STg,fEventWeightError);
         hSR[model+"_gen"].Fill(genSTg,fEventWeight);
         if (scan == GGM){
            hSR[model+"_mu2"].Fill(STg,fEventWeight);
            hSR[model+"_mu05"].Fill(STg,fEventWeight);            
         } else {
            hSR[model+"_mu2"].Fill(STg,fEventWeight*w_pdf->at(4));
            hSR[model+"_mu05"].Fill(STg,fEventWeight*w_pdf->at(8));
         }

      }
      else {
      hISRWeight[model+"_after"].Fill(-1,fEventWeight);
      hISRWeight[model+"_afterErr"].Fill(-1,fEventWeightError);
      }

      
   } // evt loop
   model = getModelName(scan, 640, 630);
   io::log<<"";
   io::log*"vetoed "*iFastSimVeto*"/"*iBeforeVeto>>"events";

   io::log<<"";
  // io::log*"rel. acceptance low PV set "*PVlowSR*"/"*PVlowAll>>"";  
   io::log<<"";
 //  io::log*"rel. acceptance high PV set "*PVhighSR*"/"*PVhighAll>>"";
 
   TString sScan="unkown_scan";
   if (scan==GGM)  sScan="GGM";
   else if (scan==TChiWg) sScan="TChiWg";
   else if (scan==TChiNg) sScan="TChiNg";
   else if (scan==T5Wg) sScan="T5Wg";
   else if (scan==T5gg) sScan="T5gg";
   else if (scan==T6Wg) sScan="T6Wg";
   else if (scan==T6gg) sScan="T6gg";

 //  for (auto &it: PVlowSR) std::cout << it.first << " : " << ((1.*it.second/PVlowAll[it.first]) - (1.*PVhighSR[it.first]/PVhighAll[it.first]))*100/2. << std::endl; //balanced deviation because of PU in percent

   TString fill_line;

   double uncert_PU = 0.;
   for (auto &it: PVlowSR) {
      uncert_PU = fabs(((1.*it.second/PVlowAll[it.first]) - (1.*PVhighSR[it.first]/PVhighAll[it.first]))*100/2.);
      fill_line = it.first +" : " + std::to_string(uncert_PU);
      fit_result_weights<< fill_line;
   }
    
   io::RootFileReader dataReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("danilo_distributions%.1f",cfg.processFraction*100));
   TH1F hData(*dataReader.read<TH1F>("pre_ph165/c_MET100/MT100/METl300vMTl300/exclusiv/absphiMETnJetPh/SinglePhoton"));
   float const nData=hData.Integral();

   std::map<TString,TGraph2D> grAcc;
   std::map<TString,TGraph2D> grScaleUnc;
   std::map<TString,TGraph2D> grCont;


   io::RootFileSaver saver_hist(TString::Format("signal_scan_exclusiv_%s.root",cfg.treeVersion.Data()),"",false);
   TString sVar;
   
   for (auto const &map: hSR) {
      std::string const model=map.first;
      if (model.find("_mu")!=std::string::npos) continue; // variations are handeled in nominal iteration
      if (model.find("_gen")!=std::string::npos) continue; // gen case handeled in nominal iteration
      if (model.find("SRErrISR")!=std::string::npos) continue; // gen case handeled in nominal iteration     
      // scale lumi
      TH1F* hcf;
      if (model.find("T5gg")!=std::string::npos) {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlowT5Wg"+model.substr(model.find('_'))).c_str());
         sScan="T5gg";
      } else if (model.find("T6gg")!=std::string::npos) {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlowT6Wg"+model.substr(model.find('_'))).c_str());
         sScan="T6gg";
      } else if (model.find("GGM")!=std::string::npos) {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlowGMSB"+model.substr(model.find('_'))).c_str());
         sScan="GGM";         
      } else {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlow"+model).c_str());
         if (scan==T5Wg) sScan="T5Wg";
     //    std::cout << model << ".root" << std::endl;
      }
      assert(hcf);
      float Ngen=hcf->GetBinContent(2);
      std::cout << model << std::endl;
 //     if (scan==T5Wg) {  //after discussion with conveners we need to assume 50% to W and 50% to gamma, so that we need the full scan
         // consider that gg, gW, and WW events were generated
  //       Ngen/=2.0;
 //     }
      if (scan==T5gg) {
         Ngen/=4.0;
      }
    //  else if (scan==T6Wg) {
   //      Ngen/=2.0;
   //   }
      else if (scan==T6gg) {
         Ngen/=4.0;
      }

      double scaleUnc,normalization, ISR_norm;

  //    debug << model << hISRWeight[model+"_before"].Integral(0,-1) << hISRWeight[model+"_after"].Integral(0,-1) << hISRWeight[model+"_before"].Integral(0,-1)/hISRWeight[model+"_after"].Integral(0,-1) ;
      ISR_norm = hISRWeight[model+"_before"].Integral(0,-1)/hISRWeight[model+"_after"].Integral(0,-1);
      
      normalization=hPresel[model].Integral()/hPresel[model+"_mu2"].Integral();
      scaleUnc=std::fabs(1-normalization*hSR[model+"_mu2"].Integral()/hSR[model].Integral());
      normalization=hPresel[model].Integral()/hPresel[model+"_mu05"].Integral();
      scaleUnc=std::max(scaleUnc,std::fabs(1-normalization*hSR[model+"_mu05"].Integral()/hSR[model].Integral()));

      hSR[model].Scale(ISR_norm); //scale after scaleUnc determination because otherwise Eventweights don't cancel
      hSR[model+"SRErrISR"].Scale(ISR_norm);   
      hSR[model+"_gen"].Scale(ISR_norm);
      hCR[model].Scale(ISR_norm);

      std::pair<int,int> const masses=getMasses(model,scan);
      int m ;
      float xs;
      
      if (scan == GGM){
         m=masses.first*100000 + masses.second;
         xs=mXsecs[m];
      }
      else {
         m=masses.first;
         xs=mXsecs[m];
      }
      float const w=xs/Ngen*cfg.lumi;

      hSR[model].Scale(w);
      hSR[model+"SRErrISR"].Scale(w);      
      hSR[model+"_gen"].Scale(w);
      hCR[model].Scale(w);

      hist::mergeOverflow(hSR[model]);
      hist::mergeOverflow(hSR[model+"SRErrISR"]);  
      hist::mergeOverflow(hSR[model+"_gen"]);
      hist::mergeOverflow(hCR[model]);
      hist::mergeOverflow(hISRWeight[model+"_before"]);
      hist::mergeOverflow(hISRWeight[model+"_after"]);
      hist::mergeOverflow(hISRWeight[model+"_afterErr"]);
      
      hSR[model].Scale(cfg.trigger_eff_Ph);
      hSR[model+"SRErrISR"].Scale(cfg.trigger_eff_Ph);   
      hSR[model+"_gen"].Scale(cfg.trigger_eff_Ph);
      hCR[model].Scale(cfg.trigger_eff_Ph);
      
      sVar = "pre_ph165/c_MET300/MT300/STg";
      saver_hist.save(hSR[model],sVar+"/"+model);
      saver_hist.save(hSR[model+"SRErrISR"],sVar+"/"+model+"SRErrISR");  
      saver_hist.save(hSR[model+"_gen"],sVar+"/"+model+"_gen");
      
      sVar = "EWKinoPairPt";      
      saver_hist.save(hISRWeight[model+"_before"],sVar+"/"+model+"_before");
      sVar = "pre_ph165/c_MET300/MT300/STg600/EWKinoPairPt";  
      saver_hist.save(hISRWeight[model+"_after"],sVar+"/"+model+"_after");
      saver_hist.save(hISRWeight[model+"_afterErr"],sVar+"/"+model+"_afterErr");

      // acceptance
      float x,y;
      if (scan == GGM) {
         x=masses.second;
         y=masses.first;
      }
      else {      
         x=masses.first;
         y=masses.second;
      }
      grAcc[sScan].SetPoint(grAcc[sScan].GetN(),x,y,float(miAcc[model])/Ngen);
      grScaleUnc[sScan].SetPoint(grScaleUnc[sScan].GetN(),x,y,scaleUnc);

      // contamination
      hCR[model].Scale(1./nData);
      sVar="pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/exclusiv";
      saver_hist.save(hCR[model],sVar+"/"+model);
      grCont[sScan].SetPoint(grCont[sScan].GetN(),x,y,hCR[model].Integral());
   }
   file.Close();

   for (auto const& p: grAcc) {
      sVar = "pre_ph165/c_MET300/MT300/STg/exclusiv";
      saver_hist.save(grAcc[p.first],sVar+"/"+p.first+"_acceptance");
      saver_hist.save(grScaleUnc[p.first],sVar+"/"+p.first+"_scaleUnc");

      sVar="pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/exclusiv";
      saver_hist.save(grCont[p.first],sVar+"/"+p.first+"_contamination");
   }
}

extern "C"
void run()
{
 //   runScan_80X(TChiWg);
 //   runScan_80X(TChiNg);
    runScan_80X(T5Wg);
 //   runScan_80X(T5gg);
 //   runScan_80X(T6Wg);
 //   runScan_80X(T6gg);
  //runScan_80X(GGM);
}
