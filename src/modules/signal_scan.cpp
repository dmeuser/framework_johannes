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
   if (signal_m1) modelName += std::to_string(signal_m1);
   if (signal_m2) modelName += "_" + std::to_string(signal_m2);   
   
   return modelName;
}


std::pair<int,int> getMasses(std::string fileName,Scan_t scan)
{
   std::smatch m;
   std::regex e;
   if (scan==GGM)  e=".*_M2_(.*)_M1_(.*)\\.root";
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
   else if (scan==TChiNg) fname+="SMS-TChiNG.root";  
   else if (scan==T5Wg) fname+="SMS-T5Wg.root";
   else if (scan==T5gg) fname+="SMS-T5gg.root";
   else if (scan==T6Wg) fname+="SMS-T6Wg.root";
   else if (scan==T6gg) fname+="SMS-T6gg.root";
   else debug<<"unsupported scan!";
   std::map<std::string,TH1F> hSR;
   std::map<std::string,TH1F> hCR;
   std::map<std::string,TH1F> hPresel;

   TFile file(fname,"read");
   if (file.IsZombie()) {
      return;
   }
   io::log * ("Processing '"+fname+"' ");

   float w_pu;
   Char_t w_mc;
   std::vector<float> *w_pdf=0;
   std::vector<tree::Photon> *photons=0;
   std::vector<tree::Jet> *jets=0;
   std::vector<tree::Particle> *genJets=0;
   std::vector<tree::GenParticle> *genParticles=0;
   UShort_t signal_m1 = 0;
   UShort_t signal_m2 = 0;
   UShort_t signal_nBinos = 0;
   tree::MET *MET=0;
   tree::MET *genMET=0;

   TTree *tree=(TTree*)file.Get(cfg.treeName);
   tree->SetBranchAddress("pu_weight", &w_pu);
   tree->SetBranchAddress("mc_weight", &w_mc);
   tree->SetBranchAddress("pdf_weights", &w_pdf);
   tree->SetBranchAddress("photons", &photons);
   tree->SetBranchAddress("jets", &jets);
   tree->SetBranchAddress("genJets", &genJets);
   tree->SetBranchAddress("genParticles", &genParticles);
 //  tree->SetBranchAddress("modelName", &modelName);
   tree->SetBranchAddress("signal_m1", &signal_m1);
   tree->SetBranchAddress("signal_m2", &signal_m2);
 //  tree->SetBranchAddress("signal_nBinos", &signal_nBinos);  
   tree->SetBranchAddress("met", &MET);
   tree->SetBranchAddress("met_gen", &genMET);

   std::string model= "";
   std::map<std::string,int> miAcc;
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
         hSR[model+"_gen"]=hist::fromWidths((model+"genSR").c_str(),";gen STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         hCR[model]=hist::fromWidths((model+"CR").c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});

         hPresel[model]=hist::fromWidths((model).c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});
         hPresel[model+"_mu2"]=hist::fromWidths((model+"_mu2").c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});
         hPresel[model+"_mu05"]=hist::fromWidths((model+"_mu05").c_str(),";absphiMETnJetPh;EventsBIN",{0,.8,3.2},{.2,.4});
         hSR[model+"_mu2"]=hist::fromWidths((model+"SR"+"_mu2").c_str(),";STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         hSR[model+"_mu05"]=hist::fromWidths((model+"SR"+"_mu05").c_str(),";STg;EventsBIN",{600,800,1000,1300,1600},{200,200,300,300});
         miAcc[model]=0;

      }

      float fEventWeight=w_pu * w_mc;

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
      if (phoPt<180) continue;
      float const MT=phys::M_T(*pho[0],*MET);
      float STg=MET->p.Pt();
      float genSTg=genMET->p.Pt();
      for (auto const &ph: pho){
         STg+=ph->p.Pt();
         genSTg+=ph->p.Pt();
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

      bool clean_MET = true;
      
      for (auto const &jet: cjets) {
         if (jet.p.Pt() < 100) continue;            
         if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
      }

      if (!clean_MET) continue;
      

      float minDR=std::numeric_limits<float>::max();
      for (tree::Jet const &jet: *jets){
         if (jet.isLoose && jet.p.Pt()>30){
            float const dr=jet.p.DeltaR(pho[0]->p);
            float const dpt=fabs(jet.p.Pt()-phoPt)/phoPt;
            if (dr<.1 && dpt<.5) continue;
            minDR=std::min(minDR,dr);
         }
      }
      if (minDR < .5) continue;

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
      hPresel[model+"_mu2"].Fill(dPhiMETnearJetPh,fEventWeight*w_pdf->at(4));
      hPresel[model+"_mu05"].Fill(dPhiMETnearJetPh,fEventWeight*w_pdf->at(8));

      if (MET->p.Pt()>100 && MT>100) {
         if (MET->p.Pt()<300 || MT<300) {
            hCR[model].Fill(dPhiMETnearJetPh,fEventWeight);
         }
      }
      //adjust genMET study for MET uncertainty -> selection with gen met
      if (MET->p.Pt()>300) {
         if (MT>300) {
            hSR[model].Fill(STg,fEventWeight);
            hSR[model+"_gen"].Fill(genSTg,fEventWeight);
            hSR[model+"_mu2"].Fill(STg,fEventWeight*w_pdf->at(4));
            hSR[model+"_mu05"].Fill(STg,fEventWeight*w_pdf->at(8));
            if (STg>600) miAcc[model]++;
         }
      }
   } // evt loop
   io::log<<"";
   io::log*"vetoed "*iFastSimVeto*"/"*iBeforeVeto>>"events";

   io::RootFileReader dataReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   TH1F hData(*dataReader.read<TH1F>("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/SinglePhoton"));
   float const nData=hData.Integral();

   std::map<TString,TGraph2D> grAcc;
   std::map<TString,TGraph2D> grScaleUnc;
   std::map<TString,TGraph2D> grCont;
   TString sScan="unkown_scan";
   if (scan==GGM)  sScan="GGM";
   else if (scan==TChiWg) sScan="TChiWg";
   else if (scan==TChiNg) sScan="TChiNg";
   else if (scan==T5Wg) sScan="T5Wg";
   else if (scan==T5gg) sScan="T5gg";
   else if (scan==T6Wg) sScan="T6Wg";
   else if (scan==T6gg) sScan="T6gg";

   io::RootFileSaver saver_hist(TString::Format("signal_scan_%s.root",cfg.treeVersion.Data()),"",false);
   TString sVar;
   for (auto const &map: hSR) {
      std::string const model=map.first;
      if (model.find("_mu")!=std::string::npos) continue; // variations are handeled in nominal iteration
      if (model.find("_gen")!=std::string::npos) continue; // gen case handeled in nominal iteration
      // scale lumi
      TH1F* hcf;
      if (model.find("T5gg")!=std::string::npos) {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlowT5Wg"+model.substr(model.find('_'))).c_str());
         sScan="T5gg";
      } else if (model.find("T6gg")!=std::string::npos) {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlowT6Wg"+model.substr(model.find('_'))).c_str());
         sScan="T6gg";
      } else {
         hcf=(TH1F*)file.Get(("TreeWriter/hCutFlow"+model).c_str());
         if (scan==T5Wg) sScan="T5Wg";
       //  std::cout << model << ".root" << std::endl;
      }
      assert(hcf);
      float Ngen=hcf->GetBinContent(2);
      if (scan==T5Wg) {
         // consider that gg, gW, and WW events were generated
       //  if (sScan=="T5gg") Ngen/=4.0;
         Ngen/=2.0;
      }
      else if (scan==T5gg) {
         // consider that gg, gW, and WW events were generated
       //  if (sScan=="T5gg") Ngen/=4.0;
         Ngen/=4.0;
      }
      else if (scan==T6Wg) {
         // consider that gg, gW, and WW events were generated
       //  if (sScan=="T5gg") Ngen/=4.0;
         Ngen/=2.0;
      }
      else if (scan==T6gg) {
         // consider that gg, gW, and WW events were generated
       //  if (sScan=="T5gg") Ngen/=4.0;
         Ngen/=4.0;
      }

      double scaleUnc,normalization;
      normalization=hPresel[model].Integral()/hPresel[model+"_mu2"].Integral();
      scaleUnc=std::fabs(1-normalization*hSR[model+"_mu2"].Integral()/hSR[model].Integral());
      normalization=hPresel[model].Integral()/hPresel[model+"_mu05"].Integral();
      scaleUnc=std::max(scaleUnc,std::fabs(1-normalization*hSR[model+"_mu05"].Integral()/hSR[model].Integral()));

      std::pair<int,int> const masses=getMasses(model,scan);
      int const m=masses.first;
      float const xs=mXsecs[m];
      float const w=xs/Ngen*cfg.lumi;
      hSR[model].Scale(w);
      hSR[model+"_gen"].Scale(w);
      hCR[model].Scale(w);
      hist::mergeOverflow(hSR[model]);
      hist::mergeOverflow(hSR[model+"_gen"]);
      hist::mergeOverflow(hCR[model]);
      hSR[model].Scale(cfg.trigger_eff_Ph);
      hSR[model+"_gen"].Scale(cfg.trigger_eff_Ph);
      hCR[model].Scale(cfg.trigger_eff_Ph);
      sVar = "pre_ph165/c_MET300/MT300/STg";
      saver_hist.save(hSR[model],sVar+"/"+model);
      saver_hist.save(hSR[model+"_gen"],sVar+"/"+model+"_gen");

      // acceptance
      float const x=masses.first;
      float const y=masses.second;
      grAcc[sScan].SetPoint(grAcc[sScan].GetN(),x,y,float(miAcc[model])/Ngen);
      grScaleUnc[sScan].SetPoint(grScaleUnc[sScan].GetN(),x,y,scaleUnc);

      // contamination
      hCR[model].Scale(1./nData);
      sVar="pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh";
      saver_hist.save(hCR[model],sVar+"/"+model);
      grCont[sScan].SetPoint(grCont[sScan].GetN(),x,y,hCR[model].Integral());
   }
   file.Close();

   for (auto const& p: grAcc) {
      sVar = "pre_ph165/c_MET300/MT300/STg";
      saver_hist.save(grAcc[p.first],sVar+"/"+p.first+"_acceptance");
      saver_hist.save(grScaleUnc[p.first],sVar+"/"+p.first+"_scaleUnc");

      sVar="pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh";
      saver_hist.save(grCont[p.first],sVar+"/"+p.first+"_contamination");
   }
}

extern "C"
void run()
{
	//runScan_80X(TChiWg);
	//runScan_80X(TChiNg);
	runScan_80X(T5Wg);
	//runScan_80X(T5gg);
	//runScan_80X(T6Wg);
	//runScan_80X(T6gg);  
}
