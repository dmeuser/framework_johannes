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

static io::Logger logFile(TString::Format("distributions%.0f.log",cfg.processFraction*100).Data());

void saveHistograms(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH1F> &hs,hist::Histograms<TH1F> &hs_pix,bool saveData)
{
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
//         for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets","T5gg","T5Wg","GGM"}){
         for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR","T5Wg","TChiWG"}){
            saver_hist.save(*hs.getHistogram(sVar,sSample),sVar+"/"+sSample);
         }
         TH1F hEFake(*hs_pix.getHistogram(sVar,"SinglePhoton"));
         hEFake.Scale(cfg.efake.f);
         hEFake.SetLineColor(cfg.efake.color);
         // hEFake.SetFillStyle(1001);
         saver_hist.save(hEFake,sVar+"/efake");
         if (saveData) saver_hist.save(*hs.getHistogram(sVar,"SinglePhoton"),sVar+"/SinglePhoton");
         //if (saveData) saver_hist.save(*hs.getHistogram(sVar,"MET"),sVar+"/MET");        
      }
   }
}

void saveHistograms(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH2F> &hs)
{
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
//         for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets","T5gg","T5Wg","GGM"}){
         for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR","T5Wg","TChiWG"}){
            saver_hist.save(*hs.getHistogram(sVar,sSample),sVar+"/"+sSample);
         }
      }
   }
}

/* sx,sy: signal region borders; cx,xy: control region borders */
//~ static void drawCR_SR(TH2F const &h,float sx,float sy,float cx,float cy){
   //~ TAxis const
      //~ &x=*h.GetXaxis(),
      //~ &y=*h.GetYaxis();
   //~ TLine l;
   //~ l.SetLineWidth(2);
   //~ l.SetLineStyle(1);
   //~ // make a light/dark gray dashed line
   //~ l.SetLineColor(kGray);
   //~ l.DrawLine(sx,sy,sx,y.GetXmax());
   //~ l.DrawLine(sx,sy,x.GetXmax(),sy);
   //~ l.DrawLine(cx,cy,cx,y.GetXmax());
   //~ l.DrawLine(cx,cy,x.GetXmax(),cy);
   //~ // superimpose dark dashed line
   //~ l.SetLineColor(kGray+2);
   //~ l.SetLineStyle(2);
   //~ l.DrawLine(sx,sy,sx,y.GetXmax());
   //~ l.DrawLine(sx,sy,x.GetXmax(),sy);
   //~ l.DrawLine(cx,cy,cx,y.GetXmax());
   //~ l.DrawLine(cx,cy,x.GetXmax(),cy);
   //~ // labels
   //~ TLatex t;
   //~ t.SetTextFont(42);
   //~ t.SetTextSize(.04);
   //~ t.SetTextAlign(31);
   //~ t.SetTextColor(kOrange+9);
   //~ float xrange=x.GetXmax()-x.GetXmin();
   //~ float yrange=y.GetXmax()-y.GetXmin();
   //~ t.DrawLatex(x.GetXmax()-.03*xrange,sy+.01*yrange,"SR");
   //~ t.DrawLatex(x.GetXmax()-.03*xrange,cy+.01*yrange,"CR");
//~ }

/* sx,sy: signal region borders; cx,xy: control region borders */
//~ static void drawCR_SR_VR(TH2F const &h,float sx,float sy,float cx,float cy, bool logx=false){
   //~ TAxis const
      //~ &x=*h.GetXaxis(),
      //~ &y=*h.GetYaxis();
   //~ if (logx) cx=x.GetXmin();
   //~ TLine l;
   //~ l.SetLineWidth(2);
   //~ l.SetLineStyle(1);
   //~ // make a light/dark gray dashed line
   //~ l.SetLineColor(kGray);
   //~ l.DrawLine(sx,sy,sx,y.GetXmax());
   //~ l.DrawLine(cx,sy,x.GetXmax(),sy);
   //~ l.DrawLine(cx,cy,x.GetXmax(),cy);
   //~ if (!logx) l.DrawLine(cx,cy,cx,y.GetXmax());
   //~ // superimpose dark dashed line
   //~ l.SetLineColor(kGray+2);
   //~ l.SetLineStyle(2);
   //~ l.DrawLine(sx,sy,sx,y.GetXmax());
   //~ l.DrawLine(cx,sy,x.GetXmax(),sy);
   //~ l.DrawLine(cx,cy,x.GetXmax(),cy);
   //~ if (!logx) l.DrawLine(cx,cy,cx,y.GetXmax());
   //~ // labels
   //~ TLatex t;
   //~ t.SetTextFont(42);
   //~ t.SetTextSize(.04);
   //~ t.SetTextAlign(31);
   //~ t.SetTextColor(kOrange+9);
   //~ float xrange=x.GetXmax()-x.GetXmin();
   //~ float yrange=y.GetXmax()-y.GetXmin();
   //~ t.DrawLatex(x.GetXmax()-.03*xrange,sy+.01*yrange,"SR");
   //~ t.DrawLatex(x.GetXmax()-.03*xrange,cy+.01*yrange,"CR");
   //~ t.SetTextAlign(13);
   //~ t.DrawLatex(cx,y.GetXmax()-.01*yrange,"VR");
//~ }

// isolation thresholds for Iso40 in HLT
static TF1 fIso40ec("fIso40ec","4.0+0.012*x",0,1000);
static TF1 fIso40hc("fIso40hc","4.0+0.005*x",0,1000);
static TF1 fIso40tr("fIso40tr","4.0+0.002*x",0,1000);

//~ static void drawFunction(TF1 &fun){
   //~ fun.SetLineWidth(2);
   //~ fun.SetLineColor(kBlack);
   //~ fun.DrawCopy("same");
//~ }

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
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/STg_SRbin"   ,";STg;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300});
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/STg"   ,";STg;EventsBIN"           ,300,0,3000);
   
   //Hists for checking overlap
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/exclusive/STg_SRbin"   ,";STg;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300});
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/exclusive/STg"   ,";STg;EventsBIN"           ,300,0,3000);
   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/diphoton/STg_SRbin"   ,";STg;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300});
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/diphoton/STg"   ,";STg;EventsBIN"           ,300,0,3000);
   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/lepton/STg_SRbin"   ,";STg;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300});
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/lepton/STg"   ,";STg;EventsBIN"           ,300,0,3000);
   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/HTG/STg_SRbin"   ,";STg;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300});
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/HTG/STg"   ,";STg;EventsBIN"           ,300,0,3000);
   hs2d.addHist("pre_ph165/c_MET300/MT300/STg600/HTG/STg_lowHTG"   ,";STg;MET;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300}, {350,450,600,800},{100,150,200});
   hs2d.addHist("pre_ph165/c_MET300/MT300/STg600/HTG/STg_highHTG"   ,";STg;MET;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300}, {350,450,600,800},{100,150,200});
   
   //Hists for checking background prediction for combination
   ADD_HIST("pre_ph165/combined/HTG"   ,";%HTG;EventsBIN" ,30,0,2500);
   ADD_HIST("pre_ph165/combined/MET"   ,";%MET;EventsBIN" ,30,0,800);
   ADD_HIST("pre_ph165/combined/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,30,0,1.7);
   
   ADD_HIST("pre_ph165/VR_SR/HTG"   ,";%HTG;EventsBIN" ,30,0,2500);
   ADD_HIST("pre_ph165/VR_SR/MET"   ,";%MET;EventsBIN" ,30,0,800);
   ADD_HIST("pre_ph165/VR_SR/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,30,0,1.7);
   ADD_HIST("pre_ph165/VR_SR/phoPt"   ,";#gamma_{1} PT;EventsBIN" ,30,170,2000);
   ADD_HIST("pre_ph165/VR_SR/MT"   ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN" ,30,300,1000);
   ADD_HIST("pre_ph165/VR_SR/phoEta" ,";#gamma #eta;EventsBIN",52,-2.6,2.6);
   ADD_HIST("pre_ph165/VR_SR/nEle"   ,";N_{e};EventsBIN" ,5,0,5);
   ADD_HIST("pre_ph165/VR_SR/nMu"   ,";N_{#mu};EventsBIN" ,5,0,5);
   ADD_HIST("pre_ph165/VR_SR/nPho"   ,";N_{#gamma mediumId};EventsBIN" ,5,0,5);
   ADD_HIST("pre_ph165/VR_SR/phoR9"   ,";r9 #gamma;EventsBIN" ,20,0,1);
   ADD_HIST("pre_ph165/VR_SR/dRPhoLep"   ,";#DeltaR(#gamma,l_{any});EventsBIN" ,20,0,2);
   ADD_HIST("pre_ph165/VR_SR/dRPhoLepLead"   ,";#DeltaR(#gamma,l_{leading});EventsBIN" ,20,0,2);
   ADD_HIST("pre_ph165/VR_SR/massdiffElePhoToZ"   ,";|M(e,#gamma)-M_{Z}| (GeV);EventsBIN" ,30,0,30);
   
   //Same as above but without overlap to emht analysis
   ADD_HIST("pre_ph165/VR_SR/noHTG/HTG"   ,";%HTG;EventsBIN" ,30,0,2500);
   ADD_HIST("pre_ph165/VR_SR/noHTG/MET"   ,";%MET;EventsBIN" ,30,0,800);
   ADD_HIST("pre_ph165/VR_SR/noHTG/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,30,0,1.7);
   ADD_HIST("pre_ph165/VR_SR/noHTG/phoPt"   ,";#gamma_{1} PT;EventsBIN" ,30,170,2000);
   ADD_HIST("pre_ph165/VR_SR/noHTG/MT"   ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN" ,30,300,1000);
   ADD_HIST("pre_ph165/VR_SR/noHTG/phoEta" ,";#gamma #eta;EventsBIN",52,-2.6,2.6);
   
   //Same as above but without overlap to lepton analysis
   ADD_HIST("pre_ph165/VR_SR/noLepton/HTG"   ,";%HTG;EventsBIN" ,30,0,2500);
   ADD_HIST("pre_ph165/VR_SR/noLepton/MET"   ,";%MET;EventsBIN" ,30,0,800);
   ADD_HIST("pre_ph165/VR_SR/noLepton/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,30,0,1.7);
   ADD_HIST("pre_ph165/VR_SR/noLepton/phoPt"   ,";#gamma_{1} PT;EventsBIN" ,30,170,2000);
   ADD_HIST("pre_ph165/VR_SR/noLepton/MT"   ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN" ,30,300,1000);
   ADD_HIST("pre_ph165/VR_SR/noLepton/phoEta" ,";#gamma #eta;EventsBIN",52,-2.6,2.6);
   
   //Same as above but without overlap to diphoton analysis
   ADD_HIST("pre_ph165/VR_SR/noDiphoton/HTG"   ,";%HTG;EventsBIN" ,30,0,2500);
   ADD_HIST("pre_ph165/VR_SR/noDiphoton/MET"   ,";%MET;EventsBIN" ,30,0,800);
   ADD_HIST("pre_ph165/VR_SR/noDiphoton/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,30,0,1.7);
   ADD_HIST("pre_ph165/VR_SR/noDiphoton/phoPt"   ,";#gamma_{1} PT;EventsBIN" ,30,170,2000);
   ADD_HIST("pre_ph165/VR_SR/noDiphoton/MT"   ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN" ,30,300,1000);
   ADD_HIST("pre_ph165/VR_SR/noDiphoton/phoEta" ,";#gamma #eta;EventsBIN",52,-2.6,2.6);
   
   //Same as above but without any overlap to analyses included in combination
   ADD_HIST("pre_ph165/VR_SR/exclusiv/HTG"   ,";%HTG;EventsBIN" ,30,0,2500);
   ADD_HIST("pre_ph165/VR_SR/exclusiv/MET"   ,";%MET;EventsBIN" ,30,0,800);
   ADD_HIST("pre_ph165/VR_SR/exclusiv/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,30,0,1.7);
   ADD_HIST("pre_ph165/VR_SR/exclusiv/phoPt"   ,";#gamma_{1} PT;EventsBIN" ,30,170,2000);
   ADD_HIST("pre_ph165/VR_SR/exclusiv/MT"   ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN" ,30,300,1000);
   ADD_HIST("pre_ph165/VR_SR/exclusiv/phoEta" ,";#gamma #eta;EventsBIN",52,-2.6,2.6);
   
   //Same as above but without any overlap to analyses included in combination only for VR
   ADD_HIST("pre_ph165/VR/exclusiv/HTG"   ,";%HTG;EventsBIN" ,50,0,2500);
   ADD_HIST("pre_ph165/VR/exclusiv/MET"   ,";%MET;EventsBIN" ,50,0,800);
   ADD_HIST("pre_ph165/VR/exclusiv/absdPhi_pmMet_Pho"   ,";min(|#Delta#phi(#pm MET,#gamma_{1})|);EventsBIN" ,50,0,1.7);
   ADD_HIST("pre_ph165/VR/exclusiv/phoPt"   ,";#gamma_{1} PT;EventsBIN" ,50,150,700);
   ADD_HIST("pre_ph165/VR/exclusiv/MT"   ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN" ,50,300,1000);
   ADD_HIST("pre_ph165/VR/exclusiv/phoEta" ,";#gamma #eta;EventsBIN",52,-2.6,2.6);
   ADD_HIST("pre_ph165/VR/exclusiv/STG"   ,";STg;EventsBIN" ,50,400,800);
   ADD_HIST("pre_ph165/VR/exclusiv/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   
   //For Limits and syst errors
   ADD_HIST("pre_ph165/c_MET300/MT300/exclusiv/STg"   ,";STg;EventsBIN"           ,2000,0,2000);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/exclusiv/absphiMETnJetPh",";|#Delta#phi(p_{T}^{miss},nearest jet/#gamma)| (radians);EventsBIN",50,0,5);  
   
   //Files for eventnumbers of overlapping regions for combination
   std::string outdir (CMAKE_SOURCE_DIR);
   outdir = outdir + "output/events_overlap/";
   std::ofstream exclusive;
   std::ofstream htg;
   std::ofstream lepton;
   std::ofstream diphoton;
   std::ofstream total;
   
   exclusive.open(outdir + TString::Format("exclusive%.1f.txt",cfg.processFraction*100));
   htg.open(outdir + TString::Format("htg%.1f.txt",cfg.processFraction*100));
   lepton.open(outdir + TString::Format("lepton%.1f.txt",cfg.processFraction*100));
   diphoton.open(outdir + TString::Format("diphoton%.1f.txt",cfg.processFraction*100));
   total.open(outdir + TString::Format("total%.1f.txt",cfg.processFraction*100));


 //  TH1F *trigger_photon_pt_n = new TH1F("numerator);

   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      hs_notPix.setCurrentSample(dss.name);
      hs_pix.setCurrentSample(dss.name);
      hs2d.setCurrentSample(dss.name);

      bool const isData=dss.isData;
      bool const isSignal=dss.isSignal;

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

      float phoPt;

      std::shared_ptr<WeightCalculator> wCalc;
      if (dss.datasetName=="ZNuNuGJets" || dss.datasetName=="ZGTo2LG"){
         wCalc=std::make_shared<WeightCalculatorZgamma>(phoPt);
      }
      else if (dss.datasetName=="WGToLNuG"){
         wCalc=std::make_shared<WeightCalculatorWgamma>();
      }

      int iEv=0;
      //~ int events =0;
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
         hs2d.setFillWeight(fEventWeight);

         std::vector<tree::Photon const *> lPho,mPho,tPho,lPho15,mPho15,tPho15,lPixPho;
         for (tree::Photon const &ph: *photons){

            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
    //        if (ph.r9 < 0.9 || ph.r9 > 1.0) continue;
            if (fabs(ph.p.Eta())>1.4442) continue;
            if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
            if (ph.hasPixelSeed){
               lPixPho.push_back(&ph);
            } else {
               //Attention check lossePhoton 15 or new
               if (ph.isLoose) lPho.push_back(&ph);
               if (ph.isMedium) mPho.push_back(&ph);
               if (ph.isTight)  tPho.push_back(&ph);
               if (ph.isLoose15) lPho15.push_back(&ph);
               if (ph.isMedium15) mPho15.push_back(&ph);
               if (ph.isTight15)  tPho15.push_back(&ph);                        
            }
         }




         // independent of "normal"/"pixel" run
         // jet related
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         //~ int const Njets=cjets.size();
         int NbL=0, NbM=0, NbT=0;
         // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
         for (auto const &j: cjets) {
            if (j.bDiscriminator > 0.5426) NbL++;
            if (j.bDiscriminator > 0.8484) NbM++;
            if (j.bDiscriminator > 0.9535) NbT++;
         }
         //~ float const HT=phys::computeHT(cjets);
         
         TVector3 vecHT, vecGammaHT;
         bool clean_MET = true;
         
         for (auto const &jet: cjets) {
            vecHT+=jet.p;
            if (jet.p.Pt() < 100) continue;            
            if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
         }

         if (!clean_MET) continue;
    
         //~ float const dPhiMETjet=Njets>0 ? MET->p.DeltaPhi(cjets[0].p) : 4;
         //~ float const METdotJet=Njets>0 ? MET->p.Dot(cjets[0].p) : 0;

         float dPhiMETnearJet=4;
         float dPhiMETnear2Jet=4; // only considering two leading jets
         float iJet=0;
         for (auto const &j: cjets){
            iJet++;
            const float dPhi=MET->p.DeltaPhi(j.p);
            if (std::abs(dPhi) < std::abs(dPhiMETnearJet))
               dPhiMETnearJet=dPhi;
            if (iJet > 2) continue;
            if (std::abs(dPhi) < std::abs(dPhiMETnear2Jet))
               dPhiMETnear2Jet=dPhi;
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
         
         // lepton related
         int Nele=0;
         for (tree::Electron const &ele: *electrons) {
     //       if (ele.p.Pt()<20) continue;
            if (ele.p.Eta()>2.5) continue;
            if (ele.p.Eta()>1.4442 && ele.p.Eta()<1.6) continue;
            if (ele.isTight) Nele++;
         }
         int Nmu=0;
         for (tree::Muon const &mu: *muons) {
    //        if (mu.p.Pt()<20) continue;
            if (mu.p.Eta()>2.4) continue;
            if (mu.isTight && mu.rIso<0.15) Nmu++;
         }
         //~ const int Nlepton=Nele+Nmu;

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

            float relPt2Jets = 0;
            float DeltaS = 0;
            float DeltaS1 = 0;

            if (cjets.size() == 1 ) {
               TVector3 vPhotonMet = pho[0]->p + MET->p;
               DeltaS1 = std::acos((vPhotonMet.Pt()*cjets[0].p.Pt()*std::cos(vPhotonMet.DeltaPhi(cjets[0].p)))/(vPhotonMet.Pt()*cjets[0].p.Pt()));
            }
            if (cjets.size() > 1 ) {
               TVector3 vleading2jets = cjets[0].p + cjets[1].p;
               TVector3 vPhotonMet = pho[0]->p + MET->p;
               relPt2Jets = vleading2jets.Pt() / (cjets[0].p.Pt() + cjets[1].p.Pt());
               DeltaS = std::acos((vPhotonMet.Pt()*vleading2jets.Pt()*std::cos(vPhotonMet.DeltaPhi(vleading2jets)))/(vPhotonMet.Pt()*vleading2jets.Pt()));
               DeltaS1 = std::acos((vPhotonMet.Pt()*vleading2jets.Pt()*std::cos(vPhotonMet.DeltaPhi(vleading2jets)))/(vPhotonMet.Pt()*vleading2jets.Pt()));              
            }
            

            fEventWeight=*w_pu * *w_mc;
            if (wCalc) fEventWeight*=wCalc->get();
            // [1-8] are muF, muR weights
            // [4] is x2, x2; [8] is x0.5, x0.5
           // if (w_pdf->size()>=9) fEventWeight*=(*w_pdf)[1];
     //       if (w_pdf->size()>=1) debug << (*w_pdf)[5];
            
            hs_notPix.setFillWeight(fEventWeight);
            hs_pix.setFillWeight(fEventWeight);
            hs2d.setFillWeight(fEventWeight);

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

               int const Nph=pho.size();                  
               //~ float const METSHT= phys::METoverSqrtHT(MET->p.Pt(),HT);
               float const MT=phys::M_T(*pho[0],*MET);
   
               float const dPhiMETph=MET->p.DeltaPhi(pho[0]->p);
               //~ float const dPhiJetPh=Njets>0 ? pho[0]->p.DeltaPhi(cjets[0].p) : 4;
               //~ float const dPhiMETHT=Njets>0 ? MET->p.DeltaPhi(vecHT) : 4;
               //~ float const dPhiPhHT=Njets>0 ? pho[0]->p.DeltaPhi(vecHT) : 4;                            
               //~ float const METdotPh=MET->p.Dot(pho[0]->p);
   
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
                  vecGammaHT = vecHT + ph->p;
               }   

            bool triggerMatch=false;
            if (isData) {
               for (tree::Particle const& to: *triggerObjects){
                  if (to.p.DeltaR(pho[0]->p)<0.1 && std::fabs(to.p.Pt()-pho[0]->p.Pt())/pho[0]->p.Pt()<0.7) {
                     triggerMatch=true;
                  }
               }
            } 
            
            
            ////Vetos for combination
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
            if (pho[0]->isMedium && pho[0]->p.Pt() > 40 && Nph > 1 && met > 100) {
               if(pho[1]->isMedium && pho[1]->p.Pt() > 40 && phys::invmass(*pho[0],*pho[1]) > 105 && pho[0]->p.DeltaR(pho[1]->p) > 0.3) {
                  diphotonVeto = true;
               }
            }
            
            /////////////
            //EMHT Veto//
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
            if (emht > 700 && met > 350 && fabs(dPhiMETph) > 0.3 && fabs(fabs(dPhiMETph)-TMath::Pi()) > 0.3 && phoPt > 100)  {
               emhtVeto = true;
            }
            
            
            if (phoPt>180 && (!isData || (*trigger_Ph && triggerMatch))){

               if (met > 300 && MT > 300){

                  if (STg>600) {
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/STg_SRbin",STg);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/STg",STg);
                     if (isData && pass==pass_normal) {total << *runNo << ":" << *lumNo << ":" << *evtNo << std::endl;}
                     
                     if (pho[0]->isMedium && Nph > 1) {
                        if(pho[1]->isMedium && pho[1]->p.Pt() > 40 && phys::invmass(*pho[0],*pho[1]) > 105 && pho[0]->p.DeltaR(pho[1]->p) > 0.3) {
                           hs.fill("pre_ph165/c_MET300/MT300/STg600/diphoton/STg_SRbin",STg);
                           hs.fill("pre_ph165/c_MET300/MT300/STg600/diphoton/STg",STg);
                           diphotonVeto = true;
                           if (isData && pass==pass_normal) {diphoton << *runNo << ":" << *lumNo << ":" << *evtNo << std::endl;}
                        }
                     }
                     
                     if (leptoVeto == 1) {
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/lepton/STg_SRbin",STg);
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/lepton/STg",STg);
                        if (isData && pass==pass_normal) {
                           lepton << *runNo << ":" << *lumNo << ":" << *evtNo << std::endl;
                           }
                     }
                     
                     if (emht > 700 && met > 350 && fabs(dPhiMETph) > 0.3 && fabs(fabs(dPhiMETph)-TMath::Pi()) > 0.3)  {
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/HTG/STg_SRbin",STg);
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/HTG/STg",STg);
                        if (emht < 2000 && pass==pass_normal){
                           hs2d.fill("pre_ph165/c_MET300/MT300/STg600/HTG/STg_lowHTG",STg,met);
                        }
                        else if(pass==pass_normal) {
                           hs2d.fill("pre_ph165/c_MET300/MT300/STg600/HTG/STg_highHTG",STg,met);
                        }
                        emhtVeto = true;
                        if (isData && pass==pass_normal) { htg << *runNo << ":" << *lumNo << ":" << *evtNo << std::endl;}
                     }
                     
                     if (emhtVeto == 0 && leptoVeto == 0 && diphotonVeto == 0) {
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/exclusive/STg_SRbin",STg);
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/exclusive/STg",STg);
                        if (isData && pass==pass_normal) {exclusive << *runNo << ":" << *lumNo << ":" << *evtNo << std::endl;}
                     }
                  }
                  if (emhtVeto == 0 && leptoVeto == 0 && diphotonVeto == 0) {
                     hs.fill("pre_ph165/c_MET300/MT300/exclusiv/STg",STg);
                  }
               }
               if (MT > 100 && met > 100 && emhtVeto == 0 && leptoVeto == 0 && diphotonVeto == 0){
                  if (met < 300 || MT < 300){
                     hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/exclusiv/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                  }
               }
                     
               
               if (met > 100 && MT > 100) {
                  //Fill hists for background check in cut variables
                  hs.fill("pre_ph165/combined/HTG",emht);
                  hs.fill("pre_ph165/combined/MET",met);
                  hs.fill("pre_ph165/combined/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                  
                  //same but only for SR and VR
                  if (met > 300 && MT > 300) {
                     hs.fill("pre_ph165/VR_SR/HTG",emht);
                     hs.fill("pre_ph165/VR_SR/MET",met);
                     hs.fill("pre_ph165/VR_SR/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                     hs.fill("pre_ph165/VR_SR/nEle",el_comb.size());
                     hs.fill("pre_ph165/VR_SR/nMu",mu_comb.size());
                     hs.fill("pre_ph165/VR_SR/nPho",mPho.size());
                     hs.fill("pre_ph165/VR_SR/phoPt",phoPt);
                     hs.fill("pre_ph165/VR_SR/MT",MT);
                     hs.fill("pre_ph165/VR_SR/phoEta",pho[0]->p.Eta());
                     hs.fill("pre_ph165/VR_SR/phoR9",pho[0]->r9);
                     for (auto const &ele: el_comb) {
                        hs.fill("pre_ph165/VR_SR/dRPhoLep",ele->p.DeltaR(pho[0]->p));
                     }
                     for (auto const &mu: mu_comb) {
                        hs.fill("pre_ph165/VR_SR/dRPhoLep",mu->p.DeltaR(pho[0]->p));
                     }
                     if (leptoVeto == true) hs.fill("pre_ph165/VR_SR/dRPhoLepLead",leadLep->p.DeltaR(pho[0]->p));
                     if (leadLep_ele == true) {
                        hs.fill("pre_ph165/VR_SR/massdiffElePhoToZ",fabs(invmassPhoLep-91.1876));
                     }
                     
                     //without overlap to emht analysis
                     if (emhtVeto == false) {
                        hs.fill("pre_ph165/VR_SR/noHTG/HTG",emht);
                        hs.fill("pre_ph165/VR_SR/noHTG/MET",met);
                        hs.fill("pre_ph165/VR_SR/noHTG/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                        hs.fill("pre_ph165/VR_SR/noHTG/phoPt",phoPt);
                        hs.fill("pre_ph165/VR_SR/noHTG/MT",MT);
                        hs.fill("pre_ph165/VR_SR/noHTG/phoEta",pho[0]->p.Eta());
                     }
                     
                     //without overlap to lepton analysis
                     if (leptoVeto == false) {
                        hs.fill("pre_ph165/VR_SR/noLepton/HTG",emht);
                        hs.fill("pre_ph165/VR_SR/noLepton/MET",met);
                        hs.fill("pre_ph165/VR_SR/noLepton/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                        hs.fill("pre_ph165/VR_SR/noLepton/phoPt",phoPt);
                        hs.fill("pre_ph165/VR_SR/noLepton/MT",MT);
                        hs.fill("pre_ph165/VR_SR/noLepton/phoEta",pho[0]->p.Eta());
                     }
                     
                     //without overlap to diphoton analysis
                     if (diphotonVeto == false) {
                        hs.fill("pre_ph165/VR_SR/noDiphoton/HTG",emht);
                        hs.fill("pre_ph165/VR_SR/noDiphoton/MET",met);
                        hs.fill("pre_ph165/VR_SR/noDiphoton/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                        hs.fill("pre_ph165/VR_SR/noDiphoton/phoPt",phoPt);
                        hs.fill("pre_ph165/VR_SR/noDiphoton/MT",MT);
                        hs.fill("pre_ph165/VR_SR/noDiphoton/phoEta",pho[0]->p.Eta());
                     }
                     
                     //without overlap to any analysis included in the combination
                     if (leptoVeto == false && diphotonVeto == false && emhtVeto == false) {
                        hs.fill("pre_ph165/VR_SR/exclusiv/HTG",emht);
                        hs.fill("pre_ph165/VR_SR/exclusiv/MET",met);
                        hs.fill("pre_ph165/VR_SR/exclusiv/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                        hs.fill("pre_ph165/VR_SR/exclusiv/phoPt",phoPt);
                        hs.fill("pre_ph165/VR_SR/exclusiv/MT",MT);
                        hs.fill("pre_ph165/VR_SR/exclusiv/phoEta",pho[0]->p.Eta());
                     }
                     
                     //without overlap to any analysis included in the combination only VR
                     if (leptoVeto == false && diphotonVeto == false && emhtVeto == false && STg < 600) {
                        hs.fill("pre_ph165/VR/exclusiv/HTG",emht);
                        hs.fill("pre_ph165/VR/exclusiv/MET",met);
                        hs.fill("pre_ph165/VR/exclusiv/absdPhi_pmMet_Pho",TMath::Min(fabs(dPhiMETph),fabs(fabs(dPhiMETph)-TMath::Pi())));
                        hs.fill("pre_ph165/VR/exclusiv/phoPt",phoPt);
                        hs.fill("pre_ph165/VR/exclusiv/MT",MT);
                        hs.fill("pre_ph165/VR/exclusiv/phoEta",pho[0]->p.Eta());
                        hs.fill("pre_ph165/VR/exclusiv/STG",STg);
                        hs.fill("pre_ph165/VR/exclusiv/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
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
      hs2d.scaleLumi();
      hs2d.mergeOverflow();
      file.Close();
   } // datasets
   
   //closing files for evtnumbers
   exclusive.close();
   diphoton.close();
   lepton.close();
   htg.close();
   total.close();

   // calling the "normal" histogram "hs" from here, since it's the most used
   hist::Histograms<TH1F> &hs = hs_notPix;
//   std::vector<TString> samplesToCombine={"GJets","QCD","ZNuNuGJets","ZGTo2LG","WGToLNuG","ZNuNuJets","WLNuJets","diboson","T5gg","T5Wg","GGM","SinglePhoton","MET"};
   //std::vector<TString> samplesToCombine={"GJets_DR","QCD","ZNuNuGJets","ZGTo2LG","WGToLNuG","ZNuNuJets","WLNuJets","diboson","T5Wg","TChiWG","SinglePhoton","MET"};
   std::vector<TString> samplesToCombine={"GJets_DR","QCD","ZNuNuGJets","ZGTo2LG","WGToLNuG","ZNuNuJets","WLNuJets","diboson","T5Wg","TChiWG","SinglePhoton"};
   hs    .combineFromSubsamples(samplesToCombine);
   hs_pix.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   io::RootFileSaver saver("plots.root",TString::Format("danilo_distributions%.1f",cfg.processFraction*100));
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("danilo_distributions%.1f",cfg.processFraction*100),false);
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
   //check overlap
      {"pre_ph165/c_MET300/MT300/STg600/",{"STg_SRbin","STg",}},          
      {"pre_ph165/c_MET300/MT300/STg600/exclusive/",{"STg_SRbin","STg",}},
      {"pre_ph165/c_MET300/MT300/STg600/diphoton/",{"STg_SRbin","STg",}},
      {"pre_ph165/c_MET300/MT300/STg600/lepton/",{"STg_SRbin","STg",}},
      {"pre_ph165/c_MET300/MT300/STg600/HTG/",{"STg_SRbin","STg",}},
   //check background prediction for combination
      {"pre_ph165/combined/",{"HTG","MET","absdPhi_pmMet_Pho"}},
      {"pre_ph165/VR_SR/",{"HTG","MET","absdPhi_pmMet_Pho","nEle","nMu","nPho","phoPt","MT","phoEta","phoR9","dRPhoLep","dRPhoLepLead","massdiffElePhoToZ"}},
      {"pre_ph165/VR_SR/noHTG/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR_SR/noLepton/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR_SR/noDiphoton/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR_SR/exclusiv/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR/exclusiv/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta","STG","absphiMETnJetPh"}},
   //for final datacards
      {"pre_ph165/c_MET300/MT300/exclusiv/",{"STg"}},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/exclusiv/",{"absphiMETnJetPh"}},
   };
   
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         THStack st_mc=hs.getStack(sVar,{"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR"});
         gfx::LegendEntries le=hs.getLegendEntries();
         TH1F hEFake(*hs_pix.getHistogram(sVar,"SinglePhoton"));
         hEFake.Scale(cfg.efake.f);
         hEFake.SetFillColor(cfg.efake.color);
         hEFake.SetFillStyle(1001);
         THStack st_mc_efake=hist::stackPrepend(st_mc,hEFake,"hist");
         THStack &st = sVar.Contains(sPresel+"Ng")
            ? st_mc // do not draw e->gamma fake (not filled anyway)
            : st_mc_efake;
         if (!sVar.Contains(sPresel+"Ng")) le.append(hEFake,cfg.efake.label,"f");
         st.Draw();
//         auto hists=hs.getHistograms(sVar,{"T5gg","T5Wg","GGM"});
         auto hists=hs.getHistograms(sVar,{"T5Wg","TChiWG"});
         for (auto const &h: hists) h->Draw("same hist");
         le+=hs.getLegendEntries();
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         saver.save(can,sVar);
      }
   }
   //SAVING HISTOGRAMS
   msPresel_vVars={
      //check overlap
      {"pre_ph165/c_MET300/MT300/STg600/",{"STg_SRbin","STg",}},          
      {"pre_ph165/c_MET300/MT300/STg600/exclusive/",{"STg_SRbin","STg",}},
      {"pre_ph165/c_MET300/MT300/STg600/diphoton/",{"STg_SRbin","STg",}},
      {"pre_ph165/c_MET300/MT300/STg600/lepton/",{"STg_SRbin","STg",}},
      {"pre_ph165/c_MET300/MT300/STg600/HTG/",{"STg_SRbin","STg",}},
      //check background prediction for combination
      {"pre_ph165/combined/",{"HTG","MET","absdPhi_pmMet_Pho"}},
      {"pre_ph165/VR_SR/",{"HTG","MET","absdPhi_pmMet_Pho","nEle","nMu","nPho","phoPt","MT","phoEta","phoR9","dRPhoLep","dRPhoLepLead","massdiffElePhoToZ"}},
      {"pre_ph165/VR_SR/noHTG/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR_SR/noLepton/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR_SR/noDiphoton/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR_SR/exclusiv/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta"}},
      {"pre_ph165/VR/exclusiv/",{"HTG","MET","absdPhi_pmMet_Pho","phoPt","MT","phoEta","STG","absphiMETnJetPh"}},
      //for final datacards
      {"pre_ph165/c_MET300/MT300/exclusiv/",{"STg"}},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/exclusiv/",{"absphiMETnJetPh"}},
   };
   saveHistograms(msPresel_vVars,saver_hist,hs,hs_pix,true);
   
   
   //2d Plot to present overlap for combination
   TCanvas can_2d;
   TH2F overlap("Overlap_2d","",4,0.5,4.5,4,0.5,4.5);
   TH2F overlap_ratio("Overlap_ratio_2d","",4,0.5,4.5,4,0.5,4.5);
   
   int i = 1;
   TH1F *histCompl = hs.getHistogram("pre_ph165/c_MET300/MT300/STg600/STg_SRbin","SinglePhoton");
   for (TString analys: {"exclusive","HTG","lepton","diphoton"}) {
      TH1F *histTemp = hs.getHistogram("pre_ph165/c_MET300/MT300/STg600/"+analys+"/STg_SRbin","SinglePhoton");
      for (int j = 1; j <= 4; j++) {
         overlap.SetBinContent(j,i,histTemp->GetBinContent(j));
         overlap_ratio.SetBinContent(j,i,roundf(histTemp->GetBinContent(j)/(1.0*histCompl->GetBinContent(j))*100)/100);
      }
      overlap.GetYaxis()->SetBinLabel(i,analys);
      overlap_ratio.GetYaxis()->SetBinLabel(i,analys);
      i++;
   }
   can_2d.cd();
   overlap.SetMarkerColor(kRed);
   overlap.GetYaxis()->SetTitleOffset(2);
   overlap.GetXaxis()->SetTickLength(0);
   overlap.GetYaxis()->SetTickLength(0);
   overlap.SetTitle(";SignalBin Photon+ST;Analysis");
   overlap.SetStats(false);
   overlap.Draw("col TEXT");
   saver.save(can_2d,"Overlap_2d",true,false);
   
   can_2d.cd();
   overlap_ratio.SetMarkerColor(kRed);
   overlap_ratio.GetYaxis()->SetTitleOffset(2);
   overlap_ratio.GetXaxis()->SetTickLength(0);
   overlap_ratio.GetYaxis()->SetTickLength(0);
   overlap_ratio.SetTitle(";SignalBin Photon+ST;Analysis");
   overlap_ratio.SetStats(false);
   overlap_ratio.Draw("col TEXT");
   saver.save(can_2d,"overlap_ratio_2d",true,false);
   
   
   //2d Plot to present overlap for combination with all signalregion binning
   TH2F overlap_HTG("Overlap_2d_HTGbinned","",4,0.5,4.5,6,0.5,6.5);
   TH2F overlap_HTG_ratio("Overlap_2d_HTGbinned_ratio","",4,0.5,4.5,6,0.5,6.5);
   
   std::cout<<"okay"<<std::endl;
   
   TH2F *temp_low = hs2d.getHistogram("pre_ph165/c_MET300/MT300/STg600/HTG/STg_lowHTG","SinglePhoton");
   TH2F *temp_high = hs2d.getHistogram("pre_ph165/c_MET300/MT300/STg600/HTG/STg_highHTG","SinglePhoton");
   
   for (int i=1; i<=4; i++){
      for(int j=1; j<=3; j++){
         overlap_HTG.SetBinContent(i,j,temp_low->GetBinContent(i,j));
         overlap_HTG.SetBinContent(i,j+3,temp_high->GetBinContent(i,j));
         TString bin1 = "lowHTG_MET"+std::to_string(j);
         TString bin2 = "highHTG_MET"+std::to_string(j);
         overlap_HTG.GetYaxis()->SetBinLabel(j,bin1);
         overlap_HTG.GetYaxis()->SetBinLabel(j+3,bin2);
         
         overlap_HTG_ratio.SetBinContent(i,j,roundf(temp_low->GetBinContent(i,j)/(1.0*histCompl->GetBinContent(i))*100)/100);
         overlap_HTG_ratio.SetBinContent(i,j+3,roundf(temp_high->GetBinContent(i,j)/(1.0*histCompl->GetBinContent(i))*100)/100);
         overlap_HTG_ratio.GetYaxis()->SetBinLabel(j,bin1);
         overlap_HTG_ratio.GetYaxis()->SetBinLabel(j+3,bin2);
      }
   }

   can_2d.cd();
   overlap_HTG.SetMarkerColor(kRed);
   overlap_HTG.GetXaxis()->SetTitleOffset(2);
   overlap_HTG.GetXaxis()->SetTickLength(0);
   overlap_HTG.GetXaxis()->SetTitleSize(0.03);
   overlap_HTG.GetYaxis()->SetTickLength(0);
   overlap_HTG.GetYaxis()->SetLabelSize(0.03);
   overlap_HTG.GetYaxis()->SetTitleSize(0.03);
   overlap_HTG.GetYaxis()->SetTitleOffset(4);
   overlap_HTG.SetTitle(";SignalBin Photon+ST;SignalBin Photon+HTG");
   overlap_HTG.SetStats(false);
   gPad->SetLeftMargin(0.22);
   overlap_HTG.Draw("col TEXT");
   saver.save(can_2d,"Overlap_2d_HTGbinned",true,false);
   
   can_2d.cd();
   overlap_HTG_ratio.SetMarkerColor(kRed);
   overlap_HTG_ratio.GetXaxis()->SetTitleOffset(2);
   overlap_HTG_ratio.GetXaxis()->SetTickLength(0);
   overlap_HTG_ratio.GetXaxis()->SetTitleSize(0.03);
   overlap_HTG_ratio.GetYaxis()->SetTickLength(0);
   overlap_HTG_ratio.GetYaxis()->SetLabelSize(0.03);
   overlap_HTG_ratio.GetYaxis()->SetTitleSize(0.03);
   overlap_HTG_ratio.GetYaxis()->SetTitleOffset(4);
   overlap_HTG_ratio.SetTitle(";SignalBin Photon+ST;SignalBin Photon+HTG");
   overlap_HTG_ratio.SetStats(false);
   gPad->SetLeftMargin(0.22);
   overlap_HTG_ratio.Draw("col TEXT");
   saver.save(can_2d,"Overlap_2d_HTGbinned_ratio",true,false);
}
