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
         if (saveData) saver_hist.save(*hs.getHistogram(sVar,"MET"),sVar+"/MET");        
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
static void drawCR_SR(TH2F const &h,float sx,float sy,float cx,float cy){
   TAxis const
      &x=*h.GetXaxis(),
      &y=*h.GetYaxis();
   TLine l;
   l.SetLineWidth(2);
   l.SetLineStyle(1);
   // make a light/dark gray dashed line
   l.SetLineColor(kGray);
   l.DrawLine(sx,sy,sx,y.GetXmax());
   l.DrawLine(sx,sy,x.GetXmax(),sy);
   l.DrawLine(cx,cy,cx,y.GetXmax());
   l.DrawLine(cx,cy,x.GetXmax(),cy);
   // superimpose dark dashed line
   l.SetLineColor(kGray+2);
   l.SetLineStyle(2);
   l.DrawLine(sx,sy,sx,y.GetXmax());
   l.DrawLine(sx,sy,x.GetXmax(),sy);
   l.DrawLine(cx,cy,cx,y.GetXmax());
   l.DrawLine(cx,cy,x.GetXmax(),cy);
   // labels
   TLatex t;
   t.SetTextFont(42);
   t.SetTextSize(.04);
   t.SetTextAlign(31);
   t.SetTextColor(kOrange+9);
   float xrange=x.GetXmax()-x.GetXmin();
   float yrange=y.GetXmax()-y.GetXmin();
   t.DrawLatex(x.GetXmax()-.03*xrange,sy+.01*yrange,"SR");
   t.DrawLatex(x.GetXmax()-.03*xrange,cy+.01*yrange,"CR");
}

/* sx,sy: signal region borders; cx,xy: control region borders */
static void drawCR_SR_VR(TH2F const &h,float sx,float sy,float cx,float cy, bool logx=false){
   TAxis const
      &x=*h.GetXaxis(),
      &y=*h.GetYaxis();
   if (logx) cx=x.GetXmin();
   TLine l;
   l.SetLineWidth(2);
   l.SetLineStyle(1);
   // make a light/dark gray dashed line
   l.SetLineColor(kGray);
   l.DrawLine(sx,sy,sx,y.GetXmax());
   l.DrawLine(cx,sy,x.GetXmax(),sy);
   l.DrawLine(cx,cy,x.GetXmax(),cy);
   if (!logx) l.DrawLine(cx,cy,cx,y.GetXmax());
   // superimpose dark dashed line
   l.SetLineColor(kGray+2);
   l.SetLineStyle(2);
   l.DrawLine(sx,sy,sx,y.GetXmax());
   l.DrawLine(cx,sy,x.GetXmax(),sy);
   l.DrawLine(cx,cy,x.GetXmax(),cy);
   if (!logx) l.DrawLine(cx,cy,cx,y.GetXmax());
   // labels
   TLatex t;
   t.SetTextFont(42);
   t.SetTextSize(.04);
   t.SetTextAlign(31);
   t.SetTextColor(kOrange+9);
   float xrange=x.GetXmax()-x.GetXmin();
   float yrange=y.GetXmax()-y.GetXmin();
   t.DrawLatex(x.GetXmax()-.03*xrange,sy+.01*yrange,"SR");
   t.DrawLatex(x.GetXmax()-.03*xrange,cy+.01*yrange,"CR");
   t.SetTextAlign(13);
   t.DrawLatex(cx,y.GetXmax()-.01*yrange,"VR");
}

// isolation thresholds for Iso40 in HLT
static TF1 fIso40ec("fIso40ec","4.0+0.012*x",0,1000);
static TF1 fIso40hc("fIso40hc","4.0+0.005*x",0,1000);
static TF1 fIso40tr("fIso40tr","4.0+0.002*x",0,1000);

static void drawFunction(TF1 &fun){
   fun.SetLineWidth(2);
   fun.SetLineColor(kBlack);
   fun.DrawCopy("same");
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
   
   ADD_HIST("effi/trigger_effi_n"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_d"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_n_cut"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_d_cut"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_n_match"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_d_match"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_n_cut_match"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("effi/trigger_effi_d_cut_match"  ,";#gamma PT;EntriesBIN"    ,50,0,1000); 
    
   ADD_HIST("pre/minDR" ,";min#DeltaR(#gamma_{1},jet);Entries / bin",{0,.5,3,5},{.1,.5,1});
   ADD_HIST("pre/phPt"  ,";#gamma PT;EntriesBIN"    ,50,0,1000);
   ADD_HIST("pre/ph1Pt" ,";#gamma_{1} PT;Events / bin" ,{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre/phEta" ,";#gamma #eta;EventsBIN"   ,52,-2.6,2.6);
   ADD_HIST("pre/jetPt" ,";jet PT;EntriesBIN"       ,50,0,1500);
   ADD_HIST("pre/jet1Pt",";jet_{1} PT;EventsBIN"    ,50,0,1500);
   ADD_HIST("pre/HT"    ,";HT;EventsBIN"            ,50,0,2500);
   ADD_HIST("pre/MET"   ,";MET;Events / bin"        ,{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre/METS"  ,";METSIG;Events / bin"     ,{0,100,250,400,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre/METS_z",";METSIG;Events / bin"        ,{0,100,300},{5,10});
   ADD_HIST("pre/HTgen" ,";HTG;EventsBIN"           ,500,0,2500);
   ADD_HIST("pre/Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre/Ngl15"   ,";N_{#gamma loose,15};EventsBIN" ,5 ,-.5, 5-.5);  
   ADD_HIST("pre/Ngm"   ,";N_{#gamma medium};EventsBIN",5 ,-.5, 5-.5);
   ADD_HIST("pre/Ngt"   ,";N_{#gamma tight};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre/Ngpix" ,";N_{#gamma pix};EventsBIN"   ,5 ,-.5, 5-.5);
   ADD_HIST("pre/Njets" ,";N_{jets};EventsBIN"         ,15,-.5,15-.5);
   ADD_HIST("pre/MT"    ,";MT(#gamma_{1},%MET);EventsBIN"   ,50,0,1500);
   ADD_HIST("pre/phiMETph" ,";#Delta#phi(%MET,#gamma_{1});EventsBIN",50,-3.2,3.2);
   ADD_HIST("pre/phiMETjet",";#Delta#phi(%MET,jet_{1});EventsBIN"   ,50,-3.2,3.2);

   // PHOTON165
   ADD_HIST("pre_ph165/Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre_ph165/Ngl15"   ,";N_{#gamma loose,15};EventsBIN" ,5 ,-.5, 5-.5); 
   ADD_HIST("pre_ph165/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/METS",";METSIG;Events / bin",{0,100,250,400,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/METS_z",";METSIG;Events / bin",{0,100,300},{5,10});
   ADD_HIST("pre_ph165/METS_long",";METSIG;Events / bin",{0,100,250,400,1000,3000,3100},{20,30,50,100,500,100});
   ADD_HIST("pre_ph165/METS_logx",";METSIG;Events / bin",{0,100,250,400,1000,4000,3500},{10,30,50,100,500,500});
   ADD_HIST("pre_ph165/METS_fine",";METSIG;Events / bin",100,0,5000);
   ADD_HIST("pre_ph165/HT"    ,";HT;EventsBIN"            ,50,0,2500);
   ADD_HIST("pre_ph165/METSHT",";METSHT;EventsBIN"        ,100,0,500);
   ADD_HIST("pre_ph165/MT",";MT(#gamma_{1},%MET);EventsBIN",50,0,2000);
   ADD_HIST("pre_ph165/relPt2Jets",";relative Pt of 2 lead. jets;EventsBIN",25,0,1);
   ADD_HIST("pre_ph165/DeltaS",";Delta S;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/DeltaS1",";Delta S1;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/MT_z",";MT(#gamma_{1},%MET);EventsBIN",50,0,1000);
   ADD_HIST("pre_ph165/MT_zz",";MT(#gamma_{1},%MET);EventsBIN",50,0,500);
   ADD_HIST("pre_ph165/phiMETph" ,";#Delta#phi(%MET,#gamma_{1});EventsBIN",50,-3.2,3.2);

   ADD_HIST("pre_ph165/DeltaPhiMETjet100",";#Delta#phi(%MET,jet_{i} PT > 100 GeV);EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/phiMETjet",";#Delta#phi(%MET,jet_{1});EventsBIN",50,-5,5);
   ADD_HIST("pre_ph165/phiMETnJet",";#Delta#phi(%MET,nearest jet);EventsBIN",50,-5,5);
   ADD_HIST("pre_ph165/phiMETn2Jet",";#Delta#phi(%MET,nearest jet_{1/2});EventsBIN",50,-5,5);
   ADD_HIST("pre_ph165/phiMETnJetPh",";#Delta#phi(%MET,nearest jet/photon);EventsBIN",50,-5,5);
   ADD_HIST("pre_ph165/phiJetPh",";#Delta#phi(#gamma_{1},jet_{1});EventsBIN",50,-5,5);
   ADD_HIST("pre_ph165/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/METdotJet",";%VMET#upoint%VPT(jet_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/absphiMETHT",";|#Delta#phi|(%MET,#vec{HT});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/absphiPhHT",";|#Delta#phi|(#gamma_{1},#vec{HT});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/vecHTPt",";HT PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});

   ADD_HIST("pre_ph165/MT100/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/MT100/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
                                 
   ADD_HIST("pre_ph165/dPhiMETphl06/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1410},{20,30,50,10});
   ADD_HIST("pre_ph165/dPhiMETphl06/MET",";MET;Events / bin",{0,300,600,1000,1410},{20,30,50,10});
   ADD_HIST("pre_ph165/dPhiMETphl06/sieie",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);
   
   //2016 studies
   ADD_HIST("pre_ph165/MT100_l300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/MT100_l300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/MT100_l300/METS",";METSIG;Events / bin",{20,100},{5});
   ADD_HIST("pre_ph165/MT100_l300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/MT100_l300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/MT100_l300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/MT100_l300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/MT100_l300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/MT100_l300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/MT100_l300/STg"   ,";STg;EventsBIN"           ,150,0,1500);

   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_MET50/MT100/METl200vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);

   ////2016 candidate? ----
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/DeltaPhiMETjet100",";#Delta#phi(%MET,jet_{i} PT > 100 GeV);EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/vecHTPt",";HT PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1500},{20,30,50,100});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/MET",";p_{T}^{miss};Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETph",";|#Delta#phi|(p_{T}^{miss},#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJet",";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh",";|#Delta#phi(p_{T}^{miss},nearest jet/#gamma)| (radians);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_phoPtl500",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_phoPtl700",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);   
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_phoPtg700",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);   
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_JESu",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_JESd",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETHT",";|#Delta#phi|(%MET,#vec{HT});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiPhHT",";|#Delta#phi|(#gamma_{1},#vec{HT});EventsBIN",50,0,5);

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETphl02/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1500},{20,30,50,100});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1500},{20,30,50,100});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/STg"   ,";STg;EventsBIN"           ,150,0,1500);   

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/r9_phoPtl500",";R_{9}(#gamma)",110,0,1.1);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/r9_phoPtl700",";R_{9}(#gamma)",110,0,1.1);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/r9_phoPtg700",";R_{9}(#gamma)",110,0,1.1);
   ADD_HIST("pre_ph165/c_MET100/MT100/r9",";R_{9}(#gamma)",110,0,1.1);  

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/relPt2Jets",";relative Pt of 2 lead. jets;EventsBIN",25,0,1);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/DeltaS",";Delta S;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/DeltaS1",";Delta S1;EventsBIN",32,0,3.2);

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/Ngl15"   ,";N_{#gamma loose,15};EventsBIN" ,5 ,-.5, 5-.5); 
   
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/Nj"    ,";N_{jets};EventsBIN"     ,8,0-0.5,8-0.5);

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/absphiMETnJetPh",";|#Delta#phi|(p_{T}^{miss},nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/STg"   ,";STg;EventsBIN"           ,150,0,1500);

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/absphiMETnJetPh",";|#Delta#phi|(p_{T}^{miss},nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/STg"   ,";STg;EventsBIN"           ,150,0,1500);

   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/absphiMETnJetPh",";|#Delta#phi|(p_{T}^{miss},nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/STg"   ,";STg;EventsBIN"           ,150,0,1500);
        
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/absphiMETnJetPh",";|#Delta#phi|(p_{T}^{miss},nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/STg"   ,";STg;EventsBIN"           ,150,0,1500);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/relPt2Jets",";relative Pt of 2 lead. jets;EventsBIN",25,0,1);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/DeltaS",";Delta S;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/DeltaS1",";Delta S1;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/MT"    ,";MT(#gamma_{1},p_{T}^{miss});EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/MET",";p_{T}^{miss};Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET100/MT100/METl400vMTl400/Nj"    ,";N_{jets};EventsBIN"     ,8,0-0.5,8-0.5);
   
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET150/MT150/METl400vMTl400/Nj"    ,";N_{jets};EventsBIN"     ,8,0-0.5,8-0.5);
   
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET150/MT100/METl300vMTl300/Nj"    ,";N_{jets};EventsBIN"     ,8,0-0.5,8-0.5);    

   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);

   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);    
   
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);   
   ADD_HIST("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);


  // 2016 SR
   ADD_HIST("pre_ph165/c_MET300/MT300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg"   ,";STg;EventsBIN"           ,2000,0,2000);  
   ADD_HIST("pre_ph165/c_MET300/MT300/relPt2Jets",";relative Pt of 2 lead. jets;EventsBIN",25,0,1);
   ADD_HIST("pre_ph165/c_MET300/MT300/DeltaS",";Delta S;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET300/MT300/DeltaS1",";Delta S1;EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET300/MT300/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/c_MET300/MT300/Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre_ph165/c_MET300/MT300/Ngl15"   ,";N_{#gamma loose,15};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre_ph165/c_MET300/MT300/ph1Pt"   ,";#gamma_{1} PT;EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET300/MT300/vecHTPt",";HT PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   
   ADD_HIST("pre_ph165/c_MET400/MT400/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_MET400/MT400/STg"   ,";STg;EventsBIN"           ,200,0,2000);
 
   
   // CR
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/METS",";METSIG;Events / bin",{20,100},{5});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETn2Jet",";|#Delta#phi|(%MET,nearest jet_{1/2});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/METdotJet",";%VMET#upoint%VPT(jet_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/STg"   ,";STg;EventsBIN"           ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/Nj"    ,";N_{jets};EventsBIN"     ,8,0-0.5,8-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/Nele"  ,";N_{electrons};EventsBIN",5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/Nmu"   ,";N_{muons};EventsBIN"    ,5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/Nl"    ,";N_{leptons};EventsBIN"  ,5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/sieie1",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/sipip1",";SIPIP(#gamma_{1});EventsBIN"      ,300,0,0.03);

   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETjetJESu",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETjetJESd",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);

   // with b veto
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METS",";METSIG;Events / bin",{20,100},{5});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETn2Jet",";|#Delta#phi|(%MET,nearest jet_{1/2});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METdotJet",";%VMET#upoint%VPT(jet_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/STg"   ,";STg;EventsBIN"           ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nj"    ,";N_{jets};EventsBIN"     ,8,0-0.5,8-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nele"  ,";N_{electrons};EventsBIN",5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nmu"   ,";N_{muons};EventsBIN"    ,5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nl"    ,";N_{leptons};EventsBIN"  ,5,0-0.5,5-0.5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/sieie1",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/sipip1",";SIPIP(#gamma_{1});EventsBIN"      ,300,0,0.03);

   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETjetJESu",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETjetJESd",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);

   // other WP
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bL/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0bM/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);

   // beam halo would be mostly at eta = +- 1.4
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/sieie1",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);

   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/METS",";METSIG;Events / bin",{20,100},{5});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/STg"   ,";STg;EventsBIN"           ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);

   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/METS",";METSIG;Events / bin",{20,100},{5});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/METSHT",";METSHT;EventsBIN"        ,50,0,50);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/STg"   ,";STg;EventsBIN"           ,150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/MT"    ,";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);

   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/ph1Pt",";#gamma_{1} PT;Events / bin",{0,300,600,1000,1010},{20,30,50,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/MET",";MET;Events / bin",{0,100,400,500,1000,1010},{20,30,50,100,10});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/METS",";METSIG;Events / bin",{20,100},{5});
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETnJet",";|#Delta#phi|(%MET,nearest jet);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETn2Jet",";|#Delta#phi|(%MET,nearest jet_{1/2});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETnJetPh",";|#Delta#phi|(%MET,nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);
   ADD_HIST("pre_ph165/c_S30/MT100/Sl100vMTl300/METdotJet",";%VMET#upoint%VPT(jet_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);


   // SR
   ADD_HIST("pre_ph165/c_S80/MT300/ph1Pt",";#gamma_{1} PT;EventsBIN",100,0,1000);
   ADD_HIST("pre_ph165/c_S80/MT300/MET",";MET;EventsBIN",100,0,1000);
   ADD_HIST("pre_ph165/c_S80/MT300/METS",";METSIG;EventsBIN",200,100,2100);
   ADD_HIST("pre_ph165/c_S80/MT300/MT",";MT(#gamma_{1},%MET);EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_S80/MT300/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S80/MT300/absphiMETjet",";|#Delta#phi|(%MET,jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S80/MT300/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_S80/MT300/METSHT",";METSHT;EventsBIN"        ,100,0,100);
   ADD_HIST("pre_ph165/c_S80/MT300/STg"   ,";STg;EventsBIN"           ,200,0,2000);
   ADD_HIST("pre_ph165/c_S80/MT300/sieie1",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S80/MT300/sipip1",";SIPIP(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S80/MT300/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);
   // beam halo would be mostly at eta = +- 1.4
   ADD_HIST("pre_ph165/c_S80/MT300/eta1/absphiMETph",";|#Delta#phi|(%MET,#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S80/MT300/eta1/sieie1",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S80/MT300/eta1/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);

   // excluding VR
   ADD_HIST("pre_ph165/c_S80/MT300/STg600/STg"   ,";STg;EventsBIN"           ,200,0,2000);
   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/STg"   ,";STg;EventsBIN"           ,300,0,3000);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/STg_SRbin"   ,";STg;EventsBIN" ,{600,800,1000,1300,1600},{200,200,300,300});
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);

   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/Ngl15"   ,";N_{#gamma loose,15};EventsBIN" ,5 ,-.5, 5-.5);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/ph1Pt"   ,";#gamma_{1} PT;EventsBIN",150,0,1500);   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/r9",";R_{9}(#gamma)",110,0,1.1);  
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/r9_Stgg1300",";R_{9}(#gamma)",110,0,1.1);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/absphiMETph",";|#Delta#phi|(p_{T}^{miss},#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/MET",";p_{T}^{miss};EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/DeltaPhiMETjet100",";#Delta#phi(p_{T}^{miss},jet_{i} PT > 100 GeV);EventsBIN",32,0,3.2);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/absphiMETnJetPh",";|#Delta#phi|(p_{T}^{miss},nearest jet/#gamma);EventsBIN",50,0,5);


   ADD_HIST("pre_ph165/c_MET300/MT300/absphiMETHT",";|#Delta#phi|(p_{T}^{miss},#vec{HT});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET300/MT300/absphiPhHT",";|#Delta#phi|(#gamma_{1},#vec{HT});EventsBIN",50,0,5);

   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/STg",";STg;EventsBIN"           ,300,0,3000);
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/ph1Pt",";#gamma_{1} PT;EventsBIN",150,0,1500);   
   ADD_HIST("pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/MET",";p_{T}^{miss};EventsBIN",150,0,1500);
   
   // validation region
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/ph1Pt",";#gamma_{1} PT;EventsBIN",100,0,1000);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/MET",";p_{T}^{miss};EventsBIN",100,0,1000);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/METS",";METSIG;EventsBIN",200,100,2100);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/MT",";MT(#gamma_{1},p_{T}^{miss});EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/absphiMETph",";|#Delta#phi|(p_{T}^{miss},#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/HT"    ,";HT;EventsBIN"            ,150,0,1500);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/METSHT",";METSHT;EventsBIN"        ,100,0,100);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/STg"   ,";STg;EventsBIN"           ,200,0,2000);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/sieie1",";SIEIE(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/sipip1",";SIPIP(#gamma_{1});EventsBIN"      ,300,0,0.03);
   ADD_HIST("pre_ph165/c_S80/MT300/STgl600/Nph"   ,";N_{#gamma};EventsBIN"   ,5,0-0.5,5-0.5);

   // validation region 2016
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/ph1Pt",";#gamma_{1} PT;EventsBIN",100,0,1000);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/MET",";p_{T}^{miss};EventsBIN",100,0,1000);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/MT",";MT(#gamma_{1},p_{T}^{miss});EventsBIN",150,0,1500);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/absphiMETph",";|#Delta#phi|(p_{T}^{miss},#gamma_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/absphiMETjet",";|#Delta#phi|(p_{T}^{miss},jet_{1});EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/absphiMETnJetPh",";|#Delta#phi|(p_{T}^{miss},nearest jet/#gamma);EventsBIN",50,0,5);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/STg"   ,";STg;EventsBIN"           ,200,0,2000);
   ADD_HIST("pre_ph165/c_MET300/MT300/STgl600/METdotPh",";%VMET#upoint%VPT(#gamma_{1}) [GeV^{2}];EventsBIN",50,-100000,100000);     

   // PHOTON135MET100

   // HLT_PFMET170

   hist::Histograms<TH2F> h2s(vsDatasubsets);
   h2s.addHist("METS MT",";METSIG;MT(#gamma_{1},p_{T}^{miss});Events / bin",40,0,2000,40,0,2000);
   h2s.addHist("METS pt",";METSIG;#gamma_{1} PT;Events / bin",40,0,2000,40,0,1000);
   h2s.addHist("METS MET",";METSIG;MET;Events / bin",40,0,2000,40,0,1000);
   h2s.addHist("METSHT MT",";METSHT;MT(#gamma_{1},p_{T}^{miss});Events / bin",40,0,80,40,0,2000);
   h2s.addHist("METSHT METS",";METSHT;METSIG;Events / bin",40,0,80,40,0,2000);
   h2s.addHist("Iso/ec",";#gamma_{1} PT;EcalPFCluster Iso (GeV);Events / bin",50,0,1000,50,0,40);
   h2s.addHist("Iso/hc",";#gamma_{1} PT;HcalPFCluster Iso (GeV);Events / bin",50,0,1000,50,0,40);
   h2s.addHist("Iso/tr",";#gamma_{1} PT;SC Track Iso (GeV);Events / bin",50,0,1000,50,0,40);
   h2s.addHist("pt MET",";#gamma_{1} PT;MET;Events / bin",50,0,1000,50,0,1000);
   h2s.addHist("pt MT",";#gamma_{1} PT;MT(#gamma_{1},p_{T}^{miss});Events / bin",50,0,1000,40,0,2000);
   h2s.addHist("Ngl Ngpix",";N_{#gamma loose};N_{#gamma pix};Events / bin",5 ,-.5, 5-.5,5 ,-.5, 5-.5);
   
   h2s.addHist("pre_ph165/MET MT",";MET;MT(#gamma_{1},p_{T}^{miss});Events / bin",50,0,1000,50,0,1000);  
   h2s.addHist("pre_ph165/dPhiMETph pt",";|#Delta#phi|(p_{T}^{miss},#gamma_{1});PT;Events / bin",32,0,3.2,50,0,1000);
   h2s.addHist("pre_ph165/METS MT",";METSIG;MT(#gamma_{1},p_{T}^{miss});Events / bin",50,0,500,50,0,1000);
   h2s.addHist("pre_ph165/c_S100/MT300/METS STg",";METSIG;STg;Events / bin",200,100,2100,200,0,2000);
   h2s.addHist("pre_ph165/c_S80/MT300/METS STg",";METSIG;STg;Events / bin",200,100,2100,200,0,2000);

   h2s.addCounter("total");
   h2s.addCounter("cr"); // control region
   h2s.addCounter("vr"); // validation region
   h2s.addCounter("!Iso40"); // not passing trigger iso



 //  TH1F *trigger_photon_pt_n = new TH1F("numerator);

   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      hs_notPix.setCurrentSample(dss.name);
      hs_pix.setCurrentSample(dss.name);
      h2s.setCurrentSample(dss.name);

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
      int events =0;
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
         h2s.setFillWeight(fEventWeight);

         std::vector<tree::Photon const *> lPho,mPho,tPho,lPho15,mPho15,tPho15,lPixPho;
         for (tree::Photon const &ph: *photons){
            (ph.hasPixelSeed?hs_pix:hs_notPix).fill("pre/phEta",ph.p.Eta());

            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
    //        if (ph.r9 < 0.9 || ph.r9 > 1.0) continue;
            if (fabs(ph.p.Eta())>1.4442) continue;
            if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
            if (ph.hasPixelSeed){
               hs_pix.fill("pre/phPt",ph.p.Pt());
               lPixPho.push_back(&ph);
            } else {
               hs_notPix.fill("pre/phPt",ph.p.Pt());
               //Attention check lossePhoton 15 or new
               if (ph.isLoose) lPho.push_back(&ph);
               if (ph.isMedium) mPho.push_back(&ph);
               if (ph.isTight)  tPho.push_back(&ph);
               if (ph.isLoose15) lPho15.push_back(&ph);
               if (ph.isMedium15) mPho15.push_back(&ph);
               if (ph.isTight15)  tPho15.push_back(&ph);                        
            }
         }

         hs_notPix.fill("pre/Ngl",lPho.size());
         hs_notPix.fill("pre/Ngl15",lPho15.size());         
         hs_notPix.fill("pre/Ngm",mPho.size());
         hs_notPix.fill("pre/Ngt",tPho.size());
         hs_notPix.fill("pre/Ngpix",lPixPho.size());
         h2s.fill("Ngl Ngpix",lPho.size(),lPixPho.size());



         // independent of "normal"/"pixel" run
         // jet related
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         int const Njets=cjets.size();
         int NbL=0, NbM=0, NbT=0;
         // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
         for (auto const &j: cjets) {
            if (j.bDiscriminator > 0.5426) NbL++;
            if (j.bDiscriminator > 0.8484) NbM++;
            if (j.bDiscriminator > 0.9535) NbT++;
         }
         float const HT=phys::computeHT(cjets);
         
         TVector3 vecHT, vecGammaHT;
         bool clean_MET = true;
         
         for (auto const &jet: cjets) {
            vecHT+=jet.p;
            if (jet.p.Pt() < 100) continue;            
            if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
         }

         if (!clean_MET) continue;
    
         float const dPhiMETjet=Njets>0 ? MET->p.DeltaPhi(cjets[0].p) : 4;
         float const METdotJet=Njets>0 ? MET->p.Dot(cjets[0].p) : 0;

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
         const int Nlepton=Nele+Nmu;

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
            h2s.setFillWeight(fEventWeight);

            float minDR=std::numeric_limits<float>::max();
            for (tree::Jet const &jet: *jets){
               if (jet.isLoose && jet.p.Pt()>30){
                  float const dr=jet.p.DeltaR(pho[0]->p);
                  float const dpt=fabs(jet.p.Pt()-phoPt)/phoPt;
                  if (dr<.1 && dpt<.5) continue;
                  minDR=std::min(minDR,dr);
               }
            }
            hs.fill("pre/minDR",minDR);
            //jet und photon separated
            if (minDR < .5) continue;

               int const Nph=pho.size();                  
               float const METSHT= phys::METoverSqrtHT(MET->p.Pt(),HT);
               float const MT=phys::M_T(*pho[0],*MET);
   
               float const dPhiMETph=MET->p.DeltaPhi(pho[0]->p);
               float const dPhiJetPh=Njets>0 ? pho[0]->p.DeltaPhi(cjets[0].p) : 4;
               float const dPhiMETHT=Njets>0 ? MET->p.DeltaPhi(vecHT) : 4;
               float const dPhiPhHT=Njets>0 ? pho[0]->p.DeltaPhi(vecHT) : 4;                            
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

            if (met > 100 && MT > 100){
               if (*baseMETTr){
                  hs.fill("effi/trigger_effi_d",phoPt);
                  if (*trigger_Ph){
                     hs.fill("effi/trigger_effi_n",phoPt);
                  }
               }
               if (phoPt > 180 && *baseMETTr){
                  hs.fill("effi/trigger_effi_d_cut",phoPt);
                  if (*trigger_Ph){
                     hs.fill("effi/trigger_effi_n_cut",phoPt);
                  }
               }
   
               
               if (*baseMETTr){
                  hs.fill("effi/trigger_effi_d_match",phoPt);
                  if (*trigger_Ph && triggerMatch){
                     hs.fill("effi/trigger_effi_n_match",phoPt);
                  }
               }
               if (phoPt > 180 && *baseMETTr){
                  hs.fill("effi/trigger_effi_d_cut_match",phoPt);
                  if (*trigger_Ph && triggerMatch){
                     hs.fill("effi/trigger_effi_n_cut_match",phoPt);
                  }
               }
            } 


            if (phoPt>180 && (!isData || (*trigger_Ph && triggerMatch))){
            
               for (auto const &jet: cjets) {
                  if (jet.p.Pt() < 100) continue;
                  hs.fill("pre_ph165/DeltaPhiMETjet100",std::fabs(MET->p.DeltaPhi(jet.p)));                 
               }
               
               hs.fill("pre_ph165/Ngl",lPho.size());
               hs.fill("pre_ph165/Ngl15",lPho15.size()); 
               hs.fill("pre_ph165/ph1Pt",phoPt);
               hs.fill("pre_ph165/MET",met);
               hs.fill("pre_ph165/METS",MET->sig);
               hs.fill("pre_ph165/METS_z",MET->sig);
               hs.fill("pre_ph165/METS_long",MET->sig);
               hs.fill("pre_ph165/METS_logx",MET->sig);
               hs.fill("pre_ph165/METS_fine",MET->sig);
               hs.fill("pre_ph165/HT",HT);
               hs.fill("pre_ph165/METSHT",METSHT);
               hs.fill("pre_ph165/MT",MT);
               hs.fill("pre_ph165/MT_z",MT);
               hs.fill("pre_ph165/MT_zz",MT);
               hs.fill("pre_ph165/phiMETph",dPhiMETph);
               hs.fill("pre_ph165/phiMETjet",dPhiMETjet);
               hs.fill("pre_ph165/phiMETnJet",dPhiMETnearJet);
               hs.fill("pre_ph165/phiMETn2Jet",dPhiMETnear2Jet);
               hs.fill("pre_ph165/phiMETnJetPh",dPhiMETnearJetPh);
               hs.fill("pre_ph165/phiJetPh",dPhiJetPh);
               hs.fill("pre_ph165/METdotPh",METdotPh);
               hs.fill("pre_ph165/METdotJet",METdotJet);
               hs.fill("pre_ph165/DeltaS",DeltaS);
               hs.fill("pre_ph165/DeltaS1",DeltaS1);              
               hs.fill("pre_ph165/relPt2Jets",relPt2Jets);
               hs.fill("pre_ph165/absphiMETHT",std::abs(dPhiMETHT));
               hs.fill("pre_ph165/absphiPhHT",std::abs(dPhiPhHT));
               hs.fill("pre_ph165/vecHTPt",vecHT.Pt());
               
               h2s.fill("pre_ph165/METS MT",MET->sig,MT);
               h2s.fill("pre_ph165/MET MT",met,MT);
               h2s.fill("pre_ph165/dPhiMETph pt",std::abs(dPhiMETph),phoPt);              
               
               if (std::abs(dPhiMETph) < 0.6 ) {
                  hs.fill("pre_ph165/dPhiMETphl06/ph1Pt",phoPt);
                  hs.fill("pre_ph165/dPhiMETphl06/sieie",pho[0]->sigmaIetaIeta);
                  hs.fill("pre_ph165/dPhiMETphl06/MET",met);                                             
               }

            //   if (METdotPh > 0) continue;

               //2016 studies

               if (MT > 100) {
                  hs.fill("pre_ph165/MT100/METdotPh",METdotPh);
                  hs.fill("pre_ph165/MT100/ph1Pt",phoPt);
                                    
                  if (MT < 300){
                     hs.fill("pre_ph165/MT100_l300/ph1Pt",phoPt);
                     hs.fill("pre_ph165/MT100_l300/MET",met);
                     hs.fill("pre_ph165/MT100_l300/absphiMETph",std::abs(dPhiMETph));
                     hs.fill("pre_ph165/MT100_l300/absphiMETjet",std::abs(dPhiMETjet));
                     hs.fill("pre_ph165/MT100_l300/absphiMETnJet",std::abs(dPhiMETnearJet));
                     hs.fill("pre_ph165/MT100_l300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                     hs.fill("pre_ph165/MT100_l300/HT",HT);
                     hs.fill("pre_ph165/MT100_l300/METSHT",METSHT);
                     hs.fill("pre_ph165/MT100_l300/METS",MET->sig);
                     hs.fill("pre_ph165/MT100_l300/STg",STg);
                  }
                  if (met > 50){
                     if (met < 200 || MT < 300){
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/MET",met);
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETph",std::abs(dPhiMETph));
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/HT",HT);
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/METSHT",METSHT);
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/STg",STg);
                        hs.fill("pre_ph165/c_MET50/MT100/METl200vMTl300/MT",MT);
                     }
                     if (met > 100){
                                        

                           hs.fill("pre_ph165/c_MET100/MT100/r9",pho[0]->r9);
                        if (met < 300 || MT < 300){
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/ph1Pt",phoPt);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/MET",met);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETph",std::abs(dPhiMETph));
                           if (std::abs(dPhiMETph) < 0.2) hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETphl02/ph1Pt",phoPt);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                           if (std::abs(dPhiMETph) > 0.3) {
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/ph1Pt",phoPt);
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/STg",STg);
                           }
                           if (phoPt < 500){
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_phoPtl500",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/r9_phoPtl500",pho[0]->r9);
                           }else if (phoPt < 700){
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_phoPtl700",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/r9_phoPtl700",pho[0]->r9);
                           }else{
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_phoPtg700",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/r9_phoPtg700",pho[0]->r9);
                           }
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETHT",std::abs(dPhiMETHT));
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiPhHT",std::abs(dPhiPhHT));
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/vecHTPt",vecHT.Pt());

                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_JESu",std::abs(dPhiMETnearJetPh_JESu));
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh_JESd",std::abs(dPhiMETnearJetPh_JESd));

                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/DeltaS",DeltaS);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/DeltaS1",DeltaS1);              
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/relPt2Jets",relPt2Jets);

                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/Ngl",lPho.size());
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/Ngl15",lPho15.size());                        
                           
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/HT",HT);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/METSHT",METSHT);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/STg",STg);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/MT",MT);
                           hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/Nj",Njets);
                           for (auto const &jet: cjets) {
                              if (jet.p.Pt() < 100) continue;
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/DeltaPhiMETjet100",std::fabs(MET->p.DeltaPhi(jet.p)));                 
                           }
                           
                           
                           if (Njets < 3){
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/absphiMETjet",std::abs(dPhiMETjet));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/STg",STg);                                                    
                           }
                           if (NbT == 0){
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/absphiMETjet",std::abs(dPhiMETjet));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/STg",STg);                                                                
                           }
                           if (Nlepton == 0){
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/absphiMETjet",std::abs(dPhiMETjet));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/STg",STg);                                                                
                           }
                           if (Nlepton == 1){
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/absphiMETjet",std::abs(dPhiMETjet));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/STg",STg);
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/MT",MT);
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/DeltaS",DeltaS);
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/DeltaS1",DeltaS1);              
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/relPt2Jets",relPt2Jets);
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/ph1Pt",phoPt);
                              hs.fill("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/MET",met);                                                            
                           }                                                                           
                        }
                        if (met < 400 || MT < 400){
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/ph1Pt",phoPt);
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/MET",met);
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETph",std::abs(dPhiMETph));
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETjet",std::abs(dPhiMETjet));
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETnJet",std::abs(dPhiMETnearJet));
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/HT",HT);
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/METSHT",METSHT);
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/STg",STg);
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/MT",MT);
                           hs.fill("pre_ph165/c_MET100/MT100/METl400vMTl400/Nj",Njets);                           
                        }
                        if (met > 150 && MT > 150){
                           if (met < 400 || MT < 400){
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/ph1Pt",phoPt);
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/MET",met);
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETph",std::abs(dPhiMETph));
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETjet",std::abs(dPhiMETjet));
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETnJet",std::abs(dPhiMETnearJet));
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/HT",HT);
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/METSHT",METSHT);
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/STg",STg);
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/MT",MT);
                              hs.fill("pre_ph165/c_MET150/MT150/METl400vMTl400/Nj",Njets);
                           }
                        }                                            
                        if( met > 150){
                            if (met < 300 || MT < 300){
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/ph1Pt",phoPt);
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/MET",met);
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETph",std::abs(dPhiMETph));
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/HT",HT);
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/METSHT",METSHT);
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/STg",STg);
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/MT",MT);
                              hs.fill("pre_ph165/c_MET150/MT100/METl300vMTl300/Nj",Njets);
                           }
                        }                      
                     }
                                         
                  }                 
                  
                  if (METSHT > 5){
                     if (METSHT < 20 || MT < 300){
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/MET",met);
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETph",std::abs(dPhiMETph));
                        if (Njets > 0) hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/HT",HT);
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/METSHT",METSHT);
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/STg",STg);
                        hs.fill("pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/MT",MT);
                     }                    
                  }
                  if (METSHT > 10){
                     if (METSHT < 20 || MT < 300){
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/MET",met);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETph",std::abs(dPhiMETph));
                        if (Njets > 0) hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/HT",HT);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/METSHT",METSHT);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/STg",STg);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/MT",MT);
                     }
                     if (METSHT < 15 || MT < 300){
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/MET",met);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETph",std::abs(dPhiMETph));
                        if (Njets > 0) hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/HT",HT);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/METSHT",METSHT);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/STg",STg);
                        hs.fill("pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/MT",MT);
                     }                                         
                  }                  
               }

               if (met > 300 && MT > 300){
                  hs.fill("pre_ph165/c_MET300/MT300/STg",STg);
                  hs.fill("pre_ph165/c_MET300/MT300/METSHT",METSHT);
                  hs.fill("pre_ph165/c_MET300/MT300/DeltaS",DeltaS);
                  hs.fill("pre_ph165/c_MET300/MT300/DeltaS1",DeltaS1);              
                  hs.fill("pre_ph165/c_MET300/MT300/relPt2Jets",relPt2Jets);
                  hs.fill("pre_ph165/c_MET300/MT300/METdotPh",METdotPh);
                  hs.fill("pre_ph165/c_MET300/MT300/Ngl",lPho.size());
                  hs.fill("pre_ph165/c_MET300/MT300/Ngl15",lPho15.size());
                  hs.fill("pre_ph165/c_MET300/MT300/ph1Pt",phoPt);
                  hs.fill("pre_ph165/c_MET300/MT300/absphiMETHT",std::abs(dPhiMETHT));
                  hs.fill("pre_ph165/c_MET300/MT300/absphiPhHT",std::abs(dPhiPhHT));
                  hs.fill("pre_ph165/c_MET300/MT300/vecHTPt",vecHT.Pt());             
                  if (met > 400 && MT > 400){
                     hs.fill("pre_ph165/c_MET400/MT400/STg",STg);
                     hs.fill("pre_ph165/c_MET400/MT400/METSHT",METSHT);                                        
                  }
                  if (STg<600) {
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/ph1Pt",phoPt);
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/MET",met);
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/MT",MT);
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/absphiMETph",std::abs(dPhiMETph));
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/absphiMETjet",std::abs(dPhiMETjet));
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/STg",STg);
                     hs.fill("pre_ph165/c_MET300/MT300/STgl600/METdotPh",METdotPh);
                  } else {
                     for (auto const &jet: cjets) {
                        if (jet.p.Pt() < 100) continue;
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/DeltaPhiMETjet100",std::fabs(MET->p.DeltaPhi(jet.p)));                 
                     }
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/STg",STg);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/STg_SRbin",STg);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/METdotPh",METdotPh);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/Ngl",lPho.size());
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/Ngl15",lPho15.size());
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/r9",pho[0]->r9);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/ph1Pt",phoPt);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/MET",met);
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                     hs.fill("pre_ph165/c_MET300/MT300/STg600/absphiMETph",std::abs(dPhiMETph));
                     if (std::abs(dPhiMETph) > 0.3) {
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/STg",STg);
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/MET",met);
                     }
                     if (STg > 1300){
                        hs.fill("pre_ph165/c_MET300/MT300/STg600/r9_Stgg1300",pho[0]->r9);
                     }
                  }
               }
               


               //2015 CR
               /*
               if (MET->sig>30){                  
                  if (MT>100) {
                     if (MET->sig<80 || MT<300) {
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/MET",met);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/METS",MET->sig);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETph",std::abs(dPhiMETph));
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETn2Jet",std::abs(dPhiMETnear2Jet));
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/METdotPh",METdotPh);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/METdotJet",METdotJet);

           //             hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETjetJESu",std::abs(Njets>0 ? MET_JESu->p.DeltaPhi(cjets[0].p) : 4));
           //             hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/absphiMETjetJESd",std::abs(Njets>0 ? MET_JESd->p.DeltaPhi(cjets[0].p) : 4));

                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/HT",HT);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/METSHT",METSHT);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/STg",STg);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/MT",MT);

                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/Nph",Nph);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/Nj",Njets);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/Nele",Nele);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/Nmu",Nmu);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/Nl",Nlepton);

                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/sieie1",pho[0]->sigmaIetaIeta);
                        hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/sipip1",pho[0]->sigmaIphiIphi);

                        if (NbL==0) {
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bL/absphiMETjet",std::abs(dPhiMETjet));
                        }
                        if (NbM==0) {
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bM/absphiMETjet",std::abs(dPhiMETjet));
                        }
                        if (NbT==0) {
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/ph1Pt",phoPt);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/MET",met);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METS",MET->sig);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETph",std::abs(dPhiMETph));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETjet",std::abs(dPhiMETjet));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETnJet",std::abs(dPhiMETnearJet));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETn2Jet",std::abs(dPhiMETnear2Jet));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METdotPh",METdotPh);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METdotJet",METdotJet);

               //            hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETjetJESu",std::abs(Njets>0 ? MET_JESu->p.DeltaPhi(cjets[0].p) : 4));
               //            hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/absphiMETjetJESd",std::abs(Njets>0 ? MET_JESd->p.DeltaPhi(cjets[0].p) : 4));

                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/HT",HT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/METSHT",METSHT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/STg",STg);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/MT",MT);

                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nph",Nph);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nj",Njets);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nele",Nele);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nmu",Nmu);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/Nl",Nlepton);

                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/sieie1",pho[0]->sigmaIetaIeta);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/sipip1",pho[0]->sigmaIphiIphi);
                        }
                        if (std::abs(pho[0]->p.Eta()) > 1) {
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/absphiMETph",std::abs(dPhiMETph));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/sieie1",pho[0]->sigmaIetaIeta);
                        }

                        if (Nlepton==0) {
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/ph1Pt",phoPt);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/MET",met);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/METS",MET->sig);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/absphiMETjet",std::abs(dPhiMETjet));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/HT",HT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/METSHT",METSHT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/STg",STg);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/MT",MT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/Nph",Nph);
                        } else if (Nlepton==1) {
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/ph1Pt",phoPt);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/MET",met);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/METS",MET->sig);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/absphiMETjet",std::abs(dPhiMETjet));
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/HT",HT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/METSHT",METSHT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/STg",STg);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/MT",MT);
                           hs.fill("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/Nph",Nph);
                        }
                     }
                     if (MET->sig<100 || MT<300) {
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/MET",met);
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/METS",MET->sig);
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETph",std::abs(dPhiMETph));
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETnJet",std::abs(dPhiMETnearJet));
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETn2Jet",std::abs(dPhiMETnear2Jet));
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/absphiMETnJetPh",std::abs(dPhiMETnearJetPh));
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/METdotPh",METdotPh);
                        hs.fill("pre_ph165/c_S30/MT100/Sl100vMTl300/METdotJet",METdotJet);
                     }
                  }
               } // S>30
               */
               /*
               if (MET->sig>80) {
                  if (MT>300) {
                     hs.fill("pre_ph165/c_S80/MT300/ph1Pt",phoPt);
                     hs.fill("pre_ph165/c_S80/MT300/MET",met);
                     hs.fill("pre_ph165/c_S80/MT300/METS",MET->sig);
                     hs.fill("pre_ph165/c_S80/MT300/MT",MT);
                     hs.fill("pre_ph165/c_S80/MT300/absphiMETph",std::abs(dPhiMETph));
                     hs.fill("pre_ph165/c_S80/MT300/absphiMETjet",std::abs(dPhiMETjet));
                     hs.fill("pre_ph165/c_S80/MT300/HT",HT);
                     hs.fill("pre_ph165/c_S80/MT300/METSHT",METSHT);
                     hs.fill("pre_ph165/c_S80/MT300/STg",STg);
                     hs.fill("pre_ph165/c_S80/MT300/sieie1",pho[0]->sigmaIetaIeta);
                     hs.fill("pre_ph165/c_S80/MT300/sipip1",pho[0]->sigmaIphiIphi);
                     hs.fill("pre_ph165/c_S80/MT300/Nph",Nph);

                     h2s.fill("pre_ph165/c_S80/MT300/METS STg",MET->sig,STg);

                     if (std::abs(pho[0]->p.Eta()) > 1) {
                        hs.fill("pre_ph165/c_S80/MT300/eta1/absphiMETph",std::abs(dPhiMETph));
                        hs.fill("pre_ph165/c_S80/MT300/eta1/sieie1",pho[0]->sigmaIetaIeta);
                        hs.fill("pre_ph165/c_S80/MT300/eta1/Nph",Nph);
                     }
                     if (STg<600) {
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/ph1Pt",phoPt);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/MET",met);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/METS",MET->sig);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/MT",MT);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/absphiMETph",std::abs(dPhiMETph));
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/absphiMETjet",std::abs(dPhiMETjet));
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/HT",HT);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/METSHT",METSHT);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/STg",STg);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/sieie1",pho[0]->sigmaIetaIeta);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/sipip1",pho[0]->sigmaIphiIphi);
                        hs.fill("pre_ph165/c_S80/MT300/STgl600/Nph",Nph);
                     } else {
                        hs.fill("pre_ph165/c_S80/MT300/STg600/STg",STg);
                     }
                  } // MT>300
               } // S>80 */
            } // phoPt>180

            // TRIGGER
            /* pay attention if definition of vairable is already implemented on this depth
            if (isData) {
               bool trigger = (phoPt>180 ? *trigger_Ph : *trigger_PhMET);
               if (!trigger) continue;
               fEventWeight=1;
            } else {
               // scale MC with trigger efficiency
               fEventWeight*=(phoPt>180 ? cfg.trigger_eff_Ph : cfg.trigger_eff_PhMET);
               if (wCalc) fEventWeight*=wCalc->get();
            }
            hs_notPix.setFillWeight(fEventWeight);
            hs_pix.setFillWeight(fEventWeight);
            h2s.setFillWeight(fEventWeight);

            // 2d-hists (no efake extimation here)
            if (pass==pass_normal){
               // Signal/Control region
               h2s.fill("METS MT",MET->sig,MT);
               h2s.fill("METS pt",MET->sig,phoPt);
               h2s.fill("METS MET",MET->sig,met);
               h2s.fill("METSHT MT",METSHT,MT);

               h2s.fill("pt MET",phoPt,met);

               // Trigger Isolation
               // h2s.fill("Iso/ec",phoPt,pho[0]->ecalPFClIso);
               // h2s.fill("Iso/hc",phoPt,pho[0]->hcalPFClIso);
               // h2s.fill("Iso/tr",phoPt,pho[0]->trackIso);

               // correlations
               h2s.fill("METSHT METS",METSHT,MET->sig);
               h2s.fill("pt MT",phoPt,MT);

               // bool const bIso40ec=pho[0]->ecalPFClIso<fIso40ec(phoPt);
               // bool const bIso40hc=pho[0]->hcalPFClIso<fIso40hc(phoPt);
               // bool const bIso40tr=pho[0]->trackIso   <fIso40tr(phoPt);
               // bool const bIso40=bIso40ec && bIso40hc && bIso40tr;
               // if (!bIso40) h2s.count("!Iso40");

               h2s.count("total");
            }*/

         } // normal/pixel pass

      } // evt loop
      io::log<<"";

      hs_notPix.scaleLumi();
      hs_notPix.mergeOverflow();
      hs_pix.scaleLumi();
      hs_pix.mergeOverflow();
      h2s.scaleLumi();
      file.Close();
   } // datasets

   // calling the "normal" histogram "hs" from here, since it's the most used
   hist::Histograms<TH1F> &hs = hs_notPix;
//   std::vector<TString> samplesToCombine={"GJets","QCD","ZNuNuGJets","ZGTo2LG","WGToLNuG","ZNuNuJets","WLNuJets","diboson","T5gg","T5Wg","GGM","SinglePhoton","MET"};
   std::vector<TString> samplesToCombine={"GJets_DR","QCD","ZNuNuGJets","ZGTo2LG","WGToLNuG","ZNuNuJets","WLNuJets","diboson","T5Wg","TChiWG","SinglePhoton","MET"};
   hs    .combineFromSubsamples(samplesToCombine);
   hs_pix.combineFromSubsamples(samplesToCombine);
   h2s   .combineFromSubsamples(samplesToCombine);
   io::RootFileSaver saver("plots.root",TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),false);
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      
      {"effi/",{"trigger_effi_n","trigger_effi_d","trigger_effi_n_cut","trigger_effi_d_cut",}},
      {"pre/",{"minDR","phPt","ph1Pt","phEta","jetPt","jet1Pt","MET","HT","METSHT","METS","METS_z","HTgen","Ngl","Ngl15","Ngm","Ngt","Ngpix","Njets","MT","phiMETph","phiMETjet"}},
      {"pre_ph165/"      ,{"Ngl","Ngl15","ph1Pt","MET","METS","METS_z","METS_long","METS_logx","METS_fine","HT","METSHT","MT","relPt2Jets", "DeltaS","DeltaS1","MT_z","MT_zz","phiMETph","phiMETjet","phiMETnJet","phiMETn2Jet","phiJetPh","METdotPh","METdotJet","phiMETnJetPh",
         }},
      {"pre_ph165/MT100/"      ,{"METdotPh","ph1Pt",
         }},         
      {"pre_ph165/dPhiMETphl06/"      ,{"ph1Pt","MET","sieie",
         }},         
      {"pre_ph165/MT100_l300/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg",
         }},
      {"pre_ph165/c_MET50/MT100/METl200vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/",{"Ngl","Ngl15","ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh","absphiMETnJetPh_JESu","absphiMETnJetPh_JESd",
                                             "HT","METSHT","STg","MT","Nj","relPt2Jets", "DeltaS","DeltaS1","absphiMETnJetPh_phoPtl500","absphiMETnJetPh_phoPtl700","absphiMETnJetPh_phoPtg700",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/",{"absphiMETjet","absphiMETnJetPh","STg",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/0b/",{"absphiMETjet","absphiMETnJetPh","STg",
         }},        
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/0l/",{"absphiMETjet","absphiMETnJetPh","STg",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/1l/",{"absphiMETjet","absphiMETnJetPh","STg","ph1Pt","MET","MT","relPt2Jets", "DeltaS","DeltaS1",
         }}, 
      {"pre_ph165/c_MET150/MT100/METl300vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT","Nj",
         }},
      {"pre_ph165/c_MET150/MT150/METl400vMTl400/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT","Nj",
         }},
      {"pre_ph165/c_MET100/MT100/METl400vMTl400/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT","Nj",
         }},   /*            
      {"pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},         
      {"pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},
      {"pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},*/
   //new SR
      {"pre_ph165/c_MET300/MT300/",{"METSHT","STg", "relPt2Jets", "DeltaS","DeltaS1","METdotPh","Ngl","Ngl15",
         }},
      {"pre_ph165/c_MET300/MT300/STg600/",{"STg","STg_SRbin","METdotPh","Ngl","Ngl15", "absphiMETnJetPh",
         }},        
      {"pre_ph165/c_MET400/MT400/",{"METSHT","STg",
         }},
   //OLD
 /*     {"pre_ph165/c_S30/MT100/Sl80vMTl300/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETn2Jet","absphiMETnJetPh","METdotPh","METdotJet",
                                             "absphiMETjetJESu","absphiMETjetJESd",
                                             "Nph","Nj","Nele","Nmu","Nl","HT","METSHT","STg","MT","sieie1","sipip1",
                                             "0bL/absphiMETjet","0bM/absphiMETjet",
         }},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETn2Jet","absphiMETnJetPh","METdotPh","METdotJet",
                                                 "absphiMETjetJESu","absphiMETjetJESd",
                                                 "Nph","Nj","Nele","Nmu","Nl","HT","METSHT","STg","MT","sieie1","sipip1",
         }},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/0l/",{"ph1Pt","MET","METS","absphiMETjet","Nph","HT","METSHT","STg","MT"}},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/1l/",{"ph1Pt","MET","METS","absphiMETjet","Nph","HT","METSHT","STg","MT"}},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/",{"absphiMETph","sieie1"}},
      {"pre_ph165/c_S30/MT100/Sl100vMTl300/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETn2Jet","absphiMETnJetPh","METdotPh","METdotJet",}},
   */   // VR
   //   {"pre_ph165/c_S80/MT300/STgl600/",{"ph1Pt","MET","MT","HT","METSHT","absphiMETjet","absphiMETph","sieie1","sipip1","Nph","STg","METS",}},
      //VR 2016
      {"pre_ph165/c_MET300/MT300/STgl600/",{"ph1Pt","MET","MT","absphiMETjet","absphiMETph","STg","absphiMETnJetPh","METdotPh",}},                  
   };

   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         THStack st_mc=hs.getStack(sVar,{"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR"});
   //      THStack st_mc=hs.getStack(sVar,{"diboson","ZNuNuJets","WLNuJets","TTJets","WGToLNuG"});
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
   // divided by bin width
   
   msPresel_vVars={
      {"pre_ph165/"      ,{"ph1Pt","MET","METS","METS_z","METS_long","METS_logx","MT"}}
   };
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         THStack st_mc=hs.getStack(sVar,{"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets_DR"},true);
         gfx::LegendEntries le=hs.getLegendEntries();
         TH1F hEFake(*hs_pix.getHistogram(sVar,"SinglePhoton"));
         hEFake.Scale(cfg.efake.f);
         hist::divideByBinWidth(hEFake);
         hEFake.SetFillColor(cfg.efake.color);
         hEFake.SetFillStyle(1001);
         THStack st_mc_efake=hist::stackPrepend(st_mc,hEFake,"hist");
         THStack &st = sVar.Contains(sPresel+"Ng")
            ? st_mc // do not draw e->gamma fake (not filled anyway)
            : st_mc_efake;
         if (!sVar.Contains(sPresel+"Ng")) le.append(hEFake,cfg.efake.label,"f");
         st.Draw();
  //       auto hists=hs.getHistograms(sVar,{"T5gg","T5Wg","GGM"},true);
         auto hists=hs.getHistograms(sVar,{"TChiWG","T5Wg"},true);
         for (auto const &h: hists) h->Draw("same hist");
         le+=hs.getLegendEntries();
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         saver.save(can,sVar+"_byBW");
      }
   }

   // binned plots
   for (std::string sVar:{"HT","HTgen","jetPt","MET"}){
      TString saveName=TString::Format("pre/binned/%s_",sVar.c_str());
      sVar="pre/"+sVar;
      for (auto sDataset:{"GJets","QCD"}){
         std::vector<TString> samples =cfg.datasets.getDataset(sDataset).getSubsetNames();
         std::reverse(samples.begin(),samples.end());
         auto hists=hs.getHistograms(sVar,samples);
         hists[0]->Draw("axis");
         double max=0;
         for (auto const &h: hists){
            h->Draw("hist same");
            max=std::max(max,h->GetMaximum());
         }
         hists[0]->SetMaximum(max);

         TLegend leg=hs.getLegendEntries().buildLegend(.7,.7);
         leg.Draw();
         saver.save(can,saveName+sDataset);
         THStack st=hs.getStack(sVar,samples);
         st.Draw();
         TLegend leg2=hs.getLegendEntries().buildLegend(.7,.7);
         leg2.Draw();
         saver.save(can,saveName+sDataset+"_st");
      }
   }
   
   // which histograms to save
   // preselection and blind Signal Region: do not save data, do not draw data
   /*
   msPresel_vVars={
      {"pre_ph165/",{"METS_logx","METS_fine","METSHT",}},
      {"pre_ph165/c_S80/MT300/",{"METSHT",}},
   };
   saveHistograms(msPresel_vVars,saver_hist,hs,hs_pix,false);
   msPresel_vVars={
      {"pre_ph165/c_S80/MT300/",{"METS STg"}},
   };
   saveHistograms(msPresel_vVars,saver_hist,h2s);
   // unblinded SR plots
   msPresel_vVars={
      {"pre_ph165/c_S80/MT300/",{"ph1Pt","MET","MT","HT","absphiMETph","sieie1","sipip1","Nph","STg","METS",}},
      {"pre_ph165/c_S80/MT300/eta1/",{"absphiMETph","sieie1","Nph",}},
      {"pre_ph165/c_S80/MT300/STg600/",{"STg",}},
   };
   saveHistograms(msPresel_vVars,saver_hist,hs,hs_pix,true);
   * */
   // control regions: save and plot data later

   //SAVING HISTOGRAMS
   msPresel_vVars={
      //2016 studies
      
      {"effi/",{"trigger_effi_n","trigger_effi_d","trigger_effi_n_cut","trigger_effi_d_cut",
               "trigger_effi_n_match","trigger_effi_d_match","trigger_effi_n_cut_match","trigger_effi_d_cut_match",
         }},
      {"pre_ph165/"      ,{"Ngl","Ngl15","ph1Pt","MET","METS","METS_z","METS_long","METS_logx","METS_fine","HT","METSHT","MT",
                           "relPt2Jets", "DeltaS", "DeltaS1","MT_z","MT_zz","phiMETph","phiMETjet","phiMETnJet","phiMETn2Jet",
                           "phiJetPh","METdotPh","METdotJet","phiMETnJetPh","absphiMETHT","absphiPhHT","vecHTPt","DeltaPhiMETjet100",
         }},
      {"pre_ph165/MT100/"      ,{"METdotPh","ph1Pt",
         }}, 
      {"pre_ph165/dPhiMETphl06/"      ,{"ph1Pt","MET","sieie",
         }}, 
      {"pre_ph165/MT100_l300/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg",
         }},
      {"pre_ph165/c_MET50/MT100/METl200vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/",{"Ngl","Ngl15","ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh","absphiMETnJetPh_JESu","absphiMETnJetPh_JESd",
                                             "HT","METSHT","STg","MT","Nj","relPt2Jets", "DeltaS","DeltaS1",
                                             "absphiMETnJetPh_phoPtl500","absphiMETnJetPh_phoPtl700","absphiMETnJetPh_phoPtg700",
                                             "r9_phoPtl500","r9_phoPtl700","r9_phoPtg700","absphiMETHT","absphiPhHT","vecHTPt","DeltaPhiMETjet100",
         }},

      {"pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETphl02/",{"ph1Pt",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/",{"absphiMETnJetPh","ph1Pt","STg",
         }},
      {"pre_ph165/c_MET100/MT100/",{"r9",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/",{"absphiMETjet","absphiMETnJetPh","STg",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/0b/",{"absphiMETjet","absphiMETnJetPh","STg",
         }},        
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/0l/",{"absphiMETjet","absphiMETnJetPh","STg",
         }},
      {"pre_ph165/c_MET100/MT100/METl300vMTl300/1l/",{"absphiMETjet","absphiMETnJetPh","STg","ph1Pt","MET","MT","relPt2Jets", "DeltaS","DeltaS1",
         }},                  
      {"pre_ph165/c_MET150/MT100/METl300vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT","Nj",
         }},
      {"pre_ph165/c_MET150/MT150/METl400vMTl400/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT","Nj",
         }},
      {"pre_ph165/c_MET100/MT100/METl400vMTl400/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT","Nj",
         }},
         /*
      {"pre_ph165/c_METSHT5/MT100/METSHTl20vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},         
      {"pre_ph165/c_METSHT10/MT100/METSHTl20vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},
      {"pre_ph165/c_METSHT10/MT100/METSHTl15vMTl300/",{"ph1Pt","MET","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETnJetPh",
                                             "HT","METSHT","STg","MT",
         }},*/
      //new SR
      {"pre_ph165/c_MET300/MT300/",{"METSHT","STg", "relPt2Jets", "DeltaS","DeltaS1","METdotPh","Ngl","Ngl15","ph1Pt","absphiMETHT","absphiPhHT","vecHTPt",
         }},
      {"pre_ph165/c_MET300/MT300/STg600/",{"STg","STg_SRbin","METdotPh","Ngl","Ngl15","ph1Pt","r9","r9_Stgg1300","MET","absphiMETph","DeltaPhiMETjet100", "absphiMETnJetPh",
         }},
      {"pre_ph165/c_MET300/MT300/STg600/dPhiMETPhg03/",{"STg","ph1Pt","MET",
         }},      
      {"pre_ph165/c_MET400/MT400/",{"METSHT","STg",
         }},
      // CR
      /*
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETn2Jet","absphiMETnJetPh","METdotPh","METdotJet",
                                             "absphiMETjetJESu","absphiMETjetJESd",
                                             "Nph","Nj","Nele","Nmu","Nl","HT","METSHT","STg","MT","sieie1","sipip1",
                                             "0bL/absphiMETjet","0bM/absphiMETjet",
         }},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETn2Jet","absphiMETnJetPh","METdotPh","METdotJet",
                                                 "absphiMETjetJESu","absphiMETjetJESd",
                                                 "Nph","Nj","Nele","Nmu","Nl","HT","METSHT","STg","MT","sieie1","sipip1",
         }},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/0l/",{"ph1Pt","MET","METS","absphiMETjet","Nph","HT","METSHT","STg","MT"}},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/1l/",{"ph1Pt","MET","METS","absphiMETjet","Nph","HT","METSHT","STg","MT"}},
      {"pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/",{"absphiMETph","sieie1"}},
      {"pre_ph165/c_S30/MT100/Sl100vMTl300/",{"ph1Pt","MET","METS","absphiMETph","absphiMETjet","absphiMETnJet","absphiMETn2Jet","absphiMETnJetPh","METdotPh","METdotJet",}},
 */     // VR
 //     {"pre_ph165/c_S80/MT300/STgl600/",{"ph1Pt","MET","MT","HT","METSHT","absphiMETjet","absphiMETph","sieie1","sipip1","Nph","STg","METS",}},
      //VR 2016
      {"pre_ph165/c_MET300/MT300/STgl600/",{"ph1Pt","MET","MT","absphiMETjet","absphiMETph","STg","absphiMETnJetPh","METdotPh",}},    
   };
   saveHistograms(msPresel_vVars,saver_hist,hs,hs_pix,true);

   // data plots: all saved histograms (which are considered for fitting) and explict CR
   gfx::SplitCan spcan;
;

   // empty bin checking
   can.cd();
   can.SetLogy(false);
   logFile<<"==== Empty bin uncertainty fill-up ====";
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         logFile<<" -- "+sPresel+sVar+" -- ";
         for (auto sDataset:{"WGToLNuG","GJets_DR"}){
            std::vector<TString> samples=cfg.datasets.getDataset(sDataset).getSubsetNames();
            std::reverse(samples.begin(),samples.end());
            auto hists=hs.getHistograms(sPresel+sVar,samples);
            hists[0]->Draw("axis");
            double max=0;
            for (auto const &h: hists){
               h->Draw("hist same");
               h->SetFillColor(h->GetLineColor());
               max=std::max(max,h->GetMaximum());
            }
            hists[0]->SetMaximum(max);

            TLegend leg=hs.getLegendEntries().buildLegend(.7,.7);
            leg.Draw();
            saver.save(can,"emptyBin/"+sPresel+sVar+"_"+sDataset);

            THStack st;
            for (unsigned i=0; i<hists.size(); i++) {
               TString sample=samples[i];
               TH1F &h=*hists[i];

               float meanWeight=h.Integral()/h.GetEntries();
               logFile<<"* "+sample;
               logFile<<TString::Format("mean w=%f/%f=%f",h.Integral(),h.GetEntries(),meanWeight);
               int nBins=h.GetNbinsX();
               int iMax=0;
               int i2Max=0;
               float vMax=h.GetBinCenter(iMax);
               int nEmptyBins=0;
               for (int i=1; i<=nBins; i++) {
                  if (h.GetBinContent(i) != 0) iMax=i;
               }
               if (iMax!=0) {
                  vMax=h.GetBinCenter(iMax);
                  i2Max=h.FindBin(2*vMax);
               } else {
                  vMax=-1;
                  i2Max=0;
               }
               for (int i=1; i<i2Max; i++) {
                  if (h.GetBinContent(i) == 0) {
                     nEmptyBins++;
                  }
               }
               float err=meanWeight/TMath::Sqrt(nEmptyBins);
               for (int i=1; i<i2Max; i++) {
                  if (h.GetBinContent(i) == 0) {
                     h.SetBinError(i,err);
                  }
               }
               logFile<<TString::Format("iMax=%d, vMax=%.1f -> i2Max=%d",iMax,vMax,i2Max);
               logFile<<TString::Format("filling empty %d bins",nEmptyBins);
               h.SetMarkerSize(0);
               h.SetFillStyle(3354);
               h.SetFillColor(kBlack);
               h.DrawClone("e2");
               h.SetFillStyle(0);
               h.Draw("same hist");
               saver.save(can,"emptyBin/"+sPresel+sVar+"_"+sample);
               h.SetFillStyle(1001);
               h.SetFillColor(h.GetLineColor());
               st.Add(&h,"hist");
            }
            TH1F hErr(*(TH1F*)st.GetStack()->Last());
            hErr.SetFillStyle(3354);
            hErr.SetFillColor(kBlack);
            hErr.Draw("axis");
            st.Draw("same");
            hErr.Draw("same e2");
            saver.save(can,"emptyBin/"+sPresel+sVar+"_"+sDataset+"_filled");
         }
      }
   }

   // VR plots
   can.cd();

   // 2d
   
   can.cd();
   can.SetRightMargin (.15);
   can.SetBottomMargin(.22);
   can.SetLogy(false);
   can.SetLogz();
   for (TString const &sVar:h2s.getVariableNames()){
      auto st=h2s.getStack(sVar,{"TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","GJets_DR","QCD"});
      TH2F hBkg(*(TH2F*)st.GetStack()->Last());
      hBkg.Draw("colz");
      gfx::setupZaxis(hBkg);
      can.SetLogx(false);
      if (sVar.Contains("MET MT")){
         drawCR_SR(hBkg,300,300,100,100);
      } else if (sVar=="METS pt"){
         drawCR_SR_VR(hBkg,300,180,50,50,true);
         can.SetLogx();
      } else if (sVar=="Iso/ec"){
         drawFunction(fIso40ec);
      } else if (sVar=="Iso/hc"){
         drawFunction(fIso40hc);
      } else if (sVar=="Iso/tr"){
         drawFunction(fIso40tr);
      }

      TLatex label=gfx::cornerLabel("",3);
      label.DrawLatex(.1,.1,"Background");

      can.RedrawAxis();
      saver.save(can,"2d/"+sVar+"_bkg");
      for (auto sDataset:{"T5Wg","TChiWG","TTGJets","GJets_DR","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","SinglePhoton"}){
         auto hists=h2s.getHistograms(sVar,{sDataset});
         TH2F &h=*hists.at(0);
         h.Draw("colz");
         gfx::setupZaxis(h);
         can.SetLogx(false);
         if (sVar.Contains("MET MT")){
            drawCR_SR(h,300,300,100,100);
         } else if (sVar=="METS pt"){
            drawCR_SR_VR(h,300,180,50,50,true);
            can.SetLogx();
         } else if (sVar=="Iso/ec"){
            drawFunction(fIso40ec);
         } else if (sVar=="Iso/hc"){
            drawFunction(fIso40hc);
         } else if (sVar=="Iso/tr"){
            drawFunction(fIso40tr);
         }
         label.DrawLatex(.1,.1,cfg.datasets.getDataset(sDataset).label);
         can.RedrawAxis();
         saver.save(can,"2d/"+sVar+"_"+sDataset);
      }
   }

   // signal fraction in CR
   float fCR_bkg=0.0;
   float fCR_bkg_e=0.0;
   for (TString const &bkg:{"TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","GJets","QCD"}){
      fCR_bkg+=h2s.getCount("cr",bkg);
      fCR_bkg_e+=TMath::Power(h2s.getCountError("cr",bkg),2);
   }
   fCR_bkg_e=TMath::Sqrt(fCR_bkg_e);
   logFile<<"Signal fraction in CR";
   logFile<<TString::Format("  Bkg %g+-%g",fCR_bkg,fCR_bkg_e);
//   for (TString const &sig:{"T5gg","GGM"}){
   for (TString const &sig:{"T5Wg","TChiWG"}){      
      float const count=h2s.getCount("cr",sig);
      logFile<<TString::Format("  %s %g+-%g (frac: %g)",sig.Data(),count,h2s.getCountError("cr",sig),count/(fCR_bkg+count));
   }

   // Trigger Isolation
   logFile<<"Not passing Iso40";
 //  for (TString const &sig:{"T5gg","GGM"}){
   for (TString const &sig:{"T5Wg","TChiWG"}){
      float const count=h2s.getCount("!Iso40",sig);
      float const total=h2s.getCount("total",sig);
      logFile<<TString::Format("  %s %g/%g (frac: %g)",sig.Data(),count,total,count/total);
   }
}
