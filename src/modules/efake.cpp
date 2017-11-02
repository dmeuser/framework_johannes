/* electron fake studies (closure test) */

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>

Config const &cfg=Config::get();

// Macro to add histograms to all, the normal Histograms instance
// and the "pixel" Histograms and "gen-electron" Histograms
#define ADD_HIST(args...) do { \
      hs_notPix.addHist(args); \
      hs_pix.addHist(args);    \
      hs_genE.addHist(args);   \
   } while (0)

enum pass_t {pass_normal,pass_pixel};

extern "C"
void run()
{
   std::vector<TString> datasetsToUse={"TTJets","diboson","WLNuJets"};
   hist::Histograms<TH1F> hs_notPix(cfg.datasets.getDatasubsetNames(datasetsToUse));
   hist::Histograms<TH1F> hs_pix(cfg.datasets.getDatasubsetNames(datasetsToUse));
   hist::Histograms<TH1F> hs_genE(cfg.datasets.getDatasubsetNames(datasetsToUse));

   ADD_HIST("pre/ph1Pt",";#gamma_{1} PT;Events / bin",{75,825},{25});
   ADD_HIST("pre/MET"  ,";MET;Events / bin"          ,{0,625},{25});
//   ADD_HIST("pre/METS" ,";METSIG;Events / bin"       ,{0,20,30,60,62},{5,10,15,2});
   ADD_HIST("pre/MT"    ,";MT(#gamma_{1},%MET);EventsBIN"   ,{0,625},{25});
   ADD_HIST("pre/STg"   ,";STg;EventsBIN"           ,{75,1025},{25});

   ADD_HIST("loose/ph1Pt",";#gamma_{1} PT;Events / bin",{75,825},{25});
   ADD_HIST("loose/MET"  ,";MET;Events / bin"         ,{25,625},{25});
//   ADD_HIST("loose/METS" ,";METSIG;Events / bin"       ,{0,20,30,60,62},{5,10,15,2});
   ADD_HIST("loose/MT"    ,";MT(#gamma_{1},%MET);EventsBIN"  ,{25,625},{25});
   ADD_HIST("loose/STg"   ,";STg;EventsBIN"           ,{100,1025},{25});

   ADD_HIST("CR/ph1Pt",";#gamma_{1} PT;Events / bin",{0,1100},{100});
   ADD_HIST("CR/MET"  ,";MET;Events / bin"          ,{0,900},{100});
//   ADD_HIST("CR/METS" ,";METSIG;Events / bin"       ,{0,20,30,60,62},{5,10,15,2});
   ADD_HIST("CR/MT"    ,";MT(#gamma_{1},%MET);EventsBIN"   ,{0,700},{100});
   ADD_HIST("CR/STg"   ,";STg;EventsBIN"           ,{100,1700},{100});

  

   for (auto const &dss: cfg.datasets.getDatasubsets(datasetsToUse)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      hs_notPix.setCurrentSample(dss.name);
      hs_pix.setCurrentSample(dss.name);
      hs_genE.setCurrentSample(dss.name);

      // bool const isData=dss.isData;

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<tree::Photon>> photons(reader, "photons");
      TTreeReaderValue<std::vector<tree::Muon>>   muons  (reader, "muons");
      TTreeReaderValue<std::vector<tree::Jet>>    jets   (reader, "jets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<float> HTgen(reader, "genHt");
      TTreeReaderValue<bool> trigger(reader, "HLT_Photon165_HE10_v");

      int iEv=0;
      while (reader.Next()){
         iEv++;
         if (iEv>cfg.processFraction*dss.entries) break;

         float fEventWeight=*w_pu * *w_mc;
         hs_notPix.setFillWeight(fEventWeight);
         hs_pix.setFillWeight(fEventWeight);
         hs_genE.setFillWeight(fEventWeight);

         std::vector<tree::Photon const *> lPho,mPho,tPho,lPixPho;
         for (tree::Photon const &ph: *photons){
            // (ph.hasPixelSeed?hs_pix:hs_notPix).fill("pre/phEta",ph.p.Eta());

            if (fabs(ph.p.Eta())>1.4442) continue;
            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
            if (ph.hasPixelSeed){
               // hs_pix.fill("pre/phPt",ph.p.Pt());
               lPixPho.push_back(&ph);
            } else {
               // hs_notPix.fill("pre/phPt",ph.p.Pt());
               lPho.push_back(&ph);
               if (ph.isMedium) mPho.push_back(&ph);
               if (ph.isTight)  tPho.push_back(&ph);
            }
         }

         for (pass_t pass:{pass_normal,pass_pixel}){
            fEventWeight=*w_pu * *w_mc;
            hs_notPix.setFillWeight(fEventWeight);
            hs_pix.setFillWeight(fEventWeight);
            hs_genE.setFillWeight(fEventWeight);

            if (pass==pass_normal){
               if (lPho.empty()) continue;
            } else if (pass==pass_pixel) {
               // if (!isData) continue; // MC pixel-pass for closure
               if (lPixPho.empty()) continue;
               // pixel-object needs to be leading EM object
               if (!lPho.empty() && lPho[0]->p.Pt()>lPixPho[0]->p.Pt()) continue;
            }
            // find e->photon fakes in MC (in data genParticles is empty, so nothing is rejected)
            bool const isGenE=(pass==pass_normal && phys::matchesGen(*lPho[0],*genParticles,11,0.1,0.5));

            std::vector<tree::Photon const*> const &pho = (pass==pass_pixel) ? lPixPho : lPho;
            hist::Histograms<TH1F> &hs = (pass==pass_pixel) ? hs_pix : (!isGenE ? hs_notPix : hs_genE);

            float phoPt=pho[0]->p.Pt();
            float met=MET->p.Pt();

            float minDR=std::numeric_limits<float>::max();
            for (tree::Jet const &jet: *jets){
               if (jet.isLoose && jet.p.Pt()>30){
                  float const dr=jet.p.DeltaR(pho[0]->p);
                  float const dpt=fabs(jet.p.Pt()-phoPt)/phoPt;
                  if (dr<.1 && dpt<.5) continue;
                  minDR=std::min(minDR,dr);
               }
            }
            // hs.fill("pre/minDR",minDR);
            if (minDR < .5) continue;
            if (phoPt < 100) continue;

            // "pre" selection

            float STg=met;
            for (auto const &ph: pho){
               STg+=ph->p.Pt();
            }
            
            float const MT=phys::M_T(*pho[0],*MET);

            hs.fill("pre/ph1Pt",phoPt);
            hs.fill("pre/MET",met);
    //        hs.fill("pre/METS",MET->sig);
            hs.fill("pre/MT",MT);
            hs.fill("pre/STg",STg);

            if (met < 50 || MT < 50) continue;

            hs.fill("loose/ph1Pt",phoPt);
            hs.fill("loose/MET",met);
      //      hs.fill("loose/METS",MET->sig);
            hs.fill("loose/MT",MT);
            hs.fill("loose/STg",STg);

            if (met < 100 || MT < 100) continue;
            if (met > 300 && MT > 300) continue;
            if (phoPt < 180) continue;           

            hs.fill("CR/ph1Pt",phoPt);
            hs.fill("CR/MET",met);
      //      hs.fill("CR/METS",MET->sig);
            hs.fill("CR/MT",MT);
            hs.fill("CR/STg",STg);                                       
         } // normal/pixel pass loop
      } // evt loop

      hs_notPix.scaleLumi();
      hs_notPix.mergeOverflow();
      hs_pix.scaleLumi();
      hs_pix.mergeOverflow();
      hs_genE.scaleLumi();
      hs_genE.mergeOverflow();
      file.Close();
   } // datasets

   hs_notPix.combineFromSubsamples({"diboson"});
   hs_pix.combineFromSubsamples({"diboson"});
   hs_genE.combineFromSubsamples({"diboson"});
   hs_notPix.combineFromSubsamples({"WLNuJets"});
   hs_pix.combineFromSubsamples({"WLNuJets"});
   hs_genE.combineFromSubsamples({"WLNuJets"});   
   io::RootFileSaver saver("plots.root",TString::Format("efake%.1f",cfg.processFraction*100));
   gfx::SplitCan can;
   can.pU_.SetLogy();
   for (std::string sVar:{"pre/MT","pre/STg","pre/ph1Pt","pre/MET","loose/MT","loose/STg","loose/ph1Pt","loose/MET","CR/MT","CR/STg","CR/ph1Pt","CR/MET"}){
      gfx::LegendEntries le;
      TH1F h_genE(*(TH1F*)hs_genE.getStack(sVar,datasetsToUse).GetStack()->Last());
      TH1F h_pix(*(TH1F*)hs_pix.getStack(sVar,datasetsToUse).GetStack()->Last());
      TH1F h_pred(h_pix);
      h_pred.Scale(cfg.efake.f_mc/(1.-cfg.efake.f_mc));

      double integral, e_integral;
      integral=h_pred.IntegralAndError(0,h_pred.GetNbinsX()+1,e_integral);
      io::log<<TString::Format("Pred.  : %.2f+-%.2f",integral, e_integral);
      integral=h_genE.IntegralAndError(0,h_genE.GetNbinsX()+1,e_integral);
      io::log<<TString::Format("Matched: %.2f+-%.2f",integral, e_integral);

      hist::divideByBinWidth(h_genE);
      hist::divideByBinWidth(h_pred);
      h_pred.SetFillStyle(1001);
      h_pred.SetFillColor(kGray+1);
      h_pred.SetLineColor(kBlack);
      h_pred.SetMarkerSize(0);
      h_pred.SetLineWidth(1);
      can.cdUp();
      TH1 *h_pred_err=h_pred.DrawCopy("e2");
      TH1F h_pred_syst_err(h_pred);

      h_pred_err->SetLineColor(kWhite);
      h_pred_err->SetFillColor(kMagenta+4);
      h_pred_err->SetFillStyle(3345);          
      h_pred.SetFillStyle(0);

      for (int i=0; i <= h_pred_syst_err.GetNbinsX(); i++){
         h_pred_syst_err.SetBinError(i, (h_pred_syst_err.GetBinContent(i)*0.3));
      }
      h_pred_syst_err.SetLineColor(kWhite);
      h_pred_syst_err.SetFillColor(kRed+1);
      h_pred_syst_err.SetFillStyle(3354);
      h_pred_syst_err.Draw("same e2");
        
      for (int i=0; i<=h_pred.GetNbinsX()+1; i++) h_pred.SetBinError(i,1e-10); // setting to zero for all bins makes histogram errorless ("E" option has no effect)
      hist::setMaximum(h_pred_syst_err,{h_genE});
      hist::setMinimum(h_pred_syst_err,{h_genE},0.9,false);
      h_pred.Draw("same");
      TGraphErrors gr_genE(&h_genE); // to show errors of points outside of range
      gr_genE.Draw("p0");
            le.append(h_genE, "gen matched e#rightarrow#gamma","pe");
      le.append(h_pred, "predicted e#rightarrow#gamma","l");
      le.append(*h_pred_err, "#sigma_{stat, pred}","f");
      le.append(h_pred_syst_err, "#sigma_{syst, pred}","f");     
      TLegend leg=le.buildLegend(.55,.7);
      leg.Draw();

      can.cdLow();
      TH1F hRatio=hist::getRatio(h_genE,h_pred,"ratio",hist::ONLY1);
      TGraphErrors grRatioStat=hist::getRatioGraph(*(TH1F*)h_pred_err,h_pred,"ratio",hist::ONLY1);
      TGraphErrors grRatioSyst=hist::getRatioGraph((TH1F)h_pred_syst_err,h_pred,"ratio",hist::ONLY1);

      hRatio.SetMaximum(1.9);
      hRatio.SetMinimum(0.1);
      hRatio.Draw("axis");
      grRatioStat.SetFillStyle(3345);
      grRatioStat.SetFillColor(kMagenta+4);
      
      grRatioSyst.SetFillStyle(3354);
      grRatioSyst.SetFillColor(kRed+1);
         
      gfx::scaleXerrors(grRatioStat,.9);
      gfx::scaleXerrors(grRatioSyst,.5);
           
      grRatioStat.Draw("2");
      grRatioSyst.Draw("2");     
      hRatio.Draw("same axig");
      hRatio.Draw("same pe0");

      saver.save(can,sVar);
   }
}
