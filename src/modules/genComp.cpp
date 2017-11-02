/* Generator comparisons */

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLine.h>
#include <TTreeReader.h>

static Config const &cfg=Config::get();

void run(TString name,std::vector<TString> datasetsToCompare,float HTgen_cut=0.0, float phPt_cut=0.0)
{
   hist::Histograms<TH1F> hs(cfg.datasets.getDatasubsetNames(datasetsToCompare));
   hs.addHist("phPt"  ,";#gamma^{loose} PT;Entries / bin" ,{0,300,500,1000},{20,40,100});
   hs.addHist("ph1Pt" ,";#gamma_{1}^{loose} PT;Events / bin" ,{0,300,500,1500},{20,40,100});
   hs.addHist("STg" ,";#S_{T}^{#gamma} PT;Events / bin" ,{0,300,500,2000},{20,40,100});
   hs.addHist("jetPt" ,";jet PT;EntriesBIN"               ,50,0,1500);
   hs.addHist("jet1Pt",";jet_{1} PT;EventsBIN"            ,50,0,1500);
   hs.addHist("HT"    ,";HT;Events / bin" ,{0,1000,1500,2500},{50,100,250});
   hs.addHist("MET"   ,";MET;Events / bin",{0,300,500,1000},{20,40,100});
   hs.addHist("METSHT",";METSHT;EventsBIN"                ,50,0,50);
 //  hs.addHist("HTgen" ,";HTG;Events / bin",{0,1000,1500,2500},{50,100,250});
   hs.addHist("Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   hs.addHist("Ngm"   ,";N_{#gamma medium};EventsBIN",5 ,-.5, 5-.5);
   hs.addHist("Ngt"   ,";N_{#gamma tight};EventsBIN" ,5 ,-.5, 5-.5);
   hs.addHist("Nggen" ,";N_{#gamma gen};EventsBIN" ,5 ,-.5, 5-.5);
   hs.addHist("Njets" ,";N_{jets};EventsBIN"         ,11,-.5,11-.5);
   hs.addHist("genphPt" ,";#gamma^{gen} PT;Entries / bin" ,{0,300,500,1000},{20,40,100});
   hs.addHist("genph1Pt",";#gamma_{1}^{gen} PT;Events / bin" ,{0,300,500,1000},{20,40,100});

   for (auto const &dss: cfg.datasets.getDatasubsets(datasetsToCompare)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      hs.setCurrentSample(dss.name);

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<tree::Photon>> photons(reader, "photons");
      TTreeReaderValue<std::vector<tree::Muon>>   muons  (reader, "muons");
      TTreeReaderValue<std::vector<tree::Jet>>    jets   (reader, "jets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<tree::MET> MET(reader, "met");
 //     TTreeReaderValue<float> HTgen(reader, "genHt");

      float phoPt;

      std::shared_ptr<WeightCalculator> wCalc;
      if (dss.datasetName=="ZNuNuGJets" || dss.datasetName=="ZGTo2LG"){
         wCalc=std::make_shared<WeightCalculatorZgamma>(phoPt);
      }
      else if (dss.datasetName=="WGToLNuG"){
         wCalc=std::make_shared<WeightCalculatorWgamma>();
      }

      int iEv=0;
      while (reader.Next()){
         iEv++;
         if (iEv>cfg.processFraction*dss.entries) break;
         hs.setFillWeight(*w_pu * *w_mc);

         std::vector<tree::Photon const *> lPho,mPho,tPho;
         for (tree::Photon const &ph: *photons){
            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
            if (ph.hasPixelSeed || fabs(ph.p.Eta())>1.4442) continue;
            lPho.push_back(&ph);
            hs.fill("phPt",ph.p.Pt());
            if (ph.isMedium){
               mPho.push_back(&ph);
               if (ph.isTight){
                  tPho.push_back(&ph);
               }
            }
         }

         if (lPho.size()==0) continue;
         float lPhoPt=lPho[0]->p.Pt();
         phoPt = lPho[0]->p.Pt();
         if (wCalc) hs.setFillWeight(*w_pu * *w_mc *(wCalc->get()));
      
         float met=MET->p.Pt();
         float const MT=phys::M_T(*lPho[0],*MET);
   //      if (lPho[0]->r9<.9) continue;

         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         const float HT=phys::computeHT(cjets);

         hs.fill("HT",HT);
   //      hs.fill("HTgen",*HTgen);
   //      if (*HTgen<HTgen_cut) continue;

         float minDR=std::numeric_limits<float>::max();
         for (tree::Jet const &jet: *jets){
            if (jet.isLoose && jet.p.Pt()>30){
               float const dr=jet.p.DeltaR(lPho[0]->p);
               float const dpt=fabs(jet.p.Pt()-lPhoPt)/lPhoPt;
               if (dr<.1 && dpt<.5) continue;
               minDR=std::min(minDR,dr);
            }
         }
         if (minDR < .5) continue;

         if (lPhoPt<180) continue;
         if (met<100 || MT < 100) continue;

         float STg=met;
         for (auto const &ph: lPho){
               STg+=ph->p.Pt();
         }

         hs.fill("ph1Pt",lPhoPt);

         hs.fill("Ngl",lPho.size());
         hs.fill("Ngm",mPho.size());
         hs.fill("Ngt",tPho.size());

         hs.fill("STg",STg);

         hs.fill("MET",met);

         for (tree::Jet const &jet: cjets){
            hs.fill("jetPt",jet.p.Pt());
         }
         if (cjets.size()>0) hs.fill("jet1Pt",cjets[0].p.Pt());
         hs.fill("Njets",cjets.size());
         hs.fill("METSHT", MET->p.Pt()/TMath::Sqrt(HT));

         int Nggen=0;
         for (auto const &genP: *genParticles){
            if (! (genP.pdgId==22)) continue;  //was always false before! because !genP.pdgId==22 
            Nggen++;
            if (Nggen==1) hs.fill("genph1Pt",genP.p.Pt());
            hs.fill("genphPt",genP.p.Pt());
         }
         hs.fill("Nggen",Nggen);
      } // evt loop

      hs.scaleLumi();
      hs.mergeOverflow();
      file.Close();
   } // datasets

   hs.combineFromSubsamples(datasetsToCompare);

   io::RootFileSaver saver("plots.root",TString::Format("generatorComp%.1f/%s",cfg.processFraction*100,name.Data()));
   gfx::SplitCan spcan;

   for (TString sVar: {"STg","ph1Pt","phPt","HT","MET","Ngl","Nggen","Njets","genphPt","genph1Pt"}){
      spcan.cdUp();
      gPad->SetLogy();
      auto hists=hs.getHistograms(sVar,datasetsToCompare);
      TH1F &hFrame=*(TH1F*)hists[0]->DrawClone("axis");
      Color::reset();
      double max=0;
      for (auto const &h: hists){
         h->Draw("same hist e");
         h->SetLineColor(Color::next());
         h->SetMarkerSize(0);
         max=std::max(max,h->GetMaximum());
      }
      hFrame.SetMaximum(max*1.5);
      gfx::LegendEntries legE = hs.getLegendEntries();
      TLegend leg=legE.buildLegend(.7,.75);
      leg.Draw();

      spcan.cdLow();
      TH1F histRatioFrame=hist::getRatio(*hists[0],*hists[0],"ratio",hist::ONLY1);
      histRatioFrame.Draw("axis");
      histRatioFrame.Draw("same axig");
      max=0;
      for (auto const &h: hists){
         TH1F histRatio=hist::getRatio(*h,*hists[0],"ratio",hist::ONLY1);
         histRatio.DrawClone("same hist e");
         max=std::max(max,histRatio.GetMaximum());
      }
      histRatioFrame.SetMinimum(0);
      histRatioFrame.SetMaximum(max);
      saver.save(spcan,sVar);
   }
}

/*
void WG()
{
   hist::Histograms<TH1F> hs(cfg.datasets.getDatasubsetNames({"WGToLNuG","WGToLNuG_madgraphMLM"}));
   hs.addHist("phPt"  ,";#gamma^{loose} PT;Entries / bin" ,{0,300,500,1000},{20,40,100});
   hs.addHist("ph1Pt" ,";#gamma_{1}^{loose} PT;Events / bin" ,{0,300,500,1000},{20,40,100});
   hs.addHist("jetPt" ,";jet PT;EntriesBIN"               ,50,0,1500);
   hs.addHist("jet1Pt",";jet_{1} PT;EventsBIN"            ,50,0,1500);
   hs.addHist("HT"    ,";HT;Events / bin" ,{0,1000,1500,2500},{50,100,250});
   hs.addHist("MET"   ,";MET;Events / bin",{0,300,500,1000},{20,40,100});
   hs.addHist("METSHT",";METSHT;EventsBIN"                ,50,0,50);
   hs.addHist("HTgen" ,";HTG;Events / bin",{0,1000,1500,2500},{50,100,250});
   hs.addHist("Ngl"   ,";N_{#gamma loose};EventsBIN" ,5 ,-.5, 5-.5);
   hs.addHist("Ngm"   ,";N_{#gamma medium};EventsBIN",5 ,-.5, 5-.5);
   hs.addHist("Ngt"   ,";N_{#gamma tight};EventsBIN" ,5 ,-.5, 5-.5);
   hs.addHist("Njets" ,";N_{jets};EventsBIN"         ,11,-.5,11-.5);
   for (auto const &dss: cfg.datasets.getDatasubsets({"WGToLNuG","WGToLNuG_madgraphMLM"})){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      hs.setCurrentSample(dss.name);

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<tree::Photon>> photons(reader, "photons");
      TTreeReaderValue<std::vector<tree::Muon>>   muons  (reader, "muons");
      TTreeReaderValue<std::vector<tree::Jet>>    jets   (reader, "jets");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<float> HTgen(reader, "genHt");

      int iEv=0;
      while (reader.Next()){
         iEv++;
         if (iEv>cfg.processFraction*dss.entries) break;
         hs.setFillWeight(*w_pu * *w_mc);

         std::vector<tree::Photon const *> lPho,mPho,tPho;
         for (tree::Photon const &ph: *photons){
            if (ph.hasPixelSeed || fabs(ph.p.Eta())>1.4442) continue;
            lPho.push_back(&ph);
            hs.fill("phPt",ph.p.Pt());
            if (ph.isMedium){
               mPho.push_back(&ph);
               if (ph.isTight){
                  tPho.push_back(&ph);
               }
            }
         }

         if (lPho.size()==0) continue;
         float lPhoPt=lPho[0]->p.Pt();
         float met=MET->p.Pt();
         if (lPho[0]->r9<.9) continue;
         if (lPhoPt<50) continue;
         if (met<100) continue;

         float minDR=std::numeric_limits<float>::max();
         for (tree::Jet const &jet: *jets){
            if (jet.isLoose && jet.p.Pt()>30){
               float const dr=jet.p.DeltaR(lPho[0]->p);
               float const dpt=fabs(jet.p.Pt()-lPhoPt)/lPhoPt;
               if (dr<.1 && dpt<.5) continue;
               minDR=std::min(minDR,dr);
            }
         }
         hs.fill("minDR",minDR);
         if (minDR < .5) continue;

         hs.fill("Ngl",lPho.size());
         hs.fill("Ngm",mPho.size());
         hs.fill("Ngt",tPho.size());
         hs.fill("ph1Pt",lPhoPt);

         hs.fill("MET",met);

         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         for (tree::Jet const &jet: cjets){
            hs.fill("jetPt",jet.p.Pt());
         }
         if (cjets.size()>0) hs.fill("jet1Pt",cjets[0].p.Pt());
         hs.fill("Njets",cjets.size());
         const float HT=phys::computeHT(cjets);
         hs.fill("HT",HT);
         hs.fill("METSHT", MET->p.Pt()/TMath::Sqrt(HT));
         hs.fill("HTgen",*HTgen);
      } // evt loop

      hs.scaleLumi();
      file.Close();
   } // datasets

   io::RootFileSaver saver("plots.root",TString::Format("generatorComp%.1f",cfg.processFraction*100));
   gfx::SplitCan spcan;

   for (TString sVar: {"ph1Pt","phPt","HT","MET","HTgen","Ngl","Njets"}){
      spcan.cdUp();
      gPad->SetLogy();
      gfx::LegendEntries legE;
      auto hists=hs.getHistograms(sVar,{"WGToLNuG-amcatnloFXFX","WGToLNuG-madgraphMLM"});
      legE.append(*hists[1], "MG MLM" ,"pe");
      hists[0]->SetFillStyle(1001);
      hists[0]->SetFillColor(kGray+1);
      hists[0]->SetLineColor(kBlack);
      hists[0]->SetMarkerSize(0);
      hists[0]->SetLineWidth(1);
      TH1 *hEbar=hists[0]->DrawCopy("e2");
      legE.append(*hEbar, "aMC@NLO","fl");
      TH1F hDraw(*hists[0]);
      hDraw.SetFillStyle(0);
      for (int i=0; i<=hDraw.GetNbinsX()+1; i++) hDraw.SetBinError(i,1e-10); // setting to zero for all bins makes histogram errorless ("E" option has no effect)
      hDraw.Draw("e same");
      hists[1]->Draw("e same");
      hists[1]->SetLineColor(kBlack);
      hists[1]->SetMarkerStyle(4);
      TLegend leg=legE.buildLegend(.7,.75);
      leg.SetHeader("W(#rightarrowl#nu)+#gamma");
      leg.Draw();

      TH1F histRatio=hist::getRatio(*hists[0],*hists[1]);
      spcan.cdLow();
      histRatio.Draw("axis"); // draw first for axes
      histRatio.Draw("same axig");
      TGraphErrors grRatio(&histRatio);
      grRatio.Draw("SAME PE");
      float ratio=hists[0]->Integral()/hists[1]->Integral();
      TLine l;
      l.DrawLine(histRatio.GetXaxis()->GetXmin(),ratio,histRatio.GetXaxis()->GetXmax(),ratio);
      gPad->RedrawAxis();
      saver.save(spcan,"WG_"+sVar);
   }
}
*/

extern "C"
void run()
{
//   run("GJets",{"GJets","GJets_DR"},0.0,0.0);
   run("ZNuG", {"ZNuNuGJets" ,"ZGTo2NuG"} ,0.0,0.0);
   run("WLNuG", {"WGToLNuG" ,"WGToLNuG_SUS"} ,0.0,0.0);   

   // WG();
}
