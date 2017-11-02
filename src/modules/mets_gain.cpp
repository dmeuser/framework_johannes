#include "Config.hpp"
#include "tools/io.hpp"

#include <TMath.h>

Config const &cfg=Config::get();

TH1F cumulative(TH1F const &h)
{
   TH1F hInt(h);
   float integral=0;
   float statErr=0;
   for (int i=h.GetNbinsX()+1; i>=0; i--) {
      integral+=h.GetBinContent(i);
      statErr+=h.GetBinError(i)*h.GetBinError(i);
      hInt.SetBinContent(i,integral);
      hInt.SetBinError(i,TMath::Sqrt(statErr));
   }
   return hInt;
}

extern "C"
void run()
{
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   std::map<TString,TH1F> hIntB;
   hIntB["METS_fine"]=TH1F(*histReader.read<TH1F>("pre_ph165/METS_fine/GJets"));
   hIntB["METS_fine"].Reset();
   hIntB["METSHT"]=TH1F(*histReader.read<TH1F>("pre_ph165/METSHT/GJets"));
   hIntB["METSHT"].Reset();
   std::map<TString,TH1F> hIntS;
   hIntS["METS_fine"]=TH1F(*histReader.read<TH1F>("pre_ph165/METS_fine/GGM"));
   hIntS["METSHT"]=TH1F(*histReader.read<TH1F>("pre_ph165/METSHT/GGM"));

   for (TString sVar: {"METS_fine","METSHT"}) {
      for (TString sSample: {"diboson","ZNuNuJets","WLNuJets","TTJets","TTGJets","ZGTo2LG","ZNuNuGJets","WGToLNuG","QCD","GJets","efake"}) {
         TH1F h(*histReader.read<TH1F>("pre_ph165/"+sVar+"/"+sSample));
         hIntB[sVar].Add(&h);
      }
      hIntB[sVar]=cumulative(hIntB[sVar]);
      hIntS[sVar]=cumulative(hIntS[sVar]);
   }


   std::map<TString,TH1F> hIntSB;
   io::RootFileSaver saver("plots.root","mets_gain");
   TCanvas can;
   for (TString sVar: {"METS_fine","METSHT"}) {
      TH1F &hB=hIntB[sVar];
      hB.Draw("hist");
      hB.GetYaxis()->SetTitle("cumulative");
      hB.SetLineColor(kBlack);
      TH1F &hS=hIntS[sVar];
      hS.Draw("hist same");
      saver.save(can,sVar+"_int");
      TH1F &hSB=hIntSB[sVar]=TH1F(hS);
      float s=0, b=0;
      for (int i=hS.GetNbinsX()+1; i>=0; i--) {
         s=hS.GetBinContent(i);
         b=hB.GetBinContent(i);
         hSB.SetBinContent(i,(b>0 ? s/b:0));
      }
      hSB.Draw("hist");
      hSB.GetYaxis()->SetTitle("S/B");
      saver.save(can,sVar+"_SB");
      hSB.SetLineColor(kBlue);
   }
   // somehow crashes. when map is deleted??
   // maybe ownership problem with root file
   // but plots are stored ok.
   debug <<"not crashed yet";
}
