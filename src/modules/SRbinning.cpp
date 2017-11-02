#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

#include <TFractionFitter.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TMath.h>
#include <TRandom3.h>

#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooEllipse.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

Config const &cfg=Config::get();
io::RootFileSaver saver("plots.root","SRbinning");

template <class H>
H combineHists(io::RootFileReader const& histReader, TString sPresel, TString sVar, std::vector<TString> sampleNames)
{
   assert(sampleNames.size()>0);
   H const *hFirst=histReader.read<H>(sPresel+sVar+"/"+sampleNames[0]);
   H hComb(*hFirst);
   hComb.SetLineColor(hFirst->GetLineColor());
   hComb.SetFillColor(hFirst->GetLineColor());
   for (unsigned i=1; i<sampleNames.size(); i++){
      hComb.Add(histReader.read<H>(sPresel+sVar+"/"+sampleNames[i]));
   }
   return hComb;
}

TH2F getCumulative(TH2F const& h) {
   TH2F hi(h);
   // overflow bins:
   int iXof=h.GetNbinsX()+1;
   int iYof=h.GetNbinsY()+1;
   for (int x=h.GetNbinsX(); x>=0; x--) {
      for (int y=h.GetNbinsY(); y>=0; y--) {
         // integrate with overflow
         hi.SetBinContent(x,y,h.Integral(x,iXof,y,iYof));
      }
   }
   return hi;
}

void setupHistForStack(TH1 &h){
   h.SetFillStyle(1001);
   h.SetFillColor(h.GetLineColor());
   h.SetLineColor(kBlack);
   h.SetLineWidth(1);
}

void optim1d(TString sSelection,TString sVar,std::vector<float> upperBinEdges)
{
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   TString saveName=sSelection+sVar;

   int stage=upperBinEdges.size();
   float upperBinEdge=-1;
   if (!upperBinEdges.empty()) {
      // the stage this call is computing
      upperBinEdge=upperBinEdges.back();
      upperBinEdges.pop_back();
   }
   if (stage!=0) {
      saveName+=TString::Format("_%d",stage);
      // compute the next stage
      optim1d(sSelection,sVar,upperBinEdges);
   }

   std::vector<TString> vsVgComponents{"WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"};
   std::vector<TString> vsGJComponents{"QCD","GJets"};

   TH1F hVg=combineHists<TH1F>(histReader,sSelection,sVar,vsVgComponents);
   TH1F hGJ=combineHists<TH1F>(histReader,sSelection,sVar,vsGJComponents);
   hVg.Scale(cfg.trigger_eff_Ph*cfg.sf.Vg);
   hGJ.Scale(cfg.trigger_eff_Ph*cfg.sf.GJ);

   hVg.SetLineColor(kOrange-3);
   hGJ.SetLineColor(kGray);
   setupHistForStack(hVg);
   setupHistForStack(hGJ);

   std::map<TString,TH1F> fixHists;
   for (TString sSample:{"efake","TTJets","TTGJets","diboson","GGM","T5gg","T5Wg"}){
      fixHists[sSample]=*(TH1F*)histReader.read<TH1F>(sSelection+sVar+"/"+sSample);
      TH1 &h=fixHists[sSample];
      h.Scale(cfg.trigger_eff_Ph);
      setupHistForStack(h);
   }
   fixHists["efake"].Scale(1./cfg.trigger_eff_Ph); // efake is data-driven

   fixHists["TTcomb"]=fixHists["TTJets"];
   fixHists["TTcomb"].Add(&fixHists["TTGJets"]);
   THStack st;
   st.Add(&fixHists["efake"],"hist");
   st.Add(&fixHists["TTcomb"],"hist");
   st.Add(&fixHists["diboson"],"hist");

   st.Add(&hGJ,"hist");
   st.Add(&hVg,"hist");

   std::map<TString,TH1F> hists;
   hists["bkg"]=TH1F(*(TH1F*)st.GetStack()->Last());
   hists["bkg"].SetFillStyle(0);
   hists["bkg"].SetLineWidth(2);

   TCanvas can;
   gfx::LegendEntries le;
   hists["bkg"].Draw("hist");
   le.append(hists["bkg"],"Background","l");
   for (TString sSample:{"GGM","T5gg","T5Wg"}){
      TH1F &h=hists[sSample]=TH1F(fixHists[sSample]);
      h.SetFillStyle(0);
      h.SetLineWidth(2);
      h.SetLineColor(h.GetFillColor());
      h.Draw("same hist");
      le.append(h,sSample,"l");
   }

   // integrals
   for (TString sSample:{"bkg","GGM","T5gg","T5Wg"}){
      TH1F &hInt=hists[sSample+"_int"]=TH1F(hists[sSample]);
      TH1F &h=hists[sSample];
      hInt.Reset();
      hInt.SetLineStyle(7);
      float integral=0;
      float statErr=0;
      int iUpperBinEdge=upperBinEdge<0 ? h.GetNbinsX()+1 : h.FindBin(upperBinEdge);
      for (int i=iUpperBinEdge; i>=0; i--) {
         integral+=h.GetBinContent(i);
         statErr+=h.GetBinError(i);
         hInt.SetBinContent(i,integral);
         hInt.SetBinError(i,TMath::Sqrt(statErr));
      }
      hInt.Draw("same hist");
   }
   le.buildLegend(.49,.25,-1,-1,2).DrawClone();
   gfx::LegendEntries le_lines;
   le_lines.append(hists["bkg"],"Events / bin","l");
   le_lines.append(hists["bkg_int"],"#sum^{}_{#it{S}>x} Events","l");
   le_lines.buildLegend(.75,.8).DrawClone();
   hists["bkg"].SetMaximum(hists["bkg_int"].GetMaximum());
   can.SetLogy(true);
   saver.save(can,saveName+"_int");

   hists["bkg"].Draw("axis");
   for (TString sSample:{"GGM","T5gg","T5Wg"}){
      TH1F &hS=hists[sSample+"_int"];
      TH1F &hB=hists["bkg_int"];

      TH1F &hSB=hists[sSample+"_s/b"]=TH1F(hS);
      TH1F &hSSB=hists[sSample+"_s/sb"]=TH1F(hS);
      TH1F &hSSBstat=hists[sSample+"_s/sbstat"]=TH1F(hS);
      hSB.Reset();
      hSSB.Reset();
      hSSB.SetLineStyle(1);
      hSSBstat.Reset();
      hSSBstat.SetLineStyle(3);
      float s=0, b=0, estat=0;
      for (int i=hS.GetNbinsX()+1; i>=0; i--) {
         s=hS.GetBinContent(i);
         b=hB.GetBinContent(i);
         estat=hB.GetBinError(i);
         hSB.SetBinContent(i,(b>0 ? s/TMath::Sqrt(b) : 0));
         hSSB.SetBinContent(i,(b>0 ? s/TMath::Sqrt(s+b) : 0));
         hSSBstat.SetBinContent(i,(b>0 ? s/TMath::Sqrt(s+b+estat*estat) : 0));
      }
      hSB.Draw("same hist");
      hSSB.Draw("same hist");
      hSSBstat.Draw("same hist");
   }
   hists["bkg"].SetMaximum(.001);
   hist::setMaximum(hists["bkg"],{hists["T5gg_s/b"],hists["GGM_s/sb"]});
   hists["bkg"].GetYaxis()->SetTitle("");
   le.pop_front(); // remove bkg
   le.buildLegend(.35,.75).DrawClone();
   le_lines.clear();
   le_lines.append(2,7,kBlack,"s/b^{1/2}","l");
   le_lines.append(2,1,kBlack,"s/(s+b)^{1/2}","l");
   le_lines.append(2,3,kBlack,"s/(s+b+#sigma_{b stat}^{2})^{1/2}","l");
   le_lines.buildLegend(.35,.75,.65).DrawClone();
   can.SetLogy(false);
   saver.save(can,saveName);
}

void optim2d(TString sSelection,TString sVar)
{
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   TString saveName=sSelection+"2d/"+sVar;
   std::vector<TString> vsVgComponents{"WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"};
   std::vector<TString> vsGJComponents{"QCD","GJets"};

   TH2F hVg=combineHists<TH2F>(histReader,sSelection,sVar,vsVgComponents);
   TH2F hGJ=combineHists<TH2F>(histReader,sSelection,sVar,vsGJComponents);
   hVg.Scale(cfg.trigger_eff_Ph*cfg.sf.Vg);
   hGJ.Scale(cfg.trigger_eff_Ph*cfg.sf.GJ);

   std::map<TString,TH2F> fixHists;
   // efake not calculated in 2d
   for (TString sSample:{/*"efake",*/"TTJets","TTGJets","diboson","GGM","T5gg","T5Wg"}){
      fixHists[sSample]=*(TH2F*)histReader.read<TH2F>(sSelection+sVar+"/"+sSample);
      TH2 &h=fixHists[sSample];
      h.Scale(cfg.trigger_eff_Ph);
      setupHistForStack(h);
   }
   fixHists["efake"].Scale(1./cfg.trigger_eff_Ph); // efake is data-driven
   THStack st;

   st.Add(&fixHists["TTGJets"],"hist");
   st.Add(&fixHists["TTJets"],"hist");
   st.Add(&fixHists["diboson"],"hist");

   st.Add(&hGJ,"hist");
   st.Add(&hVg,"hist");

   std::map<TString,TH2F> hists;
   hists["bkg"]=TH2F(*(TH2F*)st.GetStack()->Last());
   for (TString sSample:{"GGM","T5gg","T5Wg"}) {
      hists[sSample]=TH2F(fixHists[sSample]);
   }

   TCanvas can;
   can.SetRightMargin (.15);
   can.SetBottomMargin(.22);
   can.SetLogz();

   TLatex latex=gfx::cornerLabel("",3);
   for (TString sSample:{"bkg","GGM","T5gg","T5Wg"}) {
      TH2F &h=hists[sSample];
      h.Draw("colz");
      gfx::setupZaxis(h);
      latex.DrawLatex(.02,.02,sSample);
      saver.save(can,saveName+"_"+sSample);

      // Integrals
      TH2F &hi=hists[sSample+"_int"]=getCumulative(hists[sSample]);
      hi.Draw("colz");
      hi.GetZaxis()->SetTitle("#lower[-.3]{#scale[.5]{#sum}}#lower[.1]{#scale[.5]{#it{S}>x,S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}>y}} Events");
      gfx::setupZaxis(hi);
      latex.DrawLatex(.02,.02,sSample);
      saver.save(can,saveName+"_"+sSample+"_int");
   }

   can.SetLogz(false);
   for (TString sSample:{"GGM","T5gg","T5Wg"}) {
      TH2F &hS=hists[sSample+"_int"];
      TH2F &hB=hists["bkg_int"];

      TH2F &hSB=hists[sSample+"_s/b"]=TH2F(hS);
      TH2F &hSSB=hists[sSample+"_s/sb"]=TH2F(hS);
      TH2F &hSSBstat=hists[sSample+"_s/sbstat"]=TH2F(hS);

      float s,b,estat;
      for (int x=0; x<=hB.GetNbinsX(); x++){
         for (int y=0; y<=hB.GetNbinsY(); y++){
            s=hS.GetBinContent(x,y);
            b=hB.GetBinContent(x,y);
            estat=hB.GetBinError(x,y);
            hSB.SetBinContent     (x,y,(b>0 ? s/TMath::Sqrt(b) : 0));
            hSSB.SetBinContent    (x,y,(b>0 ? s/TMath::Sqrt(s+b) : 0));
            hSSBstat.SetBinContent(x,y,(b>0 ? s/TMath::Sqrt(s+b+estat*estat) : 0));
         }
      }
      hSB.Draw("colz");
      hSB.GetZaxis()->SetTitle("s/b^{1/2}");
      gfx::setupZaxis(hSB,false);
      latex.DrawLatex(.02,.02,sSample);
      saver.save(can,saveName+"_"+sSample+"_SB");

      hSSB.Draw("colz");
      hSSB.GetZaxis()->SetTitle("s/(s+b)^{1/2}");
      gfx::setupZaxis(hSSB,false);
      latex.DrawLatex(.02,.02,sSample);
      saver.save(can,saveName+"_"+sSample+"_SSB");

      hSSBstat.Draw("colz");
      hSSBstat.GetZaxis()->SetTitle("s/(s+b+#sigma_{b stat}^{2})^{1/2}");
      gfx::setupZaxis(hSSBstat,false);
      latex.DrawLatex(.02,.02,sSample);
      saver.save(can,saveName+"_"+sSample+"_SSBstat");
   }
}

extern "C"
void run()
{
   optim1d("pre_ph165/c_S80/MT300/","METS",{750,400,200});
   optim1d("pre_ph165/c_S80/MT300/","STg",{1000,800,600});

   optim2d("pre_ph165/c_S80/MT300/","METS STg");
}
