#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"


#include <TGaxis.h>
#include <TMarker.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TStyle.h>

#include "TNamed.h"
#include "RooChi2Var.h"
#include "RooExtendPdf.h"
#include "RooMinuit.h"
#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooEllipse.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"

Config const &cfg=Config::get();
io::RootFileSaver saver("plots.root","fit_weighting");
TRandom3 g_rand(1111);

Color_t color1=kOrange+7;
Color_t color2=kAzure-6;
Color_t fillcolor1=kOrange-9;
Color_t fillcolor1_2=kOrange-3;
Color_t fillcolor2=kAzure+10;
Color_t fillcolor2_2=kAzure-3;

class Rebinner
{
public:
   // if iRebin==0: use vectors
   Rebinner(int iRebin)
      : iRebin_(iRebin)
      {}
   Rebinner(std::vector<float> rebinEdges,std::vector<float> rebinWidths)
      : iRebin_(0)
      {
         newBinEdges_=hist::getBinVector(rebinEdges,rebinWidths);
      }
   TH1F operator()(TH1F const &h) const {
      if (iRebin_!=0){
         TH1F hNew(h);
         hNew.Rebin(iRebin_);
         return hNew;
      } else {
         return hist::rebinned(h,newBinEdges_);
      }
   }
   int numberOfBins(){
      if (iRebin_==0) return newBinEdges_.size()-1;
      return 0;
   }
private:
   int iRebin_;
   std::vector<double> newBinEdges_;
};

struct FitResult
{
   TString selection,distribution,distLabel;
   int numberOfBins;
   float r1,r2,e1,e2; // results and errors
   int weight_index;
   float rho; // correlation
};
bool cmpSelection   (const FitResult &a,const FitResult &b) { return a.selection.CompareTo(b.selection)      >0; }
bool cmpDistribution(const FitResult &a,const FitResult &b) { return a.distribution.CompareTo(b.distribution)>0; }
bool cmpNoBins      (const FitResult &a,const FitResult &b) { return a.numberOfBins>b.numberOfBins; }

std::vector<FitResult> g_fitResults;

TH1F combineHists(std::map<TString,TH1F> singleHists,std::vector<TString> sampleNames)
{
   assert(sampleNames.size()>0);
   TH1F hComb(singleHists[sampleNames[0]]);
   hComb.Reset();
   hComb.SetLineColor(singleHists[sampleNames[0]].GetLineColor());
   hComb.SetFillColor(singleHists[sampleNames[0]].GetLineColor());
   for (unsigned i=0; i<sampleNames.size(); i++){
      TH1F h(singleHists[sampleNames[i]]);
      hComb.Add(&h);
   }
   return hComb;
}

enum FitMode_t
{
   DATA, // normal mode, fitting to data
};

// General implementation of fit. Call on of the functions below, not this one.
FitResult fit(TString sPresel,TString sVar,int iRebin,
              std::vector<float> rebinEdges,std::vector<float> rebinWidths,
              std::vector<std::vector<TString>> fitCombinations,
              std::vector<TString> fixedBkg,
              FitMode_t fitMode,
              int weight_index,
              float f1, float f2 // "closure" scale factors
   )
{
   TH1::SetDefaultSumw2();
   TH1F::SetDefaultSumw2();   
   assert(fitCombinations.size() == 2); // maybe generalize this
   if (fitMode==DATA) assert(f1<0 && f2<0);
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("weight_studies%.1f",cfg.processFraction*100));

   bool isSim=(fitMode!=DATA);
   bool savePlots(f1<0 && f2<0); // don't save scaled closures
   TString saveName=sPresel;
   switch (fitMode) {
   default: break;
   }
   saveName+=+sVar;

   Rebinner rebinned(iRebin);
   if (iRebin==0) rebinned=Rebinner(rebinEdges,rebinWidths);

   // extract the separate backgrounds and the "fit" name
   std::vector<TString> bkgToFit;
   for (std::vector<TString> const &comb: fitCombinations){
      for (unsigned i=1; i<comb.size(); i++){
         bkgToFit.push_back(comb[i]);
      }
   }
   // collect all backgrounds for "raw" plot
   std::vector<TString> bkgAll;
   bkgAll.reserve(fixedBkg.size()+bkgToFit.size());
   bkgAll.insert(bkgAll.end(),fixedBkg.begin(),fixedBkg.end());
   bkgAll.insert(bkgAll.end(),bkgToFit.begin(),bkgToFit.end());

   // extract all bkg hists and scale with trigger efficiency
   std::map<TString,TH1F> bkgHists;
   for (TString sSample: bkgAll){
      bkgHists[sSample]=rebinned(*(TH1F*)histReader.read<TH1F>(sPresel+sVar+"/"+sSample));
      TH1F &h=bkgHists[sSample];
      h.Scale(cfg.trigger_eff_Ph);
      h.SetFillColor(h.GetLineColor());
      h.SetFillStyle(1001);
      h.SetLineColor(kBlack);
      h.SetLineWidth(1);
   }
   bkgHists["efake"].Scale(1./cfg.trigger_eff_Ph); // efake is data-driven

   THStack st;
   for (TString sSample: bkgAll) {
      TH1F *h=(TH1F*)bkgHists[sSample].Clone();
      st.Add(h,"hist");
   }

   // combine samples
   std::map<TString,TH1F> combHists;
   combHists["fixed"]=combineHists(bkgHists,fixedBkg);
   for (std::vector<TString> const& combination: fitCombinations){
      TString const combinationName=combination[0];
      std::vector<TString> toCombine(combination.begin()+1,combination.end());
      combHists[combinationName]=combineHists(bkgHists,toCombine);
      combHists[combinationName].SetFillStyle(0);
      combHists[combinationName].SetLineWidth(2);
   }
   TString sFitSample1=fitCombinations[0][0];  
   TString sFitSample2=fitCombinations[1][0];

   TH1F hData;
   TH1F hSignal=rebinned(*(TH1F*)histReader.read<TH1F>(sPresel+sVar+"/TChiWG"));
   hSignal.Scale(cfg.trigger_eff_Ph);
   TString dataLabel;
   if (fitMode==DATA) {
      // for regular fit, use data histogram
      TH1F* hData_orig=histReader.read<TH1F>(sPresel+sVar+"/SinglePhoton");
      hData=rebinned(*hData_orig);
      dataLabel="data";
   }
   
   TCanvas can;
   hData.SetLineColor(kBlack);
   hData.Draw();
   st.Draw("same");
   hData.Draw("same pe1");
   hSignal.Draw("same hist");
   float plotMax=std::max(st.GetMaximum(),hData.GetMaximum())*1.1;
   st.SetMaximum(plotMax);
   can.RedrawAxis();

   can.Clear();
   gfx::LegendEntries le;
   combHists["(#gamma+)jets"].Draw("hist");
   bkgHists["WGToLNuG"].Draw("same hist");
   bkgHists["ZNuNuGJets"].Draw("same hist");
   bkgHists["WGToLNuG"].SetLineColor(bkgHists["WGToLNuG"].GetFillColor());
   bkgHists["ZNuNuGJets"].SetLineColor(bkgHists["ZNuNuGJets"].GetFillColor());
   bkgHists["WGToLNuG"].SetLineWidth(2);
   bkgHists["ZNuNuGJets"].SetLineWidth(2);
   bkgHists["WGToLNuG"].SetFillStyle(0);
   bkgHists["ZNuNuGJets"].SetFillStyle(0);
   le.append(combHists["(#gamma+)jets"],"(#gamma+)jets","l");
   le.append(bkgHists["WGToLNuG"],"W+#gamma","l");
   le.append(bkgHists["ZNuNuGJets"],"Z+#gamma","l");
   le.buildLegend(.7,.7).DrawClone();

   combHists[sFitSample1].SetLineColor(kOrange-3);
   combHists[sFitSample1].SetFillColor(kOrange-3);
   combHists[sFitSample2].SetLineColor(kAzure-6);
   combHists[sFitSample2].SetFillColor(kAzure-6);
   combHists["fixed"].SetLineColor(kGray+2);
   combHists["fixed"].SetFillColor(kGray+2);

   float n1Before=combHists[sFitSample1].Integral();
   float n2Before=combHists[sFitSample2].Integral();
   float nFixedBefore=combHists["fixed"].Integral();

   hData.Draw("pe1");
   hData.SetMaximum(std::max(st.GetMaximum(),hData.GetMaximum())*1.1);
   combHists[sFitSample1].Draw("same hist");
   combHists[sFitSample2].Draw("same hist");
   hSignal.Draw("same hist");
   can.RedrawAxis();
   le.clear();
   le.append(hData,dataLabel,"pe");
   le.append(combHists[sFitSample1],sFitSample1,"l");
   le.append(combHists[sFitSample2],sFitSample2,"l");
   le.append(hSignal,"TChiWG","l");
   le.buildLegend(.7,.7).DrawClone();


//////////////////////////////////////////////////////////////////////////////////
 ///////////Start changing/////////////////////////////////////////

   //eigener Chi2Test
   gStyle->SetOptStat(0000);

   float Ndata = hData.Integral();
   float total_scale = (Ndata-nFixedBefore)/(n1Before+n2Before);
   float scale_EWK = 0;
   float scale_Multi = 0;
   double fitwidth = 0.15;
   float binerror= 0;

   Double_t chiSquared[250], steps[250];

   TH1F *H_fixed = (TH1F*)combHists["fixed"].Clone("H_fixed");
   for (int i = 0; i <= H_fixed->GetNbinsX(); i++){
      binerror = H_fixed->GetBinError(i)*H_fixed->GetBinError(i);
      binerror = binerror + H_fixed->GetBinContent(i)*0.3 * H_fixed->GetBinContent(i)*0.3;
      binerror = TMath::Sqrt(binerror);
      H_fixed->SetBinError(i,binerror);
   }
      
   for ( int i = 1; i < 250; i++) {
      TH1F *H_EWK = (TH1F*)combHists[sFitSample1].Clone("H_EWK");	
      TH1F *H_Multi = (TH1F*)combHists[sFitSample2].Clone("H_Multi");
  
      steps[i-1] = 0.04 + i*0.02;
        
      scale_EWK = steps[i-1];
      scale_Multi = (total_scale*(n1Before + n2Before) - scale_EWK*n1Before)/n2Before; 

      H_EWK->Scale(scale_EWK);
      H_Multi->Scale(scale_Multi);
      
      H_EWK->Add(H_Multi);
      H_EWK->Add(H_fixed);

      chiSquared[i-1] = hData.Chi2Test(H_EWK,"UW OF CHI2");   
   }
   
   TCanvas can3;
   
	TGraph *chi2_graph = new TGraph(249,steps,chiSquared);
	chi2_graph->GetXaxis()->SetRangeUser(0.4,2.);
	chi2_graph->GetXaxis()->SetTitle("scale factors");
   double y_min = chiSquared[10];    
   for(int i = 1; i<100; i++){
      if(chiSquared[i] < y_min) y_min = chiSquared[i];
   }
   int y_index = std::distance(chiSquared, std::find(chiSquared, chiSquared + 249, y_min));
   double x_of_min = steps[y_index];
	chi2_graph->GetYaxis()->SetRangeUser(0,100); //(0.5,2.7);
	chi2_graph->GetYaxis()->SetTitle("#chi^{2}");
   chi2_graph->GetXaxis()->SetTitle("scale factors");
	chi2_graph->SetMarkerStyle(21);
	chi2_graph->SetMarkerColor(2);
	chi2_graph->SetMarkerSize(.9);
	chi2_graph->GetYaxis()->SetTitleSize(0.09);
   chi2_graph->GetXaxis()->SetTitleSize(0.06);		
	chi2_graph->GetYaxis()->SetLabelSize(0.06);
	chi2_graph->GetXaxis()->SetLabelSize(0.06);
	chi2_graph->GetXaxis()->SetTitleOffset(0.85);
	chi2_graph->GetYaxis()->SetTitleOffset(0.5);		
   chi2_graph->SetTitle("");
	gPad->SetGridx(1);
	gPad->SetGridy(1);	
	gPad->SetBottomMargin(0.12);
   gStyle->SetOptStat(0000);
	chi2_graph->Draw("AP E1 same");
   TLatex *CMS_text = new TLatex(0.12,0.91,"#scale[1]{CMS}#scale[0.3]{ }#scale[0.7]{#bf{#it{private work}}}#bf{                  #scale[0.7]{33.59 fb^{-1} (13 TeV)}}");	
   CMS_text->SetNDC();
   CMS_text->SetTextSize(0.069);
	CMS_text->Draw();		

   gStyle->SetOptFit(0000);
   chi2_graph->Fit("pol2","0 R","",x_of_min - fitwidth,x_of_min + fitwidth);
   
   TF1 *fitFunc = (TF1*)chi2_graph->GetFunction("pol2");
   fitFunc->SetLineColor(kRed);
    
   float a0 = fitFunc->GetParameter(0); //constant term  
   float a1 = fitFunc->GetParameter(1); //linear coefficient, [0] would be constant term
   float a2 = fitFunc->GetParameter(2);  //quadratic coefficient

   float fitscale_Vg = 0; // x value at chi2 miniumum
   float deltaChi2 = 1; // change chi2 by this value for the uncertainty
   float fitErr_Vg = 0;
   if (a2 > 0){
      fitscale_Vg = -a1/(2*a2); 
      fitErr_Vg = sqrt(deltaChi2/a2); // x value of y value increased by one
   }

   saver.save(can3,saveName+"_chi2fit_Vg",isSim);

   for ( int i = 1; i < 250; i++) {
      
      TH1F *H_EWK = (TH1F*)combHists[sFitSample1].Clone("H_EWK");	
      TH1F *H_Multi = (TH1F*)combHists[sFitSample2].Clone("H_Multi");

      steps[i-1] = 0.04 + i*0.02;
        
      scale_Multi = steps[i-1];
      scale_EWK = (total_scale*(n1Before + n2Before) - scale_Multi*n2Before)/n1Before; 

      H_EWK->Scale(scale_EWK);
      H_Multi->Scale(scale_Multi);
      
      H_EWK->Add(H_Multi);
      H_EWK->Add(H_fixed);

      chiSquared[i-1] = hData.Chi2Test(H_EWK,"UW OF CHI2");   
   }
   
   TCanvas can4;
   
	TGraph *chi2_graph_2 = new TGraph(249,steps,chiSquared);
	chi2_graph_2->GetXaxis()->SetRangeUser(0.4,2.);
	chi2_graph_2->GetXaxis()->SetTitle("scale factors");
   y_min = chiSquared[10];    
   for(int i = 1; i<100; i++){
      if(chiSquared[i] < y_min) y_min = chiSquared[i];
   }
   y_index = std::distance(chiSquared, std::find(chiSquared, chiSquared + 249, y_min));
   x_of_min = steps[y_index];
	chi2_graph_2->GetYaxis()->SetRangeUser(0,100); //(0.5,2.7);
	chi2_graph_2->GetYaxis()->SetTitle("#chi^{2}");
   chi2_graph_2->GetXaxis()->SetTitle("scale factors");
	chi2_graph_2->SetMarkerStyle(21);
	chi2_graph_2->SetMarkerColor(2);
	chi2_graph_2->SetMarkerSize(.9);
	chi2_graph_2->GetYaxis()->SetTitleSize(0.09);
   chi2_graph_2->GetXaxis()->SetTitleSize(0.06);		
	chi2_graph_2->GetYaxis()->SetLabelSize(0.06);
	chi2_graph_2->GetXaxis()->SetLabelSize(0.06);
	chi2_graph_2->GetXaxis()->SetTitleOffset(0.85);
	chi2_graph_2->GetYaxis()->SetTitleOffset(0.5);		
   chi2_graph_2->SetTitle("");
	gPad->SetGridx(1);
	gPad->SetGridy(1);	
	gPad->SetBottomMargin(0.12);
	chi2_graph_2->Draw("AP E1 same");
	CMS_text->Draw();		

   gStyle->SetOptFit(0000); 
   chi2_graph_2->Fit("pol2","0R","",x_of_min - fitwidth,x_of_min + fitwidth);

   TF1 *fitFunc_2 = (TF1*)chi2_graph_2->GetFunction("pol2");
   fitFunc_2->SetLineColor(kRed);
    
   a0 = fitFunc_2->GetParameter(0); //constant term  
   a1 = fitFunc_2->GetParameter(1); //linear coefficient, [0] would be constant term
   a2 = fitFunc_2->GetParameter(2);  //quadratic coefficient

   float fitscale_gJ = 0; // x value at chi2 miniumum
   deltaChi2 = 1; // change chi2 by this value for the uncertainty
   float fitErr_gJ = 0;
 
   if (a2 > 0){
      fitscale_gJ = -a1/(2*a2);
      fitErr_gJ = sqrt(deltaChi2/a2); // x value of y value increased by one      
   }

   gStyle->SetOptStat(0000);
   saver.save(can4,saveName+"_chi2fit_gJ",isSim);   
   //From here only best fit values
   
   float sf1=fitscale_Vg;
   float sf2=fitscale_gJ;
   float e_sf1=fitErr_Vg;
   float e_sf2=fitErr_gJ;

   float n1after= n1Before*sf1;
   float n2after= n2Before*sf2;   
   
   TString sResult1=TString::Format("%.2f#pm%.2f",sf1,e_sf1);
   TString sResult2=TString::Format("%.2f#pm%.2f",sf2 ,e_sf2);
   io::log<<sFitSample1+TString::Format(": %.0f -> %.0f ",n1Before,n1after)+sResult1;
   io::log<<sFitSample2+TString::Format(": %.0f -> %.0f ",n2Before,n2after)+sResult2;

   combHists[sFitSample1].Scale(sf1);
   combHists[sFitSample2].Scale(sf2);

   gfx::SplitCan spcan(.7);
   spcan.cdUp();
   //add up all histos for fit result 
   TH1F hFitTotal(combHists[sFitSample1]);
   hFitTotal.Add(&combHists[sFitSample2]);
   hFitTotal.Add(&combHists["fixed"]);
 //  hFitTotal.Add(H_fixed);
   hFitTotal.SetFillStyle(1001);
   hFitTotal.SetMarkerSize(0);
   hFitTotal.SetFillColor(kWhite);
   hFitTotal.SetLineColor(kRed+1);
   hFitTotal.Draw("hist");
   hFitTotal.SetMaximum(std::max(st.GetMaximum(),hData.GetMaximum())*1.1);
   TH1F hFitTotalErr(hFitTotal);
   TH1F white(hFitTotal);
     
   hFitTotalErr.SetLineColor(kWhite);
   hFitTotalErr.SetFillColor(kRed+1);
   hFitTotalErr.SetFillStyle(3354);
   hFitTotalErr.Draw("same e2");

   TH1F hFixedErr(combHists["fixed"]);
   for (int i = 0; i <= H_fixed->GetNbinsX(); i++){
      hFixedErr.SetBinContent(i,hFitTotal.GetBinContent(i));      
      hFixedErr.SetBinError(i,H_fixed->GetBinError(i));
   }
   hFixedErr.SetFillStyle(3345);
   hFixedErr.SetFillColor(kMagenta+4);
   hFixedErr.SetLineColor(kWhite);
   hFixedErr.SetMarkerSize(0);
   hFixedErr.Draw("same e2");
   
   combHists["fixed"].SetFillStyle(1001);
   combHists["fixed"].SetFillColor(kMagenta+4);
   combHists["fixed"].Draw("same hist");
   combHists[sFitSample1].SetLineColor(kOrange-3);
   combHists[sFitSample2].SetLineColor(kAzure+10);  
   combHists[sFitSample1].Draw("same hist");
   combHists[sFitSample2].Draw("same hist");
   hData.Draw("same pe1");

   white.SetFillColor(kWhite);
   white.SetLineColor(kWhite);
   white.SetFillStyle(1001);

   le.clear();
   le.append(hData,dataLabel,"pe");
   le.append(combHists[sFitSample1],sFitSample1+": "+sResult1,"l");
   le.append(white,"","f");
   le.append(combHists[sFitSample2],sFitSample2+": "+sResult2,"l");
   le.append(combHists["fixed"],"fixed bkg","f");
   le.append(hFitTotal,"total fit","l");
   le.append(hFixedErr,"#sigma_{syst, fixed}","f");   
   le.append(hFitTotalErr,"#sigma_{stat}","f");
   TLegend leg=le.buildLegend(.36,.69,1-gPad->GetRightMargin(),-1,2);
   leg.Draw();
   spcan.pU_.RedrawAxis();

   spcan.cdLow();
   TH1F hPull=hist::getPull(hData,hFitTotal);
   hPull.SetFillColor(kRed+1);  
   hPull.Draw("hist");
   if (savePlots) saver.save(spcan,saveName+"_fit",isSim);

   FitResult fr = {
      sPresel,sVar,hData.GetXaxis()->GetTitle(),
      rebinned.numberOfBins(),
      sf1,sf2,e_sf1,e_sf2,
      weight_index,
      0 //correlation
   };
   g_fitResults.push_back(fr);
   return fr;
}

// fit functions to call
FitResult fit(TString sPresel,TString sVar,int iRebin,
              std::vector<std::vector<TString>> fitCombinations,
              std::vector<TString> fixedBkg,
              FitMode_t fitMode=DATA,
              int weight_index = -1,
              float f1=-1, float f2=-1 // "closure" scale factors
   )
{
   return fit(sPresel,sVar,iRebin,{},{},fitCombinations,fixedBkg,fitMode,weight_index,f1,f2);
}
FitResult fit(TString sPresel,TString sVar,std::vector<float> rebinEdges,std::vector<float> rebinWidths,
              std::vector<std::vector<TString>> fitCombinations,
              std::vector<TString> fixedBkg,
              FitMode_t fitMode=DATA,
              int weight_index =-1,
              float f1=-1, float f2=-1 // "closure" scale factors
   )
{
   return fit(sPresel,sVar,0,rebinEdges,rebinWidths,fitCombinations,fixedBkg,fitMode,weight_index,f1,f2);
}

std::string getSelectionString(FitResult res){
   std::string selection(res.selection);
   std::vector<std::string> parts=util::to_vector<std::string>(selection,'/');
   selection.clear();
   selection=parts[(parts.size()-1)];
   return selection;
}

enum SummaryMode_t
{
   SRT_DST, // sort by distribution
   SRT_SEL, // sort by selection
};

void resultSummary(FitResult referenceResult,SummaryMode_t summaryMode,FitMode_t fitMode)
{
   io::log<<"";
   io::log<<"##### SUMMARY #####";

   bool isSim=(fitMode!=DATA);

   float coord1=-1.;
   float coord2=2.0;
   float coordRange=coord2-coord1;
   float axRange1=0.39999;
   float axRange2=2.1001;

   int iResults=g_fitResults.size();
   debug << "hier" << iResults;
   int iResult=0;

   TCanvas can;
   can.SetLeftMargin(0.0);
   // draw histogram without axes -> specifies coordianate system
   TH1F hFrame("hFrame","",10,coord1,coord2);
   hFrame.SetLineColor(kWhite);
   hFrame.Draw("ah");
   can.SetFrameLineColor(kWhite); // no box around plot

   TGraphErrors gr1(iResults);
   TGraphErrors gr2(iResults);
   TH1F hPull1("",";pull;EntriesBIN",10,-4,4);
   TH1F hPull2("",";pull;EntriesBIN",10,-4,4);

   TLatex ltx;
   ltx.SetTextFont(42);
   ltx.SetTextAlign(21);
   ltx.SetTextSize(.04);
   TLine line;
   TBox box;
   box.SetFillStyle(1001);

   TString fitModeLabel;

   float reference1=1;
   float reference2=1;
   if (fitMode==DATA){
      reference1=referenceResult.r1;
      reference2=referenceResult.r2;

      //include scale uncertainty
      float scale_uncert_Vg = 0.;
      float scale_uncert_gJ = 0.;
      float errorVg = TMath::Sqrt(referenceResult.e1*referenceResult.e1 + scale_uncert_Vg*scale_uncert_Vg);
      float errorgJets = TMath::Sqrt(referenceResult.e2*referenceResult.e2 + scale_uncert_gJ*scale_uncert_gJ);     

   //   box.SetFillColor(fillcolor1);
  //    box.DrawBox(referenceResult.r1-referenceResult.e1,0,referenceResult.r1+referenceResult.e1,1);
  //    box.DrawBox(referenceResult.r1-errorVg,0,referenceResult.r1+errorVg,1);
      box.SetFillColor(fillcolor1_2);
      box.DrawBox(referenceResult.r1-referenceResult.e1,0,referenceResult.r1+referenceResult.e1,1);     
 //     box.SetFillColor(fillcolor2);
  //    box.DrawBox(referenceResult.r2-referenceResult.e2,0,referenceResult.r2+referenceResult.e2,1);
 //     box.DrawBox(referenceResult.r2-errorgJets,0,referenceResult.r2+errorgJets,1);
      box.SetFillColor(fillcolor2_2);
      box.DrawBox(referenceResult.r2-referenceResult.e2,0,referenceResult.r2+referenceResult.e2,1);
      
      line.SetLineColor(color1);
      line.DrawLine(referenceResult.r1,0,referenceResult.r1,1);
      line.SetLineColor(color2);
      line.DrawLine(referenceResult.r2,0,referenceResult.r2,1);
      ltx.SetTextColor(color1);
      ltx.DrawLatex(referenceResult.r1,1.02,"V(+#gamma)");
      ltx.SetTextColor(color2);
      ltx.DrawLatex(referenceResult.r2,1.02,"(#gamma+)jets");
   }
   
   line.SetLineStyle(2);
   line.SetLineColor(kBlack);
   ltx.SetTextSize(TMath::Min(.6/iResults,0.025));
   ltx.SetTextAlign(12);
   ltx.SetTextColor(kBlack);

   std::vector<FitResult> fitResults(g_fitResults);
   if (summaryMode==SRT_SEL)      std::sort(fitResults.begin(),fitResults.end(),&cmpSelection);
   else if (summaryMode==SRT_DST) std::sort(fitResults.begin(),fitResults.end(),&cmpDistribution);

   TString currentGroup;

   float min1=referenceResult.r1-referenceResult.e1;
   float max1=referenceResult.r1+referenceResult.e1;
   float min2=referenceResult.r2-referenceResult.e2;
   float max2=referenceResult.r2+referenceResult.e2;
   io::Logger fit_result_weights(("weights_fitresult.tex"));
   fit_result_weights<< "# weight_index, weight, Vgamma scale, gJets scale";
   TString fill_line;

   for (FitResult const &fr: fitResults){
      std::string selection="weight "+std::to_string(fr.weight_index);

      if (fr.distribution.Contains("JESu")) selection+=" (JES up)";
      if (fr.distribution.Contains("JESd")) selection+=" (JES down)";

      // print results
      io::log>>TString::Format("%.2f+-%.2f %.2f+-%.2f (rho=%.2f)",fr.r1,fr.e1,fr.r2,fr.e2,fr.rho);
      // vizualize results
      float const y=1-float(iResult+.5)/iResults;
      gr1.SetPoint(iResult,fr.r1,y);
      gr1.SetPointError(iResult,fr.e1,0);
      gr2.SetPoint(iResult,fr.r2,y);
      gr2.SetPointError(iResult,fr.e2,0);

      //systematic studies to check the shape variation
      fill_line = std::to_string(fr.weight_index) +" " + std::to_string(fr.r1)+ " "+ std::to_string(fr.r2) ;
      fit_result_weights<< fill_line;
   
      TString label=selection;
      ltx.DrawLatex(coord1+.01*coordRange,y,label);

      TString const &newGroup=(summaryMode==SRT_SEL ? selection : fr.distribution);
      if (currentGroup=="") currentGroup=newGroup;
      if (newGroup!=currentGroup){
         float const yy=1-float(iResult)/iResults;
         line.DrawLine(coord1+.01*coordRange,yy,coord2,yy);
         currentGroup=newGroup;
      }
      min1=TMath::Min(min1,fr.r1);
      max1=TMath::Max(max1,fr.r1);
      min2=TMath::Min(min2,fr.r2);
      max2=TMath::Max(max2,fr.r2);

      hPull1.Fill((fr.r1-reference1)/fr.e1);
      hPull2.Fill((fr.r2-reference2)/fr.e2);

      iResult++;
   }

   gr1.SetMarkerColor(color1);
   gr1.SetLineColor(color1);
   gr1.SetLineWidth(2);
   gr1.SetMarkerStyle(kOpenCircle);
   gr2.SetMarkerColor(color2);
   gr2.SetLineColor(color2);
   gr2.SetLineWidth(2);
   gr2.SetMarkerStyle(kOpenSquare);
   gr1.Draw("p");
   gr2.Draw("p");
   TGaxis axis(axRange1,0,axRange2,0,axRange1,axRange2,505,"");
   axis.SetTitle("SF fit");
   axis.CenterTitle(true);
   axis.SetTitleFont(42);
   axis.SetLabelFont(42);
   axis.SetTitleSize(.06);
   axis.SetLabelSize(.05);
   axis.Draw();

   TString saveName;
   switch (fitMode) {
   case DATA:    { saveName="summary/"; break; }
   }
   switch (summaryMode) {
   case SRT_DST: { saveName+="dist"; break; }
   case SRT_SEL: { saveName+="sel"; break; }
   }
   saver.save(can,saveName+"_1d",true,isSim);

   // Pulls
   TCanvas can2;
   hPull1.SetLineColor(color1);
   hPull2.SetLineColor(color2);
   hPull1.Draw("hist");
   hPull2.Draw("hist same");
   hist::setMaximum(hPull1,{hPull2});

   gfx::cornerLabel(fitModeLabel,1).DrawClone();
   TString text=TString::Format("#mu=%.1f#pm%.1f #sigma=%.1f",
                                hPull1.GetMean(),hPull1.GetMeanError(),hPull1.GetStdDev());
   TLatex label=gfx::cornerLabel(text);
   label.SetTextColor(color1);
   label.DrawClone();
   text=TString::Format("#mu=%.1f#pm%.1f #sigma=%.1f",
                        hPull2.GetMean(),hPull2.GetMeanError(),hPull2.GetStdDev());
   label.SetTitle(text);
   label.SetY(label.GetY()-0.05);
   label.SetTextColor(color2);
   label.DrawClone();

   saver.save(can2,saveName+"_pull");
}

void run (FitMode_t fitMode){
   TH1::SetDefaultSumw2();
   g_fitResults.clear();
   std::vector<TString> fixedBkg={"diboson","TTJets","TTGJets","efake"};
   std::vector<std::vector<TString>> fitCombinations={{"V(+#gamma)","WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"},{"#gamma+jets","GJets_DR"}};

   FitResult referenceResult;
   int nWeights= 109;

   referenceResult=fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_0",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode,0);
   for (int i = 1; i < nWeights; i++){ //i = 0 is reference fit
      fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_"+std::to_string(i),{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode,i);
   }
   resultSummary(referenceResult,SRT_SEL,fitMode);

}

extern "C"
void run()
{
   run(DATA);
}

