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
io::RootFileSaver saver("plots.root","roofit");
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
   float rho; // correlation
};
bool cmpSelection   (const FitResult &a,const FitResult &b) { return a.selection.CompareTo(b.selection)      >0; }
bool cmpDistribution(const FitResult &a,const FitResult &b) { return a.distribution.CompareTo(b.distribution)>0; }
bool cmpNoBins      (const FitResult &a,const FitResult &b) { return a.numberOfBins>b.numberOfBins; }

std::vector<FitResult> g_fitResults;

RooDataHist toRooDataHist(RooRealVar const& x, TH1 const &h)
{
   return RooDataHist(h.GetName(),h.GetTitle(),x,&h);
}

TH1F combineHists(std::map<TString,TH1F> singleHists,std::vector<TString> sampleNames)
{
   assert(sampleNames.size()>0);
   TH1F hComb(singleHists[sampleNames[0]]);
   hComb.Reset();
   hComb.SetLineColor(singleHists[sampleNames[0]].GetLineColor());
   hComb.SetFillColor(singleHists[sampleNames[0]].GetLineColor());
   for (unsigned i=0; i<sampleNames.size(); i++){
      TH1F h(singleHists[sampleNames[i]]);
      // to change the W/Z ratio and check the fit result
      // if (sampleNames[i].Contains("W")) {
      //    io::log*sampleNames[i];
      //    h.Scale(50/78.);
      //    io::log>>h.Integral();
      // }
      // if (sampleNames[i].Contains("Z")) {
      //    io::log*sampleNames[i];
      //    h.Scale(50/22.);
      //    io::log>>h.Integral();
      // }
      hComb.Add(&h);
   }
   // debug<<hComb.Integral();
   return hComb;
}

enum FitMode_t
{
   DATA, // normal mode, fitting to data
   CLOSURE, // do a closure by generating "per bin" toy data
   CONT_GGM, // same as closure, but with signal contamination
};

// General implementation of fit. Call on of the functions below, not this one.
FitResult fit(TString sPresel,TString sVar,int iRebin,
              std::vector<float> rebinEdges,std::vector<float> rebinWidths,
              std::vector<std::vector<TString>> fitCombinations,
              std::vector<TString> fixedBkg,
              FitMode_t fitMode,
              float f1, float f2 // "closure" scale factors
   )
{
   TH1::SetDefaultSumw2();
   TH1F::SetDefaultSumw2();   
   assert(fitCombinations.size() == 2); // maybe generalize this
   if (fitMode==DATA) assert(f1<0 && f2<0);
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   bool isSim=(fitMode!=DATA);
   bool savePlots(f1<0 && f2<0); // don't save scaled closures
   TString saveName=sPresel;
   switch (fitMode) {
   case CLOSURE: { saveName+="closure/";break; }
   case CONT_GGM:{ saveName+="contamin/";break; }
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
      dataLabel="Data";
   } else {
      // add up all MC
      hData=combHists["fixed"];
      TH1F hToAdd1(combHists[sFitSample1]);
      TH1F hToAdd2(combHists[sFitSample2]);
      if (f1>0) hToAdd1.Scale(f1);
      if (f2>0) hToAdd2.Scale(f2);
      hData.Add(&hToAdd1);
      hData.Add(&hToAdd2);

      if (fitMode==CLOSURE){
         dataLabel="toy data";
      } else if (fitMode==CONT_GGM){
         dataLabel="signal inj. toy";
         // inject signal
         hData.Add(&hSignal);
      }

      // dice from a poissonian for each bin
      for (int i=0; i<=hData.GetNbinsX()+1;i++){
         int const v=g_rand.Poisson(hData.GetBinContent(i));
         hData.SetBinContent(i,v);
         hData.SetBinError(i,TMath::Sqrt(v));
      }
      hData.SetMarkerStyle(kOpenCircle);
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
   if (savePlots) saver.save(can,saveName+"_rawstack",true,isSim);

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
   if (savePlots) saver.save(can,saveName+"_raw",true,isSim);

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

   if (savePlots) saver.save(can,saveName+"_combined",true,isSim);

   // RooMsgService::instance().setSilentMode(true); // no effect for current usage
   // https://root.cern.ch/doc/master/namespaceRooFit.html#a36e12ae9ea9c0d3ab48f71c4ffcdace3
   // RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

   /*
   float rangeL=hData.GetXaxis()->GetBinLowEdge(1);
   float rangeU=hData.GetXaxis()->GetBinUpEdge(hData.GetNbinsX());
   RooRealVar x("x","x",rangeL,rangeU);
  // RooDataHist dhData=toRooDataHist(x,hData);
   RooDataHist dhData=toRooDataHist(x,hData);  
   RooDataHist dhSample1=toRooDataHist(x,combHists[sFitSample1]);
   RooDataHist dhSample2=toRooDataHist(x,combHists[sFitSample2]);
   RooDataHist dhFixed  =toRooDataHist(x,combHists["fixed"]);

   RooRealVar nSample1("nSample1","sample1 count",n1Before,2,2*hData.Integral());
   RooRealVar nSample2("nSample2","sample2 count",n2Before,2,2*hData.Integral());
   RooRealVar nFixed  ("nFixed"  ,"fixed count"  ,nFixedBefore,nFixedBefore,nFixedBefore);
   nFixed.setConstant(true);

   RooHistPdf pdfSample1("pdfSample1","pdfSample1",x,dhSample1);
   RooHistPdf pdfSample2("pdfSample2","pdfSample2",x,dhSample2);
   RooHistPdf pdfFixed  ("pdfFixed","pdfFixed",x,dhFixed);


   TH1* check = pdfSample1.createHistogram("pdf",x);
   TH1* check2 = pdfSample2.createHistogram("pdf",x);
   
   double eintrag = 0;
   double error = 0;   
   std::cout << "****************************" << std::endl;
   for (int i = 1; i < check->GetNbinsX(); i++){
      eintrag = check->GetBinContent(i);
      std::cout << "binContent:   " << eintrag << std::endl;
      error = check->GetBinError(i);
      std::cout << "fehler:    " << error << std::endl;     

   }
   for (int i = 1; i < check2->GetNbinsX(); i++){
      eintrag = check2->GetBinContent(i);
      std::cout << "2 binContent:   " << eintrag << std::endl;
      error = check2->GetBinError(i);
      std::cout << "2 fehler:    " << error << std::endl;     

   }   
   std::cout << "****************************" << std::endl;
   */
   
 /* 
   RooExtendPdf pdfExtSample1("pdfSample1","pdfSample1",pdfSample1,nSample1);
   RooExtendPdf pdfExtSample2("pdfSample2","pdfSample2",pdfSample2,nSample2);
   RooExtendPdf pdfExtFixed  ("pdfFixed","pdfFixed",pdfFixed,nFixed);
   
   RooAddPdf model("model","model",RooArgList(pdfSample1,pdfSample2,pdfFixed),RooArgList(nSample1,nSample2,nFixed));
 //  RooAddPdf model("model","model",RooArgList(pdfExtSample1,pdfExtSample2,pdfExtFixed));
   RooChi2Var chi2("chi2","chi2",model,dhData,RooAbsData::SumW2);
   RooMinuit m(chi2);
   m.migrad();
   m.hesse();
   * */
//   RooFitResult *fitResult=model.chi2FitTo(dhData,RooFit::Extended(true),RooFit::Save(),RooFit::SumW2Error(true));
//   float correlation=fitResult->correlation(nSample1,nSample2);
   float correlation = 0.5;

   
//////////////////////////////////////////////////////////////////////////////////
 ///////////Start changing/////////////////////////////////////////
 /*
   float n1after=nSample1.getValV();
   float n2after=nSample2.getValV();
   float e_n1after=nSample1.getError();
   float e_n2after=nSample2.getError(); */

   
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
 //     debug << "Fehler:  " << H_fixed->GetBinError(i);
      binerror = H_fixed->GetBinError(i)*H_fixed->GetBinError(i);
      binerror = binerror + H_fixed->GetBinContent(i)*0.3 * H_fixed->GetBinContent(i)*0.3;
      binerror = TMath::Sqrt(binerror);
  //    debug << "Fehler:  " << binerror;
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
      if(chiSquared[i] < y_min)
      y_min = chiSquared[i];
     }
   int y_index = std::distance(chiSquared, std::find(chiSquared, chiSquared + 249, y_min));
    std::cout << "y min:  " << y_min << std::endl;  
   std::cout << "y index:  " << y_index << std::endl;
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

	std::cout << "  -------------------------------------------" << std::endl;
   std::cout << " Fit" << std::endl;
   std::cout << "  -------------------------------------------" << std::endl;
   gStyle->SetOptFit(0000);
   chi2_graph->Fit("pol2","0 R","",x_of_min - fitwidth,x_of_min + fitwidth);
   
   TF1 *fitFunc = (TF1*)chi2_graph->GetFunction("pol2");
   fitFunc->SetLineColor(kRed);
    
   float a0 = fitFunc->GetParameter(0); //constant term  
   float a1 = fitFunc->GetParameter(1); //linear coefficient, [0] would be constant term
   float a2 = fitFunc->GetParameter(2);  //quadratic coefficient
   double offset_y_m = 0;

   float fitscale_Vg = 0; // x value at chi2 miniumum
   float deltaChi2 = 1; // change chi2 by this value for the uncertainty
   float fitErr_Vg = 0;
   if (a2 > 0){
      fitscale_Vg = -a1/(2*a2); 
      fitErr_Vg = sqrt(deltaChi2/a2); // x value of y value increased by one
      offset_y_m = a0 - 0.25*a1*a1/a2;
   }
	double offset_x_m = fitscale_Vg;
	double mean_up_m = fitscale_Vg + fitErr_Vg;
	double mean_down_m = fitscale_Vg - fitErr_Vg;

	std::cout << std::endl;
	std::cout << "  ***Fit to scale" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "  RESULTS:" << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  offset_y = " << offset_y_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  offset_x = " << offset_x_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  sigma_x_up = " << mean_up_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  sigma_x_down = " << mean_down_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  error = " << fitErr_Vg << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;

   saver.save(can3,saveName+"_chi2fit_Vg",isSim);

   for ( int i = 1; i < 250; i++) {
      
      TH1F *H_EWK = (TH1F*)combHists[sFitSample1].Clone("H_EWK");	
      TH1F *H_Multi = (TH1F*)combHists[sFitSample2].Clone("H_Multi");
  //    TH1F *H_fixed = (TH1F*)combHists["fixed"].Clone("H_fixed");
         
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
      if(chiSquared[i] < y_min)
      y_min = chiSquared[i];
     }
   y_index = std::distance(chiSquared, std::find(chiSquared, chiSquared + 249, y_min));
   std::cout << "y min:  " << y_min << std::endl;  
   std::cout << "y index:  " << y_index << std::endl;
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

	std::cout << "  -------------------------------------------" << std::endl;
   std::cout << " Fit" << std::endl;
   std::cout << "  -------------------------------------------" << std::endl;
   gStyle->SetOptFit(0000); 
   chi2_graph_2->Fit("pol2","0R","",x_of_min - fitwidth,x_of_min + fitwidth);

   
   TF1 *fitFunc_2 = (TF1*)chi2_graph_2->GetFunction("pol2");
   fitFunc_2->SetLineColor(kRed);

     std::cout << "  -------------------HERE------------------------" << std::endl;
    
   a0 = fitFunc_2->GetParameter(0); //constant term  
   a1 = fitFunc_2->GetParameter(1); //linear coefficient, [0] would be constant term
   a2 = fitFunc_2->GetParameter(2);  //quadratic coefficient

   float fitscale_gJ = 0; // x value at chi2 miniumum
   deltaChi2 = 1; // change chi2 by this value for the uncertainty
   float fitErr_gJ = 0;


     
   if (a2 > 0){
      fitscale_gJ = -a1/(2*a2);
      fitErr_gJ = sqrt(deltaChi2/a2); // x value of y value increased by one      
      offset_y_m = a0 - 0.25*a1*a1/a2;
   }
	offset_x_m = fitscale_gJ;
	mean_up_m = fitscale_gJ + fitErr_gJ;
	mean_down_m = fitscale_gJ - fitErr_gJ;

	std::cout << std::endl;
	std::cout << "  ***Fit to scale" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "  RESULTS:" << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  offset_y = " << offset_y_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  offset_x = " << offset_x_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  sigma_x_up = " << mean_up_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  sigma_x_down = " << mean_down_m << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
	std::cout << "  error = " << fitErr_gJ << std::endl;
	std::cout << "  -------------------------------------------" << std::endl;
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
   hFitTotal.GetYaxis()->SetTitleOffset(hFitTotal.GetYaxis()->GetTitleOffset()*1.1);
   hFitTotal.Draw("hist");
   hFitTotal.SetMaximum(std::max(st.GetMaximum()*1.3,hData.GetMaximum())*1.3);
   hFitTotal.SetMinimum(0.1);   
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
   le.append(combHists[sFitSample1],"SF_{V(#gamma)}:  _{ }"+sResult1,"l");
   le.append(white,"","f");
   le.append(combHists[sFitSample2],"SF_{#gamma+jets}: "+sResult2,"l");
   le.append(combHists["fixed"],"Fixed bkg","f");
   le.append(hFixedErr,"#sigma_{syst, fixed}","f");  
   le.append(hFitTotal,"Total fit","l"); 
   le.append(hFitTotalErr,"#sigma_{stat}","f");
   TLegend leg=le.buildLegend(.32,.69,0.97-gPad->GetRightMargin(),-1,2);
   leg.Draw();
   spcan.pU_.RedrawAxis();

   spcan.cdLow();
   double error;
   for (int i = 0; i <= H_fixed->GetNbinsX(); i++){
      error = hFitTotal.GetBinError(i)*hFitTotal.GetBinError(i) + H_fixed->GetBinError(i)*H_fixed->GetBinError(i);
      hFitTotal.SetBinError(i,TMath::Sqrt(error));
   }
   TH1F hPull=hist::getPull(hData,hFitTotal);
   hPull.SetFillColor(kRed+1);
   hPull.GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma_{stat}}");
   hPull.GetYaxis()->SetTitleOffset(hFitTotal.GetYaxis()->GetTitleOffset()*0.1);
   hPull.GetYaxis()->SetTitleSize(0.17);
   hPull.Draw("hist");
   hPull.GetYaxis()->SetTitleOffset(hPull.GetYaxis()->GetTitleOffset()*0.1);
   hPull.GetYaxis()->SetTitleSize(0.17);
   hPull.Draw("hist same"); 
   if (savePlots) saver.save(spcan,saveName+"_fit",isSim);

   can.cd();
   // CL contour drawing:
   // - RooFitResult::plotOn
   // - more general: RooMinuit::contour

   // automatic drawing of CL ellipse with error+covariance lines
   // RooPlot* frame = new RooPlot(nSample1,nSample2,50,150,230,350);
   // fitResult->plotOn(frame,nSample1,nSample2,"ME12AHV");
   // frame->Draw();
   // RooEllipse* err_ell=(RooEllipse*)frame->findObject("contour");
   // manually:
   // RooEllipse ell("error-contour",n1after,n2after,e_n1after,e_n2after,correlation,1000);
   // ell.Draw("same");

   // translate to scale factors
   RooEllipse ell(TString::Format(";SF_{%s};SF_{%s}",sFitSample1.Data(),sFitSample2.Data()),sf1,sf2,e_sf1,e_sf2,correlation,1000);
   ell.SetFillColor(kGray);
   ell.SetLineColor(kGray);
   ell.Draw("af");
   TMarker bestFit(sf1,sf2,kFullStar);
   bestFit.SetMarkerSize(4);
   bestFit.Draw();
   le.clear();
   le.append(bestFit,"best fit","p");
   le.append(TString::Format("(%.2f,%.2f)",sf1,sf2));
   le.append(TString::Format("#rho=%.2f",correlation));
   le.append(ell,"68% CL","f");
   le.buildLegend(.7,.75).DrawClone();

   if (savePlots) saver.save(can,saveName+"_cont",true,isSim);

   // store the result
   FitResult fr = {
      sPresel,sVar,hData.GetXaxis()->GetTitle(),
      rebinned.numberOfBins(),
      sf1,sf2,e_sf1,e_sf2,
      correlation
   };
   g_fitResults.push_back(fr);
   return fr;
}

// fit functions to call
FitResult fit(TString sPresel,TString sVar,int iRebin,
              std::vector<std::vector<TString>> fitCombinations,
              std::vector<TString> fixedBkg,
              FitMode_t fitMode=DATA,
              float f1=-1, float f2=-1 // "closure" scale factors
   )
{
   return fit(sPresel,sVar,iRebin,{},{},fitCombinations,fixedBkg,fitMode,f1,f2);
}
FitResult fit(TString sPresel,TString sVar,std::vector<float> rebinEdges,std::vector<float> rebinWidths,
              std::vector<std::vector<TString>> fitCombinations,
              std::vector<TString> fixedBkg,
              FitMode_t fitMode=DATA,
              float f1=-1, float f2=-1 // "closure" scale factors
   )
{
   return fit(sPresel,sVar,0,rebinEdges,rebinWidths,fitCombinations,fixedBkg,fitMode,f1,f2);
}

std::string getSelectionString(FitResult res){
   std::string selection(res.selection);
   std::vector<std::string> parts=util::to_vector<std::string>(selection,'/');
   selection.clear();
   for (unsigned i=1; i<parts.size(); i++){
      selection+=parts[i]+" ";
   }
   selection=selection.substr(2,selection.size()-3); // remove "c_" and whitespace
   return selection;
}

enum SummaryMode_t
{
   SRT_DST, // sort by distribution
   SRT_SEL, // sort by selection
   BINNING, // show number of bins
};

void resultSummary(FitResult referenceResult,SummaryMode_t summaryMode,FitMode_t fitMode)
{
   io::log<<"";
   io::log<<"##### SUMMARY #####";

   bool isSim=(fitMode!=DATA);

   float coord1=-1.;
   float coord2=2.2;
   float coordRange=coord2-coord1;
   float axRange1=0.39999;
   float axRange2=2.2001;

   int iResults=g_fitResults.size();
   int iResult=0;

   TCanvas can;
   can.SetLeftMargin(0.0);
   // draw histogram without axes -> specifies coordianate system
   TH1F hFrame("hFrame","",10,coord1,coord2);
   hFrame.SetLineColor(kWhite);
   // hFrame.Draw("axis"); // for debugging (check axis range)
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
      float errorVg = TMath::Sqrt(referenceResult.e1*referenceResult.e1 + 0.0722*0.0722);
      float errorgJets = TMath::Sqrt(referenceResult.e2*referenceResult.e2 + 0.097*0.097);     

      box.SetFillColor(fillcolor1);
  //    box.DrawBox(referenceResult.r1-referenceResult.e1,0,referenceResult.r1+referenceResult.e1,1);
      box.DrawBox(referenceResult.r1-errorVg,0,referenceResult.r1+errorVg,1);
      box.SetFillColor(fillcolor1_2);
      box.DrawBox(referenceResult.r1-referenceResult.e1,0,referenceResult.r1+referenceResult.e1,1);     
      box.SetFillColor(fillcolor2);
  //    box.DrawBox(referenceResult.r2-referenceResult.e2,0,referenceResult.r2+referenceResult.e2,1);
      box.DrawBox(referenceResult.r2-errorgJets,0,referenceResult.r2+errorgJets,1);
      box.SetFillColor(fillcolor2_2);
      box.DrawBox(referenceResult.r2-referenceResult.e2,0,referenceResult.r2+referenceResult.e2,1);
      
      line.SetLineColor(color1);
      line.DrawLine(referenceResult.r1,0,referenceResult.r1,1);
      line.SetLineColor(color2);
      line.DrawLine(referenceResult.r2,0,referenceResult.r2,1);
      ltx.SetTextColor(color1);
      ltx.DrawLatex(referenceResult.r1,1.02,"SF_{V(#gamma)}");
      ltx.SetTextColor(color2);
      ltx.DrawLatex(referenceResult.r2,1.02,"SF_{(#gamma+)jets}");
   } else {
      if (fitMode==CLOSURE){
         fitModeLabel="Closure test";
      } else if (fitMode==CONT_GGM){
         fitModeLabel="Signal contamination test";
      }
      line.SetLineColor(kBlack);
      line.DrawLine(1,0,1,1);

      ltx.SetTextColor(color1);
      ltx.DrawLatex(.5*(1+axRange1),1.02,"SF_{V(#gamma)}");
      ltx.SetTextColor(color2);
      ltx.DrawLatex(.5*(1+axRange2),1.02,"SF_{(#gamma+)jets}");

      ltx.SetTextColor(kBlack);
      ltx.SetTextAlign(11);
      ltx.DrawLatex(coord1+.01*coordRange,1.02,fitModeLabel);
   }

   line.SetLineStyle(2);
   line.SetLineColor(kBlack);
   ltx.SetTextSize(TMath::Min(.6/iResults,0.025));
   ltx.SetTextAlign(12);
   ltx.SetTextColor(kBlack);

   std::vector<FitResult> fitResults(g_fitResults);
   if (summaryMode==SRT_SEL)      std::sort(fitResults.begin(),fitResults.end(),&cmpSelection);
   else if (summaryMode==SRT_DST) std::sort(fitResults.begin(),fitResults.end(),&cmpDistribution);
   else if (summaryMode==BINNING) std::sort(fitResults.begin(),fitResults.end(),&cmpNoBins);
   TString currentGroup;

   float min1=referenceResult.r1-referenceResult.e1;
   float max1=referenceResult.r1+referenceResult.e1;
   float min2=referenceResult.r2-referenceResult.e2;
   float max2=referenceResult.r2+referenceResult.e2;
   for (FitResult const &fr: fitResults){
      std::string selection=getSelectionString(fr);

      if (fr.distribution.Contains("JESu")) selection+=" (JES up)";
      if (fr.distribution.Contains("JESd")) selection+=" (JES down)";

      // print results
      TString label=(summaryMode==SRT_SEL ? selection+": "+fr.distribution : fr.distribution+": "+selection);
      io::log*label;
      io::log>>TString::Format("%.2f+-%.2f %.2f+-%.2f (rho=%.2f)",fr.r1,fr.e1,fr.r2,fr.e2,fr.rho);

      // vizualize results
      float const y=1-float(iResult+.5)/iResults;
      gr1.SetPoint(iResult,fr.r1,y);
      gr1.SetPointError(iResult,fr.e1,0);
      gr2.SetPoint(iResult,fr.r2,y);
      gr2.SetPointError(iResult,fr.e2,0);

      label=(summaryMode==SRT_SEL ? selection+": "+fr.distLabel : fr.distLabel+": "+selection);
      // if (summaryMode==BINNING) label+=TString::Format(" [%d]",fr.numberOfBins);
      if (summaryMode==BINNING) label=TString::Format("%d bins",fr.numberOfBins);
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
   case CLOSURE: { saveName="closure/"; break; }
   case CONT_GGM:{ saveName="contamin/"; break; }
   }
   switch (summaryMode) {
   case SRT_DST: { saveName+="dist"; break; }
   case SRT_SEL: { saveName+="sel"; break; }
   case BINNING: { saveName+="binning"; break; }
   }
   saver.save(can,saveName+"_1d",true,isSim);

   TCanvas can2;
   gfx::LegendEntries le;
   min1=min1-(max1-min1)*.05;
   max1=max1+(max1-min1)*.05;
   min2=min2-(max2-min2)*.05;
   max2=max2+(max2-min2)*.05;
   TH2F h2Frame("h2Frame",TString::Format(";SF_{%s};SF_{%s}","V(#gamma)","(#gamma+)jets"),10,min1,max1,10,min2,max2);
   h2Frame.Draw("axis");
   RooEllipse ell("",referenceResult.r1,referenceResult.r2,referenceResult.e1,referenceResult.e2,referenceResult.rho,1000);
   ell.SetFillColor(kGray);
   ell.SetLineColor(kGray);
   ell.Draw("f");
   le.append(ell,"68% CL","f");

   TMarker marker;
   marker.SetMarkerStyle(kFullCircle);
   marker.SetMarkerSize(2);
   marker.SetMarkerColor(kGray+1);
   marker.DrawMarker(referenceResult.r1,referenceResult.r2);
   le.prepend(*marker.Clone(),"reference","p");

   int index = 2;
   //systematic studies to check the shape variation
   io::Logger fit_result_weights((std::to_string(index)+"-fitresult.tex"));
   fit_result_weights<< "# weight_index, weight, Vgamma scale, gJets scale";
   TString fill_line = std::to_string(index) +" " + std::to_string(referenceResult.r1)+ " "+ std::to_string(referenceResult.r2) ;
   fit_result_weights<< fill_line;

   marker.SetMarkerColor(kBlack);
   marker.SetMarkerSize(1);
   int markerStyle=23;
   currentGroup="";
   for (FitResult const &fr: fitResults){
      TString const &newGroup=(summaryMode==SRT_SEL ? getSelectionString(fr) : fr.distLabel);
      if (newGroup!=currentGroup){
         markerStyle++;
         marker.SetMarkerStyle(markerStyle);
         le.append(*marker.Clone(),newGroup,"p");
         currentGroup=newGroup;
      }
      marker.DrawMarker(fr.r1,fr.r2);
   }

   TLegend leg=le.buildLegend(.55,.6);
   leg.Draw();
   saver.save(can2,saveName+"_2d",true,isSim);

   // Pulls
   can2.Clear();
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
  std::vector<std::vector<TString>> fitCombinations={{"V(#gamma)","WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"},{"#gamma+jets","GJets_DR"}}; 
//   std::vector<std::vector<TString>> fitCombinations={{"V(+#gamma)","WGToLNuG_SUS","ZGTo2NuG","ZGTo2LG","ZNuNuJets","WLNuJets"},{"#gamma+jets","GJets_DR"}};
 //  std::vector<std::vector<TString>> fitCombinations={{"V(+#gamma)","WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"},{"(#gamma+)jets","QCD","GJets"}};
//   std::vector<std::vector<TString>> fitCombinations={{"V(+#gamma)","WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"},{"(#gamma+)jets","GJets"}};
   // double/half binning
   
   FitResult referenceResult=fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);

   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.1,.2},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.4,.8},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.4,1.2,2.8,3.2},{.1,.4,.8,.4},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.4,1.2,2.8,3.2},{.2,.8,1.6,.4},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.4,1.2,2.8,3.2},{.1,.8,1.6,.4},fitCombinations,fixedBkg,fitMode);
   /*
   fit("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","absphiMETjet",{0,.4,3.2,3.4},{.1,.2,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","absphiMETjet",{0,.4,3.2,3.4},{.1,.4,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","absphiMETjet",{0,.4,.8,3.2,3.4},{.4,.4,.8,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","absphiMETjet",{0,.4,2.2,3.2,3.4},{.4,1.8,1.0,.1},fitCombinations,fixedBkg,fitMode);
*/
   // fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETnJet",{0,.4,2.,3.2,3.4},{.2,.4,.6,.1},fitCombinations,fixedBkg,fitMode);
   // fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETnJet",{0,.4,2.,3.2,3.4},{.1,.2,.3,.1},fitCombinations,fixedBkg,fitMode);
   // fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETnJet",{0,.4,2.,3.2,3.4},{.4,.8,1.2,.1},fitCombinations,fixedBkg,fitMode);

   // fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETph",{0,1,2,2.4,3.2,3.4},{1,.5,.4,.2,.1},fitCombinations,fixedBkg,fitMode);
   // fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETph",{0,1,2,2.4,3.2,3.4},{.5,.25,.2,.1,.1},fitCombinations,fixedBkg,fitMode);
   // fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETph",{0,1.6,2.4,3.2,3.4},{.8,.8,.4,.1},fitCombinations,fixedBkg,fitMode);

   resultSummary(referenceResult,BINNING,fitMode);

   // nominal binning
   g_fitResults.clear();

 //  fit("pre_ph165/c_S30/MT100/Sl80vMTl300/0bL/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
 //  fit("pre_ph165/c_S30/MT100/Sl80vMTl300/0bM/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);

   //2016 studies
  
   referenceResult=fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_phoPtl500",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_phoPtl700",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_phoPtg700",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);
   
 //  fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETph",{0,1,2,2.4,3.2,3.4},{1,.5,.4,.2,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);   
   fit("pre_ph165/c_MET100/MT100/METl400vMTl400/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_JESu",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_JESd",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET150/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET150/MT150/METl400vMTl400/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","STg"    ,{300,500,700,800},{20,40,100},fitCombinations,fixedBkg,fitMode);     
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","relPt2Jets",{0.08,1},{0.08},fitCombinations,fixedBkg,fitMode);
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","DeltaS",{0.2,2.4,3.2},{0.2,0.1},fitCombinations,fixedBkg,fitMode);
 //  fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","DeltaS1",{0,2.4,3.2},{0.2,0.1},fitCombinations,fixedBkg,fitMode);

 
 
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);
     fit("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  
//   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  
 //  fit("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/","absphiMETnJetPh",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);

   
   fit("pre_ph165/c_MET150/MT100/METl300vMTl300/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  

 //  fit("pre_ph165/c_MET150/MT100/METl300vMTl300/","absphiMETph",{0,1,2,2.4,3.2,3.4},{1,.5,.4,.2,.1},fitCombinations,fixedBkg,fitMode);
    
   fit("pre_ph165/c_MET150/MT150/METl400vMTl400/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  

 //  fit("pre_ph165/c_MET150/MT150/METl400vMTl400/","absphiMETph",{0,1,2,2.4,3.2,3.4},{1,.5,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  
   fit("pre_ph165/c_MET100/MT100/METl400vMTl400/","absphiMETjet",{0,1.2,2.8,3.2,3.4},{.2,.4,.2,.1},fitCombinations,fixedBkg,fitMode);  
 //  fit("pre_ph165/c_MET100/MT100/METl400vMTl400/","absphiMETph",{0,1,2,2.4,3.2,3.4},{1,.5,.4,.2,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","MET",{100,310,410},{30,100},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","MT"    ,{100,300,500,600},{40,20,100},fitCombinations,fixedBkg,fitMode);

 

   std::vector<FitResult> fitResults(g_fitResults);

   g_fitResults.clear();

   // fits not in summary
/*
   fit("pre_ph165/c_S30/MT100/Sl100vMTl300/","absphiMETph",{0,1,2,2.4,3.2,3.4},{1,.5,.4,.2,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl100vMTl300/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl100vMTl300/","absphiMETnJet",{0,.4,2.,3.2,3.4},{.2,.4,.6,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl100vMTl300/","absphiMETn2Jet",{0,.4,2.,3.2,3.4},{.2,.4,.6,.1},fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl100vMTl300/","absphiMETnJetPh",{0,.4,2.,2.8,3.2,3.4},{.2,.4,.8,.4,.1},fitCombinations,fixedBkg,fitMode);

 //  fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETjetJESu",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);
 //  fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETjetJESd",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);

   fit("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},fitCombinations,fixedBkg,fitMode);

   fit("pre_ph165/c_S30/MT100/Sl80vMTl300/","METS",1,fitCombinations,fixedBkg,fitMode);
   fit("pre_ph165/c_S30/MT100/Sl100vMTl300/","METS",1,fitCombinations,fixedBkg,fitMode);

*/
   g_fitResults=fitResults;
   resultSummary(referenceResult,SRT_SEL,fitMode);
   resultSummary(referenceResult,SRT_DST,fitMode);

}

struct Chi2Result
{
   double yieldBkg;
   double yieldPred,systErrPred;
   double chi2Ndf,pVal;
};

Chi2Result getSR_chi2(TString sPresel,TString sVar,
                      std::vector<float> rebinEdges,std::vector<float> rebinWidths,
                      std::vector<std::vector<TString>> fitCombinations,
                      std::vector<TString> fixedBkg,
                      float f1, float f2,
                      FitResult fr,
                      FitMode_t fitMode
   )
{
   assert(fitCombinations.size() == 2); // maybe generalize this

   TString saveName="closure_scale/";
   if (fitMode==CONT_GGM) saveName="contamin_scale/";

   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   Rebinner rebinned(rebinEdges,rebinWidths);

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

   TH1F hTotal;
   // add up all MC
   hTotal=combHists["fixed"];
   {
      TH1F hToAdd1(combHists[sFitSample1]);
      TH1F hToAdd2(combHists[sFitSample2]);
      if (f1>0) hToAdd1.Scale(f1);
      if (f2>0) hToAdd2.Scale(f2);
      hTotal.Add(&hToAdd1);
      hTotal.Add(&hToAdd2);
   }

   TH1F hBkg(hTotal);

   TString dataLabel,totalLabel,predLabel;
   if (fitMode==CLOSURE){
      dataLabel="toy data";
      totalLabel="#Sigma B";
      predLabel="prediction";
   } else if (fitMode==CONT_GGM){
      dataLabel="signal injected toy data";
      totalLabel="#Sigma B + S";
      predLabel="signal contaminated prediction";
      // inject signal
      TH1F hSignal=rebinned(*(TH1F*)histReader.read<TH1F>(sPresel+sVar+"/TChiWG"));
      hSignal.Scale(cfg.trigger_eff_Ph);
      hTotal.Add(&hSignal);
   }

   // dice from a poissonian for each bin
   TH1F hData(hTotal);
   for (int i=0; i<=hData.GetNbinsX()+1;i++){
      int const v=g_rand.Poisson(hData.GetBinContent(i));
      hData.SetBinContent(i,v);
      hData.SetBinError(i,TMath::Sqrt(v));
   }

   TH1F hPred;
   TH1F hSyst;
   float totalSystErr;
   // add up all MC
   hPred=combHists["fixed"];
   {
      TH1F hToAdd1(combHists[sFitSample1]);
      TH1F hToAdd2(combHists[sFitSample2]);
      hToAdd1.Scale(fr.r1);
      hToAdd2.Scale(fr.r2);
      hPred.Add(&hToAdd1);
      hPred.Add(&hToAdd2);
      float eVg,eGJ;
      float erri;
      eVg=hToAdd1.Integral()*fr.e1/fr.r1;
      eGJ=hToAdd2.Integral()*fr.e2/fr.r2;
      erri=util::quadSum<double>({eVg,eGJ});
      erri+=2*fr.rho*eVg*eGJ; // correlation
      totalSystErr=TMath::Sqrt(erri);
      hSyst=hPred;
      for (int i=0; i<=hSyst.GetNbinsX()+1; i++) {
         eVg=hToAdd1.GetBinContent(i)*fr.e1/fr.r1;
         eGJ=hToAdd2.GetBinContent(i)*fr.e2/fr.r2;
         erri=util::quadSum<double>({eVg,eGJ});
         erri+=2*fr.rho*eVg*eGJ; // correlation
         hSyst.SetBinError(i,TMath::Sqrt(erri));
      }
   }

   TCanvas can;
   can.SetLogy();
   gfx::LegendEntries le;
   hData.SetLineColor(kBlack);
   hData.SetMarkerStyle(kOpenCircle);
   hData.Draw();
   hSyst.Draw("same e2");
   hSyst.SetFillColor(kGray+1);
   hSyst.SetLineColor(kGray+1);
   hSyst.SetMarkerSize(0);
   if (fitMode==CONT_GGM) {
      hBkg.Draw("same hist");
      hBkg.SetFillStyle(0);
      hBkg.SetLineStyle(2);
      le.append(hBkg,"#Sigma B","l");
      hTotal.SetLineColor(kRed);
   }
   hTotal.Draw("same hist");
   hTotal.SetFillStyle(0);
   hTotal.SetLineStyle(2);
   hPred.Draw("same hist e");
   hPred.SetLineColor(kBlue);
   hPred.SetLineWidth(2);
   hPred.SetFillStyle(0);
   hPred.SetMarkerSize(0);
   hData.Draw("same pe1");

   le.append(hTotal,totalLabel,"l");
   le.append(hData,dataLabel,"pe");
   le.append(hPred,predLabel,"le");
   le.append(hSyst,"#pm#sigma_{syst} (fit)","f");
   le.buildLegend(.51,.75).DrawClone();
   can.RedrawAxis();
   saver.save(can,TString::Format(saveName+"dist/%.1f_%.1f",f1,f2),true,true);

   // add systematic uncertainty

   float erri;
   for (int i=0; i<=hSyst.GetNbinsX()+1; i++) {          
      erri=util::quadSum<double>({hPred.GetBinError(i),hSyst.GetBinError(i)});
      hPred.SetBinError(i,TMath::Sqrt(erri));
   }
 //  debug << "Chi2Result    " << hData.Chi2Test(&hPred,"UW CHI2/NDF");
   return Chi2Result{hBkg.Integral(),hPred.Integral(),totalSystErr,hData.Chi2Test(&hPred,"UW CHI2/NDF"),hData.Chi2Test(&hPred,"UW")};
}

void runClosureRange(FitMode_t fitMode, int Nrepetitions=1)
{
   assert(fitMode==CLOSURE || fitMode==CONT_GGM);

   TString saveName="closure_scale/";
   if (fitMode==CONT_GGM) saveName="contamin_scale/";

   std::vector<TString> fixedBkg={"diboson","TTJets","TTGJets","efake"};
 //  std::vector<std::vector<TString>> fitCombinations={{"V(+#gamma)","WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"},{"(#gamma+)jets","QCD","GJets"}};
   std::vector<std::vector<TString>> fitCombinations={{"V(+#gamma)","WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"},{"(#gamma+)jets","GJets_DR"}};

   TGraphErrors gr1_SFVg;
   TGraphErrors gr1_SFGJ;
   TGraph       gr1_total;
   TGraphErrors gr1_pred;
   TGraph       gr1_chi2Ndf;
   TH1F         h1_chi2Ndf("",";#chi^{2}/N_{df};Entries / bin",5,0,5);
   TH1F         h1_pVal("",";p-value;Entries / bin",5,0,1);
   gr1_SFVg.SetMarkerColor(fillcolor1);
   gr1_SFVg.SetLineColor(fillcolor1);
   gr1_SFVg.SetLineWidth(2);
   gr1_SFVg.SetMarkerStyle(kOpenCircle);
   gr1_SFGJ.SetMarkerColor(fillcolor2);
   gr1_SFGJ.SetLineColor(fillcolor2);
   gr1_SFGJ.SetLineWidth(2);
   gr1_SFGJ.SetMarkerStyle(kOpenSquare);
   gr1_chi2Ndf.SetTitle(";f_{V(+#gamma)};");
   gr1_chi2Ndf.SetMarkerStyle(kOpenCross);
   gr1_pred.SetMarkerStyle(kOpenCircle);
   TGraphErrors gr2_SFVg(gr1_SFVg);
   TGraphErrors gr2_SFGJ(gr1_SFGJ);
   TGraph       gr2_total;
   TGraphErrors gr2_pred(gr1_pred);
   TGraph       gr2_chi2Ndf(gr1_chi2Ndf);
   TH1F         h2_chi2Ndf("",";#chi^{2}/N_{df};Entries / bin",5,0,5);
   TH1F         h2_pVal("",";p-value;Entries / bin",5,0,1);
   gr2_chi2Ndf.SetTitle(";f_{(#gamma+)jets};");

   TGraphErrors gr1_SFVg_mean(gr1_SFVg);
   TGraphErrors gr1_SFGJ_mean(gr2_SFGJ);
   TGraphErrors gr2_SFVg_mean(gr1_SFVg);
   TGraphErrors gr2_SFGJ_mean(gr2_SFGJ);


   FitResult fr;
   int iPoint=-1;
   float f_start=0.5;
   float f=f_start;
   Chi2Result chi2result;
   for (; f<=2.0; f+=0.04) {
      iPoint++;
      float meanW1=0;
      float sumW1=0;
      float meanW2=0;
      float sumW2=0;
      for (int iRepeat=0; iRepeat<Nrepetitions; iRepeat++) {
         fr=fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode,f,-1);
         // last repetition will be shown in these graphs:
         gr1_SFVg.SetPoint(iPoint,f,fr.r1);
         gr1_SFVg.SetPointError(iPoint,0,fr.e1);
         gr1_SFGJ.SetPoint(iPoint,f,fr.r2);
         gr1_SFGJ.SetPointError(iPoint,0,fr.e2);
         meanW1+=1/fr.e1/fr.e1*fr.r1;
         sumW1+=1/fr.e1/fr.e1;
         meanW2+=1/fr.e2/fr.e2*fr.r2;
         sumW2+=1/fr.e2/fr.e2;
         chi2result=getSR_chi2("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},fitCombinations,fixedBkg,f,-1,fr,fitMode);
         gr1_total.SetPoint(iPoint,f,chi2result.yieldBkg);
         gr1_pred.SetPoint(iPoint,f,chi2result.yieldPred);
         gr1_pred.SetPointError(iPoint,0,chi2result.systErrPred);
         gr1_chi2Ndf.SetPoint(iPoint,f,chi2result.chi2Ndf);
         h1_chi2Ndf.Fill(chi2result.chi2Ndf);
         h1_pVal.Fill(chi2result.pVal);
      }
      meanW1/=sumW1;
      meanW2/=sumW2;
      gr1_SFVg_mean.SetPoint(iPoint,f,meanW1);
      gr1_SFVg_mean.SetPointError(iPoint,0,TMath::Sqrt(1./sumW1));
      gr1_SFGJ_mean.SetPoint(iPoint,f,meanW2);
      gr1_SFGJ_mean.SetPointError(iPoint,0,TMath::Sqrt(1./sumW2));

      meanW1=0;
      sumW1=0;
      meanW2=0;
      sumW2=0;
      for (int iRepeat=0; iRepeat<Nrepetitions; iRepeat++) {
         fr=fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode,-1,f);
         // last repetition will be shown in these graphs:
         gr2_SFVg.SetPoint(iPoint,f,fr.r1);
         gr2_SFVg.SetPointError(iPoint,0,fr.e1);
         gr2_SFGJ.SetPoint(iPoint,f,fr.r2);
         gr2_SFGJ.SetPointError(iPoint,0,fr.e2);
         meanW1+=1/fr.e1/fr.e1*fr.r1;
         sumW1+=1/fr.e1/fr.e1;
         meanW2+=1/fr.e2/fr.e2*fr.r2;
         sumW2+=1/fr.e2/fr.e2;
         chi2result=getSR_chi2("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},fitCombinations,fixedBkg,-1,f,fr,fitMode);
         gr2_total.SetPoint(iPoint,f,chi2result.yieldBkg);
         gr2_pred.SetPoint(iPoint,f,chi2result.yieldPred);
         gr2_pred.SetPointError(iPoint,0,chi2result.systErrPred);
         gr2_chi2Ndf.SetPoint(iPoint,f,chi2result.chi2Ndf);
         h2_chi2Ndf.Fill(chi2result.chi2Ndf);
         h2_pVal.Fill(chi2result.pVal);
      }
      meanW1/=sumW1;
      meanW2/=sumW2;
      gr2_SFVg_mean.SetPoint(iPoint,f,meanW1);
      gr2_SFVg_mean.SetPointError(iPoint,0,TMath::Sqrt(1./sumW1));
      gr2_SFGJ_mean.SetPoint(iPoint,f,meanW2);
      gr2_SFGJ_mean.SetPointError(iPoint,0,TMath::Sqrt(1./sumW2));
   }
   TCanvas can;
   gfx::LegendEntries le1;
   le1.append(gr1_SFVg,"SF_{V(#gamma)}","pe");
   le1.append(gr1_SFGJ,"SF_{(#gamma+)jets}","pe");
   le1.append(gr1_chi2Ndf,"SR #chi^{2}/N_{df}","p");
   gfx::LegendEntries le2;
   le2.append(gr1_total,"true bkg (SR)","l");
   le2.append(gr1_pred,"prediction","pe");
   TLine line;

   gr1_chi2Ndf.Draw("ap");
   gr1_SFVg.Draw("p");
   gr1_SFGJ.Draw("p");
   gr1_chi2Ndf.Draw("p");
   line.SetLineColor(color1);
   line.DrawLine(.9*f_start,.9*f_start,f,f);
   line.SetLineColor(color2);
   line.DrawLine(.9*f_start,1,f,1);
   le1.buildLegend(.4,.75).DrawClone();
   saver.save(can,saveName+"f1",true,true);

   gr1_chi2Ndf.Draw("ap");
   gr1_SFVg_mean.Draw("p");
   gr1_SFGJ_mean.Draw("p");
   gr1_chi2Ndf.Draw("p");
   line.SetLineColor(color1);
   line.DrawLine(.9*f_start,.9*f_start,f,f);
   line.SetLineColor(color2);
   line.DrawLine(.9*f_start,1,f,1);
   le1.buildLegend(.4,.75).DrawClone();
   saver.save(can,saveName+"f1_mean",true,true);

   h1_chi2Ndf.Draw("hist");
   h1_chi2Ndf.SetMinimum(0);
   gfx::cornerLabel(TString::Format("mean=%.2f",h1_chi2Ndf.GetMean()),2).DrawClone();
   saver.save(can,saveName+"chi2Ndf1",true,true);
   h1_pVal.Draw("hist");
   h1_pVal.SetMinimum(0);
   saver.save(can,saveName+"pVal1",true,true);

   gr1_total.Draw("al");
   gr1_total.SetTitle(";f_{V(+#gamma)};");
   gr1_total.GetYaxis()->SetRangeUser(100,450);
   gr1_pred.Draw("pe");
   le2.buildLegend(.4,.75).DrawClone();
   saver.save(can,saveName+"pred1",true,true);

   gr2_chi2Ndf.Draw("ap");
   gr2_SFVg.Draw("p");
   gr2_SFGJ.Draw("p");
   gr2_chi2Ndf.Draw("p");
   line.SetLineColor(color1);
   line.DrawLine(.9*f_start,1,f,1);
   line.SetLineColor(color2);
   line.DrawLine(.9*f_start,.9*f_start,f,f);
   le1.buildLegend(.4,.75).DrawClone();
   saver.save(can,saveName+"f2",true,true);

   gr2_chi2Ndf.Draw("ap");
   gr2_SFVg_mean.Draw("p");
   gr2_SFGJ_mean.Draw("p");
   gr2_chi2Ndf.Draw("p");
   line.SetLineColor(color1);
   line.DrawLine(.9*f_start,1,f,1);
   line.SetLineColor(color2);
   line.DrawLine(.9*f_start,.9*f_start,f,f);
   le1.buildLegend(.4,.75).DrawClone();
   saver.save(can,saveName+"f2_mean",true,true);

   h2_chi2Ndf.Draw("hist");
   h2_chi2Ndf.SetMinimum(0);
   gfx::cornerLabel(TString::Format("mean=%.2f",h2_chi2Ndf.GetMean()),2).DrawClone();
   saver.save(can,saveName+"chi2Ndf2",true,true);
   h2_pVal.Draw("hist");
   h2_pVal.SetMinimum(0);
   saver.save(can,saveName+"pVal2",true,true);

   gr2_total.Draw("al");
   gr2_total.SetTitle(";f_{(#gamma+)jets};");
   gr2_total.GetYaxis()->SetRangeUser(100,450);
   gr2_pred.Draw("pe");
   le2.buildLegend(.4,.75).DrawClone();
   saver.save(can,saveName+"pred2",true,true);

   // single example
   fr=fit("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.2,.4},fitCombinations,fixedBkg,fitMode,.5,1.5);
   chi2result=getSR_chi2("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},fitCombinations,fixedBkg,.5,1.5,fr,fitMode);
   io::log*"Chi2/Ndf=">>chi2result.chi2Ndf;
}

extern "C"
void run()
{
 //  run(CLOSURE);
 //  run(CONT_GGM);
   run(DATA);

//   g_rand.SetSeed(87234);
//   runClosureRange(CLOSURE,10);
//   g_rand.SetSeed(87234);
 //  runClosureRange(CONT_GGM);
}

// plotting directly from rooplot
   // RooPlot* frame = x.frame();
   // dhData.plotOn(frame);
   // model.plotOn(frame,RooFit::DrawOption("F"),RooFit::FillColor(kGray));
   // model.plotOn(frame,RooFit::Components(RooArgSet(pdfSample1)),RooFit::LineColor(kBlack));
   // // model.plotOn(frame,RooFit::Components(RooArgSet(pdfSample2)) ,RooFit::LineColor(kOrange-3));
   // model.paramOn(frame);
   // dhData.plotOn(frame);
   // frame->Draw();


// re-extract histogram
// TH1* hData2=dhData.createHistogram("name",x);
// hData2->Draw();
// does not look right:
// model.createHistogram("name",x,RooFit::Components(RooArgSet(pdfSample2)))->Draw();
