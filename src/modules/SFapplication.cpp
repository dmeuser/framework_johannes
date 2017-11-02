#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

#include <limits>

#include <TFractionFitter.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TStyle.h>

#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooEllipse.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

Config const &cfg=Config::get();
static io::RootFileSaver saver("plots.root","SFapplication");
static io::Logger logFile("SFapplication.log");


float Vg_JES_bin1 = 0.053425;
float Vg_JES_bin2 = 0.048425;
float Vg_JES_bin3 = 0.039086;
float Vg_JES_bin4 = 0.046930;

float GJ_JES_bin1 = 0.022336;
float GJ_JES_bin2 = 0.043804;
float GJ_JES_bin3 = 0.043238;
float GJ_JES_bin4 = 0.449069;

float Vg_mu_bin1 = 0.045273;
float Vg_mu_bin2 = 0.064142;
float Vg_mu_bin3 = 0.074932;
float Vg_mu_bin4 = 0.098386;

float Vg_pdf_bin1 = 0.0165178;
float Vg_pdf_bin2 = 0.0292389;
float Vg_pdf_bin3 = 0.0339644;
float Vg_pdf_bin4 = 0.0530947;

float GJ_mu_bin1 = 0.043224;
float GJ_mu_bin2 = 0.037132;
float GJ_mu_bin3 = 0.054357;
float GJ_mu_bin4 = 0.060463;

float GJ_pdf_bin1 = 0.0385493;
float GJ_pdf_bin2 = 0.0175165;
float GJ_pdf_bin3 = 0.0707174;
float GJ_pdf_bin4 = 0.0648986;

class Rebinner
{
public:
   // if iRebin==0: use vectors
   Rebinner(int iRebin, bool bDivideByBinWidth=false)
      : iRebin_(iRebin)
      , bDivideByBinWidth_(bDivideByBinWidth)
      {}
   Rebinner(std::vector<float> rebinEdges,std::vector<float> rebinWidths, bool bDivideByBinWidth=false)
      : iRebin_(0)
      , bDivideByBinWidth_(bDivideByBinWidth)
      {
         newBinEdges_=hist::getBinVector(rebinEdges,rebinWidths);
      }
   TH1F operator()(TH1F const &h) const {
      TH1F hNew;
      if (iRebin_!=0){
         hNew=h;
         hNew.Rebin(iRebin_);
      } else {
         hNew=TH1F(h);
         hNew=hist::rebinned(h,newBinEdges_,true,false); // merging overflow, not underflow
      }
      if (bDivideByBinWidth_) {
         hist::divideByBinWidth(hNew);
      }
      return hNew;
   }
   int numberOfBins(){
      if (iRebin_==0) return newBinEdges_.size()-1;
      return 0;
   }
   float binEdge(int i) { // using root bin numbers
      assert(iRebin_==0);
      i--; // using root bin numbers
      if (i==-1) return -std::numeric_limits<float>::infinity();
      if ((unsigned)i==newBinEdges_.size()) return std::numeric_limits<float>::infinity();
      return newBinEdges_.at(i);
   }
private:
   int iRebin_;
   bool bDivideByBinWidth_;
   std::vector<double> newBinEdges_;
};


TH1F combineHists(io::RootFileReader const& histReader, TString sPresel, TString sVar, std::vector<TString> sampleNames)
{
   assert(sampleNames.size()>0);
   TH1F const *hFirst=histReader.read<TH1F>(sPresel+sVar+"/"+sampleNames[0]);
   TH1F hComb(*hFirst);
   hComb.SetLineColor(hFirst->GetLineColor());
   hComb.SetFillColor(hFirst->GetLineColor());
   for (unsigned i=1; i<sampleNames.size(); i++){
      hComb.Add(histReader.read<TH1F>(sPresel+sVar+"/"+sampleNames[i]));
   }
   return hComb;
}

void setupHistForStack(TH1 &h){
   h.SetFillStyle(1001);
   h.SetFillColor(h.GetLineColor());
   h.SetLineColor(kBlack);
   h.SetLineWidth(1);
}

struct SignalBin
{
   struct Component
   {
      Component(){}
      float count=0;
      float estat=0;
      float esyst=0;
      TString formatLine(float nTotalBkg=-1,bool tex=false, bool data=false,float add=-1) const {
         if (nTotalBkg<0) nTotalBkg=count;
         float esy=esyst;
         if (add>0) esy=TMath::Sqrt(esy*esy+add*add*count*count);
         if (!tex) return TString::Format("%6.2f %6.2f %6.2f(%4.1f%%)",
                                          count,estat,esy,esy/nTotalBkg*100);
         TString countFormat=" & %6.2f";
         if (data) countFormat=" & %6.0f";
         if (esy>0 && estat>0) return TString::Format(countFormat+" & %6.2f & %6.2f & %4.1f\\%%\\\\",
                                                        count,estat,esy,esy/nTotalBkg*100);
         if (estat>0) return TString::Format(countFormat+" & %6.2f &  & \\\\",
                                             count,estat);
         return TString::Format(countFormat+" &  &  & \\\\",count);
      }
   };
   float lower,upper; // edges
   std::map<TString,Component> component;

   TString formatBlock(bool tex=false,bool showTotal=false) {
      float nBkgTot=getTotal().count;
      TString block;
      for (TString compName: {"Vg","GJ","TTcomb","efake","diboson","TChiWG","T5Wg"}) {
         TString line;
         if (component.count(compName)) {
            if (compName=="GJ") line=component[compName].formatLine(nBkgTot,tex,false,cfg.sf.uncert_gammaJ);
            if (compName=="Vg") line=component[compName].formatLine(nBkgTot,tex,false,cfg.sf.uncert_Vgamma);           
            else line=component[compName].formatLine(nBkgTot,tex);
         } else {
            continue;
         }
         if (tex) {
            compName.ReplaceAll("5","five");
            compName="\\"+compName;
         }
         block+=TString::Format("   %-10s ",compName.Data());
         block+=line;
         block+="\n";
      }
      if (showTotal) {
         block+=TString::Format("   %-10s ",tex?"\\total":"total");
         block+=getTotal().formatLine(-1,tex);
      }
      if (component.count("data")) {
         block+=TString::Format("\n   %-10s ","data");
         block+=component["data"].formatLine(-1,tex,true);
      }
      return block;
   }
   Component getTotal() const {
      Component tot;
      for (auto const &b: component) {
         if (b.first=="data") continue;
         tot.count+=b.second.count;
         tot.estat+=b.second.estat*b.second.estat;
         tot.esyst+=b.second.esyst*b.second.esyst;
      }
      if (component.count("Vg") && component.count("GJ")) {
         // background: consider correlation
         tot.esyst+=2*cfg.sf.rho
            *component.at("Vg").esyst
            *component.at("GJ").esyst;
         tot.esyst+=(cfg.sf.uncert_gammaJ*component.at("GJ").count*cfg.sf.uncert_gammaJ*component.at("GJ").count);
         tot.esyst+=(cfg.sf.uncert_Vgamma*component.at("Vg").count*cfg.sf.uncert_Vgamma*component.at("Vg").count);         
      }
      tot.estat=TMath::Sqrt(tot.estat);
      tot.esyst=TMath::Sqrt(tot.esyst);
      return tot;
   }
};

enum PlotMode_t
{
   PRE, // preselection
   CR, // control region
   SR_BLIND, // signal region without data
   SR_CHECK, // signal region control distributions with data
   SR, // signal region unblinded
   VR,
};

// General implementation of plot. Call on of the functions below, not this one.
void plot(TString sSelection,TString sVar,int iRebin,
          std::vector<float> rebinEdges,std::vector<float> rebinWidths,
          PlotMode_t plotMode
   )
{
   bool const showData=(plotMode!=SR_BLIND && plotMode!=PRE);
   bool const showSignal=(plotMode!=CR && plotMode!=VR);
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   TString saveName=sSelection+sVar;

   bool divideByBinWidth=(plotMode==CR || plotMode==VR || plotMode==PRE);
   Rebinner rebinned(iRebin,divideByBinWidth);
   if (iRebin==0) rebinned=Rebinner(rebinEdges,rebinWidths,divideByBinWidth);

   bool const customBinning = (iRebin==0 && plotMode!=CR && plotMode!=VR);
   if (customBinning) saveName+="_binning";
   if (!showData)     saveName+="_blind";

   std::vector<SignalBin> signalBins;
   std::vector<SignalBin> signalBins_signal;
   SignalBin bkgIntegral; // for the whole signal region
   SignalBin signalIntegral;
   if (plotMode==SR) {
      for (int i=0; i<rebinned.numberOfBins()+2; i++) {
         SignalBin bin;
         bin.lower=rebinned.binEdge(i);
         bin.upper=rebinned.binEdge(i+1);
         signalBins.push_back(bin);
         signalBins_signal.push_back(bin);
      }
      bkgIntegral.lower=bkgIntegral.upper=0;
      signalIntegral.lower=bkgIntegral.upper=0;
   }

//   std::vector<TString> vsVgComponents{"WGToLNuG_SUS","ZGTo2NuG","ZGTo2LG","ZNuNuJets","WLNuJets"};
   std::vector<TString> vsVgComponents{"WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"};
   std::vector<TString> vsGJComponents{"GJets_DR"};
  // std::vector<TString> vsGJComponents{"GJets_DR"};

   TH1F hVg=rebinned(combineHists(histReader,sSelection,sVar,vsVgComponents));
   TH1F hGJ=rebinned(combineHists(histReader,sSelection,sVar,vsGJComponents));
   hVg.Scale(cfg.trigger_eff_Ph*cfg.sf.Vg);
   hGJ.Scale(cfg.trigger_eff_Ph*cfg.sf.GJ);

   io::log<<TString::Format("N_Vg=%.2f",hVg.Integral(0,hVg.GetNbinsX()+1,divideByBinWidth ? "width" : ""));
   io::log<<TString::Format("N_GJ=%.2f",hGJ.Integral(0,hGJ.GetNbinsX()+1,divideByBinWidth ? "width" : ""));

   hVg.SetLineColor(kOrange-3);
   hGJ.SetLineColor(kAzure+10);
   setupHistForStack(hVg);
   setupHistForStack(hGJ);

   gfx::LegendEntries le;
   THStack st;
   std::map<TString,TH1F> fixHists;
   for (TString sSample:{"efake","TTJets","TTGJets","diboson","TChiWG","T5Wg"}){
      fixHists[sSample]=rebinned(*(TH1F*)histReader.read<TH1F>(sSelection+sVar+"/"+sSample));
      TH1 &h=fixHists[sSample];
      h.Scale(cfg.trigger_eff_Ph);
      setupHistForStack(h);
      if (sSample == "efake"){
         h.SetFillColor(kMagenta+2);
      }
      else if (sSample == "TTGJets"){
         h.SetFillColor(kBlue+2);
      }
      else if (sSample == "diboson"){
         h.SetFillColor(kMagenta+4);
      }
   }
   fixHists["efake"].Scale(1./cfg.trigger_eff_Ph); // efake is data-driven

   fixHists["TTcomb"]=fixHists["TTGJets"];
   fixHists["TTcomb"].Add(&fixHists["TTJets"]);
   st.Add(&fixHists["diboson"],"hist");
   st.Add(&fixHists["efake"],"hist");
   st.Add(&fixHists["TTcomb"],"hist");
   le.prepend(fixHists["diboson"],cfg.datasets.getLabel("diboson"),"f");
   le.prepend(fixHists["efake"],cfg.datasets.getLabel("efake"),"f");
   le.prepend(fixHists["TTcomb"],"t#bar{t}(+#gamma)","f");

   if (plotMode==CR) {
      // CR: more GJet
      st.Add(&hVg,"hist");
      st.Add(&hGJ,"hist");
      le.prepend(hVg,"V(+#gamma)","f");
      le.prepend(hGJ,"#gamma+jets","f");
   } else {
      // SR,VR: more Vg
      st.Add(&hGJ,"hist");
      st.Add(&hVg,"hist");
      le.prepend(hGJ,"#gamma+jets","f");
      le.prepend(hVg,"V(+#gamma)","f");
   }
   TH1F hStackSum(*(TH1F*)st.GetStack()->Last());
   for (TString sSample:{"TChiWG","T5Wg"}) {
      fixHists[sSample+"_stacked"]=fixHists[sSample];
      TH1F &h=fixHists[sSample+"_stacked"];
      h.Add(&hStackSum);
      h.SetLineColor(h.GetFillColor());
      h.SetLineWidth(2);
      h.SetFillStyle(0);
   }

   TH1F hData;
   if (showData){
      hData=rebinned(*histReader.read<TH1F>(sSelection+sVar+"/SinglePhoton"));
      le.prepend(hData,"data","pe");
   }
   TH1F hStatErr(hStackSum);
   hStatErr.SetMarkerSize(0);
   TH1F hSystErr(hStatErr);
   float erri, error_scale_Vg, error_scale_GJ, error_corr;
   float eVg,eGJ,eVg_mu,eGJ_mu,eVg_pdf,eGJ_pdf;
   float eVg_JES,eGJ_JES;
   float comb_uncert;
   // systematic uncert. from fit
   for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
      // total uncertainties:
      eVg=hVg.GetBinContent(i)*cfg.sf.e_Vg/cfg.sf.Vg;
      eGJ=hGJ.GetBinContent(i)*cfg.sf.e_GJ/cfg.sf.GJ;
      error_corr=2*cfg.sf.rho*eVg*eGJ;
      if (plotMode == SR  && sVar=="STg"){         
         if (i == 1){
            eVg_mu=hVg.GetBinContent(i)*Vg_mu_bin1;
            eVg_pdf=hVg.GetBinContent(i)*Vg_pdf_bin1;
            eVg_JES=hVg.GetBinContent(i)*Vg_JES_bin1;            
            error_scale_Vg=util::sqrtQuadSum<double>({eVg_mu,eVg_pdf});
            error_scale_Vg=util::sqrtQuadSum<double>({error_scale_Vg,eVg_JES});
            eVg = util::quadSum<double>({eVg,error_scale_Vg});
   
            eGJ_mu=hGJ.GetBinContent(i)*GJ_mu_bin1;
            eGJ_pdf=hGJ.GetBinContent(i)*GJ_pdf_bin1;
            eGJ_JES=hGJ.GetBinContent(i)*GJ_JES_bin1; 
            error_scale_GJ=util::sqrtQuadSum<double>({eGJ_mu,eGJ_pdf});
            error_scale_GJ=util::sqrtQuadSum<double>({error_scale_GJ,eGJ_JES});
            eGJ = util::quadSum<double>({eGJ,error_scale_GJ});
         }else if (i ==2) {
            eVg_mu=hVg.GetBinContent(i)*Vg_mu_bin2;
            eVg_pdf=hVg.GetBinContent(i)*Vg_pdf_bin2;
            eVg_JES=hVg.GetBinContent(i)*Vg_JES_bin2; 
            error_scale_Vg=util::sqrtQuadSum<double>({eVg_mu,eVg_pdf});
            error_scale_Vg=util::sqrtQuadSum<double>({error_scale_Vg,eVg_JES});
            eVg = util::quadSum<double>({eVg,error_scale_Vg});
   
            eGJ_mu=hGJ.GetBinContent(i)*GJ_mu_bin2;
            eGJ_pdf=hGJ.GetBinContent(i)*GJ_pdf_bin2;
            eGJ_JES=hGJ.GetBinContent(i)*GJ_JES_bin2;
            error_scale_GJ=util::sqrtQuadSum<double>({eGJ_mu,eGJ_pdf});
            error_scale_GJ=util::sqrtQuadSum<double>({error_scale_GJ,eGJ_JES});
            eGJ = util::quadSum<double>({eGJ,error_scale_GJ});
         }else if (i ==3) {
            eVg_mu=hVg.GetBinContent(i)*Vg_mu_bin3;
            eVg_pdf=hVg.GetBinContent(i)*Vg_pdf_bin3;
            eVg_JES=hVg.GetBinContent(i)*Vg_JES_bin3;
            error_scale_Vg=util::sqrtQuadSum<double>({eVg_mu,eVg_pdf});
            error_scale_Vg=util::sqrtQuadSum<double>({error_scale_Vg,eVg_JES});
            eVg = util::quadSum<double>({eVg,error_scale_Vg});
   
            eGJ_mu=hGJ.GetBinContent(i)*GJ_mu_bin3;
            eGJ_pdf=hGJ.GetBinContent(i)*GJ_pdf_bin3;
            eGJ_JES=hGJ.GetBinContent(i)*GJ_JES_bin3;
            error_scale_GJ=util::sqrtQuadSum<double>({eGJ_mu,eGJ_pdf});
            error_scale_GJ=util::sqrtQuadSum<double>({error_scale_GJ,eGJ_JES});
            eGJ = util::quadSum<double>({eGJ,error_scale_GJ});
         }else if (i ==4) {
            eVg_mu=hVg.GetBinContent(i)*Vg_mu_bin4;
            eVg_pdf=hVg.GetBinContent(i)*Vg_pdf_bin4;
            eVg_JES=hVg.GetBinContent(i)*Vg_JES_bin4;
            error_scale_Vg=util::sqrtQuadSum<double>({eVg_mu,eVg_pdf});
            error_scale_Vg=util::sqrtQuadSum<double>({error_scale_Vg,eVg_JES});
            eVg = util::quadSum<double>({eVg,error_scale_Vg});
   
            eGJ_mu=hGJ.GetBinContent(i)*GJ_mu_bin4;
            eGJ_pdf=hGJ.GetBinContent(i)*GJ_pdf_bin4;
            eGJ_JES=hGJ.GetBinContent(i)*GJ_JES_bin4;
            error_scale_GJ=util::sqrtQuadSum<double>({eGJ_mu,eGJ_pdf});
            error_scale_GJ=util::sqrtQuadSum<double>({error_scale_GJ,eGJ_JES});
            eGJ = util::quadSum<double>({eGJ,error_scale_GJ});
         }
      }else{         
         error_scale_Vg=hVg.GetBinContent(i)*cfg.sf.uncert_Vgamma;
         eVg = util::quadSum<double>({eVg,error_scale_Vg});
         error_scale_GJ=hGJ.GetBinContent(i)*cfg.sf.uncert_gammaJ; 
         eGJ = util::quadSum<double>({eGJ,error_scale_GJ});
      }         
      erri=eVg + eGJ; //components are already squared
      std::cout << "bin number:  " << i << std::endl;
      std::cout << "yield and syst uncert Vgamma:  " << hVg.GetBinContent(i) << "$\%$" << TMath::Sqrt(eVg) << std::endl;
      std::cout << "yield and syst uncert gammaJets:  " << hGJ.GetBinContent(i) << "$\%$" << TMath::Sqrt(eGJ) << std::endl;
      std::cout << "syst error squared:  " << erri << "   correlation error squared:  " << error_corr << std::endl;
      std::cout << "correlation error fraction:  " << error_corr/erri << std::endl;
      erri+=error_corr;
      hSystErr.SetBinError(i,TMath::Sqrt(erri));
   }
   // syst. uncert. from fixed backgrounds
   for (TString sSample:{"efake","TTcomb","diboson"}){
      TH1F const &h=fixHists[sSample];
      float const uncert=cfg.datasets.getSystUncert(sSample);
      for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
         erri=h.GetBinContent(i)*uncert;
         std::cout << "bin number:  " << i << std::endl;
         if (sSample == "efake"){
            std::cout << "yield and syst uncert efake:  " << h.GetBinContent(i) << "$\%$" << erri << std::endl;
            }
         else if (sSample == "TTcomb"){
            std::cout << "yield and syst uncert TTcomb:  " << h.GetBinContent(i) << "$\%$" << erri << std::endl;
            }
         else if (sSample == "diboson"){
            std::cout << "yield and syst uncert diboson:  " << h.GetBinContent(i) << "$\%$" << erri << std::endl;
            }
         erri=util::sqrtQuadSum<double>({erri,hSystErr.GetBinError(i)});
         hSystErr.SetBinError(i,erri);
         comb_uncert = util::sqrtQuadSum<double>({hStatErr.GetBinError(i),hSystErr.GetBinError(i)});
         hStatErr.SetBinError(i,comb_uncert);
      }
   }

   // store separate background uncertainties
   if (plotMode==SR) {
      for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
         std::map<TString,SignalBin::Component> &bini=signalBins[i].component;
         bini["Vg"].count=hVg.GetBinContent(i);
         bini["Vg"].estat=hVg.GetBinError(i);
         bini["Vg"].esyst=hVg.GetBinContent(i)*cfg.sf.e_Vg/cfg.sf.Vg;

         bini["GJ"].count=hGJ.GetBinContent(i);
         bini["GJ"].estat=hGJ.GetBinError(i);
         bini["GJ"].esyst=hGJ.GetBinContent(i)*cfg.sf.e_GJ/cfg.sf.GJ;

         bini["data"].count=hData.GetBinContent(i);
         bini["data"].estat=bini["data"].esyst=0;

         for (TString sSample:{"efake","TTcomb","diboson"}){
            bini[sSample].count=fixHists[sSample].GetBinContent(i);
            bini[sSample].estat=fixHists[sSample].GetBinError(i);
            bini[sSample].esyst=fixHists[sSample].GetBinContent(i)*cfg.datasets.getSystUncert(sSample);
         }
         for (TString sSample:{"TChiWG","T5Wg"}) {
            signalBins_signal[i].component[sSample].count=fixHists[sSample].GetBinContent(i);
            signalBins_signal[i].component[sSample].estat=fixHists[sSample].GetBinError(i);
            signalBins_signal[i].component[sSample].esyst=0;
         }
      }
      double integr,integrErr;
      integr=hVg.IntegralAndError(1,hVg.GetNbinsX()+1,integrErr);
      bkgIntegral.component["Vg"].count=integr;
      bkgIntegral.component["Vg"].estat=integrErr;
      bkgIntegral.component["Vg"].esyst=integr*cfg.sf.e_Vg/cfg.sf.Vg;

      integr=hGJ.IntegralAndError(1,hGJ.GetNbinsX()+1,integrErr);
      bkgIntegral.component["GJ"].count=integr;
      bkgIntegral.component["GJ"].estat=integrErr;
      bkgIntegral.component["GJ"].esyst=integr*cfg.sf.e_GJ/cfg.sf.GJ;

      integr=hData.IntegralAndError(1,hData.GetNbinsX()+1,integrErr);
      bkgIntegral.component["data"].count=integr;
      bkgIntegral.component["data"].estat=bkgIntegral.component["data"].esyst=0;

      for (TString sSample:{"efake","TTcomb","diboson"}){
         integr=fixHists[sSample].IntegralAndError(1,fixHists[sSample].GetNbinsX()+1,integrErr);
         bkgIntegral.component[sSample].count=integr;
         bkgIntegral.component[sSample].estat=integrErr;
         bkgIntegral.component[sSample].esyst=integr*cfg.datasets.getSystUncert(sSample);
      }
      for (TString sSample:{"TChiWG","T5Wg"}) {
         integr=fixHists[sSample].IntegralAndError(1,fixHists[sSample].GetNbinsX()+1,integrErr);
         signalIntegral.component[sSample].count=integr;
         signalIntegral.component[sSample].estat=integrErr;
         signalIntegral.component[sSample].esyst=0;
         debug << sSample << signalIntegral.component[sSample].count;
      }
   }

   gfx::SplitCan spcan;
   TCanvas can;
   if (plotMode==PRE) {
      can.cd();
      can.SetLogy();
   } else {
      spcan.cdUp();
      spcan.pU_.SetLogy();
   }
   hStackSum.Draw("axis");

   if (plotMode==SR && sVar=="METS") hStackSum.GetYaxis()->SetRangeUser(0.2,300);
   if (plotMode==SR && sVar=="STg")  hStackSum.GetYaxis()->SetRangeUser(0.4,500);
   if (plotMode==CR && sVar=="METS") hStackSum.GetYaxis()->SetRangeUser(2e-2,200);
   if (plotMode==PRE && sVar=="METS_logx") can.SetLogx();
   st.Draw("same");
   st.GetYaxis()->SetTitle(hStatErr.GetYaxis()->GetTitle());
   st.GetXaxis()->SetTitle(hStatErr.GetXaxis()->GetTitle());
   if (plotMode==CR) logFile<<"-- "+saveName;

   for (TString sSample:{"TChiWG","T5Wg"}) {
      TH1F &h=fixHists[sSample+"_stacked"];
      if (showSignal) {
         if (plotMode==PRE) {
            TH1F &hr=fixHists[sSample];
            hr.SetLineColor(hr.GetFillColor());
            hr.SetLineWidth(2);
            hr.SetFillStyle(0);
            hr.Draw("hist same");
            hr.Draw("hist same");
            le.append(hr,sSample,"l");
         } else {
            h.Draw("hist same");
            le.append(h,sSample,"l");
         }
      }
      if (plotMode==CR) {
         float nS=fixHists[sSample].Integral(0,h.GetNbinsX()+1,"width");
         float nSB=h.Integral(0,h.GetNbinsX()+1,"width");
         logFile<<sSample+TString::Format(": s/s+b = %f/%f = %f",nS,nSB,nS/nSB);
      }
   }
   // hStatErr.SetFillStyle(3354);
   // hStatErr.SetFillColor(kBlack);
   // hStatErr.Draw("same e2");
    hSystErr.SetFillStyle(3354);
    hSystErr.SetFillColor(kGray);
    hSystErr.SetLineColor(kWhite);
    hSystErr.Draw("same e2");
    le.append(hSystErr,"#sigma_{syst}","f");
   if (showData) hData.Draw("same pe1");
   spcan.pU_.RedrawAxis();
   can.RedrawAxis();
   le.buildLegend(.55,.71,-1,-1,2).DrawClone();
   if (plotMode==CR) {
      TString txt="CR";
      if (sSelection.Contains("/0l")) txt+=", 0 leptons";
      if (sSelection.Contains("/1l")) txt+=", 1 lepton";
      if (sSelection.Contains("/0b")) txt+=", 0 b";
      gfx::cornerLabel(txt,1).DrawClone();
   }else if (plotMode==VR) {
      TString txt="VR";
      gfx::cornerLabel(txt,1).DrawClone();
   }
   if (plotMode==PRE) {
      saver.save(can,saveName,!showData);
   } else {
      spcan.cdLow();
      TGraphErrors grRatioStat=hist::getRatioGraph(hStatErr,st,"ratio",hist::ONLY1);
      TH1F hRatioSyst=hist::getRatio(hSystErr,st,"ratio",hist::ONLY1);
      TH1F hRatio;
      if (showData) {
         hRatio=hist::getRatio(hData,st,"ratio",hist::ONLY1);
         hRatio.Draw("axis e0");
         hRatio.SetMaximum(1.9);
         hRatio.SetMinimum(0.1);
      } else {
         hRatioSyst.Draw("axis e0");
         hRatioSyst.SetMaximum(1.9);
         hRatioSyst.SetMinimum(0.1);
      }
      hRatioSyst.SetFillStyle(1001);
      hRatioSyst.SetFillColor(kGray);
      hRatioSyst.SetLineColor(kGray);
      hRatioSyst.Draw("same e2");
      grRatioStat.SetFillStyle(1001);
      grRatioStat.SetFillColor(kGray+1);
      grRatioStat.SetLineColor(kGray+1);
      gfx::scaleXerrors(grRatioStat,.3);
      grRatioStat.Draw("2");
      hRatio.Draw("same e0 axig");
      if (showData) hRatio.Draw("same pe0");
      spcan.pL_.RedrawAxis();
      le.clear();
      le.append(grRatioStat,"#sigma_{tot}","f");
      le.append(hRatioSyst,"#sigma_{syst}","f");
      le.buildLegend(.4,.38,-1,.55,2).DrawClone();
      saver.save(spcan,saveName,!showData);
   }

   if (plotMode==SR && sVar=="STg" && !sSelection.Contains("STg600")) {
      bool tex=false;
      for (unsigned iBin=0; iBin<signalBins.size(); iBin++) {
         logFile*sVar/signalBins[iBin].lower>>signalBins[iBin].upper;
         logFile<<signalBins[iBin].formatBlock(tex,true);
         logFile<<signalBins_signal[iBin].formatBlock(tex);
      }
      logFile*sVar>>"integral";
      logFile<<bkgIntegral.formatBlock(tex,true);
      logFile<<signalIntegral.formatBlock(tex);

      io::Logger texSummary((sVar+"-summary.tex").Data());
      texSummary<<"\\begin{tabular}{LRRRR}";
      texSummary<<"\\ST [\\text{GeV}] &  \\text{prediction} & \\sigma_\\text{stat} &  \\sigma_\\text{syst} & \\text{data}\\\\\\hline";
      tex=true;
      unsigned nBins=signalBins.size();
      for (unsigned iBin=1; iBin<nBins-1; iBin++) {
         float l=signalBins[iBin].lower;
         float u=signalBins[iBin].upper;
         // actually contains overflow:
         if (iBin==nBins-2) u=std::numeric_limits<float>::infinity();

         io::Logger texFile((sVar+TString::Format("-%.0f-%.0f.tex",l,u)).Data());

         texFile<<"\\begin{tabular}{lrrrr}";
         TString binString=sVar+" bin: ";
         binString+=TString::Format("$%.0f-%.0f\\GeV$",l,u);
         binString.ReplaceAll("inf","\\infty");
         binString.ReplaceAll("STg","\\ST{}");
         texFile<<"   \\multicolumn{5}{c}{"+binString+"}\\\\";
         texFile<<"               &  yield & $\\sigma_\\text{stat}$ &  $\\sigma_\\text{syst}$ & $(/$tot. bkg.)\\\\\\hline";
         texFile<<signalBins[iBin].formatBlock(tex,true);
         texFile*signalBins_signal[iBin].formatBlock(tex);
         texFile<<"\\end{tabular}";

         binString=TString::Format("%.0f-%.0f     ",l,u);
         binString.ReplaceAll("inf","\\infty");
         binString.ReplaceAll("STg","\\ST{}");
         texSummary*TString::Format("%-16s &",binString.Data());
         SignalBin::Component tot=signalBins[iBin].getTotal();
    //     float esy=TMath::Sqrt(tot.esyst*tot.esyst+.13*.13*signalBins[iBin].component["GJ"].count*signalBins[iBin].component["GJ"].count);
         float esy=TMath::Sqrt(tot.esyst*tot.esyst);
         texSummary*TString::Format(" %6.1f & \\pm%6.1f & \\pm%6.1f &",tot.count,tot.estat,esy);
         texSummary<<TString::Format("%6.0f \\\\",signalBins[iBin].component["data"].count);
      }
      texSummary<<"\\end{tabular}";
      io::Logger texFile((sVar+"-integral.tex").Data());
      texFile<<"\\begin{tabular}{lrrrr}";
      texFile<<"               &  yield & $\\sigma_\\text{stat}$ &  $\\sigma_\\text{syst}$ & $(/$tot. bkg.)\\\\\\hline";
      texFile<<bkgIntegral.formatBlock(tex,true);
      texFile*signalIntegral.formatBlock(tex);
      texFile<<"\\end{tabular}";

      static io::RootFileSaver saver_yields("yields.root","");
      for (auto comp: signalBins[0].component) {
         TString sBkg=comp.first;
         TH1F h("","",nBins-2,0,1);
         TH1F hSys("","",nBins-2,0,1);
         for (unsigned iBin=1; iBin<nBins-1; iBin++) {
            h.SetBinContent(iBin,signalBins[iBin].component[sBkg].count);
            h.SetBinError(iBin,signalBins[iBin].component[sBkg].estat);
            hSys.SetBinContent(iBin,signalBins[iBin].component[sBkg].esyst);
         }
         saver_yields.save(h,sSelection+sVar+"/"+sBkg);
         saver_yields.save(hSys,sSelection+sVar+"/"+sBkg+"_esyst");
      }
   }
}

// plot functions to call
void plot(TString sSelection,TString sVar,int iRebin,PlotMode_t plotMode)
{
   return plot(sSelection,sVar,iRebin,{},{},plotMode);
}
void plot(TString sSelection,TString sVar,
          std::vector<float> rebinEdges,std::vector<float> rebinWidths,
          PlotMode_t plotMode
   )
{
   return plot(sSelection,sVar,0,rebinEdges,rebinWidths,plotMode);
}

extern "C"
void run()
{
   
   plot("pre_ph165/","METdotPh",{-100000,100000},{4000},CR);
   plot("pre_ph165/","METdotJet",{-100000,100000},{4000},CR);
   plot("pre_ph165/","absphiMETHT",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/","absphiPhHT",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/","DeltaPhiMETjet100",{0,3.2},{.1},CR);
   plot("pre_ph165/","vecHTPt",{0,300,600,1000,1010},{20,30,50,10},CR);
   
   plot("pre_ph165/MT100/","METdotPh",{-100000,100000},{4000},CR);
 //  plot("pre_ph165/MT100/","ph1Pt",{0,300,600,1400,1410},{20,30,50,10},CR);    
   //2016 studies
   plot("pre_ph165/dPhiMETphl06/","ph1Pt",{0,300,600,1400,1410},{20,30,50,10},CR);
   plot("pre_ph165/dPhiMETphl06/","MET",{0,300,600,1000,1410},{20,30,50,10},CR);
   plot("pre_ph165/dPhiMETphl06/","sieie"      ,{0,0.03},{0.0001},CR);
 
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
 //  plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","ph1Pt",{0,500,750,1000,1300},{500,250,250,300},CR);
   
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","MET",{80,100,250,350},{20,30,100},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETph",{0,1.8,2.6,3.3},{0.1,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh",{0,.8,3.2},{.2,.4},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETjet",{0,1.8,2.6,3.3},{.3,.2,.1},CR);   
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_MET100/MT100/","r9",{0,1.1},{0.05},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","r9_phoPtl500",{0,1.1},{0.05},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","r9_phoPtl700",{0,1.1},{0.05},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","r9_phoPtg700",{0,1.1},{0.05},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_phoPtl500",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_phoPtl700",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETnJetPh_phoPtg700",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiMETHT",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","absphiPhHT",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","vecHTPt",{0,300,600,1000,1010},{20,30,50,10},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","DeltaPhiMETjet100",{0,3.2},{.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETphl02/","ph1Pt",{0,500,750,1000,1300},{500,250,250,300},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/","absphiMETnJetPh",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/","ph1Pt",{0,500,750,1000,1300},{500,250,250,300},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/dPhiMETPhg03/","STg"   ,{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   

 //  plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","STg",{600,800,1000,1300,1600},{200,200,300,300},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","Nj",1,CR);

   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/","absphiMETnJetPh",{0,1.8,2.6,3.3},{.3,.2,.1},CR);   
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/Njl3/","absphiMETjet",{0,1.8,2.6,3.3},{.3,.2,.1},CR);

   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/","absphiMETnJetPh",{0,1.8,2.6,3.3},{.3,.2,.1},CR);   
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/0b/","absphiMETjet",{0,1.8,2.6,3.3},{.3,.2,.1},CR);

   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/","absphiMETnJetPh",{0,1.8,2.6,3.3},{.3,.2,.1},CR);   
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/0l/","absphiMETjet",{0,1.8,2.6,3.3},{.3,.2,.1},CR);

   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","absphiMETnJetPh",{0,1.8,2.6,3.3},{.3,.2,.1},CR);   
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","absphiMETjet",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","relPt2Jets",{0,1},{0.04},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","DeltaS",{0,3.2},{0.2},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","DeltaS1",{0,3.2},{0.2},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","MT",{100,300,500,600},{40,20,100},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","ph1Pt",{0,300,600,1000,1010},{20,30,50,10},CR);
   plot("pre_ph165/c_MET100/MT100/METl300vMTl300/1l/","MET",{80,100,250,350},{20,30,100},CR);
       
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","MET",{80,100,250,350},{20,30,100},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_MET100/MT100/METl400vMTl400/","Nj",1,CR);
   
   plot("pre_ph165/","MT",{0,400,800,2000},{40,80,200},CR);  
   plot("pre_ph165/","ph1Pt",{0,400,800,1200},{40,80,200},CR);
      
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","MET",{80,100,250,350},{20,30,100},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_MET150/MT100/METl300vMTl300/","Nj",1,CR);
       
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","ph1Pt",{0,300,500,600,1000},{20,50,100,200},CR);
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","MET",{80,100,250,350},{20,30,100},CR);
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);     
   plot("pre_ph165/c_MET150/MT150/METl400vMTl400/","Nj",1,CR);   
   
   // Control region
   /*
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","MET",{80,100,220,400,500},{20,30,180,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","METS",{20,80,85},{5,5},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","sieie1",{0.0,0.008,0.011},{0.0005,0.0002},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","sipip1",{0.0,0.009,0.011,0.015,0.03},{0.001,0.0002,0.0005,0.005},CR);

   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","Nph",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","Nj",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","Nele",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","Nmu",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/","Nl",1,CR);

   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/","sieie1",{0.0,0.008,0.011},{0.0005,0.0002},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/eta1/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},CR);

   // with b-veto
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","MET",{80,100,220,400,500},{20,30,180,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","METS",{20,80,85},{5,5},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","sieie1",{0.0,0.008,0.011},{0.0005,0.0002},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","sipip1",{0.0,0.009,0.011,0.015,0.03},{0.001,0.0002,0.0005,0.005},CR);

   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","Nph",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","Nj",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","Nele",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","Nmu",1,CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0bT/","Nl",1,CR);

   // Split by number of leptons
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","MET",{80,100,220,400,500},{20,30,180,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","METS",{20,80,85},{5,5},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","absphiMETjet",{0,1.8,2.6,3.3},{.3,.2,.1},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/0l/","Nph",1,CR);

   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","ph1Pt",{0,300,500,600,1000,1100},{20,50,100,200,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","MET",{80,100,220,400,500},{20,30,180,100},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","METS",{20,80,85},{5,5},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","HT",{0,300,500,1200,1550},{100,200,350,350},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","STg",{260,330,450,750,1000,1250},{90,120,150,250,250},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","MT",{100,250,400,600,700,900},{150,20,50,100,200},CR);
   plot("pre_ph165/c_S30/MT100/Sl80vMTl300/1l/","Nph",1,CR);


   // Validation region
   plot("pre_ph165/c_S80/MT300/STgl600/","ph1Pt",{140,300,400},{20,50},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","MET",{100,280,380},{30,50},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","METS",{80,200,500},{40,50},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","HT",{0,300,500,1200,1550},{100,200,700,350},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","STg",{260,620},{20},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","MT",{100,250,400,600,700,900},{150,20,50,100,200},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","sieie1",{0.0,0.008,0.011},{0.0005,0.0002},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","sipip1",{0.0,0.009,0.011,0.015,0.03},{0.001,0.0002,0.0005,0.005},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},VR);
   plot("pre_ph165/c_S80/MT300/STgl600/","Nph",1,VR);
   */

   //Validation region 2016
   plot("pre_ph165/c_MET300/MT300/STgl600/","STg",{260,620},{20},VR);   
   plot("pre_ph165/c_MET300/MT300/STgl600/","ph1Pt",{140,300,400},{20,50},VR);
   plot("pre_ph165/c_MET300/MT300/STgl600/","MET",{100,280,380},{30,50},VR);
   plot("pre_ph165/c_MET300/MT300/STgl600/","MT",{100,280,400,600,700,900},{180,120,50,100,200},VR);  
   plot("pre_ph165/c_MET300/MT300/STgl600/","absphiMETjet",{0,.4,3.2,3.4},{.2,.4,.1},VR);
   plot("pre_ph165/c_MET300/MT300/STgl600/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},VR);
   plot("pre_ph165/c_MET300/MT300/STgl600/","absphiMETnJetPh",{0,1.8,2.6,3.3},{.3,.2,.1},VR);
   plot("pre_ph165/c_MET300/MT300/STgl600/","METdotPh",{-100000,100000},{10000},VR);
   
   // preselection
//   plot("pre_ph165/","METS_logx",1,PRE);

   // Signal region
   /*
   plot("pre_ph165/c_S80/MT300/","ph1Pt",5,SR_BLIND);
   plot("pre_ph165/c_S80/MT300/","MET",5,SR_BLIND);
   plot("pre_ph165/c_S80/MT300/","METS",5,SR_BLIND);
   plot("pre_ph165/c_S80/MT300/","HT",5,SR_BLIND);
   plot("pre_ph165/c_S80/MT300/","METSHT",5,SR_BLIND);
   plot("pre_ph165/c_S80/MT300/","STg",5,SR_BLIND);*/
   plot("pre_ph165/c_MET300/MT300/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR);
   plot("pre_ph165/c_MET300/MT300/","STg",{480,600,800,1000,1300,1600},{120,200,200,300,300},SR);  
   plot("pre_ph165/c_MET300/MT300/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},SR);
   plot("pre_ph165/c_MET300/MT300/","absphiMETHT",{0,1.8,2.6,3.3},{.3,.2,.1},SR);
   plot("pre_ph165/c_MET300/MT300/","absphiPhHT",{0,1.8,2.6,3.3},{.3,.2,.1},SR);
   plot("pre_ph165/c_MET300/MT300/","vecHTPt",{0,300,600,1000,1010},{20,30,50,10},SR);
   plot("pre_ph165/c_MET300/MT300/","METdotPh",{-100000,-90000,100000},{5000,19000},SR);
   plot("pre_ph165/c_MET300/MT300/","METSHT",{0,4,12,20,40,100},{2,1,2,20,60},SR);
   plot("pre_ph165/c_MET300/MT300/","ph1Pt",{0,1500},{100},SR);
   plot("pre_ph165/c_MET300/MT300/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR);
   
   // halo/spike check
   /*
   plot("pre_ph165/c_S80/MT300/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},SR_CHECK);
   plot("pre_ph165/c_S80/MT300/","sieie1",{0.0,0.008,0.011},{0.0005,0.0002},SR_CHECK);
   plot("pre_ph165/c_S80/MT300/","sipip1",{0.0,0.009,0.011,0.015,0.03},{0.001,0.0002,0.0005,0.005},CR);
   plot("pre_ph165/c_S80/MT300/","Nph",1,SR_CHECK);
   plot("pre_ph165/c_S80/MT300/eta1/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},SR_CHECK);
   plot("pre_ph165/c_S80/MT300/eta1/","sieie1",{0.0,0.008,0.011},{0.0005,0.0002},SR_CHECK);
   plot("pre_ph165/c_S80/MT300/eta1/","Nph",1,SR_CHECK);
   

   // Interpretation binning
   plot("pre_ph165/c_S80/MT300/","METS",{80,200,400,750,800},{120,200,350,50},SR_BLIND);
   plot("pre_ph165/c_S80/MT300/","STg",{260,600,800,1000,1200},{340,200,200,200},SR_BLIND);
   // unblinded
   plot("pre_ph165/c_S80/MT300/","METS",{80,200,400,750,800},{120,200,350,50},SR);
   plot("pre_ph165/c_S80/MT300/","STg",{260,600,800,1000,1200},{340,200,200,200},SR);

   // excluding VR
   plot("pre_ph165/c_S80/MT300/STg600/","STg",{600,800,1000,1200},{200,200,200},SR);
*/
   // further SR distributions with data

   //2016
   plot("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR);  
   plot("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","ph1Pt",{0,700,1500},{100,200},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","MET",{0,700,1500},{100,200},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","absphiMETph",{0,1.8,2.6,3.3},{.3,.2,.1},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","r9",{0,1.1},{0.05},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","r9_Stgg1300",{0,1.1},{0.05},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","METdotPh",{-100000,-90000,100000},{5000,19000},SR);
   plot("pre_ph165/c_MET300/MT300/STg600/","DeltaPhiMETjet100",{0,3.2},{.1},SR);
  
   
   //2015
   /*
   plot("pre_ph165/c_S80/MT300/","ph1Pt",{0,400,700,1000},{100,150,300},SR);
   plot("pre_ph165/c_S80/MT300/","MET"  ,{100,400,700,1000},{100,150,300},SR);
   plot("pre_ph165/c_S80/MT300/","MT"   ,{0,600,900,1500},{100,150,600},SR);
   plot("pre_ph165/c_S80/MT300/","HT"   ,{0,300,600,1000,1500},{100,150,400,500},SR);
    */
}
