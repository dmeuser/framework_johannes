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
static io::RootFileSaver saver("plots.root",TString::Format("danilo_SFapplication%.1f",cfg.processFraction*100));
static io::Logger logFile("SFapplication.log");


float Vg_JES_bin1 = 0.059295;
float Vg_JES_bin2 = 0.055648;
float Vg_JES_bin3 = 0.053110;
float Vg_JES_bin4 = 0.050381;

float GJ_JES_bin1 = 0.009304;
float GJ_JES_bin2 = 0.070263;
float GJ_JES_bin3 = 0.070263;
float GJ_JES_bin4 = 0.323459;

float Vg_mu_bin1 = 0.038008;
float Vg_mu_bin2 = 0.054549;
float Vg_mu_bin3 = 0.065654;
float Vg_mu_bin4 = 0.090400;

float Vg_pdf_bin1 = 0.0159707;
float Vg_pdf_bin2 = 0.0200078;
float Vg_pdf_bin3 = 0.0244478;
float Vg_pdf_bin4 = 0.0382614;

float GJ_mu_bin1 = 0.051947;
float GJ_mu_bin2 = 0.027569;
float GJ_mu_bin3 = 0.054835;
float GJ_mu_bin4 = 0.071058;

float GJ_pdf_bin1 = 0.0815288;
float GJ_pdf_bin2 = 0.0255214;
float GJ_pdf_bin3 = 0.0192825;
float GJ_pdf_bin4 = 0.0512977;

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
   io::RootFileReader histReader(TString::Format("danilo_histograms_%s.root",cfg.treeVersion.Data()),TString::Format("danilo_distributions%.1f",cfg.processFraction*100));

   TString saveName=sSelection+sVar;

   bool divideByBinWidth=(plotMode==!CR || plotMode==!VR || plotMode==PRE);
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
   le.prepend(fixHists["TTcomb"],"t#bar{t}(#gamma)","f");

   if (plotMode==CR) {
      // CR: more GJet
      st.Add(&hVg,"hist");
      st.Add(&hGJ,"hist");
      le.prepend(hVg,"V(#gamma)","f");
      le.prepend(hGJ,"#gamma+jets","f");
   } else {
      // SR,VR: more Vg
      st.Add(&hGJ,"hist");
      st.Add(&hVg,"hist");
      le.prepend(hGJ,"#gamma+jets","f");
      le.prepend(hVg,"V(#gamma)","f");
   }
   TH1F hStackSum(*(TH1F*)st.GetStack()->Last());
   for (TString sSample:{"TChiWG","T5Wg"}) {
      fixHists[sSample+"_stacked"]=fixHists[sSample];
      TH1F &h=fixHists[sSample+"_stacked"];
      h.Add(&hStackSum);
      h.SetLineColor(h.GetFillColor());
      h.SetLineWidth(3);
      h.SetFillStyle(0);
   }

   TH1F hData;
   if (showData){
      hData=rebinned(*histReader.read<TH1F>(sSelection+sVar+"/SinglePhoton"));
      le.prepend(hData,"Data","pe");
   }
   TH1F hStatErr(hStackSum);
   hStatErr.SetMarkerSize(0);
   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << hStatErr.GetBinError(1) << std::endl;
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
      std::cout << "yield and syst uncert Vgamma:  " << hVg.GetBinContent(i) << "$\\pm$" << TMath::Sqrt(eVg) << std::endl;
      std::cout << "yield and syst uncert gammaJets:  " << hGJ.GetBinContent(i) << "$\\pm$" << TMath::Sqrt(eGJ) << std::endl;
      //std::cout << "syst error squared:  " << erri << "   correlation error squared:  " << error_corr << std::endl;
     // std::cout << "correlation error fraction:  " << error_corr/erri << std::endl;
      erri+=error_corr;
      hSystErr.SetBinError(i,TMath::Sqrt(erri));
   }

   // syst. uncert. from fixed backgrounds
   for (TString sSample:{"efake","TTcomb","diboson"}){
      TH1F const &h=fixHists[sSample];
      float const uncert=cfg.datasets.getSystUncert(sSample);
      for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
         erri=h.GetBinContent(i)*uncert;
         std::cout << "bin number:  " << i;
         if (sSample == "efake"){
            std::cout << "  yield and syst uncert efake:  " << h.GetBinContent(i) << "$\\pm$" << erri << std::endl;
            }
         else if (sSample == "TTcomb"){
            std::cout << "  yield and syst uncert TTcomb:  " << h.GetBinContent(i) << "$\\pm$" << erri << std::endl;
            }
         else if (sSample == "diboson"){
            std::cout << "  yield and syst uncert diboson:  " << h.GetBinContent(i) << "$\\pm$" << erri << std::endl;
            }
         erri=util::sqrtQuadSum<double>({erri,hSystErr.GetBinError(i)});
         hSystErr.SetBinError(i,erri);
         std::cout << "bin syst error total: " << erri << std::endl;
         std::cout << "bin stat error total: " << hStatErr.GetBinError(i) << std::endl;         
         comb_uncert = util::sqrtQuadSum<double>({hStatErr.GetBinError(i),hSystErr.GetBinError(i)});
         //hStatErr is now total error!
       //  hStatErr.SetBinError(i,comb_uncert);
         std::cout << "bin total error total: " << comb_uncert << std::endl;
         
      }
   }

   for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
      std::cout << "!!!!!!!!!BIN!!!!!!!!!!!" << i << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << hStatErr.GetBinError(i) << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << hSystErr.GetBinError(i) << std::endl;
      hStatErr.SetBinError(i, util::sqrtQuadSum<double>({hStatErr.GetBinError(i),hSystErr.GetBinError(i)}));
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << hStatErr.GetBinError(i) << std::endl;
   }

   // store separate background uncertainties
   if (plotMode==SR  && sVar=="STg") {
      for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
         std::map<TString,SignalBin::Component> &bini=signalBins[i].component;
         bini["Vg"].count=hVg.GetBinContent(i);
         std::cout << "bin:  " << i << "  Vgamma yield and stat uncert:  " << hVg.GetBinContent(i) << "$\\pm$" << hVg.GetBinError(i)<< std::endl;
         bini["Vg"].estat=hVg.GetBinError(i);
         bini["Vg"].esyst=hVg.GetBinContent(i)*cfg.sf.e_Vg/cfg.sf.Vg;

         bini["GJ"].count=hGJ.GetBinContent(i);
         std::cout << "bin:  " << i << "  gammaJets yield and stat uncert:  " << hGJ.GetBinContent(i) << "$\\pm$" << hGJ.GetBinError(i)<< std::endl;
         bini["GJ"].estat=hGJ.GetBinError(i);
         bini["GJ"].esyst=hGJ.GetBinContent(i)*cfg.sf.e_GJ/cfg.sf.GJ;

         bini["data"].count=hData.GetBinContent(i);
         bini["data"].estat=bini["data"].esyst=0;

         for (TString sSample:{"efake","TTcomb","diboson"}){
            bini[sSample].count=fixHists[sSample].GetBinContent(i);
            bini[sSample].estat=fixHists[sSample].GetBinError(i);
            if (sSample == "efake"){
               std::cout << "bin:  " << i << "  efake yield and stat uncert:  " << fixHists[sSample].GetBinContent(i) << "$\\pm$" << fixHists[sSample].GetBinError(i) << std::endl;
               }
            else if (sSample == "TTcomb"){
               std::cout << "bin:  " << i << "  TTcomb yield and stat uncert:  " << fixHists[sSample].GetBinContent(i) << "$\\pm$" << fixHists[sSample].GetBinError(i)<< std::endl;
               }
            else if (sSample == "diboson"){
               std::cout << "bin:  " << i << "  diboson yield and stat uncert:  " << fixHists[sSample].GetBinContent(i) << "$\\pm$" << fixHists[sSample].GetBinError(i)<< std::endl;
               }
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
   hStackSum.SetStats(false);
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
            hr.SetLineWidth(3);
            hr.SetFillStyle(0);
            hr.Draw("hist same");
            hr.Draw("hist same");
            le.append(hr,sSample,"l");
         } else {
            if (sSample == "T5Wg"){
               h.SetLineColor(kGreen+1);
            }
            h.Draw("hist same");
            h.SetLineWidth(3);
            le.append(h,sSample,"l");
         }
      }
      if (plotMode==CR) {
         float nS=fixHists[sSample].Integral(0,h.GetNbinsX()+1,"width");
         float nSB=h.Integral(0,h.GetNbinsX()+1,"width");
         logFile<<sSample+TString::Format(": s/s+b = %f/%f = %f",nS,nSB,nS/nSB);
      }
   }
   //~ hStatErr.SetFillStyle(3354);
   //~ hStatErr.SetFillColor(kBlack);
   //~ hStatErr.Draw("same e2");
   hSystErr.SetFillStyle(3354);
   hSystErr.SetFillColor(kGray);
   hSystErr.SetLineColor(kWhite);
   hSystErr.Draw("same 0e2");

   std::cout << bkgIntegral.formatBlock(false,true) << std::endl;

    le.append(hSystErr,"#sigma_{syst}","f");
   if (showData) hData.Draw("same pe1");
   spcan.pU_.RedrawAxis();
   can.RedrawAxis();
   //le.buildLegend(.55,.71,-1,-1,2).DrawClone();
   le.buildLegend(.7,.75,-1,-1,2).DrawClone();
   if (plotMode==CR) {
      TString txt=""; //add "CR" in plot if wished
      if (sSelection.Contains("/0l")) txt+=", 0 leptons";
      if (sSelection.Contains("/1l")) txt+=", 1 lepton";
      if (sSelection.Contains("/0b")) txt+=", 0 b";
      gfx::cornerLabel(txt,1).DrawClone();
   }else if (plotMode==VR) {
      TString txt=""; // here set "VR" if you want to indicate the region
      gfx::cornerLabel(txt,1).DrawClone();
   }
   if (plotMode==PRE) {
      saver.save(can,saveName,!showData);
   } else {
      spcan.cdLow();
      TGraphErrors grRatioStat=hist::getRatioGraph(hStatErr,st,"Ratio",hist::ONLY1);
      TH1F hRatioSyst=hist::getRatio(hSystErr,st,"Ratio",hist::ONLY1);
      TH1F hRatio;
      if (showData) {
         hRatio=hist::getRatio(hData,st,"Data/Pred.",hist::ONLY1);
         hRatio.SetStats(false);
         hRatio.Draw("axis e0");
         hRatio.SetMaximum(1.9);
         hRatio.SetMinimum(0.1);
      } else {
         hRatioSyst.SetStats(false);
         hRatioSyst.Draw("axis e0");
         hRatioSyst.SetMaximum(1.9);
         hRatioSyst.SetMinimum(0.1);
      }
      hRatioSyst.SetFillStyle(1001);
      hRatioSyst.SetFillColor(kGray);
      hRatioSyst.SetLineColor(kGray);
      hRatioSyst.Draw("same 0e2");
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
      le.buildLegend(.4,.43,-1,.6,2).DrawClone();
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
      for (int i = 0; i < hSystErr.GetNbinsX(); i++){
         std::cout << signalBins[i].formatBlock(false,true) << std::endl;
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
   plot("pre_ph165/combined/","HTG",{0,2400},{240},VR);
   plot("pre_ph165/combined/","MET",{300,800},{50},VR);
   plot("pre_ph165/combined/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   
   plot("pre_ph165/VR_SR/","HTG",{0,2400},{240},VR);
   plot("pre_ph165/VR_SR/","MET",{300,800},{50},VR);
   plot("pre_ph165/VR_SR/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   plot("pre_ph165/VR_SR/","phoPt",{200,1000},{80},VR);
   plot("pre_ph165/VR_SR/","MT",{300,1000},{70},VR);
   plot("pre_ph165/VR_SR/","phoEta",{-2.6,2.6},{0.57},VR);
   plot("pre_ph165/VR_SR/","nEle",{0,5},{1},VR);
   plot("pre_ph165/VR_SR/","nMu",{0,5},{1},VR);
   plot("pre_ph165/VR_SR/","lepPt",{0,800},{80},VR);
   plot("pre_ph165/VR_SR/","nPho",{0,5},{1},VR);
   
   plot("pre_ph165/VR_SR/noHTG/","HTG",{0,2400},{240},VR);
   plot("pre_ph165/VR_SR/noHTG/","MET",{300,800},{50},VR);
   plot("pre_ph165/VR_SR/noHTG/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   plot("pre_ph165/VR_SR/noHTG/","phoPt",{200,1000},{100},VR);
   plot("pre_ph165/VR_SR/noHTG/","MT",{300,1000},{100},VR);
   plot("pre_ph165/VR_SR/noHTG/","phoEta",{-2.6,2.6},{0.57},VR);
   
   plot("pre_ph165/VR_SR/noLepton/","HTG",{0,2400},{240},VR);
   plot("pre_ph165/VR_SR/noLepton/","MET",{300,800},{50},VR);
   plot("pre_ph165/VR_SR/noLepton/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   plot("pre_ph165/VR_SR/noLepton/","phoPt",{200,1000},{80},VR);
   plot("pre_ph165/VR_SR/noLepton/","MT",{300,1000},{70},VR);
   plot("pre_ph165/VR_SR/noLepton/","phoEta",{-2.6,2.6},{0.57},VR);
   
   plot("pre_ph165/VR_SR/noDiphoton/","HTG",{0,2400},{240},VR);
   plot("pre_ph165/VR_SR/noDiphoton/","MET",{300,800},{50},VR);
   plot("pre_ph165/VR_SR/noDiphoton/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   plot("pre_ph165/VR_SR/noDiphoton/","phoPt",{200,1000},{80},VR);
   plot("pre_ph165/VR_SR/noDiphoton/","MT",{300,1000},{70},VR);
   plot("pre_ph165/VR_SR/noDiphoton/","phoEta",{-2.6,2.6},{0.57},VR);
   
   plot("pre_ph165/VR_SR/exclusiv/","HTG",{0,2400},{240},VR);
   plot("pre_ph165/VR_SR/exclusiv/","MET",{300,800},{50},VR);
   plot("pre_ph165/VR_SR/exclusiv/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   plot("pre_ph165/VR_SR/exclusiv/","phoPt",{200,1000},{100},VR);
   plot("pre_ph165/VR_SR/exclusiv/","MT",{300,1000},{100},VR);
   plot("pre_ph165/VR_SR/exclusiv/","phoEta",{-2.6,2.6},{0.57},VR);
   
   plot("pre_ph165/VR/exclusiv/","HTG",{0,1200},{120},VR);
   plot("pre_ph165/VR/exclusiv/","MET",{200,500},{50},VR);
   plot("pre_ph165/VR/exclusiv/","absdPhi_pmMet_Pho",{0,1.6},{0.16},VR);
   plot("pre_ph165/VR/exclusiv/","phoPt",{150,400},{22},VR);
   plot("pre_ph165/VR/exclusiv/","MT",{200,700},{50},VR);
   plot("pre_ph165/VR/exclusiv/","phoEta",{-2.6,2.6},{0.57},VR);
   plot("pre_ph165/VR/exclusiv/","STG",{400,620},{20},VR);
}
