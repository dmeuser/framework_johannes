#ifndef HIST_HPP__
#define HIST_HPP__

#include "Config.hpp"
#include "Dataset.hpp"
#include "tools/gfx.hpp"

#include <map>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>

namespace hist
{
   enum ErrorType
   {
      STAT,
      SYST,
      COMB,
      // special for ratios, residuals, ...
      ONLY1,
      ONLY2
   };

   template <class HIST>
   class Histograms
   {
   public:
      Histograms(std::vector<TString> const &samples,DatasetCollection const &datasets=Config::get().datasets);
      Histograms(std::vector<std::string> const &samples,DatasetCollection const &datasets=Config::get().datasets);

      // 1d
      void addHist(TString const &varName,TString const &title, Int_t nbinsx, Double_t xlow, Double_t xup);
      void addHist(TString const &varName,TString const &title, std::vector<float> edges, std::vector<float> widths);
      // 2d
      void addHist(TString const &varName,TString const &title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);
      void addHist(TString const &varName,TString const &title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y);
      // simple counters
      void addCounter(TString const &varName);

      void setCurrentSample(TString const &current);
      void setFillWeight(float w);
      void fill(TString const &varName,float x); // 1d
      void fill(TString const &varName,float x,float y); // 2d   
      void fillweight(TString const &varName,float x,float w);
      void fillbin(TString const &varName,TString const &binName);
      void fillbinFake(TString const &varName,TString const &binName);
      void count(TString const &varName);
      void scaleLumi(); // scales MC with lumi weight and trigger efficiency. Data is ignored.
      void mergeOverflow(bool includeUnderflow=true); // add the overflow to the last bin (and underflow to first)
      void combineFromSubsamples(std::vector<TString> const &samples);
      std::vector<TString> getVariableNames();
      std::vector<HIST*> getHistograms(TString const &varName,std::vector<TString> const &samples,bool divideByBinWidth=false);
      HIST* getHistogram(TString const &varName,TString const &sample,bool divideByBinWidth=false);
      THStack getStack(TString const &varName,std::vector<TString> const& samples,bool divideByBinWidth=false);
      float getCount(TString const &varName, TString const &sample);
      float getCountError(TString const &varName, TString const &sample);

      gfx::LegendEntries getLegendEntries();
   private:
      static int s_iInstances; // count all instances
      int const iInstance_;    // store which instance this is

      DatasetCollection const &datasets;

      std::vector<TString> vsSamples_;
      std::string sCurrentSample_;
      std::map<TString,std::map<TString,HIST>> mmH_; // [var][sample]
      std::map<TString,std::map<TString,float>> mCount_; // [var][sample]
      std::map<TString,std::map<TString,float>> mCountError2_; // [var][sample]

      float fWeight_;

      gfx::LegendEntries le_;
   };
   template class Histograms<TH1F>;
   template class Histograms<TH2F>;

   // forbid:
   template <>
   void Histograms<TH2F>::addHist(TString const&,TString const&,Int_t,Double_t,Double_t)=delete;
   template <>
   void Histograms<TH2F>::addHist(TString const&,TString const&,std::vector<float>,std::vector<float>)=delete;
   template <>
   void Histograms<TH1F>::addHist(TString const&,TString const&,Int_t,Double_t,Double_t,Int_t,Double_t,Double_t)=delete;
   template <>
   void Histograms<TH1F>::addHist(TString const &varName,TString const &title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y)=delete;
   template <>
   void Histograms<TH2F>::fill(TString const&,float)=delete;
   template <>
   void Histograms<TH2F>::fillbin(TString const &varName,TString const &binName)=delete;
   template <>
   void Histograms<TH2F>::fillbinFake(TString const &varName,TString const &binName)=delete;
   template <>
   void Histograms<TH1F>::fill(TString const&,float,float)=delete;
   //~ template <>
   //~ void Histograms<TH2F>::mergeOverflow(bool)=delete;

   std::vector<double> getBinVector(std::vector<float> edges, std::vector<float> widths);
   TH1F fromWidths(const char *name, const char *title,std::vector<float> edges, std::vector<float> widths);
   TH2F fromWidths_2d(const char *name, const char *title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y);


   TH1F rebinned(TH1F const &h, std::vector<float> const &edges, std::vector<float> const &widths,bool mergeOverflow=true,bool mergeUnderflow=true);
   TH1F rebinned(TH1F const &h, std::vector<double> const &binedges,bool mergeOverflow=true,bool mergeUnderflow=true);
   void divideByBinWidth(TH1& h,bool divideLastBin=true);
   void mergeOverflow(TH1& h, bool includeUnderflow=true);
   void mergeOverflow(TH2& h, bool includeUnderflow=true);

   void setMaximum(TH1& h,std::vector<TH1F> hists,float multiplier=1.1);
   void setMinimum(TH1& h,std::vector<TH1F> hists,float multiplier=0.9,bool allowNegative=true);

   TH1F getRatio   (TH1F const &h1,TH1F const &h2,TString title="ratio",ErrorType et=COMB);
   TH1F getResidual(TH1F const &h1,TH1F const &h2,TString title="residual",ErrorType et=COMB);
   TH1F getPull    (TH1F const &h1,TH1F const &h2,TString title="pull",ErrorType et=COMB);
   TH1F getRatio   (TH1F const &h1,THStack &h2,TString title="ratio",ErrorType et=COMB);
   TH1F getPull    (TH1F const &h1,THStack &h2,TString title="pull",ErrorType et=COMB);

   TGraphErrors getRatioGraph(TH1F const &h1,TH1F const &h2,TString title="ratio",ErrorType et=COMB);
   TGraphErrors getRatioGraph(TH1F const &h1,THStack &h2,TString title="ratio",ErrorType et=COMB);

   THStack stackPrepend(THStack const& stOld, TH1F &h, Option_t *option="");
}

#endif /* HIST_HPP__ */
