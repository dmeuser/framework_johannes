#ifndef EFFICIENCY_HPP
#define EFFICIENCY_HPP

#include "tools/io.hpp"
#include "tools/hist.hpp"

#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

namespace eff
{
   class Efficiency
   {
   public:
      Efficiency(){}
      Efficiency(TString name,TString var,Int_t nbins, Double_t xlow, Double_t xup, std::vector<float> const plateauCut={});
      Efficiency(TString name,TString var,std::vector<float> edges,std::vector<float> widths,std::vector<float> const plateauCut={});
      void fill(Bool_t bPassed, Double_t x);
      void save(io::RootFileSaver const &saver,TString datasetName,bool sim);
   private:
      TEfficiency eff_;
      TString name_;
      std::vector<float> plateauCut_;
      bool varBins_;
   };

   class Efficiencies
   {
   public:
      Efficiencies(){}
      void add(TString name,TString var,Int_t nbins,Double_t xlow,Double_t xup,std::vector<float> const plateauCut={});
      void add(TString name,TString var,std::vector<float> edges,std::vector<float> widths,std::vector<float> const plateauCut={});
      void fill(TString name,Bool_t bPassed, Double_t x);
      void saveAll(io::RootFileSaver const &saver,TString datasetName,bool sim);
   private:
      std::map<TString,Efficiency> effs_;
   };
}
#endif /* EFFICIENCY_HPP */
