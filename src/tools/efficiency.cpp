#include "efficiency.hpp"


void eff::Efficiency::save(io::RootFileSaver const &saver,TString datasetName,bool sim)
{
   TCanvas can;
   // histograms
   TH1 *htot=eff_.GetCopyTotalHisto();
   htot->SetLineColor(kGray+1);
   htot->SetFillColor(kGray+1);
   htot->SetFillStyle(1001);
   TH1 *hpass=eff_.GetCopyPassedHisto();
   hpass->SetLineColor(kBlack);
   htot->Draw();
   if (varBins_) htot->GetYaxis()->SetTitle("Events / bin");
   else          htot->GetYaxis()->SetTitle("EventsBIN");
   hpass->Draw("same");
   TLatex label=gfx::cornerLabel(name_+" ("+datasetName+")");
   label.DrawClone();
   // can.SetLogy();
   saver.save(can,name_+"_hists",true,sim);

   // efficiency
   TGraphAsymmErrors *effGr=eff_.CreateGraph();
   effGr->Draw("AP");
   effGr->GetYaxis()->SetNdivisions(507);
   effGr->GetXaxis()->SetLimits(htot->GetXaxis()->GetXmin(),htot->GetXaxis()->GetXmax());
   effGr->GetYaxis()->SetRangeUser(0,1.1);
   if (sim) effGr->SetMarkerStyle(24);

   gfx::LegendEntries le;
   le.append(*effGr,sim ? "MC" : "data","pe");
   TBox box;
   box.SetFillColor(kGray);
   box.SetFillStyle(1001);
   box.SetLineWidth(2);
   box.SetLineStyle(2);
   box.SetLineColor(kGray+2);
   TLine line;
   line.SetLineWidth(2);
   line.SetLineStyle(2);
   line.SetLineColor(kGray+2);
   for (unsigned iLow=0; iLow<plateauCut_.size();iLow++){
      TH1 *htot =eff_.GetCopyTotalHisto();
      TH1 *hpass=eff_.GetCopyPassedHisto();
      int lowerBin=htot->FindBin(plateauCut_[iLow]);
      int upperBin=iLow+1<plateauCut_.size() ? htot->FindBin(plateauCut_[iLow+1])-1 : -1;
      int tot=htot ->Integral(lowerBin,upperBin);
      int pas=hpass->Integral(lowerBin,upperBin);
      float efficiency=pas/float(tot);
      float e_u= TEfficiency::ClopperPearson(tot,pas,0.682689492137,true) -efficiency;
      float e_d=-TEfficiency::ClopperPearson(tot,pas,0.682689492137,false)+efficiency;

      float lowEdge=htot->GetBinLowEdge(lowerBin);
      float uppEdge=(upperBin==-1)?effGr->GetXaxis()->GetXmax():htot->GetBinLowEdge(upperBin+1);
      TString sEff=TString::Format("%.1f^{#plus%.1f}_{#minus%.1f}%% ",efficiency*100,e_u*100,e_d*100);

      box.DrawBox(lowEdge,efficiency-e_d,uppEdge,efficiency+e_u);
      line.DrawLine(lowEdge,efficiency,uppEdge,efficiency);
      le.append(box,sEff,"fl");
      le.buildLegend(lowEdge+.3*(uppEdge-lowEdge),.03,lowEdge+.7*(uppEdge-lowEdge),.25,1,"").DrawClone();
      le.pop_back(); // don't draw this in the next plateau
   }
   label.DrawClone();
   effGr->Draw("same P");
   can.RedrawAxis();

   can.SetLogy(false);
   saver.save(can,name_,true,sim);
}

eff::Efficiency::Efficiency(TString name,TString var,Int_t nbins, Double_t xlow, Double_t xup, std::vector<float> const plateauCut)
   : eff_(name,";"+var+";TREFF",nbins,xlow,xup)
   , name_(name)
   , plateauCut_(plateauCut)
   , varBins_(false)
{}

eff::Efficiency::Efficiency(TString name,TString var,std::vector<float> edges,std::vector<float> widths,std::vector<float> const plateauCut)
   : name_(name)
   , plateauCut_(plateauCut)
   , varBins_(true)
{
   std::vector<double> xbins=hist::getBinVector(edges, widths);
   eff_=TEfficiency(name,";"+var+";TREFF",xbins.size()-1,&xbins[0]);
}

void eff::Efficiency::fill(Bool_t bPassed, Double_t x)
{
   eff_.Fill(bPassed,x);
}

////////////////
// Efficiencies
////////////////

void eff::Efficiencies::add(TString name,TString var,Int_t nbins,Double_t xlow,Double_t xup,std::vector<float> const plateauCut)
{
   effs_[name]=Efficiency(name,var,nbins,xlow,xup,plateauCut);
}

void eff::Efficiencies::add(TString name,TString var,std::vector<float> edges,std::vector<float> widths,std::vector<float> const plateauCut)
{
   effs_[name]=Efficiency(name,var,edges,widths,plateauCut);
}

void eff::Efficiencies::fill(TString name,Bool_t bPassed, Double_t x)
{
   effs_[name].fill(bPassed,x);
}

void eff::Efficiencies::saveAll(io::RootFileSaver const &saver,TString datasetName,bool sim)
{
   // TCanvas can;
   for (auto &it: effs_){
      it.second.save(saver,datasetName,sim);
   }
}
