#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

#include <TGraph2D.h>
#include <TGraphSmooth.h>

Config const &cfg=Config::get();

enum Scan_t
{
   TChiWg,
   TChiNg,
   T5gg,
   T5Wg,
   T6gg,
   T6Wg,
   GGM,
};
std::map<int,float> getXsecs(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==TChiWg) path += "/xsec_N2C1_wino.csv";
   if (scan==TChiNg) path += "/xsec_comb_wino.csv";  
   unsigned const nCol=3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==TChiWg || scan==TChiNg) {
            m_M_XS[(int)values[0]]=values[1]/1000.; // convert to pb
         } else {
            debug<<"unkown scan";
         }
      }
   }
   return m_M_XS;
}

std::map<int,float> getXsec_unc(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==TChiWg) path += "/xsec_N2C1_wino.csv";
   if (scan==TChiNg) path += "/xsec_comb_wino.csv";
   unsigned const nCol=3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==TChiWg || scan==TChiNg) {
            m_M_XS[(int)values[0]]=values[2]/1000.; // convert to pb
         } else {
            debug<<"unkown scan";
         }
      }
   }
   return m_M_XS;
}

void plot_1d(Scan_t scan)
{
   assert(scan==TChiWg || scan==TChiNg);
   TString sScan;
   if (scan==TChiWg) sScan="TChiWg";
   if (scan==TChiNg) sScan="TChiNg";

   TString title, title2;
   if (scan==TChiWg)  title="pp #rightarrow"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, ";
   if (scan==TChiWg)  title2=
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma/W^{#pm}#tilde{G}";
   if (scan==TChiNg)  title="pp #rightarrow"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}/"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, ";
   if (scan==TChiNg)  title2=
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}#rightarrow"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}+soft, "                        
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma/(H,Z)#tilde{G}";
   std::map<int,float> mXsecs=getXsecs(scan);
   std::map<int,float> mXsec_unc=getXsec_unc(scan);

   io::RootFileReader fileReader(TString::Format("limits_%s.root",sScan.Data()));
   io::RootFileSaver aux_saver(TString::Format("aux_limits_%s.root",sScan.Data()),"");
   std::map<TString,TGraph> gr;
   for (TString lvl: {"exp","exp+1","exp-1","exp+2","exp-2","obs"}) {
      TGraph2D gr2(*fileReader.read<TGraph2D>("gr_"+lvl));    
      gr[lvl]=TGraph(gr2.GetN());
      gr[lvl].SetLineWidth(2);
      double *x=gr2.GetX();
      double *z=gr2.GetZ();
      for (int i=0; i<gr2.GetN(); i++) {
         gr[lvl].SetPoint(i,x[i],z[i]*mXsecs[int(x[i])]);
      }
      gr[lvl].Sort();
      aux_saver.save(gr[lvl],lvl);
   }
   TGraph2D grSig(*fileReader.read<TGraph2D>("gr_significance"));
   TGraph gr1DSig=TGraph(grSig.GetN());
   double *xSig=grSig.GetX();
   double *zSig=grSig.GetZ();
   for (int i=0; i<grSig.GetN(); i++) {
      gr1DSig.SetPoint(i,xSig[i],zSig[i]);
   }
   gr1DSig.Sort();

   for (TString lvl: {"1","2"}) {
      gr[lvl]=TGraph(gr["exp+"+lvl]);
      double x, y;
      for (int i=gr["exp-"+lvl].GetN()-1; i>=0; i--) {
         gr["exp-"+lvl].GetPoint(i,x,y);
         gr[lvl].SetPoint(gr[lvl].GetN(),x,y);
      }
      gr[lvl].GetPoint(0,x,y);
      gr[lvl].SetPoint(gr[lvl].GetN(),x,y);
   }
   gr["xs"]=TGraph(mXsecs.size());
   gr["xs"].SetLineWidth(2);
   gr["xs+1"]=TGraph(gr["xs"]);
   gr["xs-1"]=TGraph(gr["xs"]);
   int iP=0;
   for (auto m_xs :mXsecs) {
      int const m=m_xs.first;
      gr["xs"].SetPoint(iP,m,m_xs.second);
      gr["xs+1"].SetPoint(iP,m,m_xs.second+mXsec_unc[m]);
      gr["xs-1"].SetPoint(iP,m,m_xs.second-mXsec_unc[m]);
      iP++;
   }

   aux_saver.save(gr["xs"],"xs");
   aux_saver.save(gr["xs+1"],"xs+1");
   aux_saver.save(gr["xs-1"],"xs-1");

   io::RootFileSaver saver("plots.root","limit_1d");
   TCanvas can;
   can.SetLogy();

   gr["2"].Draw("afl");
   gr["2"].SetFillColor(kYellow-7);
   gr["2"].SetLineColor(kYellow-7);
   gr["1"].Draw("fl");
   gr["1"].SetFillColor(kGreen-7);
   gr["1"].SetLineColor(kGreen-7);
   gr["exp"].Draw("l");
   gr["exp"].SetLineStyle(2);
   gr["exp"].SetLineColor(kRed);
   gr["obs"].Draw("l");
   gr["obs"].SetLineColor(kBlack);
   gr["xs"].Draw("l");
   gr["xs"].SetLineColor(kBlue);
   gr["xs"].SetLineStyle(2);
   gr["xs+1"].Draw("l");
   gr["xs+1"].SetLineColor(kBlue);
   gr["xs-1"].Draw("l");
   gr["xs-1"].SetLineColor(kBlue);

   if (scan == TChiNg) gr["2"].GetXaxis()->SetRangeUser(300,1150);
   if (scan == TChiWg) gr["2"].GetXaxis()->SetRangeUser(300,1000);
   if (scan == TChiNg) gr["2"].GetYaxis()->SetRangeUser(0.0004,0.8);
   if (scan == TChiWg) gr["2"].GetYaxis()->SetRangeUser(0.0008,0.4);
   gr["2"].GetYaxis()->SetTitleSize(.055);
   gr["2"].GetYaxis()->SetTitleOffset(1.5);
   gr["2"].SetTitle(";m_{NLSP} [GeV];95% CL upper limit / cross section [pb]   ");

   gfx::LegendEntries le;
   le.append(gr["exp"],"exp.","l");
   le.append(gr["1"],"#pm 1 #sigma","f");
   le.append(gr["2"],"#pm 2 #sigma","f");
   TLegend l1=le.buildLegend(.5,.85,.97,.92,3);
   l1.SetColumnSeparation(-0.2);
   l1.DrawClone();
   le.clear();
   le.append(gr["obs"],"obs. limit","l");
   le.buildLegend(.5,.78,.97,.85,2).DrawClone();
   le.clear();
   le.append(gr["xs"],"cross section","l");
   le.append(gr["xs+1"],"#pm 1 #sigma","l");
   TLegend l2=le.buildLegend(.5,.71,.97,.78,2);
   l2.SetColumnSeparation(-0.2);
   l2.DrawClone();


   TLatex l=gfx::cornerLabel(title,3);
   l.SetY(0.22);
   l.DrawClone();
   TLatex l3=gfx::cornerLabel(title2,3);
   l3.SetY(0.17);
   l3.DrawClone();

   can.RedrawAxis();
   saver.save(can,sScan,true,false);

   TCanvas can2;
   gr1DSig.Draw("al");
   gr1DSig.SetLineColor(kBlue);
   gr1DSig.SetLineWidth(2);
   gr1DSig.SetTitle(";m_{NLSP} [GeV];Significance (s.d.)");
   gr1DSig.GetYaxis()->SetTitleOffset(0.8);
   if (scan == TChiNg) gr1DSig.GetXaxis()->SetRangeUser(300,1150);
   if (scan == TChiWg) gr1DSig.GetXaxis()->SetRangeUser(300,1000);
   if (scan == TChiNg) gr1DSig.GetYaxis()->SetRangeUser(-3,3);
   if (scan == TChiWg) gr1DSig.GetYaxis()->SetRangeUser(-3,3);
   l.DrawClone();
   l3.DrawClone();
   can2.RedrawAxis();
   saver.save(can2,sScan+"_significance",true,false);
}

// save significance for 2D
void plot_significance2D(Scan_t scan)
{
   assert(scan!=TChiWg);
   TString sScan;
   if (scan==T5gg) sScan="T5gg";
   if (scan==T5Wg) sScan="T5Wg";
   if (scan==T6gg) sScan="T6gg";
   if (scan==T6Wg) sScan="T6Wg";
   if (scan==GGM) sScan="GGM";

   io::RootFileReader fileReader(TString::Format("limits_%s.root",sScan.Data()));
   io::RootFileSaver saver("plots.root","limit_2d");

   TH2F grSig(*fileReader.read<TH2F>("h_significance"));
   
   TCanvas can;
  // can.SetLogy();

   grSig.Draw("colz");
   if (scan==T5gg) grSig.SetTitle("m_{#tilde{g}};m_{#tilde{#chi_{1}^{0}}}");
   if (scan==T5Wg) grSig.SetTitle("m_{#tilde{g}};m_{#tilde{#chi_{1}^{0/#pm}}}");
   if (scan==T6gg) grSig.SetTitle("m_{#tilde{q}};m_{#tilde{#chi_{1}^{0}}}");
   if (scan==T6Wg) grSig.SetTitle("m_{#tilde{q}};m_{#tilde{#chi_{1}^{0/#pm}}}");
   grSig.GetZaxis()->SetTitle("Significance (s.d.)");
   grSig.GetZaxis()->SetRangeUser(-3,3);
   grSig.GetXaxis()->SetRangeUser(800,2500);
   grSig.GetYaxis()->SetTitleOffset(0.8);

   can.RedrawAxis();

   saver.save(can,sScan+"_significance",true,false);


}

// save limits for auxiliary material
void save_aux(Scan_t scan)
{
   assert(scan!=TChiWg);
   TString sScan;
   if (scan==T5gg) sScan="T5gg";
   if (scan==T5Wg) sScan="T5Wg";
   if (scan==T6gg) sScan="T6gg";
   if (scan==T6Wg) sScan="T6Wg";
   if (scan==GGM) sScan="GGM";

   io::RootFileReader fileReader(TString::Format("limits_%s.root",sScan.Data()));
   io::RootFileSaver aux_saver(TString::Format("aux_limits_%s.root",sScan.Data()),"");

   std::map<TString,TGraph> gr;
   for (TString lvl: {"exp","exp+1","exp-1","obs","obs+1","obs-1"}) {
      TGraph gr(*fileReader.read<TGraph>("gr_"+lvl+"C_sm"));
      aux_saver.save(gr,lvl);
   }
   if (scan==T5gg || scan==T5Wg) {
      TH2F h2(*fileReader.read<TH2F>("h_obs_xs_redone"));
      aux_saver.save(h2,"h_obs_xs");
   } else {
      TH2F h2(*fileReader.read<TH2F>("h_obs_xs"));
      aux_saver.save(h2,"h_obs_xs");
   }
}

extern "C"
void run()
{
  // plot_1d(TChiWg);
 //  plot_1d(TChiNg);
   plot_significance2D(T5gg);
   
 //  save_aux(T5gg);
 //  save_aux(T5Wg);
 //  save_aux(GGM);
}
