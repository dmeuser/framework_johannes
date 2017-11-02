/* plot signal cross sections */

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>

#include <regex>

Config const &cfg=Config::get();

enum Scan_t
{
   gluglu,
   GGM,
   CN,
};


std::map<int,float> getXsecs(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==gluglu) path += "/xsec_gluglu.csv";
   if (scan==GGM)  path += "/xsec_wino-bino.csv";
   if (scan==CN)  path += "/xsec_N2C1_wino.csv";
   unsigned const nCol=scan==GGM ? 4 : 3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==GGM) {
            // for key, combine M2 and M1 with a 0 in between
            // debug<<((int)values[0]*100000+(int)values[1]);
            m_M_XS[(int)values[0]*100000+(int)values[1]]=values[2]/1000.; // convert to pb
         } else if (scan==CN) {
            // if (values[0]<350 || values[0]>1000) continue;
            m_M_XS[(int)values[0]]=values[1]/1000.; // convert to pb
         } else {
            if (values[0]<1000) continue;
            m_M_XS[(int)values[0]]=values[1];
         }
      }
   }
   return m_M_XS;
}

std::map<int,float> getXsec_errors(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==gluglu) path += "/xsec_gluglu.csv";
   if (scan==GGM)  path += "/xsec_wino-bino.csv";
   if (scan==CN)  path += "/xsec_N2C1_wino.csv";
   unsigned const nCol=scan==GGM ? 4 : 3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==GGM) {
            // for key, combine M2 and M1 with a 0 in between
            // debug<<((int)values[0]*100000+(int)values[1]);
            m_M_XS[(int)values[0]*100000+(int)values[1]]=values[2]/1000.; // convert to pb
         } else if (scan==CN) {
            // if (values[0]<350 || values[0]>1000) continue;
            m_M_XS[(int)values[0]]=values[2]/values[1]; // make relative
         } else {
            if (values[0]<1000) continue;
            m_M_XS[(int)values[0]]=values[2]/100; // file contains percentual errors
         }
      }
   }
   return m_M_XS;
}

void GluGlu()
{
   TString sScan="gluglu";
   std::map<int,float> xsecs=getXsecs(gluglu);
   std::map<int,float> e_xsecs=getXsec_errors(gluglu);
   TGraph gr(xsecs.size());
   TGraphErrors gr_e(xsecs.size());
   int iP=0;
   for (auto xsec: xsecs) {
      int const m=xsec.first;
      float const xs=xsec.second;
      float const e_xs=e_xsecs[m]*xs;
      gr.SetPoint(iP,m,xs);
      gr_e.SetPoint(iP,m,xs);
      gr_e.SetPointError(iP,0,e_xs);
      // debug*m>>xs;
      // debug<<e_xsecs[m];
      iP++;
   }

   io::RootFileSaver saver("plots.root","signal_xsec");
   TCanvas can;
   can.SetLogy();
   gr_e.SetFillColor(kGray);
   gr_e.SetFillStyle(1001);
   gr_e.Draw("a3");
   gr_e.SetTitle(";m(#tilde{g}) [GeV]; #sigma [pb]");
   gr_e.SetLineWidth(2);
   gr_e.SetLineStyle(7);
   gr_e.GetXaxis()->SetRangeUser(1000,3000);
   gr.Draw("c");
   gr.SetLineWidth(2);
   gr.SetLineStyle(7);

   gfx::LegendEntries le;
   le.append(gr_e,"pp#rightarrow#tilde{g}#tilde{g}","fl");
   le.buildLegend(.6,.86).DrawClone();
   can.RedrawAxis();
   saver.save(can,sScan);
}

void CharginoNeutralino()
{
   TString sScan="CN";
   std::map<int,float> xsecs=getXsecs(CN);
   std::map<int,float> e_xsecs=getXsec_errors(CN);
   TGraph gr(xsecs.size());
   TGraphErrors gr_e(xsecs.size());
   int iP=0;
   for (auto xsec: xsecs) {
      int const m=xsec.first;
      float const xs=xsec.second;
      float const e_xs=e_xsecs[m]*xs;
      gr.SetPoint(iP,m,xs);
      gr_e.SetPoint(iP,m,xs);
      gr_e.SetPointError(iP,0,e_xs);
      // debug*m>>xs;
      // debug<<e_xsecs[m];
      iP++;
   }

   io::RootFileSaver saver("plots.root","signal_xsec");
   TCanvas can;
   can.SetLogy();
   gr_e.SetFillColor(kGray);
   gr_e.SetFillStyle(1001);
   gr_e.Draw("a3");
   gr_e.SetTitle(";m(#tilde{#chi}^{#pm}, #tilde{#chi}^{0}) [GeV]; #sigma [pb]");
   gr_e.SetLineWidth(2);
   gr_e.SetLineStyle(7);
   gr_e.GetXaxis()->SetRangeUser(400,1000);
   gr_e.GetYaxis()->SetRangeUser(8e-4,2e-1);
   gr.Draw("c");
   gr.SetLineWidth(2);
   gr.SetLineStyle(7);

   gfx::LegendEntries le;
   le.append(gr_e,"pp#rightarrow#tilde{#chi}^{#pm}#tilde{#chi}^{0}","fl");
   le.buildLegend(.6,.86).DrawClone();
   can.RedrawAxis();
   saver.save(can,sScan);
}

void ggm()
{
   TString sScan="GGM";
   std::map<int,float> xsecs=getXsecs(GGM);

   TH2F h("",";m_{#tilde{B}} [GeV];m_{#tilde{W}} [GeV];#sigma [pb]", 33, 205-12.5, 1005+12.5, 33, 215-12.5, 1015+12.5); // 32 points in each direction
   for (auto xsec: xsecs) {
      int const key=xsec.first;
      int const m2=key/100000;
      int const m1=key%100000;
      float const xs=xsec.second;
      h.SetBinContent(h.FindBin(m1,m2),xs);
      // debug*m2>>m1;
      // debug<<xs;
   }
   io::RootFileSaver saver("plots.root","signal_xsec");
   TCanvas can;
   can.SetLogz();
   can.SetRightMargin (.15);
   can.SetBottomMargin(.22);
   h.Draw("colz");
   h.SetMinimum(1e-3);
   gfx::setupZaxis(h);
   gfx::cornerLabel("GGM",4).DrawClone();
   saver.save(can,sScan);
}

extern "C"
void run()
{
   GluGlu();
   CharginoNeutralino();
   ggm();
}
