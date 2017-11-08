#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

#include <TStyle.h>
#include <TGraph2D.h>
#include <TGraphSmooth.h>

Config const &cfg=Config::get();

enum Scan_t
{
   T5gg,
   T5Wg,
   TChiWg,
   T6gg,
   T6Wg,
   TChiNg,
   GGM
};

void transpose(TGraph &gr)
{
   double x,y;
   for (int i=0; i<gr.GetN(); i++) {
      gr.GetPoint(i,x,y);
      gr.SetPoint(i,y,x);
   }
}

void removeErrors(TGraphErrors &gr)
{
   for (int i=0; i<gr.GetN(); i++) {
      gr.SetPointError(i,0,0);
   }
}

void drawContours(TGraph2D &grAcc, Scan_t scan)
{
   double min=grAcc.GetZmin();
   double max=grAcc.GetZmax();
   double xRange=grAcc.GetXmax()-grAcc.GetXmin();
   double yRange=grAcc.GetYmax()-grAcc.GetYmin();

   TLatex label=gfx::cornerLabel("",1);
   label.SetNDC(false);
   label.SetTextAlign(22);
   label.SetTextSize(.03);
   for (float lvl=0.1; lvl < 1.0; lvl+=0.1) {
      if (lvl<min || lvl>max) continue;
      TList *conts=grAcc.GetContourList(lvl);
      if (!conts) continue;
      TGraphSmooth gs;
      TGraph *gr=(TGraph*)conts->First();
      if (scan==GGM) transpose(*gr);
      // gr = gs.Approx(gr,"linear");
      // gr = gs.SmoothSuper(gr);
      gr = gs.SmoothLowess(gr);
      removeErrors(*(TGraphErrors*)gr);
      if (scan==GGM) transpose(*gr);
      gr->SetLineWidth(2);
      gr->DrawClone("same l");
      double x,y;
      gr->GetPoint(gr->GetN()*2/3,x,y);
      if (scan==GGM) {
         x+=.05*xRange;
      } else {
         y+=.02*yRange;
      }
      label.DrawLatex(x,y,TString::Format("%.0f%%",lvl*100));
   }
}

void drawContour(TGraph2D &gr, Scan_t scan, float lvl)
{
   double min=gr.GetZmin();
   double max=gr.GetZmax();
   double xRange=gr.GetXmax()-gr.GetXmin();
   double yRange=gr.GetYmax()-gr.GetYmin();

   TLatex label=gfx::cornerLabel("",1);
   label.SetNDC(false);
   label.SetTextAlign(22);
   label.SetTextSize(.03);

   if (lvl<min || lvl>max) return;
   TList *conts=gr.GetContourList(lvl);
   if (!conts) return;
   TGraphSmooth gs;
   TGraph *cont=(TGraph*)conts->First();
   // if (scan==GGM) transpose(*cont);
   // cont = gs.Approx(cont,"linear");
   // cont = gs.SmoothSuper(cont);
   // cont = gs.SmoothLowess(cont);
   // removeErrors(*(TGraphErrors*)cont);
   // if (scan==GGM) transpose(*cont);
   cont->SetLineWidth(2);
   cont->DrawClone("same l");
   double x,y;
   cont->GetPoint(cont->GetN()*2/3,x,y);
   if (scan==GGM) {
      x+=.05*xRange;
   } else {
      y+=.02*yRange;
   }
   label.DrawLatex(x,y,TString::Format("%.0f%%",lvl*100));

}

void drawLabels(Scan_t scan,TString text)
{
   TString sScan;
   if (scan==T5gg) sScan="T5gg";
   if (scan==T5Wg) sScan="T5Wg";
   if (scan==TChiWg) sScan="TChiWg";
   if (scan==T6gg) sScan="T6gg";
   if (scan==T6Wg) sScan="T6Wg";
   if (scan==TChiNg) sScan="TChiNg";

   TLatex label= scan==GGM ? gfx::cornerLabel(text,4) : gfx::cornerLabel(text,1);
   if (scan==GGM) label.SetY(label.GetY()+.05);
   label.DrawClone();
   label.SetTitle(sScan);
   label.SetY(label.GetY()-.05);
   label.DrawClone();
}

void run(Scan_t scan)
{
   TString sScan;
   if (scan==T5gg) sScan="T5gg";
   if (scan==T5Wg) sScan="T5Wg";
   if (scan==TChiWg) sScan="TChiWg";
   if (scan==T6gg) sScan="T6gg";
   if (scan==T6Wg) sScan="T6Wg";
   if (scan==TChiNg) sScan="TChiNg";
   if (scan==GGM) sScan="GGM";  
   TString title;
   if (scan==T5gg) title=";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV];m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} [GeV];";
   if (scan==T5Wg) title=";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV];m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} [GeV];";
   if (scan==TChiWg) title=";m_{NLSP} [GeV];";
   if (scan==T6gg) title=";m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} [GeV];m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} [GeV];";
   if (scan==T6Wg) title=";m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} [GeV];m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} [GeV];";
   if (scan==TChiNg) title=";m_{NLSP} [GeV];";
   if (scan==GGM) title=";m_{#tilde{B}} (GeV);m_{#tilde{W}} (GeV)";  
   io::RootFileReader fileReader(TString::Format("signal_scan_%s.root",cfg.treeVersion.Data()),"pre_ph165");

   io::RootFileSaver saver("plots.root","");
   TCanvas can;

   if (scan==TChiWg || scan == TChiNg) {
      TGraph2D grAcc(*fileReader.read<TGraph2D>("c_MET300/MT300/STg/"+sScan+"_acceptance"));
      TGraph grA(grAcc.GetN());
      grA.SetTitle(title+"Acceptance");
      double *x=grAcc.GetX();
      double *z=grAcc.GetZ();
      for (int i=0; i<grAcc.GetN(); i++) {
         grA.SetPoint(i,x[i],z[i]);
      }
      grA.Sort();
      grA.Draw("apl");
      drawLabels(scan,"Acceptance");
      saver.save(can,"acceptance/"+sScan);

      TGraph2D grScaleUnc(*fileReader.read<TGraph2D>("c_MET300/MT300/STg/"+sScan+"_scaleUnc"));
      TGraph grSU(grScaleUnc.GetN());
      grSU.SetTitle(title+"Scale uncertainty");
      x=grScaleUnc.GetX();
      z=grScaleUnc.GetZ();
      for (int i=0; i<grScaleUnc.GetN(); i++) {
         grSU.SetPoint(i,x[i],100*z[i]);
      }
      grSU.Sort();
      grSU.Draw("apl");
      grSU.GetYaxis()->SetRangeUser(0,20);
      grSU.GetYaxis()->SetNoExponent();
      drawLabels(scan,"Scale uncertainty in % on acceptance");
      saver.save(can,"scaleUnc/"+sScan);

      TGraph2D grCont(*fileReader.read<TGraph2D>("c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/"+sScan+"_contamination"));
      TGraph grC(grCont.GetN());
      grC.SetTitle(title+"Signal fraction");
      x=grCont.GetX();
      z=grCont.GetZ();
      for (int i=0; i<grCont.GetN(); i++) {
         grC.SetPoint(i,x[i],z[i]);
      }
      grC.Sort();
      grC.Draw("apl");
      drawLabels(scan,"CR signal contamination");
      saver.save(can,"contamination/"+sScan);
   } else {
      // Acceptance
      TGraph2D grAcc(*fileReader.read<TGraph2D>("c_MET300/MT300/STg/"+sScan+"_acceptance"));
      grAcc.SetTitle(title+"Acceptance");
      TH2D hAcc(*grAcc.GetHistogram());

      grAcc.Draw("tri1");
      saver.save(grAcc,"acceptance/"+sScan);

      hAcc.Draw("axis");
      drawContours(grAcc,scan);
      drawLabels(scan,"Acceptance");
      saver.save(can,"acceptance/"+sScan+"_cont");

      can.SetRightMargin (.15);
      can.SetBottomMargin(.22);
      hAcc.Draw("colz");
      gfx::setupZaxis(hAcc,false);
      // drawContours(grAcc,scan);
      drawLabels(scan,"Acceptance");
      if (scan==GGM) {
         hAcc.GetYaxis()->SetRangeUser(265,1050);
      } else if (scan==T5gg || scan==T5Wg ) {
         hAcc.GetXaxis()->SetRangeUser(1275,2500);
      }
      else if (scan==T6gg || scan==T6Wg ) {
         hAcc.GetXaxis()->SetRangeUser(1100,2150);
      }
      saver.save(can,"acceptance/"+sScan+"_colz");

      if (scan!=GGM) {
         TGraph2D grScaleUnc(*fileReader.read<TGraph2D>("c_MET300/MT300/STg/"+sScan+"_scaleUnc"));
         grScaleUnc.SetTitle(title);
         TH2D hScaleUnc(*grScaleUnc.GetHistogram());
         hScaleUnc.Scale(100.);
         hScaleUnc.Draw("colz");
         gfx::setupZaxis(hScaleUnc,false);
         drawLabels(scan,"Scale uncertainty in % on acceptance");
         saver.save(can,"scaleUnc/"+sScan+"_colz");
      }

      // Signal Contamination
      TGraph2D grCont(*fileReader.read<TGraph2D>("c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/"+sScan+"_contamination"));
      grCont.SetTitle(title+"Signal fraction");
      TH2D hCont(*grCont.GetHistogram());

      hCont.Draw("colz");
      gfx::setupZaxis(hCont,false);
      drawContour(grCont,scan,.05);
      drawLabels(scan,"CR signal contamination");
      saver.save(can,"contamination/"+sScan+"_colz");
   }
}

extern "C"
void run()
{
   run(GGM);
 //  run(T5gg);
 //  run(T5Wg);
 //  run(TChiWg);
 //  run(T6gg);
 //  run(T6Wg);
 //  run(TChiNg);
}
