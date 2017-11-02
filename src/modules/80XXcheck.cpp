#include "Config.hpp"
#include "tools/io.hpp"
#include "tools/gfx.hpp"

#include <TFile.h>
#include <TTree.h>


void run(TString sVar, TString label)
{
   Config cfg=Config::get();

   TFile fNew("/.automount/home/home__home4/institut_1b/lange/run2/dl/80X/fudge/T5Wg_1350_1000.root");
   TFile fOld("/.automount/home/home__home4/institut_1b/lange/run2/v13/T5Wg_1350_1000.root");
   if (fOld.IsZombie()||fNew.IsZombie()) {
      return;
   }
   TTree *tNew=(TTree*)fNew.Get("TreeWriter/eventTree");
   TTree *tOld=(TTree*)fOld.Get("TreeWriter/eventTree");

   TCanvas can;
   can.SetLogy();
   io::RootFileSaver saver("plots.root","80XXcheck",true);

   TString presel="photons[0].p.Pt()>180";

   tNew->Draw(sVar,presel);
   TH1F hNew(*(TH1F*)gPad->GetPrimitive("htemp"));
   TH1F hOld(hNew);
   hOld.SetName("hOld");
   tOld->Draw(sVar+">>hOld",presel);

   hNew.Scale(1./hNew.Integral());
   hOld.Scale(1./hOld.Integral());
   hNew.SetLineColor(kRed);
   hOld.SetLineColor(kBlack);
   hNew.SetLineWidth(2);
   hOld.SetLineWidth(2);

   hNew.SetTitle(";"+label+";normalized events");
   hNew.Draw();
   hOld.Draw("same");

   gfx::LegendEntries le;
   le.append(hOld,"74X","l");
   le.append(hNew,"80X","l");
   le.buildLegend(.7,.75).DrawClone();

   saver.save(can,label,true,false);

   fNew.Close();
   fOld.Close();
}

extern "C"
void run()
{
   run("met.p.Pt()","MET");
   run("met.sig"   ,"METSIG");
   run("photons[0].p.Pt()","#gamma_{1} PT");
}
