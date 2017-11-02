#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/efficiency.hpp"

#include <TChain.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTreeReader.h>

static Config const &cfg=Config::get();

class PhotonTrigger
{
public:
   PhotonTrigger(TTreeReader &reader,int threshold,TString name="")
      : threshold(threshold)
      , branch(name!="" ? name : TString::Format("HLT_Photon%d_v",threshold))
      , hasPrescale(reader.GetTree()->GetBranch(branch+"_pre"))
      , trig_(new TTreeReaderValue<bool>(reader,branch))
      , pres_(0)
      {
         if (hasPrescale) pres_=new TTreeReaderValue<int>(reader,branch+"_pre");
      }
   ~PhotonTrigger(){
      if (pres_) delete pres_;
      if (trig_) delete trig_;
   }
   bool fired() {return **trig_;}
   int  pres()  {return (hasPrescale ? **pres_ : 1);}

   int const threshold;
   TString const branch;
   bool const hasPrescale;
private:
   TTreeReaderValue<bool>* trig_;
   TTreeReaderValue<int>*  pres_;
};

static void process(TString const &datasetName,bool sim)
{
   io::log<<"Processing '"+datasetName+"'";
   Dataset const &dataset=cfg.datasets.getDataset(datasetName);
   TChain chain(cfg.treeName);
   for (auto const &dss: dataset.subsets){
      chain.Add(dss.getPath());
   }

   TTreeReader reader(&chain);
   TTreeReaderValue<std::vector<tree::Photon>> photons(reader, "photons");
   TTreeReaderValue<std::vector<tree::Jet>>    jets   (reader, "jets"   );
   TTreeReaderValue<tree::MET> MET(reader, "met");

   TTreeReaderValue<bool> baseMu(reader, "HLT_Mu50_v");
   TTreeReaderValue<bool> baseMET(reader, "HLT_PFMET170_v");
   TTreeReaderValue<bool> Photon165_HE10(reader, "HLT_Photon165_HE10_v");
   TTreeReaderValue<bool> Photon175(reader, "HLT_Photon175_v");

   std::vector<int> thresholds{22,30,36,50,75,90,120};
   std::map<int,PhotonTrigger*> triggers;
   for (int thr: thresholds){
      triggers[thr]=new PhotonTrigger(reader,thr);
   }
   triggers[165]=new PhotonTrigger(reader,165,"HLT_Photon165_HE10_v");
   triggers[175]=new PhotonTrigger(reader,175);
   thresholds.push_back(165);
   thresholds.push_back(175);

   std::map<int, TH1F> hists,hists_pre;
   for (int thr: thresholds){
      TString const name=TString::Format("h%d",thr);
      hists[thr]=TH1F(name,";#gamma_{1} PT;EventsBIN",100,0,250);
      hists_pre[thr]=TH1F(name+"_pre",";#gamma_{1} PT;EventsBIN",100,0,250);
   }

   eff::Efficiencies effs;
   effs.add("bMET/pt","PT",{0,100,300,800},{4,10,100});
   effs.add("bMu/pt" ,"PT",{0,100,300,800},{4,10,100});

   std::map<int,int> nTriggered;

   while (reader.Next()){
      for (int thr: thresholds){
         if (triggers.at(thr)->fired()) nTriggered[thr]++;
      }

      std::vector<tree::Photon const*> lPho;
      for (tree::Photon const& ph: *photons){
         if (!ph.hasPixelSeed and fabs(ph.p.Eta())<1.4442){
            lPho.push_back(&ph);
         }
      }
      if (lPho.size()==0) continue;

      float const pt=lPho[0]->p.Pt();

      float minDR=std::numeric_limits<float>::max();
      for (tree::Jet const &jet: *jets){
         if (jet.isLoose && jet.p.Pt()>30){
            float const dr=jet.p.DeltaR(lPho[0]->p);
            float const dpt=fabs(jet.p.Pt()-pt)/pt;
            if (dr<.1 && dpt<.5) continue;
            minDR=std::min(minDR,dr);
         }
      }

      bool const c_minDR=minDR>.5;
      if (!c_minDR) continue;

      // combine triggers
      bool const allPhotonTriggers=*Photon165_HE10 || triggers.at(22)->fired() || triggers.at(30)->fired() || triggers.at(36)->fired() || triggers.at(50)->fired() || triggers.at(75)->fired() || triggers.at(90)->fired() || triggers.at(120)->fired() || *Photon175;

      // actual efficiency measurements:

      if (*baseMET){
         effs.fill("bMET/pt",allPhotonTriggers,pt);
      } // base met

      if (*baseMu){
         effs.fill("bMu/pt",allPhotonTriggers,pt);
      } // base muon

      // general histograms
      for (int thr: thresholds){
         if (triggers.at(thr)->fired()){
            hists[thr].Fill(pt);
            hists_pre[thr].Fill(pt,triggers.at(thr)->pres());
         }
      }
   } // evt loop

   io::RootFileSaver saver("plots.root",TString::Format("prescales/%s",dataset.name.c_str()));
   effs.saveAll(saver,datasetName,sim);

   for (auto thr_n: nTriggered){
      io::log<<TString::Format("Threshold=%d  N=%d",thr_n.first,thr_n.second);
   }

   TCanvas can;
   double max=0;
   TH1F* hDrawn=(TH1F*)hists[thresholds.back()].DrawClone("axis");
   gfx::LegendEntries le;
   Color::reset();
   for (int thr: thresholds){
      TH1F& h=hists[thr];
      h.SetLineColor(Color::next());
      h.Draw("hist same");
      max=std::max(max,h.GetMaximum());
      le.append(h,TString::Format("%d",thr),"l");
   }
   hDrawn->SetMaximum(max);
   TLegend leg=le.buildLegend(.7,.6,-1,-1,2);
   leg.Draw();
   saver.save(can,"all");

   hDrawn=(TH1F*)hists_pre[thresholds.back()].DrawClone("axis");
   le.clear();
   Color::reset();
   for (int thr: thresholds){
      TH1F& h=hists_pre[thr];
      h.SetLineColor(Color::next());
      h.Draw("hist same");
      max=std::max(max,h.GetMaximum());
      le.append(h,TString::Format("%d",thr),"l");
   }
   hDrawn->SetMaximum(max);
   TLegend leg2=le.buildLegend(.7,.6,-1,-1,2);
   leg2.Draw();
   saver.save(can,"all_pre");

   for (auto m: triggers){
      delete m.second;
   }
}

extern "C"
void run()
{
   process("SingleMuon",false);
   process("MET",false);
   process("SinglePhoton",false);
}
