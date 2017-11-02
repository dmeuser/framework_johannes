#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/efficiency.hpp"

#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTreeReader.h>

static Config const &cfg=Config::get();
static io::Logger logFile("trigger.log");

// isolation thresholds for Iso40 in HLT
static TF1 fIso40ec("fIso40ec","4.0+0.012*x",0,1000);
static TF1 fIso40hc("fIso40hc","4.0+0.005*x",0,1000);
static TF1 fIso40tr("fIso40tr","4.0+0.002*x",0,1000);

static void process(TString const &datasetName,bool sim)
{
   Dataset const &dataset=cfg.datasets.getDataset(datasetName);
   TChain chain(cfg.treeName);
   for (auto const &dss: dataset.subsets){
      chain.Add(dss.getPath());
   }

   TTreeReader reader(&chain);
   TTreeReaderValue<std::vector<tree::Photon>> photons(reader, "photons");
   TTreeReaderValue<std::vector<tree::Jet>>    jets   (reader, "jets"   );
   TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
   TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
   TTreeReaderValue<tree::MET> MET(reader, "met");

   TTreeReaderValue<bool> sigTr(reader, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v");
   TTreeReaderValue<bool> basePhTr (reader, "HLT_Photon36_R9Id90_HE10_IsoM_v");
   TTreeReaderValue<bool> baseMETTr(reader, "HLT_PFMET170_HBHECleaned_v");
   TTreeReaderValue<bool> baseMuTr (reader, "HLT_Mu50_v");
   TTreeReaderValue<bool> Photon165_HE10(reader, "HLT_Photon165_HE10_v");
   TTreeReaderValue<bool> Photon175(reader, "HLT_Photon175_v");
   TTreeReaderValue<bool> Photon135MET100_JetIdCleaned(reader,"HLT_Photon135_PFMET100_JetIdCleaned_v"); // early
   TTreeReaderValue<bool> Photon135MET100_normal(reader,"HLT_Photon135_PFMET100_v");                    // late
   TTreeReaderValue<bool> Photon135MET100_NoiseCleaned(reader,"HLT_Photon135_PFMET100_NoiseCleaned_v"); // MC
   TTreeReaderValue<bool> HT200(reader,"HLT_PFHT200_v");
   TTreeReaderValue<bool> HT250(reader,"HLT_PFHT250_v");
   TTreeReaderValue<bool> HT300(reader,"HLT_PFHT300_v");
   TTreeReaderValue<bool> HT350(reader,"HLT_PFHT350_v");
   TTreeReaderValue<bool> HT400(reader,"HLT_PFHT400_v");
   TTreeReaderValue<bool> HT475(reader,"HLT_PFHT475_v");
   TTreeReaderValue<bool> HT600(reader,"HLT_PFHT600_v");
   TTreeReaderValue<bool> HT650(reader,"HLT_PFHT650_v");
   TTreeReaderValue<bool> HT800(reader,"HLT_PFHT800_v");


   TTreeReaderValue<UInt_t> runNo(reader, "runNo");

   eff::Efficiencies effs;
   effs.add("pt","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("pt c_minDR","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("pt c_minDR Iso40","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("pt_tight","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("pt_noR9","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("minDR","min#DeltaR(#gamma_{1},jet)",{0,.5,3,5},{.1,.5,1});
   effs.add("met","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("S","METSIG",{0,20,50,200},{5,10,50});
   effs.add("met c_minDR","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("met_pt50-180 c_minDR","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("met_pt>50","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("met_tight","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("sieie","SIEIE",30,0,.02);

   effs.add("Photon175 bMu","PT",{0,500},{20},{200});
   effs.add("Photon165_HE10 bMu","PT",{0,500},{20},{180});
   effs.add("OR_Photon165_HE10","PT",{0,50,120,300,800},{5,20,30,100});
   effs.add("+Photon165_HE10","PT",{0,50,120,300,800},{5,20,30,100},{50,180});
   effs.add("pt bMu","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("pt bMu c_minDR","PT",{0,50,300,800},{5,50,100},{50});
   effs.add("met bMu","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("met bMu c_minDR","MET",{0,120,200,500},{10,20,150},{100});
   effs.add("sieie bMu","SIEIE",30,0,.02);

   effs.add("g135m100 pt" ,"PT",{0,500},{20},{140});
   effs.add("g135m100 met","MET",{0,500},{20},{240});
   effs.add("g135m100 S"  ,"METSIG",{0,1000},{20});

   // the final efficiencies in relevant selections
   effs.add("final/PhMet_pt","PT",{20,50,130,180},{5,20,50},{50});
   effs.add("final/PhMet_met","MET",{0,120,200},{10,20},{100});
   effs.add("final/Ph","PT",{0,1000},{20},{180});

   // dependencies in plateau
   effs.add("ph165/met","MET",{0,500},{50});
   effs.add("ph165/METS","METSIG",{0,100,200,600},{20,50,100});
   effs.add("ph165/STg","STg",{0,800,1500},{50,100});
   effs.add("ph165/MT","MT",{0,1000},{50});
   effs.add("ph165/pt bHt","PT",{0,500,1000,1500},{20,50,100},{180});

   if (sim){
      effs.add("pt bMet100cut c_minDR","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("pt bMet100cut c_genM c_minDR","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("pt bMet100-120cut c_minDR","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("pt bMet120-140cut c_minDR","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("pt bMet140-160cut c_minDR","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("pt bMet160cut c_minDR","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("pt bMet100cut c_minDR c_Iso40","PT",{0,50,300,2000},{5,50,100},{50});
      effs.add("sieie bMet100cut","SIEIE",30,0,.02);

      effs.add("OR_Photon165_HE10 bMet100cut","PT",{0,50,120,300,800},{5,20,30,100});
      effs.add("+Photon165_HE10 bMet100cut","PT",{0,50,120,300,800},{5,20,30,100});

      effs.add("Photon36_R9Id90_HE10_IsoM","PT",{0,50,300,2000},{5,50,100});
      effs.add("Photon36_R9Id90_HE10_IsoM bMet100cut","PT",{0,50,300,2000},{5,50,100});

      // efficiency of trigger isolation
      effs.add("Iso40ec","PT",{0,50,300,2000},{5,50,100});
      effs.add("Iso40hc","PT",{0,50,300,2000},{5,50,100});
      effs.add("Iso40tr","PT",{0,50,300,2000},{5,50,100});
   }

   TEfficiency eDpTdR("eDpTdR",";dpt;dr;TREFF",20,0,2,20,0,.6);
    
   
   std::map<int,int> mRunTotal, mRunPass;

   // for the muon base trigger
   int iTotal=0;
   int iPass =0;

   while (reader.Next()){
      std::vector<tree::Photon const*> lPho,tPho;
      for (tree::Photon const& ph: *photons){
         if (!ph.hasPixelSeed and fabs(ph.p.Eta())<1.4442){
            if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
            if (ph.isLoose) lPho.push_back(&ph);
            if (ph.isTight) tPho.push_back(&ph);
         }
      }
      if (lPho.size()==0) continue;
      
      std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);

      float const pt=lPho[0]->p.Pt();
      float const r9=lPho[0]->r9;
      float const met=MET->p.Pt();
      float STg=met;
      float const MT=phys::M_T(*lPho[0],*MET);
      for (tree::Photon const *ph: lPho){
         STg+=ph->p.Pt();
      }

      float minDR=std::numeric_limits<float>::max();
      for (tree::Jet const &jet: *jets){
         if (jet.isLoose && jet.p.Pt()>30){
            float const dr=jet.p.DeltaR(lPho[0]->p);
            float const dpt=fabs(jet.p.Pt()-pt)/pt;
            eDpTdR.Fill(*sigTr,dpt,dr);
            if (dr<.1 && dpt<.5) continue;
            minDR=std::min(minDR,dr);
         }
      }


      bool clean_MET = true;
      
      for (auto const &jet: cjets) {
         if (jet.p.Pt() < 100) continue;            
         if (std::fabs(MET->p.DeltaPhi(jet.p)) < 0.3) clean_MET = false;
      }

      if (!clean_MET) continue;

      // require photon to have fired the trigger
      bool triggerMatch=false; //false
      for (tree::Particle const& to: *triggerObjects){
         if (to.p.DeltaR(lPho[0]->p)<0.1 && std::fabs(to.p.Pt()-lPho[0]->p.Pt())/lPho[0]->p.Pt()<0.7) {
            triggerMatch=true;
         }
      }
      
      bool const ph165_he10=*Photon165_HE10 && triggerMatch;

      // combine triggers
      bool const Photon135MET100= *Photon135MET100_JetIdCleaned || *Photon135MET100_normal || *Photon135MET100_NoiseCleaned;
      bool const c_minDR=minDR>.5;
      bool const c_pt   =pt>50;
      bool const c_r9   =r9>.9;
      bool const c_met  =met>100;
      // bool const c_Iso40ec=lPho[0]->ecalPFClIso<fIso40ec(pt);
      // bool const c_Iso40hc=lPho[0]->hcalPFClIso<fIso40hc(pt);
      // bool const c_Iso40tr=lPho[0]->trackIso   <fIso40tr(pt);
      // bool const c_Iso40=c_Iso40ec && c_Iso40hc && c_Iso40tr;

      // actual efficiency measurements:

      bool const phMet_or_165=*sigTr||ph165_he10;
      bool const phMet_plus_165= (pt<180) ? *sigTr : ph165_he10;
      bool const allHT = *HT200 ||*HT250 ||*HT300 ||*HT350 ||*HT400 ||*HT475 ||*HT600 ||*HT650 ||*HT800;

      if (*basePhTr){
         if (c_r9){
            effs.fill("met",*sigTr,met);
            if (c_minDR){
               effs.fill("met c_minDR",*sigTr,met);
               if (pt>50 && pt<180){
                  effs.fill("met_pt50-180 c_minDR",*sigTr,met);
                  effs.fill("final/PhMet_met",*sigTr,met);
                  effs.fill("S",*sigTr,MET->sig);
               }
            }
            if (pt>50) effs.fill("met_pt>50",*sigTr,met);
         }
         if (tPho.size()!=0 && tPho[0]->r9>.9) effs.fill("met_tight",*sigTr,met);
      } // base photon

      // Met trigger not in simulation:
      bool baseMET=sim ? (met>170) : *baseMETTr;
      if (baseMET && met > 100 && MT > 100){
         if (tPho.size()!=0 && tPho[0]->r9>.9){
            effs.fill("pt_tight",*sigTr,tPho[0]->p.Pt());
         }

         if (c_minDR) effs.fill("final/Ph",ph165_he10,pt);
                   
         effs.fill("pt_noR9",*sigTr,pt);

         if (c_r9){
            effs.fill("pt",*sigTr,pt);
            if (c_minDR){
               effs.fill("pt c_minDR",*sigTr,pt);
               effs.fill("OR_Photon165_HE10",phMet_or_165,pt);
               effs.fill("+Photon165_HE10",phMet_plus_165,pt);
               if (c_met) {
                  if (pt<=180) effs.fill("final/PhMet_pt",*sigTr,pt);
               }
               // if (c_Iso40) effs.fill("pt c_minDR Iso40",*sigTr,pt);
            }
            if (c_pt){
               effs.fill("minDR",*sigTr,minDR);
               effs.fill("sieie",*sigTr,lPho[0]->sigmaIetaIeta);
            }
         }
         if (met>240 && c_minDR) effs.fill("g135m100 pt",Photon135MET100,pt);
      } // base met

      if (allHT){
         if (c_minDR) {
            effs.fill("ph165/pt bHt",ph165_he10,pt);
            if (pt>180) {
               effs.fill("ph165/met",ph165_he10,met);
               effs.fill("ph165/METS",ph165_he10,MET->sig);
               effs.fill("ph165/STg",ph165_he10,STg);
               effs.fill("ph165/MT",ph165_he10,MT);
            }
         }
      } // base HT

      if (sim && met>100 && c_r9) {
         if (c_minDR){
            effs.fill("pt bMet100cut c_minDR",*sigTr,pt);
            effs.fill("OR_Photon165_HE10 bMet100cut",phMet_or_165,pt);
            effs.fill("+Photon165_HE10 bMet100cut",phMet_plus_165,pt);
            // if (c_Iso40) effs.fill("pt bMet100cut c_minDR c_Iso40",*sigTr,pt);

            bool c_genM=false;
            float dR,dPt;
            for (auto const &genp:*genParticles){
               if (genp.pdgId!=22) continue; // no photon
               dR=lPho[0]->p.DeltaR(genp.p);
               dPt=fabs(lPho[0]->p.Pt()-genp.p.Pt())/lPho[0]->p.Pt();
               if (dR<.1 && dPt<.1){
                  c_genM=true;
                  break;
               }
            }
            if (c_genM) effs.fill("pt bMet100cut c_genM c_minDR",*sigTr,pt);

            if      (met<120) effs.fill("pt bMet100-120cut c_minDR",*sigTr,pt);
            else if (met<140) effs.fill("pt bMet120-140cut c_minDR",*sigTr,pt);
            else if (met<160) effs.fill("pt bMet140-160cut c_minDR",*sigTr,pt);
            else              effs.fill("pt bMet160cut c_minDR",    *sigTr,pt);
         }
         if (c_pt) effs.fill("sieie bMet100cut", *sigTr, lPho[0]->sigmaIetaIeta);
      } // base met100cut

      if (*Photon175 && pt>240 && c_minDR){
         effs.fill("g135m100 met",Photon135MET100,met);
         effs.fill("g135m100 S",Photon135MET100,MET->sig);
      }

      if (*baseMuTr && c_r9){
         effs.fill("Photon175 bMu",*Photon175,pt);
         effs.fill("Photon165_HE10 bMu",ph165_he10,pt);
         if (c_pt){
            effs.fill("met bMu",*sigTr,met);
            if (c_minDR) effs.fill("met bMu c_minDR",*sigTr,met);
         }
         if (c_met){
            effs.fill("pt bMu",*sigTr,pt);
            if (c_minDR) effs.fill("pt bMu c_minDR",*sigTr,pt);
         }
         if (c_met && c_pt) effs.fill("sieie bMu",*sigTr,lPho[0]->sigmaIetaIeta);

         if (c_met && c_pt && c_minDR){
            iTotal++;
            mRunTotal[*runNo]++;
            if (*sigTr){
               iPass++;
               mRunPass[*runNo]++;
            }
         }
      } // base mu

      if (sim && c_r9 && c_minDR){
         effs.fill("Photon36_R9Id90_HE10_IsoM",*basePhTr,pt);
         if (c_met) effs.fill("Photon36_R9Id90_HE10_IsoM bMet100cut",*basePhTr,pt);

         // efficiency of trigger isolation
         // effs.fill("Iso40ec",c_Iso40ec,pt);
         // effs.fill("Iso40hc",c_Iso40hc,pt);
         // effs.fill("Iso40tr",c_Iso40tr,pt);
      } // sim only
   } // evt loop

   io::RootFileSaver saver("plots.root",TString::Format("trigger/%s",dataset.name.c_str()));
   effs.saveAll(saver,datasetName,sim);
   TCanvas can;
   can.SetRightMargin (.13);
   can.SetBottomMargin(.20);

   TH2* eH2=eDpTdR.CreateHistogram();
   eH2->Draw("colz");
   gfx::setupZaxis(*eH2);
   saver.save(can,"dR_dPT",true,sim);

   if (datasetName!="SinglePhoton" && datasetName!="MET"){
      float eff=iPass/float(iTotal)*100;
      logFile<<datasetName+TString::Format(
         ": ε=%.2f+%.2f-%.2f (Muon base trigger)",
         eff,
         TEfficiency::ClopperPearson(iTotal,iPass,0.682689492137,true)*100-eff,
         eff-TEfficiency::ClopperPearson(iTotal,iPass,0.682689492137,false)*100
         );
   }
   if (datasetName=="SingleMuon"){
      TEfficiency eff("effRun",";Run Number;Trigger Efficiency",mRunTotal.size(),0,1);
      int i=0;
      for (auto pTot: mRunTotal){
         i++;
         eff.SetTotalEvents (i,pTot.second);
         eff.SetPassedEvents(i,mRunPass[pTot.first]);
      }
      TH1F hDummy(*(TH1F*)eff.GetPassedHistogram());
      hDummy.Divide(eff.GetTotalHistogram());
      hDummy.Draw("axis");
      hDummy.GetYaxis()->SetNdivisions(507);
      hDummy.GetXaxis()->SetNdivisions(-1);
      hDummy.GetXaxis()->SetTitleOffset(1.5);
      hDummy.GetXaxis()->CenterLabels();
      i=0;
      std::set<int> strangeRuns={259626,259636,259637,259681,259682,259683,259685};
      for (auto pTot: mRunTotal){
         i++;
         const int iRun=pTot.first;
         TString label=TString::Format("%d",iRun);
         if (strangeRuns.count(iRun)!=0) label="#color[633]{"+label+"}";
         hDummy.GetXaxis()->SetBinLabel(i,label);
      }
      TGraphAsymmErrors *effGr=eff.CreateGraph();
      effGr->Draw("P same");
      saver.save(can,"runNo",true,sim);
   }
}

extern "C"
void run()
{
   process("MET",false);
  // process("SinglePhoton",false);
  // process("GJets",true);
  // process("WGToLNuG",true);
 //  process("ZNuNuGJets",true);
 //  process("T5gg",true);
}
