/*
Extract Run:Lumi:Event for events passing some criteria.
Can then be passed to edmPickEvents.py
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents
*/

#include "Config.hpp"
#include "tools/io.hpp"
#include "tree/TreeParticles.hpp"

#include "tools/physics.hpp"

#include <TChain.h>
#include <TTreeReader.h>

static Config const &cfg=Config::get();

void printEdmFormat(std::vector<std::tuple<UInt_t,UInt_t,ULong64_t> > eventList)
{
   for (auto ev :eventList)
      io::log<<TString::Format("%i:%i:%i",std::get<0>(ev),std::get<1>(ev),(int)std::get<2>(ev));
}

void qcd100jetPt_MET()
{
   Datasubset const &dss=cfg.datasets.getDatasubset("QCD_HT100to200");
   TChain chain(cfg.treeName);
   chain.Add(dss.getPath());

   TTreeReader reader(&chain);
   TTreeReaderValue<std::vector<tree::Jet>> jets(reader, "jets");
   TTreeReaderValue<tree::MET>              MET (reader, "met");
   TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
   TTreeReaderValue<UInt_t>    runNo(reader, "runNo");
   TTreeReaderValue<UInt_t>    lumNo(reader, "lumNo");

   std::vector<std::tuple<UInt_t,UInt_t,ULong64_t> > pickedJET;
   std::vector<std::tuple<UInt_t,UInt_t,ULong64_t> > pickedMET;
   while (reader.Next()){
      bool strangeJet=false;
      for (tree::Jet const &j:*jets){
         if (j.isLoose and j.p.Pt()>400) strangeJet=true;
      }
      if (strangeJet)
         pickedJET.push_back(std::make_tuple(*runNo,*lumNo,*evtNo));
      if (MET->p.Pt() > 400)
         pickedMET.push_back(std::make_tuple(*runNo,*lumNo,*evtNo));
   }
   io::log<<TString::Format("picked %i probematic jet-pt events",(int)pickedJET.size());
   printEdmFormat(pickedJET);
   io::log<<TString::Format("picked %i probematic MET events",(int)pickedMET.size());
   printEdmFormat(pickedMET);
}

/* find ecal spike events */
void showerShape()
{
   Dataset const &ds=cfg.datasets.getDataset("SinglePhoton");
   assert(ds.subsets.size()==1);
   Datasubset const &dss=ds.subsets[0];
   TFile file(dss.getPath(),"read");
   io::log<<dss.getPath();
   TTreeReader reader(cfg.treeName, &file);
   TTreeReaderValue<std::vector<tree::Photon>>   photons  (reader, "photons");
   TTreeReaderValue<tree::MET> MET(reader, "met");
   TTreeReaderValue<bool> trigger_Ph   (reader, "HLT_Photon165_HE10_v");

   TString inSR="SR\n";
   TString inCR="CR\n";
   while (reader.Next()){
      if (!*trigger_Ph) continue;

      std::vector<tree::Photon const *> lPho;
      for (tree::Photon const &ph: *photons){
         if (ph.hasPixelSeed) continue;
         if (fabs(ph.p.Eta())>1.4442) continue;
         lPho.push_back(&ph);
      }

      if (lPho.empty()) continue;
      tree::Photon const &pho=*lPho[0];
      float const phoPt=pho.p.Pt();
      float const met=MET->p.Pt();
      float const MT=phys::M_T(pho,*MET);

      if (phoPt<180) continue;
      bool const SR = ((MET->sig>80) && (MT>300));
      bool const CR = ((MET->sig>30) && (MT>100) && !SR);

      if (pho.sigmaIetaIeta<0.001 || pho.sigmaIetaIeta<0.001) {
         TString line=TString::Format("|ph1 pt %.0f sieie %f sipip %f | met %.0f | phi(met,ph) %f|\n",
                               phoPt,pho.sigmaIetaIeta,pho.sigmaIphiIphi,met,pho.p.DeltaPhi(MET->p));
         if      (SR) inSR+=line;
         else if (CR) inCR+=line;
      }
   }
   io::log<<inSR;
   io::log<<inCR;
}

extern "C"
void run()
{
   // qcd100jetPt_MET();
   showerShape();
}
