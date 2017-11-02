/* Z contribution studies */

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"

#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLine.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>

static Config const &cfg=Config::get();

//useful variables, TO BE exported

float DeltaPhi(float a, float b) {
  float temp = std::abs(a-b);
  if (temp <= TMath::Pi())
    return temp;
  else
    return  2.*TMath::Pi() - temp;
} 

void run(TString name,std::vector<TString> datasetsToCompare)
{
   hist::Histograms<TH1F> hs(cfg.datasets.getDatasubsetNames(datasetsToCompare));
   hs.addHist("loose_selection/HT",";HT;Events / bin"                                        ,{0,1000,1500,2500},{50,100,250});
   hs.addHist("loose_selection/MET"   ,";MET;Events / bin"                                   ,{0,300,500,1000},{20,40,100}); 
   hs.addHist("loose_selection/lPhoton1Pt"  ,";#gamma^{loose} PT;Entries / bin"              ,{0,300,500,1000},{20,40,100}); 
   hs.addHist("loose_selection/lPhoton1Pt_ISR"  ,";#gamma^{loose} PT;Entries / bin"          ,{0,300,500,1000},{20,40,100});
   hs.addHist("loose_selection/lPhoton1Pt_ISR_norm"  ,";#gamma^{loose} PT;Entries / bin"     ,1,150,250);   
   hs.addHist("loose_selection/lPhoton1Pt_weighted"  ,";#gamma^{loose} PT;Entries / bin"     ,{0,300,500,1000},{20,40,100});
   hs.addHist("loose_selection/photon_pt_genWeighted"  ,";#gamma^{loose} PT;Entries / bin"   ,50,0,500);    
   hs.addHist("loose_selection/lPhoton1Pt_ISR_weighted"  ,";#gamma^{loose} PT;Entries / bin" ,{0,300,500,1000},{20,40,100});           
   hs.addHist("loose_selection/jetAllPt" ,";jet PT;EntriesBIN"                  ,50,0,1500);
   hs.addHist("loose_selection/jet1Pt",";jet_{1} PT;EventsBIN"                  ,50,0,1500);   
   hs.addHist("loose_selection/NlPhoton",";N_{#gamma loose};EventsBIN"          ,5 ,-.5, 5-.5); 
   hs.addHist("loose_selection/NlElectron",";N_{e loose};EventsBIN"             ,5 ,-.5, 5-.5);    
   hs.addHist("loose_selection/NtMuon",";N_{#mu tight};EventsBIN"               ,5 ,-.5, 5-.5);                    
   hs.addHist("loose_selection/N_reco_ele_onlyISR_genE",";N_{e loose};EventsBIN"    ,5 ,-.5, 5-.5);  
   hs.addHist("loose_selection/N_reco_mu_onlyISR_genE",";N_{#mu tight};EventsBIN"   ,5 ,-.5, 5-.5);    
   hs.addHist("loose_selection/N_reco_ele_onlyISR_genMu",";N_{e loose};EventsBIN"   ,5 ,-.5, 5-.5);  
   hs.addHist("loose_selection/N_reco_mu_onlyISR_genMu",";N_{#mu tight};EventsBIN"  ,5 ,-.5, 5-.5);     
   hs.addHist("loose_selection/Njets",";N_{#gamma loose};EventsBIN"                 ,10 ,-.5, 10-.5);  
   hs.addHist("loose_selection/DR_ph_lep1", ";#DeltaR(#gamma,l_{1});EventsBIN"      ,64,0,6.4);
   hs.addHist("loose_selection/DR_ph_lep1_oneL", ";#DeltaR(#gamma,l_{1});EventsBIN" ,64,0,6.4);   
   hs.addHist("loose_selection/DR_ph_lep2", ";#DeltaR(#gamma,l_{2});EventsBIN"      ,64,0,6.4);   
   hs.addHist("loose_selection/DR_ph_lep1_onZ", ";#DeltaR(#gamma,l_{1});EventsBIN"  ,64,0,6.4);   
   hs.addHist("loose_selection/DR_ph_lep2_onZ", ";#DeltaR(#gamma,l_{2});EventsBIN"  ,64,0,6.4);
   hs.addHist("loose_selection/Dphi_ph_lep1", ";#Delta#Phi(#gamma,l_{1});EventsBIN" ,64,0,6.4);
   hs.addHist("loose_selection/Dphi_ph_lep1_oneL", ";#Delta#Phi(#gamma,l_{1});EventsBIN" ,64,0,6.4);   
   hs.addHist("loose_selection/Dphi_ph_lep2", ";#Delta#Phi(#gamma,l_{2});EventsBIN" ,64,0,6.4);   
   hs.addHist("loose_selection/Dphi_ph_lep1_onZ", ";#Delta#Phi(#gamma,l_{1});EventsBIN"    ,64,0,6.4);   
   hs.addHist("loose_selection/Dphi_ph_lep2_onZ", ";#Delta#Phi(#gamma,l_{2});EventsBIN"    ,64,0,6.4);         
   hs.addHist("loose_selection/photonLepton_pt" ,";e#gamma_{1}^{loose} PT;Events / bin"    ,60,0,300);  
   hs.addHist("loose_selection/trans_mass_1D",";M_{T, e MET};EventsBIN"                    ,50,0,200); 
   hs.addHist("loose_selection/trans_mass_1D_photon",";M_{T, e MET #gamma};EventsBIN"      ,50,0,200);    
   hs.addHist("loose_selection/phLepMetInvMass",";Invariant mass of photon and lepton(s)"  ,50,0,200);
   hs.addHist("loose_selection/LepInvMass",";Invariant mass of leptons"                    ,30,0,300);      
   hs.addHist("loose_selection/leptoMet"   ,";lepto MET;Events / bin"                      ,{0,400,600,1000},{20,50,100}); 
   hs.addHist("loose_selection/met_genMother"   ,";lepto MET;Events / bin"                 ,{0,400,600,1000},{20,50,100});  
   hs.addHist("loose_selection/fake_met"   ,";fake MET;Events / bin"                       ,{0,400,600,1000},{20,50,100});
   hs.addHist("loose_selection/mother_pt"   ,";PT(mother part.);Events / bin"              ,{0,400,600,1000},{20,50,100});         
   hs.addHist("loose_selection/leptoMetGen"   ,";lepto MET;Events / bin"                   ,{0,400,600,1000},{20,50,100});    
   hs.addHist("loose_selection/leptoMet_ISR"   ,";lepto MET;Events / bin"                  ,{0,400,600,1000},{20,50,100}); 
   hs.addHist("loose_selection/leptoMetStack"   ,";lepto MET;Events / bin"                 ,{0,400,600,1000},{20,50,100});      
   
   hs.addHist("trigger_selection/HT",";HT;Events / bin"                           ,{0,1000,1500,2500},{50,100,250});
   hs.addHist("trigger_selection/MET"   ,";MET;Events / bin"                      ,{0,300,500,1000},{20,40,100});   
   hs.addHist("trigger_selection/lPhoton1Pt"  ,";#gamma^{loose} PT;Entries / bin" ,{0,300,500,1000},{20,40,100});   
   hs.addHist("trigger_selection/jetAllPt" ,";jet PT;EntriesBIN"                  ,50,0,1500);
   hs.addHist("trigger_selection/jet1Pt",";jet_{1} PT;EventsBIN"                  ,50,0,1500);   
   hs.addHist("trigger_selection/NlPhoton",";N_{#gamma loose};EventsBIN"          ,5 ,-.5, 5-.5); 
   hs.addHist("trigger_selection/NlElectron",";N_{#gamma loose};EventsBIN"        ,5 ,-.5, 5-.5);    
   hs.addHist("trigger_selection/NtMuon",";N_{#gamma loose};EventsBIN"            ,5 ,-.5, 5-.5);      
   hs.addHist("trigger_selection/Njets",";N_{#gamma loose};EventsBIN"             ,10 ,-.5, 10-.5); 
   
   hs.addHist("zInv_selection/HT",";HT;Events / bin"                           ,{0,1000,1500,2500},{50,100,250});
   hs.addHist("zInv_selection/MET"   ,";MET;Events / bin"                      ,{0,300,500,1000},{20,40,100});   
   hs.addHist("zInv_selection/lPhoton1Pt"  ,";#gamma^{loose} PT;Entries / bin" ,{0,300,500,1000},{20,40,100});   
   hs.addHist("zInv_selection/jetAllPt" ,";jet PT;EntriesBIN"                  ,50,0,1500);
   hs.addHist("zInv_selection/jet1Pt",";jet_{1} PT;EventsBIN"                  ,50,0,1500);   
   hs.addHist("zInv_selection/NlPhoton",";N_{#gamma loose};EventsBIN"          ,5 ,-.5, 5-.5); 
   hs.addHist("zInv_selection/NlElectron",";N_{#gamma loose};EventsBIN"        ,5 ,-.5, 5-.5);    
   hs.addHist("zInv_selection/NtMuon",";N_{#gamma loose};EventsBIN"            ,5 ,-.5, 5-.5);      
   hs.addHist("zInv_selection/Njets",";N_{#gamma loose};EventsBIN"             ,10 ,-.5, 10-.5);             
   hs.addHist("zInv_selection/leptoMet"   ,";lepto MET;Events / bin"           ,{0,400,600,1000},{20,50,100}); 
   hs.addHist("zInv_selection/leptoMet_ISR"   ,";lepto MET;Events / bin"       ,{0,400,600,1000},{20,50,100}); 
     
            
   hs.addHist("phPt"  ,";#gamma^{loose} PT;Entries / bin"                  ,{0,300,500,1000},{20,40,100});
   hs.addHist("ph1Pt" ,";#gamma_{1}^{loose} PT;Events / bin"               ,{0,300,500,1000},{20,40,100});
   hs.addHist("jetPt" ,";jet PT;EntriesBIN"                                ,50,0,1500);
   hs.addHist("jet1Pt",";jet_{1} PT;EventsBIN"                             ,50,0,1500);
   hs.addHist("METSHT",";METSHT;EventsBIN"                                 ,50,0,50);
   hs.addHist("Ngl"   ,";N_{#gamma loose};EventsBIN"                       ,5 ,-.5, 5-.5);
   hs.addHist("Ngm"   ,";N_{#gamma medium};EventsBIN"                      ,5 ,-.5, 5-.5);
   hs.addHist("Ngt"   ,";N_{#gamma tight};EventsBIN"                       ,5 ,-.5, 5-.5);
   hs.addHist("Njets" ,";N_{jets};EventsBIN"                               ,11,-.5,11-.5);
   
   hist::Histograms<TH2F> h2s(cfg.datasets.getDatasubsetNames(datasetsToCompare));
   h2s.addHist("2D_inv_masses2L",";M_{ee} [GeV];M_{ee#gamma} [GeV];Events / bin" ,100,0,150,80,0,400);
   h2s.addHist("2D_inv_masses2L_ISR",";M_{ll} [GeV];M_{ll#gamma} [GeV];Events / bin" ,100,0,150,80,0,400);   
   h2s.addHist("2D_inv_masses2L_FSR",";M_{ll} [GeV];M_{ll#gamma} [GeV];Events / bin" ,100,0,150,80,0,400); 
   h2s.addHist("2D_inv_masses2L_ISR_W",";M_{ll} [GeV];M_{ll#gamma} [GeV];Events / bin" ,100,0,150,80,0,400);   
   h2s.addHist("2D_inv_masses2L_FSR_W",";M_{ll} [GeV];M_{ll#gamma} [GeV];Events / bin" ,100,0,150,80,0,400); 
   h2s.addHist("2D_inv_masses2L_ZNu",";M_{ll} [GeV];M_{ll#gamma} [GeV];Events / bin" ,100,0,150,80,0,400);        
   h2s.addHist("2D_inv_masses1L",";M_{e MET};M_{e MET#gamma};Events / bin"       ,100,0,150,80,0,400);   
   h2s.addHist("2D_trans_mass",";M_{T, e MET};M_{T, e MET#gamma};Events / bin"   ,100,0,150,80,0,400);      
   h2s.addHist("2D_inv_masses0L",";M_{MET};M_{MET#gamma};Events / bin"           ,100,0,150,80,0,400);


//Initalizing Histograms for one distribution filled by different samples   
   TH1F *h1 = new TH1F("leptoMetStack", "Zll reco lepto MET;EventsBIN", 50, 0, 500);
   TH1F *h2 = new TH1F("leptoMetStack_Wlu", "W reco lepto MET;EventsBIN", 50, 0, 500);
   TH1F *h_leptoMet_gen_Zll = new TH1F("leptoMet_gen", "leptoMet", 100, 0, 500);
   // TH1F *h_leptoMet_gen_Wlnu = new TH1F("leptoMetgen", "leptoMet", 100, 0, 500);
   TH1F *h_met_Znunu = new TH1F("MET reco", ";Znunu reco MET; EventsBIN", 100, 0, 500); 
   TH1F *h_met_genMother = new TH1F("MET genMother", ";Zll gen lepto MET;EventsBIN", 50, 0, 500); 
   TH1F *h_photon_pt_genMother = new TH1F("photon pt genMother", "Zll #gamma^{loose} PT; EventsBIN", 50,0,500);
   TH1F *h_met_genMother_W = new TH1F("MET genMother W", ";W gen lepto MET;EventsBIN", 50, 0, 500); 
   TH1F *h_photon_pt_genMother_W = new TH1F("photon pt genMother W", ";W #gamma^{loose} PT; EventsBIN", 50,0,500);           
   TH1F *h_met_Znunu_gen = new TH1F("MET gen - reco", ";Znunu gen MET; EventsBIN", 100, 0, 500);      
   TH1F *h3 = new TH1F("Photon_pt_trigger_selection_Zll", "Zll #gamma^{loose} PT; EventsBIN", 100, 0, 1000); 
   TH1F *h_photon_pt_trigger_selection_gen = new TH1F("Photon_pt_trigger_selection_Zll_gen", "Zll #gamma^{loose} PT; EventsBIN", 50, 0, 500); 
   TH1F *h_photon_pt_trigger_selection_gen_W = new TH1F("Photon_pt_trigger_selection_W_gen", "W #gamma^{loose} PT; EventsBIN", 50, 0, 500);           
   TH1F *h4 = new TH1F("Photon_pt_trigger_selection_Wlu", "W #gamma^{loose} PT; EventsBIN", 100, 0, 1000);  
   
   
              
   for (auto const &dss: cfg.datasets.getDatasubsets(datasetsToCompare)){
      
      float w=dss.xsec/float(dss.Ngen)*cfg.lumi;
      debug << w;
      
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      hs.setCurrentSample(dss.name);
      h2s.setCurrentSample(dss.name);
      
      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<float> HTgen(reader, "genHt");      
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<tree::MET> MET(reader, "met");      
      TTreeReaderValue<std::vector<tree::Photon>> photons(reader, "photons");
      TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
      TTreeReaderValue<std::vector<tree::Muon>>   muons  (reader, "muons");
      TTreeReaderValue<std::vector<tree::Jet>>    jets   (reader, "jets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles (reader, "genParticles");      
      TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles (reader, "intermediateGenParticles");      

      //only MC, otherwise only weight if !isData
      float counter = 0;
      float counter_onZ = 0;
      int iEv=0;
      while (reader.Next()){
         iEv++;
         counter++;
         
         if (iEv>cfg.processFraction*dss.entries) break;
         hs.setFillWeight(*w_pu * *w_mc);
         h2s.setFillWeight(*w_pu * *w_mc);
         
         //Object IDs
                    
         std::vector<tree::Photon const *> lPho,mPho,tPho;
         for (tree::Photon const &ph: *photons){
            // Pixelseedveto and eta cut FOR ALL photons
            if (ph.hasPixelSeed || fabs(ph.p.Eta())>1.4442) continue;
            lPho.push_back(&ph);
            hs.fill("phPt",ph.p.Pt());
            if (ph.isMedium){
               mPho.push_back(&ph);
               if (ph.isTight){
                  tPho.push_back(&ph);
               }
            }
         }
         
         std::vector<tree::Electron const *> lEle,mEle,tEle;
         for (tree::Electron const &ele: *electrons){
            if (ele.isLoose){           
                lEle.push_back(&ele);
                if (ele.isMedium){
                    mEle.push_back(&ele);
                    if (ele.isTight){
                        tEle.push_back(&ele);
                    }
                }
            }
         } 
             
         std::vector<tree::Muon const *> tMu;
         for (tree::Muon const &mu: *muons){
            if (mu.isTight){         
                tMu.push_back(&mu);
            }
         }
         
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);        
         
         //loose selection         
         if (lPho.size()==0) continue;
         if (lPho[0]->r9<.9) continue;
         
         //define variables and further objects after loose selection      
         float met=MET->p.Pt();     
         const float HT=phys::computeHT(cjets);              
         float lPhoPt=lPho[0]->p.Pt();
         float jet1Pt = 0;         
         if (cjets.size() > 0){ 
            jet1Pt=cjets[0].p.Pt();            
         }
         float leptoMet = 0; 
         float leptoMetGen = 0;
         // float leptoMetGenDaughters = 0;
         // float PtDaughters = 0;
         float fake_met = 0;              
         TVector3 PhotonLepton;
         TVector3 PhotonLepton2;         
         TVector3 Zleptons;   
         TVector3 Daughters;
         TVector3 PhotonDaughters;
         TVector3 NewMet;  
         TVector3 RecoNewMet;
         TVector3 fakeMet;                                         
         float invmass = 0;
         float invmassLep = 0;
         float invmassGen = 0;
         float invmassGenLep = 0;         
         float DR_ph_lep1 = 0;
         float DR_ph_lep2 = 0;         
         float Dphi_ph_lep1 = 0;
         float Dphi_ph_lep2 = 0;
         float M1 = 0;
         float M2 = 0;
  //       float dileptonWeight = 1./0.1568;
  //       float WlnuWeight = 1./1.6144;
  //       float dileptonWeight = 1./0.2735;  //both
         float dileptonWeight = 1./0.1669; //only mu
     //    float dileptonGenWeight = 1./0.5495; //140 cut
 //        float dileptonGenWeight = 1./0.60592;  //150 cut 
         float dileptonGenWeight = 1./0.58262;  //150 cut, without taus                       
         float WlnuWeight = 1./3.29;  
  //       float WlnuGenWeight = 1./5.115;
         float WlnuGenWeight = 1./4.80208; // without taus             
                                 
         float minDR=std::numeric_limits<float>::max();
         for (tree::Jet const &jet: *jets){
            if (jet.isLoose){
               float const dr=jet.p.DeltaR(lPho[0]->p);
               float const dpt=fabs(jet.p.Pt()-lPhoPt)/lPhoPt;
               if (dr<.1 && dpt<.5) continue;
               minDR=std::min(minDR,dr);
            }
         }        
         if (minDR < .4) continue;   
                           
         TString selection = "loose_selection/";         
          
         //standard plots         
         hs.fill(selection+"HT",HT);   
         hs.fill(selection+"MET",met);                 
         hs.fill(selection+"lPhoton1Pt",lPhoPt);
         hs.fillweight(selection+"lPhoton1Pt_weighted",lPhoPt,dileptonWeight);            
         hs.fill(selection+"Njets",cjets.size());  
         hs.fill(selection+"NlPhoton",lPho.size()); 
         hs.fill(selection+"NlElectron",lEle.size());          
         hs.fill(selection+"NtMuon",tMu.size());  
         hs.fill(selection+"jet1Pt",jet1Pt);            
         
         if (cjets.size() > 0){                        
            for (tree::Jet const &jet: cjets){
               hs.fill(selection+"jetAllPt",jet.p.Pt());
            }
         }
         
         //further analysis
         if (lEle.size() > 0) {
            DR_ph_lep1 = lPho[0]->p.DeltaR(lEle[0]->p);
            Dphi_ph_lep1 = lPho[0]->p.DeltaPhi(lEle[0]->p);                      
            hs.fill(selection+"DR_ph_lep1_oneL",DR_ph_lep1);
            hs.fill(selection+"Dphi_ph_lep1_oneL",Dphi_ph_lep1);          
         }
         if (tMu.size() > 0) {
            DR_ph_lep1 = lPho[0]->p.DeltaR(tMu[0]->p);
            Dphi_ph_lep1 = lPho[0]->p.DeltaPhi(tMu[0]->p);                                   
            hs.fill(selection+"DR_ph_lep1_oneL",DR_ph_lep1);
            hs.fill(selection+"Dphi_ph_lep1_oneL",Dphi_ph_lep1);                                   
         }
         
         for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){
            if (Igen.daughters.size() == 2 ) {
               if (Igen.daughters[0].pdgId == 13 || Igen.daughters[0].pdgId == -13 ||Igen.daughters[1].pdgId == 13 || Igen.daughters[1].pdgId == -13){
                  hs.fill(selection+"N_reco_ele_onlyISR_genMu",lEle.size());
                  hs.fill(selection+"N_reco_mu_onlyISR_genMu",tMu.size()); 
               }
               if (Igen.daughters[0].pdgId == 11 || Igen.daughters[0].pdgId == -11 ||Igen.daughters[1].pdgId == 11 || Igen.daughters[1].pdgId == -11){
                  hs.fill(selection+"N_reco_ele_onlyISR_genE",lEle.size());
                  hs.fill(selection+"N_reco_mu_onlyISR_genE",tMu.size()); 
               } 
               
            }
         }                 
         
         if (TString(dss.name).Contains("ZGTo2LG")){
            for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){              
               if ( abs(Igen.daughters[0].pdgId) != 15 ) {
                  M1 = Igen.daughters[0].p.Mag(); 
                  M2 = Igen.daughters[1].p.Mag();                      
                  Daughters = Igen.daughters[0].p + Igen.daughters[1].p;
                  PhotonDaughters = Daughters + lPho[0]->p;  
                  invmassGen = sqrt((lPho[0]->p.Mag() + M1 + M2)*(lPho[0]->p.Mag() + M1 + M2) - PhotonDaughters.Mag()*PhotonDaughters.Mag());
                  invmassGenLep = sqrt((M1 + M2)*(M1 + M2) - Daughters.Mag()*Daughters.Mag());                                                                   
                  if (Igen.daughters.size() == 2 ) {                           
                     // PtDaughters = Daughters.Pt();
                     NewMet = Igen.p + MET->p;
                     // leptoMetGenDaughters = met + PtDaughters;
                     leptoMetGen = NewMet.Pt();    
                           
                     hs.fillweight(selection+"photon_pt_genWeighted",lPhoPt,dileptonGenWeight);
                     hs.fillweight(selection+"fake_met",met,dileptonGenWeight);                     
                     hs.fillweight(selection+"mother_pt",Igen.p.Pt(),dileptonGenWeight);
                     h2s.fill("2D_inv_masses2L_ISR",invmassGenLep,invmassGen);                    
                     
                     if (lPhoPt < 150){                      
                        hs.fill(selection+"met_genMother", leptoMetGen);               
                        h_met_genMother->Fill(leptoMetGen, dileptonGenWeight* *w_mc * *w_pu * w); 
                        h_photon_pt_genMother->Fill(lPhoPt, dileptonGenWeight* *w_mc * *w_pu * w); 
                        if (leptoMetGen > 100 ) {
                           h_photon_pt_trigger_selection_gen->Fill(lPhoPt,dileptonGenWeight * *w_mc * *w_pu * w); 
                        }
                     }
                     else if ( lPhoPt > 150 && lPhoPt < 250){
                      hs.fill(selection+"lPhoton1Pt_ISR_norm",lPhoPt);
                     }
                  }
                  else if (Igen.daughters.size() == 3 ) {
                     h2s.fill("2D_inv_masses2L_FSR",invmassLep,invmass);
                  }                  
               }
            }
            
            /*
            if (lEle.size() > 1 ) {
               leptoMet = met + lEle[0]->p.Pt() + lEle[1]->p.Pt();
               PhotonLepton = lPho[0]->p + lEle[0]->p + lEle[1]->p;
               Zleptons = lEle[0]->p + lEle[1]->p;
               invmass = sqrt((lPho[0]->p.Mag() + lEle[0]->p.Mag() + lEle[1]->p.Mag())*(lPho[0]->p.Mag() + lEle[0]->p.Mag() + lEle[1]->p.Mag()) - PhotonLepton.Mag()*PhotonLepton.Mag());
               invmassLep = sqrt((lEle[0]->p.Mag() + lEle[1]->p.Mag())*(lEle[0]->p.Mag() + lEle[1]->p.Mag()) - Zleptons.Mag()*Zleptons.Mag());
               DR_ph_lep2 = lPho[0]->p.DeltaR(lEle[1]->p);
               Dphi_ph_lep2 = lPho[0]->p.DeltaPhi(lEle[1]->p); 
               
               hs.fill(selection+"leptoMet", leptoMet);
               for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){
                  if (Igen.daughters.size() == 2 ) {                                               
              //    if( Dphi_ph_lep2 > 0.5) {
                     hs.fill(selection+"leptoMet_ISR", leptoMet);
                     hs.fill(selection+"lPhoton1Pt_ISR",lPhoPt);                                     
                     hs.fillweight(selection+"lPhoton1Pt_ISR_weighted",lPhoPt,dileptonWeight); 
                     if (lPhoPt < 150){                 
                        h1->Fill(leptoMet,dileptonWeight* *w_mc * *w_pu * w);
                        hs.fillweight(selection+"leptoMetStack", leptoMet,dileptonWeight);
                        if (leptoMet > 100){
                           h3->Fill(lPhoPt,dileptonWeight * *w_mc * *w_pu * w); 
                        }
                     }                  
                  }
                  
               }
               hs.fill(selection+"phLepMetInvMass", invmass);
               hs.fill(selection+"LepInvMass", invmassLep);               
               hs.fill(selection+"photonLepton_pt",lPho[0]->p.Pt());           
               for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){ 
                  if (Igen.daughters.size() == 3 ) {
                     h2s.fill("2D_inv_masses2L",invmassLep,invmass);
                  }
               }
                       
               hs.fill(selection+"DR_ph_lep1",DR_ph_lep1);
               hs.fill(selection+"DR_ph_lep2",DR_ph_lep2);
               hs.fill(selection+"Dphi_ph_lep1",Dphi_ph_lep1);
               hs.fill(selection+"Dphi_ph_lep2",Dphi_ph_lep2);               
               
               if (invmass > 80 && invmass < 100) {
                  counter_onZ++;
                  hs.fill(selection+"DR_ph_lep1_onZ", DR_ph_lep1);
                  hs.fill(selection+"DR_ph_lep2_onZ", DR_ph_lep2); 
                  hs.fill(selection+"Dphi_ph_lep1_onZ", Dphi_ph_lep1);
                  hs.fill(selection+"Dphi_ph_lep2_onZ", Dphi_ph_lep2);                                   
               } 
               
            } 
            */
            if (tMu.size() > 1 ) {
            //   leptoMet = met + tMu[0]->p.Pt() + tMu[1]->p.Pt();  
               RecoNewMet = MET->p + tMu[0]->p + tMu[1]->p;
               leptoMet = RecoNewMet.Pt();
               PhotonLepton = lPho[0]->p + tMu[0]->p + tMu[1]->p;
               Zleptons = tMu[0]->p + tMu[1]->p;
               invmass = sqrt((lPho[0]->p.Mag() + tMu[0]->p.Mag() + tMu[1]->p.Mag())*(lPho[0]->p.Mag() + tMu[0]->p.Mag() + tMu[1]->p.Mag()) - PhotonLepton.Mag()*PhotonLepton.Mag());
               invmassLep = sqrt((tMu[0]->p.Mag() + tMu[1]->p.Mag())*(tMu[0]->p.Mag() + tMu[1]->p.Mag()) - Zleptons.Mag()*Zleptons.Mag());
               DR_ph_lep2 = lPho[0]->p.DeltaR(tMu[1]->p);
               Dphi_ph_lep2 = lPho[0]->p.DeltaPhi(tMu[1]->p);                
                              
               hs.fill(selection+"leptoMet", leptoMet);
               hs.fill(selection+"leptoMetGen", leptoMetGen);               
               h2s.fill("2D_inv_masses2L",invmassLep,invmass);               
       //        if( Dphi_ph_lep2 > 0.5) {               
               for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){
                  if (Igen.daughters.size() == 2 ) {                                    
                     hs.fill(selection+"leptoMet_ISR", leptoMet);
                     hs.fill(selection+"lPhoton1Pt_ISR",lPhoPt);                                       
                     hs.fillweight(selection+"lPhoton1Pt_ISR_weighted",lPhoPt,dileptonWeight);                 
                     if (lPhoPt < 150){                 
                        h1->Fill(leptoMet,dileptonWeight* *w_mc * *w_pu * w);  
                        h_leptoMet_gen_Zll->Fill(leptoMetGen,dileptonWeight* *w_mc * *w_pu * w);                                          
                        hs.fillweight(selection+"leptoMetStack", leptoMet,dileptonWeight);
                        if (leptoMet > 100){
                           h3->Fill(lPhoPt,dileptonWeight * *w_mc * *w_pu * w); 
                        } 
                     }
                     h2s.fill("2D_inv_masses2L_ISR",invmassLep,invmass);                    
                  }
                               
               }
                           
               hs.fill(selection+"phLepMetInvMass", invmass);
               hs.fill(selection+"LepInvMass", invmassLep);               
               hs.fill(selection+"photonLepton_pt",lPho[0]->p.Pt());  
               for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){ 
                  if (Igen.daughters.size() == 3 ) {
                     h2s.fill("2D_inv_masses2L_FSR",invmassLep,invmass);
                  }
               }
                                                               
               hs.fill(selection+"DR_ph_lep1",DR_ph_lep1);
               hs.fill(selection+"DR_ph_lep2",DR_ph_lep2);
               hs.fill(selection+"Dphi_ph_lep1",Dphi_ph_lep1);
               hs.fill(selection+"Dphi_ph_lep2",Dphi_ph_lep2);                
               
               if (invmass > 80 && invmass < 100) {
                  counter_onZ++;
                  hs.fill(selection+"DR_ph_lep1_onZ", DR_ph_lep1);
                  hs.fill(selection+"DR_ph_lep2_onZ", DR_ph_lep2);
                  hs.fill(selection+"Dphi_ph_lep1_onZ", Dphi_ph_lep1);
                  hs.fill(selection+"Dphi_ph_lep2_onZ", Dphi_ph_lep2);                                    
               }               
            }           
            
         }
        
         else if (TString(dss.name).Contains("WGToLNuG")){            
            for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){
               if ( abs(Igen.daughters[0].pdgId) != 15 ) {
                  M1 = Igen.daughters[0].p.Mag(); 
                  M2 = Igen.daughters[1].p.Mag();                      
                  Daughters = Igen.daughters[0].p + Igen.daughters[1].p;
                  PhotonDaughters = Daughters + lPho[0]->p;  
                  invmassGen = sqrt((lPho[0]->p.Mag() + M1 + M2)*(lPho[0]->p.Mag() + M1 + M2) - PhotonDaughters.Mag()*PhotonDaughters.Mag());
                  invmassGenLep = sqrt((M1 + M2)*(M1 + M2) - Daughters.Mag()*Daughters.Mag());
                  if (Igen.daughters.size() == 2 ) {                                    
                     // PtDaughters = Daughters.Pt();
                     NewMet = Igen.daughters[0].p + MET->p;
                     // leptoMetGenDaughters = met + PtDaughters;
                     leptoMetGen = NewMet.Pt();
                     fakeMet = MET->p - Igen.daughters[1].p;
                     fake_met = fakeMet.Pt();                          
                     hs.fillweight(selection+"photon_pt_genWeighted",lPhoPt,WlnuGenWeight);                     
                     hs.fillweight(selection+"fake_met",fake_met,WlnuGenWeight);
                     hs.fillweight(selection+"mother_pt",Igen.p.Pt(),WlnuGenWeight);                     
                     h2s.fill("2D_inv_masses2L_ISR_W",invmassGenLep,invmassGen);    
                     if (lPhoPt < 150){                      
                        hs.fill(selection+"met_genMother", leptoMetGen);               
                        h_met_genMother_W->Fill(leptoMetGen, WlnuGenWeight* *w_mc * *w_pu * w); 
                        h_photon_pt_genMother_W->Fill(lPhoPt, WlnuGenWeight* *w_mc * *w_pu * w); 
                        if (leptoMetGen > 100 ) {
                           h_photon_pt_trigger_selection_gen_W->Fill(lPhoPt,WlnuGenWeight * *w_mc * *w_pu * w); 
                        }
                     }
                     else if ( lPhoPt > 150 && lPhoPt < 250){
                      hs.fill(selection+"lPhoton1Pt_ISR_norm",lPhoPt);
                     }
                  }
                  else if (Igen.daughters.size() == 3 ) { 
                     h2s.fill("2D_inv_masses2L_FSR_W",invmassGenLep,invmassGen);  
                  }                   
               }
            }            
            if (lEle.size() > 0 ) {
          //     leptoMet = met + lEle[0]->p.Pt();
               RecoNewMet = MET->p + lEle[0]->p;
               leptoMet = RecoNewMet.Pt();          
               PhotonLepton = lPho[0]->p + lEle[0]->p + MET->p;
               PhotonLepton2 = lPho[0]->p + lEle[0]->p;               
               Zleptons = lEle[0]->p + MET->p;           
               invmass = sqrt((lPho[0]->p.Mag() + lEle[0]->p.Mag() + met)*(lPho[0]->p.Mag() + lEle[0]->p.Mag() + met) - PhotonLepton.Mag()*PhotonLepton.Mag());
               invmassLep = sqrt((lEle[0]->p.Mag() + met)*(lEle[0]->p.Mag() + met) - Zleptons.Mag()*Zleptons.Mag());
               float trans_mass_combi = TMath::Sqrt(2*PhotonLepton2.Pt()*met*(1.- TMath::Cos(DeltaPhi(MET->p.Phi(),PhotonLepton2.Phi()))));
               float trans_mass = TMath::Sqrt(2*lEle[0]->p.Pt()*met*(1.- TMath::Cos(DeltaPhi(MET->p.Phi(),lEle[0]->p.Phi()))));               

               hs.fill(selection+"leptoMet", leptoMet);               
  //             if( Dphi_ph_lep1 > 0.5) {               
               for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){                  
                  if (Igen.daughters.size() == 2 ) {
                     hs.fill(selection+"leptoMet_ISR", leptoMet); 
                     hs.fill(selection+"lPhoton1Pt_ISR",lPhoPt);                  
                     hs.fillweight(selection+"lPhoton1Pt_ISR_weighted",lPhoPt,WlnuWeight);             
                     if (lPhoPt < 150){                       
                        h2->Fill(leptoMet,WlnuWeight * *w_mc * *w_pu * w);                                   
                        hs.fillweight(selection+"leptoMetStack", leptoMet,WlnuWeight);
                        if (leptoMet > 100){
                           h4->Fill(lPhoPt,WlnuWeight * *w_mc * *w_pu * w); 
                        } 
                     }                                                  
                  }
               }
               
               hs.fill(selection+"phLepMetInvMass", invmass);
               hs.fill(selection+"LepInvMass", invmassLep); 
               hs.fill(selection+"trans_mass_1D", trans_mass);   
               hs.fill(selection+"trans_mass_1D_photon", trans_mass_combi);      
               hs.fill(selection+"photonLepton_pt",PhotonLepton2.Pt());                     
               h2s.fill("2D_inv_masses1L",invmassLep,invmass);
               h2s.fill("2D_trans_mass",trans_mass,trans_mass_combi);                                          
            }
            else if (tMu.size() > 0 ) {
        //       leptoMet = met + tMu[0]->p.Pt();               
               RecoNewMet = MET->p + tMu[0]->p;
               leptoMet = RecoNewMet.Pt();  
               PhotonLepton = lPho[0]->p + tMu[0]->p + MET->p;
               PhotonLepton2 = lPho[0]->p + tMu[0]->p;               
               Zleptons = tMu[0]->p + MET->p;           
               invmass = sqrt((lPho[0]->p.Mag() + tMu[0]->p.Mag() + met)*(lPho[0]->p.Mag() + tMu[0]->p.Mag() + met) - PhotonLepton.Mag()*PhotonLepton.Mag());
               invmassLep = sqrt((tMu[0]->p.Mag() + met)*(tMu[0]->p.Mag() + met) - Zleptons.Mag()*Zleptons.Mag());
               float trans_mass_combi = TMath::Sqrt(2*PhotonLepton2.Pt()*met*(1.- TMath::Cos(DeltaPhi(MET->p.Phi(),PhotonLepton2.Phi()))));
               float trans_mass = TMath::Sqrt(2*tMu[0]->p.Pt()*met*(1.- TMath::Cos(DeltaPhi(MET->p.Phi(),tMu[0]->p.Phi()))));               

               hs.fill(selection+"leptoMet", leptoMet);               
  //             if( Dphi_ph_lep1 > 0.5) {               
               for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){
                  if (Igen.daughters.size() == 2 ) {
                     hs.fill(selection+"leptoMet_ISR", leptoMet); 
                     hs.fill(selection+"lPhoton1Pt_ISR",lPhoPt);
                     hs.fillweight(selection+"lPhoton1Pt_ISR_weighted",lPhoPt,WlnuWeight);               
                     if (lPhoPt < 150){                       
                        h2->Fill(leptoMet,WlnuWeight* *w_mc * *w_pu * w);                                   
                        hs.fillweight(selection+"leptoMetStack", leptoMet,WlnuWeight);
                        if (leptoMet > 100){
                           h4->Fill(lPhoPt,WlnuWeight * *w_mc * *w_pu * w); 
                        } 
                     }                                                  
                  }
               }
               hs.fill(selection+"phLepMetInvMass", invmass);
               hs.fill(selection+"LepInvMass", invmassLep); 
               hs.fill(selection+"trans_mass_1D", trans_mass);   
               hs.fill(selection+"trans_mass_1D_photon", trans_mass_combi);      
               hs.fill(selection+"photonLepton_pt",PhotonLepton2.Pt());                     
               h2s.fill("2D_inv_masses1L",invmassLep,invmass);
               h2s.fill("2D_trans_mass",trans_mass,trans_mass_combi);                                          
            }            
            
         }        
                
         
         if (TString(dss.name).Contains("ZNuNuGJets")){    
            
            for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){                 
               M1 = Igen.daughters[0].p.Mag(); 
               M2 = Igen.daughters[1].p.Mag();                      
               Daughters = Igen.daughters[0].p + Igen.daughters[1].p;
               PhotonDaughters = Daughters + lPho[0]->p;  
               invmassGen = sqrt((lPho[0]->p.Mag() + M1 + M2)*(lPho[0]->p.Mag() + M1 + M2) - PhotonDaughters.Mag()*PhotonDaughters.Mag());
               invmassGenLep = sqrt((M1 + M2)*(M1 + M2) - Daughters.Mag()*Daughters.Mag());
               h_met_Znunu->Fill(met, *w_mc * *w_pu * w);
               h_met_Znunu_gen->Fill(Igen.p.Pt(), *w_mc * *w_pu * w);
               fakeMet = MET->p - Igen.p;
               fake_met = fakeMet.Pt();
               hs.fill(selection+"photon_pt_genWeighted",lPhoPt);
               hs.fill(selection+"fake_met",fake_met);                     
               hs.fill(selection+"mother_pt",Igen.p.Pt());                     
               h2s.fill("2D_inv_masses2L_ZNu",invmassGenLep,invmassGen);
               if (lPhoPt > 150 && lPhoPt < 250){                 
                  hs.fill(selection+"lPhoton1Pt_ISR_norm",lPhoPt);
               }                 
            }                   
            hs.fill(selection+"leptoMet", met);                   
            hs.fill(selection+"leptoMet_ISR", met);                  
            hs.fill(selection+"lPhoton1Pt_ISR",lPhoPt);                   
            hs.fill(selection+"lPhoton1Pt_ISR_weighted",lPhoPt);                  
            if (lPhoPt >= 150){                          
               hs.fill(selection+"met_genMother", met);       
               h1->Fill(met, *w_mc * *w_pu * w);
               h2->Fill(met, *w_mc * *w_pu * w);
               h_met_genMother->Fill(met, *w_mc * *w_pu * w);
               h_photon_pt_genMother->Fill(lPhoPt, *w_mc * *w_pu * w);
               h_met_genMother_W->Fill(met, *w_mc * *w_pu * w);
               h_photon_pt_genMother_W->Fill(lPhoPt, *w_mc * *w_pu * w);                                            
               if (met > 100){
                  h3->Fill(lPhoPt, *w_mc * *w_pu * w);
                  h4->Fill(lPhoPt, *w_mc * *w_pu * w); 
                  h_photon_pt_trigger_selection_gen->Fill(lPhoPt, *w_mc * *w_pu * w); 
                  h_photon_pt_trigger_selection_gen_W->Fill(lPhoPt, *w_mc * *w_pu * w);                                    
               }                                    
               hs.fill(selection+"leptoMetStack", met); 
            }                             
            PhotonLepton = lPho[0]->p + MET->p;
            invmass = sqrt((lPho[0]->p.Mag() + met)*(lPho[0]->p.Mag() + met) - PhotonLepton.Mag()*PhotonLepton.Mag());            
            invmassLep = sqrt(met*met - MET->p.Mag()*MET->p.Mag());
            hs.fill(selection+"phLepMetInvMass", invmass); 
            hs.fill(selection+"LepInvMass", invmass);  
            h2s.fill("2D_inv_masses0L",invmassLep,invmass);                             
         }
        
                                       
         //trigger selection
         if ( (lPhoPt > 50) && (met>100) ){
            selection = "trigger_selection/";
            
            hs.fill(selection+"HT",HT);   
            hs.fill(selection+"MET",met);                 
            hs.fill(selection+"lPhoton1Pt",lPhoPt);   
            hs.fill(selection+"Njets",cjets.size());  
            hs.fill(selection+"NlPhoton",lPho.size()); 
            hs.fill(selection+"NlElectron",lEle.size());          
            hs.fill(selection+"NtMuon",tMu.size());           
            hs.fill(selection+"jet1Pt",jet1Pt);          
            
            if (cjets.size() > 0){                        
               for (tree::Jet const &jet: cjets){
                  hs.fill(selection+"jetAllPt",jet.p.Pt());
               }
            }
         }         
         
         //zInv selection
         if (lPhoPt > 150) {      

            selection = "zInv_selection/";
            
            hs.fill(selection+"HT",HT);   
            hs.fill(selection+"MET",met);                 
            hs.fill(selection+"lPhoton1Pt",lPhoPt);   
            hs.fill(selection+"Njets",cjets.size());  
            hs.fill(selection+"NlPhoton",lPho.size()); 
            hs.fill(selection+"NlElectron",lEle.size());          
            hs.fill(selection+"NtMuon",tMu.size());  
            hs.fill(selection+"jet1Pt",jet1Pt);          
            if (cjets.size() > 0){                        
               for (tree::Jet const &jet: cjets){
                  hs.fill(selection+"jetAllPt",jet.p.Pt());
               }
            }            
            for (tree::IntermediateGenParticle const &Igen: *intermediateGenParticles){
               if (Igen.daughters.size() == 2 ) {
                  if ( abs(Igen.daughters[0].pdgId) != 15 ) {                   
                     if (TString(dss.name).Contains("ZGTo2LG")){                        
                        hs.fillweight(selection+"leptoMet", leptoMet,dileptonWeight);              
                        hs.fillweight(selection+"leptoMet_ISR", leptoMetGen,dileptonGenWeight);
                     }
                     else if (TString(dss.name).Contains("WGToLNuG")){                    
                        if ( abs(Igen.daughters[0].pdgId) != 15 ) {                      
                           hs.fillweight(selection+"leptoMet", leptoMet,WlnuWeight);              
                           hs.fillweight(selection+"leptoMet_ISR", leptoMetGen,WlnuGenWeight);    
                        }              
                     }
                  }                        
                  if (TString(dss.name).Contains("ZNuNuGJets")){                        
                     hs.fill(selection+"leptoMet", met);              
                     hs.fill(selection+"leptoMet_ISR", met);
                  }   
               }
            }         
         }
                 
         hs.fill("Ngl",lPho.size());
         hs.fill("Ngm",mPho.size());
         hs.fill("Ngt",tPho.size());
                  
         for (tree::Jet const &jet: cjets){
            hs.fill("jetPt",jet.p.Pt());
         }
         if (cjets.size()>0) hs.fill("jet1Pt",cjets[0].p.Pt());
         hs.fill("Njets",cjets.size());
         hs.fill("METSHT", MET->p.Pt()/TMath::Sqrt(HT));
      } // evt loop
      
      float ratio = 0;
      if (counter > 0 && counter_onZ > 0) {
         ratio = counter_onZ/counter;
      }
      std::cout << ratio << std::endl;
      
      hs.scaleLumi();
      h2s.scaleLumi();      
      hs.mergeOverflow();
               
      file.Close();
   } // datasets

   hs.combineFromSubsamples(datasetsToCompare);
   h2s.combineFromSubsamples(datasetsToCompare);  

   io::RootFileSaver saver("plots.root",TString::Format("zInvisible%.1f/%s",cfg.processFraction*100,name.Data()));
   gfx::SplitCan spcan;            
   TCanvas can;            
   //plot loose selection
   for (TString sVar: {"photon_pt_genWeighted","HT", "MET","met_genMother","fake_met","mother_pt","lPhoton1Pt","lPhoton1Pt_ISR","lPhoton1Pt_ISR_norm", "lPhoton1Pt_weighted","lPhoton1Pt_ISR_weighted", "Njets", "NlPhoton", "NlElectron", "NtMuon", "jet1Pt", "jetAllPt", 
      "DR_ph_lep1", "DR_ph_lep1_oneL","DR_ph_lep2","DR_ph_lep1_onZ","DR_ph_lep2_onZ","Dphi_ph_lep1","Dphi_ph_lep1_oneL", "Dphi_ph_lep2","Dphi_ph_lep1_onZ","Dphi_ph_lep2_onZ",
      "photonLepton_pt","trans_mass_1D","trans_mass_1D_photon", "phLepMetInvMass","LepInvMass","leptoMet","leptoMetGen", "leptoMet_ISR", "leptoMetStack",
      "N_reco_ele_onlyISR_genMu","N_reco_mu_onlyISR_genMu","N_reco_ele_onlyISR_genE","N_reco_mu_onlyISR_genE"}){
      sVar = "loose_selection/"+sVar;

      spcan.cdUp();
      gPad->SetLogy();
      auto hists=hs.getHistograms(sVar,datasetsToCompare);
      TH1F &hFrame=*(TH1F*)hists[0]->DrawClone("axis");
      Color::reset();
      double max=0;
      for (auto &h: hists){
     //    integral = h->Integral();
     //    h->Scale(scale);
         h->Draw("same hist e");
         h->SetLineColor(Color::next());
         h->SetMarkerSize(0);
         max=std::max(max,h->GetMaximum());        
      }
      hFrame.SetMaximum(max*1.5);
      hFrame.SetMinimum(0.05);      
      gfx::LegendEntries legE = hs.getLegendEntries();
      TLegend leg=legE.buildLegend(.7,.75);
      leg.Draw();

      spcan.cdLow();
      TH1F histRatioFrame=hist::getRatio(*hists[0],*hists[0],"ratio",hist::ONLY1);
      histRatioFrame.Draw("axis");
      histRatioFrame.Draw("same axig");
      max=0;
      for (auto const &h: hists){
         TH1F histRatio=hist::getRatio(*h,*hists[0],"ratio",hist::ONLY1);
         histRatio.DrawClone("same hist e");
         max=std::max(max,histRatio.GetMaximum());
      }
      histRatioFrame.SetMinimum(0);
      histRatioFrame.SetMaximum(10);
      saver.save(spcan,sVar);   
   }
   //plot trigger selection
   for (TString sVar: {"HT", "MET", "lPhoton1Pt", "Njets", "NlPhoton", "NlElectron", "NtMuon", "jet1Pt", "jetAllPt"}){
      sVar = "trigger_selection/"+sVar;      
      spcan.cdUp();
      gPad->SetLogy();      
      auto hists=hs.getHistograms(sVar,datasetsToCompare);
      TH1F &hFrame=*(TH1F*)hists[0]->DrawClone("axis");
      Color::reset();
      double max=0;
      for (auto &h: hists){
     //    integral = h->Integral();
     //    h->Scale(scale);
         h->Draw("same hist e");
         h->SetLineColor(Color::next());
         h->SetMarkerSize(0);
         max=std::max(max,h->GetMaximum());        
      }
      hFrame.SetMaximum(max*1.5);
      hFrame.SetMinimum(0.05);      
      gfx::LegendEntries legE = hs.getLegendEntries();
      TLegend leg=legE.buildLegend(.7,.75);
      leg.Draw();

      spcan.cdLow();
      TH1F histRatioFrame=hist::getRatio(*hists[0],*hists[0],"ratio",hist::ONLY1);
      histRatioFrame.Draw("axis");
      histRatioFrame.Draw("same axig");
      max=0;
      for (auto const &h: hists){
         TH1F histRatio=hist::getRatio(*h,*hists[0],"ratio",hist::ONLY1);
         histRatio.DrawClone("same hist e");
         max=std::max(max,histRatio.GetMaximum());
      }
      histRatioFrame.SetMinimum(0);
      histRatioFrame.SetMaximum(10);
      saver.save(spcan,sVar);   
   }   
   //plot zInv selection
   for (TString sVar: {"HT", "MET", "lPhoton1Pt", "Njets", "NlPhoton", "NlElectron", "NtMuon", "jet1Pt", "jetAllPt","leptoMet", "leptoMet_ISR"}){
      sVar = "zInv_selection/"+sVar;     
      spcan.cdUp();
      gPad->SetLogy();
      auto hists=hs.getHistograms(sVar,datasetsToCompare);
      TH1F &hFrame=*(TH1F*)hists[0]->DrawClone("axis");
      Color::reset();
      double max=0;
      for (auto &h: hists){
     //    integral = h->Integral();
     //    h->Scale(scale);
         h->Draw("same hist e");
         h->SetLineColor(Color::next());
         h->SetMarkerSize(0);
         max=std::max(max,h->GetMaximum());        
      }
      hFrame.SetMaximum(max*1.5);
      hFrame.SetMinimum(0.05);      
      gfx::LegendEntries legE = hs.getLegendEntries();
      TLegend leg=legE.buildLegend(.7,.75);
      leg.Draw();

      spcan.cdLow();
      TH1F histRatioFrame=hist::getRatio(*hists[0],*hists[0],"ratio",hist::ONLY1);
      histRatioFrame.Draw("axis");
      histRatioFrame.Draw("same axig");
      max=0;
      for (auto const &h: hists){
         TH1F histRatio=hist::getRatio(*h,*hists[0],"ratio",hist::ONLY1);
         histRatio.DrawClone("same hist e");
         max=std::max(max,histRatio.GetMaximum());
      }
      histRatioFrame.SetMinimum(0);
      histRatioFrame.SetMaximum(5);
      saver.save(spcan,sVar);   
   }    
   //plot rest
   for (TString sVar: {"ph1Pt","phPt","Ngl","Ngm","Ngt","Njets"}){   
      spcan.cdUp();
      gPad->SetLogy();
      auto hists=hs.getHistograms(sVar,datasetsToCompare);
      TH1F &hFrame=*(TH1F*)hists[0]->DrawClone("axis");
      Color::reset();
      double max=0;
      for (auto &h: hists){
     //    integral = h->Integral();
     //    h->Scale(scale);
         h->Draw("same hist e");
         h->SetLineColor(Color::next());
         h->SetMarkerSize(0);
         max=std::max(max,h->GetMaximum());        
      }
      hFrame.SetMaximum(max*1.5);
      hFrame.SetMinimum(0.05);      
      gfx::LegendEntries legE = hs.getLegendEntries();
      TLegend leg=legE.buildLegend(.7,.75);
      leg.Draw();

      spcan.cdLow();
      TH1F histRatioFrame=hist::getRatio(*hists[0],*hists[0],"ratio",hist::ONLY1);
      histRatioFrame.Draw("axis");
      histRatioFrame.Draw("same axig");
      max=0;
      for (auto const &h: hists){
         TH1F histRatio=hist::getRatio(*h,*hists[0],"ratio",hist::ONLY1);
         histRatio.DrawClone("same hist e");
         max=std::max(max,histRatio.GetMaximum());
      }
      histRatioFrame.SetMinimum(0);
      histRatioFrame.SetMaximum(10);
      saver.save(spcan,sVar);   
   }   
   ////////
   // 2d //
   ////////
      
   can.cd();
   can.SetRightMargin (.14);
   can.SetBottomMargin(.14);
   can.SetLogy(false);
   can.SetLogz();
   h1->Draw("hist e");
   saver.save(can,"leptoMetStack");   
   h_leptoMet_gen_Zll->Draw("hist e");
   saver.save(can,"leptoMet_gen_Zll");   
   can.Clear();
   h_met_Znunu_gen->SetLineColor(kRed);
   h_met_Znunu_gen->Draw("hist e");   
   h_met_Znunu->Draw("hist same e");   
   saver.save(can,"Met_reco_gen_Znunu");   
   can.Clear();  
   h_met_genMother->Draw("hist e");
   saver.save(can,"met_genMother");
   can.Clear();    
   h_photon_pt_genMother->Draw("hist e");
   saver.save(can,"photon_pt_genMother"); 
   h_met_genMother_W->Draw("hist e");
   saver.save(can,"met_genMother_W");
   can.Clear();    
   h_photon_pt_genMother_W->Draw("hist e");
   saver.save(can,"photon_pt_genMother_W");   
   can.Clear();
   h_met_genMother->Draw("hist e");
   h_met_genMother_W->SetLineColor(kRed);      
   h_met_genMother_W->Draw("hist same e");  
   saver.save(can,"full_met_both"); 
   can.Clear();            
   h2->Draw("hist e");
   saver.save(can,"leptoMetStack_Wlu");   
   can.Clear();   
   h3->Draw("hist e");
   saver.save(can,"Photon_pt_trigger_selection_Zll");   
   can.Clear(); 
   h_photon_pt_trigger_selection_gen->Draw("hist e");
   saver.save(can,"Photon_pt_trigger_selection_Zll_gen"); 
   h_photon_pt_trigger_selection_gen_W->Draw("hist e");
   saver.save(can,"Photon_pt_trigger_selection_W_gen");        
   can.Clear(); 
   h4->Draw("hist e");
   saver.save(can,"Photon_pt_trigger_selection_Wlu");   
   can.Clear(); 
   h1->Draw("hist e");
   h2->SetLineColor(kRed);
   h2->Draw("hist e same");  
   saver.save(can,"leptoMet_both");   
   can.Clear();             
   for (auto const &dss: cfg.datasets.getDatasubsets(datasetsToCompare)){   
      if (TString(dss.name).Contains("ZGTo2LG")){   
         for (TString sVar: {"2D_inv_masses2L","2D_inv_masses2L_ISR","2D_inv_masses2L_FSR"}){
            auto hists=h2s.getHistograms(sVar,datasetsToCompare);       
            hists[0]->DrawClone("axis");
            Color::reset();
            for (auto &h: hists){    
               h->Draw("colz");
               can.RedrawAxis();
               saver.save(can,sVar);
            }
         }
      }
      if (TString(dss.name).Contains("WGToLNuG")){   
         for (TString sVar: {"2D_inv_masses1L", "2D_trans_mass","2D_inv_masses2L_ISR_W","2D_inv_masses2L_FSR_W"}){
            auto hists=h2s.getHistograms(sVar,datasetsToCompare);       
            hists[0]->DrawClone("axis");
            Color::reset();
            for (auto &h: hists){    
               h->Draw("colz");
               can.RedrawAxis();
               saver.save(can,sVar);
            } 
         }
      }
         
      for (TString sVar: {"2D_inv_masses0L","2D_inv_masses2L_ZNu"}){
         auto hists=h2s.getHistograms(sVar,datasetsToCompare);       
         hists[0]->DrawClone("axis");
         Color::reset();
         for (auto &h: hists){    
            h->Draw("colz");
            can.RedrawAxis();
            saver.save(can,sVar);
         }  
      }          
   }     
}

  
extern "C"
void run()
{
   run("Zbkg",{"ZNuNuGJets","ZGTo2LG","WGToLNuG"});
 //  run("Z_2l",{"ZGTo2LG"});
 //  run("Z_0l",{"ZNuNuGJets"});
 //  run("Z_1L",{"WGToLNuG"});          
}
