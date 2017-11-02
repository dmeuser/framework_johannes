#include "physics.hpp"

#include <limits>

std::vector<tree::Jet> phys::getCleanedJets(std::vector<tree::Jet> const &jets)
{
   std::vector<tree::Jet> cjets;
   for (tree::Jet j: jets){
      // vary jet energy scale up or down (resort after loop!)
      // j.p*=(1.-j.uncert);
      if (!j.isLoose || j.p.Pt()<30 || fabs(j.p.Eta())>3.0) continue;
      if (j.hasPhotonMatch || j.hasElectronMatch || j.hasMuonMatch) continue;
      cjets.push_back(j);
   }
   // sort(cjets.begin(), cjets.end(), tree::PtGreater);
   return cjets;
}


/*
 * - compute HT using the given jet collection
 * - jets have to be cleaned of photons/electrons and
 *   quality criteria have to be applied before, if wanted
 * - any kinematic selection (pT, eta) has to be done before
 */
float phys::computeHT(std::vector<tree::Jet> const &jets)
{
   float HT=0.0;
   for (auto const &jet: jets) HT+=jet.p.Pt();
   return HT;
}

float phys::METoverSqrtHT(float MET, float HT)
{
   return (HT==0.0
           ? std::numeric_limits<float>::infinity()
           : MET/TMath::Sqrt(HT));
}

float phys::M_T2(TVector3 const &v1, TVector3 const &v2)
{
   return 2.0*v1.Pt()*v2.Pt()*(1-TMath::Cos(v1.DeltaPhi(v2)));
}

float phys::M_T(TVector3 const &v1, TVector3 const &v2)
{
   return TMath::Sqrt(M_T2(v1,v2));
}

float phys::M_T(tree::Particle const &p1, tree::Particle const &p2)
{
   return M_T(p1.p,p2.p);
}

bool phys::matchesGen(tree::Particle const &p,std::vector<tree::GenParticle> const &genP,int pdgId,float dR,float rel_dp)
{
   for (tree::GenParticle const &gP: genP){
      float const rDP=fabs(p.p.Pt()-gP.p.Pt())/p.p.Pt();
      if (abs(gP.pdgId)==pdgId && gP.p.DeltaR(p.p)<dR && rDP<rel_dp) {
         return true;
      }
   }
   return false;
}
