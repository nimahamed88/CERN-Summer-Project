#ifndef __nanogen_helpers_h__ 
#define __nanogen_helpers_h__ 

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <vector>

using namespace ROOT::VecOps; 

using rvec_f = const RVec<float>;
using rvec_i = const RVec<int>;
using rvec_b = const RVec<bool>;

/**
    @short checks for o.p charge op. flavour lepton pairs compatible with the J/Psi mass
*/
rvec_b dileptonCands(const rvec_i &pdgId, const rvec_f &pt,const rvec_i &eta, const rvec_i &phi,float m0=3.096,float deltam=1.2)
{
    std::vector<bool> isValidCand(pdgId.size(),false);
    
    //loop over the list of leptons
    for(size_t i=0; i<pdgId.size(); i++) {
    
        //try to make an op. charge op. flavour pair
        for(size_t j=i+1; i<pdgId.size(); i++) {
            
            //require same flavour
            if(abs(pdgId[i])!=abs(pdgId[j])) continue;
            if(pdgId[i]*pdgId[j]>0) continue;

            //assign lepton mass
            float mass(-1);
            if(abs(pdgId[i])==11) mass=0.000511;
            if(abs(pdgId[i])==13) mass=0.105658;
            if(mass<0) continue;
            
            //compute mass of the system
            ROOT::Math::PtEtaPhiMVector pi(pt[i], eta[i], phi[i], mass);
            ROOT::Math::PtEtaPhiMVector pj(pt[j], eta[j], phi[j], mass);
            float mll((pi+pj).M());

            //check compatibility with the required mass window
            if(fabs(mll-m0)>deltam) continue;
            isValidCand[i]=true;
            isValidCand[j]=true;        
        }
    }

    return rvec_b(isValidCand.begin(), isValidCand.end());
}


/**
    @short checks if a set of objects are isolated with respect to a reference in the eta-phi plane
*/
rvec_b crossClean(const rvec_f &eta,const rvec_f &phi, const rvec_f &eta_ref, const rvec_f &phi_ref,float cone=0.4)
{
    std::vector<bool> isIso;
    
    //loop over the first list of objects
    for(size_t i=0; i<eta.size(); i++) {

        float minDR(9999.);
        for(size_t j=0; i<eta_ref.size(); i++) {
            minDR = min(minDR,DeltaR(eta[i],eta_ref[j],phi[i],phi_ref[j]));
        }
        isIso.push_back( minDR>cone );
    }

    return rvec_b(isIso.begin(), isIso.end());
}

/**
   @returns a kinematics feature of a two body system
*/
float kinematics2l(const int &pdgId1, const float &pt1, const float &eta1, const float &phi1,
                   const int &pdgId2, const float &pt2, const float &eta2, const float &phi2,
                   std::string kin="mass")
{
    float m1 = abs(pdgId1)==11 ? 0.000511 : 0.105658;
    float m2 = abs(pdgId2)==11 ? 0.000511 : 0.105658;
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, m1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, m2);
    if(kin=="pt") return(p1+p2).Pt();
    if(kin=="eta") return (p1+p2).Eta();
    if(kin=="phi") return (p1+p2).Phi();
    return (p1+p2).M();
}


#endif
