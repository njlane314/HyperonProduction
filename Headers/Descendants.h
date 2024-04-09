#ifndef _Descendents_h_
#define _Descendents_h_

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace hyperon{

std::unordered_map<int, art::Ptr<recob::PFParticle>> GetPFParticleMap(const std::vector<art::Ptr<recob::PFParticle>> &allPFParticles){
    std::unordered_map<int, art::Ptr<recob::PFParticle>> pfParticleMap;

    for (const auto &pfParticle : allPFParticles)
    {
        if (!pfParticleMap.emplace(pfParticle->Self(), pfParticle).second)
            throw cet::exception("RecoHelper::GetPFParticleMap") << " - Found repeated PFParticle with Self = " << pfParticle->Self() << "." << std::endl;
    }

    return pfParticleMap;
}

std::vector<art::Ptr<recob::PFParticle>> GetDaughters(const art::Ptr<recob::PFParticle> &particle, const std::unordered_map<int, art::Ptr<recob::PFParticle>> &pfParticleMap){
   std::vector<art::Ptr<recob::PFParticle>> daughters;

   for (int i = 0; i < particle->NumDaughters(); ++i)
   {
      const auto daughterIter = pfParticleMap.find(particle->Daughter(i));
      if (daughterIter == pfParticleMap.end())
         throw cet::exception("RecoHelper::GetDaughter") << " - Couldn't find daughter PFParticle in hierarchy." << std::endl;

      daughters.push_back(daughterIter->second);
   }

   return daughters;
}

void GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const std::unordered_map<int, art::Ptr<recob::PFParticle>> &pfParticleMap, std::vector<art::Ptr<recob::PFParticle>> &downstreamParticles){

   downstreamParticles.push_back(particle);
   for (const auto &daughter : GetDaughters(particle, pfParticleMap))
      GetDownstreamParticles(daughter, pfParticleMap, downstreamParticles);
   
}

unsigned int GetNDescendents(const art::Ptr<recob::PFParticle> &particle, const std::unordered_map<int, art::Ptr<recob::PFParticle>> &pfParticleMap){
   
   unsigned int nDescendents = 0;
   std::vector<art::Ptr<recob::PFParticle>> downstreamParticles;
   GetDownstreamParticles(particle, pfParticleMap, downstreamParticles);
   for (const auto &downstreamParticle : downstreamParticles)
      if (downstreamParticle != particle) nDescendents++;
   return nDescendents;

}

}

#endif