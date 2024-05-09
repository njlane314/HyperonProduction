#ifndef _ReconstructionAnalyser_h_
#define _ReconstructionAnalyser_h_

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Headers/LLRPID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_ProtonMuonLookup.h"

#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/RecoHelper.h"
#include "ubana/HyperonProduction/Objects/SimHelper.h"
#include "ubana/HyperonProduction/Alg/PIDManager.h"
#include "ubana/HyperonProduction/Modules/SubModules/ParticleTrackerAnalyser.h"

#include "TVector3.h"

using std::string;

namespace hyperon {

struct RecoData 
{
   TVector3 RecoPrimaryVertex = TVector3(-1000, -1000, -1000);

   int NPrimaryDaughters; 
   int NPrimaryTrackDaughters;
   int NPrimaryShowerDaughters;

   std::vector<RecoParticle> TrackPrimaryDaughters;
   std::vector<RecoParticle> ShowerPrimaryDaughters;

   std::vector<TVector3> TrackStarts;

   size_t TrueMuonIndex = -1;
   size_t TrueDecayPionPlusIndex = -1;
   size_t TrueDecayPionMinusIndex = -1;

   bool GoodReco = false;
};

class ReconstructionAnalyser {

   public:

      ReconstructionAnalyser(art::Event const& E, 
                        bool IsData, 
                        string PFParticleLabel, 
                        string TrackLabel, 
                        string ShowerLabel, 
                        string VertexLabel, 
                        string PIDLabel, 
                        string CaloLabel, 
                        string HitLabel,
                        string HitTruthAssnLabel, 
                        string TrackHitAssnLabel, 
                        string MetaDataLabel, 
                        string GenLabel,
                        string G4Label, 
                        bool DoGetPIDs, 
                        bool IncludeCosmics, 
                        bool ParticleGunMode = false);

      ReconstructionAnalyser(art::Event const& E, 
                        bool IsData, 
                        fhicl::ParameterSet ParamSet, 
                        bool ParticleGunMode = false);

      void PrepareInfo(); 
      TVector3 GetPrimaryVertex();
      void SetIndices(std::vector<bool> IsSignal);

      RecoData GetInfo();
      void SetResRangeCutoff(double cutoff){ resRangeCutoff = cutoff; }

   private:

      art::Handle<std::vector<recob::PFParticle>> handlePFParticle;
      std::vector<art::Ptr<recob::PFParticle>> vectPFParticle;

      art::Handle<std::vector<recob::Track>> handleTrack;
      std::vector<art::Ptr<recob::Track>> vectTrack;

      art::Handle<std::vector<recob::Shower>> handleShower;
      std::vector<art::Ptr<recob::Shower>> vectShower;

      art::Handle<std::vector<recob::Hit>> handleHit;
      std::vector<art::Ptr<recob::Hit>> vectHit;

      RecoParticle MakeRecoParticle(const art::Ptr<recob::PFParticle> &PFP);

      art::FindManyP<recob::Vertex>* assocPFParticleVertex;
      art::FindManyP<recob::Track>* assocPFParticleTrack;
      art::FindManyP<recob::Shower>* assocPFParticleShower;
      art::FindManyP<larpandoraobj::PFParticleMetadata>* assocPFParticleMetadata;
      art::FindManyP<recob::Hit>* assocTrackHit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* particlesPerHit;
      art::FindManyP<anab::Calorimetry>* assocTrackCalo;
      art::FindManyP<anab::ParticleID>* assocTrackPID;

      ParticleTrackerAnalyser* particleTrackerAnalyser = nullptr;
      PIDManager PIDCalc;      

      RecoData eventReco;
      size_t neutrinoIndex = 99999;
      std::map<size_t, int> particleMap;

      void SetPFPMetaData(const art::Ptr<recob::PFParticle> &PFP, RecoParticle &Particle);
      void GetTrackData(const art::Ptr<recob::PFParticle> &PFP, RecoParticle &Particle);
      void TruthMatch(const art::Ptr<recob::Track> &Trk, RecoParticle &Particle);
      void GetPIDs(const art::Ptr<recob::Track> &Trk, RecoParticle &Particle);
      void GetVertexData(const art::Ptr<recob::PFParticle> &PFP, RecoParticle &Particle);

      bool isData;
      bool doGetPIDs = true;
      double resRangeCutoff = 5; 
      const bool includeCosmics;
      const bool particleGunMode;
};

}

#endif