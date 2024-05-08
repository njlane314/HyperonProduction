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
#include "ubana/HyperonProduction/Headers/LLR_PID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_proton_muon_lookup.h"
#include "ubana/HyperonProduction/Headers/LLR_PID_K.h"
#include "ubana/HyperonProduction/Headers/LLRPID_kaon_proton_lookup.h"
//#include "ubana/HyperonProduction/Headers/Descendants.h"

#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"
#include "ubana/HyperonProduction/Alg/PIDManager.h"
#include "ubana/HyperonProduction/Modules/SubModules/ParticleTrackerAnalyser.h"
#include "ubana/HyperonProduction/Alg/BDTHandle.h"

#include "TVector3.h"

using std::string;

namespace hyperon {

struct RecoData {

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

      ReconstructionAnalyser(art::Event const& e,bool isdata,string pfparticlelabel,string tracklabel,
                        string showerlabel,string vertexlabel,string pidlabel,string calolabel,string hitlabel,
                        string hittruthassnlabel,string trackhitassnlabel,string metadatalabel,string genlabel,
                        string g4label,bool dogetpids,bool includecosmics,bool particlegunmode=false);

      ReconstructionAnalyser(art::Event const& e,bool isdata,fhicl::ParameterSet pset,bool particlegunmode=false);

      void PrepareInfo(); 
      TVector3 GetPrimaryVertex();
      void SetIndices(std::vector<bool> IsSignal);

      RecoData GetInfo();
      void SetResRangeCutoff(double cutoff){ ResRangeCutoff = cutoff; }

   private:

      art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
      std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;

      art::Handle<std::vector<recob::Track>> Handle_Track;
      std::vector<art::Ptr<recob::Track>> Vect_Track;

      art::Handle<std::vector<recob::Shower>> Handle_Shower;
      std::vector<art::Ptr<recob::Shower>> Vect_Shower;

      art::Handle<std::vector<recob::Hit>> Handle_Hit;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit;

      RecoParticle MakeRecoParticle(const art::Ptr<recob::PFParticle> &pfp);

      art::FindManyP<recob::Vertex>* Assoc_PFParticleVertex;
      art::FindManyP<recob::Track>* Assoc_PFParticleTrack;
      art::FindManyP<recob::Shower>* Assoc_PFParticleShower;
      art::FindManyP<larpandoraobj::PFParticleMetadata>* Assoc_PFParticleMetadata;
      art::FindManyP<recob::Hit>* Assoc_TrackHit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* Assoc_MCParticleBacktracker;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo;
      art::FindManyP<anab::ParticleID>* Assoc_TrackPID;

      searchingfornues::LLRPID llr_pid_calculator;
      searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;

      searchingfornuesk::LLRPIDK llr_pid_calculator_kaon;
      searchingfornuesk::KaonProtonLookUpParameters kaonproton_parameters;

      ParticleTrackerAnalyser* G4T = nullptr;
      PIDManager PIDCalc;      

      RecoData eventReco;
      size_t neutrinoID = 99999;
      std::map<size_t,int> m_PFPID_TrackIndex;

      void GetPFPMetadata(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);
      void GetTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);
      void TruthMatch(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void GetVertexData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);

      bool IsData;
      bool DoGetPIDs=true;
      double ResRangeCutoff=5; 
      const bool IncludeCosmics;
      const bool ParticleGunMode;
};

}

#endif