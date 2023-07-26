#ifndef _SubModuleG4Truth_h_
#define _SubModuleG4Truth_h_

#include <string>
#include <vector>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "cetlib_except/exception.h"

#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Headers/FV.h"
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"

namespace hyperon {

// Used for comparing vertex positions
const double _EPSILON_ = 0.0001;

struct G4Truth {

   // Flags applying to the entire event 
   // Use for sample orthogonality
   bool EventHasHyperon = false;
   bool EventHasNeutronScatter = false;

   // Flags for each MCTruth
   std::vector<bool> InActiveTPC;
   std::vector<bool> IsHyperon;
   std::vector<bool> IsLambda;
   std::vector<bool> IsLambdaCharged;
   std::vector<bool> IsSigmaZero;
   std::vector<bool> IsSigmaZeroCharged;
   std::vector<bool> IsAssociatedHyperon;

   double Weight = 1.0;

   std::vector<SimParticle> Lepton;
   std::vector<SimParticle> Hyperon;
   std::vector<SimParticle> PrimaryNucleon;
   std::vector<SimParticle> PrimaryPion;
   std::vector<SimParticle> PrimaryKaon;
   std::vector<SimParticle> PrimaryNucleus;
   std::vector<SimParticle> Decay;
   std::vector<SimParticle> KaonDecay;
   std::vector<SimParticle> SigmaZeroDecayPhoton;
   std::vector<SimParticle> SigmaZeroDecayLambda;

   //TVector3 DecayVertex;

   std::vector<double> TruePrimaryVertex_X;
   std::vector<double> TruePrimaryVertex_Y;
   std::vector<double> TruePrimaryVertex_Z;
   
   std::vector<double> DecayVertex_X;
   std::vector<double> DecayVertex_Y;
   std::vector<double> DecayVertex_Z;

};

class SubModuleG4Truth {

   public:

      //SubModuleG4Truth();
      SubModuleG4Truth(art::Event const& e,std::string genlabel,std::string g4label,bool particlegunmode=false);
      SubModuleG4Truth(art::Event const& e,fhicl::ParameterSet pset,bool particlegunmode=false);

      void GetParticleLists();
      G4Truth GetG4Info();

      void GetPrimaryParticles();
      void GetHyperonDecay();
      void GetKaonDecay();
      void GetSigmaZeroDecay();
      bool FindNeutronScatter();
      int GetOrigin(int trackid);
      void MCTruthMatch(SimParticle &P);
      void SetFlags();
      
      void SetNeutronScatterThresholds(double neutronscatterprotonthresh,double neutronscatterpionthresh);
      void SetDecayThresholds(double decayprotonthresh,double decaypionthresh);

   private:
      
      art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
      std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;

      art::Handle<std::vector<simb::MCParticle>> Handle_G4;
      std::vector<art::Ptr<simb::MCParticle>> Vect_G4;

      std::vector<int> Daughter_IDs;        // IDs of Lambda,SigmaP,SigmaM decay products
      std::vector<int> SigmaZero_Daughter_IDs; // IDs of SigmaZero decay products
      std::vector<int> Primary_IDs;         // IDs of particles produced at primary vertex
      std::vector<int> Kaon_Daughter_IDs;   // IDs of Kaon decay products

      std::vector<TVector3> PrimaryVertices;
      bool PosMatch(TVector3 Pos1,TVector3 Pos2);

      std::map<int,art::Ptr<simb::MCParticle>> partByID;

      std::vector<int> GetChildIDs(const art::Ptr<simb::MCParticle> &g4p,bool IsNeutron=false);

      double NeutronScatterProtonThresh = 0.15;
      double NeutronScatterPionThresh = 0.05;

      double DecayProtonThresh = 0.0;
      double DecayPionThresh = 0.0;

      G4Truth theTruth;

      int NMCTruths;
       
      const bool ParticleGunMode;

};

}

#endif
