#ifndef _ParticleTrackerAnalyser_h_
#define _ParticleTrackerAnalyser_h_

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

const double _EPSILON_ = 0.0001;

struct EventTruth 
{
   double Weight = 1.0;

   std::vector<bool> inActiveTPC;

   std::vector<int> nHyperons;
   std::vector<int> nKaonShorts;
   std::vector<int> nPions;

   std::vector<double> PrimaryVerticesX;
   std::vector<double> PrimaryVerticesY;
   std::vector<double> PrimaryVerticesZ;

   std::vector<bool> hasKaonShort;
   std::vector<bool> hasKaonShortCharged;

   std::vector<SimParticle> PrimaryLeptons;
   std::vector<SimParticle> PrimaryPions;
   std::vector<SimParticle> PrimaryHyperons;

   std::vector<SimParticle> PrimaryKaonShorts; 
   std::vector<SimParticle> KaonShortChargedDaughters;

   bool eventHasKaonShort;
   bool eventHasHyperon;
   bool eventHasNeutronScatter;
};

class ParticleTrackerAnalyser 
{
   public:

      ParticleTrackerAnalyser(art::Event const& e, std::string genlabel, std::string g4label, bool particlegunmode = false);
      ParticleTrackerAnalyser(art::Event const& e, fhicl::ParameterSet pset, bool particlegunmode = false);

      void GetParticleLists();
      EventTruth GetEventTruth();

      void GetPrimaryParticles();
      void ConstructKaonShorts();
      bool FindNeutronScatter();

      void MCTruthMatch(SimParticle &MCParticle);
      void MCTruthMatch(SimParticle &MCParticle, int TrackID);
   
      void SetFlags();
      
      void SetNeutronScatterThresholds(double ProtonThreshold, double PionThreshold);
      void SetDecayThresholds(double ProtonThreshold, double PionThreshold);

      int GetOrigin(int TrackIdentifier);

   private:
      
      art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
      std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;

      art::Handle<std::vector<simb::MCParticle>> Handle_G4;
      std::vector<art::Ptr<simb::MCParticle>> Vect_G4;

      std::vector<int> primaryParticles;
      std::vector<int> primaryKaonShorts;
      std::vector<int> primaryKaonShortDaughters;

      std::vector<TVector3> primaryVertices;
      bool ProximityMatch(TVector3 FirstVector, TVector3 SecondVector);

      std::map<int,art::Ptr<simb::MCParticle>> particleMap;

      std::vector<int> GetDaughters(const art::Ptr<simb::MCParticle> &MotherParticle);
      std::vector<int> GetDescendants(const art::Ptr<simb::MCParticle> &MotherParticle);
      
      double neutronScatterProtonThresh = 0.15;
      double neutronScatterPionThresh = 0.05;

      double decayProtonThresh = 0.0;
      double decayPionThresh = 0.0;

      EventTruth eventTruth;

      int nMCTruths;

      const bool particleGunMode;
};

}

#endif
