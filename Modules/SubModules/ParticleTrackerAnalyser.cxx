#ifndef _ParticleTrackerAnalyser_cxx_
#define _ParticleTrackerAnalyser_cxx_

#include "ParticleTrackerAnalyser.h"

using namespace hyperon;

ParticleTrackerAnalyser::ParticleTrackerAnalyser(art::Event const& e, std::string genlabel, std::string g4label, bool particlegunmode) :
particleGunMode(particlegunmode)
{
   if(!e.getByLabel(genlabel, Handle_MCTruth))  
      throw cet::exception("ParticleTrackerAnalyser") << "No MCTruth data product!" << std::endl;

   if(!e.getByLabel(g4label, Handle_G4)) 
      throw cet::exception("ParticleTrackerAnalyser") << "No Geant4 data products found!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth, Handle_MCTruth);  
   art::fill_ptr_vector(Vect_G4, Handle_G4);

   for(const art::Ptr<simb::MCParticle> &Particle : Vect_G4){
      particleMap.insert(std::make_pair(Particle->TrackId(), Particle));
   }
}

ParticleTrackerAnalyser::ParticleTrackerAnalyser(art::Event const& e, fhicl::ParameterSet pset, bool particlegunmode) :
ParticleTrackerAnalyser(e,
   pset.get<std::string>("GeneratorModuleLabel", "generator"),
   pset.get<std::string>("G4ModuleLabel", "largeant"),
   particlegunmode)
{
   SetNeutronScatterThresholds(pset.get<double>("NeutronScatterProtonThresh", 0.15), pset.get<double>("NeutronScatterPionThresh", 0.05));
   SetDecayThresholds(pset.get<double>("DecayProtonThresh", 0.0), pset.get<double>("DecayPionThresh", 0.0));
}

void ParticleTrackerAnalyser::GetParticleLists()
{
   primaryParticles.clear();

   primaryKaonShorts.clear();
   primaryKaonShortDaughters.clear();

   for(const art::Ptr<simb::MCParticle> &Particle : Vect_G4){
      if(Particle->Mother() != 0) continue;
      primaryParticles.push_back(Particle->TrackId());

      std::vector<int> Daughters = GetDaughters(Particle);

      int PdgCode = Particle->PdgCode(); 
      if(isNeutralKaon(PdgCode) && Particle->EndProcess() == "Decay"){
         for(size_t d = 0; d < Daughters.size(); d++){
            if(particleMap.find(Daughters[d]) == particleMap.end()) continue;

            art::Ptr<simb::MCParticle> DaughterParticle = particleMap[Daughters.at(d)];
            if(DaughterParticle->PdgCode() != 130 && DaughterParticle->PdgCode() != 310)
               throw std::invalid_argument("Neutral kaon daughter with pdg code " + std::to_string(DaughterParticle->PdgCode()) + " expected 130 or 310!");
            
            if(isKaonShort(DaughterParticle->PdgCode())){
               primaryKaonShorts.push_back(DaughterParticle->TrackId());
            
               if(DaughterParticle->EndProcess() == "Decay"){
                  std::vector<int> KaonShortDaughters = GetDaughters(DaughterParticle);   
                  primaryKaonShortDaughters.insert(primaryKaonShortDaughters.begin(), KaonShortDaughters.begin(), KaonShortDaughters.end());  
               }
            }
         }
      }
   }

   for(const art::Ptr<simb::MCTruth> &MCTruth : Vect_MCTruth){         
      if(!particleGunMode){
         for(int p = 0; p < MCTruth->NParticles(); p++){
            simb::MCParticle Particle = MCTruth->GetParticle(p);
            if(isNeutrino(Particle.PdgCode()) && Particle.StatusCode() == 0){ 
               eventTruth.PrimaryVerticesX.push_back(Particle.Vx());
               eventTruth.PrimaryVerticesY.push_back(Particle.Vy());
               eventTruth.PrimaryVerticesZ.push_back(Particle.Vz());
            }
         }
      }       
      else if(MCTruth->NParticles()){
         eventTruth.PrimaryVerticesX.push_back(MCTruth->GetParticle(0).Vx());
         eventTruth.PrimaryVerticesY.push_back(MCTruth->GetParticle(0).Vy());
         eventTruth.PrimaryVerticesZ.push_back(MCTruth->GetParticle(0).Vz());
      }
      else if(MCTruth->NParticles() == 0)
         throw cet::exception("ParticleTrackerAnalyser") << "MCTruth made with particle gun contains no particles!" << std::endl;
   }

   nMCTruths = Vect_MCTruth.size();

   if(primaryVertices.size() != Vect_MCTruth.size()) 
      throw cet::exception("ParticleTrackerAnalyser") << "Primary vertices and MCTruth vector-size mismatch!" << std::endl;
}

std::vector<int> ParticleTrackerAnalyser::GetDaughters(const art::Ptr<simb::MCParticle> &MotherParticle)
{
   std::vector<int> DecayDaughters;

   for(int i = 0; i < MotherParticle->NumberDaughters(); i++){
      if(particleMap.find(MotherParticle->Daughter(i)) == particleMap.end()) continue;

      art::Ptr<simb::MCParticle> DaughterParticle = particleMap[MotherParticle->Daughter(i)];

      double DaughterStartX = DaughterParticle->Position().X();
      double DaughterStartY = DaughterParticle->Position().Y();
      double DaughterStartZ = DaughterParticle->Position().Z();

      double MotherEndX = MotherParticle->EndPosition().X();
      double MotherEndY = MotherParticle->EndPosition().Y();
      double MotherEndZ = MotherParticle->EndPosition().Z();

      if(!ProximityMatch(TVector3(DaughterStartX, DaughterStartY, DaughterStartZ), TVector3(MotherEndX, MotherEndY, MotherEndZ))) continue;

      DecayDaughters.push_back(MotherParticle->Daughter(i));
   }

   return DecayDaughters;
}

std::vector<int> ParticleTrackerAnalyser::GetDescendants(const art::Ptr<simb::MCParticle> &MotherParticle)
{
   std::vector<int> Descendants;
   std::vector<int> Daughters = GetDaughters(MotherParticle);

   Descendants.insert(Descendants.end(), Daughters.begin(), Daughters.end());

   for(int Daughter : Daughters){
      art::Ptr<simb::MCParticle> DaughterParticle = particleMap.at(Daughter);
      std::vector<int> DaughterDescendants = GetDescendants(DaughterParticle);

      Descendants.insert(Descendants.end(), DaughterDescendants.begin(), DaughterDescendants.end());
   }

   return Descendants;
}

bool ParticleTrackerAnalyser::ProximityMatch(TVector3 FirstVector, TVector3 SecondVector)
{
   return (FirstVector - SecondVector).Mag() < _EPSILON_;
}

EventTruth ParticleTrackerAnalyser::GetEventTruth()
{
   GetParticleLists();

   eventTruth.inActiveTPC.resize(nMCTruths);

   eventTruth.PrimaryVerticesX.resize(nMCTruths);
   eventTruth.PrimaryVerticesY.resize(nMCTruths);
   eventTruth.PrimaryVerticesZ.resize(nMCTruths);
   
   eventTruth.hasKaonShort.resize(nMCTruths);
   eventTruth.hasKaonShortCharged.resize(nMCTruths);

   eventTruth.PrimaryLeptons.clear();
   eventTruth.PrimaryKaonShorts.clear();

   eventTruth.KaonShortChargedDaughters.clear();

   eventTruth.eventHasKaonShort = false;
   eventTruth.eventHasHyperon = false;
   eventTruth.eventHasNeutronScatter = false;

   GetPrimaryParticles();

   if(primaryKaonShorts.size()) ConstructKaonShorts(); 

   SetFlags(); 

   return eventTruth;
}

void ParticleTrackerAnalyser::GetPrimaryParticles()
{
   for(size_t p = 0; p < primaryParticles.size(); p++){

      if(particleMap.find(primaryParticles[p]) == particleMap.end()) continue;

      art::Ptr<simb::MCParticle> TruthParticle = particleMap[primaryParticles.at(p)];
      int PdgCode = TruthParticle->PdgCode();

      if(PdgCode > 10000) continue;

      SimParticle MCParticle = MakeSimParticle(*TruthParticle);
      MCParticle.Origin = 1;
      MCTruthMatch(MCParticle);

      if(isLepton(PdgCode)) eventTruth.PrimaryLeptons.push_back(MCParticle);
      else if(isPion(PdgCode)) eventTruth.PrimaryPions.push_back(MCParticle);
      else if(isHyperon(PdgCode)) eventTruth.PrimaryHyperons.push_back(MCParticle); 
   }
}

void ParticleTrackerAnalyser::ConstructKaonShorts()
{
   for(size_t d = 0; d < primaryKaonShorts.size(); d++){
      if(particleMap.find(primaryKaonShorts[d]) == particleMap.end()) continue;

      art::Ptr<simb::MCParticle> TruthParticle = particleMap[primaryKaonShorts[d]];
       
      SimParticle MCParticle = MakeSimParticle(*TruthParticle);
      MCParticle.Origin = 7;

      MCTruthMatch(MCParticle, primaryKaonShorts[d]);
  
      eventTruth.PrimaryKaonShorts.push_back(MCParticle);     
   }

    for(size_t d = 0; d < primaryKaonShortDaughters.size(); d++){
      if(particleMap.find(primaryKaonShortDaughters[d]) == particleMap.end()) continue;

      art::Ptr<simb::MCParticle> TruthParticle = particleMap[primaryKaonShortDaughters[d]];
       
      SimParticle MCParticle = MakeSimParticle(*TruthParticle);
      MCParticle.Origin = 8;

      MCTruthMatch(MCParticle, primaryKaonShortDaughters[d]);
  
      eventTruth.KaonShortChargedDaughters.push_back(MCParticle);     
   }
}

bool ParticleTrackerAnalyser::FindNeutronScatter()
{
   std::vector<int> ScatterNeutrons;

   for(size_t i = 0; i < primaryParticles.size(); i++){
      if(particleMap.find(primaryParticles[i]) == particleMap.end()) continue;

      art::Ptr<simb::MCParticle> Particle = particleMap[primaryParticles[i]];

      if(Particle->PdgCode() != 2112) continue;
      if(Particle->EndProcess() == "neutronInelastic") ScatterNeutrons.push_back(Particle->TrackId()); 
   }

   for(size_t i = 0; i < ScatterNeutrons.size(); i++){
      art::Ptr<simb::MCParticle> particle = particleMap[ScatterNeutrons[i]];
      std::vector<int> ScatterNeutronDaughters = GetDaughters(particle);      

      int nProtons = 0, nPions = 0;

      for(size_t j = 0; j < ScatterNeutronDaughters.size(); j++){
         art::Ptr<simb::MCParticle> ScatterDaughter = particleMap[ScatterNeutronDaughters[i]];
         double Momentum = sqrt(ScatterDaughter->Px()*ScatterDaughter->Px() + ScatterDaughter->Py()*ScatterDaughter->Py() + ScatterDaughter->Pz()*ScatterDaughter->Pz());
         if(ScatterDaughter->PdgCode() == 2212 && Momentum > neutronScatterProtonThresh) nProtons++; 
         if(abs(ScatterDaughter->PdgCode()) == 211 && Momentum > neutronScatterPionThresh) nPions++; 
      }
       
      if(nProtons >= 2 || nPions >= 2 || (nProtons >= 1 && nPions >= 1)) return true;
   }

   return false;
}

void ParticleTrackerAnalyser::SetNeutronScatterThresholds(double ProtonThreshold, double PionThreshold)
{
   neutronScatterProtonThresh = ProtonThreshold;
   neutronScatterPionThresh = PionThreshold;
}

void ParticleTrackerAnalyser::SetDecayThresholds(double ProtonThreshold, double PionThreshold)
{
   decayProtonThresh = ProtonThreshold;
   decayPionThresh = PionThreshold;
}

void ParticleTrackerAnalyser::MCTruthMatch(SimParticle &MCParticle)
{ 
   for(size_t i = 0; i < primaryVertices.size(); i++){
      if(ProximityMatch(TVector3(MCParticle.StartX, MCParticle.StartY, MCParticle.StartZ), primaryVertices.at(i))) MCParticle.MCTruthIndex = i;
   }
}

void ParticleTrackerAnalyser::MCTruthMatch(SimParticle &MCParticle, int TrackID)
{
   if(particleMap.find(TrackID) == particleMap.end()) {
      MCParticle.MCTruthIndex = -1;
      return;
   }
  
  art::Ptr<simb::MCParticle> Particle = particleMap.at(TrackID);
   while(true){
      if(particleMap.find(Particle->Mother()) == particleMap.end() || Particle->Mother() == 0) break;
      Particle = particleMap.at(Particle->Mother());
   }

   for(size_t i = 0; i < primaryVertices.size(); i++){
      if(ProximityMatch(TVector3(Particle->Position().X(), Particle->Position().Y(), Particle->Position().Z()), primaryVertices.at(i))) MCParticle.MCTruthIndex = i;
   }
}

void ParticleTrackerAnalyser::SetFlags()
{
   for(int i = 0; i < nMCTruths; i++){

      eventTruth.inActiveTPC[i] = inActiveTPC(primaryVertices.at(i));

      int NHyperons = 0, NKaonShorts = 0, NPions = 0;

      for(size_t h = 0; h < eventTruth.PrimaryHyperons.size(); h++){
         if(eventTruth.PrimaryHyperons.at(h).MCTruthIndex == i){
            NHyperons++;
         }
      }

      for(size_t k = 0; k < eventTruth.PrimaryKaonShorts.size(); k++){
         if(eventTruth.PrimaryKaonShorts.at(k).MCTruthIndex == i){
            NKaonShorts++;
            eventTruth.hasKaonShort[k] = true;
         }
      }

      for(size_t p = 0; p < eventTruth.PrimaryPions.size(); p++){
         if(eventTruth.PrimaryPions.at(p).MCTruthIndex == i && isPion(eventTruth.PrimaryPions.at(p).PDG)){
            NPions++;
         }
      }

      bool HasPionPlus = false, HasPionMinus = false;

      for(size_t d = 0; d < eventTruth.KaonShortChargedDaughters.size(); d++){
         if(eventTruth.KaonShortChargedDaughters.size() != 2) continue; 
         if(eventTruth.KaonShortChargedDaughters.at(d).MCTruthIndex == i){
            if(eventTruth.KaonShortChargedDaughters.at(d).PDG == 221) HasPionPlus = true;  
            if(eventTruth.KaonShortChargedDaughters.at(d).PDG == -211) HasPionMinus = true;  
         }

         if(HasPionPlus && HasPionMinus) eventTruth.hasKaonShortCharged[i] = true;
      }

      eventTruth.nHyperons.push_back(NHyperons);
      eventTruth.nKaonShorts.push_back(NKaonShorts);
      eventTruth.nPions.push_back(NPions);
   }

   eventTruth.eventHasNeutronScatter = FindNeutronScatter();   
   eventTruth.eventHasHyperon = (eventTruth.PrimaryHyperons.size() > 0);
   eventTruth.eventHasKaonShort = (eventTruth.PrimaryKaonShorts.size() > 0);
}

int ParticleTrackerAnalyser::GetOrigin(int TrackIdentifier){
   if(std::find(primaryParticles.begin(), primaryParticles.end(), TrackIdentifier) != primaryParticles.end()){
      return 1;
   }
   else if(std::find(primaryKaonShorts.begin(), primaryKaonShorts.end(), TrackIdentifier) != primaryKaonShorts.end()){
      return 2;
   }
   else if(std::find(primaryKaonShortDaughters.begin(), primaryKaonShortDaughters.end(), TrackIdentifier) != primaryKaonShortDaughters.end()){
      return 3;
   }
   else{
      return 4;
   }
}

#endif