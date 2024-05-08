#ifndef ReconstructionAnalyser_cxx_
#define _ReconstructionAnalyser_cxx_

#include "ReconstructionAnalyser.h"

using namespace hyperon;

ReconstructionAnalyser::ReconstructionAnalyser(art::Event const& E, bool IsData, fhicl::ParameterSet ParamSet, bool ParticleGunMode) :
ReconstructionAnalyser(E, IsData,
                  ParamSet.get<std::string>("PFParticleModuleLabel"),
                  ParamSet.get<std::string>("TrackModuleLabel"),
                  ParamSet.get<std::string>("ShowerModuleLabel"),
                  ParamSet.get<std::string>("VertexModuleLabel"),
                  ParamSet.get<std::string>("PIDModuleLabel"),
                  ParamSet.get<std::string>("CaloModuleLabel"),
                  ParamSet.get<std::string>("HitModuleLabel"),
                  ParamSet.get<std::string>("HitTruthAssnLabel"),
                  ParamSet.get<std::string>("TrackHitAssnLabel"),
                  ParamSet.get<std::string>("MetadataModuleLabel"),
                  ParamSet.get<std::string>("GeneratorModuleLabel"),
                  ParamSet.get<std::string>("G4ModuleLabel"),
                  ParamSet.get<bool>("DoGetPIDs", true),
                  ParamSet.get<bool>("includeCosmics", false),
                  ParticleGunMode)
{}

ReconstructionAnalyser::ReconstructionAnalyser(art::Event const& E, bool IsData, 
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
                  bool includeCosmics, 
                  bool ParticleGunMode):
PIDCalc(),
doGetPIDs(DoGetPIDs),
includeCosmics(includeCosmics),
particleGunMode(ParticleGunMode)
{
   isData = IsData;

   if(!E.getByLabel(PFParticleLabel, handlePFParticle)) 
      throw cet::exception("ReconstructionAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   if(!E.getByLabel(TrackLabel, handleTrack)) 
      throw cet::exception("ReconstructionAnalyser") << "No Track Data Products Found!" << std::endl;

   if(!E.getByLabel(ShowerLabel, handleShower)) 
      throw cet::exception("ReconstructionAnalyser") << "No Shower Data Products Found!" << std::endl;

   if(!E.getByLabel(HitLabel, handleHit)) 
      throw cet::exception("ReconstructionAnalyser") << "No Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(vectPFParticle, handlePFParticle);
   art::fill_ptr_vector(vectTrack, handleTrack);
   art::fill_ptr_vector(vectShower, handleShower);
   art::fill_ptr_vector(vectHit, handleHit);

   assocPFParticleVertex = new art::FindManyP<recob::Vertex>(vectPFParticle, E, VertexLabel);    
   assocPFParticleTrack = new art::FindManyP<recob::Track>(vectPFParticle, E, TrackLabel);    
   assocPFParticleShower = new art::FindManyP<recob::Shower>(vectPFParticle, E, ShowerLabel);    
   assocPFParticleMetadata = new art::FindManyP<larpandoraobj::PFParticleMetadata>(vectPFParticle, E, MetaDataLabel);   
   assocTrackHit = new  art::FindManyP<recob::Hit>(vectTrack, E, TrackHitAssnLabel);
   particlesPerHit = new art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>(handleHit, E, HitTruthAssnLabel);

   if(DoGetPIDs){
      assocTrackCalo = new art::FindManyP<anab::Calorimetry>(vectTrack, E, CaloLabel);
      assocTrackPID = new art::FindManyP<anab::ParticleID>(vectTrack, E, PIDLabel);
   }

   if(!isData){
      particleTrackerAnalyser = new ParticleTrackerAnalyser(E, GenLabel, G4Label, ParticleGunMode);
      particleTrackerAnalyser->GetParticleLists();
   }

   particleMap.clear(); 
}

void ReconstructionAnalyser::PrepareInfo()
{
   eventReco.RecoPrimaryVertex = GetPrimaryVertex();

   for(const art::Ptr<recob::PFParticle> &PFP : vectPFParticle){
      if(!includeCosmics && PFP->Parent() != neutrinoIndex && particleMap.find(PFP->Parent()) == particleMap.end()) continue; 
      RecoParticle Particle = MakeRecoParticle(PFP);
      
      if(PFP->Parent() == neutrinoIndex){
         Particle.Parentage = 1;
         Particle.InNuSlice = true;         
      }
      else if(particleMap.find(PFP->Parent()) != particleMap.end()){         
         Particle.Parentage = 2; 
         Particle.ParentIndex = particleMap[PFP->Parent()]; 
      }

      // Pandora stores tracks as a muon, and showers as an electron
      if(Particle.PDG == 13){
         eventReco.TrackPrimaryDaughters.push_back(Particle);
         if(Particle.InNuSlice) particleMap[PFP->Self()] = Particle.Index;
      }
      else if(Particle.PDG == 11) eventReco.ShowerPrimaryDaughters.push_back(Particle);      
   }

   eventReco.NPrimaryDaughters = eventReco.TrackPrimaryDaughters.size() + eventReco.ShowerPrimaryDaughters.size();
   eventReco.NPrimaryTrackDaughters = eventReco.TrackPrimaryDaughters.size();
   eventReco.NPrimaryShowerDaughters = eventReco.ShowerPrimaryDaughters.size();
}

RecoData ReconstructionAnalyser::GetInfo()
{
   return eventReco;
}

TVector3 ReconstructionAnalyser::GetPrimaryVertex()
{
   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   for(const art::Ptr<recob::PFParticle> &PFP : vectPFParticle){

      if(PFP->IsPrimary() && isNeutrino(PFP->PdgCode())){
         neutrinoIndex = PFP->Self();

         std::vector<art::Ptr<recob::Vertex>> PFPVertex = assocPFParticleVertex->at(PFP.key());

         for(const art::Ptr<recob::Vertex> &Vertex : PFPVertex){

            geo::Point_t Point = { Vertex->position().X() , Vertex->position().Y() , Vertex->position().Z() };
            geo::Vector_t SpaceChargeCorrections = SCE->GetPosOffsets(Point);

            return TVector3(Vertex->position().X() + SpaceChargeCorrections.X(), 
                           Vertex->position().Y() - SpaceChargeCorrections.Y(), 
                           Vertex->position().Z() - SpaceChargeCorrections.Z());
         }
      }
   }

   return TVector3(-1000,-1000,-1000);
}


RecoParticle ReconstructionAnalyser::MakeRecoParticle(const art::Ptr<recob::PFParticle> &PFP)
{
   RecoParticle Particle;

   Particle.PDG = PFP->PdgCode();

   std::vector<art::Ptr<recob::Track>> VecTracks = assocPFParticleTrack->at(PFP.key());
   std::vector<art::Ptr<recob::Shower>> VecShowers = assocPFParticleShower->at(PFP.key());

   if(PFP->PdgCode() == 13 && VecTracks.size() != 1) Particle.PDG = 0;
   if(PFP->PdgCode() == 11 && VecShowers.size() != 1) Particle.PDG = 0;

   SetPFPMetaData(PFP, Particle);

   if(VecTracks.size() == 1){
      GetTrackData(PFP, Particle);
      GetVertexData(PFP, Particle);
   }

   return Particle;
}

void ReconstructionAnalyser::SetPFPMetaData(const art::Ptr<recob::PFParticle> &PFP, RecoParticle &Particle)
{
   std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> PFPMeta = assocPFParticleMetadata->at(PFP.key());

   for(const art::Ptr<larpandoraobj::PFParticleMetadata> &Meta : PFPMeta){

      const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(Meta->GetPropertiesMap());

      if (!pfParticlePropertiesMap.empty()){
         for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){
            if(it->first == "TrackScore"){
               Particle.TrackShowerScore = it->second;
            }
            else if(it->first == "NuScore"){
               Particle.NuScore = it->second;
            }
         }
      }	
   }

}

void ReconstructionAnalyser::GetTrackData(const art::Ptr<recob::PFParticle> &PFP, RecoParticle &Particle)
{
   std::vector<art::Ptr<recob::Track>> VecTracks = assocPFParticleTrack->at(PFP.key());

   if(VecTracks.size() != 1) return;

   art::Ptr<recob::Track> Trk = VecTracks.at(0);

   SetTrackVariables(Particle, Trk);

   if(!isData) TruthMatch(Trk, Particle);
   if(doGetPIDs) GetPIDs(Trk, Particle);
   
   eventReco.TrackStarts.push_back(TVector3(Trk->Start().X(), Trk->Start().Y(), Trk->Start().Z()));
   Particle.Index = eventReco.TrackStarts.size() - 1;
}

void ReconstructionAnalyser::TruthMatch(const art::Ptr<recob::Track> &Trk, RecoParticle &Particle)
{
   std::vector<art::Ptr<recob::Hit>> VecHits = assocTrackHit->at(Trk.key());

   std::unordered_map<int, double> TrackMap;
   int MaxHits = -1;

   simb::MCParticle const* MatchedParticle = NULL;

   std::vector<simb::MCParticle const*> VecMCParticle;
   std::vector<anab::BackTrackerHitMatchingData const*> VecMathingHits;

   for(size_t h = 0; h < VecHits.size(); ++h){

      VecMCParticle.clear();
      VecMathingHits.clear();

      particlesPerHit->get(VecHits[h].key(), VecMCParticle, VecMathingHits);

      for(size_t p = 0; p < VecMCParticle.size(); ++p){
         TrackMap[VecMCParticle[p]->TrackId()]++; 

         //Choose particle depositing energy in the most VecHits
         if(TrackMap[VecMCParticle[p]->TrackId()] > MaxHits){
            MaxHits = TrackMap[VecMCParticle[p]->TrackId()];
            MatchedParticle = VecMCParticle[p];
         }
      }
   }

   if(MatchedParticle != NULL){ 
      SimParticle SP = MakeSimParticle(*MatchedParticle);
            
      SP.Origin = particleTrackerAnalyser->GetOrigin(MatchedParticle->TrackId());
      particleTrackerAnalyser->MCTruthMatch(SP, MatchedParticle->TrackId());
 
      Particle.HasTruth = true;
      Particle.MCTruthIndex = SP.MCTruthIndex;
      Particle.TrackTruePDG = SP.PDG;
      Particle.TrackTrueE = SP.E;
      Particle.TrackTruePx = SP.Px;
      Particle.TrackTruePy = SP.Py;
      Particle.TrackTruePz = SP.Pz;
      Particle.TrackTrueEndE = SP.E;
      Particle.TrackTrueEndPx = SP.EndPx;
      Particle.TrackTrueEndPy = SP.EndPy;
      Particle.TrackTrueEndPz = SP.EndPz;
      Particle.TrackTrueModMomentum = SP.ModMomentum;
      Particle.TrackTrueEndModMomentum = SP.EndModMomentum;
      Particle.TrackTrueKE = SP.KE;
      Particle.TrackTrueEndKE = SP.EndKE;
      Particle.TrackTrueLength = SP.Travel;
      Particle.TrackTrueOrigin = SP.Origin;
      Particle.TrackTruthPurity = (double)MaxHits/VecHits.size();
   }
   else Particle.HasTruth = false;
}

void ReconstructionAnalyser::GetPIDs(const art::Ptr<recob::Track> &Trk, RecoParticle &Particle)
{
   std::vector<art::Ptr<anab::Calorimetry>> CaloFromTrack = assocTrackCalo->at(Trk.key());
   std::vector<art::Ptr<anab::ParticleID>> TrackPID = assocTrackPID->at(Trk.key());
   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = TrackPID.at(0)->ParticleIDAlgScores();

   PIDStore Store = PIDCalc.GetPIDScores(Trk, CaloFromTrack, AlgScoresVec);

   Particle.Track_LLR_PID = Store.LLR;

   Particle.MeandEdX_Plane0 = Store.MeandEdX_Plane0;
   Particle.MeandEdX_Plane1 = Store.MeandEdX_Plane1;
   Particle.MeandEdX_Plane2 = Store.MeandEdX_Plane2;
   Particle.MeandEdX_ThreePlane = Store.MeandEdX_3Plane;
}

void ReconstructionAnalyser::GetVertexData(const art::Ptr<recob::PFParticle> &PFP, RecoParticle &Particle)
{
   std::vector<art::Ptr<recob::Vertex>> PFPVertex = assocPFParticleVertex->at(PFP.key());

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   for(const art::Ptr<recob::Vertex> &Vertex : PFPVertex){

      geo::Point_t Point = {Vertex->position().X(), Vertex->position().Y(), Vertex->position().Z()};                
      geo::Vector_t SpaceChargeCorrections = SCE->GetPosOffsets(Point);

      TVector3 Pos(Vertex->position().X() + SpaceChargeCorrections.X(), Vertex->position().Y() - SpaceChargeCorrections.Y(), Vertex->position().Z() - SpaceChargeCorrections.Z());

      Particle.SetVertex(Pos);
      Particle.Displacement = (Pos - eventReco.RecoPrimaryVertex).Mag();
   }
}

void ReconstructionAnalyser::SetIndices(std::vector<bool> IsSignal)
{      
   bool ContainsSignal = std::find(IsSignal.begin(), IsSignal.end(),true) == IsSignal.end();

   bool FoundMuon = false, FoundDecayPionPlus = false, FoundDecayPionMinus = false;

   for(size_t i = 0; i < eventReco.TrackPrimaryDaughters.size(); i++){

      RecoParticle Particle = eventReco.TrackPrimaryDaughters.at(i);

      if(!FoundMuon && abs(Particle.TrackTruePDG) == 13 && Particle.TrackTrueOrigin == 1){ 
         eventReco.TrueMuonIndex = Particle.Index;
         FoundMuon = true; 
      }

      if(ContainsSignal && !FoundDecayPionPlus && isPionPlus(Particle.TrackTruePDG) && Particle.TrackTrueOrigin == 3){
         eventReco.TrueDecayPionPlusIndex = Particle.Index;
         FoundDecayPionPlus = true;
      }

      if(ContainsSignal &&  FoundDecayPionMinus && isPionMinus(Particle.TrackTruePDG) && Particle.TrackTrueOrigin == 3){
         eventReco.TrueDecayPionMinusIndex = Particle.Index;
       FoundDecayPionMinus = true;
      }
   }

   for(size_t i = 0; i < eventReco.ShowerPrimaryDaughters.size(); i++){

      RecoParticle Particle = eventReco.ShowerPrimaryDaughters.at(i);

      if(!FoundMuon && isMuon(Particle.TrackTruePDG) && Particle.TrackTrueOrigin == 1){ 
         eventReco.TrueMuonIndex = Particle.Index;
         FoundMuon = true; 
      }

      if(ContainsSignal && !FoundDecayPionPlus && isPionPlus(Particle.TrackTruePDG) && Particle.TrackTrueOrigin == 3){
         eventReco.TrueDecayPionPlusIndex = Particle.Index;
         FoundDecayPionPlus = true;
      }

      if(ContainsSignal &&  FoundDecayPionMinus && isPionMinus(Particle.TrackTruePDG) && Particle.TrackTrueOrigin == 3){
         eventReco.TrueDecayPionMinusIndex = Particle.Index;
       FoundDecayPionMinus = true;
      }
   }

   if(ContainsSignal && FoundDecayPionPlus && FoundDecayPionMinus) eventReco.GoodReco = true;
}

#endif