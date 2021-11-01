#ifndef _SubModuleG4Truth_cxx_
#define _SubModuleG4Truth_cxx_

#include "SubModuleG4Truth.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(art::Event const& e,std::string genlabel,std::string g4label){

   if(!e.getByLabel(genlabel,Handle_MCTruth))  
      throw cet::exception("SubModuleG4Truth") << "No MC Truth data product!" << std::endl;

   if(!e.getByLabel(g4label,Handle_G4)) 
      throw cet::exception("SubModuleG4Truth") << "No Geant4 Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  
   art::fill_ptr_vector(Vect_G4,Handle_G4);

   // Create map between particle ID and Particles
   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
      partByID.insert(std::make_pair(g4p->TrackId(),g4p));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(art::Event const& e,fhicl::ParameterSet pset) :
   SubModuleG4Truth(e,
         pset.get<std::string>("GeneratorModuleLabel","generator"),
         pset.get<std::string>("G4ModuleLabel","largeant"))
{
   SetNeutronScatterThresholds(pset.get<double>("NeutronScatterProtonThresh",0.15),pset.get<double>("NeutronScatterPionThresh",0.05));
   SetDecayThresholds(pset.get<double>("DecayProtonThresh",0.0),pset.get<double>("DecayPionThresh",0.0));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetParticleLists(){

   Primary_IDs.clear();
   Daughter_IDs.clear();
   SigmaZero_Daughter_IDs.clear();
   Kaon_Daughter_IDs.clear();

   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){

      // Skip any particles not produced at PV
      if(g4p->Mother() != 0)  continue;

      Primary_IDs.push_back(g4p->TrackId());

      std::vector<int> IDs = GetChildIDs(g4p);

      if(isHyperon(g4p->PdgCode())){

         if(g4p->PdgCode() == 3212 && g4p->EndProcess() == "Decay"){
            SigmaZero_Daughter_IDs.insert(SigmaZero_Daughter_IDs.begin(),IDs.begin(),IDs.end());                         

         }

         else if(g4p->EndProcess() == "Decay")
            Daughter_IDs.insert(Daughter_IDs.begin(),IDs.begin(),IDs.end());                

      }		

      else if(isKaon(g4p->PdgCode()) && g4p->EndProcess() == "Decay")
         Kaon_Daughter_IDs.insert(Kaon_Daughter_IDs.begin(),IDs.begin(),IDs.end());                     

   }


   // Get the SigmaZero daughters  
   for(size_t i_d=0;i_d<SigmaZero_Daughter_IDs.size();i_d++){

      if(partByID.find(SigmaZero_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[SigmaZero_Daughter_IDs.at(i_d)];

      if(part->PdgCode() == 3122){
         std::vector<int> IDs = GetChildIDs(part);
         Daughter_IDs.insert(Daughter_IDs.begin(),IDs.begin(),IDs.end());                     
      }     

   } 

   // Set list of Primary vertices for matching to multiple MCTruths

   for(const art::Ptr<simb::MCTruth> &theMCTruth : Vect_MCTruth){        
      for(int k_particles=0;k_particles<theMCTruth->NParticles();k_particles++){
         simb::MCParticle Part = theMCTruth->GetParticle(k_particles);
         if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) 
            PrimaryVertices.push_back(TVector3(Part.Vx(),Part.Vy(),Part.Vz()));

      }       
   }

   NMCTruths = Vect_MCTruth.size();

   if(PrimaryVertices.size() != Vect_MCTruth.size())         
      throw cet::exception("SubModuleG4Truth") << "Vertex/MCTruth vector size mismatch" << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

G4Truth SubModuleG4Truth::GetG4Info(){

   GetParticleLists();

   // Clear everything
   /*
      theTruth.IsHyperon = false;
      theTruth.IsLambda = false;
      theTruth.IsLambdaCharged = false;
      */

   theTruth.InActiveTPC.resize(NMCTruths);
   theTruth.IsHyperon.resize(NMCTruths);
   theTruth.IsLambda.resize(NMCTruths);
   theTruth.IsLambdaCharged.resize(NMCTruths);
   theTruth.IsSigmaZero.resize(NMCTruths);
   theTruth.IsSigmaZeroCharged.resize(NMCTruths);
   theTruth.IsAssociatedHyperon.resize(NMCTruths);

   theTruth.DecayVertex_X.resize(NMCTruths);
   theTruth.DecayVertex_Y.resize(NMCTruths);
   theTruth.DecayVertex_Z.resize(NMCTruths);

   for(int i_t=0;i_t<NMCTruths;i_t++){
      theTruth.DecayVertex_X[i_t] = -1000;
      theTruth.DecayVertex_Y[i_t] = -1000;
      theTruth.DecayVertex_Z[i_t] = -1000;
   }

   theTruth.Lepton.clear();
   theTruth.Hyperon.clear();
   theTruth.PrimaryNucleon.clear();
   theTruth.PrimaryPion.clear();
   theTruth.PrimaryKaon.clear();
   theTruth.Decay.clear();
   theTruth.KaonDecay.clear();
   theTruth.SigmaZeroDecayPhoton.clear();
   theTruth.SigmaZeroDecayLambda.clear();

   GetPrimaryParticles();

   if(SigmaZero_Daughter_IDs.size()) GetSigmaZeroDecay(); 
   if(theTruth.Hyperon.size() || theTruth.SigmaZeroDecayLambda.size()) GetHyperonDecay();
   if(theTruth.PrimaryKaon.size()) GetKaonDecay();

   SetFlags(); 

   /*    
         if(theTruth.Hyperon.size()) theTruth.IsHyperon = true;
         if(theTruth.Hyperon.size() == 1 && theTruth.Hyperon.at(0).PDG == 3122) theTruth.IsLambda = true;
         if(theTruth.Hyperon.size() == 1 && theTruth.Hyperon.at(0).PDG == 3212) theTruth.IsSigmaZero = true;
         if(theTruth.Hyperon.size() && theTruth.PrimaryKaon.size()) theTruth.IsAssociatedHyperon = true;


         if(theTruth.Decay.size() == 2){
         if(theTruth.IsLambda && theTruth.Decay.at(0).PDG == 2212 && theTruth.Decay.at(1).PDG == -211) theTruth.IsLambdaCharged = true;
         if(theTruth.IsLambda && theTruth.Decay.at(1).PDG == 2212 && theTruth.Decay.at(0).PDG == -211) theTruth.IsLambdaCharged = true; 
         }
         */


   return theTruth;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetPrimaryParticles(){

   for(size_t i_p=0;i_p<Primary_IDs.size();i_p++){

      // Geant does not always keep everything it simulates, make sure particle is in list of IDs (crashes otherwise!)
      if(partByID.find(Primary_IDs[i_p]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs.at(i_p)];

      // Anything with very large pdg code is a nucleus, skip these
      if(part->PdgCode() > 10000) continue;

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 1;
      MCTruthMatch(P);

      //hyperon produced at primary vertex
      if(isHyperon(part->PdgCode())){
         theTruth.Hyperon.push_back(P);
         //theTruth.IsHyperon = true;	
      }

      if(isLepton(part->PdgCode()) || isNeutrino(part->PdgCode())) theTruth.Lepton.push_back(P);
      if(isNucleon(part->PdgCode())) theTruth.PrimaryNucleon.push_back(P);
      if(isPion(part->PdgCode())) theTruth.PrimaryPion.push_back(P);
      if(isKaon(part->PdgCode())) theTruth.PrimaryKaon.push_back(P);

   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetHyperonDecay(){

   if(!Daughter_IDs.size()) return;

   for(size_t i_d=0;i_d<Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 2;

      // Check which MCTruth this decay belongs to
      for(size_t i_h=0;i_h<theTruth.Hyperon.size();i_h++){
         SimParticle H = theTruth.Hyperon.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))){
            P.MCTruthIndex = H.MCTruthIndex;            
            theTruth.DecayVertex_X[P.MCTruthIndex] = P.StartX;
            theTruth.DecayVertex_Y[P.MCTruthIndex] = P.StartY;
            theTruth.DecayVertex_Z[P.MCTruthIndex] = P.StartZ;
         }
      } 

      for(size_t i_h=0;i_h<theTruth.SigmaZeroDecayLambda.size();i_h++){
         SimParticle H = theTruth.SigmaZeroDecayLambda.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))) P.MCTruthIndex = H.MCTruthIndex;            
      } 


      theTruth.Decay.push_back(P);     
   }


   //theTruth.DecayVertex = TVector3(theTruth.Decay.at(0).StartX,theTruth.Decay.at(0).StartY,theTruth.Decay.at(0).StartZ);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetKaonDecay(){

   for(size_t i_d=0;i_d<Kaon_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Kaon_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Kaon_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 4;

      // Check which MCTruth this decay belongs to
      for(size_t i_k=0;i_k<theTruth.PrimaryKaon.size();i_k++){
         SimParticle K = theTruth.PrimaryKaon.at(i_k);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(K.EndX,K.EndY,K.EndZ))) P.MCTruthIndex = K.MCTruthIndex;
      } 

      theTruth.Decay.push_back(P);     
   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetSigmaZeroDecay(){

   for(size_t i_d=0;i_d<SigmaZero_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(SigmaZero_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[SigmaZero_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 5;

      // Check which MCTruth this decay belongs to
      for(size_t i_h=0;i_h<theTruth.Hyperon.size();i_h++){
         SimParticle H = theTruth.Hyperon.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))) P.MCTruthIndex = H.MCTruthIndex;
      } 

      if(part->PdgCode() == 3122)
         theTruth.SigmaZeroDecayLambda.push_back(P);     
      else if(part->PdgCode() == 22)
         theTruth.SigmaZeroDecayPhoton.push_back(P);     

   }



}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> SubModuleG4Truth::GetChildIDs(art::Ptr<simb::MCParticle> g4p,bool IsNeutron){

   std::vector<int> DecayProduct_IDs;

   if(g4p->EndProcess() != "Decay" && !IsNeutron) return DecayProduct_IDs;

   for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[g4p->Daughter(i_d)];

      if(part->PdgCode() > 10000) continue;

      double X = part->Position().X();
      double Y = part->Position().Y();
      double Z = part->Position().Z();

      double EndX = g4p->EndPosition().X();
      double EndY = g4p->EndPosition().Y();
      double EndZ = g4p->EndPosition().Z();

      //if(X != EndX || Y != EndY || Z != EndZ) continue;

      if(!PosMatch(TVector3(X,Y,Z),TVector3(EndX,EndY,EndZ))) continue;

      DecayProduct_IDs.push_back(g4p->Daughter(i_d));
   }

   return DecayProduct_IDs;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int SubModuleG4Truth::GetOrigin(int trackid){

   if(std::find(Primary_IDs.begin(),Primary_IDs.end(),trackid) != Primary_IDs.end()) return 1;
   else if(std::find(Daughter_IDs.begin(),Daughter_IDs.end(),trackid) != Daughter_IDs.end()) return 2;
   else if(std::find(Kaon_Daughter_IDs.begin(),Kaon_Daughter_IDs.end(),trackid) != Kaon_Daughter_IDs.end()) return 4;
   else if(std::find(SigmaZero_Daughter_IDs.begin(),SigmaZero_Daughter_IDs.end(),trackid) != SigmaZero_Daughter_IDs.end()) return 5;
   else return 3;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleG4Truth::FindNeutronScatter(){

   std::vector<int> Neutron_IDs;

   for(size_t i=0;i<Primary_IDs.size();i++){

      if(partByID.find(Primary_IDs[i]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs[i]];

      if(part->PdgCode() != 2112) continue;

      if(part->EndProcess() == "neutronInelastic") Neutron_IDs.push_back(part->TrackId()); 

   }

   for(size_t i=0;i<Neutron_IDs.size();i++){
      art::Ptr<simb::MCParticle> part = partByID[Neutron_IDs[i]];
      std::vector<int> NeutronDaughter_IDs = GetChildIDs(part,true);      

      int nProtons = 0;
      int nPions = 0;

      for(size_t i_d=0;i_d<NeutronDaughter_IDs.size();i_d++){

         // Look for protons and charged pions
         art::Ptr<simb::MCParticle> part2 = partByID[NeutronDaughter_IDs[i_d]];


         double P = sqrt(part2->Px()*part2->Px() + part2->Py()*part2->Py() + part2->Pz()*part2->Pz());

         if(part2->PdgCode() == 2212 && P > NeutronScatterProtonThresh) nProtons++; 
         if(abs(part2->PdgCode()) == 211 && P > NeutronScatterPionThresh) nPions++; 

      }

      // Flag event as containing neutron scatter if there are either 2+ protons, 2+ charged pions or 1+ of each
      if(nProtons >= 2 || nPions >= 2 || (nProtons >= 1 && nPions >= 1)) return true;

   }

   return false;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetNeutronScatterThresholds(double neutronscatterprotonthresh,double neutronscatterpionthresh){

   NeutronScatterProtonThresh = neutronscatterprotonthresh;
   NeutronScatterPionThresh = neutronscatterpionthresh;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetDecayThresholds(double decayprotonthresh,double decaypionthresh){

   DecayProtonThresh = decayprotonthresh;
   DecayPionThresh = decaypionthresh;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleG4Truth::PosMatch(TVector3 Pos1,TVector3 Pos2){

   return (Pos1-Pos2).Mag() < _EPSILON_;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::MCTruthMatch(SimParticle &P){

   for(size_t i_mct=0;i_mct<PrimaryVertices.size();i_mct++)
      if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),PrimaryVertices.at(i_mct)))
         P.MCTruthIndex = i_mct;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetFlags(){

   /*
      if(theTruth.Hyperon.size()) theTruth.IsHyperon = true;
      if(theTruth.Hyperon.size() == 1 && theTruth.Hyperon.at(0).PDG == 3122) theTruth.IsLambda = true;
      if(theTruth.Hyperon.size() == 1 && theTruth.Hyperon.at(0).PDG == 3212) theTruth.IsSigmaZero = true;

      if(theTruth.Decay.size() == 2){
      if(theTruth.IsLambda && theTruth.Decay.at(0).PDG == 2212 && theTruth.Decay.at(1).PDG == -211) theTruth.IsLambdaCharged = true;
      if(theTruth.IsLambda && theTruth.Decay.at(1).PDG == 2212 && theTruth.Decay.at(0).PDG == -211) theTruth.IsLambdaCharged = true; 
      }


   // Set Assoc hyperon to true if there is a hyperon and kaon from the same PV
   for(size_t i_h=0;i_h<theTruth.Hyperon.size();i_h++)
   for(size_t i_k=0;i_k<theTruth.PrimaryKaon.size();i_k++)
   if(theTruth.Hyperon.at(i_h).MCTruthIndex == theTruth.PrimaryKaon.at(i_k).MCTruthIndex)
   theTruth.IsAssociatedHyperon = true;
   */

   // Iterate over the different MCTruths in the event, set their flags
   for(int i_t=0;i_t<NMCTruths;i_t++){

      theTruth.InActiveTPC[i_t] = inActiveTPC(PrimaryVertices.at(i_t));

      theTruth.IsHyperon[i_t] = false;
      theTruth.IsLambda[i_t] = false;
      theTruth.IsLambdaCharged[i_t] = false;
      theTruth.IsSigmaZero[i_t] = false;
      theTruth.IsAssociatedHyperon[i_t] = false;

      int nHyperons=0,nKaons=0;

      for(size_t i_h=0;i_h<theTruth.Hyperon.size();i_h++){
         if(theTruth.Hyperon.at(i_h).MCTruthIndex == i_t){
            nHyperons++;
            theTruth.IsHyperon[i_t] = true;
         }
         if(theTruth.Hyperon.at(i_h).MCTruthIndex == i_t && theTruth.Hyperon.at(i_h).PDG == 3122) theTruth.IsLambda[i_t] = true;
         if(theTruth.Hyperon.at(i_h).MCTruthIndex == i_t && theTruth.Hyperon.at(i_h).PDG == 3212) theTruth.IsSigmaZero[i_t] = true;
      }

      int nProducts=0;
      bool hasProton=false,hasPion=false;

      for(size_t i_d=0;i_d<theTruth.Decay.size();i_d++){
         if(theTruth.Decay.at(i_d).MCTruthIndex == i_t){
            nProducts++;
            if(theTruth.Decay.at(i_d).PDG == 2212) hasProton = true;  
            if(theTruth.Decay.at(i_d).PDG == -211) hasPion = true;  
         }
      }
     
      if(nHyperons == 1 && theTruth.IsLambda.at(i_t) && nProducts == 2 && hasProton && hasPion) theTruth.IsLambdaCharged[i_t] = true;
      if(nHyperons == 1 && theTruth.IsSigmaZero.at(i_t) && nProducts == 2 && hasProton && hasPion) theTruth.IsSigmaZeroCharged[i_t] = true;

      for(size_t i_k=0;i_k<theTruth.PrimaryKaon.size();i_k++)
         if(theTruth.PrimaryKaon.at(i_k).MCTruthIndex == i_t) nKaons++;

      // If there are multiple hyperons/antihyperons or hyperons and kaons present together
      // flag as associated hyperon
      if(nHyperons > 1 || (nHyperons >= 1 && nKaons >= 1)) theTruth.IsAssociatedHyperon[i_t] = true;

   }

   theTruth.EventHasNeutronScatter = FindNeutronScatter();   
   theTruth.EventHasHyperon = std::find(theTruth.IsHyperon.begin(), theTruth.IsHyperon.end(), true) != theTruth.IsHyperon.end();


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
