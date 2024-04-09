#ifndef RecoVariables_cxx_
#define _RecoVariables_cxx_

#include "SubModuleReco.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//SubModuleReco::SubModuleReco(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleReco::SubModuleReco(art::Event const& e,bool isdata,fhicl::ParameterSet pset,bool particlegunmode) :
SubModuleReco(e,isdata,
                  pset.get<std::string>("PFParticleModuleLabel"),
                  pset.get<std::string>("TrackModuleLabel"),
                  pset.get<std::string>("ShowerModuleLabel"),
                  pset.get<std::string>("VertexModuleLabel"),
                  pset.get<std::string>("PIDModuleLabel"),
                  pset.get<std::string>("CaloModuleLabel"),
                  pset.get<std::string>("HitModuleLabel"),
                  pset.get<std::string>("HitTruthAssnLabel"),
                  pset.get<std::string>("TrackHitAssnLabel"),
                  pset.get<std::string>("MetadataModuleLabel"),
                  pset.get<std::string>("GeneratorModuleLabel"),
                  pset.get<std::string>("G4ModuleLabel"),
                  pset.get<bool>("DoGetPIDs",true),
                  pset.get<bool>("IncludeCosmics",false),
                  particlegunmode)
{

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleReco::SubModuleReco(art::Event const& e,bool isdata,string pfparticlelabel,string tracklabel,
                                     string showerlabel,string vertexlabel,string pidlabel,string calolabel,string hitlabel,
                                     string hittruthassnlabel,string trackhitassnlabel,string metadatalabel,string genlabel,
                                     string g4label,bool dogetpids,bool includecosmics,bool particlegunmode) :
PIDCalc(),
DoGetPIDs(dogetpids),
IncludeCosmics(includecosmics),
ParticleGunMode(particlegunmode)
{

   IsData = isdata;

   if(!e.getByLabel(pfparticlelabel,Handle_PFParticle)) 
      throw cet::exception("SubModuleReco") << "No PFParticle Data Products Found!" << std::endl;

   if(!e.getByLabel(tracklabel,Handle_Track)) 
      throw cet::exception("SubModuleReco") << "No Track Data Products Found!" << std::endl;

   if(!e.getByLabel(showerlabel,Handle_Shower)) 
      throw cet::exception("SubModuleReco") << "No Shower Data Products Found!" << std::endl;

   if(!e.getByLabel(hitlabel,Handle_Hit)) 
      throw cet::exception("SubModuleReco") << "No Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);
   art::fill_ptr_vector(Vect_Track,Handle_Track);
   art::fill_ptr_vector(Vect_Shower,Handle_Shower);
   art::fill_ptr_vector(Vect_Hit,Handle_Hit);

   Assoc_PFParticleVertex = new art::FindManyP<recob::Vertex>(Vect_PFParticle,e,vertexlabel);    
   Assoc_PFParticleTrack = new art::FindManyP<recob::Track>(Vect_PFParticle,e,tracklabel);    
   Assoc_PFParticleShower = new art::FindManyP<recob::Shower>(Vect_PFParticle,e,showerlabel);    
   Assoc_PFParticleMetadata = new art::FindManyP<larpandoraobj::PFParticleMetadata>(Vect_PFParticle,e,metadatalabel);   
   Assoc_TrackHit = new  art::FindManyP<recob::Hit>(Vect_Track,e,trackhitassnlabel);
   ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hit,e,hittruthassnlabel);

   if(DoGetPIDs){
      Assoc_TrackCalo = new art::FindManyP<anab::Calorimetry>(Vect_Track,e,calolabel);
      Assoc_TrackPID = new art::FindManyP<anab::ParticleID>(Vect_Track,e,pidlabel);
   }

   llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
   llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
   llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

   llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
   llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
   llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

   llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
   llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
   llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

   llr_pid_calculator_kaon.set_dedx_binning(0, kaonproton_parameters.dedx_edges_pl_0);
   llr_pid_calculator_kaon.set_par_binning(0, kaonproton_parameters.parameters_edges_pl_0);
   llr_pid_calculator_kaon.set_lookup_tables(0, kaonproton_parameters.dedx_pdf_pl_0);

   llr_pid_calculator_kaon.set_dedx_binning(1, kaonproton_parameters.dedx_edges_pl_1);
   llr_pid_calculator_kaon.set_par_binning(1, kaonproton_parameters.parameters_edges_pl_1);
   llr_pid_calculator_kaon.set_lookup_tables(1, kaonproton_parameters.dedx_pdf_pl_1);

   llr_pid_calculator_kaon.set_dedx_binning(2, kaonproton_parameters.dedx_edges_pl_2);
   llr_pid_calculator_kaon.set_par_binning(2, kaonproton_parameters.parameters_edges_pl_2);
   llr_pid_calculator_kaon.set_lookup_tables(2, kaonproton_parameters.dedx_pdf_pl_2);

   if(!IsData){
      G4T = new SubModuleG4Truth(e,genlabel,g4label,ParticleGunMode);
      G4T->GetParticleLists();
   }

   m_PFPID_TrackIndex.clear(); 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////


void SubModuleReco::PrepareInfo(){

   theData.RecoPrimaryVertex = GetPrimaryVertex();

   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
      if(!IncludeCosmics && pfp->Parent() != neutrinoID && m_PFPID_TrackIndex.find(pfp->Parent()) == m_PFPID_TrackIndex.end()) continue; 
      RecoParticle P = MakeRecoParticle(pfp);
      
      if(pfp->Parent() == neutrinoID){
         P.Parentage = 1;
         P.InNuSlice = true;         
      }
      else if(m_PFPID_TrackIndex.find(pfp->Parent()) != m_PFPID_TrackIndex.end()){         
         P.Parentage = 2; 
         P.ParentIndex = m_PFPID_TrackIndex[pfp->Parent()]; 
      }

      if(P.PDG == 13){
         theData.TrackPrimaryDaughters.push_back(P);
         if(P.InNuSlice) m_PFPID_TrackIndex[pfp->Self()] = P.Index;
      }
      else if(P.PDG == 11) theData.ShowerPrimaryDaughters.push_back(P);      
   }

   theData.NPrimaryDaughters = theData.TrackPrimaryDaughters.size() + theData.ShowerPrimaryDaughters.size();
   theData.NPrimaryTrackDaughters = theData.TrackPrimaryDaughters.size();
   theData.NPrimaryShowerDaughters = theData.ShowerPrimaryDaughters.size();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoData SubModuleReco::GetInfo(){
   return theData;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 SubModuleReco::GetPrimaryVertex(){

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){

      if(pfp->IsPrimary() && isNeutrino(pfp->PdgCode())){

         neutrinoID = pfp->Self();

         std::vector<art::Ptr<recob::Vertex>> pfpVertex = Assoc_PFParticleVertex->at(pfp.key());

         for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

            geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };
            geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

            return TVector3(vtx->position().X() + sce_corr.X(),vtx->position().Y()-sce_corr.Y(),vtx->position().Z()-sce_corr.Z());
         }
      }
   }

   // If there is no neutrino candidate
   return TVector3(-1000,-1000,-1000);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoParticle SubModuleReco::MakeRecoParticle(const art::Ptr<recob::PFParticle> &pfp){

   RecoParticle P;

   P.PDG = pfp->PdgCode();

   std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack->at(pfp.key());
   std::vector<art::Ptr<recob::Shower>> pfpShowers = Assoc_PFParticleShower->at(pfp.key());

   if(pfp->PdgCode() == 13 && pfpTracks.size() != 1) P.PDG = 0;
   if(pfp->PdgCode() == 11 && pfpShowers.size() != 1) P.PDG = 0;

   GetPFPMetadata(pfp,P);

   if(pfpTracks.size() == 1){
      GetTrackData(pfp,P);
      GetVertexData(pfp,P);
   }

   return P;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetPFPMetadata(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

   std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMeta = Assoc_PFParticleMetadata->at(pfp.key());

   for(const art::Ptr<larpandoraobj::PFParticleMetadata> &meta : pfpMeta){

      const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(meta->GetPropertiesMap());

      if (!pfParticlePropertiesMap.empty()){
         for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){
            if(it->first == "TrackScore"){
               P.TrackShowerScore = it->second;
            }
         }
      }	
   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

   std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack->at(pfp.key());

   if(pfpTracks.size() != 1) return;

   art::Ptr<recob::Track> trk = pfpTracks.at(0);

   // Sets track length/position related variables
   int ndesc = GetNDescendents(pfp,GetPFParticleMap(Vect_PFParticle)); 
   SetTrackVariables(P,trk,ndesc);

   if(!IsData) TruthMatch(trk,P);
   if(DoGetPIDs) GetPIDs(trk,P);
   
   theData.TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z()));
   P.Index = theData.TrackStarts.size() - 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::TruthMatch(const art::Ptr<recob::Track> &trk,RecoParticle &P){

   std::vector<art::Ptr<recob::Hit>> hits = Assoc_TrackHit->at(trk.key());

   std::unordered_map<int,double>  trkide;
   int maxhits=-1;

   simb::MCParticle const* matchedParticle = NULL;

   std::vector<simb::MCParticle const*> particleVec;
   std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

   for(size_t i_hit=0;i_hit<hits.size();++i_hit){

      particleVec.clear();
      matchVec.clear();
      ParticlesPerHit->get(hits[i_hit].key(),particleVec,matchVec);

      for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){

         trkide[particleVec[i_particle]->TrackId()]++; 

         //new method - choose particle depositing energy in the most hits
         if(trkide[particleVec[i_particle]->TrackId()] > maxhits){
            maxhits = trkide[particleVec[i_particle]->TrackId()];
            matchedParticle = particleVec[i_particle];
         }
      }
   }

   if(matchedParticle != NULL){ 

      SimParticle SP = MakeSimParticle(*matchedParticle);
            
      SP.Origin = G4T->GetOrigin(matchedParticle->TrackId());
      G4T->MCTruthMatch(SP,matchedParticle->TrackId());
 
      P.HasTruth = true;
      P.MCTruthIndex = SP.MCTruthIndex;
      P.TrackTruePDG = SP.PDG;
      P.TrackTrueE = SP.E;
      P.TrackTruePx = SP.Px;
      P.TrackTruePy = SP.Py;
      P.TrackTruePz = SP.Pz;
      P.TrackTrueEndE = SP.E;
      P.TrackTrueEndPx = SP.EndPx;
      P.TrackTrueEndPy = SP.EndPy;
      P.TrackTrueEndPz = SP.EndPz;
      P.TrackTrueModMomentum = SP.ModMomentum;
      P.TrackTrueEndModMomentum = SP.EndModMomentum;
      P.TrackTrueKE = SP.KE;
      P.TrackTrueEndKE = SP.EndKE;
      P.TrackTrueLength = SP.Travel;
      P.TrackTrueOrigin = SP.Origin;
      P.TrackTruthPurity = (double)maxhits/hits.size();
   }
   else P.HasTruth = false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P){

   std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = Assoc_TrackCalo->at(trk.key());
   std::vector<art::Ptr<anab::ParticleID>> trackPID = Assoc_TrackPID->at(trk.key());
   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

   PIDStore store = PIDCalc.GetPIDScores(trk, caloFromTrack, AlgScoresVec);
   P.Track_LLR_PID = store.LLR;
   P.Track_LLR_PID_Kaon = store.LLR_Kaon;
   P.Track_LLR_PID_Kaon_Partial = store.LLR_Kaon_Partial;
   P.MeandEdX_Plane0 = store.MeandEdX_Plane0;
   P.MeandEdX_Plane1 = store.MeandEdX_Plane1;
   P.MeandEdX_Plane2 = store.MeandEdX_Plane2;
   P.MeandEdX_ThreePlane = store.MeandEdX_3Plane;

   P.Track_Bragg_Pion = store.BraggWeighted_Pion;
   P.Track_Bragg_Muon = store.BraggWeighted_Muon;
   P.Track_Bragg_Proton = store.BraggWeighted_Proton;
   P.Track_Bragg_Kaon = store.BraggWeighted_Kaon;
   P.Track_Bragg_Sigma = store.BraggWeighted_Sigma;
   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetVertexData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

   std::vector<art::Ptr<recob::Vertex>> pfpVertex = Assoc_PFParticleVertex->at(pfp.key());

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

      geo::Point_t point = {vtx->position().X(),vtx->position().Y(),vtx->position().Z()};                
      geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

      TVector3 pos(vtx->position().X()+sce_corr.X(),vtx->position().Y()-sce_corr.Y(),vtx->position().Z()-sce_corr.Z());

      P.SetVertex(pos);
      P.Displacement = (pos-theData.RecoPrimaryVertex).Mag();

   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::SetIndices(std::vector<bool> IsSignal,std::vector<bool> IsSignalSigmaZero){
      
   bool ContainsSignal = std::find(IsSignal.begin(),IsSignal.end(),true) == IsSignal.end() 
                      || std::find(IsSignalSigmaZero.begin(),IsSignalSigmaZero.end(),true) == IsSignalSigmaZero.end();

   bool found_muon=false,found_decayproton=false,found_decaypion=false;

   for(size_t i_tr=0;i_tr<theData.TrackPrimaryDaughters.size();i_tr++){

      RecoParticle P = theData.TrackPrimaryDaughters.at(i_tr);

      if(!found_muon && abs(P.TrackTruePDG) == 13 && P.TrackTrueOrigin == 1){ 
         theData.TrueMuonIndex = P.Index;
         found_muon = true; 
      }

      if(ContainsSignal && !found_decayproton && P.TrackTruePDG == 2212 && P.TrackTrueOrigin == 2){
         theData.TrueDecayProtonIndex = P.Index;
         found_decayproton = true;
      }

      if(ContainsSignal && !found_decaypion && P.TrackTruePDG == -211 && P.TrackTrueOrigin == 2){
         theData.TrueDecayPionIndex = P.Index;
         found_decaypion = true;
      }
   }

   for(size_t i_sh=0;i_sh<theData.ShowerPrimaryDaughters.size();i_sh++){

      RecoParticle P = theData.ShowerPrimaryDaughters.at(i_sh);

      if(!found_muon && abs(P.TrackTruePDG) == 13 && P.TrackTrueOrigin == 1){ 
         theData.TrueMuonIndex = P.Index;
         found_muon = true; 
      }

      if(ContainsSignal && !found_decayproton && P.TrackTruePDG == 2212 && P.TrackTrueOrigin == 2){
         theData.TrueDecayProtonIndex = P.Index;
         found_decayproton = true;
      }

      if(ContainsSignal && !found_decaypion && P.TrackTruePDG == -211 && P.TrackTrueOrigin == 2){
         theData.TrueDecayPionIndex = P.Index;
         found_decaypion = true;
      }
   }

if(ContainsSignal && found_decayproton && found_decaypion) theData.GoodReco = true;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
