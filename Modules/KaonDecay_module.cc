////////////////////////////////////////////////////////////////////////
// Class:       KaonDecay
// Plugin Type: analyzer (art v3_03_01)
// File:        KaonDecay_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "cetlib_except/exception.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Local includes
#include "ubana/HyperonProduction/Alg/PIDManager.h"
//#include "ubana/HyperonProduction/Alg/LLRPIDHelper.h"
//#include "ubana/HyperonProduction/Alg/MeandEdXCalculator.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

//root includes
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace hyperon {
   class KaonDecay;
}


class hyperon::KaonDecay : public art::EDAnalyzer {
   public:
      explicit KaonDecay(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      KaonDecay(KaonDecay const&) = delete;
      KaonDecay(KaonDecay&&) = delete;
      KaonDecay& operator=(KaonDecay const&) = delete;
      KaonDecay& operator=(KaonDecay&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;

      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:

      const double _EPSILON_ = 0.00001;

      // Output trees
      TTree * OutputTree;

      // Basic event info
      unsigned int t_EventID;
      int t_run,t_subrun,t_event;

      std::vector<double> t_MCPDG;
      std::vector<double> t_MCMomentum;
      std::vector<double> t_MCEndMomentum;
      std::vector<double> t_MCStartX;
      std::vector<double> t_MCStartY;
      std::vector<double> t_MCStartZ;
      std::vector<double> t_MCEndX;
      std::vector<double> t_MCEndY;
      std::vector<double> t_MCEndZ;
      std::vector<double> t_MCRange;
      std::vector<double> t_MCStartPX;
      std::vector<double> t_MCStartPY;
      std::vector<double> t_MCStartPZ;
      std::vector<double> t_MCEndPX;
      std::vector<double> t_MCEndPY;
      std::vector<double> t_MCEndPZ;
      std::vector<std::string> t_MCEndProcess;
      std::vector<std::vector<int>> t_MCDaughterPDG_v;
      std::vector<std::vector<double>> t_MCDaughterMomentum_v; 

      std::vector<double> t_TrackLength;
      std::vector<double> t_TrackKaonMomentum;
      std::vector<double> t_TrackMuonMomentum;
      std::vector<double> t_TrackMeandEdX;
      std::vector<int> t_TrackTruePDG;
      std::vector<int> t_TrackTrueOrigin;      
      std::vector<double> t_TrackStart_X;
      std::vector<double> t_TrackStart_Y;
      std::vector<double> t_TrackStart_Z;
      std::vector<double> t_TrackEnd_X;
      std::vector<double> t_TrackEnd_Y;
      std::vector<double> t_TrackEnd_Z;
      std::vector<std::vector<float>> t_TrackdEdX_Plane0_v;
      std::vector<std::vector<float>> t_TrackResidualRange_Plane0_v;
      std::vector<std::vector<float>> t_TrackPitch_Plane0_v;
      std::vector<std::vector<float>> t_TrackdEdX_Corrected_Plane0_v;
      std::vector<std::vector<float>> t_TrackdEdX_Plane1_v;
      std::vector<std::vector<float>> t_TrackResidualRange_Plane1_v;
      std::vector<std::vector<float>> t_TrackPitch_Plane1_v;
      std::vector<std::vector<float>> t_TrackdEdX_Corrected_Plane1_v;
      std::vector<std::vector<float>> t_TrackdEdX_Plane2_v;
      std::vector<std::vector<float>> t_TrackResidualRange_Plane2_v;
      std::vector<std::vector<float>> t_TrackPitch_Plane2_v;
      std::vector<std::vector<float>> t_TrackdEdX_Corrected_Plane2_v;
      std::vector<double> t_TrackBraggPID_Plane0;
      std::vector<double> t_TrackBraggPID_Plane1;
      std::vector<double> t_TrackBraggPID_Plane2;
      std::vector<double> t_TrackBraggPID_3Plane;
      std::vector<double> t_TrackBraggShiftPID_Plane0;
      std::vector<double> t_TrackBraggShiftPID_Plane1;
      std::vector<double> t_TrackBraggShiftPID_Plane2;
      std::vector<double> t_TrackLLRPID;
      std::vector<double> t_TrackLLRPID_Kaon;
      std::vector<double> t_TrackLLRPID_Kaon_Partial;
 
      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      std::string f_G4Label;
      std::string f_TrackLabel;
      std::string f_HitLabel;
      std::string f_TrackHitAssn;
      std::string f_TrackCaloAssn;
      std::string f_HitTruthAssn;
      std::string f_PIDLabel;     
 
      bool f_Debug = false;

      ///////////////////////
      //      Objects      //
      ///////////////////////

      //LLRPIDHelper LLRPIDCalc;
      //MeandEdXCalculator dEdXCalc;
      PIDManager PIDCalc;
      trkf::TrackMomentumCalculator trkm{0};
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::KaonDecay::KaonDecay(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_G4Label(p.get<std::string>("G4Label","largeant")),
   f_TrackLabel(p.get<std::string>("TrackLabel","pandora")),
   f_HitLabel(p.get<std::string>("HitLabel","gaushit")),
   f_TrackHitAssn(p.get<std::string>("TrackHitAssnLabel","pandora")),
   f_TrackCaloAssn(p.get<std::string>("TrackCaloAssnLabel","pandoracali")),
   f_HitTruthAssn(p.get<std::string>("HitTruthAssnLabel","gaushitTruthMatch")),
   f_PIDLabel(p.get<std::string>("PIDLabel","pandorapid")),
   f_Debug(p.get<bool>("Debug",false)),
   //LLRPIDCalc(),
   //dEdXCalc()
   PIDCalc()
{
}

void hyperon::KaonDecay::analyze(art::Event const& e)
{
   if(f_Debug) std::cout << "New Event" << std::endl;

   t_MCPDG.clear();
   t_MCMomentum.clear();
   t_MCEndMomentum.clear();
   t_MCStartX.clear();
   t_MCStartY.clear();
   t_MCStartZ.clear();
   t_MCEndX.clear();
   t_MCEndY.clear();
   t_MCEndZ.clear();
   t_MCRange.clear();
   t_MCStartPX.clear();
   t_MCStartPY.clear();
   t_MCStartPZ.clear();
   t_MCEndPX.clear();
   t_MCEndPY.clear();
   t_MCEndPZ.clear();
   t_MCEndProcess.clear();
   t_MCDaughterPDG_v.clear();
   t_MCDaughterMomentum_v.clear(); 
   t_TrackLength.clear();
   t_TrackKaonMomentum.clear();
   t_TrackMuonMomentum.clear();
   t_TrackMeandEdX.clear();
   t_TrackTruePDG.clear();
   t_TrackTrueOrigin.clear();     
   t_TrackStart_X.clear();
   t_TrackStart_Y.clear();
   t_TrackStart_Z.clear();
   t_TrackEnd_X.clear();
   t_TrackEnd_Y.clear();
   t_TrackEnd_Z.clear();
   t_TrackdEdX_Plane0_v.clear();
   t_TrackResidualRange_Plane0_v.clear();
   t_TrackPitch_Plane0_v.clear();
   t_TrackdEdX_Corrected_Plane0_v.clear();
   t_TrackdEdX_Plane1_v.clear();
   t_TrackResidualRange_Plane1_v.clear();
   t_TrackPitch_Plane1_v.clear();
   t_TrackdEdX_Corrected_Plane1_v.clear();
   t_TrackdEdX_Plane2_v.clear();
   t_TrackResidualRange_Plane2_v.clear();
   t_TrackPitch_Plane2_v.clear();
   t_TrackdEdX_Corrected_Plane2_v.clear();
   t_TrackBraggPID_Plane0.clear();
   t_TrackBraggPID_Plane1.clear();
   t_TrackBraggPID_Plane2.clear();
   t_TrackBraggShiftPID_Plane0.clear();
   t_TrackBraggShiftPID_Plane1.clear();
   t_TrackBraggShiftPID_Plane2.clear();
   t_TrackBraggPID_3Plane.clear();
   t_TrackLLRPID.clear();
   t_TrackLLRPID_Kaon.clear();
   t_TrackLLRPID_Kaon_Partial.clear();
   
   // General Event Info

   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   // Get G4 info
   art::Handle<std::vector<simb::MCParticle>> Handle_G4;
   std::vector<art::Ptr<simb::MCParticle>> Vect_G4;

   if(!e.getByLabel(f_G4Label,Handle_G4))  
      throw cet::exception("KaonDecay") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_G4,Handle_G4);

   // Create map between particle ID and Particles
   std::map<int,art::Ptr<simb::MCParticle>> partByID;
   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
      partByID.insert(std::make_pair(g4p->TrackId(),g4p));

   // Store the IDs of the decay products
   std::vector<int> daughter_IDs;

   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){

      if(g4p->Mother() != 0)  continue;

      t_MCPDG.push_back(g4p->PdgCode());        
      t_MCMomentum.push_back(sqrt(g4p->Px()*g4p->Px()+g4p->Py()*g4p->Py()+g4p->Pz()*g4p->Pz()));
      t_MCEndMomentum.push_back(sqrt(g4p->EndPx()*g4p->EndPx()+g4p->EndPy()*g4p->EndPy()+g4p->EndPz()*g4p->EndPz())); 
      t_MCStartX.push_back(g4p->Position().X());
      t_MCStartY.push_back(g4p->Position().Y());
      t_MCStartZ.push_back(g4p->Position().Z());
      t_MCEndX.push_back(g4p->EndPosition().X());
      t_MCEndY.push_back(g4p->EndPosition().Y());
      t_MCEndZ.push_back(g4p->EndPosition().Z());                
      t_MCStartPX.push_back(g4p->Momentum().X());
      t_MCStartPY.push_back(g4p->Momentum().Y());
      t_MCStartPZ.push_back(g4p->Momentum().Z());
      t_MCEndPX.push_back(g4p->EndMomentum().X());
      t_MCEndPY.push_back(g4p->EndMomentum().Y());
      t_MCEndPZ.push_back(g4p->EndMomentum().Z());                
      t_MCEndProcess.push_back(g4p->EndProcess());

      double range = 0.0;
      for(size_t i_p=1;i_p<g4p->NumberTrajectoryPoints();i_p++)
         range += TVector3(g4p->Position(i_p).X()-g4p->Position(i_p-1).X(),g4p->Position(i_p).Y()-g4p->Position(i_p-1).Y(),g4p->Position(i_p).Z()-g4p->Position(i_p-1).Z()).Mag();

      t_MCRange.push_back(range);

      std::vector<int> daughter_pdgs;
      std::vector<double> daughter_momenta;
      for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

         if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;
         art::Ptr<simb::MCParticle> part = partByID[g4p->Daughter(i_d)];

         // calculate the distance from the end position of the kaon to the start position of this particle
         TLorentzVector d = g4p->EndPosition() - part->Position();
         double dist = (TVector3(d.X(),d.Y(),d.Z())).Mag(); 
         if(dist > _EPSILON_) continue;

         daughter_pdgs.push_back(part->PdgCode());
         daughter_momenta.push_back(sqrt(part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz()));
         daughter_IDs.push_back(part->TrackId());
      }
      t_MCDaughterPDG_v.push_back(daughter_pdgs);
      t_MCDaughterMomentum_v.push_back(daughter_momenta);      
   }

   // Get Reco info
   art::Handle<std::vector<recob::Track>> Handle_Tracks;
   std::vector<art::Ptr<recob::Track>> Vect_Tracks;
   art::Handle<std::vector<recob::Hit>> Handle_Hits;
   std::vector<art::Ptr<recob::Hit>> Vect_Hits;
     
   if(!e.getByLabel(f_TrackLabel,Handle_Tracks))  
      throw cet::exception("KaonDecay") << "No Track data product!" << std::endl;

   if(!e.getByLabel(f_HitLabel,Handle_Hits)) 
      throw cet::exception("KaonDecay") << "Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_Tracks,Handle_Tracks);
   art::fill_ptr_vector(Vect_Hits,Handle_Hits);

   art::FindManyP<recob::Hit>* Assoc_TrackHit = new art::FindManyP<recob::Hit>(Vect_Tracks,e,f_TrackHitAssn);
   art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo = new art::FindManyP<anab::Calorimetry>(Vect_Tracks,e,f_TrackCaloAssn);
   art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hits,e,f_HitTruthAssn);
   art::FindManyP<anab::ParticleID>* Assoc_TrackPID = new art::FindManyP<anab::ParticleID>(Vect_Tracks,e,f_PIDLabel);

   for(const art::Ptr<recob::Track> &trk : Vect_Tracks){

        // Truth match
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

        int origin = 3;
        int pdg = 0;

        if(matchedParticle != NULL){
        pdg = matchedParticle->PdgCode();
           if(matchedParticle->Mother() == 0) origin = 1;
           if(std::find(daughter_IDs.begin(),daughter_IDs.end(),matchedParticle->TrackId()) != daughter_IDs.end()) origin = 2;        
        }

        t_TrackLength.push_back(trk->Length());
        t_TrackTruePDG.push_back(pdg);
        t_TrackTrueOrigin.push_back(origin);      
        t_TrackKaonMomentum.push_back((156.222*pow(trk->Length()-0.0426198,0.274777) + 1.5323*trk->Length())/1e3);
        t_TrackMuonMomentum.push_back(trkm.GetTrackMomentum(trk->Length(),13));
        t_TrackStart_X.push_back(trk->Start().X());
        t_TrackStart_Y.push_back(trk->Start().Y());
        t_TrackStart_Z.push_back(trk->Start().Z());
        t_TrackEnd_X.push_back(trk->End().X());
        t_TrackEnd_Y.push_back(trk->End().Y());
        t_TrackEnd_Z.push_back(trk->End().Z());

        std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = Assoc_TrackCalo->at(trk.key());

        // Get kaon bragg likelihood
        std::vector<art::Ptr<anab::ParticleID>> trackPID = Assoc_TrackPID->at(trk.key());
        std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
        double BraggPID_Plane0 = -1;
        double BraggPID_Plane1 = -1;
        double BraggPID_Plane2 = -1;
        double BraggShiftPID_Plane0 = -1;
        double BraggShiftPID_Plane1 = -1;
        double BraggShiftPID_Plane2 = -1;

        for(size_t i_algscore=0;i_algscore<AlgScoresVec.size();i_algscore++){
           anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
           if(AlgScore.fAssumedPdg == 321 && AlgScore.fAlgName=="BraggPeakLLH" && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
              if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 0) BraggPID_Plane0 = AlgScore.fValue;
              if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 1) BraggPID_Plane1 = AlgScore.fValue;
              if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 2) BraggPID_Plane2 = AlgScore.fValue;
           }
           if(AlgScore.fAssumedPdg == 321 && AlgScore.fAlgName=="BraggPeakLLH_shift" && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
              if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 0) BraggShiftPID_Plane0 = AlgScore.fValue;
              if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 1) BraggShiftPID_Plane1 = AlgScore.fValue;
              if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 2) BraggShiftPID_Plane2 = AlgScore.fValue;
           }
        }

        PIDStore store = PIDCalc.GetPIDs(trk,caloFromTrack,AlgScoresVec);
        t_TrackMeandEdX.push_back(store.MeandEdX_3Plane);

        t_TrackdEdX_Plane0_v.push_back(store.dEdX_Plane0);
        t_TrackPitch_Plane0_v.push_back(store.Pitch_Plane0);
        t_TrackResidualRange_Plane0_v.push_back(store.ResidualRange_Plane0);
        t_TrackdEdX_Corrected_Plane0_v.push_back(store.dEdX_Corrected_Plane0);
        t_TrackdEdX_Plane1_v.push_back(store.dEdX_Plane1);
        t_TrackPitch_Plane1_v.push_back(store.Pitch_Plane1);
        t_TrackResidualRange_Plane1_v.push_back(store.ResidualRange_Plane1);
        t_TrackdEdX_Corrected_Plane1_v.push_back(store.dEdX_Corrected_Plane1);
        t_TrackdEdX_Plane2_v.push_back(store.dEdX_Plane2);
        t_TrackPitch_Plane2_v.push_back(store.Pitch_Plane2);
        t_TrackResidualRange_Plane2_v.push_back(store.ResidualRange_Plane2);
        t_TrackdEdX_Corrected_Plane2_v.push_back(store.dEdX_Corrected_Plane2);

        t_TrackBraggPID_Plane0.push_back(BraggPID_Plane0);
        t_TrackBraggPID_Plane1.push_back(BraggPID_Plane1);
        t_TrackBraggPID_Plane2.push_back(BraggPID_Plane2);
        t_TrackBraggShiftPID_Plane0.push_back(BraggShiftPID_Plane0);
        t_TrackBraggShiftPID_Plane1.push_back(BraggShiftPID_Plane1);
        t_TrackBraggShiftPID_Plane2.push_back(BraggShiftPID_Plane2);

        double BraggPID_3Plane = BraggPID_Plane0*PIDCalc.PlaneWeight(trk,0) + BraggPID_Plane1*PIDCalc.PlaneWeight(trk,1) + BraggPID_Plane2*PIDCalc.PlaneWeight(trk,2);
        BraggPID_3Plane /= (PIDCalc.PlaneWeight(trk,0) + PIDCalc.PlaneWeight(trk,1) + PIDCalc.PlaneWeight(trk,2));

        t_TrackBraggPID_3Plane.push_back(BraggPID_3Plane);
        t_TrackLLRPID.push_back(store.LLR);
        t_TrackLLRPID_Kaon.push_back(store.LLR_Kaon);
        t_TrackLLRPID_Kaon_Partial.push_back(store.LLR_Kaon_Partial);
   }  
 
   OutputTree->Fill();
}

///////////////////////////////////////////////////////////////	

void hyperon::KaonDecay::beginJob(){

   if(f_Debug) std::cout << "Begin job" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //             Output Tree	           //
   //////////////////////////////////////////

   OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");

   OutputTree->Branch("EventID",&t_EventID);
   OutputTree->Branch("run",&t_run);
   OutputTree->Branch("subrun",&t_subrun);
   OutputTree->Branch("event",&t_event);

   OutputTree->Branch("MCPDG",&t_MCPDG);
   OutputTree->Branch("MCMomentum",&t_MCMomentum);
   OutputTree->Branch("MCEndMomentum",&t_MCEndMomentum);
   OutputTree->Branch("MCStartX",&t_MCStartX);
   OutputTree->Branch("MCStartY",&t_MCStartY);
   OutputTree->Branch("MCStartZ",&t_MCStartZ);
   OutputTree->Branch("MCEndX",&t_MCEndX);
   OutputTree->Branch("MCEndY",&t_MCEndY);
   OutputTree->Branch("MCEndZ",&t_MCEndZ);
   OutputTree->Branch("MCRange",&t_MCRange);
   OutputTree->Branch("MCStartPX",&t_MCStartPX);
   OutputTree->Branch("MCStartPY",&t_MCStartPY);
   OutputTree->Branch("MCStartPZ",&t_MCStartPZ);
   OutputTree->Branch("MCEndPX",&t_MCEndPX);
   OutputTree->Branch("MCEndPY",&t_MCEndPY);
   OutputTree->Branch("MCEndPZ",&t_MCEndPZ);
   OutputTree->Branch("MCEndProcess",&t_MCEndProcess);
   OutputTree->Branch("MCDaughterPDG_v",&t_MCDaughterPDG_v);
   OutputTree->Branch("MCDaughterMomentum_v",&t_MCDaughterMomentum_v);

   OutputTree->Branch("TrackLength",&t_TrackLength);
   OutputTree->Branch("TrackKaonMomentum",&t_TrackKaonMomentum);
   OutputTree->Branch("TrackMuonMomentum",&t_TrackMuonMomentum);
   OutputTree->Branch("TrackMeandEdX",&t_TrackMeandEdX);
   OutputTree->Branch("TrackTruePDG",&t_TrackTruePDG);
   OutputTree->Branch("TrackTrueOrigin",&t_TrackTrueOrigin);
   OutputTree->Branch("TrackStart_X",&t_TrackStart_X);
   OutputTree->Branch("TrackStart_Y",&t_TrackStart_Y);
   OutputTree->Branch("TrackStart_Z",&t_TrackStart_Z);
   OutputTree->Branch("TrackEnd_X",&t_TrackEnd_X);
   OutputTree->Branch("TrackEnd_Y",&t_TrackEnd_Y);
   OutputTree->Branch("TrackEnd_Z",&t_TrackEnd_Z);
/*
   OutputTree->Branch("TrackdEdX_Plane0_v",&t_TrackdEdX_Plane0_v);
   OutputTree->Branch("TrackPitch_Plane0_v",&t_TrackPitch_Plane0_v);
   OutputTree->Branch("TrackResidualRange_Plane0_v",&t_TrackResidualRange_Plane0_v);
   OutputTree->Branch("TrackdEdX_Corrected_Plane0_v",&t_TrackdEdX_Corrected_Plane0_v);
   OutputTree->Branch("TrackdEdX_Plane1_v",&t_TrackdEdX_Plane1_v);
   OutputTree->Branch("TrackPitch_Plane1_v",&t_TrackPitch_Plane1_v);
   OutputTree->Branch("TrackResidualRange_Plane1_v",&t_TrackResidualRange_Plane1_v);
   OutputTree->Branch("TrackdEdX_Corrected_Plane1_v",&t_TrackdEdX_Corrected_Plane1_v);
   OutputTree->Branch("TrackdEdX_Plane2_v",&t_TrackdEdX_Plane2_v);
   OutputTree->Branch("TrackPitch_Plane2_v",&t_TrackPitch_Plane2_v);
   OutputTree->Branch("TrackResidualRange_Plane2_v",&t_TrackResidualRange_Plane2_v);
   OutputTree->Branch("TrackdEdX_Corrected_Plane2_v",&t_TrackdEdX_Corrected_Plane2_v);
*/
   OutputTree->Branch("TrackBraggPID_Plane0",&t_TrackBraggPID_Plane0);   
   OutputTree->Branch("TrackBraggPID_Plane1",&t_TrackBraggPID_Plane1);
   OutputTree->Branch("TrackBraggPID_Plane2",&t_TrackBraggPID_Plane2);
   OutputTree->Branch("TrackBraggShiftPID_Plane0",&t_TrackBraggShiftPID_Plane0);   
   OutputTree->Branch("TrackBraggShiftPID_Plane1",&t_TrackBraggShiftPID_Plane1);
   OutputTree->Branch("TrackBraggShiftPID_Plane2",&t_TrackBraggShiftPID_Plane2);
   OutputTree->Branch("TrackBraggPID_3Plane",&t_TrackBraggPID_3Plane);
   OutputTree->Branch("TrackLLRPID",&t_TrackLLRPID);
   OutputTree->Branch("TrackLLRPID_Kaon",&t_TrackLLRPID_Kaon);
   OutputTree->Branch("TrackLLRPID_Kaon_Partial",&t_TrackLLRPID_Kaon_Partial); 

   if(f_Debug) std::cout << "Finished begin job" << std::endl;
}

void hyperon::KaonDecay::endJob()
{
}

void hyperon::KaonDecay::beginSubRun(const art::SubRun& sr)
{
}

void hyperon::KaonDecay::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(hyperon::KaonDecay)
