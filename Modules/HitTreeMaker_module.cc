////////////////////////////////////////////////////////////////////////
// Class:       HitTreeMaker
// Plugin Type: analyzer (art v3_03_01)
// File:        HitTreeMaker_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <string>

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// Root includes
#include "TTree.h"
#include "TVector3.h"

// Local includes
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"

namespace hyperon {
   class HitTreeMaker;
}

class hyperon::HitTreeMaker : public art::EDAnalyzer {
   public:
      explicit HitTreeMaker(fhicl::ParameterSet const& p);

      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      HitTreeMaker(HitTreeMaker const&) = delete;
      HitTreeMaker(HitTreeMaker&&) = delete;
      HitTreeMaker& operator=(HitTreeMaker const&) = delete;
      HitTreeMaker& operator=(HitTreeMaker&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;

      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:

      // General event info
      unsigned int t_EventID;
      int t_run,t_subrun,t_event;

      ///////////////////////////////
      //     Wire Signal Info      //
      ///////////////////////////////

      TTree *t_HitTree;

      // Track Start (in xyz and tick/channel space)
      std::vector<int> t_TrackStart_Channel_Plane0;
      std::vector<int> t_TrackStart_Time_Plane0;
      std::vector<int> t_TrackEnd_Channel_Plane0;
      std::vector<int> t_TrackEnd_Time_Plane0;

      std::vector<int> t_TrackStart_Channel_Plane1;
      std::vector<int> t_TrackStart_Time_Plane1;
      std::vector<int> t_TrackEnd_Channel_Plane1;
      std::vector<int> t_TrackEnd_Time_Plane1;

      std::vector<int> t_TrackStart_Channel_Plane2;
      std::vector<int> t_TrackStart_Time_Plane2;
      std::vector<int> t_TrackEnd_Channel_Plane2;
      std::vector<int> t_TrackEnd_Time_Plane2;

      std::vector<double> t_TrackStart_X;
      std::vector<double> t_TrackStart_Y;
      std::vector<double> t_TrackStart_Z;
      std::vector<double> t_TrackEnd_X;
      std::vector<double> t_TrackEnd_Y;
      std::vector<double> t_TrackEnd_Z;

      std::vector<double> t_TrackDir_X;
      std::vector<double> t_TrackDir_Y;
      std::vector<double> t_TrackDir_Z;
      std::vector<double> t_TrackEndDir_X;
      std::vector<double> t_TrackEndDir_Y;
      std::vector<double> t_TrackEndDir_Z;

      std::vector<double> t_TrackGrad_Plane0;
      std::vector<double> t_TrackGrad_Plane1;
      std::vector<double> t_TrackGrad_Plane2;
      std::vector<double> t_TrackEndGrad_Plane0;
      std::vector<double> t_TrackEndGrad_Plane1;
      std::vector<double> t_TrackEndGrad_Plane2;
      std::vector<double> t_TrackAngle_Plane0;
      std::vector<double> t_TrackAngle_Plane1;
      std::vector<double> t_TrackAngle_Plane2;
      std::vector<double> t_TrackEndAngle_Plane0;
      std::vector<double> t_TrackEndAngle_Plane1;
      std::vector<double> t_TrackEndAngle_Plane2;

      // Shower Start (in xyz and tick/channel space)
      std::vector<int> t_ShowerStart_Channel_Plane0;
      std::vector<int> t_ShowerStart_Time_Plane0;

      std::vector<int> t_ShowerStart_Channel_Plane1;
      std::vector<int> t_ShowerStart_Time_Plane1;

      std::vector<int> t_ShowerStart_Channel_Plane2;
      std::vector<int> t_ShowerStart_Time_Plane2;

      std::vector<double> t_ShowerStart_X;
      std::vector<double> t_ShowerStart_Y;
      std::vector<double> t_ShowerStart_Z;

      std::vector<double> t_ShowerDir_X;
      std::vector<double> t_ShowerDir_Y;
      std::vector<double> t_ShowerDir_Z;

      // Track Calo Start (in xyz and tick/channel space)
      std::vector<int> t_CaloStart_Channel_Plane0;
      std::vector<int> t_CaloStart_Time_Plane0;
      std::vector<double> t_CaloStart_X_Plane0;
      std::vector<double> t_CaloStart_Y_Plane0;
      std::vector<double> t_CaloStart_Z_Plane0;

      std::vector<int> t_CaloStart_Channel_Plane1;
      std::vector<int> t_CaloStart_Time_Plane1;
      std::vector<double> t_CaloStart_X_Plane1;
      std::vector<double> t_CaloStart_Y_Plane1;
      std::vector<double> t_CaloStart_Z_Plane1;

      std::vector<int> t_CaloStart_Channel_Plane2;
      std::vector<int> t_CaloStart_Time_Plane2;
      std::vector<double> t_CaloStart_X_Plane2;
      std::vector<double> t_CaloStart_Y_Plane2;
      std::vector<double> t_CaloStart_Z_Plane2;

      // All hits in the event
      std::vector<int> t_Hit_Channel_Plane0;
      std::vector<int> t_Hit_Tick_Plane0;
      std::vector<double> t_Hit_Amplitude_Plane0;     
      std::vector<int> t_Hit_StartTick_Plane0;
      std::vector<int> t_Hit_EndTick_Plane0;

      std::vector<int> t_Hit_Channel_Plane1;
      std::vector<int> t_Hit_Tick_Plane1;
      std::vector<double> t_Hit_Amplitude_Plane1;     
      std::vector<int> t_Hit_StartTick_Plane1;
      std::vector<int> t_Hit_EndTick_Plane1;

      std::vector<int> t_Hit_Channel_Plane2;
      std::vector<int> t_Hit_Tick_Plane2;
      std::vector<double> t_Hit_Amplitude_Plane2;     
      std::vector<int> t_Hit_StartTick_Plane2;
      std::vector<int> t_Hit_EndTick_Plane2;

      // Hits belonging to tracks
      std::vector<int> t_TrackHit_TrackIndex_Plane0;
      std::vector<int> t_TrackHit_Channel_Plane0;
      std::vector<int> t_TrackHit_Tick_Plane0;
      std::vector<double> t_TrackHit_Amplitude_Plane0;     
      std::vector<int> t_TrackHit_StartTick_Plane0;
      std::vector<int> t_TrackHit_EndTick_Plane0;

      std::vector<int> t_TrackHit_TrackIndex_Plane1;
      std::vector<int> t_TrackHit_Channel_Plane1;
      std::vector<int> t_TrackHit_Tick_Plane1;
      std::vector<double> t_TrackHit_Amplitude_Plane1;     
      std::vector<int> t_TrackHit_StartTick_Plane1;
      std::vector<int> t_TrackHit_EndTick_Plane1;

      std::vector<int> t_TrackHit_TrackIndex_Plane2;
      std::vector<int> t_TrackHit_Channel_Plane2;
      std::vector<int> t_TrackHit_Tick_Plane2;
      std::vector<double> t_TrackHit_Amplitude_Plane2;     
      std::vector<int> t_TrackHit_StartTick_Plane2;
      std::vector<int> t_TrackHit_EndTick_Plane2;

      // Hits belonging to showers
      std::vector<int> t_ShowerHit_ShowerIndex_Plane0;
      std::vector<int> t_ShowerHit_Channel_Plane0;
      std::vector<int> t_ShowerHit_Tick_Plane0;
      std::vector<double> t_ShowerHit_Amplitude_Plane0;     
      std::vector<int> t_ShowerHit_StartTick_Plane0;
      std::vector<int> t_ShowerHit_EndTick_Plane0;

      std::vector<int> t_ShowerHit_ShowerIndex_Plane1;
      std::vector<int> t_ShowerHit_Channel_Plane1;
      std::vector<int> t_ShowerHit_Tick_Plane1;
      std::vector<double> t_ShowerHit_Amplitude_Plane1;     
      std::vector<int> t_ShowerHit_StartTick_Plane1;
      std::vector<int> t_ShowerHit_EndTick_Plane1;

      std::vector<int> t_ShowerHit_ShowerIndex_Plane2;
      std::vector<int> t_ShowerHit_Channel_Plane2;
      std::vector<int> t_ShowerHit_Tick_Plane2;
      std::vector<double> t_ShowerHit_Amplitude_Plane2;     
      std::vector<int> t_ShowerHit_StartTick_Plane2;
      std::vector<int> t_ShowerHit_EndTick_Plane2;

      // All hits in the slice
      std::vector<int> t_SliceHit_Channel_Plane0;
      std::vector<int> t_SliceHit_Tick_Plane0;
      std::vector<double> t_SliceHit_Amplitude_Plane0;     
      std::vector<int> t_SliceHit_StartTick_Plane0;
      std::vector<int> t_SliceHit_EndTick_Plane0;

      std::vector<int> t_SliceHit_Channel_Plane1;
      std::vector<int> t_SliceHit_Tick_Plane1;
      std::vector<double> t_SliceHit_Amplitude_Plane1;     
      std::vector<int> t_SliceHit_StartTick_Plane1;
      std::vector<int> t_SliceHit_EndTick_Plane1;

      std::vector<int> t_SliceHit_Channel_Plane2;
      std::vector<int> t_SliceHit_Tick_Plane2;
      std::vector<double> t_SliceHit_Amplitude_Plane2;     
      std::vector<int> t_SliceHit_StartTick_Plane2;
      std::vector<int> t_SliceHit_EndTick_Plane2;

      // Hits in slice but no belonging to a track or shower
      std::vector<int> t_SpareSliceHit_Channel_Plane0;
      std::vector<int> t_SpareSliceHit_Tick_Plane0;
      std::vector<double> t_SpareSliceHit_Amplitude_Plane0;     
      std::vector<int> t_SpareSliceHit_StartTick_Plane0;
      std::vector<int> t_SpareSliceHit_EndTick_Plane0;

      std::vector<int> t_SpareSliceHit_Channel_Plane1;
      std::vector<int> t_SpareSliceHit_Tick_Plane1;
      std::vector<double> t_SpareSliceHit_Amplitude_Plane1;     
      std::vector<int> t_SpareSliceHit_StartTick_Plane1;
      std::vector<int> t_SpareSliceHit_EndTick_Plane1;

      std::vector<int> t_SpareSliceHit_Channel_Plane2;
      std::vector<int> t_SpareSliceHit_Tick_Plane2;
      std::vector<double> t_SpareSliceHit_Amplitude_Plane2;     
      std::vector<int> t_SpareSliceHit_StartTick_Plane2;
      std::vector<int> t_SpareSliceHit_EndTick_Plane2;

      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      bool fDebug;

      // Producer module labels
      std::string fSliceLabel;
      std::string fTrackLabel;
      std::string fShowerLabel;
      std::string fPFParticleLabel;
      std::string fCaloLabel;
      std::string fClusterLabel;
      std::string fClusterHitAssocLabel;
      std::string fHitLabel;
      std::string fSliceHitAssocLabel;
      std::string fPFParticleSliceAssoc;
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::HitTreeMaker::HitTreeMaker(fhicl::ParameterSet const& p)
   : EDAnalyzer{p}   // ,
   // More initializers here.
{
   fDebug = p.get<bool>("Debug","false");

   // Module labels
   fSliceLabel = p.get<std::string>("SliceLabel","pandora");
   fPFParticleLabel = p.get<std::string>("PFParticleLabel","pandora");
   fTrackLabel = p.get<std::string>("TrackLabel","pandora");
   fShowerLabel = p.get<std::string>("ShowerLabel","pandora");
   fCaloLabel = p.get<std::string>("CaloLabel","pandoracali");
   fClusterLabel = p.get<std::string>("ClusterLabel","pandora");
   fClusterHitAssocLabel = p.get<std::string>("ClusterHitAssocLabel","pandora");
   fHitLabel = p.get<std::string>("HitLabel","gaushit");
   fSliceHitAssocLabel = p.get<std::string>("SliceHitAssocLabel","pandora");
   fPFParticleSliceAssoc = p.get<std::string>("PFParticleSliceAssoc","pandora");
}

void hyperon::HitTreeMaker::analyze(art::Event const& e)
{
   if(fDebug) std::cout << "New Event" << std::endl;

   ////////////////////////////////
   //  Reset All of the Vectors  //
   ////////////////////////////////

   t_TrackStart_Channel_Plane0.clear();
   t_TrackStart_Time_Plane0.clear();
   t_TrackEnd_Channel_Plane0.clear();
   t_TrackEnd_Time_Plane0.clear();

   t_TrackStart_Channel_Plane1.clear();
   t_TrackStart_Time_Plane1.clear();
   t_TrackEnd_Channel_Plane1.clear();
   t_TrackEnd_Time_Plane1.clear();

   t_TrackStart_Channel_Plane2.clear();
   t_TrackStart_Time_Plane2.clear();
   t_TrackEnd_Channel_Plane2.clear();
   t_TrackEnd_Time_Plane2.clear();

   t_TrackStart_X.clear();
   t_TrackStart_Y.clear();
   t_TrackStart_Z.clear();
   t_TrackEnd_X.clear();
   t_TrackEnd_Y.clear();
   t_TrackEnd_Z.clear();

   t_TrackDir_X.clear();
   t_TrackDir_Y.clear();
   t_TrackDir_Z.clear();
   t_TrackEndDir_X.clear();
   t_TrackEndDir_Y.clear();
   t_TrackEndDir_Z.clear();

   t_TrackGrad_Plane0.clear();
   t_TrackGrad_Plane1.clear();
   t_TrackGrad_Plane2.clear();
   t_TrackEndGrad_Plane0.clear();
   t_TrackEndGrad_Plane1.clear();
   t_TrackEndGrad_Plane2.clear();
   t_TrackAngle_Plane0.clear();
   t_TrackAngle_Plane1.clear();
   t_TrackAngle_Plane2.clear();
   t_TrackEndAngle_Plane0.clear();
   t_TrackEndAngle_Plane1.clear();
   t_TrackEndAngle_Plane2.clear();

   t_CaloStart_Channel_Plane0.clear();
   t_CaloStart_Time_Plane0.clear();
   t_CaloStart_X_Plane0.clear();
   t_CaloStart_Y_Plane0.clear();
   t_CaloStart_Z_Plane0.clear();

   t_CaloStart_Channel_Plane1.clear();
   t_CaloStart_Time_Plane1.clear();
   t_CaloStart_X_Plane1.clear();
   t_CaloStart_Y_Plane1.clear();
   t_CaloStart_Z_Plane1.clear();

   t_CaloStart_Channel_Plane2.clear();
   t_CaloStart_Time_Plane2.clear();
   t_CaloStart_X_Plane2.clear();
   t_CaloStart_Y_Plane2.clear();
   t_CaloStart_Z_Plane2.clear();

   t_ShowerStart_Channel_Plane0.clear();
   t_ShowerStart_Time_Plane0.clear();

   t_ShowerStart_Channel_Plane1.clear();
   t_ShowerStart_Time_Plane1.clear();

   t_ShowerStart_Channel_Plane2.clear();
   t_ShowerStart_Time_Plane2.clear();

   t_ShowerStart_X.clear();
   t_ShowerStart_Y.clear();
   t_ShowerStart_Z.clear();

   t_ShowerDir_X.clear();
   t_ShowerDir_Y.clear();
   t_ShowerDir_Z.clear();

   t_Hit_Channel_Plane0.clear();
   t_Hit_Tick_Plane0.clear();
   t_Hit_Amplitude_Plane0.clear();
   t_Hit_StartTick_Plane0.clear();
   t_Hit_EndTick_Plane0.clear();

   t_Hit_Channel_Plane1.clear();
   t_Hit_Tick_Plane1.clear();
   t_Hit_Amplitude_Plane1.clear();
   t_Hit_StartTick_Plane1.clear();
   t_Hit_EndTick_Plane1.clear();

   t_Hit_Channel_Plane2.clear();
   t_Hit_Tick_Plane2.clear();
   t_Hit_Amplitude_Plane2.clear();
   t_Hit_StartTick_Plane2.clear();
   t_Hit_EndTick_Plane2.clear();

   t_TrackHit_TrackIndex_Plane0.clear();
   t_TrackHit_Channel_Plane0.clear();
   t_TrackHit_Tick_Plane0.clear();
   t_TrackHit_Amplitude_Plane0.clear();
   t_TrackHit_StartTick_Plane0.clear();
   t_TrackHit_EndTick_Plane0.clear();

   t_TrackHit_TrackIndex_Plane1.clear();
   t_TrackHit_Channel_Plane1.clear();
   t_TrackHit_Tick_Plane1.clear();
   t_TrackHit_Amplitude_Plane1.clear();
   t_TrackHit_StartTick_Plane1.clear();
   t_TrackHit_EndTick_Plane1.clear();

   t_TrackHit_TrackIndex_Plane2.clear();
   t_TrackHit_Channel_Plane2.clear();
   t_TrackHit_Tick_Plane2.clear();
   t_TrackHit_Amplitude_Plane2.clear();
   t_TrackHit_StartTick_Plane2.clear();
   t_TrackHit_EndTick_Plane2.clear();

   t_ShowerHit_ShowerIndex_Plane0.clear();
   t_ShowerHit_Channel_Plane0.clear();
   t_ShowerHit_Tick_Plane0.clear();
   t_ShowerHit_Amplitude_Plane0.clear();
   t_ShowerHit_StartTick_Plane0.clear();
   t_ShowerHit_EndTick_Plane0.clear();

   t_ShowerHit_ShowerIndex_Plane1.clear();
   t_ShowerHit_Channel_Plane1.clear();
   t_ShowerHit_Tick_Plane1.clear();
   t_ShowerHit_Amplitude_Plane1.clear();
   t_ShowerHit_StartTick_Plane1.clear();
   t_ShowerHit_EndTick_Plane1.clear();

   t_ShowerHit_ShowerIndex_Plane2.clear();
   t_ShowerHit_Channel_Plane2.clear();
   t_ShowerHit_Tick_Plane2.clear();
   t_ShowerHit_Amplitude_Plane2.clear();
   t_ShowerHit_StartTick_Plane2.clear();
   t_ShowerHit_EndTick_Plane2.clear();

   t_SliceHit_Channel_Plane0.clear();
   t_SliceHit_Tick_Plane0.clear();
   t_SliceHit_Amplitude_Plane0.clear();
   t_SliceHit_StartTick_Plane0.clear();
   t_SliceHit_EndTick_Plane0.clear();

   t_SliceHit_Channel_Plane1.clear();
   t_SliceHit_Tick_Plane1.clear();
   t_SliceHit_Amplitude_Plane1.clear();
   t_SliceHit_StartTick_Plane1.clear();
   t_SliceHit_EndTick_Plane1.clear();

   t_SliceHit_Channel_Plane2.clear();
   t_SliceHit_Tick_Plane2.clear();
   t_SliceHit_Amplitude_Plane2.clear();
   t_SliceHit_StartTick_Plane2.clear();
   t_SliceHit_EndTick_Plane2.clear();

   t_SpareSliceHit_Channel_Plane0.clear();
   t_SpareSliceHit_Tick_Plane0.clear();
   t_SpareSliceHit_Amplitude_Plane0.clear();
   t_SpareSliceHit_StartTick_Plane0.clear();
   t_SpareSliceHit_EndTick_Plane0.clear();

   t_SpareSliceHit_Channel_Plane1.clear();
   t_SpareSliceHit_Tick_Plane1.clear();
   t_SpareSliceHit_Amplitude_Plane1.clear();
   t_SpareSliceHit_StartTick_Plane1.clear();
   t_SpareSliceHit_EndTick_Plane1.clear();

   t_SpareSliceHit_Channel_Plane2.clear();
   t_SpareSliceHit_Tick_Plane2.clear();
   t_SpareSliceHit_Amplitude_Plane2.clear();
   t_SpareSliceHit_StartTick_Plane2.clear();
   t_SpareSliceHit_EndTick_Plane2.clear();

   ////////////////////////////
   //  Event ID Information  //
   ////////////////////////////

   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   /////////////////////////////////////////////////////////
   //  Obtain Start Positions of Tracks in the Hierarchy  //
   /////////////////////////////////////////////////////////

   if(fDebug) std::cout << "Getting Reco'd Particles" << std::endl;

   // Setup handles
   art::Handle<std::vector<recob::Slice>> sliceHandle;
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
   art::Handle<std::vector<recob::Track>> trackHandle;
   art::Handle<std::vector<recob::Shower>> showerHandle;
   art::Handle<std::vector<anab::Calorimetry>> caloHandle;
   art::Handle<std::vector<recob::Cluster>> clusterHandle;
   art::Handle<std::vector<recob::Hit>> hitHandle;

   std::vector<art::Ptr<recob::Slice>> sliceVect;
   std::vector<art::Ptr<recob::PFParticle>> pfparticleVect;
   std::vector<art::Ptr<recob::Track>> trackVect;
   std::vector<art::Ptr<recob::Shower>> showerVect;
   std::vector<art::Ptr<anab::Calorimetry>> caloVect;
   std::vector<art::Ptr<recob::Cluster>> clusterVect;
   std::vector<art::Ptr<recob::Hit>> hitVect;

   // Fill Slice vector
   if(e.getByLabel(fSliceLabel,sliceHandle))
      art::fill_ptr_vector(sliceVect,sliceHandle);

   // Fill PFP vector
   if(e.getByLabel(fPFParticleLabel,pfparticleHandle))
      art::fill_ptr_vector(pfparticleVect,pfparticleHandle);

   // If PFP vector is empty, add blank entry to tree and go to next event
   if(!pfparticleVect.size()) {
      t_HitTree->Fill();
      return;
   }

   // Fill track vector
   if(e.getByLabel(fTrackLabel,trackHandle)) art::fill_ptr_vector(trackVect,trackHandle);
   else std::cout << "Track handle not setup" << std::endl;

   // Fill shower vector
   if(e.getByLabel(fShowerLabel,showerHandle)) art::fill_ptr_vector(showerVect,showerHandle);
   else std::cout << "Shower handle not setup" << std::endl;

   // Fill cluster vector
   if(e.getByLabel(fClusterLabel,clusterHandle)) art::fill_ptr_vector(clusterVect,clusterHandle);
   else std::cout << "Cluster handle not setup" << std::endl;

   // Fill hit vector
   if(e.getByLabel(fHitLabel,hitHandle)) art::fill_ptr_vector(hitVect,hitHandle);
   else std::cout << "Hit handle not setup" << std::endl;

   // Setup Assns

   // Tracks Assoc with PFPs
   art::FindManyP<recob::Track> trackAssoc(pfparticleVect,e,fTrackLabel);

   // Showers Assoc with PFPs
   art::FindManyP<recob::Shower> showerAssoc(pfparticleVect,e,fShowerLabel);

   // Calo Assoc to Tracks
   art::FindManyP<anab::Calorimetry> caloTrackAssoc(trackVect,e,fCaloLabel);

   // Cluster Assoc to PFPs
   art::FindManyP<recob::Cluster> clusterPFParticleAssoc(pfparticleVect,e,fClusterLabel);

   // Hit Assoc to Clusters
   art::FindManyP<recob::Hit> hitClusterAssoc(clusterVect,e,fClusterHitAssocLabel);

   // Slice Assoc to PFPs
   art::FindManyP<recob::Slice> slicePFParticleAssoc(pfparticleVect,e,fPFParticleSliceAssoc);

   // Hit Assoc to Slices
   art::FindManyP<recob::Hit> hitSliceAssoc(sliceVect,e,fSliceHitAssocLabel);

   // Go through the list of pandora PFP's, find the reconstructed neutrino
   size_t neutrinoID = 99999;
   std::vector<art::Ptr<recob::Slice>> nu_slice;
   for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){
      if((pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12))){
         neutrinoID = pfp->Self();
         nu_slice = slicePFParticleAssoc.at(pfp.key()); 
      }
   }

   std::vector<int> trackhit_keys;
   std::vector<int> showerhit_keys;
   std::vector<size_t> PFP_IDs;

   for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){

      // Only store information from particles in the hierarchy
      if(pfp->Parent() != neutrinoID && std::find(PFP_IDs.begin(),PFP_IDs.end(),pfp->Parent()) == PFP_IDs.end()) continue;
      //if(pfp->Parent() != neutrinoID) continue;

      std::vector<art::Ptr<recob::Track>> pfpTracks = trackAssoc.at(pfp.key());
      std::vector<art::Ptr<recob::Shower>> pfpShowers = showerAssoc.at(pfp.key());
      std::vector<art::Ptr<recob::Cluster>> clusterFromPFP = clusterPFParticleAssoc.at(pfp.key());
      PFP_IDs.push_back(pfp->Self());

      if(pfpTracks.size() == 1){

         art::Ptr<recob::Track> trk = pfpTracks.at(0);

         // Setup Calo Assn for next couple of steps	
         std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = caloTrackAssoc.at(trk.key());

         // Get Track Start point
         TVector3 TrackStart(trk->Start().X(),trk->Start().Y(),trk->Start().Z());
         TVector3 TrackEnd(trk->End().X(),trk->End().Y(),trk->End().Z());
         TVector3 TrackDir(trk->StartDirection().X(),trk->StartDirection().Y(),trk->StartDirection().Z());
         TVector3 TrackEndDir(trk->EndDirection().X(),trk->EndDirection().Y(),trk->EndDirection().Z());

         t_TrackStart_Channel_Plane0.push_back(U_wire(TrackStart)); 
         t_TrackStart_Time_Plane0.push_back(tick(TrackStart)); 
         t_TrackEnd_Channel_Plane0.push_back(U_wire(TrackEnd)); 
         t_TrackEnd_Time_Plane0.push_back(tick(TrackEnd)); 

         t_TrackStart_Channel_Plane1.push_back(V_wire(TrackStart)); 
         t_TrackStart_Time_Plane1.push_back(tick(TrackStart)); 
         t_TrackEnd_Channel_Plane1.push_back(V_wire(TrackEnd)); 
         t_TrackEnd_Time_Plane1.push_back(tick(TrackEnd)); 

         t_TrackStart_Channel_Plane2.push_back(Y_wire(TrackStart)); 
         t_TrackStart_Time_Plane2.push_back(tick(TrackStart)); 
         t_TrackEnd_Channel_Plane2.push_back(Y_wire(TrackEnd)); 
         t_TrackEnd_Time_Plane2.push_back(tick(TrackEnd)); 

         t_TrackStart_X.push_back(TrackStart.X());
         t_TrackStart_Y.push_back(TrackStart.Y());
         t_TrackStart_Z.push_back(TrackStart.Z());
         t_TrackEnd_X.push_back(TrackEnd.X());
         t_TrackEnd_Y.push_back(TrackEnd.Y());
         t_TrackEnd_Z.push_back(TrackEnd.Z());

         t_TrackDir_X.push_back(TrackDir.X());
         t_TrackDir_Y.push_back(TrackDir.Y());
         t_TrackDir_Z.push_back(TrackDir.Z());
         t_TrackEndDir_X.push_back(TrackEndDir.X());
         t_TrackEndDir_Y.push_back(TrackEndDir.Y());
         t_TrackEndDir_Z.push_back(TrackEndDir.Z());
        
         t_TrackGrad_Plane0.push_back(1.0/dUdt(TrackDir));
         t_TrackGrad_Plane1.push_back(1.0/dVdt(TrackDir));
         t_TrackGrad_Plane2.push_back(1.0/dYdt(TrackDir));
         t_TrackEndGrad_Plane0.push_back(1.0/dUdt(TrackEndDir));
         t_TrackEndGrad_Plane1.push_back(1.0/dVdt(TrackEndDir));
         t_TrackEndGrad_Plane2.push_back(1.0/dYdt(TrackEndDir));
         t_TrackAngle_Plane0.push_back(AngleU(TrackDir));
         t_TrackAngle_Plane1.push_back(AngleV(TrackDir));
         t_TrackAngle_Plane2.push_back(AngleY(TrackDir));
         t_TrackEndAngle_Plane0.push_back(AngleU(TrackEndDir));
         t_TrackEndAngle_Plane1.push_back(AngleV(TrackEndDir));
         t_TrackEndAngle_Plane2.push_back(AngleY(TrackEndDir));

         // Get Calo Start Point

         int U = -1000;
         int V = -1000;
         int Y = -1000;
         int tick_Plane0 = -1000;                                        
         int tick_Plane1 = -1000;                                        
         int tick_Plane2 = -1000;                                        

         TVector3 CaloStart_Plane0(-1000,-1000,-1000);
         TVector3 CaloStart_Plane1(-1000,-1000,-1000);
         TVector3 CaloStart_Plane2(-1000,-1000,-1000);

         for(size_t i_plane=0;i_plane<caloFromTrack.size();i_plane++){

            auto thisPlaneCalo = caloFromTrack.at(i_plane);
            int planeno = thisPlaneCalo->PlaneID().Plane;

            // Get the last point in the calo object
            if(thisPlaneCalo->XYZ().size()){

               // Start of the track is the last point in the Calo point list
               TVector3 CaloStart(thisPlaneCalo->XYZ().at(thisPlaneCalo->XYZ().size()-1).X(),thisPlaneCalo->XYZ().at(thisPlaneCalo->XYZ().size()-1).Y(),thisPlaneCalo->XYZ().at(thisPlaneCalo->XYZ().size()-1).Z());

               if(planeno == 0){
                  U = U_wire(CaloStart);
                  tick_Plane0 = tick(CaloStart);
                  CaloStart_Plane0 = CaloStart;
               } 
               if(planeno == 1){
                  V = V_wire(CaloStart);
                  tick_Plane1 = tick(CaloStart);
                  CaloStart_Plane1 = CaloStart;
               }
               if(planeno == 2){
                  Y = Y_wire(CaloStart);
                  tick_Plane2 = tick(CaloStart);
                  CaloStart_Plane2 = CaloStart;
               } 
            }
         }              

         t_CaloStart_Channel_Plane0.push_back(U);
         t_CaloStart_Time_Plane0.push_back(tick_Plane0);
         t_CaloStart_X_Plane0.push_back(CaloStart_Plane0.X());
         t_CaloStart_Y_Plane0.push_back(CaloStart_Plane0.Y());
         t_CaloStart_Z_Plane0.push_back(CaloStart_Plane0.Z());

         t_CaloStart_Channel_Plane1.push_back(V);
         t_CaloStart_Time_Plane1.push_back(tick_Plane1);
         t_CaloStart_X_Plane1.push_back(CaloStart_Plane1.X());
         t_CaloStart_Y_Plane1.push_back(CaloStart_Plane1.Y());
         t_CaloStart_Z_Plane1.push_back(CaloStart_Plane1.Z());

         t_CaloStart_Channel_Plane2.push_back(Y);
         t_CaloStart_Time_Plane2.push_back(tick_Plane2);
         t_CaloStart_X_Plane2.push_back(CaloStart_Plane2.X());
         t_CaloStart_Y_Plane2.push_back(CaloStart_Plane2.Y());
         t_CaloStart_Z_Plane2.push_back(CaloStart_Plane2.Z());

         for(size_t i_c=0;i_c<clusterFromPFP.size();i_c++){

            std::vector<art::Ptr<recob::Hit>> hitFromCluster = hitClusterAssoc.at(clusterFromPFP.at(i_c).key());

            for(size_t i_h=0;i_h<hitFromCluster.size();i_h++){
               art::Ptr<recob::Hit> hit = hitFromCluster.at(i_h);

               trackhit_keys.push_back(hit.key());

               if(hit->View() == 0){
                  t_TrackHit_TrackIndex_Plane0.push_back(t_TrackStart_Channel_Plane0.size()-1);
                  t_TrackHit_Channel_Plane0.push_back(hit->Channel());
                  t_TrackHit_Tick_Plane0.push_back(hit->PeakTime());
                  t_TrackHit_Amplitude_Plane0.push_back(hit->PeakAmplitude());
                  t_TrackHit_StartTick_Plane0.push_back(hit->StartTick());
                  t_TrackHit_EndTick_Plane0.push_back(hit->EndTick());
               }
               if(hit->View() == 1){
                  t_TrackHit_TrackIndex_Plane1.push_back(t_TrackStart_Channel_Plane1.size()-1);
                  t_TrackHit_Channel_Plane1.push_back(hit->Channel());
                  t_TrackHit_Tick_Plane1.push_back(hit->PeakTime());
                  t_TrackHit_Amplitude_Plane1.push_back(hit->PeakAmplitude());
                  t_TrackHit_StartTick_Plane1.push_back(hit->StartTick());
                  t_TrackHit_EndTick_Plane1.push_back(hit->EndTick());
               }
               if(hit->View() == 2){
                  t_TrackHit_TrackIndex_Plane2.push_back(t_TrackStart_Channel_Plane2.size()-1);
                  t_TrackHit_Channel_Plane2.push_back(hit->Channel());
                  t_TrackHit_Tick_Plane2.push_back(hit->PeakTime());
                  t_TrackHit_Amplitude_Plane2.push_back(hit->PeakAmplitude());
                  t_TrackHit_StartTick_Plane2.push_back(hit->StartTick());
                  t_TrackHit_EndTick_Plane2.push_back(hit->EndTick());
               }
            }
         }
      } //pfpTracks.size() == 1 

      if(pfpShowers.size() == 1){

         art::Ptr<recob::Shower> shwr = pfpShowers.at(0);

         // Get Shower Start point
         TVector3 ShowerStart(shwr->ShowerStart().X(),shwr->ShowerStart().Y(),shwr->ShowerStart().Z());
         TVector3 ShowerDir(shwr->Direction().X(),shwr->Direction().Y(),shwr->Direction().Z());

         t_ShowerStart_Channel_Plane0.push_back(U_wire(ShowerStart)); 
         t_ShowerStart_Time_Plane0.push_back(tick(ShowerStart)); 

         t_ShowerStart_Channel_Plane1.push_back(V_wire(ShowerStart)); 
         t_ShowerStart_Time_Plane1.push_back(tick(ShowerStart)); 

         t_ShowerStart_Channel_Plane2.push_back(Y_wire(ShowerStart)); 
         t_ShowerStart_Time_Plane2.push_back(tick(ShowerStart)); 

         t_ShowerStart_X.push_back(ShowerStart.X());
         t_ShowerStart_Y.push_back(ShowerStart.Y());
         t_ShowerStart_Z.push_back(ShowerStart.Z());

         t_ShowerDir_X.push_back(ShowerDir.X());
         t_ShowerDir_Y.push_back(ShowerDir.Y());
         t_ShowerDir_Z.push_back(ShowerDir.Z());

         for(size_t i_c=0;i_c<clusterFromPFP.size();i_c++){

            std::vector<art::Ptr<recob::Hit>> hitFromCluster = hitClusterAssoc.at(clusterFromPFP.at(i_c).key());

            for(size_t i_h=0;i_h<hitFromCluster.size();i_h++){
               art::Ptr<recob::Hit> hit = hitFromCluster.at(i_h);

               showerhit_keys.push_back(hit.key());

               if(hit->View() == 0){
                  t_ShowerHit_ShowerIndex_Plane2.push_back(t_ShowerStart_Channel_Plane0.size()-1);
                  t_ShowerHit_Channel_Plane0.push_back(hit->Channel());
                  t_ShowerHit_Tick_Plane0.push_back(hit->PeakTime());
                  t_ShowerHit_Amplitude_Plane0.push_back(hit->PeakAmplitude());
                  t_ShowerHit_StartTick_Plane0.push_back(hit->StartTick());
                  t_ShowerHit_EndTick_Plane0.push_back(hit->EndTick());
               }
               if(hit->View() == 1){
                  t_ShowerHit_ShowerIndex_Plane2.push_back(t_ShowerStart_Channel_Plane1.size()-1);
                  t_ShowerHit_Channel_Plane1.push_back(hit->Channel());
                  t_ShowerHit_Tick_Plane1.push_back(hit->PeakTime());
                  t_ShowerHit_Amplitude_Plane1.push_back(hit->PeakAmplitude());
                  t_ShowerHit_StartTick_Plane1.push_back(hit->StartTick());
                  t_ShowerHit_EndTick_Plane1.push_back(hit->EndTick());
               }
               if(hit->View() == 2){
                  t_ShowerHit_ShowerIndex_Plane2.push_back(t_ShowerStart_Channel_Plane2.size()-1);
                  t_ShowerHit_Channel_Plane2.push_back(hit->Channel());
                  t_ShowerHit_Tick_Plane2.push_back(hit->PeakTime());
                  t_ShowerHit_Amplitude_Plane2.push_back(hit->PeakAmplitude());
                  t_ShowerHit_StartTick_Plane2.push_back(hit->StartTick());
                  t_ShowerHit_EndTick_Plane2.push_back(hit->EndTick());
               }
            }
         }
      } //pfpShowers.size() == 1
   } //end of PFP loop

   // Get the hits in the neutrino slice
   for(size_t i_s=0;i_s<nu_slice.size();i_s++){    
      std::vector<art::Ptr<recob::Hit>> sliceHits = hitSliceAssoc.at(nu_slice.at(i_s).key());
      for(size_t i_h=0;i_h<sliceHits.size();i_h++){

         art::Ptr<recob::Hit> hit = sliceHits.at(i_h);
        
         bool spare_hit = std::find(trackhit_keys.begin(),trackhit_keys.end(),hit.key()) == trackhit_keys.end() && std::find(showerhit_keys.begin(),showerhit_keys.end(),hit.key()) == showerhit_keys.end();

         if(hit->View() == 0){
            t_SliceHit_Channel_Plane0.push_back(hit->Channel());
            t_SliceHit_Tick_Plane0.push_back(hit->PeakTime());
            t_SliceHit_Amplitude_Plane0.push_back(hit->PeakAmplitude());
            t_SliceHit_StartTick_Plane0.push_back(hit->StartTick());
            t_SliceHit_EndTick_Plane0.push_back(hit->EndTick());
            if(spare_hit){
               t_SpareSliceHit_Channel_Plane0.push_back(hit->Channel());
               t_SpareSliceHit_Tick_Plane0.push_back(hit->PeakTime());
               t_SpareSliceHit_Amplitude_Plane0.push_back(hit->PeakAmplitude());
               t_SpareSliceHit_StartTick_Plane0.push_back(hit->StartTick());
               t_SpareSliceHit_EndTick_Plane0.push_back(hit->EndTick());
            }
         }
         if(hit->View() == 1){
            t_SliceHit_Channel_Plane1.push_back(hit->Channel());
            t_SliceHit_Tick_Plane1.push_back(hit->PeakTime());
            t_SliceHit_Amplitude_Plane1.push_back(hit->PeakAmplitude());
            t_SliceHit_StartTick_Plane1.push_back(hit->StartTick());
            t_SliceHit_EndTick_Plane1.push_back(hit->EndTick());
            if(spare_hit){
               t_SpareSliceHit_Channel_Plane1.push_back(hit->Channel());
               t_SpareSliceHit_Tick_Plane1.push_back(hit->PeakTime());
               t_SpareSliceHit_Amplitude_Plane1.push_back(hit->PeakAmplitude());
               t_SpareSliceHit_StartTick_Plane1.push_back(hit->StartTick());
               t_SpareSliceHit_EndTick_Plane1.push_back(hit->EndTick());
            }
         }
         if(hit->View() == 2){
            t_SliceHit_Channel_Plane2.push_back(hit->Channel());
            t_SliceHit_Tick_Plane2.push_back(hit->PeakTime());
            t_SliceHit_Amplitude_Plane2.push_back(hit->PeakAmplitude());
            t_SliceHit_StartTick_Plane2.push_back(hit->StartTick());
            t_SliceHit_EndTick_Plane2.push_back(hit->EndTick());
            if(spare_hit){
               t_SpareSliceHit_Channel_Plane2.push_back(hit->Channel());
               t_SpareSliceHit_Tick_Plane2.push_back(hit->PeakTime());
               t_SpareSliceHit_Amplitude_Plane2.push_back(hit->PeakAmplitude());
               t_SpareSliceHit_StartTick_Plane2.push_back(hit->StartTick());
               t_SpareSliceHit_EndTick_Plane2.push_back(hit->EndTick());
            }
         }                     
      }
   }

   /////////////////////////////
   //     Obtain Hit Info     //
   /////////////////////////////

   // Iterate through all of the wires, record signal at every tick with nonzero signal
   for(const art::Ptr<recob::Hit> &hit : hitVect){

      if(hit->View() == 0){
         t_Hit_Channel_Plane0.push_back(hit->Channel());
         t_Hit_Tick_Plane0.push_back(hit->PeakTime());
         t_Hit_Amplitude_Plane0.push_back(hit->PeakAmplitude());
         t_Hit_StartTick_Plane0.push_back(hit->StartTick());
         t_Hit_EndTick_Plane0.push_back(hit->EndTick());
      }

      if(hit->View() == 1){
         t_Hit_Channel_Plane1.push_back(hit->Channel());
         t_Hit_Tick_Plane1.push_back(hit->PeakTime());
         t_Hit_Amplitude_Plane1.push_back(hit->PeakAmplitude());
         t_Hit_StartTick_Plane1.push_back(hit->StartTick());
         t_Hit_EndTick_Plane1.push_back(hit->EndTick());
      }

      if(hit->View() == 2){
         t_Hit_Channel_Plane2.push_back(hit->Channel());
         t_Hit_Tick_Plane2.push_back(hit->PeakTime());
         t_Hit_Amplitude_Plane2.push_back(hit->PeakAmplitude());
         t_Hit_StartTick_Plane2.push_back(hit->StartTick());
         t_Hit_EndTick_Plane2.push_back(hit->EndTick());
      }
   } // loop over hits

   // Fill tree! 
   t_HitTree->Fill();
}

//////////////////////////////////////////////////////////////////

void hyperon::HitTreeMaker::beginJob(){

   if(fDebug) std::cout << "Begin job" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //              Wire Tree		   //
   //////////////////////////////////////////

   t_HitTree=tfs->make<TTree>("HitTree","Wire Tree");

   t_HitTree->Branch("EventID",&t_EventID);
   t_HitTree->Branch("run",&t_run);
   t_HitTree->Branch("subrun",&t_subrun);
   t_HitTree->Branch("event",&t_event);

   // All the hits in the event
   t_HitTree->Branch("Hit_Channel_Plane0",&t_Hit_Channel_Plane0);
   t_HitTree->Branch("Hit_Tick_Plane0",&t_Hit_Tick_Plane0);
   t_HitTree->Branch("Hit_Amplitude_Plane0",&t_Hit_Amplitude_Plane0);
   t_HitTree->Branch("Hit_StartTick_Plane0",&t_Hit_StartTick_Plane0);
   t_HitTree->Branch("Hit_EndTick_Plane0",&t_Hit_EndTick_Plane0);

   t_HitTree->Branch("Hit_Channel_Plane1",&t_Hit_Channel_Plane1);
   t_HitTree->Branch("Hit_Tick_Plane1",&t_Hit_Tick_Plane1);
   t_HitTree->Branch("Hit_Amplitude_Plane1",&t_Hit_Amplitude_Plane1);
   t_HitTree->Branch("Hit_StartTick_Plane1",&t_Hit_StartTick_Plane1);
   t_HitTree->Branch("Hit_EndTick_Plane1",&t_Hit_EndTick_Plane1);

   t_HitTree->Branch("Hit_Channel_Plane2",&t_Hit_Channel_Plane2);
   t_HitTree->Branch("Hit_Tick_Plane2",&t_Hit_Tick_Plane2);
   t_HitTree->Branch("Hit_Amplitude_Plane2",&t_Hit_Amplitude_Plane2);
   t_HitTree->Branch("Hit_StartTick_Plane2",&t_Hit_StartTick_Plane2);
   t_HitTree->Branch("Hit_EndTick_Plane2",&t_Hit_EndTick_Plane2);

   // Hits associated with tracks in the neutrino slice
   t_HitTree->Branch("TrackHit_TrackIndex_Plane0",&t_TrackHit_TrackIndex_Plane0);
   t_HitTree->Branch("TrackHit_Channel_Plane0",&t_TrackHit_Channel_Plane0);
   t_HitTree->Branch("TrackHit_Tick_Plane0",&t_TrackHit_Tick_Plane0);
   t_HitTree->Branch("TrackHit_Amplitude_Plane0",&t_TrackHit_Amplitude_Plane0);
   t_HitTree->Branch("TrackHit_StartTick_Plane0",&t_TrackHit_StartTick_Plane0);
   t_HitTree->Branch("TrackHit_EndTick_Plane0",&t_TrackHit_EndTick_Plane0);

   t_HitTree->Branch("TrackHit_TrackIndex_Plane1",&t_TrackHit_TrackIndex_Plane1);
   t_HitTree->Branch("TrackHit_Channel_Plane1",&t_TrackHit_Channel_Plane1);
   t_HitTree->Branch("TrackHit_Tick_Plane1",&t_TrackHit_Tick_Plane1);
   t_HitTree->Branch("TrackHit_Amplitude_Plane1",&t_TrackHit_Amplitude_Plane1);
   t_HitTree->Branch("TrackHit_StartTick_Plane1",&t_TrackHit_StartTick_Plane1);
   t_HitTree->Branch("TrackHit_EndTick_Plane1",&t_TrackHit_EndTick_Plane1);

   t_HitTree->Branch("TrackHit_TrackIndex_Plane2",&t_TrackHit_TrackIndex_Plane2);
   t_HitTree->Branch("TrackHit_Channel_Plane2",&t_TrackHit_Channel_Plane2);
   t_HitTree->Branch("TrackHit_Tick_Plane2",&t_TrackHit_Tick_Plane2);
   t_HitTree->Branch("TrackHit_Amplitude_Plane2",&t_TrackHit_Amplitude_Plane2);
   t_HitTree->Branch("TrackHit_StartTick_Plane2",&t_TrackHit_StartTick_Plane2);
   t_HitTree->Branch("TrackHit_EndTick_Plane2",&t_TrackHit_EndTick_Plane2);

   // Hits associated with tracks in the neutrino slice
   t_HitTree->Branch("ShowerHit_ShowerIndex_Plane0",&t_ShowerHit_ShowerIndex_Plane0);
   t_HitTree->Branch("ShowerHit_Channel_Plane0",&t_ShowerHit_Channel_Plane0);
   t_HitTree->Branch("ShowerHit_Tick_Plane0",&t_ShowerHit_Tick_Plane0);
   t_HitTree->Branch("ShowerHit_Amplitude_Plane0",&t_ShowerHit_Amplitude_Plane0);
   t_HitTree->Branch("ShowerHit_StartTick_Plane0",&t_ShowerHit_StartTick_Plane0);
   t_HitTree->Branch("ShowerHit_EndTick_Plane0",&t_ShowerHit_EndTick_Plane0);

   t_HitTree->Branch("ShowerHit_ShowerIndex_Plane1",&t_ShowerHit_ShowerIndex_Plane1);
   t_HitTree->Branch("ShowerHit_Channel_Plane1",&t_ShowerHit_Channel_Plane1);
   t_HitTree->Branch("ShowerHit_Tick_Plane1",&t_ShowerHit_Tick_Plane1);
   t_HitTree->Branch("ShowerHit_Amplitude_Plane1",&t_ShowerHit_Amplitude_Plane1);
   t_HitTree->Branch("ShowerHit_StartTick_Plane1",&t_ShowerHit_StartTick_Plane1);
   t_HitTree->Branch("ShowerHit_EndTick_Plane1",&t_ShowerHit_EndTick_Plane1);

   t_HitTree->Branch("ShowerHit_ShowerIndex_Plane2",&t_ShowerHit_ShowerIndex_Plane2);
   t_HitTree->Branch("ShowerHit_Channel_Plane2",&t_ShowerHit_Channel_Plane2);
   t_HitTree->Branch("ShowerHit_Tick_Plane2",&t_ShowerHit_Tick_Plane2);
   t_HitTree->Branch("ShowerHit_Amplitude_Plane2",&t_ShowerHit_Amplitude_Plane2);
   t_HitTree->Branch("ShowerHit_StartTick_Plane2",&t_ShowerHit_StartTick_Plane2);
   t_HitTree->Branch("ShowerHit_EndTick_Plane2",&t_ShowerHit_EndTick_Plane2);

   // All hits in the neutrino slice
   t_HitTree->Branch("SliceHit_Channel_Plane0",&t_SliceHit_Channel_Plane0);
   t_HitTree->Branch("SliceHit_Tick_Plane0",&t_SliceHit_Tick_Plane0);
   t_HitTree->Branch("SliceHit_Amplitude_Plane0",&t_SliceHit_Amplitude_Plane0);
   t_HitTree->Branch("SliceHit_StartTick_Plane0",&t_SliceHit_StartTick_Plane0);
   t_HitTree->Branch("SliceHit_EndTick_Plane0",&t_SliceHit_EndTick_Plane0);

   t_HitTree->Branch("SliceHit_Channel_Plane1",&t_SliceHit_Channel_Plane1);
   t_HitTree->Branch("SliceHit_Tick_Plane1",&t_SliceHit_Tick_Plane1);
   t_HitTree->Branch("SliceHit_Amplitude_Plane1",&t_SliceHit_Amplitude_Plane1);
   t_HitTree->Branch("SliceHit_StartTick_Plane1",&t_SliceHit_StartTick_Plane1);
   t_HitTree->Branch("SliceHit_EndTick_Plane1",&t_SliceHit_EndTick_Plane1);

   t_HitTree->Branch("SliceHit_Channel_Plane2",&t_SliceHit_Channel_Plane2);
   t_HitTree->Branch("SliceHit_Tick_Plane2",&t_SliceHit_Tick_Plane2);
   t_HitTree->Branch("SliceHit_Amplitude_Plane2",&t_SliceHit_Amplitude_Plane2);
   t_HitTree->Branch("SliceHit_StartTick_Plane2",&t_SliceHit_StartTick_Plane2);
   t_HitTree->Branch("SliceHit_EndTick_Plane2",&t_SliceHit_EndTick_Plane2);

   // Hits in neutrino slice not associated with any track or shower
   t_HitTree->Branch("SpareSliceHit_Channel_Plane0",&t_SpareSliceHit_Channel_Plane0);
   t_HitTree->Branch("SpareSliceHit_Tick_Plane0",&t_SpareSliceHit_Tick_Plane0);
   t_HitTree->Branch("SpareSliceHit_Amplitude_Plane0",&t_SpareSliceHit_Amplitude_Plane0);
   t_HitTree->Branch("SpareSliceHit_StartTick_Plane0",&t_SpareSliceHit_StartTick_Plane0);
   t_HitTree->Branch("SpareSliceHit_EndTick_Plane0",&t_SpareSliceHit_EndTick_Plane0);

   t_HitTree->Branch("SpareSliceHit_Channel_Plane1",&t_SpareSliceHit_Channel_Plane1);
   t_HitTree->Branch("SpareSliceHit_Tick_Plane1",&t_SpareSliceHit_Tick_Plane1);
   t_HitTree->Branch("SpareSliceHit_Amplitude_Plane1",&t_SpareSliceHit_Amplitude_Plane1);
   t_HitTree->Branch("SpareSliceHit_StartTick_Plane1",&t_SpareSliceHit_StartTick_Plane1);
   t_HitTree->Branch("SpareSliceHit_EndTick_Plane1",&t_SpareSliceHit_EndTick_Plane1);

   t_HitTree->Branch("SpareSliceHit_Channel_Plane2",&t_SpareSliceHit_Channel_Plane2);
   t_HitTree->Branch("SpareSliceHit_Tick_Plane2",&t_SpareSliceHit_Tick_Plane2);
   t_HitTree->Branch("SpareSliceHit_Amplitude_Plane2",&t_SpareSliceHit_Amplitude_Plane2);
   t_HitTree->Branch("SpareSliceHit_StartTick_Plane2",&t_SpareSliceHit_StartTick_Plane2);
   t_HitTree->Branch("SpareSliceHit_EndTick_Plane2",&t_SpareSliceHit_EndTick_Plane2);

   // Track starts/ends
   t_HitTree->Branch("TrackStart_Channel_Plane0",&t_TrackStart_Channel_Plane0);
   t_HitTree->Branch("TrackStart_Time_Plane0",&t_TrackStart_Time_Plane0);
   t_HitTree->Branch("TrackEnd_Channel_Plane0",&t_TrackEnd_Channel_Plane0);
   t_HitTree->Branch("TrackEnd_Time_Plane0",&t_TrackEnd_Time_Plane0);

   t_HitTree->Branch("TrackStart_Channel_Plane1",&t_TrackStart_Channel_Plane1);
   t_HitTree->Branch("TrackStart_Time_Plane1",&t_TrackStart_Time_Plane1);
   t_HitTree->Branch("TrackEnd_Channel_Plane1",&t_TrackEnd_Channel_Plane1);
   t_HitTree->Branch("TrackEnd_Time_Plane1",&t_TrackEnd_Time_Plane1);

   t_HitTree->Branch("TrackStart_Channel_Plane2",&t_TrackStart_Channel_Plane2);
   t_HitTree->Branch("TrackStart_Time_Plane2",&t_TrackStart_Time_Plane2);
   t_HitTree->Branch("TrackEnd_Channel_Plane2",&t_TrackEnd_Channel_Plane2);
   t_HitTree->Branch("TrackEnd_Time_Plane2",&t_TrackEnd_Time_Plane2);

   t_HitTree->Branch("TrackStart_X",&t_TrackStart_X);
   t_HitTree->Branch("TrackStart_Y",&t_TrackStart_Y);
   t_HitTree->Branch("TrackStart_Z",&t_TrackStart_Z);
   t_HitTree->Branch("TrackEnd_X",&t_TrackEnd_X);
   t_HitTree->Branch("TrackEnd_Y",&t_TrackEnd_Y);
   t_HitTree->Branch("TrackEnd_Z",&t_TrackEnd_Z);

   t_HitTree->Branch("TrackDir_X",&t_TrackDir_X);
   t_HitTree->Branch("TrackDir_Y",&t_TrackDir_Y);
   t_HitTree->Branch("TrackDir_Z",&t_TrackDir_Z);
   t_HitTree->Branch("TrackEndDir_X",&t_TrackEndDir_X);
   t_HitTree->Branch("TrackEndDir_Y",&t_TrackEndDir_Y);
   t_HitTree->Branch("TrackEndDir_Z",&t_TrackEndDir_Z);

   t_HitTree->Branch("TrackGrad_Plane0",&t_TrackGrad_Plane0);
   t_HitTree->Branch("TrackGrad_Plane1",&t_TrackGrad_Plane1);
   t_HitTree->Branch("TrackGrad_Plane2",&t_TrackGrad_Plane2);
   t_HitTree->Branch("TrackEndGrad_Plane0",&t_TrackEndGrad_Plane0);
   t_HitTree->Branch("TrackEndGrad_Plane1",&t_TrackEndGrad_Plane1);
   t_HitTree->Branch("TrackEndGrad_Plane2",&t_TrackEndGrad_Plane2);
   t_HitTree->Branch("TrackAngle_Plane0",&t_TrackAngle_Plane0);
   t_HitTree->Branch("TrackAngle_Plane1",&t_TrackAngle_Plane1);
   t_HitTree->Branch("TrackAngle_Plane2",&t_TrackAngle_Plane2);
   t_HitTree->Branch("TrackEndAngle_Plane0",&t_TrackEndAngle_Plane0);
   t_HitTree->Branch("TrackEndAngle_Plane1",&t_TrackEndAngle_Plane1);
   t_HitTree->Branch("TrackEndAngle_Plane2",&t_TrackEndAngle_Plane2);

   t_HitTree->Branch("CaloStart_Channel_Plane0",&t_CaloStart_Channel_Plane0);
   t_HitTree->Branch("CaloStart_Time_Plane0",&t_CaloStart_Time_Plane0);        
   t_HitTree->Branch("CaloStart_X_Plane0",&t_CaloStart_X_Plane0);
   t_HitTree->Branch("CaloStart_Y_Plane0",&t_CaloStart_Y_Plane0);
   t_HitTree->Branch("CaloStart_Z_Plane0",&t_CaloStart_Z_Plane0);

   t_HitTree->Branch("CaloStart_Channel_Plane1",&t_CaloStart_Channel_Plane1);
   t_HitTree->Branch("CaloStart_Time_Plane1",&t_CaloStart_Time_Plane1);
   t_HitTree->Branch("CaloStart_X_Plane1",&t_CaloStart_X_Plane1);
   t_HitTree->Branch("CaloStart_Y_Plane1",&t_CaloStart_Y_Plane1);
   t_HitTree->Branch("CaloStart_Z_Plane1",&t_CaloStart_Z_Plane1);

   t_HitTree->Branch("CaloStart_Channel_Plane2",&t_CaloStart_Channel_Plane2);
   t_HitTree->Branch("CaloStart_Time_Plane2",&t_CaloStart_Time_Plane2);
   t_HitTree->Branch("CaloStart_X_Plane2",&t_CaloStart_X_Plane2);
   t_HitTree->Branch("CaloStart_Y_Plane2",&t_CaloStart_Y_Plane2);
   t_HitTree->Branch("CaloStart_Z_Plane2",&t_CaloStart_Z_Plane2);

   // Shower starts/ends
   t_HitTree->Branch("ShowerStart_Channel_Plane0",&t_ShowerStart_Channel_Plane0);
   t_HitTree->Branch("ShowerStart_Time_Plane0",&t_ShowerStart_Time_Plane0);

   t_HitTree->Branch("ShowerStart_Channel_Plane1",&t_ShowerStart_Channel_Plane1);
   t_HitTree->Branch("ShowerStart_Time_Plane1",&t_ShowerStart_Time_Plane1);

   t_HitTree->Branch("ShowerStart_Channel_Plane2",&t_ShowerStart_Channel_Plane2);
   t_HitTree->Branch("ShowerStart_Time_Plane2",&t_ShowerStart_Time_Plane2);

   t_HitTree->Branch("ShowerStart_X",&t_ShowerStart_X);
   t_HitTree->Branch("ShowerStart_Y",&t_ShowerStart_Y);
   t_HitTree->Branch("ShowerStart_Z",&t_ShowerStart_Z);

   t_HitTree->Branch("ShowerDir_X",&t_ShowerDir_X);
   t_HitTree->Branch("ShowerDir_Y",&t_ShowerDir_Y);
   t_HitTree->Branch("ShowerDir_Z",&t_ShowerDir_Z);

   if(fDebug) std::cout << "Finished begin job" << std::endl;
}

void hyperon::HitTreeMaker::endJob()
{

}

void hyperon::HitTreeMaker::beginSubRun(const art::SubRun& sr)
{

}

void hyperon::HitTreeMaker::endSubRun(const art::SubRun& sr)
{

}

DEFINE_ART_MODULE(hyperon::HitTreeMaker)
