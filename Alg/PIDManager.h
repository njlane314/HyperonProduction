#ifndef _PIDManager_h_
#define _PIDManager_h_

#include <string>
#include <vector>
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "cetlib_except/exception.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "ubana/HyperonProduction/Headers/LLR_PID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_proton_muon_lookup.h"
#include "ubana/HyperonProduction/Headers/LLR_PID_K.h"
#include "ubana/HyperonProduction/Headers/LLRPID_kaon_proton_lookup.h"
#include "ubana/HyperonProduction/Headers/LLRPID_correction_lookup.h"
#include "TVector3.h"

namespace hyperon {

   //Define wire plane angles
   const double plane0_wireangle = 30*6.28/360.0;
   const double plane1_wireangle = -30*6.28/360.0;
   const double plane2_wireangle = 90*6.28/360.0;

   struct PIDStore {

      double Weight_Plane0 = 0;
      double Weight_Plane1 = 0;
      double Weight_Plane2 = 0;

      double MeandEdX_Plane0 = -1;
      double MeandEdX_Plane1 = -1;
      double MeandEdX_Plane2 = -1;
      double MeandEdX_3Plane = -1;

      double LLR;
      double LLR_Kaon;
      double LLR_Kaon_Partial; 

      double Bragg_Kaon_Plane0;
      double Bragg_Kaon_Plane1;
      double Bragg_Kaon_Plane2;    
      double Bragg_Kaon_3Plane;

      std::vector<float> dEdX_Plane0;
      std::vector<float> ResidualRange_Plane0;
      std::vector<float> Pitch_Plane0;
      std::vector<float> dEdX_Plane1;
      std::vector<float> ResidualRange_Plane1;
      std::vector<float> Pitch_Plane1;
      std::vector<float> dEdX_Plane2;
      std::vector<float> ResidualRange_Plane2;
      std::vector<float> Pitch_Plane2;

      std::vector<float> dEdX_Corrected_Plane0;
      std::vector<float> dEdX_Corrected_Plane1;
      std::vector<float> dEdX_Corrected_Plane2;

   };

   class PIDManager {

      public: 

         PIDManager();

         double GetMeandEdX(art::Ptr<anab::Calorimetry> calo);
         void ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store);
         void LLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store);
         void BraggPID(art::Ptr<recob::Track> track,std::vector<anab::sParticleIDAlgScores> algscores_v,PIDStore& store);
         PIDStore GetPIDs(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::vector<anab::sParticleIDAlgScores> algscores_v);   
               
         double PlaneWeight(TVector3 dir,int i_pl);
         double PlaneWeight(art::Ptr<recob::Track> track,int i_pl);

         void SetTophatThresh(double thresh){ TophatThresh = thresh; }

      private:

         searchingfornues::LLRPID llr_pid_calculator;
         searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
         searchingfornuesk::LLRPIDK llr_pid_calculator_kaon;
         searchingfornuesk::KaonProtonLookUpParameters kaonproton_parameters;
         searchingfornues::CorrectionLookUpParameters correction_parameters;
         

         // Miniumum value of sin2(angle between track and wires)
         double TophatThresh = 0.175;
         double ResRangeCutoff=5; 
   };
}

#endif
