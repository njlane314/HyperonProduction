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
#include "ubana/HyperonProduction/Headers/LLRPID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_ProtonMuonLookup.h"
#include "ubana/HyperonProduction/Headers/LLRPID_CorrectionLookup.h"
#include "TVector3.h"

namespace hyperon {

   // Define wire plane angles

   const double wireAnglePlane0 = 30 * (6.28 / 360.0);
   const double wireAnglePlane1 = -30 * (6.28 / 360.0);
   const double wireAnglePlane2 = 90 * (6.28 / 360.0);

   struct ParticleIdentifierStore 
   {
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

      double BraggWeighted_Pion; 
      double BraggWeighted_Muon;
      double BraggWeighted_Proton;
      double BraggWeighted_Kaon;
      double BraggWeighted_Sigma; 

      std::vector<float> dEdX_Plane0;
      std::vector<float> dEdX_Plane1;
      std::vector<float> dEdX_Plane2;

      std::vector<float> ResidualRange_Plane0;
      std::vector<float> ResidualRange_Plane1;
      std::vector<float> ResidualRange_Plane2;

      std::vector<float> Pitch_Plane0;
      std::vector<float> Pitch_Plane1;
      std::vector<float> Pitch_Plane2;

      std::vector<float> dEdX_Corrected_Plane0;
      std::vector<float> dEdX_Corrected_Plane1;
      std::vector<float> dEdX_Corrected_Plane2;
   };

   class PIDManager {

      public: 

         PIDManager();

         double GetMeandEdX(art::Ptr<anab::Calorimetry> Calo);
         void SetThreePlaneMeandEdX(art::Ptr<recob::Track> Trk, std::vector<art::Ptr<anab::Calorimetry>> VectCalo, ParticleIdentifierStore &Store);
         void LLRPID(std::vector<art::Ptr<anab::Calorimetry>> VectCalo, ParticleIdentifierStore& Store);
         double GetBraggLikelihood(art::Ptr<recob::Track> Trk, std::vector<anab::sParticleIDAlgScores> VectorAlgScores, int PDG, anab::kTrackDir Dir);
         void SetBraggScores(art::Ptr<recob::Track> Trk, std::vector<anab::sParticleIDAlgScores> VectorAlgScores, ParticleIdentifierStore &Store);

         ParticleIdentifierStore GetPIDScores(art::Ptr<recob::Track> Trk, std::vector<art::Ptr<anab::Calorimetry>> VectCalo, std::vector<anab::sParticleIDAlgScores> VectorAlgScores);   
               
         double GetPlaneWeight(art::Ptr<recob::Track> Trk, int Plane);

         void SetTophatThresh(double thresh){ tophatThresh = thresh; }

      private:

         searchingfornues::LLRPID particleIdentificationCalculator;
         searchingfornues::ProtonMuonLookUpParameters parametersProtonMuon;
         searchingfornues::CorrectionLookUpParameters parametersCorrections;

         // Miniumum value of sin2(angle between track and wires)
         double tophatThresh = 0.175;
         double resRangeCutoff = 5; 
   };
}

#endif
