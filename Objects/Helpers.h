#ifndef _Helpers_h_
#define _Helpers_h_

//local includes
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Headers/TrackWiggliness.h"

//larsoft objects
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

//uboonecode objects
#include "nusimdata/SimulationBase/MCParticle.h"

//root objects
#include "TVector3.h"

// Helper functions used to transform larsoft objects into SimParticle and RecoParticle

namespace hyperon {

inline SimParticle MakeSimParticle(simb::MCParticle Part){

   SimParticle S;
   S.SetKinematics(Part.Momentum(),Part.EndMomentum(),Part.Mass());
   S.PDG = Part.PdgCode();
   S.SetPositions(Part.Position(),Part.EndPosition());
   return S;
}

// Helper function for setting track variables in Reco Particle
inline void SetTrackVariables(RecoParticle &P , art::Ptr<recob::Track> trk){

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   P.TrackLength = trk->Length();
   P.TrackDirectionX = trk->StartDirection().X();
   P.TrackDirectionY = trk->StartDirection().Y();
   P.TrackDirectionZ = trk->StartDirection().Z();

   trkf::TrackMomentumCalculator trkm{0};
   P.ProtonMomentum = trkm.GetTrackMomentum(trk->Length(),2212);	
   P.MuonMomentum = trkm.GetTrackMomentum(trk->Length(),13);
   // Kaon momentum estimator - see docdb # 38619
   P.KaonMomentum = (156.222*pow(P.TrackLength-0.0426198,0.274777) + 1.5323*P.TrackLength)/1e3;

   geo::Point_t point = {trk->Start().X(),trk->Start().Y(),trk->Start().Z()};
   geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
   TVector3 Start(trk->Start().X()+sce_corr.X(),trk->Start().Y()-sce_corr.Y(),trk->Start().Z()-sce_corr.Z());
   point = {trk->End().X(),trk->End().Y(),trk->End().Z()};
   sce_corr = SCE->GetPosOffsets(point);
   TVector3 End(trk->End().X()+sce_corr.X(),trk->End().Y()-sce_corr.Y(),trk->End().Z()-sce_corr.Z());

   P.SetTrackPositions(Start,End);

   P.TrackWiggliness = GetTrackWiggliness(trk);
}

}

#endif
