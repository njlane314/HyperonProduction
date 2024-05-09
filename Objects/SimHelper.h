#ifndef _SimHelper_h_
#define _SimHelper_h_

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

namespace hyperon {

   inline SimParticle MakeSimParticle(simb::MCParticle Part)
   {
      SimParticle S;
      S.SetKinematics(Part.Momentum(), Part.EndMomentum(), Part.Mass());
      S.PDG = Part.PdgCode();
      S.SetPositions(Part.Position(), Part.EndPosition());

      return S;
   }

}

#endif