#ifndef _TrackWiggliness_h_
#define _TrackWiggliness_h_

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "lardataobj/RecoBase/Track.h"

// Functions for calculating track wiggliness, created by Andy Smith (U of Cambridge)
// https://github.com/a-d-smith/ubcc1pi/blob/3a0636df2e2d11581608857e862c17b5bcbfbc10/ubcc1pi/Helpers/RecoHelper.cxx

namespace hyperon {

inline std::vector<size_t> GetValidPoints(const art::Ptr<recob::Track> &track){

   std::vector<size_t> validPoints;

   const auto firstValidPoint = track->FirstValidPoint();
   validPoints.push_back(firstValidPoint);

   auto nextValidPoint = track->NextValidPoint(firstValidPoint + 1);
   while(nextValidPoint != recob::TrackTrajectory::InvalidIndex){
      validPoints.push_back(nextValidPoint);
      nextValidPoint = track->NextValidPoint(nextValidPoint + 1);
   }

   return validPoints;
}

inline float GetTrackWiggliness(const art::Ptr<recob::Track> &track){

   const auto validPoints = GetValidPoints(track);
   if (validPoints.size() < 3)
      return 0.f;

   std::vector<float> thetaVector;
   float thetaSum = 0.f;
   for (unsigned int i = 1; i < validPoints.size(); ++i){
      const auto dir = track->DirectionAtPoint(validPoints.at(i));
      const auto dirPrev = track->DirectionAtPoint(validPoints.at(i - 1));

      // Bind between -1 and 1 at floating precision to avoid issues with cast from double
      const auto cosTheta = std::min(1.f, std::max(-1.f, static_cast<float>(dir.Dot(dirPrev))));
      const auto theta = std::acos(cosTheta);

      thetaSum += theta;
      thetaVector.push_back(theta);
   }

   const auto thetaMean = thetaSum / static_cast<float>(thetaVector.size());

   float thetaDiffSum = 0.f;
   for (const auto &theta : thetaVector)
      thetaDiffSum += std::pow(theta - thetaMean, 2);
   
   const auto variance = thetaDiffSum / static_cast<float>(thetaVector.size() - 1);
   return std::sqrt(variance);
}

}

#endif
