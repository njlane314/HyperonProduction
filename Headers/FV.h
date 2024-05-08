#ifndef _FV_h_
#define _FV_h_

#include "TVector3.h"

namespace hyperon {

const std::vector<double> TPCCenter = { 126.625 , 0.97 , 518.5 }; //center of active TPC
const std::vector<double> TPCSideLengths = { 236.35 , 233.0 , 1036.8 }; //side lengths of active TPC

// Use whole TPC - decide which volumes to cut later
const double FVxmin = 0.0;
const double FVxmax = 256.35;
const double FVymin = -115.53;
const double FVymax = 117.47;
const double FVzmin = 0.1;
const double FVzmax = 1036.9;

inline bool inActiveTPC(TVector3 pos){
   if(pos.X() > FVxmax || pos.X() < FVxmin) return false;
   if(pos.Y() > FVymax || pos.Y() < FVymin) return false;
   if(pos.Z() > FVzmax || pos.Z() < FVzmin) return false;
   return true;
}

}

#endif