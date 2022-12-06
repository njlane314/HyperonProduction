#ifndef _BDTHandle_h_
#define _BDTHandle_h_

#include <string>
#include <iostream>
#include <cstdlib>

#include "cetlib_except/exception.h"

#include "TString.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

namespace hyperon {

class BDTHandle {

   public:

      BDTHandle(std::string weightsdir);
      double GetScore(double kaontrackpID,double kaontrackbraggpid);

   private:

      TMVA::Reader *reader=nullptr;

      Float_t v_KaonTrackPID;
      Float_t v_KaonTrackBraggPID;

      const std::string WeightsDir;
};

}

#endif
