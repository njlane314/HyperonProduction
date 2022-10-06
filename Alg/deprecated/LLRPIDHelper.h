#ifndef _LLRPIDHelper_h_
#define _LLRPIDHelper_h_

#include <string>
#include <vector>

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "cetlib_except/exception.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "ubana/HyperonProduction/Headers/LLR_PID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_proton_muon_lookup.h"
#include "ubana/HyperonProduction/Headers/LLR_PID_K.h"
#include "ubana/HyperonProduction/Headers/LLRPID_kaon_proton_lookup.h"

struct LLRPID_Result {
   double Score;
   double Score_Kaon;
   double Score_Kaon_Partial;
};

class LLRPIDHelper {

   public:

      LLRPIDHelper();
      LLRPID_Result GetScores(std::vector<art::Ptr<anab::Calorimetry>> calo_v);

   private:

      searchingfornues::LLRPID llr_pid_calculator;
      searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
      searchingfornuesk::LLRPIDK llr_pid_calculator_kaon;
      searchingfornuesk::KaonProtonLookUpParameters kaonproton_parameters;

      double ResRangeCutoff=5; 
};

#endif
