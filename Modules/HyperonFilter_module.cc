////////////////////////////////////////////////////////////////////////
// Class:       HyperonFilter
// Plugin Type: filter (art v3_01_02)
// File:        HyperonFilter_module.cc
//
// Purpose: Finds events containing hyperon production (direct and associated)
//
// Generated at Thu Sep 16 09:04:20 2021 by Christopher Thorpe using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"

namespace hyperon {
   class HyperonFilter;
}


class hyperon::HyperonFilter : public art::EDFilter {
   public:
      explicit HyperonFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      HyperonFilter(HyperonFilter const&) = delete;
      HyperonFilter(HyperonFilter&&) = delete;
      HyperonFilter& operator=(HyperonFilter const&) = delete;
      HyperonFilter& operator=(HyperonFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      fhicl::ParameterSet f_G4;
};


hyperon::HyperonFilter::HyperonFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   f_G4(p.get<fhicl::ParameterSet>("Geant4"))
{
}

bool hyperon::HyperonFilter::filter(art::Event& e)
{

      SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4);
      G4Truth G4T = G4_SM->GetG4Info();

      bool EventHasHyperon = G4T.EventHasHyperon;       
      
      delete G4_SM;
  
      return EventHasHyperon;
}

DEFINE_ART_MODULE(hyperon::HyperonFilter)
