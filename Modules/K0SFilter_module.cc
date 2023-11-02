////////////////////////////////////////////////////////////////////////
// Class:       K0SFilter
// Plugin Type: filter (art v3_01_02)
// File:        K0SFilter_module.cc
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

#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"

namespace hyperon {
   class K0SFilter;
}


class hyperon::K0SFilter : public art::EDFilter {
   public:
      explicit K0SFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      K0SFilter(K0SFilter const&) = delete;
      K0SFilter(K0SFilter&&) = delete;
      K0SFilter& operator=(K0SFilter const&) = delete;
      K0SFilter& operator=(K0SFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      bool f_GetGeneratorInfo;
      bool f_GetG4Info;

      fhicl::ParameterSet f_Generator;
      fhicl::ParameterSet f_G4;

      bool f_Debug = false;
};


hyperon::K0SFilter::K0SFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   f_GetGeneratorInfo(p.get<bool>("GetGeneratorInfo",true)),   
   f_GetG4Info(p.get<bool>("GetG4Info",true)),   
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Debug(p.get<bool>("Debug",false))
{
}

bool hyperon::K0SFilter::filter(art::Event& e)
{

   bool pass = false;

   if(f_GetGeneratorInfo){
     if(f_Debug) std::cout << "Getting EG Info" << std::endl;
     SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e,f_Generator);
     GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();
     pass = GenT.EventHasKaon;
   }

   if(f_GetG4Info){
     if(f_Debug) std::cout << "Getting G4 Info" << std::endl;
     SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4);
     G4Truth G4T = G4_SM->GetG4Info();
     pass = G4T.EventHasK0S;       
     delete G4_SM;
   }

return pass;

}

DEFINE_ART_MODULE(hyperon::K0SFilter)
