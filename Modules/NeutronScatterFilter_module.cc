////////////////////////////////////////////////////////////////////////
// Class:       NeutronScatterFilter
// Plugin Type: filter (art v3_01_02)
// File:        NeutronScatterFilter_module.cc
//
// Purpose: Finds events containing neutron scatters
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
   class NeutronScatterFilter;
}


class hyperon::NeutronScatterFilter : public art::EDFilter {
   public:
      explicit NeutronScatterFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      NeutronScatterFilter(NeutronScatterFilter const&) = delete;
      NeutronScatterFilter(NeutronScatterFilter&&) = delete;
      NeutronScatterFilter& operator=(NeutronScatterFilter const&) = delete;
      NeutronScatterFilter& operator=(NeutronScatterFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:


      bool f_GetGeneratorInfo;
      bool f_GetG4Info;

      fhicl::ParameterSet f_Generator;
      fhicl::ParameterSet f_G4;

      bool f_Debug = false;
};


hyperon::NeutronScatterFilter::NeutronScatterFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   f_GetGeneratorInfo(p.get<bool>("GetGeneratorInfo",true)),   
   f_GetG4Info(p.get<bool>("GetG4Info",true)),   
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Debug(p.get<bool>("Debug",false))
{
}

bool hyperon::NeutronScatterFilter::filter(art::Event& e)
{
   bool pass = true;

   if(f_GetGeneratorInfo){
      if(f_Debug) std::cout << "Getting EG Info" << std::endl;
      SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e,f_Generator);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();
      pass = GenT.EventHasFinalStateNeutron;
      delete Generator_SM;
   }

   if(f_GetG4Info){
      if(f_Debug) std::cout << "Getting G4 Info" << std::endl;
      SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4);
      G4Truth G4T = G4_SM->GetG4Info();
      pass = G4T.EventHasNeutronScatter;       
      delete G4_SM;
   }

   return pass;
}

DEFINE_ART_MODULE(hyperon::NeutronScatterFilter)
