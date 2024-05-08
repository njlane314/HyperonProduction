#ifndef _GeneratorAnalyser_h_
#define _GeneratorAnalyser_h_

#include <string>
#include <vector>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "cetlib_except/exception.h"

#include "ubana/HyperonProduction/Headers/FV.h"
#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"

namespace hyperon {

struct GeneratorTruth 
{
   double Weight = 1.0;

   int nMCTruths = 0;
   int nMCTruthsInTPC = 0;

   std::vector<SimParticle> Neutrinos;

   std::vector<std::string> CCNC;
   std::vector<std::string> Mode;

   std::vector<double> W;
   std::vector<double> X;
   std::vector<double> Y;
   std::vector<double> QSqr; 
   std::vector<double> Pt;
   std::vector<double> Theta;

   bool eventHasNeutralKaon;
   bool eventHasHyperon;
   bool eventHasNeutron;
};

class GeneratorAnalyser {

public:

   GeneratorAnalyser(art::Event const& e, fhicl::ParameterSet pset, bool particlegunmode = false);
   GeneratorTruth GetGeneratorTruth();

private:

   art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
   std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;

   GeneratorTruth generatorTruth;

   const bool particleGunMode;
};

}

#endif