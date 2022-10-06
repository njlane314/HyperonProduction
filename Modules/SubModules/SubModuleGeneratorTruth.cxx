#ifndef _SubModuleGeneratorTruth_cxx_
#define _SubModuleGeneratorTruth_cxx_

#include "SubModuleGeneratorTruth.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleGeneratorTruth::SubModuleGeneratorTruth(art::Event const& e,fhicl::ParameterSet pset){

   if(!e.getByLabel(pset.get<std::string>("GeneratorModuleLabel","generator"),Handle_MCTruth))  
      throw cet::exception("SubModuleGeneratorTruth") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  

   HyperonPDGs = pset.get<std::vector<int>>("HyperonPDGs",{3122,3212,3112,3222});
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

GeneratorTruth SubModuleGeneratorTruth::GetGeneratorTruth(){

   if(!Vect_MCTruth.size()){
      std::cout << "MCTruth vector is empty" << std::endl;
      return theTruth;
   }

   theTruth.NMCTruths = Vect_MCTruth.size();

   int i_truth=0;

   for(const art::Ptr<simb::MCTruth> &theMCTruth : Vect_MCTruth){

      simb::MCNeutrino Nu = theMCTruth->GetNeutrino();

      int mode = Nu.Mode();
      int ccnc = Nu.CCNC();

      if(ccnc == 0) theTruth.CCNC.push_back("CC");
      else theTruth.CCNC.push_back("NC");

      if(mode == 0) theTruth.Mode.push_back("QEL");
      else if(mode == 1) theTruth.Mode.push_back("RES");
      else if(mode == 2) theTruth.Mode.push_back("DIS");
      else if(mode == 3) theTruth.Mode.push_back("COH");
      else if(mode == 5) theTruth.Mode.push_back("ElectronScattering");
      else if(mode == 10) theTruth.Mode.push_back("MEC");
      else if(mode == 11) theTruth.Mode.push_back("Diffractive");
      else if(mode == 1095) theTruth.Mode.push_back("HYP");
      else theTruth.Mode.push_back("Other");	

      for(int k_particles=0;k_particles<theMCTruth->NParticles();k_particles++){

         simb::MCParticle Part = theMCTruth->GetParticle(k_particles);

         //if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) 
         // theTruth.TruePrimaryVertex.SetXYZ(Part.Vx(),Part.Vy(),Part.Vz());

         if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) {
            theTruth.TruePrimaryVertex_X.push_back(Part.Vx());
            theTruth.TruePrimaryVertex_Y.push_back(Part.Vy());
            theTruth.TruePrimaryVertex_Z.push_back(Part.Vz());
            if(inActiveTPC(TVector3(Part.Vx(),Part.Vy(),Part.Vz()))) theTruth.NMCTruthsInTPC++;
         }

         if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){
            SimParticle P = MakeSimParticle(Part);
            P.Origin = 0;
            P.MCTruthIndex = i_truth;
            theTruth.Neutrino.push_back(P);
         }

         // If there is a hyperon in the final state in a QEL event, change mode to HYP
         if(isHyperon(Part.PdgCode()) && Part.StatusCode() == 1 && mode == 0) theTruth.Mode.back() = "HYP";
          
         if(Part.StatusCode() == 1 && Part.PdgCode() == 2112) theTruth.EventHasFinalStateNeutron = true;
         if(Part.StatusCode() == 1 && isHyperon(Part.PdgCode()) && std::find(HyperonPDGs.begin(),HyperonPDGs.end(),abs(Part.PdgCode())) != HyperonPDGs.end()){
            theTruth.EventHasHyperon = true;
         }
      }

      i_truth++;
   }

      if(theTruth.Neutrino.size() != Vect_MCTruth.size())         
         throw cet::exception("SubModuleGeneratorTruth") << "Sim Neutrino/MCTruth vector size mismatch" << std::endl;
    
   return theTruth;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
