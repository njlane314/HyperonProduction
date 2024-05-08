#ifndef _GeneratorAnalyser_cxx_
#define _GeneratorAnalyser_cxx_

#include "GeneratorAnalyser.h"

using namespace hyperon;

GeneratorAnalyser::GeneratorAnalyser(art::Event const& e, fhicl::ParameterSet pset, bool particlegunmode) :
particleGunMode(particlegunmode)
{
   if(!e.getByLabel(pset.get<std::string>("GeneratorModuleLabel", "generator"), Handle_MCTruth))  
      throw cet::exception("GeneratorAnalyser") << "No MCTruth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth, Handle_MCTruth);  
}

GeneratorTruth GeneratorAnalyser::GetGeneratorTruth()
{
   if(!Vect_MCTruth.size()){
      std::cout << "MCTruth vector is empty" << std::endl;
      return generatorTruth;
   }

   generatorTruth.nMCTruths = Vect_MCTruth.size();

   int TruthIndex = 0;

   bool EventHasNeutralKaon = false;
   bool EventHasHyperon = false;
   bool EventHasNeutron = false;

   for(const art::Ptr<simb::MCTruth> &MCTruth : Vect_MCTruth){

      simb::MCNeutrino Nu = MCTruth->GetNeutrino();
      int NParticles = MCTruth->NParticles();

      int CCNC = Nu.CCNC();
      int Mode = Nu.Mode();
      int InteractionType = Nu.InteractionType();

      double W = Nu.W();
      double X = Nu.X();
      double Y = Nu.Y();
      double QSqr = Nu.QSqr();
      double Pt = Nu.Pt();
      double Theta = Nu.Theta();

      std::cout << "Neutrinos interaction mode: " << Mode << std::endl;
      std::cout << "Neutrinos interaction type: " << InteractionType << std::endl;

      if(CCNC == 0){
         generatorTruth.CCNC.push_back("CC");
      }
      else{
         generatorTruth.CCNC.push_back("NC");
      }

      generatorTruth.W.push_back(W);
      generatorTruth.X.push_back(X);
      generatorTruth.Y.push_back(Y);
      generatorTruth.QSqr.push_back(QSqr);
      generatorTruth.Pt.push_back(Pt);
      generatorTruth.Theta.push_back(Theta);

      if(Mode == 0){
         generatorTruth.Mode.push_back("QEL");
      }
      else if(Mode == 1){
         generatorTruth.Mode.push_back("RES");
      }
      else if(Mode == 2){
         generatorTruth.Mode.push_back("DIS");
      }
      else if(Mode == 3){
         generatorTruth.Mode.push_back("COH");
      }
      else if(Mode == 10){
         generatorTruth.Mode.push_back("MEC");
      }
      else if(Mode == 1095){
         generatorTruth.Mode.push_back("HYP");
      }
      else{
         generatorTruth.Mode.push_back("Other");	
      }

      for(int i = 0; i < NParticles; i++){

         simb::MCParticle Particle = MCTruth->GetParticle(i);

         if(Particle.StatusCode() == 0){
            if(isNeutrino(Particle.PdgCode())){
               SimParticle P = MakeSimParticle(Particle);
               P.Origin = 0;
               P.MCTruthIndex = TruthIndex; 

               generatorTruth.Neutrinos.push_back(P);
            }
         }
         else if(Particle.StatusCode() == 1){
            if(isHyperon(Particle.PdgCode())){
               EventHasHyperon = true;
            }
            else if(isNeutralKaon(Particle.PdgCode())){
               EventHasNeutralKaon = true;
            }
            else if(isNeutron(Particle.PdgCode())){
               EventHasNeutron = true;
            }
         }
      }

      TruthIndex++;
   }

   generatorTruth.eventHasNeutralKaon = EventHasNeutralKaon;
   generatorTruth.eventHasHyperon = EventHasHyperon;
   generatorTruth.eventHasNeutron = EventHasNeutron;

   if(!particleGunMode && generatorTruth.Neutrinos.size() != Vect_MCTruth.size())         
      throw cet::exception("GeneratorAnalyser") << "Sim Neutrinos/MCTruth vector size mismatch" << std::endl;

   return generatorTruth;
}

#endif