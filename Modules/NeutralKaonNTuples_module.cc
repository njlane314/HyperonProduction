#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"

#include "ubana/HyperonProduction/Alg/ConnectednessHelper.h"

#include "ubana/HyperonProduction/Modules/SubModules/GeneratorAnalyser.h"
#include "ubana/HyperonProduction/Modules/SubModules/ParticleTrackerAnalyser.h"
#include "ubana/HyperonProduction/Modules/SubModules/ReconstructionAnalyser.h"

namespace hyperon 
{
   class NeutralKaonNTuples;
}

class hyperon::NeutralKaonNTuples : public art::EDAnalyzer 
{
   public:
      explicit NeutralKaonNTuples(fhicl::ParameterSet const& p);

      NeutralKaonNTuples(NeutralKaonNTuples const&) = delete;
      NeutralKaonNTuples(NeutralKaonNTuples&&) = delete;
      NeutralKaonNTuples& operator=(NeutralKaonNTuples const&) = delete;
      NeutralKaonNTuples& operator=(NeutralKaonNTuples&&) = delete;

      void analyze(art::Event const& e) override;

      void beginJob() override;
      void endJob() override;

      void FinishEvent();

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:
      TTree *OutputTree;
      TTree *MetaTree;

      unsigned int t_EventIdentifier;
      int t_Run, t_Subrun, t_Event;

      // Generator event variables

      double t_gen_Weight = 1.0;

      int t_gen_nMCTruths = 0;	
      int t_gen_nMCTruthsInTPC = 0;	

      std::vector<SimParticle> t_gen_Neutrinos;

      std::vector<std::string> t_gen_CCNC;
      std::vector<std::string> t_gen_Mode;

      std::vector<double> t_gen_W;
      std::vector<double> t_gen_X;
      std::vector<double> t_gen_Y;
      std::vector<double> t_gen_QSqr; 
      std::vector<double> t_gen_Pt;
      std::vector<double> t_gen_Theta;

      // Particle tracker variables

      std::vector<bool> t_event_inActiveTPC;

      std::vector<int> t_event_nHyperons;
      std::vector<int> t_event_nKaonShorts ;
      std::vector<int> t_event_nPions;

      std::vector<double> t_event_PrimaryVerticesX;
      std::vector<double> t_event_PrimaryVerticesY;
      std::vector<double> t_event_PrimaryVerticesZ;
      
      std::vector<bool> t_event_hasKaonShort;
      std::vector<bool> t_event_hasKaonShortCharged;
      std::vector<bool> t_event_hasKaonShortChargedSignal;

      std::vector<SimParticle> t_event_PrimaryLeptons; 

      std::vector<SimParticle> t_event_PrimaryKaonShorts;
      std::vector<SimParticle> t_event_KaonShortChargedDaughters;

      bool t_event_eventHasKaonShort;
      bool t_event_eventHasHyperon;
      bool t_event_eventHasNeutronScatter;

      // Track reconstruction variables

      int t_reco_nPrimaryDaughters;
      int t_reco_nPrimaryTrackDaughters;
      int t_reco_nPrimaryShowerDaughters;

      std::vector<RecoParticle> t_reco_TrackPrimaryDaughters;
      std::vector<RecoParticle> t_reco_ShowerPrimaryDaughters;   

      TVector3 t_reco_InteractionVertex;
      bool t_reco_KoanShortChargedGoodReco;

      std::vector<RecoParticle> t_reco_RepassTrackPrimaryDaughters;
      std::vector<RecoParticle> t_reco_RepassShowerPrimaryDaughters;

      // Connectedness test variables

      std::vector<std::vector<int>> t_conn_SeedIndexesPlane0;
      std::vector<std::vector<int>> t_conn_OutputIndexesPlane0;
      std::vector<std::vector<int>> t_conn_OutputSizesPlane0;
      std::vector<std::vector<int>> t_conn_SeedChannelsPlane0;
      std::vector<std::vector<int>> t_conn_SeedTicksPlane0;

      std::vector<std::vector<int>> t_conn_SeedIndexesPlane1;
      std::vector<std::vector<int>> t_conn_OutputIndexesPlane1;
      std::vector<std::vector<int>> t_conn_OutputSizesPlane1;
      std::vector<std::vector<int>> t_conn_SeedChannelsPlane1;
      std::vector<std::vector<int>> t_conn_SeedTicksPlane1;

      std::vector<std::vector<int>> t_conn_SeedIndexesPlane2;
      std::vector<std::vector<int>> t_conn_OutputIndexesPlane2;
      std::vector<std::vector<int>> t_conn_OutputSizesPlane2;
      std::vector<std::vector<int>> t_conn_SeedChannelsPlane2;
      std::vector<std::vector<int>> t_conn_SeedTicksPlane2;

      // Data variables

      std::vector<std::string> t_data_SysDials;
      std::vector<std::vector<double>> t_data_SysWeights;

      // Meta variables

      int m_nEvents;

      int m_nKaonShortChargedSignal;         
      int m_nKaonShortChargedGoodReco;

      double m_POT = 0; // total POT of the sample

      // Configuration parameters

      bool f_GetGeneratorInfo;
      bool f_GetG4Info;
      bool f_GetRecoInfo;
      bool f_GetConnInfo;

      fhicl::ParameterSet f_Generator;
      fhicl::ParameterSet f_G4;
      fhicl::ParameterSet f_Reco;

      std::vector<fhicl::ParameterSet> f_RecoRepass;
      std::string f_WireLabel;
      std::vector<art::InputTag> f_WeightLabels;
      std::string f_POTSummaryLabel;

      bool f_ParticleGun = false;
      bool f_IsData;
      bool f_Debug = false;

      ConnectednessHelper Conn_Helper;
};

hyperon::NeutralKaonNTuples::NeutralKaonNTuples(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_GetGeneratorInfo(p.get<bool>("GetGeneratorInfo", true)),   
   f_GetG4Info(p.get<bool>("GetG4Info", true)),   
   f_GetRecoInfo(p.get<bool>("GetRecoInfo", true)),   
   f_GetConnInfo(p.get<bool>("GetConnInfo", true)),   
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Reco(p.get<fhicl::ParameterSet>("Reco")),
   f_RecoRepass(p.get<std::vector<fhicl::ParameterSet>>("RecoRepass", {})),
   f_WireLabel(p.get<std::string>("WireLabel")),
   f_WeightLabels(p.get<std::vector<art::InputTag>>("WeightCalculators", {})),
   f_POTSummaryLabel(p.get<std::string>("POTSummaryLabel")),
   f_ParticleGun(p.get<bool>("ParticleGun", false)),
   f_IsData(p.get<bool>("IsData")),
   f_Debug(p.get<bool>("Debug", false)),
   Conn_Helper(p.get<bool>("DrawConnectedness", false))   // ,
{
   if(f_WeightLabels.size()){
      std::cout << "Getting weights from data products with tags:" << std::endl;
      for(size_t i = 0; i < f_WeightLabels.size(); i++) std::cout << f_WeightLabels.at(i) << std::endl;
   }

   if(f_Reco.get<bool>("IncludeCosmics", false) && f_GetConnInfo)
      std::cout << std::endl << "HyperonNTuples WARNING: Requesting connectedness information with cosmics included. This will take a very long time." << std::endl 
                             << "Set GetConnInfo fhicl parameter to false disable connectedness" << std::endl << std::endl;
      
}

void hyperon::NeutralKaonNTuples::analyze(art::Event const& e)
{
   if(f_Debug) std::cout << "-------- New Event --------" << std::endl;

   t_EventIdentifier = e.id().event();
   t_Run = e.run();
   t_Subrun = e.subRun();
   t_Event = e.event();

   // Initialising generator variables

   t_gen_Weight = 1.0;
   
   t_gen_nMCTruths = 0;
   t_gen_nMCTruthsInTPC = 0;

   t_gen_Neutrinos.clear();

   t_gen_CCNC.clear();
   t_gen_Mode.clear();

   t_gen_W.clear();
   t_gen_X.clear();
   t_gen_Y.clear();
   t_gen_QSqr.clear();
   t_gen_Pt.clear();
   t_gen_Theta.clear();

   // Initialising particle tracker variables

   t_event_inActiveTPC.clear();

   t_event_nHyperons.clear();
   t_event_nKaonShorts .clear();
   t_event_nPions.clear();

   t_event_PrimaryVerticesX.clear();
   t_event_PrimaryVerticesY.clear();
   t_event_PrimaryVerticesZ.clear();

   t_event_hasKaonShort.clear();
   t_event_hasKaonShortCharged.clear();

   t_event_PrimaryLeptons.clear();

   t_event_PrimaryKaonShorts.clear();
   t_event_KaonShortChargedDaughters.clear();

   // Initialising track reconstruction variables

   t_reco_nPrimaryDaughters = 0;
   t_reco_nPrimaryTrackDaughters = 0;
   t_reco_nPrimaryShowerDaughters = 0;

   t_reco_TrackPrimaryDaughters.clear();
   t_reco_ShowerPrimaryDaughters.clear();

   t_reco_InteractionVertex.SetXYZ(-1000, -1000, -1000);
   t_reco_KoanShortChargedGoodReco = false;

   t_reco_RepassTrackPrimaryDaughters.clear();
   t_reco_RepassShowerPrimaryDaughters.clear();

   // Initialising connectedness test variables

   t_conn_SeedIndexesPlane0.clear();
   t_conn_OutputIndexesPlane0.clear();
   t_conn_OutputSizesPlane0.clear();
   t_conn_SeedChannelsPlane0.clear();
   t_conn_SeedTicksPlane0.clear();

   t_conn_SeedIndexesPlane1.clear();
   t_conn_OutputIndexesPlane1.clear();
   t_conn_OutputSizesPlane1.clear();
   t_conn_SeedChannelsPlane1.clear();
   t_conn_SeedTicksPlane1.clear();

   t_conn_SeedIndexesPlane2.clear();
   t_conn_OutputIndexesPlane2.clear();
   t_conn_OutputSizesPlane2.clear();
   t_conn_SeedChannelsPlane2.clear();
   t_conn_SeedTicksPlane2.clear();

   // Initialising data variables

   t_data_SysDials.clear();
   t_data_SysWeights.clear();

   // Assigining variables

   if(!f_IsData && f_GetGeneratorInfo){

      if(f_Debug) std::cout << "Getting generator information..." << std::endl;

      GeneratorAnalyser* Generator_SM = new GeneratorAnalyser(e, f_Generator, f_ParticleGun);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();

      t_gen_Weight *= GenT.Weight;

      t_gen_nMCTruths = GenT.nMCTruths;
      t_gen_nMCTruthsInTPC = GenT.nMCTruthsInTPC;

      t_gen_Neutrinos = GenT.Neutrinos;

      t_gen_CCNC = GenT.CCNC;
      t_gen_Mode = GenT.Mode;

      t_gen_W = GenT.W;
      t_gen_X = GenT.X;
      t_gen_Y = GenT.Y;
      t_gen_QSqr = GenT.QSqr;
      t_gen_Pt = GenT.Pt;
      t_gen_Theta = GenT.Theta;

      delete Generator_SM;
   }

   if(!f_IsData && f_GetG4Info){

      if(f_Debug) std::cout << "Getting particle tracker information..." << std::endl;

      ParticleTrackerAnalyser* G4_SM = new ParticleTrackerAnalyser(e, f_G4, f_ParticleGun);
      EventTruth G4T = G4_SM->GetEventTruth();

      t_event_inActiveTPC = G4T.inActiveTPC;

      t_event_nHyperons = G4T.nHyperons;
      t_event_nKaonShorts  = G4T.nKaonShorts;
      t_event_nPions = G4T.nPions;

      t_event_hasKaonShort = G4T.hasKaonShort;
      t_event_hasKaonShortCharged = G4T.hasKaonShortCharged;

      t_event_PrimaryLeptons = G4T.PrimaryLeptons;

      t_event_PrimaryKaonShorts = G4T.PrimaryKaonShorts;
      t_event_KaonShortChargedDaughters = G4T.KaonShortChargedDaughters;

      t_event_eventHasKaonShort = G4T.eventHasKaonShort;
      t_event_eventHasHyperon = G4T.eventHasHyperon;
      t_event_eventHasNeutronScatter = G4T.eventHasNeutronScatter;

      t_event_hasKaonShortChargedSignal.resize(t_gen_nMCTruths);

      for(int i = 0; i < t_gen_nMCTruths; i++){
         t_event_hasKaonShortChargedSignal[i] = t_event_inActiveTPC.at(i) && t_event_hasKaonShortCharged.at(i);
      }

      delete G4_SM;
   }

   if(f_GetRecoInfo){

      if(f_Debug) std::cout << "Getting reconstruction information..." << std::endl;

      ReconstructionAnalyser* Reco_SM = new ReconstructionAnalyser(e, f_IsData, f_Reco, f_ParticleGun);
      Reco_SM->PrepareInfo();
      Reco_SM->SetIndices(t_event_hasKaonShortChargedSignal);
      RecoData RecoD =  Reco_SM->GetInfo();   

      t_reco_nPrimaryDaughters = RecoD.NPrimaryDaughters;
      t_reco_nPrimaryTrackDaughters = RecoD.NPrimaryTrackDaughters;
      t_reco_nPrimaryShowerDaughters = RecoD.NPrimaryShowerDaughters;

      t_reco_TrackPrimaryDaughters = RecoD.TrackPrimaryDaughters;
      t_reco_ShowerPrimaryDaughters = RecoD.ShowerPrimaryDaughters;

      t_reco_InteractionVertex = RecoD.RecoPrimaryVertex;
      t_reco_KoanShortChargedGoodReco = RecoD.GoodReco;

      // Results of connectedness test on different combinations of tracks

      if(f_Debug) std::cout << "Performing connectedness tests..." << std::endl;

      if(f_GetConnInfo){

         CTOutcome ConnData = Conn_Helper.PrepareAndTestEvent(e, f_WireLabel, RecoD.TrackStarts);   

         t_conn_SeedIndexesPlane0 = ConnData.SeedIndexesPlane0;
         t_conn_OutputIndexesPlane0 = ConnData.OutputIndexesPlane0;
         t_conn_OutputSizesPlane0 = ConnData.OutputSizesPlane0;
         t_conn_SeedChannelsPlane0 = ConnData.SeedChannelsPlane0;
         t_conn_SeedTicksPlane0 = ConnData.SeedTicksPlane0;
         t_conn_SeedIndexesPlane1 = ConnData.SeedIndexesPlane1;
         t_conn_OutputIndexesPlane1 = ConnData.OutputIndexesPlane1;
         t_conn_OutputSizesPlane1 = ConnData.OutputSizesPlane1;
         t_conn_SeedChannelsPlane1 = ConnData.SeedChannelsPlane1;
         t_conn_SeedTicksPlane1 = ConnData.SeedTicksPlane1;
         t_conn_SeedIndexesPlane2 = ConnData.SeedIndexesPlane2;
         t_conn_OutputIndexesPlane2 = ConnData.OutputIndexesPlane2;
         t_conn_OutputSizesPlane2 = ConnData.OutputSizesPlane2;
         t_conn_SeedChannelsPlane2 = ConnData.SeedChannelsPlane2;
         t_conn_SeedTicksPlane2 = ConnData.SeedTicksPlane2;
      }

      delete Reco_SM;

      if(f_RecoRepass.size()){
        std::cout << "Getting repassed information..." << std::endl;
         for(size_t i = 0; i < f_RecoRepass.size(); i++){
            ReconstructionAnalyser* Reco_SM_Repass = new ReconstructionAnalyser(e, f_IsData, f_RecoRepass.at(i));
            Reco_SM_Repass->PrepareInfo();
         
            RecoData RecoD_Repass =  Reco_SM_Repass->GetInfo();   
            t_reco_RepassTrackPrimaryDaughters.insert(t_reco_RepassTrackPrimaryDaughters.end(), RecoD_Repass.TrackPrimaryDaughters.begin(), RecoD_Repass.TrackPrimaryDaughters.end());
            t_reco_RepassShowerPrimaryDaughters.insert(t_reco_RepassShowerPrimaryDaughters.end(), RecoD_Repass.ShowerPrimaryDaughters.begin(), RecoD_Repass.ShowerPrimaryDaughters.end());

            delete Reco_SM_Repass;
         }
      }
   }

   if(!f_IsData){

      std::vector<std::map<std::string, std::vector<double>>> theweightmap(t_gen_nMCTruths); 

      for(size_t i = 0; i < f_WeightLabels.size(); i++){

         std::cout << "Getting new weight products with label " << f_WeightLabels.at(i) << std::endl;

         art::Handle<std::vector<evwgh::MCEventWeight>> Handle_EventWeight;
         std::vector<art::Ptr<evwgh::MCEventWeight>> Vect_EventWeight;

         if(!e.getByLabel(f_WeightLabels.at(i), Handle_EventWeight)) 
            throw cet::exception("NeutralKaonNTuples") << "No EventWeight Found!" << std::endl;

         art::fill_ptr_vector(Vect_EventWeight, Handle_EventWeight);

         if(!Vect_EventWeight.size())
            throw cet::exception("NeutralKaonNTuples") << "Weight vector empty!" << std::endl;

         if(Vect_EventWeight.size() != (size_t)t_gen_nMCTruths)
            throw cet::exception("NeutralKaonNTuples") << "Weight vector size != nMCTruths" << std::endl;

         for(size_t j = 0; j < Vect_EventWeight.size(); j++){       

            std::cout << "Getting weights for truth " << j << std::endl;

            std::map<std::string, std::vector<double>> theWeights = Vect_EventWeight.at(j)->fWeight;
            std::map<std::string, std::vector<double>>::iterator it;

            for(it = theWeights.begin(); it != theWeights.end(); it++){

               if(it->first ==  "empty") continue;

               bool dial_found=false;

               std::map<std::string, std::vector<double>>::iterator it2;
               for(it2 = theweightmap.at(j).begin(); it2 != theweightmap.at(j).end(); it2++){
                  if(it->first == it2->first){
                     dial_found = true;
                     theweightmap.at(j)[it->first].insert(theweightmap.at(j)[it->first].end(), it->second.begin(), it->second.end());
                  }
               }

               if(!dial_found)
                  theweightmap.at(j)[it->first] = it->second;

            }
         } 
      }

      if(theweightmap.size()){
         std::map<std::string,std::vector<double>>::iterator it;
         for(it = theweightmap.at(0).begin(); it != theweightmap.at(0).end(); it++){
            std::cout << "Organising weights for dial " << it->first << std::endl;

            t_data_SysDials.push_back(it->first);
            t_data_SysWeights.push_back(it->second);

            for(size_t i = 1; i < theweightmap.size(); i++){
               if(theweightmap.at(i).find(it->first) == theweightmap.at(i).end()) 
                  throw cet::exception("NeutralKaonNTuples") << "Dial " << it->first << " not found in weights for MC truth " << i << std::endl;

               if(theweightmap.at(i)[it->first].size() != t_data_SysWeights.back().size())
                  throw cet::exception("NeutralKaonNTuples") << "Dial " << it->first << " weight vector mismatch" << std::endl;   

               for(size_t j = 0; j < t_data_SysWeights.back().size(); j++)
                  t_data_SysWeights.back().at(j) *= theweightmap.at(i)[it->first].at(j);
            }                  
         }
      }
   }

   FinishEvent();
}

void hyperon::NeutralKaonNTuples::FinishEvent()
{
   if(f_Debug) std::cout << "Finishing event..." << std::endl;

   OutputTree->Fill();

   m_nEvents++;

   if(std::find(t_event_hasKaonShortChargedSignal.begin(), t_event_hasKaonShortChargedSignal.end(), true) != t_event_hasKaonShortChargedSignal.end()) m_nKaonShortChargedSignal++;
   if(t_reco_KoanShortChargedGoodReco) m_nKaonShortChargedGoodReco++;

   if(f_Debug) std::cout << "Finished event!" << std::endl;
}	

void hyperon::NeutralKaonNTuples::beginJob()
{
   if(f_Debug) std::cout << "Beginning job..." << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   // Output tree

   OutputTree = tfs->make<TTree>("OutputTree", "Truth Info Tree");

   OutputTree->Branch("IsData", &f_IsData);

   OutputTree->Branch("t_EventIdentifier", &t_EventIdentifier);
   OutputTree->Branch("t_Run", &t_Run);
   OutputTree->Branch("t_Subrun", &t_Subrun);
   OutputTree->Branch("t_Event", &t_Event);

   // Assign generator variable trees

   OutputTree->Branch("Weight", &t_gen_Weight);

   OutputTree->Branch("nMCTruths", &t_gen_nMCTruths);
   OutputTree->Branch("nMCTruthsInTPC", &t_gen_nMCTruthsInTPC);

   OutputTree->Branch("CCNC", "vector<string>", &t_gen_CCNC);
   OutputTree->Branch("Mode", "vector<string>", &t_gen_Mode);

   OutputTree->Branch("W", &t_gen_W);
   OutputTree->Branch("X", &t_gen_X);
   OutputTree->Branch("Y", &t_gen_Y);
   OutputTree->Branch("QSqr", &t_gen_QSqr);
   OutputTree->Branch("Pt", &t_gen_Pt);
   OutputTree->Branch("Theta", &t_gen_Theta);

   OutputTree->Branch("Neutrinos", "vector<SimParticle>", &t_gen_Neutrinos);

   // Assign particle tracker variable trees

   OutputTree->Branch("inActiveTPC", "vector<bool>", &t_event_inActiveTPC);

   OutputTree->Branch("nHyperons", "vector<int>", &t_event_nHyperons);
   OutputTree->Branch("nKaonShorts", "vector<int>", &t_event_nKaonShorts);
   OutputTree->Branch("nPions", "vector<int>", &t_event_nPions);

   OutputTree->Branch("PrimaryVerticesX", "vector<double>", &t_event_PrimaryVerticesX);
   OutputTree->Branch("PrimaryVerticesY", "vector<double>", &t_event_PrimaryVerticesY);
   OutputTree->Branch("PrimaryVerticesZ", "vector<double>", &t_event_PrimaryVerticesZ);

   OutputTree->Branch("hasKaonShort", "vector<bool>", &t_event_hasKaonShort);
   OutputTree->Branch("hasKaonShortCharged", "vector<bool>", &t_event_hasKaonShortCharged);

   OutputTree->Branch("PrimaryLeptons", "vector<SimParticle>", &t_event_PrimaryLeptons);
   OutputTree->Branch("PrimaryKaonShorts", "vector<SimParticle>", &t_event_PrimaryKaonShorts);
   OutputTree->Branch("KaonShortChargedDaughters", "vector<SimParticle>", &t_event_KaonShortChargedDaughters);

   // Assign reconstruced variable trees

   OutputTree->Branch("nPrimaryDaughters", &t_reco_nPrimaryDaughters);
   OutputTree->Branch("nPrimaryTrackDaughters", &t_reco_nPrimaryTrackDaughters);
   OutputTree->Branch("nPrimaryShowerDaughters", &t_reco_nPrimaryShowerDaughters);

   OutputTree->Branch("TrackPrimaryDaughters", "vector<RecoParticle>", &t_reco_TrackPrimaryDaughters);
   OutputTree->Branch("ShowerPrimaryDaughters", "vector<RecoParticle>", &t_reco_ShowerPrimaryDaughters);

   OutputTree->Branch("RepassTracklikePrimaryDaughters", "vector<RecoParticle>", &t_reco_RepassTrackPrimaryDaughters);
   OutputTree->Branch("RepassShowerlikePrimaryDaughters", "vector<RecoParticle>", &t_reco_RepassShowerPrimaryDaughters);

   OutputTree->Branch("RecoInteractionVretex", "TVector3", &t_reco_InteractionVertex);

   OutputTree->Branch("KaonShortChargedGoodReco", &t_reco_KoanShortChargedGoodReco);

   // Assign connectedness test variable trees

   OutputTree->Branch("ConnSeedIndexesPlane0", &t_conn_SeedIndexesPlane0);
   OutputTree->Branch("ConnOutputIndexesPlane0", &t_conn_OutputIndexesPlane0);
   OutputTree->Branch("ConnOutputSizesPlane0", &t_conn_OutputSizesPlane0);
   OutputTree->Branch("ConnSeedChannelsPlane0", &t_conn_SeedChannelsPlane0);
   OutputTree->Branch("ConnSeedTicksPlane0", &t_conn_SeedTicksPlane0);
   OutputTree->Branch("ConnSeedIndexesPlane1", &t_conn_SeedIndexesPlane1);
   OutputTree->Branch("ConnOutputIndexesPlane1", &t_conn_OutputIndexesPlane1);
   OutputTree->Branch("ConnOutputSizesPlane1", &t_conn_OutputSizesPlane1);
   OutputTree->Branch("ConnSeedChannelsPlane1", &t_conn_SeedChannelsPlane1);
   OutputTree->Branch("ConnSeedTicksPlane1", &t_conn_SeedTicksPlane1);
   OutputTree->Branch("ConnSeedIndexesPlane2", &t_conn_SeedIndexesPlane2);
   OutputTree->Branch("ConnOutputIndexesPlane2", &t_conn_OutputIndexesPlane2);
   OutputTree->Branch("ConnOutputSizesPlane2", &t_conn_OutputSizesPlane2);
   OutputTree->Branch("ConnSeedChannelsPlane2", &t_conn_SeedChannelsPlane2);
   OutputTree->Branch("ConnSeedTicksPlane2", &t_conn_SeedTicksPlane2);

   OutputTree->Branch("SysDials", &t_data_SysDials);
   OutputTree->Branch("SysWeights", "vector<vector<double>>", &t_data_SysWeights);

   // Assign metadata variable trees
   
   m_nEvents = 0;

   m_nKaonShortChargedSignal = 0;
   m_nKaonShortChargedGoodReco = 0;

   m_POT=0;

   MetaTree=tfs->make<TTree>("MetaTree", "Metadata Info Tree");

   MetaTree->Branch("nEvents", &m_nEvents);

   MetaTree->Branch("NSignal", &m_nKaonShortChargedSignal);
   MetaTree->Branch("NGoodReco", &m_nKaonShortChargedGoodReco);

   MetaTree->Branch("POT", &m_POT);

   if(f_Debug) std::cout << "Finished begin job!" << std::endl;
}

void hyperon::NeutralKaonNTuples::endJob()
{
   MetaTree->Fill();
}

void hyperon::NeutralKaonNTuples::beginSubRun(const art::SubRun& sr)
{
   if(f_Debug) std::cout << "Getting subrun POT information..." << std::endl;

   art::Handle<sumdata::POTSummary> POTHandle;
   if(sr.getByLabel(f_POTSummaryLabel, POTHandle)) m_POT += POTHandle->totpot;	
}

void hyperon::NeutralKaonNTuples::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(hyperon::NeutralKaonNTuples)