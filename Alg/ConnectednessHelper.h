#ifndef _ConnectednessHelper_h_
#define _ConnectednessHelper_h_

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/Wire.h"

#include "ubana/HyperonProduction/Connectedness/Alg/ClusterBuilder.h"
#include "PositionToWire.h"

#include "TVector3.h"

namespace hyperon {

struct Wiremap {

   int View;

   std::vector<int> Channel;
   std::vector<int> Tick;
   std::vector<double> Signal;

   void LoadActivity(std::vector<art::Ptr<recob::Wire>> wires){

      Channel.clear();
      Tick.clear();
      Signal.clear();

      // Iterate through all of the wires, record signal at every tick with nonzero signal
      for(const art::Ptr<recob::Wire> &wire : wires){

         if(wire->View() == View){

            // Get regions of interest        
            unsigned int NROI = wire->SignalROI().n_ranges();
            for(size_t i_roi=0; i_roi<NROI; ++i_roi){

               // Region of tick space with nonzero activity
               recob::Wire::RegionsOfInterest_t::datarange_t const& range = wire->SignalROI().range(i_roi);

               // Iterate through the ticks in this ROI, record signal
               unsigned int thisTick = range.begin_index();

               while(thisTick < range.end_index()){

                  Channel.push_back(wire->Channel());
                  Tick.push_back(thisTick);
                  Signal.push_back(wire->Signal().at(thisTick));

                  thisTick++;

               } // while(thisTick < range.end_index

            } // loop over ROI

         } // if(view == View)

      } // loop over wires

   }

};

// Result performing the test on a single set of three tracks

struct CTSingleOutcome {

   int Plane;

   std::vector<int> SeedIndexes;
   std::vector<int> OutputIndexes;
   std::vector<int> OutputSizes;

   std::vector<int> SeedChannels;
   std::vector<int> SeedTicks;

};

// List of results from all the combinations of three tracks in event
// for all three planes

struct CTOutcome {

   std::vector<std::vector<int>> SeedIndexesPlane0;
   std::vector<std::vector<int>> OutputIndexesPlane0;
   std::vector<std::vector<int>> OutputSizesPlane0;
   std::vector<std::vector<int>> SeedChannelsPlane0;
   std::vector<std::vector<int>> SeedTicksPlane0;

   std::vector<std::vector<int>> SeedIndexesPlane1;
   std::vector<std::vector<int>> OutputIndexesPlane1;
   std::vector<std::vector<int>> OutputSizesPlane1;
   std::vector<std::vector<int>> SeedChannelsPlane1;
   std::vector<std::vector<int>> SeedTicksPlane1;

   std::vector<std::vector<int>> SeedIndexesPlane2;
   std::vector<std::vector<int>> OutputIndexesPlane2;
   std::vector<std::vector<int>> OutputSizesPlane2;
   std::vector<std::vector<int>> SeedChannelsPlane2;
   std::vector<std::vector<int>> SeedTicksPlane2;

};

class ConnectednessHelper {

   public:

      ConnectednessHelper(bool Draw);
      void LoadWireActivity(std::vector<art::Ptr<recob::Wire>> wires);
      void AddStartPositions(std::vector<TVector3> positions);
  
      std::vector<CTSingleOutcome> RunTest();

      std::vector<CTSingleOutcome> RunClustering(std::vector<int> indexes);        

      CTOutcome PrepareAndTestEvent(art::Event const& e,std::string wirelabel,std::vector<TVector3> trackstarts);

   private: 
     
      Connectedness::ClusterBuilder CPlane0;
      Connectedness::ClusterBuilder CPlane1;
      Connectedness::ClusterBuilder CPlane2;

      bool draw = false;
      std::string rse;

      Wiremap WMPlane0;
      Wiremap WMPlane1;
      Wiremap WMPlane2;

      // Starting positions of tracks in 3D (do not correct for SC)
      std::vector<TVector3> Positions3D;   
   
      // Starting positions of tracks in channel/tick space  
      std::vector<std::pair<int,int>> PositionsPlane0;
      std::vector<std::pair<int,int>> PositionsPlane1;
      std::vector<std::pair<int,int>> PositionsPlane2;
};

}

#endif