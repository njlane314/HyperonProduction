#ifndef _ConnectednessHelper_cxx_
#define _ConnectednessHelper_cxx_

#include "ConnectednessHelper.h"

using namespace hyperon;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

ConnectednessHelper::ConnectednessHelper(bool draw) :
   CPlane0(draw,"Plane0"),
   CPlane1(draw,"Plane1"),
   CPlane2(draw,"Plane2")
{
   Draw = draw;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConnectednessHelper::LoadWireActivity(std::vector<art::Ptr<recob::Wire>> wires){

   WMPlane0.View = 0;
   WMPlane0.LoadActivity(wires);
   CPlane0.Reset();
   CPlane0.ReadData(WMPlane0.Channel,WMPlane0.Tick,WMPlane0.Signal);

   WMPlane1.View = 1;
   WMPlane1.LoadActivity(wires);
   CPlane1.Reset();
   CPlane1.ReadData(WMPlane1.Channel,WMPlane1.Tick,WMPlane1.Signal);

   WMPlane2.View = 2;
   WMPlane2.LoadActivity(wires);
   CPlane2.Reset();
   CPlane2.ReadData(WMPlane2.Channel,WMPlane2.Tick,WMPlane2.Signal);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConnectednessHelper::AddStartPositions(std::vector<TVector3> positions){

   Positions3D.clear();
   PositionsPlane0.clear();
   PositionsPlane1.clear();
   PositionsPlane2.clear();

   Positions3D = positions;

   for(size_t i=0;i<Positions3D.size();i++){

      PositionsPlane0.push_back(std::make_pair(U_wire(Positions3D.at(i)),tick(Positions3D.at(i))));
      PositionsPlane1.push_back(std::make_pair(V_wire(Positions3D.at(i)),tick(Positions3D.at(i))));
      PositionsPlane2.push_back(std::make_pair(Y_wire(Positions3D.at(i)),tick(Positions3D.at(i))));

   }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<CTSingleOutcome> ConnectednessHelper::RunClustering(std::vector<int> indexes){


   // Setup the output structs   

   CTSingleOutcome OutcomePlane0;
   OutcomePlane0.Plane = 0;
   OutcomePlane0.SeedIndexes = indexes;  

   CTSingleOutcome OutcomePlane1;
   OutcomePlane1.Plane = 1;
   OutcomePlane1.SeedIndexes = indexes;  

   CTSingleOutcome OutcomePlane2;
   OutcomePlane2.Plane = 2;
   OutcomePlane2.SeedIndexes = indexes;  

   for(size_t i=0;i<indexes.size();i++){

      int index = indexes.at(i);

      OutcomePlane0.SeedChannels.push_back(PositionsPlane0.at(index).first);
      OutcomePlane0.SeedTicks.push_back(PositionsPlane0.at(index).second);

      OutcomePlane1.SeedChannels.push_back(PositionsPlane1.at(index).first);
      OutcomePlane1.SeedTicks.push_back(PositionsPlane1.at(index).second);

      OutcomePlane2.SeedChannels.push_back(PositionsPlane2.at(index).first);
      OutcomePlane2.SeedTicks.push_back(PositionsPlane2.at(index).second);

   }

   // Try generating clusters

   // Plane0 //

   std::string labelPlane0 = rse + "Plane0_";

   for(size_t i=0;i<indexes.size();i++){
      int index = indexes.at(i);

      std::pair<int,int> id_and_size =  CPlane0.MakeCluster(PositionsPlane0.at(index).first,PositionsPlane0.at(index).second,index);

      if(i < indexes.size()-1) labelPlane0 += std::to_string(indexes.at(i)) + "_";
      else labelPlane0 += std::to_string(indexes.at(i));

      OutcomePlane0.OutputIndexes.push_back(id_and_size.first);
      OutcomePlane0.OutputSizes.push_back(id_and_size.second);

   }


   if(Draw) CPlane0.DrawClustered(labelPlane0,0,-1);

   CPlane0.ClearClusters();

   // Plane1 //

   std::string labelPlane1 = rse + "Plane1_";

   for(size_t i=0;i<indexes.size();i++){
      int index = indexes.at(i);

      std::pair<int,int> id_and_size =  CPlane1.MakeCluster(PositionsPlane1.at(index).first,PositionsPlane1.at(index).second,index);

      if(i < indexes.size()-1) labelPlane1 += std::to_string(indexes.at(i)) + "_";
      else labelPlane1 += std::to_string(indexes.at(i));

      OutcomePlane1.OutputIndexes.push_back(id_and_size.first);
      OutcomePlane1.OutputSizes.push_back(id_and_size.second);

   }

   if(Draw) CPlane1.DrawClustered(labelPlane1,1,-1);

   CPlane1.ClearClusters();

   // Plane2 //

   std::string labelPlane2 = rse + "Plane2_";

   for(size_t i=0;i<indexes.size();i++){
      int index = indexes.at(i);

      std::pair<int,int> id_and_size =  CPlane2.MakeCluster(PositionsPlane2.at(index).first,PositionsPlane2.at(index).second,index);

      if(i < indexes.size()-1) labelPlane2 += std::to_string(indexes.at(i)) + "_";
      else labelPlane2 += std::to_string(indexes.at(i));

      OutcomePlane2.OutputIndexes.push_back(id_and_size.first);
      OutcomePlane2.OutputSizes.push_back(id_and_size.second);

   }

   if(Draw) CPlane2.DrawClustered(labelPlane2,2,-1);

   CPlane2.ClearClusters();

   return {OutcomePlane0,OutcomePlane1,OutcomePlane2};

}        

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<CTSingleOutcome> ConnectednessHelper::RunTest(){

   std::vector<CTSingleOutcome> AllOutcomes;

   // Need at least three tracks for the test
   if(Positions3D.size() < 3) return AllOutcomes;

   // Try all 3 track combinations possible and record outcomes

   int nseeds = static_cast<int>(Positions3D.size());

   for(int i=0;i<nseeds;i++){
      for(int j=i+1;j<nseeds;j++){
         for(int k=j+1;k<nseeds;k++){
            //std::cout << "Running clustering on combination " << i << j << k << std::endl;
            std::vector<CTSingleOutcome> Outcomes = RunClustering({i,j,k});

            AllOutcomes.push_back(Outcomes.at(0));
            AllOutcomes.push_back(Outcomes.at(1));
            AllOutcomes.push_back(Outcomes.at(2));

         }
      }
   }

   return AllOutcomes;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

CTOutcome ConnectednessHelper::PrepareAndTestEvent(art::Event const& e,std::string wirelabel,std::vector<TVector3> trackstarts){
 
   CTOutcome theOutcome;
   
   rse = std::to_string(e.run()) + "_" + std::to_string(e.subRun()) + "_" + std::to_string(e.event());

   if(trackstarts.size() < 3) return theOutcome;

   //std::cout << "Running CT Test on event with " << trackstarts.size() << " tracks" << std::endl;

   art::Handle<std::vector<recob::Wire>> Handle_Wire;
   std::vector<art::Ptr<recob::Wire>> Vect_Wire;

   if(!e.getByLabel(wirelabel,Handle_Wire)) 
      throw cet::exception("ConnectednessHelper") << "Wire data product not found!" << std::endl;

   art::fill_ptr_vector(Vect_Wire,Handle_Wire);

   //std::cout << "Loading wire activity" << std::endl;
   LoadWireActivity(Vect_Wire);
   //std::cout << "Adding track start positions" << std::endl;
   AddStartPositions(trackstarts);

   // Vector containing the result for each combintation of three tracks
   //std::cout << "Running test" << std::endl;
   std::vector<CTSingleOutcome> Outcomes = RunTest(); 

   // iterate over each combination of three tracks, store the result
   for(size_t i=0;i<Outcomes.size();i++){

      CTSingleOutcome this_Outcome = Outcomes.at(i);

      std::vector<int> this_SeedIndexes = this_Outcome.SeedIndexes;
      std::vector<int> this_OutputIndexes = this_Outcome.OutputIndexes;
      std::vector<int> this_OutputSizes = this_Outcome.OutputSizes;
      std::vector<int> this_SeedChannels = this_Outcome.SeedChannels;
      std::vector<int> this_SeedTicks = this_Outcome.SeedTicks;

      if(this_Outcome.Plane == 0){
         theOutcome.SeedIndexesPlane0.push_back(this_SeedIndexes);
         theOutcome.OutputIndexesPlane0.push_back(this_OutputIndexes);
         theOutcome.OutputSizesPlane0.push_back(this_OutputSizes);
         theOutcome.SeedChannelsPlane0.push_back(this_SeedChannels);
         theOutcome.SeedTicksPlane0.push_back(this_SeedTicks);
      } 
      else if(this_Outcome.Plane == 1){
         theOutcome.SeedIndexesPlane1.push_back(this_SeedIndexes);
         theOutcome.OutputIndexesPlane1.push_back(this_OutputIndexes);
         theOutcome.OutputSizesPlane1.push_back(this_OutputSizes);
         theOutcome.SeedChannelsPlane1.push_back(this_SeedChannels);
         theOutcome.SeedTicksPlane1.push_back(this_SeedTicks);
      } 
      else if(this_Outcome.Plane == 2){
         theOutcome.SeedIndexesPlane2.push_back(this_SeedIndexes);
         theOutcome.OutputIndexesPlane2.push_back(this_OutputIndexes);
         theOutcome.OutputSizesPlane2.push_back(this_OutputSizes);
         theOutcome.SeedChannelsPlane2.push_back(this_SeedChannels);
         theOutcome.SeedTicksPlane2.push_back(this_SeedTicks);
      } 
   }

   return theOutcome;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
