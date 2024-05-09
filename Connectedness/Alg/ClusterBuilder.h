#ifndef _ClusterBuilder_h_
#define _ClusterBuilder_h_

#include <iostream>
#include <fstream>
#include <algorithm>

#include "TH2D.h"
#include "TCanvas.h"
#include <TROOT.h>
#include <TStyle.h>

namespace Connectedness {

struct Cluster {

int ID;
   // Lists of bins in this cluster
   std::vector<int> bins_x;
   std::vector<int> bins_y;
};

class ClusterBuilder {

   public: 

      ClusterBuilder();
      ClusterBuilder(bool Draw, std::string DisplayDir = "");
      ~ClusterBuilder();

      // Setup functions

      void LoadDeadWireMaps(std::string Dir = "../");

      void SetThreshold(double Threshold = 1.8);
      void SetOffsets(int XOffset, int YOffset);
      void SetSeachArea(int XMax, int YMax);

      // Event processing functions

      void ReadData(std::vector<int> Channel, std::vector<int> Tick, std::vector<double> Signal, std::string RSE = "");

      // Returns the ID of the cluster (check this against list of existing cluster IDs to
      // see if some merging has taken place). Returns -1,-1 if seed landed on empty bin 
      std::pair<int,int> MakeCluster(int SeedChannel, int SeedTick, int Index);

      std::vector<Cluster> GetClusters();

      // Plane quality checks

      // Check seeds are not separated by dead wires
      bool SeedDeadWireCheck(std::vector<int> SeedChannel, std::vector<int> SeedTick, int Plane);


      // Empty the list of clusters, do before running different set of clusters for same event
      void ClearClusters();

      // Empty the list of clusters and delete histograms, do before reading a new event
      void Reset();
        
      void SetDisplayDir(std::string Dir);

   private:

      bool drawEverything = false;
      int xOffset = 0;
      int yOffset = 20;

      std::vector<int> DeadChannels_Plane0;
      std::vector<int> DeadChannels_Plane1;
      std::vector<int> DeadChannels_Plane2;

      double threshold = 1.8;

      TH2D *hRaw = nullptr;
      TH2D *hBinary = nullptr;
      TH2D *hClustered = nullptr;

      TCanvas *c = nullptr;

      std::vector<Cluster> clusters;

      std::pair<int,int> FindNearestOccupiedBin(TH2D *Hist, int X, int Y);
      int maxSearchX = 2;
      int maxSearchY = 15;

       void DeadWireFill(int Plane);

       void Focus();

       std::string displayDir = "";

   public:

      // pass = -1 (selection not applicable) , pass = 0 = not selected plane , pass = 1 = selected plane
      void DrawRaw(std::string RSE = "", int Pass = -1);
      void DrawBinary(std::string RSE = "", int Pass = -1);
      void DrawClustered(std::string RSE = "", int Plane = -1, int Pass = -1);   
};

}

#endif