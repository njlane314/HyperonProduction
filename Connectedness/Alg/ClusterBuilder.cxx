#ifndef _ClusterBuilder_cxx_
#define _ClusterBuilder_cxx_

#include "ClusterBuilder.h"

using namespace Connectedness;

ClusterBuilder::ClusterBuilder()
{
   drawEverything = false;
}

ClusterBuilder::ClusterBuilder(bool Draw, std::string displayDir)
{
   drawEverything = Draw;

   system("mkdir -p Displays/");

   if(drawEverything){

      std::cout << "Will Draw clusters" << std::endl;

      c = new TCanvas(("c" + displayDir).c_str(), ("c" + displayDir).c_str());

      if(displayDir != "") displayDir = displayDir;

      system(("mkdir -p Displays/" + displayDir).c_str());
      system(("mkdir -p Displays/" + displayDir + "/Pass/").c_str());
      system(("mkdir -p Displays/" + displayDir + "/Fail/").c_str());
   }
}

ClusterBuilder::~ClusterBuilder()
{
}

void ClusterBuilder::LoadDeadWireMaps(std::string Dir)
{
   std::cout << "Getting dead Channel lists" << std::endl;

   std::ifstream input_DeadChannels_Plane0(Dir + "Plane0.txt");
   std::ifstream input_DeadChannels_Plane1(Dir + "Plane1.txt");
   std::ifstream input_DeadChannels_Plane2(Dir + "Plane2.txt");

   int dead;

   while(!input_DeadChannels_Plane0.eof()){ input_DeadChannels_Plane0 >> dead; DeadChannels_Plane0.push_back(dead); }
   while(!input_DeadChannels_Plane1.eof()){ input_DeadChannels_Plane1 >> dead; DeadChannels_Plane1.push_back(dead); }
   while(!input_DeadChannels_Plane2.eof()){ input_DeadChannels_Plane2 >> dead; DeadChannels_Plane2.push_back(dead); }


   std::cout << "Found " << DeadChannels_Plane0.size() << " dead channels in Plane0" << std::endl;
   std::cout << "Found " << DeadChannels_Plane1.size() << " dead channels in Plane1" << std::endl;
   std::cout << "Found " << DeadChannels_Plane2.size() << " dead channels in Plane2" << std::endl;
}


void ClusterBuilder::SetThreshold(double Threshold)
{
   threshold = Threshold;
}


void ClusterBuilder::SetOffsets(int XOffset, int YOffset)
{
   xOffset = XOffset;
   xOffset = YOffset;
}

void ClusterBuilder::SetSeachArea(int XMax, int YMax)
{
   maxSearchX = XMax;
   maxSearchY = YMax;
}

void ClusterBuilder::ClearClusters()
{
   clusters.clear();
}

void ClusterBuilder::Reset()
{
   clusters.clear();

   if(hRaw != nullptr)   delete hRaw;
   if(hBinary != nullptr)   delete hBinary;
   if(hClustered != nullptr)   delete hClustered;
}

void ClusterBuilder::ReadData(std::vector<int> Channel, std::vector<int> Tick, std::vector<double> Signal, std::string RSE)
{
   // Channel, Tick and Signal vectors should all be of the same size
   if(Channel.size() != Tick.size() || Tick.size() != Signal.size()){ 
      std::cout << "Input Channel/Tick/Signal vectors of different sizes! Exiting" << std::endl;
      return;
   }

   int MaxChn = -10000, MinChn = 10000000;
   double MaxTk = -1e10, MinTk = 1e10;

   for(size_t i = 0; i < Channel.size(); i++){
      if(Channel.at(i) > MaxChn) MaxChn = Channel.at(i);
      if(Tick.at(i) > MaxTk) MaxTk = Tick.at(i);
      if(Channel.at(i) < MinChn) MinChn = Channel.at(i);
      if(Tick.at(i) < MinTk) MinTk = Tick.at(i);
   }

   int NChannels = MaxChn - MinChn;
   int NTicks = MaxTk - MinTk;

   // Setup the histograms
   hRaw = new TH2D(("h_Channel_vs_Tick_Raw_"+ RSE + displayDir).c_str(), "Raw Activity;Channel;Tick", NChannels, MinChn, MaxChn, NTicks, MinTk, MaxTk);
   hBinary = new TH2D(("h_Channel_vs_Tick_Binary_"+ RSE + displayDir).c_str(), "Binary Activity;Channel;Tick", NChannels, MinChn, MaxChn, NTicks, MinTk, MaxTk);

   // Fill the histograms
   for(size_t i = 0; i < Channel.size(); i++){
      hRaw->Fill(Channel.at(i), Tick.at(i), Signal.at(i));
      if(Signal.at(i) > threshold) hBinary->Fill(Channel.at(i), Tick.at(i), 1);
   }
}

std::pair<int,int> ClusterBuilder::MakeCluster(int SeedChannel, int SeedTick, int Index)
{
   // Get seed location in bin space
   int SeedChannelBin = hRaw->GetXaxis()->FindBin(SeedChannel - xOffset);
   int SeedTickBin = hRaw->GetYaxis()->FindBin(SeedTick - yOffset);

   // Check seed bin is occupied
   if(hBinary->GetBinContent(SeedChannelBin, SeedTickBin) < 1){
      std::pair<int,int> bins = FindNearestOccupiedBin(hBinary, SeedChannelBin, SeedTickBin);
      if(bins.first == -1 && bins.second == -1) return std::make_pair(-1, -1);
      else {
         SeedChannelBin = bins.first;
         SeedTickBin = bins.second;
      }
   }

   // Check if this is a bin in the list of clusters already produced, if it is
   // return the Index of that cluster

   for(size_t i_cl = 0; i_cl < clusters.size(); i_cl++){

      Cluster thisCluster = clusters.at(i_cl);

      for(size_t i_b = 0; i_b < thisCluster.bins_x.size(); i_b++){
         if(SeedChannelBin == thisCluster.bins_x.at(i_b) && SeedTickBin == thisCluster.bins_y.at(i_b)){ 
            //std::cout << "Seed already belongs to cluster " << thisCluster.Index << std::endl;
            return std::make_pair(thisCluster.ID, thisCluster.bins_x.size());
         }
      }
   }

   // CLUSTER GROWING ALG //

   // Create empty cluster struct
   Cluster C;
   C.ID = Index;

   // Add the seed bin
   C.bins_x.push_back(SeedChannelBin);
   C.bins_y.push_back(SeedTickBin);

   TH2D h_tmp = *hBinary; 

   h_tmp.SetBinContent(SeedChannelBin, SeedTickBin, 3);

   // for debugging only
   //h_tmp.Draw("colz");
   //C->Print("tmp.pdf");   

   int NeighbourX;
   int NeighbourY;

   std::vector<std::pair<int,int>> BinsAddedLastPass;
   std::vector<std::pair<int,int>> BinsAddedThisPass;

   BinsAddedLastPass.push_back(std::make_pair(SeedChannelBin, SeedTickBin));

   int NFillsThisPass = 1;

   while(NFillsThisPass > 0){

      NFillsThisPass = 0;
      BinsAddedThisPass.clear();

      // iterate over the bins added in the last Pass, check the bins that neighbour those
      for(size_t i_b = 0; i_b < BinsAddedLastPass.size(); i_b++){

         int CurrentBinX = BinsAddedLastPass.at(i_b).first;
         int CurrentBinY = BinsAddedLastPass.at(i_b).second;

         // look at each of the eight bins surrounding the current one

         NeighbourX = CurrentBinX + 1;
         NeighbourY = CurrentBinY;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX, NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX, NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX, NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX, NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }

         NeighbourX = CurrentBinX - 1;
         NeighbourY = CurrentBinY;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX, NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }

         NeighbourX = CurrentBinX;
         NeighbourY = CurrentBinY + 1;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX,NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }

         NeighbourX = CurrentBinX;
         NeighbourY = CurrentBinY - 1;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX,NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }


         NeighbourX = CurrentBinX + 1;
         NeighbourY = CurrentBinY + 1;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX,NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }

         NeighbourX = CurrentBinX - 1;
         NeighbourY = CurrentBinY + 1;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX,NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }

         NeighbourX = CurrentBinX + 1;
         NeighbourY = CurrentBinY - 1;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX,NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }

         NeighbourX = CurrentBinX - 1;
         NeighbourY = CurrentBinY - 1;                                
         // if bin is occupied, and not already part of the cluster, add it
         if(h_tmp.GetBinContent(NeighbourX,NeighbourY) > 0 && h_tmp.GetBinContent(NeighbourX,NeighbourY) < 2){ 
            h_tmp.SetBinContent(NeighbourX,NeighbourY,3);
            BinsAddedThisPass.push_back(std::make_pair(NeighbourX,NeighbourY));               
            C.bins_x.push_back(NeighbourX);   
            C.bins_y.push_back(NeighbourY);
            NFillsThisPass++;   
         }


      }//i_b

      BinsAddedLastPass = BinsAddedThisPass;

   } //while(NFillsThisPass > 0)


   clusters.push_back(C);

   return std::make_pair(Index, C.bins_x.size());

   // for debugging only
   // h_tmp.Draw("colz");
   // C->Print("tmp.pdf");   

}


std::pair<int,int> ClusterBuilder::FindNearestOccupiedBin(TH2D *Hist, int X, int Y)
{
   std::pair<int,int> NearestBin = { -1, -1};

   // Search neighbouring bins until you find an occupied one, within the ranges set by MaxX and MaxY

   for(int i_x = 1; i_x < maxSearchX + 1; i_x++){
      for(int i_y = 1; i_y < maxSearchY + 1; i_y++){
         if(Hist->GetBinContent(X + i_x, Y) > 0) { NearestBin = {X + i_x, Y}; return NearestBin; }
         if(Hist->GetBinContent(X - i_x, Y) > 0) { NearestBin = {X - i_x, Y}; return NearestBin; }
         if(Hist->GetBinContent(X , Y + i_y) > 0) { NearestBin = {X, Y + i_y}; return NearestBin; }
         if(Hist->GetBinContent(X, Y - i_y) > 0) { NearestBin = {X, Y - i_y}; return NearestBin; }
         if(Hist->GetBinContent(X + i_x, Y + i_y) > 0) { NearestBin = {X + i_x, Y + i_y}; return NearestBin; }
         if(Hist->GetBinContent(X - i_x, Y + i_y) > 0) { NearestBin = {X - i_x, Y + i_y}; return NearestBin; }
         if(Hist->GetBinContent(X + i_x, Y - i_y) > 0) { NearestBin = {X + i_x, Y - i_y}; return NearestBin; }
         if(Hist->GetBinContent(X - i_x, Y - i_y) > 0) { NearestBin = {X - i_x, Y - i_y}; return NearestBin; }
      }
   }

   // Return -1,-1 if no nearby bins are occupied
   return NearestBin;
}

std::vector<Cluster> ClusterBuilder::GetClusters()
{
   return clusters;
}

bool ClusterBuilder::SeedDeadWireCheck(std::vector<int> SeedsChannel, std::vector<int> SeedsTick, int Plane)
{
   //std::cout << "Checking if seeds are separated by dead channels" << std::endl;

   if(SeedsChannel.size() != SeedsTick.size()){
      std::cout << "Input Channel/Tick/Signal vectors of different sizes! Exiting" << std::endl;
      return true;
   }

   // if there is zero or 1 seed, this check is meaningless
   if(SeedsChannel.size() < 2) return false;

   // find lowest/highest Channel seeds
   int MinChn = SeedsChannel.at(0) - xOffset;
   int MaxChn = SeedsChannel.at(0) - xOffset;

   for(size_t i_s = 1; i_s < SeedsChannel.size(); i_s++){
      if(SeedsChannel.at(i_s) - xOffset < MinChn) MinChn = SeedsChannel.at(i_s) - xOffset;
      if(SeedsChannel.at(i_s) - xOffset > MaxChn) MaxChn = SeedsChannel.at(i_s) - xOffset;
   }

   int Ch = MinChn;

   while(Ch <= MaxChn){

      // Search the list of dead channels
      if(Plane == 0 && (std::find(DeadChannels_Plane0.begin(), DeadChannels_Plane0.end(), Ch) != DeadChannels_Plane0.end())){ return true; }
      if(Plane == 1 && (std::find(DeadChannels_Plane1.begin(), DeadChannels_Plane1.end(), Ch) != DeadChannels_Plane1.end())){ return true; }
      if(Plane == 2 && (std::find(DeadChannels_Plane2.begin(), DeadChannels_Plane2.end(), Ch) != DeadChannels_Plane2.end())){ return true; }

      Ch++;
   }

   return false;
}

void ClusterBuilder::DrawRaw(std::string RSE, int Pass)
{
   if(!drawEverything) return;

   if(Pass == -1) system(("mkdir -p Displays/" + displayDir + "/" + RSE + "/").c_str());
   else if(Pass == 0) system(("mkdir -p Displays/" + displayDir + "/Fail/" + RSE + "/").c_str());
   else if(Pass == 1) system(("mkdir -p Displays/" + displayDir + "/Pass/" + RSE + "/").c_str());

   hRaw->SetContour(100);
   hRaw->SetStats(0);
   hRaw->SetTitle(RSE.c_str());
   hRaw->Draw("colz");
   c->Print("Raw.pdf");

   if(Pass == -1) system(("mv Raw.pdf Displays/" + displayDir + "/" + RSE + "/").c_str());
   else if(Pass == 0) system(("mv Raw.pdf Displays/" + displayDir + "/Fail/" + RSE + "/").c_str());
   else if(Pass == 1) system(("mv Raw.pdf Displays/" + displayDir + "/Pass/" + RSE + "/").c_str());

   c->Clear();
}

void ClusterBuilder::DrawBinary(std::string RSE, int Pass)
{
   if(!drawEverything) return;

   std::string Cmd;
        
   if(Pass == -1) system(("mkdir -p Displays/" + displayDir + "/" + RSE + "/").c_str());
   else if(Pass == 0) system(("mkdir -p Displays/" + displayDir + "/Fail/" + RSE + "/").c_str());
   else if(Pass == 1) system(("mkdir -p Displays/" + displayDir + "/Pass/" + RSE + "/").c_str());

   Int_t colors[] = {0,4}; // #colors >= #levels - 1
   gStyle->SetPalette((sizeof(colors) / sizeof(Int_t)), colors);

   hBinary->SetStats(0);
   hBinary->SetTitle(RSE.c_str());
   hBinary->Draw("colz");
   c->Print("Binary.pdf");

   c->Clear();

   if(Pass == -1) system(("mv Binary.pdf Displays/" + displayDir + "/" + RSE + "/").c_str());
   else if(Pass == 0) system(("mv Binary.pdf Displays/" + displayDir + "/Fail/" + RSE + "/").c_str());
   else if(Pass == 1) system(("mv Binary.pdf Displays/" + displayDir + "/Pass/" + RSE + "/").c_str());

   gStyle->SetPalette();

}

void ClusterBuilder::DrawClustered(std::string RSE, int Plane, int Pass)
{
   if(!drawEverything) return;

   if(Pass == -1) system(("mkdir -p Displays/" + displayDir + "/" + RSE + "/").c_str());
   else if(Pass == 0) system(("mkdir -p Displays/" + displayDir + "/Fail/" + RSE + "/").c_str());
   else if(Pass == 1) system(("mkdir -p Displays/" + displayDir + "/Pass/" + RSE + "/").c_str());

   hClustered = (TH2D*)hBinary->Clone("hClustered");

   // Draw on dead wires - makes pdfs much bigger so use only if needed
   //if(Plane != -1) DeadWireFill(Plane);

   // Set color palette
   int NClusters = clusters.size();

   // hacky way to force color transitions to happen at the right bin height
   //hClustered->SetBinContent(1,1,-1.5);
   //hClustered->SetBinContent(1,2,NClusters+1.5);
     
   if(NClusters == 0){
      Int_t colors[] = {2,3};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }
   else if(NClusters == 1){
      Int_t colors[] = {2,3,4};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }
   else if(NClusters == 2){
      Int_t colors[] = {2,3,4,6};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }
   else if(NClusters == 3){
      Int_t colors[] = {2,3,4,6,7};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }
   else if(NClusters == 4){
      Int_t colors[] = {2,3,4,6,7,8};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }
   else if(NClusters == 5){
      Int_t colors[] = {2,3,4,6,7,40,46};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }
   else{
      Int_t colors[] = {2,3,4,6,7,40,46,12};
      gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   }

   // set bin heights equal to cluster IDs
   for(size_t i_cl = 0; i_cl < clusters.size(); i_cl++){
      for(size_t i_b = 0; i_b < clusters.at(i_cl).bins_x.size(); i_b++){
         hClustered->SetBinContent(clusters.at(i_cl).bins_x.at(i_b), clusters.at(i_cl).bins_y.at(i_b), i_cl + 2);
      }//i_b
   }//i_cl

   // Set X and Y ranges
   Focus();
 
   hClustered->SetStats(0);
   hClustered->SetTitle(RSE.c_str());
   hClustered->Draw("colz");
   c->Print("Clustered.pdf");

   if(Pass == -1) system(("mv Clustered.pdf Displays/" + displayDir + "/" + RSE + "/").c_str());
   else if(Pass == 0) system(("mv Clustered.pdf Displays/" + displayDir + "/Fail/" + RSE + "/").c_str());
   else if(Pass == 1) system(("mv Clustered.pdf Displays/" + displayDir + "/Pass/" + RSE + "/").c_str());

   c->Clear();

   DrawBinary(RSE, Pass);
   DrawRaw(RSE, Pass); 

   gStyle->SetPalette();
}

void ClusterBuilder::DeadWireFill(int Plane)
{
   if(Plane == 0){
      for(size_t i = 0; i < DeadChannels_Plane0.size(); i++){
         if(DeadChannels_Plane0.at(i) > hClustered->GetXaxis()->GetBinLowEdge(1) && DeadChannels_Plane0.at(i) < hClustered->GetXaxis()->GetBinLowEdge(hClustered->GetNbinsX())){
            for(int j = 0; j < hClustered->GetNbinsY(); j++){
               hClustered->Fill(DeadChannels_Plane0.at(i), hClustered->GetYaxis()->GetBinCenter(j), -1);
            }
         }
      }
   }

   if(Plane == 1){
      for(size_t i = 0; i < DeadChannels_Plane1.size(); i++){
         if(DeadChannels_Plane1.at(i) > hClustered->GetXaxis()->GetBinLowEdge(1) && DeadChannels_Plane1.at(i) < hClustered->GetXaxis()->GetBinLowEdge(hClustered->GetNbinsX())){
            for(int j = 0; j < hClustered->GetNbinsY(); j++){
               hClustered->Fill(DeadChannels_Plane1.at(i), hClustered->GetYaxis()->GetBinCenter(j), -1);
            }
         }
      }
   }

   if(Plane == 2){
      for(size_t i = 0; i < DeadChannels_Plane2.size(); i++){
         if(DeadChannels_Plane2.at(i) > hClustered->GetXaxis()->GetBinLowEdge(1) && DeadChannels_Plane2.at(i) < hClustered->GetXaxis()->GetBinLowEdge(hClustered->GetNbinsX())){
            for(int j = 0; j < hClustered->GetNbinsY(); j++){
               hClustered->Fill(DeadChannels_Plane2.at(i), hClustered->GetYaxis()->GetBinCenter(j), -1);
            }
         }
      }
   }
}

void ClusterBuilder::Focus()
{
   // Find range of X and Y spanned by clusters, only Draw this range + a bit of padding

   int MinX = 10000, MaxX = -10000, MinY = 10000, MaxY = -10000;

   for(size_t i_c = 0; i_c < clusters.size(); i_c++){
      Cluster C = clusters.at(i_c);

      for(size_t i_b = 0; i_b < C.bins_x.size(); i_b++){
         if(C.bins_x.at(i_b) > MaxX) MaxX = C.bins_x.at(i_b);
         if(C.bins_x.at(i_b) < MinX) MinX = C.bins_x.at(i_b);
         if(C.bins_y.at(i_b) > MaxY) MaxY = C.bins_y.at(i_b);
         if(C.bins_y.at(i_b) < MinY) MinY = C.bins_y.at(i_b);
      }
   }

   double  RangeX = MaxX - MinX;
   double  RangeY = MaxY - MinY;

   hClustered->GetXaxis()->SetRange(MinX - (RangeX * 0.2), MaxX + (RangeX * 0.2));
   hClustered->GetYaxis()->SetRange(MinY - (RangeY * 0.2), MaxY + (RangeY * 0.2));

   hRaw->GetXaxis()->SetRange(MinX - (RangeX * 0.2), MaxX + (RangeX * 0.2));
   hRaw->GetYaxis()->SetRange(MinY - (RangeY * 0.2), MaxY + (RangeY * 0.2));

   hBinary->GetXaxis()->SetRange(MinX - (RangeX * 0.2), MaxX + (RangeX * 0.2));
   hBinary->GetYaxis()->SetRange(MinY - (RangeY * 0.2), MaxY + (RangeY * 0.2));

   //hClustered->SetBinContent(MinX-RangeX*0.2+1,MinY-RangeY*0.2+1,);
   hClustered->SetBinContent(MinX - ((RangeX * 0.2) + 1), MinY - ((RangeY * 0.2) + 2), clusters.size() + 1.0);
}

void ClusterBuilder::SetDisplayDir(std::string Dir)
{
   displayDir = Dir;

   system(("mkdir -p Displays/" + displayDir).c_str());
   system(("mkdir -p Displays/" + displayDir + "/Pass/").c_str());
   system(("mkdir -p Displays/" + displayDir + "/Fail/").c_str());
}

#endif