#ifndef _PIDManager_cxx_
#define _PIDManager_cxx_

#include "ubana/HyperonProduction/Alg/PIDManager.h"

using namespace hyperon;

PIDManager::PIDManager()
{
   particleIdentificationCalculator.SetdEdxBinning(0, parametersProtonMuon.dedx_edges_pl_0);
   particleIdentificationCalculator.SetParBinning(0, parametersProtonMuon.parameters_edges_pl_0);
   particleIdentificationCalculator.SetLookupTables(0, parametersProtonMuon.dedx_pdf_pl_0);
   particleIdentificationCalculator.SetCorrParBinning(0, parametersCorrections.parameter_correction_edges_pl_0);
   particleIdentificationCalculator.SetCorrectionTables(0, parametersCorrections.correction_table_pl_0);
      
   particleIdentificationCalculator.SetdEdxBinning(1, parametersProtonMuon.dedx_edges_pl_1);
   particleIdentificationCalculator.SetParBinning(1, parametersProtonMuon.parameters_edges_pl_1);
   particleIdentificationCalculator.SetLookupTables(1, parametersProtonMuon.dedx_pdf_pl_1);
   particleIdentificationCalculator.SetCorrParBinning(1, parametersCorrections.parameter_correction_edges_pl_1);
   particleIdentificationCalculator.SetCorrectionTables(1, parametersCorrections.correction_table_pl_1);

   particleIdentificationCalculator.SetdEdxBinning(2, parametersProtonMuon.dedx_edges_pl_2);
   particleIdentificationCalculator.SetParBinning(2, parametersProtonMuon.parameters_edges_pl_2);
   particleIdentificationCalculator.SetLookupTables(2, parametersProtonMuon.dedx_pdf_pl_2);
   particleIdentificationCalculator.SetCorrParBinning(2, parametersCorrections.parameter_correction_edges_pl_2);
   particleIdentificationCalculator.SetCorrectionTables(2, parametersCorrections.correction_table_pl_2);
}

double PIDManager::GetMeandEdX(art::Ptr<anab::Calorimetry> Calo)
{
   double TotalE = 0;
   double TotalX = 0;

   if(Calo->XYZ().size() < 2) return -1;

   for(size_t i_point = 0; i_point < Calo->XYZ().size()-1; i_point++){

      anab::Point_t Pos = Calo->XYZ().at(i_point);
      anab::Point_t NextPos = Calo->XYZ().at(i_point+1);

      TVector3 D(Pos.X() - NextPos.X(), Pos.X() - NextPos.X(), Pos.X() - NextPos.X());

      TotalE += Calo->dEdx().at(i_point) * D.Mag();
      TotalX += D.Mag();
   }

   return TotalE / TotalX;
}

void PIDManager::SetThreePlaneMeandEdX(art::Ptr<recob::Track> Trk, std::vector<art::Ptr<anab::Calorimetry>> VectCalo, ParticleIdentifierStore &Store)
{
   double TotaldEdX = 0;
   double TotalWeight = 0;

   for(size_t i_plane = 0; i_plane < VectCalo.size(); i_plane++){

      int Plane = VectCalo.at(i_plane)->PlaneID().Plane;

      if(Plane != 0 && Plane != 1 && Plane != 2) continue;        

      double dEdX = GetMeandEdX(VectCalo.at(i_plane));

      if(dEdX < 0) continue;

      double PlaneWeight = GetPlaneWeight(Trk, Plane);

      if(Plane == 0){
         Store.Weight_Plane0 = PlaneWeight;       
         Store.MeandEdX_Plane0 = dEdX;
         Store.dEdX_Plane0 = VectCalo.at(i_plane)->dEdx();
         Store.ResidualRange_Plane0 = VectCalo.at(i_plane)->ResidualRange();
         Store.Pitch_Plane0 = VectCalo.at(i_plane)->TrkPitchVec();
         Store.dEdX_Corrected_Plane0 = particleIdentificationCalculator.CorrectManyHitsOnePlane(VectCalo.at(i_plane), *Trk, true, true);
      }
      if(Plane == 1){
         Store.Weight_Plane1 = PlaneWeight;       
         Store.MeandEdX_Plane1 = dEdX;
         Store.dEdX_Plane1 = VectCalo.at(i_plane)->dEdx();
         Store.ResidualRange_Plane1 = VectCalo.at(i_plane)->ResidualRange();
         Store.Pitch_Plane1 = VectCalo.at(i_plane)->TrkPitchVec();
         Store.dEdX_Corrected_Plane1 = particleIdentificationCalculator.CorrectManyHitsOnePlane(VectCalo.at(i_plane), *Trk, true, true);
      }
      if(Plane == 2){
         Store.Weight_Plane2 = PlaneWeight;       
         Store.MeandEdX_Plane2 = dEdX;
         Store.dEdX_Plane2 = VectCalo.at(i_plane)->dEdx();
         Store.ResidualRange_Plane2 = VectCalo.at(i_plane)->ResidualRange();
         Store.Pitch_Plane2 = VectCalo.at(i_plane)->TrkPitchVec();
         Store.dEdX_Corrected_Plane2 = particleIdentificationCalculator.CorrectManyHitsOnePlane(VectCalo.at(i_plane), *Trk, true, true);
      }

      TotaldEdX += dEdX * PlaneWeight;
      TotalWeight += PlaneWeight;
   }

   if(TotalWeight > 0) Store.MeandEdX_3Plane = TotaldEdX / TotalWeight;
}

void PIDManager::LLRPID(std::vector<art::Ptr<anab::Calorimetry>> VectCalo, ParticleIdentifierStore & Store)
{
   double LLRPID = 0;
   double LLRPIDScore = 0;

   for(auto const &Calo : VectCalo){

      auto const &Plane = Calo->PlaneID().Plane;
      auto const &dEdxValues = Calo->dEdx();
      auto const &RR = Calo->ResidualRange();
      auto const &Pitch = Calo->TrkPitchVec();

      std::vector<std::vector<float>> ParValues;
      ParValues.push_back(RR);
      ParValues.push_back(Pitch);

      std::vector<std::vector<float>> ParValuesPartial;
      std::vector<float> dEdxValuesPartial, RRPartial, PitchPartial;      
      
      if(Calo->dEdx().size() != Calo->ResidualRange().size() || Calo->ResidualRange().size() != Calo->TrkPitchVec().size())
         throw cet::exception("SubModuleReco") << "Track Calo point list size mismatch" << std::endl;

      for(size_t Plane = 0; Plane < Calo->dEdx().size(); Plane++){
         if(RR.at(Plane) > resRangeCutoff) continue;

         dEdxValuesPartial.push_back(Calo->dEdx().at(Plane));
         RRPartial.push_back(Calo->ResidualRange().at(Plane));
         PitchPartial.push_back(Calo->TrkPitchVec().at(Plane));        
      }

      ParValuesPartial.push_back(RRPartial);
      ParValuesPartial.push_back(PitchPartial);

      if(Calo->ResidualRange().size() == 0) continue;

      float CaloEnergy = 0;
      for(size_t i = 0; i < dEdxValues.size(); i++)
         CaloEnergy += dEdxValues[i] * Pitch[i];

      LLRPID += particleIdentificationCalculator.LLRManyHitsOnePlane(dEdxValues,ParValues,Plane);

      float CaloEnergyPartial = 0;
      for(size_t i = 0; i < dEdxValuesPartial.size(); i++)
         CaloEnergyPartial += dEdxValuesPartial[i] * PitchPartial[i];

   }

   LLRPIDScore = atan(LLRPID / 100.) * (2 / 3.14159266);

   Store.LLR = LLRPIDScore;
}

double PIDManager::GetBraggLikelihood(art::Ptr<recob::Track> Trk, std::vector<anab::sParticleIDAlgScores> VectorAlgScores, int PDG, anab::kTrackDir Dir)
{
   double ScorePlane0 = 0.0, ScorePlane1 = 0.0, ScorePlane2 = 0.0, ScoreWeighted = 0.0;

   for(size_t i_algscore = 0; i_algscore < VectorAlgScores.size(); i_algscore++){
      anab::sParticleIDAlgScores AlgScore = VectorAlgScores.at(i_algscore);

      if(AlgScore.fAssumedPdg == PDG && AlgScore.fAlgName == "BraggPeakLLH" && anab::kTrackDir(AlgScore.fTrackDir) == Dir){
         if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 0) ScorePlane0 = AlgScore.fValue;
         if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 1) ScorePlane1 = AlgScore.fValue;
         if(UBPID::uB_getSinglePlane(AlgScore.fPlaneMask) == 2) ScorePlane2 = AlgScore.fValue;
      }
   }

   ScoreWeighted = ScorePlane0 * GetPlaneWeight(Trk, 0) + ScorePlane1 * GetPlaneWeight(Trk, 1) + ScorePlane2 * GetPlaneWeight(Trk, 2);
   ScoreWeighted /= (GetPlaneWeight(Trk, 0) + GetPlaneWeight(Trk, 1) + GetPlaneWeight(Trk, 2)); 

   return ScoreWeighted;
}

void PIDManager::SetBraggScores(art::Ptr<recob::Track> Trk, std::vector<anab::sParticleIDAlgScores> VectorAlgScores, ParticleIdentifierStore & Store)
{
   Store.BraggWeighted_Pion = GetBraggLikelihood(Trk, VectorAlgScores, 211, anab::kForward); 
   Store.BraggWeighted_Muon = GetBraggLikelihood(Trk, VectorAlgScores, 13, anab::kForward);
   Store.BraggWeighted_Proton = GetBraggLikelihood(Trk, VectorAlgScores, 2212, anab::kForward);
   Store.BraggWeighted_Kaon = GetBraggLikelihood(Trk, VectorAlgScores, 321, anab::kForward);
   Store.BraggWeighted_Sigma = GetBraggLikelihood(Trk, VectorAlgScores, 3222, anab::kForward);
}

ParticleIdentifierStore  PIDManager::GetPIDScores(art::Ptr<recob::Track> Trk, std::vector<art::Ptr<anab::Calorimetry>> VectCalo, std::vector<anab::sParticleIDAlgScores> VectorAlgScores)
{
   ParticleIdentifierStore  Store;

   SetThreePlaneMeandEdX(Trk, VectCalo, Store);
   LLRPID(VectCalo, Store);
   SetBraggScores(Trk, VectorAlgScores, Store);

   return Store;
}

double PIDManager::GetPlaneWeight(art::Ptr<recob::Track> Trk, int Plane)
{
   // Returns a weight based on the angle between the Trk direction and wire Plane
   // Perpendicular tracks are more favourably weighted

   TVector3 Dir(Trk->End().x() - Trk->Start().x(), Trk->End().y() - Trk->Start().y(), Trk->End().z() - Trk->Start().z());

   TVector3 TrackVector(0, Dir.Y(), Dir.Z());
   TrackVector = TrackVector.Unit();

   TVector3 ZAxis(0, 0, 1);

   double CosThetaYZ = TrackVector.Dot(ZAxis);
   double ThetaYZ = TMath::ACos(CosThetaYZ);

   if ((Dir.Y() < 0) && (ThetaYZ > 0)) ThetaYZ *= -1;

   double thetaToWires = 0;
   if (Plane == 0) thetaToWires = std::min(std::abs(wireAnglePlane0 - ThetaYZ), std::abs((-1*(6.28-wireAnglePlane0) - ThetaYZ)));
   if (Plane == 1) thetaToWires = std::min(std::abs(wireAnglePlane1 - ThetaYZ), std::abs((-1*(6.28-wireAnglePlane1) - ThetaYZ)));
   if (Plane == 2) thetaToWires = std::min(std::abs(wireAnglePlane2 - ThetaYZ), std::abs((-1*(6.28-wireAnglePlane2) - ThetaYZ)));

   double AnglePlaneWeihght = sin(thetaToWires) * sin(thetaToWires);
   if (AnglePlaneWeihght < tophatThresh) AnglePlaneWeihght = 0;
   if (AnglePlaneWeihght != 0) AnglePlaneWeihght = 1;

   return AnglePlaneWeihght;
}

#endif