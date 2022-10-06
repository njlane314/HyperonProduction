#ifndef _PIDManager_cxx_
#define _PIDManager_cxx_

#include "ubana/HyperonProduction/Alg/PIDManager.h"

using namespace hyperon;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

PIDManager::PIDManager(){

   llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
   llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
   llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);
   llr_pid_calculator.set_corr_par_binning(0,correction_parameters.parameter_correction_edges_pl_0);
   llr_pid_calculator.set_correction_tables(0,correction_parameters.correction_table_pl_0);
      
   llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
   llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
   llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);
   llr_pid_calculator.set_corr_par_binning(1,correction_parameters.parameter_correction_edges_pl_1);
   llr_pid_calculator.set_correction_tables(1,correction_parameters.correction_table_pl_1);

   llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
   llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
   llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);
   llr_pid_calculator.set_corr_par_binning(2,correction_parameters.parameter_correction_edges_pl_2);
   llr_pid_calculator.set_correction_tables(2,correction_parameters.correction_table_pl_2);

   llr_pid_calculator_kaon.set_dedx_binning(0, kaonproton_parameters.dedx_edges_pl_0);
   llr_pid_calculator_kaon.set_par_binning(0, kaonproton_parameters.parameters_edges_pl_0);
   llr_pid_calculator_kaon.set_lookup_tables(0, kaonproton_parameters.dedx_pdf_pl_0);

   llr_pid_calculator_kaon.set_dedx_binning(1, kaonproton_parameters.dedx_edges_pl_1);
   llr_pid_calculator_kaon.set_par_binning(1, kaonproton_parameters.parameters_edges_pl_1);
   llr_pid_calculator_kaon.set_lookup_tables(1, kaonproton_parameters.dedx_pdf_pl_1);

   llr_pid_calculator_kaon.set_dedx_binning(2, kaonproton_parameters.dedx_edges_pl_2);
   llr_pid_calculator_kaon.set_par_binning(2, kaonproton_parameters.parameters_edges_pl_2);
   llr_pid_calculator_kaon.set_lookup_tables(2, kaonproton_parameters.dedx_pdf_pl_2);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::GetMeandEdX(art::Ptr<anab::Calorimetry> calo){

   double totalE=0;
   double totalX=0;

   // Sometimes this vector is empty, causes crash below, skip plane if it is
   if(calo->XYZ().size() < 2) return -1;

   for(size_t i_point = 0;i_point < calo->XYZ().size()-1;i_point++){

      anab::Point_t thisPos = calo->XYZ().at(i_point);
      anab::Point_t nextPos = calo->XYZ().at(i_point+1);

      // Step vector
      TVector3 D(thisPos.X()-nextPos.X(),thisPos.X()-nextPos.X(),thisPos.X()-nextPos.X());

      totalX += D.Mag();
      totalE += calo->dEdx().at(i_point)*D.Mag();
   }

   return totalE/totalX;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PIDManager::ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store){

   double TotaldEdX=0;
   double TotalWeight=0;

   for(size_t i_pl=0;i_pl<calo_v.size();i_pl++){

      int plane = calo_v.at(i_pl)->PlaneID().Plane;

      if(plane != 0 && plane != 1 && plane != 2) continue;        

      double dEdX = GetMeandEdX(calo_v.at(i_pl));

      // Catch default fills
      if(dEdX < 0) continue;

      double thisPlaneWeight = PlaneWeight(track,plane);

      if(plane == 0){
         store.Weight_Plane0 = thisPlaneWeight;       
         store.MeandEdX_Plane0 = dEdX;
         store.dEdX_Plane0 = calo_v.at(i_pl)->dEdx();
         store.ResidualRange_Plane0 = calo_v.at(i_pl)->ResidualRange();
         store.Pitch_Plane0 = calo_v.at(i_pl)->TrkPitchVec();
         store.dEdX_Corrected_Plane0 = llr_pid_calculator.correct_many_hits_one_plane(calo_v.at(i_pl),*track,true,true);
      }
      if(plane == 1){
         store.Weight_Plane1 = thisPlaneWeight;       
         store.MeandEdX_Plane1 = dEdX;
         store.dEdX_Plane1 = calo_v.at(i_pl)->dEdx();
         store.ResidualRange_Plane1 = calo_v.at(i_pl)->ResidualRange();
         store.Pitch_Plane1 = calo_v.at(i_pl)->TrkPitchVec();
         store.dEdX_Corrected_Plane1 = llr_pid_calculator.correct_many_hits_one_plane(calo_v.at(i_pl),*track,true,true);
      }
      if(plane == 2){
         store.Weight_Plane2 = thisPlaneWeight;       
         store.MeandEdX_Plane2 = dEdX;
         store.dEdX_Plane2 = calo_v.at(i_pl)->dEdx();
         store.ResidualRange_Plane2 = calo_v.at(i_pl)->ResidualRange();
         store.Pitch_Plane2 = calo_v.at(i_pl)->TrkPitchVec();
         store.dEdX_Corrected_Plane2 = llr_pid_calculator.correct_many_hits_one_plane(calo_v.at(i_pl),*track,true,true);
      }

      TotaldEdX += dEdX*thisPlaneWeight;
      TotalWeight += thisPlaneWeight;
   }

   if(TotalWeight > 0) store.MeandEdX_3Plane = TotaldEdX/TotalWeight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PIDManager::LLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,PIDStore& store){

   double this_llr_pid=0;
   double this_llr_pid_score=0;
   double this_llr_pid_kaon=0;
   double this_llr_pid_score_kaon=0;

   double this_llr_pid_kaon_partial = 0;
   double this_llr_pid_score_kaon_partial = 0;

   for(auto const &calo : calo_v){

      auto const &plane = calo->PlaneID().Plane;
      auto const &dedx_values = calo->dEdx();
      auto const &rr = calo->ResidualRange();
      auto const &pitch = calo->TrkPitchVec();
      std::vector<std::vector<float>> par_values;
      par_values.push_back(rr);
      par_values.push_back(pitch);

      // Get parital length PIDs
      std::vector<std::vector<float>> par_values_partial;
      std::vector<float> dedx_values_partial,rr_partial,pitch_partial;      
      if(calo->dEdx().size() != calo->ResidualRange().size() || calo->ResidualRange().size() != calo->TrkPitchVec().size())
         throw cet::exception("SubModuleReco") << "Track calo point list size mismatch" << std::endl;
      for(size_t i_p=0;i_p<calo->dEdx().size();i_p++){
         if(rr.at(i_p) > ResRangeCutoff) continue;
         dedx_values_partial.push_back(calo->dEdx().at(i_p));
         rr_partial.push_back(calo->ResidualRange().at(i_p));
         pitch_partial.push_back(calo->TrkPitchVec().at(i_p));        
      }
      par_values_partial.push_back(rr_partial);
      par_values_partial.push_back(pitch_partial);

      if(calo->ResidualRange().size() == 0) continue;

      float calo_energy = 0;
      for(size_t i=0;i<dedx_values.size();i++)
         calo_energy += dedx_values[i] * pitch[i];

      float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values,par_values,plane);
      float llr_pid_kaon = llr_pid_calculator_kaon.LLR_many_hits_one_plane(dedx_values,par_values,plane);
      this_llr_pid += llr_pid;
      this_llr_pid_kaon += llr_pid_kaon;

      // Partial length calculation
      float calo_energy_partial = 0;
      for(size_t i=0;i<dedx_values_partial.size();i++)
         calo_energy_partial += dedx_values_partial[i] * pitch_partial[i];

      float llr_pid_kaon_partial = llr_pid_calculator_kaon.LLR_many_hits_one_plane(dedx_values_partial,par_values_partial,plane);
      this_llr_pid_kaon_partial += llr_pid_kaon_partial;     
   }

   this_llr_pid_score = atan(this_llr_pid/100.)*2/3.14159266;
   this_llr_pid_score_kaon = atan(this_llr_pid_kaon/100.)*2/3.14159266;
   this_llr_pid_score_kaon_partial = atan(this_llr_pid_kaon_partial/100.)*2/3.14159266;

   store.LLR = this_llr_pid_score;
   store.LLR_Kaon = this_llr_pid_score_kaon;
   store.LLR_Kaon_Partial = this_llr_pid_score_kaon_partial;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PIDManager::BraggPID(art::Ptr<recob::Track> track,std::vector<anab::sParticleIDAlgScores> algscores_v,PIDStore& store){

   for(size_t i_algscore=0;i_algscore<algscores_v.size();i_algscore++){
      anab::sParticleIDAlgScores algscore = algscores_v.at(i_algscore);
      if(algscore.fAssumedPdg == 321 && algscore.fAlgName=="BraggPeakLLH" && anab::kTrackDir(algscore.fTrackDir) == anab::kForward){
         if(UBPID::uB_getSinglePlane(algscore.fPlaneMask) == 0) store.Bragg_Kaon_Plane0 = algscore.fValue;
         if(UBPID::uB_getSinglePlane(algscore.fPlaneMask) == 1) store.Bragg_Kaon_Plane1 = algscore.fValue;
         if(UBPID::uB_getSinglePlane(algscore.fPlaneMask) == 2) store.Bragg_Kaon_Plane2 = algscore.fValue;
      }
   }

   store.Bragg_Kaon_3Plane = store.Bragg_Kaon_Plane0*PlaneWeight(track,0) + store.Bragg_Kaon_Plane1*PlaneWeight(track,1) + store.Bragg_Kaon_Plane2*PlaneWeight(track,2);
   store.Bragg_Kaon_3Plane /= (PlaneWeight(track,0) + PlaneWeight(track,1) + PlaneWeight(track,2));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

PIDStore PIDManager::GetPIDs(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::vector<anab::sParticleIDAlgScores> algscores_v){

   PIDStore theStore;
   ThreePlaneMeandEdX(track,calo_v,theStore);
   LLRPID(calo_v,theStore);
   BraggPID(track,algscores_v,theStore);

   return theStore;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::PlaneWeight(TVector3 dir,int i_pl){

   TVector3 trackvec(0, dir.Y(), dir.Z());
   trackvec = trackvec.Unit();
   TVector3 zaxis(0, 0, 1);
   double costhetayz = trackvec.Dot(zaxis);
   double thetayz = TMath::ACos(costhetayz);
   if ((dir.Y() < 0) && (thetayz > 0)) thetayz *= -1;

   double theta_towires = 0;
   if (i_pl == 0) theta_towires = std::min(std::abs(plane0_wireangle - thetayz), std::abs((-1*(6.28-plane0_wireangle) - thetayz)));
   if (i_pl == 1) theta_towires = std::min(std::abs(plane1_wireangle - thetayz), std::abs((-1*(6.28-plane1_wireangle) - thetayz)));
   if (i_pl == 2) theta_towires = std::min(std::abs(plane2_wireangle - thetayz), std::abs((-1*(6.28-plane2_wireangle) - thetayz)));

   double angle_planeweight = sin(theta_towires)*sin(theta_towires);
   if (angle_planeweight < TophatThresh) angle_planeweight = 0;
   if (angle_planeweight != 0) angle_planeweight = 1;

   return angle_planeweight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDManager::PlaneWeight(art::Ptr<recob::Track> track,int i_pl){

   TVector3 trackdir(track->End().x()-track->Start().x(),track->End().y()-track->Start().y(),track->End().z()-track->Start().z());
   return PlaneWeight(trackdir,i_pl);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
