#ifndef _LLRPIDHelper_cxx_
#define _LLRPIDHelper_cxx_

#include "ubana/HyperonProduction/Alg/LLRPIDHelper.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

LLRPIDHelper::LLRPIDHelper(){

   llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
   llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
   llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

   llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
   llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
   llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

   llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
   llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
   llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

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

LLRPID_Result LLRPIDHelper::GetScores(std::vector<art::Ptr<anab::Calorimetry>> calo_v){

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

   LLRPID_Result theResult;
   theResult.Score = this_llr_pid_score;
   theResult.Score_Kaon = this_llr_pid_score_kaon;
   theResult.Score_Kaon_Partial = this_llr_pid_score_kaon_partial;

   return theResult;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
