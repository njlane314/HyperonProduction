#ifndef _BDTHandle_cxx_
#define _BDTHandle_cxx_

#include "ubana/HyperonProduction/Alg/BDTHandle.h"

using namespace hyperon;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

BDTHandle::BDTHandle(std::string weightsdir):
   WeightsDir(weightsdir)
{  
   TMVA::Tools::Instance();

   TString methodName = TString("BDT method");

   char const* tmp = getenv("UBOONEDATA_DIR");
   if(tmp == NULL) throw cet::exception("BDTHandle") << "UBOONEDATA_DIR not set!" << std::endl;
   TString weightfile = std::string(tmp) + "/" + WeightsDir + "/TMVAClassification_BDT.weights.xml";

   std::cout << std::endl << "BDTHandle: Opening weight file " << weightfile << std::endl << std::endl;  
   
   reader = new TMVA::Reader( "!Color:!Silent" );
   reader->AddVariable("KaonTrackPID",&v_KaonTrackPID);
   reader->AddVariable("KaonTrackBraggPID",&v_KaonTrackBraggPID);
   reader->BookMVA(methodName,weightfile);   
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BDTHandle::GetScore(double kaontrackpID,double kaontrackbraggpid){

   v_KaonTrackPID = kaontrackpID;
   v_KaonTrackBraggPID = kaontrackbraggpid;
   
   return reader->EvaluateMVA("BDT method");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
