#ifndef _LLRPID_h_
#define _LLRPID_h_

#include "Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace searchingfornues
{
   class LLRPID
   {
      private:

         size_t dEdxNumBins[3];
         std::vector<float> dEdxBinEdges[3];

         std::vector<size_t> parametersNumBins[3];
         std::vector<std::vector<float>> parametersBinEdges[3];

         std::vector<size_t> corrParametersNumBins[3];
         std::vector<std::vector<float>> corrParametersBinEdges[3];

         std::vector<float> lookupTables[3];
         std::vector<float> correctionTables[3];

      public:

         LLRPID(){}

         void SetdEdxBinning(size_t Plane, std::vector<float> BinEdges)
         {
            dEdxBinEdges[Plane] = BinEdges;
            dEdxNumBins[Plane] = (BinEdges.size() - 1);
         }

         void SetParBinning(size_t Plane, std::vector<std::vector<float>> BinEdges)
         {
            parametersBinEdges[Plane] = BinEdges;
            for (size_t i = 0; i < BinEdges.size(); i++){
               parametersNumBins[Plane].push_back(parametersBinEdges[Plane][i].size() - 1);
            }
         }
         
         void SetCorrParBinning(size_t Plane, std::vector<std::vector<float>> BinEdges)
         {
            corrParametersBinEdges[Plane] = BinEdges;
            parametersBinEdges[Plane] = BinEdges;
            for (size_t i = 0; i < BinEdges.size(); i++){
               corrParametersNumBins[Plane].push_back(corrParametersBinEdges[Plane][i].size() - 1);
            }
         }

         void SetLookupTables(size_t Plane, std::vector<float> Tables)
         {
            lookupTables[Plane] = Tables;
         }

         void SetCorrectionTables(size_t Plane, std::vector<float> Tables)
         {
            correctionTables[Plane] = Tables;
         }

         size_t Digitize(float Value, std::vector<float> BinEdges)
         {
            if (Value <= BinEdges[0])
               return 0;

            for(size_t i = 0; i < BinEdges.size(); i++){
               if (Value >= BinEdges[i])
                  continue;
               else
                  return i-1;
            }

            return BinEdges.size()-1;
         }

         size_t FindLookupTable(float dEdxValue, std::vector<float> ParValue, size_t Plane)
         {
            std::vector<size_t> ParametersBins;
            for(size_t i = 0; i < ParValue.size(); i++){
               size_t AuxIndex = Digitize(ParValue[i], parametersBinEdges[Plane][i]);
               ParametersBins.push_back(AuxIndex);
            }

            size_t LookupRow = 0, AccumulatorParBins = 1;
            for(size_t i = ParametersBins.size(); i-- > 0; ){
               LookupRow += (AccumulatorParBins * ParametersBins[i]);
               AccumulatorParBins *= parametersNumBins[Plane][i];
            }

            size_t LookupRowIndex;
            LookupRowIndex = LookupRow * dEdxNumBins[Plane];

            size_t LookupIndex = LookupRowIndex;
            LookupIndex += Digitize(dEdxValue, dEdxBinEdges[Plane]);

            return LookupIndex;
         }

         size_t FindLookupCorrParameterIndex(std::vector<float> CorrParameterValue, size_t Plane)
         {
            std::vector<size_t> CorrParametersBins;

            for(size_t i = 0; i < CorrParameterValue.size(); i++){
               size_t AuxIndex = Digitize(CorrParameterValue[i], corrParametersBinEdges[Plane][i]);
               CorrParametersBins.push_back(AuxIndex);
            }
            
            size_t LookupIndex = 0, AccumulatorParBins = 1;
            for(size_t i = CorrParametersBins.size(); i-- > 0; ){
               LookupIndex += (AccumulatorParBins * CorrParametersBins[i]);
               AccumulatorParBins *= corrParametersNumBins[Plane][i];
            }

            return LookupIndex;
         }

         float LLROneHitOnePlane(float dEdxValue, std::vector<float> ParValue, size_t Plane)
         {
            size_t Index = FindLookupTable(dEdxValue, ParValue, Plane);

            return lookupTables[Plane][Index];
         }

         float CorrectionHitOnePlane(std::vector<float> CorrParameterValue, size_t Plane)
         {
            size_t Index = FindLookupCorrParameterIndex(CorrParameterValue, Plane);

            return correctionTables[Plane][Index];
         }

         float LLRManyHitsOnePlane(std::vector<float> dEdxValues, std::vector<std::vector<float>> ParValues, size_t Plane)
         {
            float LLOut = 0;

            for(size_t i = 0; i < dEdxValues.size(); i++){
               std::vector<float> AuxPar;

               for(std::vector<float> ParValue: ParValues){
                  AuxPar.push_back(ParValue[i]);
               }

               LLOut += LLROneHitOnePlane(dEdxValues[i], AuxPar, Plane);
            }

            return LLOut;
         }

         std::vector<float> CorrectManyHitsOnePlane(std::vector<float> dEdxValues, std::vector<std::vector<float>> CorrParValues, std::vector<bool> IsToCorrect, size_t Plane)
         {
            std::vector<float> dEdxValuesCorrected;

            for(size_t i = 0; i < dEdxValues.size(); i++){
               float AuxdEdX = dEdxValues[i];

               if (IsToCorrect[i]){
                  AuxdEdX *= CorrectionHitOnePlane(CorrParValues[i], Plane);
               }

               dEdxValuesCorrected.push_back(AuxdEdX);
            }

            return dEdxValuesCorrected;
         }

         std::vector<float> CorrectManyHitsOnePlane(const art::Ptr<anab::Calorimetry> TrkCalo, const recob::Track Trk, const bool RecalibrateHits, const bool LocaldEdx)
         {
            int Plane = TrkCalo->PlaneID().Plane;
            std::vector<float> dQdxValues, dQdxValuesCorrected;

            if (!LocaldEdx)
               dQdxValues = TrkCalo->dQdx();
            else
               dQdxValues = TrkCalo->dEdx();
            if (!RecalibrateHits)
               dQdxValuesCorrected = dQdxValues;
            else
            {
               auto const &Pitch = TrkCalo->TrkPitchVec();
               auto const& VectXYZ = TrkCalo->XYZ();

               std::vector<std::vector<float>> CorrParValues;
               for (size_t i = 0; i < VectXYZ.size(); i++){
                  auto XYZ = VectXYZ[i];
                  float Dir[3];
                  searchingfornues::TrkDirectionAtXYZ(Trk, XYZ.X(), XYZ.Y(), XYZ.Z(), Dir);
                  std::vector<float> temp = searchingfornues::polarAngles(Dir[0], Dir[1], Dir[2], 2, Plane);
                  temp[0] = Pitch[i];
                  CorrParValues.push_back(temp);
               }
                     
               std::vector<bool> IsHitMC;
               const std::vector<size_t> &TpIndices = TrkCalo->TpIndices();
               for (size_t i = 0; i < TpIndices.size(); i++){
                  IsHitMC.push_back(true);
               }
          
               dQdxValuesCorrected = CorrectManyHitsOnePlane(dQdxValues, CorrParValues, IsHitMC, Plane);
            }

            return dQdxValuesCorrected;
         }
   };
}

#endif