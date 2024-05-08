#ifndef _RecoParticle_h_
#define _RecoParticle_h_

#include <iostream>

#include "TLorentzVector.h"
#include "TVector3.h"

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class RecoParticle : public TObject{
#else
class RecoParticle {
#endif

public:

RecoParticle(){}
~RecoParticle(){}

int Index;
bool InNuSlice = false;

// General reco info
int PDG; // Pandora PDG code (11 or 13)
int Parentage; // 1 - neutrino daughter, 2 - neutrino granddaughter, 3 - other
int ParentIndex=-1; // -1 - neutrino candidate or no parent
double TrackShowerScore;
double NuScore;
double X,Y,Z;
double Displacement; // Distance from PV

// Track variables
double TrackLength = 0;
double TrackDirectionX = 0, TrackDirectionY = 0, TrackDirectionZ = 0;
double TrackStartX = 0, TrackStartY = 0, TrackStartZ = 0;
double TrackEndX = 0, TrackEndY = 0, TrackEndZ = 0;
double TrackPID; 

double MeandEdX_Plane0, MeandEdX_Plane1, MeandEdX_Plane2, MeandEdX_ThreePlane;  
 
double Track_LLR_PID; 
double Track_LLR_PID_Kaon; 
double Track_LLR_PID_Kaon_Partial; 
double Track_Bragg_Pion;
double Track_Bragg_Muon;
double Track_Bragg_Proton; 
double Track_Bragg_Kaon;
double Track_Bragg_Sigma;
 
double ProtonMomentum, MuonMomentum, KaonMomentum;
double TrackWiggliness;

bool HasTruth; 
int MCTruthIndex= - 1;

int TrackTruePDG;
double TrackTrueE,TrackTruePx,TrackTruePy,TrackTruePz;
double TrackTrueEndE,TrackTrueEndPx,TrackTrueEndPy,TrackTrueEndPz;
double TrackTrueModMomentum;
double TrackTrueEndModMomentum;
double TrackTrueKE;
double TrackTrueEndKE;
double TrackTrueLength;
int TrackTrueOrigin;
double TrackTruthPurity;

inline void SetVertex(TVector3 V);
inline void SetTrackPositions(TVector3 Start,TVector3 End);
inline void Print();

#ifdef __MAKE_ROOT_DICT__
ClassDef(RecoParticle,1);
#endif

};

inline void RecoParticle::SetVertex(TVector3 V)
{
   X = V.X();
   Y = V.Y();
   Z = V.Z();
}

inline void RecoParticle::SetTrackPositions(TVector3 Start,TVector3 End)
{
   TrackStartX = Start.X();
   TrackStartY = Start.Y();
   TrackStartZ = Start.Z();

   TrackEndX = End.X();
   TrackEndY = End.Y();
   TrackEndZ = End.Z();
}

inline void RecoParticle::Print()
{
   std::cout << "Reco Info:" << std::endl;
   std::cout << "PDG Code: " << PDG << "  Track/Shower score: " << TrackShowerScore << std::endl;
   std::cout << "Track length: " << TrackLength << "  PID score: " << TrackPID <<  std::endl;
   std::cout << "Truth Info:" << std::endl;
   std::cout << "PDG: " << TrackTruePDG << "  Origin: " << TrackTrueOrigin << std::endl;
   std::cout << "Length: " << TrackTrueLength << "  KE: " << TrackTrueKE << std::endl;
}

#endif