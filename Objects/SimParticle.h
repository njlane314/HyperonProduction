#ifndef _SimParticle_h_
#define _SimParticle_h_

#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class SimParticle : public TObject{
#else
class SimParticle {
#endif
   public:

   SimParticle() {}
   ~SimParticle() {}

      int MCTruthIndex = -1;

      int PDG = 0;
      double E = 0, Px = 0, Py = 0, Pz = 0;
      double ModMomentum = 0;
      double EndE = 0, EndPx = 0, EndPy = 0, EndPz = 0;
      double EndModMomentum = 0;
      double KE = 0;
      double EndKE = 0;
      double StartX = 0, StartY = 0, StartZ = 0;
      double EndX = 0, EndY = 0, EndZ = 0;
      double Travel = 0; 
      double Theta = 0;
      double Phi = 0;
      int Origin = 0;

      double ScatterLength = 0;

      inline void SetKinematics(TLorentzVector P, TLorentzVector EndP, double Mass);
      inline void SetPositions(TLorentzVector Start, TLorentzVector End);
      inline void Print();

   #ifdef __MAKE_ROOT_DICT__
   ClassDef(SimParticle, 1);
   #endif
};

inline void SimParticle::SetKinematics(TLorentzVector P, TLorentzVector EndP, double Mass)
{
   E = P.E();
   Px = P.X();
   Py = P.Y();
   Pz = P.Z();
   ModMomentum = sqrt(Px*Px + Py*Py + Pz*Pz);
   KE = E - Mass;

   Theta = (180/3.1416)*TMath::ACos(Pz / ModMomentum);
   Phi = (180/3.1416)*TMath::ASin(Py / sqrt(Px*Px+Py*Py));

   EndE = EndP.E();
   EndPx = EndP.X();
   EndPy = EndP.Y();
   EndPz = EndP.Z();
   EndModMomentum = sqrt(EndPx*EndPx + EndPy*EndPy + EndPz*EndPz);
   EndKE = EndE - Mass;
}

inline void SimParticle::SetPositions(TLorentzVector Start, TLorentzVector End)
{
   StartX = Start.X();
   StartY = Start.Y();
   StartZ = Start.Z();

   EndX = End.X();
   EndY = End.Y();
   EndZ = End.Z();

   Travel = sqrt( (StartX - EndX)*(StartX - EndX) + (StartY - EndY)*(StartY - EndY) + (StartZ - EndZ)*(StartZ - EndZ) );
}

inline void SimParticle::Print()
{
   std::cout << "PDG: " << PDG << "  Origin: " << Origin << std::endl;
   std::cout << "Length: " << Travel << "  KE: " << KE << std::endl;
}

#endif