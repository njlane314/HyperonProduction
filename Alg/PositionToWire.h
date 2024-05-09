#ifndef _Position_To_Wire_h_
#define _Position_To_Wire_h_

#include "TVector3.h"

// Convert 3D positions to channels/ticks
// See docdb 25505-v2 for explanation

const double A_w = 3.33328;
const double C_U = 338.140;
const double C_V = 2732.53;
const double C_Y = 4799.19;
const double A_t = 18.2148;
const double C_t = 818.351;

// cos(60) and sin(60)
const double cos60 = 0.5;
const double sin60 = sqrt(3)/2.0;

int U_wire(TVector3 Pos) { return A_w * (-sin60 * Pos.Y() + cos60 * Pos.Z()) + C_U; }
int V_wire(TVector3 Pos) { return A_w * (sin60 * Pos.Y() + cos60 * Pos.Z()) + C_V; }
int Y_wire(TVector3 Pos) { return A_w * Pos.Z() + C_Y; }
int tick(TVector3 Pos) { return A_t * Pos.X() + C_t; }

double dUdt(TVector3 Dir){ return (A_w / A_t) * (-sin60 * (Dir.Y() / Dir.X()) + cos60 * (Dir.Z()/Dir.X())); }
double dVdt(TVector3 Dir){ return (A_w / A_t) * (sin60 * (Dir.Y() / Dir.X()) + cos60 * (Dir.Z()/Dir.X())); }
double dYdt(TVector3 Dir){ return (A_w / A_t) * (Dir.Z() / Dir.X()); }

double AngleU(TVector3 Dir)
{
   bool Invert = Dir.X() < 0;

   double Angle = (180 / 3.141) * atan(dUdt(Dir));
   Angle -= 2 * (Angle - 45.0);

   if(Invert && Angle < 0) Angle += 180;
   if(Invert && Angle > 0) Angle -= 180;

   return Angle;
}

double AngleV(TVector3 Dir)
{
   bool Invert = Dir.X() < 0;

   double Angle = (180 / 3.141) * atan(dVdt(Dir));
   Angle -= 2 * (Angle - 45.0);

   if(Invert && Angle < 0) Angle += 180;
   if(Invert && Angle > 0) Angle -= 180;

   return Angle;
}

double AngleY(TVector3 Dir)
{
   bool Invert = Dir.X() < 0;

   double Angle = (180 / 3.141) * atan(dYdt(Dir));
   Angle -= 2 * (Angle - 45.0);

   if(Invert && Angle < 0) Angle += 180;
   if(Invert && Angle > 0) Angle -= 180;
   
   return Angle;
}

#endif