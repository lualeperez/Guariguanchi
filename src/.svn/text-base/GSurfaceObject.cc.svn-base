#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2D.h>
#include <TText.h>
#include <TSystem.h>
#include <TLine.h>
#include <TFile.h>
#include <TEllipse.h>
#include <TVector2.h>
#include <TVector3.h>
#include "include/GSurfaceObject.h"

//C++, C
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <cstring>

using namespace std;

//====================================================================
GSurfaceObject::GSurfaceObject(TString   aName,
			       TVector3  aPosition,
			       TMatrixD  aRot,
			       bool      aIsTrackingLayer,
			       GGlobalTools* aglobal)
{
  
  Name            = aName;
  Position        = aPosition;
  IsTrackingLayer = aIsTrackingLayer;
  
  Rot.ResizeTo(3,3);
  InvRot.ResizeTo(3,3);
  Rot    = aRot;
  InvRot = Rot;
  InvRot.Invert();
  
  global = aglobal;
  
  UVector = TVector3(1.0,0.0,0.0);
  VVector = TVector3(0.0,1.0,0.0);
  WVector = TVector3(0.0,0.0,1.0);
  global->RotateVector(InvRot,UVector);
  global->RotateVector(InvRot,VVector);
  global->RotateVector(InvRot,WVector);
  
  LadderType = TString("");
  
  verbose = false;
    
}
//====================================================================
GSurfaceObject::GSurfaceObject(const GSurfaceObject& other,TString aName)
{
  
  Name     = aName;
  Position = other.Position;
  
  Rot.ResizeTo(3,3);
  InvRot.ResizeTo(3,3);
  Rot    = other.Rot;
  InvRot = other.Rot;
  InvRot.Invert();
  
  global = other.global;
  
  UVector = TVector3(1.0,0.0,0.0);
  VVector = TVector3(0.0,1.0,0.0);
  WVector = TVector3(0.0,0.0,1.0);
  global->RotateVector(InvRot,UVector);
  global->RotateVector(InvRot,VVector);
  global->RotateVector(InvRot,WVector);
  
  LadderType = other.LadderType;
  
  verbose = other.verbose;
  
}
//====================================================================
GSurfaceObject::~GSurfaceObject() 
{
  
}
//====================================================================
GSurfaceObject* GSurfaceObject::clone(TString aName) const
{
 
  return new GSurfaceObject(*this,aName);
  
}
//====================================================================
void  GSurfaceObject::SetRotMatrix(TMatrixD R)
{
  
  //Set internal rotation matrix
  
  Rot    = R;
  InvRot = R;
  InvRot.Invert(); 
  
  return;
  
}
//====================================================================
void  GSurfaceObject::CheckIndex(int idx)
{
  
  //Check index limits
  
  if(idx < 0 || idx > 2) {
    cout << endl;
    cout << "GSurfaceObject:: index is outside range (0,2). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
TVector3  GSurfaceObject::GetDerUVWFromXYZ_Numerical(TVector3 PositionXYZ,int idx)
{

  //Return the derivative of the local (U,V,W) coordinates w.r.t to global XYZ coordinates
  
  CheckIndex(idx);
  
  double h = 1.0e-2*global->GetDistanceUnit("mm");
  
  TVector3 Prev_der(-999.9,-999.9,-999.9);
  TVector3 Current_der(-999.9,-999.9,-999.9);
  TVector3 Current_der1,Current_der2;
  
  int counter = 0;
  TVector3 diff(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;
    
    TVector3 PositionXYZ_p  = PositionXYZ;
    TVector3 PositionXYZ_m  = PositionXYZ;
    TVector3 PositionXYZ_p2 = PositionXYZ;
    TVector3 PositionXYZ_m2 = PositionXYZ;

    if(counter == 1 && DoRichardson) {
      PositionXYZ_p(idx) += PositionXYZ(idx) + 0.5*h1;
      PositionXYZ_m(idx) += PositionXYZ(idx) - 0.5*h1;
    }
    
    PositionXYZ_p2(idx) = PositionXYZ(idx) + 0.5*h2;
    PositionXYZ_m2(idx) = PositionXYZ(idx) - 0.5*h2;
    
    if(DoRichardson) { 
      if(counter == 1) {
        Current_der1  = GetUVWFromXYZ(PositionXYZ_p);
        Current_der1 -= GetUVWFromXYZ(PositionXYZ_m);
        Current_der1 *= (1.0/h1);
      }
      else Current_der1 = Prev_der;
    }
    
    Current_der2  = GetUVWFromXYZ(PositionXYZ_p2);
    Current_der2 -= GetUVWFromXYZ(PositionXYZ_m2);
    Current_der2 *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der = (1.0/3.0)*(4.0*Current_der2 - Current_der1);
      if(counter == 1) Prev_der = Current_der2;
    }
    else Current_der = Current_der2;
    
    diff = Current_der - Prev_der;
    if(Current_der.Mag() > 1.0e-8) diff *= (1.0/Current_der.Mag());
    
    h /= 2.0;
    Prev_der = Current_der;
  }
  while(counter <= Nmax_iterations && diff.Mag() > epsilon_derivatives);
  
  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GSurfaceObject::GetDerUVFromXYZ_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << diff.Mag() << endl;
  }
  
  return Current_der;
  
}
//====================================================================
TVector3  GSurfaceObject::GetNormVector(TVector3 PositionXYZ)
{
  
  //Get normal vector to surface at global position PositionXYZ
  
  TVector3 PosAboveSurfUVW = GetUVWFromXYZ(PositionXYZ) + TVector3(0.0,0.0,1.0*global->GetUnit("um"));
  TVector3 PosAboveSurfXYZ = GetXYZFromUVW(PosAboveSurfUVW);
  
  TVector3  NormVector = (PosAboveSurfXYZ - PositionXYZ).Unit();
  
  return NormVector;
  
}
//====================================================================
bool  GSurfaceObject::IsInVoxel(Voxel_t aVoxel)
{
  
  //Check if surface is inside voxel
  
  std::vector<TVector3>  TestPositionList;
  TestPositionList.clear();
  
  double RU[2];
  double RV[2];
  GetVRange(RV[0],RV[1]);
  for(int iv=0;iv<NpointsV;iv++) {
    double V = RV[0] + (iv+0.5)*(RV[1] - RV[0])/NpointsV;
    GetURange(V,RU[0],RU[1]);
    for(int iu=0;iu<NpointsU;iu++) {
      double U = RU[0] + (iu+0.5)*(RU[1] - RU[0])/NpointsU;
      
      TestPositionList.push_back(GetXYZFromUVW(TVector3(U,V,0)));
    }
  }
  
  for(int ipoint=0;ipoint<int(TestPositionList.size());ipoint++) {
    if(global->IsPointInVoxel(aVoxel,TestPositionList[ipoint])) return true;
  }
  TestPositionList.clear();
  
  return false;
  
}
//====================================================================


