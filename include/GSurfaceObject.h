/***************************************************************************//**
 * @brief:      
 * @Description: Base class for surface objects
 *
 *
 * @createdby:  PEREZ PEREZ Alejandro <luis_alejandro.perez_perez@iphc.cnrs.fr> at 2017-10-01 14:07:38
 * @copyright:  (c)2017 IPHC - CNRS - Université de Strasbourg. All Rights Reserved.
 * 
 * @License: You are free to use this source files for your own development as long
 *           as it stays in a public research context. You are not allowed to use it
 *           for commercial purpose. You must put this header with laboratory and
 *           authors names in all development based on this library.
 *           When results obtained with this package (Guariguanchi) are communicated, 
 *           you should quote the package name.
 *
 * @lastchange: $Revision$
 *              $Author$
 *              $Date$
 *
 *******************************************************************************/
 

#ifndef GSurfaceObject_h
#define GSurfaceObject_h

#include "TMath.h"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TCutG.h>
#include "TStopwatch.h"
#include "include/GGlobalTools.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

const int NpointsU = 10;
const int NpointsV = 10;

using namespace std;

class GSurfaceObject {
  
private:
  
public:
  
  TString   Name;                // Surface name
  TString   Type;                // Surface type: plane, cyinders, disk, cone, ...
  TString   LadderType;          // LadderType: Spiral or alternating
  TString   GeoObjectType;       // Type of the geo object owning this surface

  bool      IsTrackingLayer;     // Is surface a tracking element
  TVector3  Position;            // Surface object position
  TVector3  UVector;             // unit vector in the tangent U direction
  TVector3  VVector;             // unit vector in the tangent V direction
  TVector3  WVector;             // unit vector in the normal  W direction
  TMatrixD  Rot;                 // Surface rotation matrix
  TMatrixD  InvRot;              // Surface inverse  rotation matrix
  
  GGlobalTools* global;
  
  // Constructor
  GSurfaceObject(TString   aName,
		 TVector3  aPosition,
		 TMatrixD  aRot,
		 bool      aIsTrackingLayer,
		 GGlobalTools* aglobal);
  
  GSurfaceObject(const GSurfaceObject& other, TString Name = TString(""));

  // Destructor
  virtual ~GSurfaceObject();
  
  //Set of generic functions
  
  //Functions to access internal variables
  TString   GetName()               const { return  Name; }
  TString   GetType()               const { return  Type; }
  
  TVector3  GetPosition()           const { return  Position; }
  TVector3  GetUVector()            const { return  UVector; }
  TVector3  GetVVector()            const { return  VVector; }
  TVector3  GetWVector()            const { return  WVector; }
  TMatrixD  GetRotMatrix()          const { return  Rot; }
  TMatrixD  GetInvRotMatrix()       const { return  InvRot; }
  
  TString   GetLadderType()         const { return  LadderType; }
  
  TString   GetGeoObjectType()      const { return  GeoObjectType; }
  
  //Functions to set internal variables
  void      SetName(TString aName)                   { Name              = aName; }
  void      SetType(TString aType)                   { Type              = aType; }
  
  void      SetPosition(TVector3 aPosition)          { Position          = aPosition; }
  void      SetUVector(TVector3 uVect)               { UVector           = uVect; }
  void      SetVVector(TVector3 vVect)               { VVector           = vVect; }
  void      SetWVector(TVector3 wVect)               { WVector           = wVect; }
  
  void      SetRotMatrix(TMatrixD R);
  
  void      SetVerbose(bool aVerbose) { verbose = aVerbose; }
  
  // Function for the numerical calculation of U,V,W coordinates w.r.t X,Y,Z
  TVector3 GetDerUVWFromXYZ_Numerical(TVector3  PositionXYZ, int idx);
  
  void     CheckIndex(int idx);
  
  void     SetLadderType(TString aLadderType)       { LadderType        = aLadderType; }
  
  void     SetGeoObjectType(TString aGeoObjectType) { GeoObjectType     = aGeoObjectType; }
  
  bool     GetVerbose() { return verbose; }
  
  bool     IsInVoxel(Voxel_t aVoxel);
  
  //Set of functions to be defined for the daughter classes
  virtual  GSurfaceObject* clone(TString aName) const;
  virtual  void      CheckInputs() {;}
  virtual  TVector3  GetInitPositionForIntersection() { return TVector3(0,0,0); }
  virtual  TVector3  GetUVWFromXYZ(TVector3 PositionXYZ) { return TVector3(-999.9,-999.9,-999.9); }
  virtual  TVector3  GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ,
						 int idx) { return TVector3(-999.9,-999.9,-999.9); }
  virtual  TVector3  GetXYZFromUVW(TVector3 PositionUVW)  { return TVector3(-999.9,-999.9,-999.9); }
  
  virtual  bool      IsInMaterial(TVector3 PositionUVW)  { return false; }
  virtual  bool      IsInSensitiveMaterial(TVector3 PositionUVW) { return false; }
  virtual  TVector3  GetNormVector(TVector3 PositionXYZ);
  virtual  void      VarySurface(double w1, double w2) {;}
  virtual  void      GetVRange(double& Vmin, double& Vmax) {;}
  virtual  void      GetURange(double V, double& Umin, double& Umax) {;}
  virtual  void      FillSurfaceRepresentation(int color,
					       std::vector<TGraph>   &GraphXY,
					       std::vector<TGraph>   &GraphZY,
					       std::vector<TGraph>   &GraphZX,
					       std::vector<TGraph>   &GraphZR,
					       int linestyle = 1) {;}

protected:

  bool verbose;
  
};

#endif //~ GSurfaceObject_h

