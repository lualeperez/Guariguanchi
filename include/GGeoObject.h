/***************************************************************************//**
 * @brief:      
 * @Description: Base class for geometrical objects
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
 

#ifndef GGeoObject_h
#define GGeoObject_h

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
#include "include/GSurfaceObject.h"
#include "include/GGlobalTools.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GGeoObject {
  
private:

public:
  
  int       GeometryID;          // Index of geometry containing this geometry object
  int       GeoObjectIdx;        // Index of this object inside geometry
  TString   GeoObjectName;       // Name
  TString   GeoObjectType;       // Type: plane, cylinder, disk, cone, petal, ...
  TString   LadderType;          // LadderType: Spiral or alternating
  TString   LayerName;           // Layer  to which this geometry object belongs
  TString   SystemName;          // System to which this geometry object belongs
  
  bool      IsSensitive;         // bool specifying if element is sensitive
  float     Thickness;           // Thickness
  TString   Material;            // Material
  float     X0;                  // Radiation length
  float     XOX0;                // Thickness in radiation length units
  float     ResolutionU;         // spatial resolution in U direction
  float     ResolutionV;         // spatial resolution in V direction
  TString   ResolutionModel;     // Resolution model name, mainly used for special cases
  int       ResolutionModelID;   // Resolution model index, mainly used for special cases
  TString   EfficiencyModel;     // Efficiency model name, mainly used for special cases
  int       EfficiencyModelID;   // Efficiency model index, mainly used for special cases
  float     DetEffic;            // detection efficiency in case of sensitive element
  float     ROtime;              // readout time in case of sensitive element
  float     BkgRate;             // background rate density (average over whole sensitive element)
  
  TVector3  Position;            // geometry element position
  TMatrixD  Rot;                 // Surface rotation matrix
  
  GSurfaceObject*  MainSurface;
  std::vector<GSurfaceObject*> BoundarySurfacesList;
  
  GGlobalTools* global;
  
  // Constructor
  GGeoObject(TString   aName,
	     int       aGeometryID,
	     int       aGeoObjectIdx,
             TVector3  aPosition,
             TMatrixD  aRot,
	     float     aThickness,
	     TString   aMaterial,
	     bool      aIsSensitive,
	     GGlobalTools* aglobal,
	     float     aResolutionU = -1.0,
	     float     aResolutionV = -1.0,
	     float     aDetEffic    = 0.0,
	     float     aROtime      = -1,
	     float     aBkgRate     = 0.0);
  
  GGeoObject(const GGeoObject& other, TString Name = TString(""));

  // Destructor
  virtual ~GGeoObject();
  
  //Functions to access internal variables
  int       GetMotherGeometryID()   const { return  GeometryID; }
  int       GetIndex()              const { return  GeoObjectIdx; }
  TString   GetName()               const { return  GeoObjectName; }
  TString   GetType()               const { return  GeoObjectType; }
  TString   GetSystemName()         const { return  SystemName; }
  TString   GetLayerName()          const { return  LayerName; }
  
  bool      GetIsSensitive()        const { return  IsSensitive; }
  TString   GetMaterial()           const { return  Material; }
  float     GetThickness()          const { return  Thickness; }
  float     GetX0()                 const { return  X0; }
  float     GetXOX0()               const { return  XOX0; }
  float     GetResolutionU()        const;
  float     GetResolutionV()        const;
  TVector2  GetResolutionUV()       const;
  TString   GetResolutionModel()    const { return  ResolutionModel; }
  int       GetResolutionModelID()  const { return  ResolutionModelID; }
  TString   GetEfficiencyModel()    const { return  EfficiencyModel; }
  int       GetEfficiencyModelID()  const { return  EfficiencyModelID; }
  float     GetDetEfficiency()      const { return  DetEffic; }
  float     GetROtime()             const { return  ROtime; }
  float     GetBkgRate()            const { return  BkgRate; }
  int       GetNBoundarySurfaces()  const { return  BoundarySurfacesList.size();}
  
  TString   GetLadderType()         const { return  LadderType; }
  
  GSurfaceObject*  GetMainSurface() { return  MainSurface; }
  GSurfaceObject*  GetBoundarySurface(int idx);
  
  double    GetBkgDensity()         const { return BkgRate*ROtime; }
  
  double    GetMass();
  
  bool      GetVerbose() { return verbose; }
  
  //Functions to set internal variables
  void      SetMotherGeometryID(int aGeometryID)     { GeometryID        = aGeometryID; }
  void      SetIndex(int aIndex)                     { GeoObjectIdx      = aIndex; }
  void      SetName(TString aName)                   { GeoObjectName     = aName; }
  void      SetType(TString aType)                   { GeoObjectType     = aType; }
  void      SetSystemName(TString aSystemName)       { SystemName        = aSystemName; }
  void      SetLayerName(TString aLayerName)         { LayerName         = aLayerName; }
  
  void      SetIsSensitive(bool aIsSensitive)        { IsSensitive       = aIsSensitive; }
  void      SetMaterial(TString aMaterial)           { Material          = aMaterial; }
  void      SetThickness(float aThickness)           { Thickness         = aThickness; }
  void      SetX0(float aX0)                         { X0                = aX0; }
  void      SetXOX0(float aXOX0)                     { XOX0              = aXOX0; }
  void      SetResolutionU(float aResolU)            { ResolutionU       = aResolU; }
  void      SetResolutionV(float aResolV)            { ResolutionV       = aResolV; }
  void      SetResolutionModel(TString aResolModel)  { ResolutionModel   = aResolModel; }
  void      SetResolutionModelID(int aResolModelID)  { ResolutionModelID = aResolModelID; }
  void      SetEfficiencyModel(TString aEfficModel)  { EfficiencyModel   = aEfficModel; }
  void      SetEfficiencyModelID(int aEfficModelID)  { EfficiencyModelID = aEfficModelID; }
  void      SetDetEfficiency(float aDetEffic)        { DetEffic          = aDetEffic; }
  void      SetROtime(float aROtime)                 { ROtime            = aROtime; }
  void      SetBkgRate(float aBkgRate)               { BkgRate           = aBkgRate; }
  
  void      SetLadderType(TString aLadderType);
  
  void      SetVerbose(bool aVerbose) { verbose = aVerbose; }
    
  bool      IsInVoxel(Voxel_t aVoxel);
  
  void      SetAllSurfacesGeoObjectType(void);
  
  //Set of functions to be defined for the daughter classes
  virtual  GGeoObject* clone(TString aName) const;
  virtual  double    GetVolume();
  virtual  void      GetThicknesses(double& w1, double& w2) {w1 = Thickness; w2 = 0.0; }
  virtual  void      CheckInputs();
  virtual  void      FillSurfaces() {;}
  virtual  bool      IsPointInsideGeometry(TVector3 PosXYZ) { return false;}
  virtual  TVector3  GetBoundaryNormVector(int idx, TVector3 PosXYZ);
  virtual  void      Print() {;}

protected:

  bool verbose;
  
};

#endif //~ GGeoObject_h

