/***************************************************************************//**
 * @brief:      
 * @Description: Class for full geometry handling
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
 

#ifndef GGeometry_h
#define GGeometry_h

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
#include "include/GGeoObject.h"
#include "include/GGeoPlane.h"
#include "include/GGeoCylinder.h"
#include "include/GGeoCylinderSection.h"
#include "include/GGeoDisk.h"
#include "include/GGeoDiskSection.h"
#include "include/GGeoPetal.h"
#include "include/GGeoCone.h"
#include "include/GGeoConeSection.h"
#include "include/GSurfaceObject.h"
#include "include/GSurfacePlane.h"
#include "include/GSurfaceCylinder.h"
#include "include/GSurfaceDisk.h"
#include "include/GSurfaceCone.h"
#include "include/GSurfacePetal.h"
#include "include/GGlobalTools.h"
#include "include/GBField.h"
#include "include/GResolutionModel.h"
#include "include/GEfficiencyModel.h"
#include "include/GUtiliratyFunctions.h"
#include "include/GTrackFinderAlgo.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GGeometry {
  
private:

public:
  
  int       GeometryID;          // Index of geometry in geometry list
  TString   GeometryName;        // Name
  
  GGeoObject* WorldVolume;       // World volume
  std::vector<GGeoObject*>  GeometryElementsList; //List with all the geometry elements of geometry
  
  GGlobalTools* global;
  
  GBField* Bfield;               // Bfield
  
  std::vector<GResolutionModel*>  ResolutionModelList; //List of resolution models for this geometry
  std::vector<GEfficiencyModel*>  EfficiencyModelList; //List of efficiency models for this geometry
  
  std::vector<TrackCuts_t>  TrackCutsList;  // List of track cuts (depending on track params and number of hits)
  
  std::vector<int>          VoxeledGeoElementsList;
  
  //Only for the case of a beam telescope analysis
  std::vector<TString>      TelescopePlanesLayerList;   //List of geometry of layers names used as Telescope
  std::vector<TString>      DUTPlanesLayerList;         //List of geometry of layers names considered as DUTs
  std::vector<int>          TelescopePlanesIndexList;   //List of geometry elements which are used as Telescope
  std::vector<int>          DUTPlanesIndexList;         //List of geometry elements which are DUTs
  
  //Systems Definition within geometry
  std::vector<System_t>  SystemsList;      //List of system objects
  std::vector<TString>   SystemNamesList;  //List of system names
  
  std::vector<GTrackFinderAlgo*> TrackFinderAlgoList; //TrackFinder algorithms list
  
  double GeoCheckPrecision;
  double BkgRateScaling;
  
  // Constructor
  GGeometry(TString   aName,
	    int       aGeometryID,
	    GGlobalTools* aglobal);

  // Destructor
  virtual ~GGeometry();
  
  //Functions to access internal variables
  int       GetGeometryID()          const { return  GeometryID; }
  TString   GetName()                const { return  GeometryName; }
  double    GetGeoCheckPrecision()   const { return  GeoCheckPrecision; }
  double    GetBkgRateScaling()      const { return  BkgRateScaling; }
  
  int       GetNGeoElements()        const { return  GeometryElementsList.size(); }
  int       GetNVoxelesGeoElements() const { return  VoxeledGeoElementsList.size(); }
  int       GetNTrackCuts()          const { return  TrackCutsList.size(); }
  
  //Beam-telescope analysis
  int       GetNTelescopePlanes()    const { return  TelescopePlanesIndexList.size(); }
  int       GetNDUTPlanes()          const { return  DUTPlanesIndexList.size(); }
  int       GetTelescopePlaneIndexInGeometry(int idx);
  int       GetDUTPlaneIndexInGeometry(int idx);
  
  //System list
  int       GetNSystems()             const { return  SystemsList.size(); }
  System_t  GetASystemElement(int idx);
  void      PushASystemIntoGeometry(System_t aSystem)  { SystemsList.push_back(aSystem); }
  
  int       GetNSystemNames()         const { return  SystemNamesList.size(); }
  TString   GetASystemName(int idx);
  
  GGeoObject*  GetWorldVolume()  { return  WorldVolume; }
  GGeoObject*  GetGeometryElement(int idx);
  GGeoObject*  GetVoxeledGeometryElement(int idx);
  GBField*     GetBField()       { return  Bfield; }
  
  int               GetNResolutionModels()  const { return ResolutionModelList.size(); }
  int               GetNEfficiencyModels()  const { return EfficiencyModelList.size(); }
  int               GetNTrackFinderAlgos()  const { return TrackFinderAlgoList.size(); }
  GResolutionModel* GetResolutionModel(int idx);
  GEfficiencyModel* GetEfficiencyModel(int idx);
  GTrackFinderAlgo* GetTrackFinderAlgo(int idx);
  TrackCuts_t       GetTrackCut(int idx);
  
  //Functions to set internal variables
  void      SetGeometryID(int aGeometryID)                    { GeometryID        = aGeometryID; }
  void      SetName(TString aName)                            { GeometryName      = aName; }
  void      SetGeoCheckPrecision(double aPrecision)           { GeoCheckPrecision = aPrecision; }
  void      SetBkgRateScaling(double scaling)                 { BkgRateScaling    = scaling; }
  
  void      PushGeoElement(GGeoObject* aGeoObject);
  void      SetWorldVolume(GGeoObject* aGeoObject);
  void      SetBField(GBField* aBfield);
  void      PushResolutionModelIntoGeometry(GResolutionModel* aResoModel);
  void      PushEfficiencyModelIntoGeometry(GEfficiencyModel* aEfficModel);
  void      PushTrackFinderAlgoIntoGeometry(GTrackFinderAlgo* aTrackFinderAlgo);
  void      PushTrackCutIntoGeometry(TrackCuts_t aTrackCut)   { TrackCutsList.push_back(aTrackCut); }
  
  void      ApplyBkgScaling();
  
  void      FillVoxeledGeoElementsList(std::vector<Voxel_t> VoxelList);
  
  void      FillResolutionModelIndexes();
  
  void      FillEfficiencyModelIndexes();
  
  void      FillBeamTestConfigLayers(std::vector<TString>  aTelescopeLayersList,
				     std::vector<TString>  aDUTLayersList);
  void      FillBeamTestConfigPlanesIndexes();
  
  void      FillSystemsList(void);
  
  int       FindTrackFinderAlgo(TVector3 x0, TVector3 p0);
  
  void      CheckGeometry();
  void      Print();
  void      PrintWeight();
  
  TVector2  GetResolutionUV(int GeoElementID,TVector3 mom, TVector3 pos);
  
  double    GetEfficiency(int GeoElementID,TVector3 mom, TVector3 pos, TString aParticle);
  
  void      ApplyTrackCuts(TVector3 InitMomentum, std::vector<IntersectionHit_t>& ItersectionHitList);
  void      ApplyTrackCuts2(TVector3 InitMomentum, std::vector<IntersectionHit_t>& ItersectionHitList);
  
protected:

  bool verbose;
  
};

#endif //~ GGeometry_h

