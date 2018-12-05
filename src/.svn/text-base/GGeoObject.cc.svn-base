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
#include "include/GGeoObject.h"

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
GGeoObject::GGeoObject(TString   aName,
		       int       aGeometryID,
		       int       aGeoObjectIdx,
		       TVector3  aPosition,
		       TMatrixD  aRot,
		       float     aThickness,
		       TString   aMaterial,
		       bool      aIsSensitive,
		       GGlobalTools* aglobal,
		       float     aResolutionU,
		       float     aResolutionV,
		       float     aDetEffic,
		       float     aROtime,
		       float     aBkgRate)
{
  
  global = aglobal;
  
  GeoObjectName     = aName;
  GeometryID        = aGeometryID;
  GeoObjectIdx      = aGeoObjectIdx;
  Position          = aPosition;
  Thickness         = aThickness;
  Material          = aMaterial;
  IsSensitive       = aIsSensitive;
  
  ResolutionU       = aResolutionU;
  ResolutionV       = aResolutionV;
  DetEffic          = aDetEffic;
  ROtime            = aROtime;
  BkgRate           = aBkgRate;
  ResolutionModel   = TString("");
  ResolutionModelID = -1;
  EfficiencyModel   = TString("");
  EfficiencyModelID = -1;
  
  Rot.ResizeTo(3,3);
  Rot    = aRot;
  
  X0   = global->GetX0FromMaterial(Material);
  XOX0 = Thickness/X0;
  
  CheckInputs();
  
  LadderType = TString("");
  SystemName = TString("");
  LayerName  = TString("");
  
  MainSurface = NULL;
  BoundarySurfacesList.clear();
  
  FillSurfaces();
  
  verbose = false;
    
}
//====================================================================
GGeoObject::GGeoObject(const GGeoObject& other,TString aName)
{
  
  global            = other.global;
  
  GeoObjectName     = aName;
  GeometryID        = other.GeometryID;
  GeoObjectIdx      = other.GeoObjectIdx;
  Position          = other.Position;
  Thickness         = other.Thickness;
  Material          = other.Material;
  IsSensitive       = other.IsSensitive;
  
  ResolutionU       = other.ResolutionU;
  ResolutionV       = other.ResolutionV;
  DetEffic          = other.DetEffic;
  ROtime            = other.ROtime;
  BkgRate           = other.BkgRate;
  ResolutionModel   = other.ResolutionModel;
  ResolutionModelID = other.ResolutionModelID;
  
  Rot.ResizeTo(3,3);
  Rot    = other.Rot;
  
  X0   = global->GetX0FromMaterial(Material);
  XOX0 = Thickness/X0;
  
  CheckInputs();
  
  LadderType = TString("");
  SystemName = TString("");
  LayerName  = TString("");
  
  MainSurface = NULL;
  BoundarySurfacesList.clear();
  
  FillSurfaces();
  
  verbose = other.verbose;
    
}
//====================================================================
GGeoObject::~GGeoObject() 
{

  //Freeing memory space
  if(MainSurface != NULL) delete MainSurface;
  for(int i=0;i<int(BoundarySurfacesList.size());i++) delete BoundarySurfacesList[i];
  BoundarySurfacesList.clear();
  
}
//====================================================================
GGeoObject* GGeoObject::clone(TString aName) const
{
 
  return new GGeoObject(*this,aName);
  
}
//====================================================================
void  GGeoObject::SetLadderType(TString aLadderType)
{
  
  LadderType = aLadderType;
  
  if(MainSurface != NULL) MainSurface->SetLadderType(aLadderType);
  for(int i=0;i<int(BoundarySurfacesList.size());i++) BoundarySurfacesList[i]->SetLadderType(aLadderType);
  
  return;
  
}
//====================================================================
void  GGeoObject::SetAllSurfacesGeoObjectType(void)
{
  
  MainSurface->SetGeoObjectType(GeoObjectType);
  for(int i=0;i<int(BoundarySurfacesList.size());i++) BoundarySurfacesList[i]->SetGeoObjectType(GeoObjectType);
  
}
//====================================================================
void  GGeoObject::CheckInputs()
{
  
  if(GeometryID < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: GeometryID < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(GeoObjectIdx < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: GeoObjectIdx < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Thickness < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: Thickness < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(!IsSensitive) return;
  
  if(ResolutionU < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: ResolutionU < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(ResolutionV < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: ResolutionV < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(DetEffic < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: DetEffic < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(BkgRate < 0) {
    cout << endl;
    cout << "GGeoObject::CheckInputs: BkgRate < 0 for GeoObject " << GeoObjectName.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
    
  return;
  
}
//====================================================================
float GGeoObject::GetResolutionU() const
{ 
  
  if(!IsSensitive) return 0.0;
  else             return  ResolutionU;
  
}
//====================================================================
float  GGeoObject::GetResolutionV() const
{
 
  if(!IsSensitive) return 0.0;
  else             return  ResolutionV;
  
}
//====================================================================
TVector2  GGeoObject::GetResolutionUV() const 
{
  
  return  TVector2(GetResolutionU(),GetResolutionV());
  
}
//====================================================================
GSurfaceObject* GGeoObject::GetBoundarySurface(int idx)
{
  
  if(idx < 0 || idx > int(BoundarySurfacesList.size()-1)) {
    cout << endl;
    cout << "GGeoObject::GetBoundarySurface:: index " << idx << " is out of boundaries for (0," << BoundarySurfacesList.size()-1 << ") for " 
         << GeoObjectName.Data() << " geometry objet. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  BoundarySurfacesList[idx];
  
}
//====================================================================
bool  GGeoObject::IsInVoxel(Voxel_t aVoxel)
{
  
  if(MainSurface->IsInVoxel(aVoxel)) return true;
  else {
    for(int iboundary=0;iboundary<int(BoundarySurfacesList.size());iboundary++) {
      if(BoundarySurfacesList[iboundary]->IsInVoxel(aVoxel)) return true;
    }
  }
  
  return  false;
  
}
//====================================================================
double  GGeoObject::GetVolume()
{
  
  return  -999.9;
  
}
//====================================================================
TVector3  GGeoObject::GetBoundaryNormVector(int idx, TVector3 PosXYZ)
{
  
  if(idx < 0 || idx > int(BoundarySurfacesList.size()-1)) {
    cout << endl;
    cout << "GGeoObject::GetBoundaryNormVector:: index " << idx << " is out of boundaries for (0," << BoundarySurfacesList.size()-1 << ") for " 
         << GeoObjectName.Data() << " geometry objet. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  BoundarySurfacesList[idx]->GetWVector();
  
}
//====================================================================
double  GGeoObject::GetMass() {
  
  double density = global->GetDensityFromMaterial(Material)/global->GetUnit("gr/cm3");
  double volume  = GetVolume()/global->GetUnit("cm3");
  double mass    = density*volume;
  
//   cout << endl;
//   cout << "GeoObject " << GeoObjectName.Data() << endl;
//   cout << "  Volume   = " << volume << " cm3" << endl;
//   cout << "  Material = " << Material.Data() << endl;
//   cout << "  Density  = " << density << " gr/cm3" << endl;
//   cout << "  mass     = " << mass << " gr" << endl;
//   cout << endl;
//   cout << endl;
  
  return  mass;
  
}
//====================================================================