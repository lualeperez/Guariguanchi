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
#include "include/GGeometry.h"

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
GGeometry::GGeometry(TString   aName,
		     int       aGeometryID,
		     GGlobalTools* aglobal)
{
  
  global = aglobal;
  
  GeometryName     = aName;
  GeometryID       = aGeometryID;
  
  WorldVolume = NULL;
  GeometryElementsList.clear();
  Bfield      = NULL;
  
  EfficiencyModelList.clear();
  ResolutionModelList.clear();
  TrackCutsList.clear();
  VoxeledGeoElementsList.clear();
  
  TelescopePlanesLayerList.clear();
  DUTPlanesLayerList.clear();
  TelescopePlanesIndexList.clear();
  DUTPlanesIndexList.clear();  
  
  SystemsList.clear();
  SystemNamesList.clear();
  
  TrackFinderAlgoList.clear();
  
  GeoCheckPrecision = 1.0*global->GetDistanceUnit("mm");
  
  BkgRateScaling    = 1.0;
  
}
//====================================================================
GGeometry::~GGeometry() 
{
   
  //Freeing memory space
  if(WorldVolume != NULL) delete WorldVolume;
  for(int i=0;i<int(GeometryElementsList.size());i++) delete GeometryElementsList[i];
  GeometryElementsList.clear();
  
  if(Bfield != NULL) delete Bfield;
  
}
//====================================================================
void   GGeometry::PushGeoElement(GGeoObject* aGeoObject)
{
 
  //Push a GGeoObject inside geometry elements list
  
  if(aGeoObject == NULL) return;
  
  int index = GeometryElementsList.size();
  GeometryElementsList.push_back(aGeoObject);
  GeometryElementsList[index]->SetIndex(index);
  
  VoxeledGeoElementsList.push_back(index);
  
  return;
  
}
//====================================================================
void   GGeometry::PushResolutionModelIntoGeometry(GResolutionModel* aResoModel)
{ 
  
  //Push a GResolutionModel inside geometry's resolution model list
  
  if(aResoModel == NULL) return;
  
  ResolutionModelList.push_back(aResoModel);
  
  return;
  
}
//====================================================================
void   GGeometry::PushEfficiencyModelIntoGeometry(GEfficiencyModel* aEfficModel)
{ 
  
  //Push a GEfficiencyModel inside geometry's resolution model list
  
  if(aEfficModel == NULL) return;
  
  EfficiencyModelList.push_back(aEfficModel);
  
  return;
  
}
//====================================================================
void  GGeometry::PushTrackFinderAlgoIntoGeometry(GTrackFinderAlgo* aTrackFinderAlgo)
{
  
  //Push a GTrackFinderAlgo inside geometry's resolution model list
  
  if(aTrackFinderAlgo == NULL) return;
  
  TrackFinderAlgoList.push_back(aTrackFinderAlgo);
  
  return;
  
}
//====================================================================
void  GGeometry::ApplyBkgScaling()
{
  
  if(BkgRateScaling < 0)    return;
  if(BkgRateScaling == 1.0) return;
  
  for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) {
    if(!GeometryElementsList[igeoElm]->GetIsSensitive()) continue;
    
    GeometryElementsList[igeoElm]->SetBkgRate(BkgRateScaling*GeometryElementsList[igeoElm]->GetBkgRate());
  }
  
  return; 
  
}
//====================================================================
void   GGeometry::FillVoxeledGeoElementsList(std::vector<Voxel_t> VoxelList)
{
  
  //Fill list of geometry elements inside the VolexList
  
  if(VoxelList.size() == 0) return;
  
  for(int i=0;i<int(GeometryElementsList.size());i++) {    
    bool IsInVoxelList = false;
    for(int ivoxel=0;ivoxel<int(VoxelList.size());ivoxel++) {
      if(GeometryElementsList[i]->IsInVoxel(VoxelList[ivoxel])) {
	IsInVoxelList = true;
	break;
      }
    }
    
    if(!IsInVoxelList) {
      int index = -999;
      for(int j=0;j<int(VoxeledGeoElementsList.size());j++) {
	if(VoxeledGeoElementsList[j] == i) {
	  index = j;
	  break;
	}
      }
      VoxeledGeoElementsList.erase(VoxeledGeoElementsList.begin() + index);
    }
  }
  
  return;
  
}
//====================================================================
void  GGeometry::FillResolutionModelIndexes()
{
  
  //Set the geometry elements resolution model indexes
  
  if(ResolutionModelList.size() == 0) return;
  
  for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) {
    int idx = -1;
    for(int iresModel=0;iresModel<int(ResolutionModelList.size());iresModel++) {
      if(GeometryElementsList[igeoElm]->GetResolutionModel() == ResolutionModelList[iresModel]->GetName()) {
	idx = iresModel;
	break;
      }
    }
    GeometryElementsList[igeoElm]->SetResolutionModelID(idx);
  }
  
  return;
  
}
//====================================================================
void  GGeometry::FillEfficiencyModelIndexes()
{
  
  //Set the geometry elements efficiency model indexes
  
  if(EfficiencyModelList.size() == 0) return;
  
  for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) {
    int idx = -1;
    for(int iefficModel=0;iefficModel<int(EfficiencyModelList.size());iefficModel++) {
      if(GeometryElementsList[igeoElm]->GetEfficiencyModel() == EfficiencyModelList[iefficModel]->GetName()) {
	idx = iefficModel;
	break;
      }
    }
    GeometryElementsList[igeoElm]->SetEfficiencyModelID(idx);
  }
  
  return;
  
}
//====================================================================
GGeoObject*  GGeometry::GetGeometryElement(int idx)
{
  
  //Get geometry element with index idx
  
  if(idx < 0 || idx > int(GeometryElementsList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetGeometryElement:: index " << idx << " is out of boundaries for (0," << GeometryElementsList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  GeometryElementsList[idx];
  
}
//====================================================================
GGeoObject*  GGeometry::GetVoxeledGeometryElement(int idx)
{
  
  //Get voxeled geometry element with index idx
  
  if(idx < 0 || idx > int(VoxeledGeoElementsList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetVoxeledGeometryElement:: index " << idx << " is out of boundaries for (0," << VoxeledGeoElementsList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  GeometryElementsList[VoxeledGeoElementsList[idx]];
  
}
//====================================================================
GResolutionModel*  GGeometry::GetResolutionModel(int idx)
{
  
  //Get resolution model with index idx
  
  if(idx < 0 || idx > int(ResolutionModelList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetResolutionModel:: index " << idx << " is out of boundaries for (0," << ResolutionModelList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return ResolutionModelList[idx];
  
}
//====================================================================
GEfficiencyModel*  GGeometry::GetEfficiencyModel(int idx)
{
  
  //Get efficiency model with index idx
  
  if(idx < 0 || idx > int(EfficiencyModelList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetEfficiencyModel:: index " << idx << " is out of boundaries for (0," << EfficiencyModelList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  EfficiencyModelList[idx];
  
}
//====================================================================
GTrackFinderAlgo*  GGeometry::GetTrackFinderAlgo(int idx)
{
  
  //Get efficiency model with index idx
  
  if(idx < 0 || idx > int(TrackFinderAlgoList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetTrackFinderAlgo:: index " << idx << " is out of boundaries for (0," << TrackFinderAlgoList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  TrackFinderAlgoList[idx];
  
}
//====================================================================
TrackCuts_t  GGeometry::GetTrackCut(int idx)
{
  
  //Get track cut with index idx
  
  if(idx < 0 || idx > int(TrackCutsList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetTrackCut:: index " << idx << " is out of boundaries for (0," << TrackCutsList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }

  return TrackCutsList[idx];
  
}
//====================================================================
void  GGeometry::FillBeamTestConfigLayers(std::vector<TString>  aTelescopeLayersList,
					  std::vector<TString>  aDUTLayersList)
{
  
  //Fill the Telescope-DUT configuration layers
  
  if(aTelescopeLayersList.size() == 0) return;
  if(aDUTLayersList.size()       == 0) return;
  
  int N = 0;

  TelescopePlanesLayerList.clear();
  N = aTelescopeLayersList.size();
  for(int kkk=0;kkk<N;kkk++) TelescopePlanesLayerList.push_back(aTelescopeLayersList[kkk]);
  
  DUTPlanesLayerList.clear();
  N = aDUTLayersList.size();
  for(int kkk=0;kkk<N;kkk++) DUTPlanesLayerList.push_back(aDUTLayersList[kkk]);

  return;
  
}
//====================================================================
void  GGeometry::FillBeamTestConfigPlanesIndexes()
{
 
  //Fill telescope-DUT configuration planes indexes
  
  if(TelescopePlanesLayerList.size() == 0 || DUTPlanesLayerList.size() == 0) return;
  
  TelescopePlanesIndexList.clear();
  DUTPlanesIndexList.clear();
  for(int igeo=0;igeo<int(VoxeledGeoElementsList.size());igeo++) { //begin loop over voxeled geoElements
    int idx = VoxeledGeoElementsList[igeo];
    TString  LayerName  = GeometryElementsList[idx]->GetLayerName();
    
    bool IsInTelescopeList = false;
    bool IsInDUTList       = false;
    
    if(GeometryElementsList[idx]->GetIsSensitive()) {
      for(int kkk=0;kkk<int(TelescopePlanesLayerList.size());kkk++) {
        if(LayerName == TelescopePlanesLayerList[kkk]) {
	  IsInTelescopeList = true;
	  TelescopePlanesIndexList.push_back(idx);
	  break;
        }
      }
    }
    
    for(int kkk=0;kkk<int(DUTPlanesLayerList.size());kkk++) {
      if(LayerName == DUTPlanesLayerList[kkk]) {
	IsInDUTList = true;
	DUTPlanesIndexList.push_back(idx);
	break;
      }
    }
    
    if(IsInTelescopeList && IsInDUTList) {
      cout << endl;
      cout << "For geometry " << GeometryName.Data() << ", the geo-element " << idx << " is both in telescop and DUT lists. Check your inputs. Exiting now!!!" << endl;
      cout << endl;
      assert(false);
    }
  } // end of loop over voxeled geoElements

  if(TelescopePlanesIndexList.size() == 0 && DUTPlanesIndexList.size() > 0) {
    cout << endl;
    cout << "For geometry " << GeometryName.Data() << ", the telecope list is empty but the DUT list is filled. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  else if(TelescopePlanesIndexList.size() > 0 && DUTPlanesIndexList.size() == 0) {
    cout << endl;
    cout << "For geometry " << GeometryName.Data() << ", the DUT list is empty but the telescope list is filled. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  TelescopePlanesLayerList.clear();
  DUTPlanesLayerList.clear();

  return;
  
}
//====================================================================
int  GGeometry::GetTelescopePlaneIndexInGeometry(int idx)
{
 
  //Get telescope plane with index idx
  
  if(idx < 0 || idx > int(TelescopePlanesIndexList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetTelescopePlaneIndexInGeometry:: index " << idx << " is out of boundaries for (0," << TelescopePlanesIndexList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  TelescopePlanesIndexList[idx];
  
}
//====================================================================
int  GGeometry::GetDUTPlaneIndexInGeometry(int idx)
{
  
  //Get DUT plane with index idx
  
  if(idx < 0 || idx > int(DUTPlanesIndexList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetDUTPlaneIndexInGeometry:: index " << idx << " is out of boundaries for (0," << DUTPlanesIndexList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  DUTPlanesIndexList[idx];
  
}
//====================================================================
void  GGeometry::FillSystemsList()
{
  
  //Fille the geometry's systems name list
  
  //Fill the systems list
  if(SystemsList.size() > 0) { //begin if systems list is not empty
    for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) { //begin loop over geo-elements
      for(int isys=0;isys<int(SystemsList.size());isys++) { //begin loop over systems list
	
        bool IsFound = false;
        for(int ilayer=0;ilayer<int(SystemsList[isys].LayersList.size());ilayer++) { //begin loop over layer list 
	  if(SystemsList[isys].LayersList[ilayer] == GeometryElementsList[igeoElm]->GetLayerName()) {
	    GeometryElementsList[igeoElm]->SetSystemName(SystemsList[isys].Name);
	    IsFound = true;
	    break;
	  }
        } //end loop over layer list
        if(IsFound) break;
	
      } //end of loop over systems list
    } // end loop over geo-elements
  } //end if systems list is not empty
  SystemsList.clear();
  
  //Now filling the geometry's system name list
  //Filling list with materials in geometries
  SystemNamesList.clear();
  for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) { //begin loop over geo-elements
    TString system = GeometryElementsList[igeoElm]->GetSystemName();
    bool IsInList = false;
    if(system == TString("")) IsInList = true;
    for(int k=0;k<int(SystemNamesList.size());k++) {
      if(SystemNamesList[k] == system) {
	IsInList = true;
	break;
      }
    }
    if(!IsInList) SystemNamesList.push_back(system);
  } // end loop over geo-elements
  
  return;
  
}
//====================================================================
System_t  GGeometry::GetASystemElement(int idx)
{
  
  //Get system element with index idx
  
  if(idx < 0 || idx > int(SystemsList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetASystemElement:: index " << idx << " is out of boundaries for (0," << SystemsList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return SystemsList[idx];
  
}
//====================================================================
TString  GGeometry::GetASystemName(int idx)
{
  
  //Get system name with index idx
  
  if(idx < 0 || idx > int(SystemNamesList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetASystemName:: index " << idx << " is out of boundaries for (0," << SystemNamesList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  return  SystemNamesList[idx];
    
}
//====================================================================
TVector2  GGeometry::GetResolutionUV(int GeoElementID,TVector3 mom, TVector3 pos)
{
 
  //Get measured space point resolution for geometry element GeoElementID intersected at position pos with momentum mom
  
  if(GeoElementID < 0 || GeoElementID > int(GeometryElementsList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetResolutionUV:: index " << GeoElementID << " is out of boundaries for (0," << GeometryElementsList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  if(!GeometryElementsList[GeoElementID]->GetIsSensitive()) return TVector2(0.0,0.0);
  
  TVector2 ResolutionUV = GeometryElementsList[GeoElementID]->GetResolutionUV();
  
  int model_idx = GeometryElementsList[GeoElementID]->GetResolutionModelID();
  if(model_idx < 0) return  ResolutionUV;
  
  return GetResolutionModel(model_idx)->GetResolution(ResolutionUV,mom,pos);
  
}
//====================================================================
double  GGeometry::GetEfficiency(int GeoElementID,TVector3 mom, TVector3 pos, TString aParticle)
{
  
  //Get measured space point detection efficiency for geometry element GeoElementID intersected at position pos with momentum mom
  
  if(GeoElementID < 0 || GeoElementID > int(GeometryElementsList.size()-1)) {
    cout << endl;
    cout << "GGeometry::GetEfficiency:: index " << GeoElementID << " is out of boundaries for (0," << GeometryElementsList.size()-1 << ") for " 
         << GeometryName.Data() << " geometry. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  if(!GeometryElementsList[GeoElementID]->GetIsSensitive()) return 0.0;
  
  double DetEffic = GeometryElementsList[GeoElementID]->GetDetEfficiency();
  
  int model_idx = GeometryElementsList[GeoElementID]->GetEfficiencyModelID();
  if(model_idx < 0) return  DetEffic;
  
  return GetEfficiencyModel(model_idx)->GetEfficiency(DetEffic,mom,pos,aParticle);
  
}
//====================================================================
void   GGeometry::ApplyTrackCuts2(TVector3 InitMomentum,
				  std::vector<IntersectionHit_t>& ItersectionHitList)
{
  
  //Apply track cuts
  
  if(TrackCutsList.size() == 0) return;
  
  double p     = InitMomentum.Mag();
  double phi   = atan2(InitMomentum.Y(),InitMomentum.X()); //in (-pi,+pi) region
  double theta = atan2(sqrt(pow(InitMomentum.X(),2) + pow(InitMomentum.Y(),2)),InitMomentum.Z()); //in (0,+pi) region
    
  int idx_cut = -999;
  for(int icut=0;icut<int(TrackCutsList[0].CutList.size());icut++) {
    if((p     >= TrackCutsList[0].CutList[icut].pMin     && p     <= TrackCutsList[0].CutList[icut].pMax)   &&
       (phi   >= TrackCutsList[0].CutList[icut].phiMin   && phi   <= TrackCutsList[0].CutList[icut].phiMax) && 
       (theta >= TrackCutsList[0].CutList[icut].thetaMin && theta <= TrackCutsList[0].CutList[icut].thetaMax)) {
	idx_cut = icut;
      break;
    }
  }
  
  if(idx_cut < 0) return; 
  
  int Nhits    = 0;
  int min_ihit = +9999999;
  int max_ihit = -9999999;
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(!ItersectionHitList[ihit].IsSensitivePoint) continue;
    if(GetGeometryElement(ItersectionHitList[ihit].geoElement_idx)->GetSystemName() != TrackCutsList[0].SystemName) continue;
      
    Nhits++;
    if(min_ihit > ihit) min_ihit = ihit;
    if(max_ihit < ihit) max_ihit = ihit;
  }
  
  if(Nhits == 0) return;
  if(Nhits >= TrackCutsList[0].CutList[idx_cut].NhitsMin && Nhits <= TrackCutsList[0].CutList[idx_cut].NhitsMax) return;
  
  std::vector<IntersectionHit_t> ItersectionHitList_ppp;
  ItersectionHitList_ppp.clear();
      
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ihit >= min_ihit && ihit <= max_ihit) continue;
    ItersectionHitList_ppp.push_back(ItersectionHitList[ihit]);
  }
  ItersectionHitList.clear();
  for(int ihit=0;ihit<int(ItersectionHitList_ppp.size());ihit++) ItersectionHitList.push_back(ItersectionHitList_ppp[ihit]);
    
  return;
  
}
//====================================================================
void   GGeometry::ApplyTrackCuts(TVector3 InitMomentum,
				 std::vector<IntersectionHit_t>& ItersectionHitList)
{
  
  //Apply track cuts
  
  if(TrackCutsList.size() == 0) return;
  
  double p     = InitMomentum.Mag();
  double phi   = atan2(InitMomentum.Y(),InitMomentum.X()); //in (-pi,+pi) region
  double theta = atan2(sqrt(pow(InitMomentum.X(),2) + pow(InitMomentum.Y(),2)),InitMomentum.Z()); //in (0,+pi) region
    
  int idx_cut = -999;
  for(int icut=0;icut<int(TrackCutsList[0].CutList.size());icut++) {
    if((p     >= TrackCutsList[0].CutList[icut].pMin     && p     <= TrackCutsList[0].CutList[icut].pMax)   &&
       (phi   >= TrackCutsList[0].CutList[icut].phiMin   && phi   <= TrackCutsList[0].CutList[icut].phiMax) && 
       (theta >= TrackCutsList[0].CutList[icut].thetaMin && theta <= TrackCutsList[0].CutList[icut].thetaMax)) {
	idx_cut = icut;
      break;
    }
  }

  if(idx_cut < 0) return; 
  
  int Nhits = 0;
  int max_ihit = -999;
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(!ItersectionHitList[ihit].IsSensitivePoint) continue;
    if(GetGeometryElement(ItersectionHitList[ihit].geoElement_idx)->GetSystemName() != TrackCutsList[0].SystemName) continue;
      
    Nhits++;
    max_ihit = ihit;
      
    if(Nhits >= TrackCutsList[0].CutList[idx_cut].NhitsMin) break;
  }
  
  //cout << "Nhits = " << Nhits << endl;
  
  if(Nhits == 0) return;
  
  std::vector<IntersectionHit_t> ItersectionHitList_ppp;
  ItersectionHitList_ppp.clear();

  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ItersectionHitList[ihit].s > ItersectionHitList[max_ihit-1].s) continue;
    ItersectionHitList_ppp.push_back(ItersectionHitList[ihit]);
  }
  ItersectionHitList.clear();
  for(int ihit=0;ihit<int(ItersectionHitList_ppp.size());ihit++) ItersectionHitList.push_back(ItersectionHitList_ppp[ihit]);
    
  return;
  
}
//====================================================================
void  GGeometry::SetWorldVolume(GGeoObject* aGeoObject)
{
  
  //Set geometry's world volume
  
  if(WorldVolume != NULL) return;
  
  TString WorldVolName = TString("Geometry ") + GeometryName + TString(" World Volume");
  WorldVolume = aGeoObject->clone(WorldVolName);
  
  return;
  
}
//====================================================================
void  GGeometry::SetBField(GBField* aBfield) {
  
  //Set geometry's magnetic field
  
  if(Bfield != NULL) return;
  
  TString BfieldName = TString("Geometry ") + GeometryName + TString(" B-field");
  Bfield = aBfield->clone(BfieldName);
  
  return;
  
}
//====================================================================
int  GGeometry::FindTrackFinderAlgo(TVector3 x0, TVector3 p0)
{
  
  std::vector<int> index_list;
  index_list.clear();
  for(int i=0;i<int(TrackFinderAlgoList.size());i++) {
    if(global->IsInitCondsInTrackFindingRegion(TrackFinderAlgoList[i]->GetRegion(),x0,p0)) {
      index_list.push_back(i);
    }
  }
  
  if(index_list.size() > 1) {
    TVector3 x0_ppp = (1.0/global->GetUnit("cm"))*x0;
    TVector3 p0_ppp = (1.0/global->GetUnit("GeV/c"))*p0;
    
    cout << endl;
    cout << "ERROR on GGeometry::FindTrackFinderAlgo:: initial conditions x0 = (" << x0_ppp.X() << "," << x0_ppp.Y() << "," << x0_ppp.Z() << ") cm and p0 = ("
         << p0_ppp.X() << "," << p0_ppp.Y() << "," << p0_ppp.Z() << ") GeV/c are inside the regions of more than one TrackFinder algorithms for geometry " << GeometryName.Data() << " :"
         << endl;
    for(int i=0;i<int(index_list.size());i++) {
      TrackFinderAlgoList[index_list[i]]->Print();
      cout << endl;
    }
    cout << " Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  else if(index_list.size() == 1) return index_list[0];
  
  return  -1;
  
}
//====================================================================
void   GGeometry::CheckGeometry()
{
 
  //Check geometry elements overlaps
  
  int Ngeo_checks = GeometryElementsList.size()*(GeometryElementsList.size() - 1)/2;
  
  cout << endl;
  cout << "Begin of overlap check for geometry " << GeometryName.Data() << endl;
  cout << "  Geometry contains " << GeometryElementsList.size() << " geometry elements. Will perform " << Ngeo_checks << " overlaps checks!" << endl;
  
  int counter = 0;
  for(int i=0;i<int(GeometryElementsList.size());i++) {
    for(int j=0;j<int(GeometryElementsList.size());j++) {
      if(i >= j) continue;
      
      counter++;
      
      cout << "  Checking overlap between Geometry elements " << i+1 << " (" << GeometryElementsList[i]->GetName().Data() << ") and "
           << j+1 << "(" << GeometryElementsList[j]->GetName().Data() << "). Check number " << counter << " out of " << Ngeo_checks << "."
           << endl;

      if(DoGeometryElementsOverlap(GeometryElementsList[i],GeometryElementsList[j],global,GeoCheckPrecision)) {
	cout << "    Overlap found. Check your inputs. Exiting now!!!" << endl;
	assert(false);
      }
      
    }
  }
  cout << "  No overlaps found!!!" << endl;
  cout << "End of overlap check for geometry " << GeometryName.Data() << endl;
  cout << endl;
  
  return;
  
}
//====================================================================
void   GGeometry::Print()
{
  
  //Print geometry
  
  int counter = 0;
  
  cout << endl;
  cout << "Being of " << GeometryName.Data() << " geometry:" << endl;
  cout << endl;
  
  if(WorldVolume != NULL) {
    cout << "  Begin Wold Geometry" << endl;
    WorldVolume->Print();
    cout << "  End   Wold Geometry" << endl;
    cout << endl;
  }
  
  if(Bfield != NULL) {
    Bfield->Print();
    cout << endl;
  }
  
  //Beam-test configuration
  if(TelescopePlanesIndexList.size() > 0 && DUTPlanesIndexList.size() > 0) {
    cout << endl;
    cout << "  Begin Beam-test configuration" << endl;
    cout << "  List of Telescope planes:" << endl;
    for(int itel=0;itel<int(TelescopePlanesIndexList.size());itel++) {
      int idx = TelescopePlanesIndexList[itel];
      cout << "   - index = " << idx << ", "
           << "geoElm-Name = " << GeometryElementsList[idx]->GetName().Data() << ", "
	   << "layer-Name = "  << GeometryElementsList[idx]->GetLayerName().Data()
           << endl;
    }
    cout << "  List of DUT planes:" << endl;
    for(int itel=0;itel<int(DUTPlanesIndexList.size());itel++) {
      int idx = DUTPlanesIndexList[itel];
      cout << "   - index = " << idx << ", "
           << "geoElm-Name = " << GeometryElementsList[idx]->GetName().Data() << ", "
	   << "layer-Name = "  << GeometryElementsList[idx]->GetLayerName().Data()
           << endl;
    }
    cout << "  End   Beam-test configuration" << endl;
    cout << endl;
  }
  
  for(int irm=0;irm<int(ResolutionModelList.size());irm++) {
    cout << "  ResolutionModel Element " << irm+1 << endl;
    if(ResolutionModelList[irm] != NULL) {
      ResolutionModelList[irm]->Print();
      cout << endl;
    }
  }
  
  for(int irm=0;irm<int(EfficiencyModelList.size());irm++) {
    cout << "  Efficiency Element " << irm+1 << endl;
    if(EfficiencyModelList[irm] != NULL) {
      EfficiencyModelList[irm]->Print();
      cout << endl;
    }
  }
  
  for(int irm=0;irm<int(TrackFinderAlgoList.size());irm++) {
    cout << "  Track Finder Algorithm Element " << irm+1 << endl;
    if(TrackFinderAlgoList[irm] != NULL) {
      TrackFinderAlgoList[irm]->Print();
      cout << endl;
    }
  }
  
  for(int irm=0;irm<int(TrackCutsList.size());irm++) {
    cout << "  Begin Track cut" << endl;
    cout << "    SystemName    "  << TrackCutsList[irm].SystemName.Data() << endl;
      
    cout << endl;
    for(int icut=0;icut<int(TrackCutsList[irm].CutList.size());icut++) {
      cout << "    Begin Cut" << endl;
      cout << "      pRange        (" << TrackCutsList[irm].CutList[icut].pMin/global->GetMomentumUnit("GeV/c") << "," << TrackCutsList[irm].CutList[icut].pMax/global->GetMomentumUnit("GeV/c") << ") GeV/c" << endl;
      cout << "      phiRange      (" << TrackCutsList[irm].CutList[icut].phiMin/global->GetAngleUnit("deg")    << "," << TrackCutsList[irm].CutList[icut].phiMax/global->GetAngleUnit("deg")    << ") deg"   << endl;
      cout << "      thetaRange    (" << TrackCutsList[irm].CutList[icut].thetaMin/global->GetAngleUnit("deg")  << "," << TrackCutsList[irm].CutList[icut].thetaMax/global->GetAngleUnit("deg")  << ") deg"   << endl;
      cout << "      NhitsRange    (" << TrackCutsList[irm].CutList[icut].NhitsMin                              << "," << TrackCutsList[irm].CutList[icut].NhitsMax                              << ")"       << endl;
      cout << "    End   Cut" << endl;
    }
    cout << endl;
    cout << "  End   Track cut" << endl;
    cout << endl;
  }
  
  counter = 0;
  for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) {
    if(GeometryElementsList[igeoElm] == NULL) continue;
    cout << "  Geometry Element " << counter+1 << endl;
    GeometryElementsList[igeoElm]->Print();
    cout << endl;
    counter++;
  }
  cout << endl;
  
  cout << "End   of " << GeometryName.Data() << " geometry:" << endl;
  cout << endl;
  
  return;
  
}
//====================================================================
void  GGeometry::PrintWeight()
{
  
  //Print geometry weight
  
  double Total_geometry_weight        = 0.0;
  double NonAsignedGeoElements_weight = 0.0;
  std::vector<double>  SystemWeightList;
  SystemWeightList.clear();
  for(int isys=0;isys<int(SystemNamesList.size());isys++) SystemWeightList.push_back(0.0);
  
  for(int igeoElm=0;igeoElm<int(GeometryElementsList.size());igeoElm++) { //begin loop over geo-elements
    if(GeometryElementsList[igeoElm]->GetSystemName() == TString("")) {
      NonAsignedGeoElements_weight += GeometryElementsList[igeoElm]->GetMass();
      continue;
    }
    
    for(int isys=0;isys<int(SystemNamesList.size());isys++) { //begin over system names
      if(SystemNamesList[isys] == GeometryElementsList[igeoElm]->GetSystemName()) {
	SystemWeightList[isys] += GeometryElementsList[igeoElm]->GetMass();
	break;
      }
    } //end of loop over system names
    
  } //end of loop over geo-elements
  
  cout << endl;
  cout << endl;
  cout << "=====================================================" << endl;
  cout << "       Weight of geometry " << GeometryName.Data() << endl;
  cout << "-----------------------------------------------------" << endl;
  if(NonAsignedGeoElements_weight > 1.0e-6*global->GetUnit("gr")) {
    Total_geometry_weight += NonAsignedGeoElements_weight;
    cout << " - total weight of non-asigned geoElements (System == \"\") = " << NonAsignedGeoElements_weight/global->GetUnit("Kg") << " Kg" << endl;
  }
  for(int isys=0;isys<int(SystemNamesList.size());isys++) {
    Total_geometry_weight += SystemWeightList[isys];
    cout << " - total " << SystemNamesList[isys].Data() << " system weight = " << SystemWeightList[isys]/global->GetUnit("Kg") << " Kg" << endl;
  }
  cout << "-----------------------------------------------------" << endl;
  cout << " - total geometry weight = " << Total_geometry_weight/global->GetUnit("Kg") << " Kg" << endl;
  cout << "=====================================================" << endl;
  cout << endl;
  cout << endl;
  
  return;
  
}
//====================================================================
