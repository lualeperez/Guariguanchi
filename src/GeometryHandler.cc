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
#include <TImage.h>
#include <TLatex.h>
#include "include/Guariguanchi.h"

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

bool Plot_TPCGAS_Material = false;
//bool Plot_TPCGAS_Material = true;

//====================================================================
void Guariguanchi::FillGeometries(void)
{
  
  // Fill geometry list from geometry datacard list
  
  GeometryList.clear();
  for(int k=0;k<int(GeometryDataCardList.size());k++) FillGeometryFromDataCard(GeometryDataCardList[k].Data(),GeometryList.size());

  bool Voxeled = false;
  std::vector<TString> MaterialList;
  MaterialList.clear();
  for(int igeo=0;igeo<int(GeometryList.size());igeo++) { //geometries loop
    GeometryList[igeo]->SetBField(GlobalBfield);                  // If geometry B-field is not specified it will be set to the GlobalBfield
    GeometryList[igeo]->SetGeoCheckPrecision(GeoCheckPrecision);  // Set geometry check (overlaps) precision to parameter
    GeometryList[igeo]->FillVoxeledGeoElementsList(VoxelList);    // Fill the geometry's list of geometry-elements inside the specified voxels
    GeometryList[igeo]->FillResolutionModelIndexes();             // Fill the geometry's resolution model indexes
    GeometryList[igeo]->FillEfficiencyModelIndexes();             // Fill the geometry's efficiency model indexes
    GeometryList[igeo]->FillBeamTestConfigPlanesIndexes();        // Fill the list of telescope and DUT planes. Only in case of BeamTelescope analysis
    GeometryList[igeo]->FillSystemsList();                        // Fill the list with the geometry systems
    GeometryList[igeo]->ApplyBkgScaling();                        // Apply global bkg scaling if it is different from zero
    
    if(!Voxeled && GeometryList[igeo]->GetNGeoElements() != GeometryList[igeo]->GetNVoxelesGeoElements()) Voxeled = true;
    
    //Filling list with materials in geometries
    for(int igeoElement=0;igeoElement<GeometryList[igeo]->GetNGeoElements();igeoElement++) { // loop over geometry elements
      TString material = GeometryList[igeo]->GetGeometryElement(igeoElement)->GetMaterial();
      bool IsInList = false;
      for(int k=0;k<int(MaterialList.size());k++) {
	if(MaterialList[k] == material) {
	  IsInList = true;
	  break;
	}
      }
      if(!IsInList) MaterialList.push_back(material);
    } // end of loop over geometry elements
    
  } //end of geometry loop
  
  for(int igeo=0;igeo<int(GeometryList.size());igeo++) {
    if(GeometryList[igeo]->GetNGeoElements() == 0) {
      cout << endl;
      cout << "ERROR Guariguanchi::FillGeometries:: Geometry " << GeometryList[igeo]->GetName().Data() << " has zero geoElements. Check your inputs. Exiting now!!!" << endl;
      cout << endl;
      assert(false);
    }
    if(GeometryList[igeo]->GetNVoxelesGeoElements() == 0) {
      cout << endl;
      cout << "ERROR Guariguanchi::FillGeometries:: Geometry " << GeometryList[igeo]->GetName().Data() << " has zero voxeled geoElements. Check your inputs. Exiting now!!!" << endl;
      cout << endl;
      assert(false);
    }
  }
  
  if(Voxeled) {
    cout << endl;
    cout << "Voxel lists:" << endl;
    for(int igeo=0;igeo<int(GeometryList.size());igeo++) {
      cout << "Geo: " << GeometryList[igeo]->GetName().Data() << " has " << GeometryList[igeo]->GetNVoxelesGeoElements() << " voxeled GeoElements ( out of " << GeometryList[igeo]->GetNGeoElements() << ")!!!" << endl;
    }
    cout << endl;
  }
  
  cout << endl;
  cout << "Material list:" << endl;
  for(int k=0;k<int(MaterialList.size());k++) {
    cout << " - material = " << MaterialList[k].Data() << endl;
  }
  cout << endl;

  return;
  
}
//====================================================================
void  Guariguanchi::FillTrackers(void)
{
  
  // Fill the tracker list
  // Each geometry has a corresponding tracker for geometry navigation and tracking related calculations

  for(int igeo=0;igeo<int(GeometryList.size());igeo++) { //geometries loop
    TString  Name = TString("Tracker for geometry") + GeometryList[igeo]->GetName();
    TVector3 InitMomentum = global->GetMomentum(global->GetMomentumFromMomVar(momArr[0],particle,momVariable),thetaArr[0],phiArr[0]);
    GTracker* aTracker = new GTracker(Name,
				      GeometryList[igeo],
				      particle,
				      ParticleOrigin,
				      InitMomentum,
				      ReferencePoint,
				      global);
    if(IncludeEloss) aTracker->SetIncludeEloss(IncludeEloss);
    
    TrackerList.push_back(aTracker);
  }  //end of geometry loop

  return;
  
}
//====================================================================
void  Guariguanchi::CheckGeometries(void)
{

  for(int igeo=0;igeo<int(GeometryList.size());igeo++) { //Begin of loop over geometries
    GeometryList[igeo]->CheckGeometry();
  } //End of loop over geometries

  return;
  
}
//====================================================================
void Guariguanchi::PrintGeometries(void)
{

  cout << endl;
  cout << "Printing out the specified geometries:" << endl;
  for(int igeo=0;igeo<int(GeometryList.size());igeo++) { //geometries loop
    cout << endl;
    GeometryList[igeo]->Print();
    cout << endl;
  } //end of geometry loop
  cout << endl;

  return;
  
}
//====================================================================
void  Guariguanchi::PrintGeometryWeights(void)
{
  
  cout << endl;
  cout << "Printing out geometries weight:" << endl;
  for(int igeo=0;igeo<int(GeometryList.size());igeo++) { //geometries loop
    cout << endl;
    GeometryList[igeo]->PrintWeight();
    cout << endl;
  } //end of geometry loop
  cout << endl;
  
  return;
  
}
//====================================================================
void  Guariguanchi::CheckGeometriesComparability(void)
{
  
  //This function checks if the geometries can be compared
  //in terms of tracking performances => same number and type of track parameters
  
  for(int itrker1=0;itrker1<int(TrackerList.size());itrker1++) {
    int Nfitparams1 = TrackerList[itrker1]->GetTrajectory()->GetNParameters();
    for(int itrker2=0;itrker2<int(TrackerList.size());itrker2++) {
      if(itrker1 >= itrker2) continue;
      int Nfitparams2 = TrackerList[itrker2]->GetTrajectory()->GetNParameters();
            
      if(Nfitparams1 != Nfitparams2) {
	cout << endl;
	cout << "Tracker for geometries " << GeometryList[itrker1]->GetName().Data() << " and " << GeometryList[itrker2]->GetName().Data() << " have different number of track parameters." << endl;
	cout << " ==> Tracking performances cannot be compared." << endl;
	cout << endl;
	assert(false);
      }
      
      for(int ipar=0;ipar<Nfitparams1;ipar++) {
	if(TrackerList[itrker1]->GetTrajectory()->GetParameterName(ipar) != TrackerList[itrker2]->GetTrajectory()->GetParameterName(ipar)) {
	  cout << endl;
	  cout << "Tracker for geometries " << GeometryList[itrker1]->GetName().Data() << " and " << GeometryList[itrker2]->GetName().Data() << " have different set of track parameters." << endl;
	  cout << " ==> Tracking performances cannot be compared." << endl;
	  cout << endl;
	  assert(false);
	}
      }
      
    }
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::CheckTelescopeConfigurations(void)
{
  
  //Check if the beam-test configurations has been well defined for all geometries
  
  if(!DoTelescopeAnalysis) return;
  
  for(int igeo=0;igeo<int(GeometryList.size());igeo++) { //Begin of loop over geometries
    int Ntelescope_planes = GeometryList[igeo]->GetNTelescopePlanes();
    int NDUT_planes       = GeometryList[igeo]->GetNDUTPlanes();
    int MinHits           = TrackerList[igeo]->GetTrajectory()->GetMinHits();
    
    if(Ntelescope_planes < MinHits && NDUT_planes < 1) {
      cout << endl;
      cout << "Gemetry " << GeometryList[igeo]->GetName().Data() 
           << " has a number of telescope (" << Ntelescope_planes << ") and DUT (" << NDUT_planes << ") planes which are smaller than the required ones " << MinHits << " and 1, respectively!!!"
	   << endl;
      cout << "Check your inputs. Exiting now!!!" << endl;
      cout << endl;
      assert(false);
    }    
  } //End of loop over geometries
  
  return;
  
}
//====================================================================
void Guariguanchi::PlotGeometries(void)
{

  GlobalFileCounter++;
  TString HistName,HistTitle;
  
  double RX[2];
  double RY[2];
  double RZ[2];
  double RR[2];
  
  const int NGeometries(GeometryList.size());
  const int NthetaPoints(thetaArr.size());
  const int NphiPoints(phiArr.size());
    
  std::vector<TGraph>*   GraphXY = new std::vector<TGraph>[NGeometries];
  std::vector<TGraph>*   GraphZY = new std::vector<TGraph>[NGeometries];
  std::vector<TGraph>*   GraphZX = new std::vector<TGraph>[NGeometries];
  std::vector<TGraph>*   GraphZR = new std::vector<TGraph>[NGeometries];
  
  RX[0] = +1.0e+20;
  RX[1] = -1.0e+20;
  RY[0] = +1.0e+20;
  RY[1] = -1.0e+20;
  RZ[0] = +1.0e+20;
  RZ[1] = -1.0e+20;
  RR[0] = +1.0e+20;
  RR[1] = -1.0e+20;

  for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
    GraphXY[igeo].clear();
    GraphZY[igeo].clear();
    GraphZX[igeo].clear();
    GraphZR[igeo].clear();
    
    for(int igeoElm=0;igeoElm<GeometryList[igeo]->GetNGeoElements();igeoElm++) { //Begin loop over geometry elements
      GGeoObject* aGeoObject = GeometryList[igeo]->GetGeometryElement(igeoElm);
    //for(int igeoElm=0;igeoElm<GeometryList[igeo]->GetNVoxelesGeoElements();igeoElm++) { //Begin loop over geometry elements
      //GGeoObject* aGeoObject = GeometryList[igeo]->GetVoxeledGeometryElement(igeoElm);
      
      int color = 0;
      TString material = aGeoObject->GetMaterial();
      color = global->GetMaterialColor(material);
      
      bool PlotThisOne = true;
      if(!Plot_TPCGAS_Material) {
	if(material == TString("TPCGAS")) PlotThisOne = false;
      }
      
      for(int iboundary=0;iboundary<aGeoObject->GetNBoundarySurfaces();iboundary++) {
	GSurfaceObject* aSurface = aGeoObject->GetBoundarySurface(iboundary);
	if(PlotThisOne) aSurface->FillSurfaceRepresentation(color,GraphXY[igeo],GraphZY[igeo],GraphZX[igeo],GraphZR[igeo]);
      }

    } // end of loop over geometry elements
    
    if(DoPlotWorldVolume) {
      //visualising the world volume
      GGeoObject* aGeoObject =  GeometryList[igeo]->GetWorldVolume();
      
      int color = kRed;
      int style = 2;
      for(int iboundary=0;iboundary<aGeoObject->GetNBoundarySurfaces();iboundary++) {
	GSurfaceObject* aSurface = aGeoObject->GetBoundarySurface(iboundary);
	aSurface->FillSurfaceRepresentation(color,GraphXY[igeo],GraphZY[igeo],GraphZX[igeo],GraphZR[igeo],style);
      }
    }
    
    if(DoPlotStepBfieldVolume) {
      //visualising Step B-field volume
      if(GeometryList[igeo]->GetBField()->GetType() == TString("MultipleSteps")) {
	GBFieldMultipleSteps* aBfield_tmp = dynamic_cast<GBFieldMultipleSteps*>(GeometryList[igeo]->GetBField());
	
	int color = kBlue;
        int style = 2;
	
	for(int ivol=0;ivol<aBfield_tmp->GetNVolumes();ivol++) {
	  GGeoObject* aGeoObject = aBfield_tmp->GetVolume(ivol);
          for(int iboundary=0;iboundary<aGeoObject->GetNBoundarySurfaces();iboundary++) {
	    GSurfaceObject* aSurface = aGeoObject->GetBoundarySurface(iboundary);
	    aSurface->FillSurfaceRepresentation(color,GraphXY[igeo],GraphZY[igeo],GraphZX[igeo],GraphZR[igeo],style);
          }
	}
      }
    }

    for(int i=0;i<int(GraphXY[igeo].size());i++) {
      for(int j=0;j<GraphXY[igeo][i].GetN();j++) {
	double x,y;
	GraphXY[igeo][i].GetPoint(j,x,y);
	if(RX[0] > x) RX[0] = x;
	if(RX[1] < x) RX[1] = x;
	if(RY[0] > y) RY[0] = y;
	if(RY[1] < y) RY[1] = y;
      }
    }
    
    for(int i=0;i<int(GraphZY[igeo].size());i++) {
      for(int j=0;j<GraphZY[igeo][i].GetN();j++) {
	double z,y;
	GraphZY[igeo][i].GetPoint(j,z,y);
	if(RZ[0] > z) RZ[0] = z;
	if(RZ[1] < z) RZ[1] = z;
	if(RY[0] > y) RY[0] = y;
	if(RY[1] < y) RY[1] = y;
      }
    }
    
    for(int i=0;i<int(GraphZX[igeo].size());i++) {
      for(int j=0;j<GraphZX[igeo][i].GetN();j++) {
	double z,x;
	GraphZX[igeo][i].GetPoint(j,z,x);
	if(RZ[0] > z) RZ[0] = z;
	if(RZ[1] < z) RZ[1] = z;
	if(RX[0] > x) RX[0] = x;
	if(RX[1] < x) RX[1] = x;
      }
    }
    
    for(int i=0;i<int(GraphZR[igeo].size());i++) {
      for(int j=0;j<GraphZR[igeo][i].GetN();j++) {
	double z,r;
	GraphZR[igeo][i].GetPoint(j,z,r);
	if(RZ[0] > z) RZ[0] = z;
	if(RZ[1] < z) RZ[1] = z;
	if(RR[0] > r) RR[0] = r;
	if(RR[1] < r) RR[1] = r;
      }
    }
  } // End of loop over geometries
  
  double xref = ParticleOrigin.X()/global->GetDistanceUnit("cm");
  double yref = ParticleOrigin.Y()/global->GetDistanceUnit("cm");
  double zref = ParticleOrigin.Z()/global->GetDistanceUnit("cm");
  double rref = sqrt(pow(xref,2) + pow(yref,2));
  if(RX[0] > xref) RX[0] = xref;
  if(RX[1] < xref) RX[1] = xref;
  if(RY[0] > yref) RY[0] = yref;
  if(RY[1] < yref) RY[1] = yref;
  if(RZ[0] > zref) RZ[0] = zref;
  if(RZ[1] < zref) RZ[1] = zref;
  if(RR[0] > rref) RR[0] = rref;
  if(RR[1] < rref) RR[1] = rref;
  
  xref = ReferencePoint.X()/global->GetDistanceUnit("cm");
  yref = ReferencePoint.Y()/global->GetDistanceUnit("cm");
  zref = ReferencePoint.Z()/global->GetDistanceUnit("cm");
  rref = sqrt(pow(xref,2) + pow(yref,2));
  if(RX[0] > xref) RX[0] = xref;
  if(RX[1] < xref) RX[1] = xref;
  if(RY[0] > yref) RY[0] = yref;
  if(RY[1] < yref) RY[1] = yref;
  if(RZ[0] > zref) RZ[0] = zref;
  if(RZ[1] < zref) RZ[1] = zref;
  if(RR[0] > rref) RR[0] = rref;
  if(RR[1] < rref) RR[1] = rref;

  double delta,porcent;
  delta   = RX[1] - RX[0];
  porcent = 0.10;
  RX[0] -= delta*porcent;
  RX[1] += delta*porcent;
  delta   = RY[1] - RY[0];
  porcent = 0.10;
  RY[0] -= delta*porcent;
  RY[1] += delta*porcent;
  delta   = RZ[1] - RZ[0];
  porcent = 0.10;
  RZ[0] -= delta*porcent;
  porcent = 0.10;
  RZ[1] += delta*porcent;
  delta   = RR[1] - RR[0];
  porcent = 0.10;
  RR[0] -= delta*porcent;
  porcent = 0.10;
  RR[1] += delta*porcent;


  int Nbins_ref = 500;
  TH2F hrefZY("hrefZY","",Nbins_ref,RZ[0],RZ[1],Nbins_ref,RY[0],RY[1]);
  hrefZY.SetXTitle("Z (cm)");
  hrefZY.GetXaxis()->CenterTitle(true);
  hrefZY.SetYTitle("Y (cm)");
  hrefZY.GetYaxis()->CenterTitle(true);
  hrefZY.SetLineColor(1);
  hrefZY.SetLineWidth(1);
  hrefZY.GetXaxis()->SetTitleSize(TheSize);
  hrefZY.GetXaxis()->SetLabelSize(TheSize);
  hrefZY.GetYaxis()->SetTitleSize(TheSize);
  hrefZY.GetYaxis()->SetLabelSize(TheSize);
  hrefZY.SetStats(false);
  
  TH2F hrefZX("hrefZX","",Nbins_ref,RZ[0],RZ[1],Nbins_ref,RX[0],RX[1]);
  hrefZX.SetXTitle("Z (cm)");
  hrefZX.GetXaxis()->CenterTitle(true);
  hrefZX.SetYTitle("X (cm)");
  hrefZX.GetYaxis()->CenterTitle(true);
  hrefZX.SetLineColor(1);
  hrefZX.SetLineWidth(1);
  hrefZX.GetXaxis()->SetTitleSize(TheSize);
  hrefZX.GetXaxis()->SetLabelSize(TheSize);
  hrefZX.GetYaxis()->SetTitleSize(TheSize);
  hrefZX.GetYaxis()->SetLabelSize(TheSize);
  hrefZX.SetStats(false);
  
  TH2F hrefXY("hrefXY","",Nbins_ref,RX[0],RX[1],Nbins_ref,RY[0],RY[1]);
  hrefXY.SetXTitle("X (cm)");
  hrefXY.GetXaxis()->CenterTitle(true);
  hrefXY.SetYTitle("Y (cm)");
  hrefXY.GetYaxis()->CenterTitle(true);
  hrefXY.SetLineColor(1);
  hrefXY.SetLineWidth(1);
  hrefXY.GetXaxis()->SetTitleSize(TheSize);
  hrefXY.GetXaxis()->SetLabelSize(TheSize);
  hrefXY.GetYaxis()->SetTitleSize(TheSize);
  hrefXY.GetYaxis()->SetLabelSize(TheSize);
  hrefXY.SetStats(false);
  
  TH2F hrefZR("hrefZR","",Nbins_ref,RZ[0],RZ[1],Nbins_ref,RR[0],RR[1]);
  hrefZR.SetXTitle("Z (cm)");
  hrefZR.GetXaxis()->CenterTitle(true);
  hrefZR.SetYTitle("R (cm)");
  hrefZR.GetYaxis()->CenterTitle(true);
  hrefZR.SetLineColor(1);
  hrefZR.SetLineWidth(1);
  hrefZR.GetXaxis()->SetTitleSize(TheSize);
  hrefZR.GetXaxis()->SetLabelSize(TheSize);
  hrefZR.GetYaxis()->SetTitleSize(TheSize);
  hrefZR.GetYaxis()->SetLabelSize(TheSize);
  hrefZR.SetStats(false);

  TH2F** hrefZY_geoRange;
  TH2F** hrefZX_geoRange;
  TH2F** hrefXY_geoRange;
  TH2F** hrefZR_geoRange;
  const int NgeoRanges(GeoRanges.size());
  if(NgeoRanges > 0) {
    hrefZY_geoRange = new TH2F*[NgeoRanges];
    hrefZX_geoRange = new TH2F*[NgeoRanges];
    hrefXY_geoRange = new TH2F*[NgeoRanges];
    hrefZR_geoRange = new TH2F*[NgeoRanges];
    
    for(int irange=0;irange<NgeoRanges;irange++) {
      double RX_geoRange[2];
      double RY_geoRange[2];
      double RZ_geoRange[2];
      double RR_geoRange[2];
      RX_geoRange[0] = GeoRanges[irange].Xmin/global->GetDistanceUnit("cm");
      RX_geoRange[1] = GeoRanges[irange].Xmax/global->GetDistanceUnit("cm");
      RY_geoRange[0] = GeoRanges[irange].Ymin/global->GetDistanceUnit("cm");
      RY_geoRange[1] = GeoRanges[irange].Ymax/global->GetDistanceUnit("cm");
      RZ_geoRange[0] = GeoRanges[irange].Zmin/global->GetDistanceUnit("cm");
      RZ_geoRange[1] = GeoRanges[irange].Zmax/global->GetDistanceUnit("cm");
      
      RR_geoRange[0] = +1.0e+20;
      RR_geoRange[1] = -1.0e+20;
      if(RX_geoRange[0] < 0 || RX_geoRange[1] < 0 || RY_geoRange[0] < 0 || RY_geoRange[1] < 0) RR_geoRange[0] = 0.0;
      RR_geoRange[0] = TMath::Min(RR_geoRange[0],sqrt(pow(RX_geoRange[0],2) + pow(RY_geoRange[0],2)));
      RR_geoRange[0] = TMath::Min(RR_geoRange[0],sqrt(pow(RX_geoRange[0],2) + pow(RY_geoRange[1],2)));
      RR_geoRange[0] = TMath::Min(RR_geoRange[0],sqrt(pow(RX_geoRange[1],2) + pow(RY_geoRange[0],2)));
      RR_geoRange[0] = TMath::Min(RR_geoRange[0],sqrt(pow(RX_geoRange[1],2) + pow(RY_geoRange[0],2)));
      RR_geoRange[1] = TMath::Max(RR_geoRange[1],sqrt(pow(RX_geoRange[0],2) + pow(RY_geoRange[0],2)));
      RR_geoRange[1] = TMath::Max(RR_geoRange[1],sqrt(pow(RX_geoRange[0],2) + pow(RY_geoRange[1],2)));
      RR_geoRange[1] = TMath::Max(RR_geoRange[1],sqrt(pow(RX_geoRange[1],2) + pow(RY_geoRange[0],2)));
      RR_geoRange[1] = TMath::Max(RR_geoRange[1],sqrt(pow(RX_geoRange[1],2) + pow(RY_geoRange[0],2)));
      
      HistName = TString("hrefZY_geoRange") + long(irange+1);
      hrefZY_geoRange[irange] = new TH2F(HistName.Data(),"",
					 100,RZ_geoRange[0],RZ_geoRange[1],
					 100,RY_geoRange[0],RY_geoRange[1]);
      hrefZY_geoRange[irange]->SetXTitle("Z (cm)");
      hrefZY_geoRange[irange]->GetXaxis()->CenterTitle(true);
      hrefZY_geoRange[irange]->SetYTitle("Y (cm)");
      hrefZY_geoRange[irange]->GetYaxis()->CenterTitle(true);
      hrefZY_geoRange[irange]->SetLineColor(1);
      hrefZY_geoRange[irange]->SetLineWidth(1);
      hrefZY_geoRange[irange]->GetXaxis()->SetTitleSize(TheSize);
      hrefZY_geoRange[irange]->GetXaxis()->SetLabelSize(TheSize);
      hrefZY_geoRange[irange]->GetYaxis()->SetTitleSize(TheSize);
      hrefZY_geoRange[irange]->GetYaxis()->SetLabelSize(TheSize);
      hrefZY_geoRange[irange]->SetStats(false);
      
      HistName = TString("hrefZX_geoRange") + long(irange+1);
      hrefZX_geoRange[irange] = new TH2F(HistName.Data(),"",
					 100,RZ_geoRange[0],RZ_geoRange[1],
					 100,RX_geoRange[0],RX_geoRange[1]);
      hrefZX_geoRange[irange]->SetXTitle("Z (cm)");
      hrefZX_geoRange[irange]->GetXaxis()->CenterTitle(true);
      hrefZX_geoRange[irange]->SetYTitle("X (cm)");
      hrefZX_geoRange[irange]->GetYaxis()->CenterTitle(true);
      hrefZX_geoRange[irange]->SetLineColor(1);
      hrefZX_geoRange[irange]->SetLineWidth(1);
      hrefZX_geoRange[irange]->GetXaxis()->SetTitleSize(TheSize);
      hrefZX_geoRange[irange]->GetXaxis()->SetLabelSize(TheSize);
      hrefZX_geoRange[irange]->GetYaxis()->SetTitleSize(TheSize);
      hrefZX_geoRange[irange]->GetYaxis()->SetLabelSize(TheSize);
      hrefZX_geoRange[irange]->SetStats(false);
      
      HistName = TString("hrefXY_geoRange") + long(irange+1);
      hrefXY_geoRange[irange] = new TH2F(HistName.Data(),"",
					 100,RX_geoRange[0],RX_geoRange[1],
					 100,RY_geoRange[0],RY_geoRange[1]);
      hrefXY_geoRange[irange]->SetXTitle("X (cm)");
      hrefXY_geoRange[irange]->GetXaxis()->CenterTitle(true);
      hrefXY_geoRange[irange]->SetYTitle("Y (cm)");
      hrefXY_geoRange[irange]->GetYaxis()->CenterTitle(true);
      hrefXY_geoRange[irange]->SetLineColor(1);
      hrefXY_geoRange[irange]->SetLineWidth(1);
      hrefXY_geoRange[irange]->GetXaxis()->SetTitleSize(TheSize);
      hrefXY_geoRange[irange]->GetXaxis()->SetLabelSize(TheSize);
      hrefXY_geoRange[irange]->GetYaxis()->SetTitleSize(TheSize);
      hrefXY_geoRange[irange]->GetYaxis()->SetLabelSize(TheSize);
      hrefXY_geoRange[irange]->SetStats(false);
      
      HistName = TString("hrefZR_geoRange") + long(irange+1);
      hrefZR_geoRange[irange] = new TH2F(HistName.Data(),"",
					 100,RZ_geoRange[0],RZ_geoRange[1],
					 100,RR_geoRange[0],RR_geoRange[1]);
      hrefZR_geoRange[irange]->SetXTitle("Z (cm)");
      hrefZR_geoRange[irange]->GetXaxis()->CenterTitle(true);
      hrefZR_geoRange[irange]->SetYTitle("R (cm)");
      hrefZR_geoRange[irange]->GetYaxis()->CenterTitle(true);
      hrefZR_geoRange[irange]->SetLineColor(1);
      hrefZR_geoRange[irange]->SetLineWidth(1);
      hrefZR_geoRange[irange]->GetXaxis()->SetTitleSize(TheSize);
      hrefZR_geoRange[irange]->GetXaxis()->SetLabelSize(TheSize);
      hrefZR_geoRange[irange]->GetYaxis()->SetTitleSize(TheSize);
      hrefZR_geoRange[irange]->GetYaxis()->SetLabelSize(TheSize);
      hrefZR_geoRange[irange]->SetStats(false);
    }
  } 
  
  int ref_line_color = 1;
  TGraph* lhZY = new TGraph();
  lhZY->SetPoint(0,hrefZY.GetXaxis()->GetXmin(),0.0);
  lhZY->SetPoint(1,hrefZY.GetXaxis()->GetXmax(),0.0);
  lhZY->SetLineColor(ref_line_color);
  lhZY->SetMarkerColor(ref_line_color);
  lhZY->SetLineWidth(1);
  lhZY->SetLineStyle(2);
  TGraph* lvZY = new TGraph();
  lvZY->SetPoint(0,0.0,hrefZY.GetYaxis()->GetXmin());
  lvZY->SetPoint(1,0.0,hrefZY.GetYaxis()->GetXmax());
  lvZY->SetLineColor(ref_line_color);
  lvZY->SetMarkerColor(ref_line_color);
  lvZY->SetLineWidth(1);
  lvZY->SetLineStyle(2);
  
  TGraph* lhZX = new TGraph();
  lhZX->SetPoint(0,hrefZX.GetXaxis()->GetXmin(),0.0);
  lhZX->SetPoint(1,hrefZX.GetXaxis()->GetXmax(),0.0);
  lhZX->SetLineColor(ref_line_color);
  lhZX->SetMarkerColor(ref_line_color);
  lhZX->SetLineWidth(1);
  lhZX->SetLineStyle(2);
  TGraph* lvZX = new TGraph();
  lvZX->SetPoint(0,0.0,hrefZX.GetYaxis()->GetXmin());
  lvZX->SetPoint(1,0.0,hrefZX.GetYaxis()->GetXmax());
  lvZX->SetLineColor(ref_line_color);
  lvZX->SetMarkerColor(ref_line_color);
  lvZX->SetLineWidth(1);
  lvZX->SetLineStyle(2);
  
  TGraph* lhXY = new TGraph();
  lhXY->SetPoint(0,hrefXY.GetXaxis()->GetXmin(),0.0);
  lhXY->SetPoint(1,hrefXY.GetXaxis()->GetXmax(),0.0);
  lhXY->SetLineColor(ref_line_color);
  lhXY->SetMarkerColor(ref_line_color);
  lhXY->SetLineWidth(1);
  lhXY->SetLineStyle(2);
  TGraph* lvXY = new TGraph();
  lvXY->SetPoint(0,0.0,hrefXY.GetYaxis()->GetXmin());
  lvXY->SetPoint(1,0.0,hrefXY.GetYaxis()->GetXmax());
  lvXY->SetLineColor(ref_line_color);
  lvXY->SetMarkerColor(ref_line_color);
  lvXY->SetLineWidth(1);
  lvXY->SetLineStyle(2);
  
  TGraph* lhZR = new TGraph();
  lhZR->SetPoint(0,hrefZR.GetXaxis()->GetXmin(),0.0);
  lhZR->SetPoint(1,hrefZR.GetXaxis()->GetXmax(),0.0);
  lhZR->SetLineColor(ref_line_color);
  lhZR->SetMarkerColor(ref_line_color);
  lhZR->SetLineWidth(1);
  lhZR->SetLineStyle(2);
  TGraph* lvZR = new TGraph();
  lvZR->SetPoint(0,0.0,hrefZR.GetYaxis()->GetXmin());
  lvZR->SetPoint(1,0.0,hrefZR.GetYaxis()->GetXmax());
  lvZR->SetLineColor(ref_line_color);
  lvZR->SetMarkerColor(ref_line_color);
  lvZR->SetLineWidth(1);
  lvZR->SetLineStyle(2);
  

  TString command;
  TString EPSName  = TheOutputFile + TString("_") + long(GlobalFileCounter) + TString(".eps");
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");

  TCanvas* c1 = new TCanvas("c1","c1",LargeCanvasX,LargeCanvasY);
  //TCanvas* c1 = new TCanvas("c1","c1");
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  //c1->SetRightMargin(0.15);

  c1->Print(EPSNameO.Data());
  
  for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZY");
    hrefZY.SetTitle(HistName.Data());
    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZX");
    hrefZX.SetTitle(HistName.Data());
    HistName = GeometryList[igeo]->GetName() + TString(" geometry XY");
    hrefXY.SetTitle(HistName.Data());
    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZR");
    hrefZR.SetTitle(HistName.Data());

    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hrefZY.Draw();
    for(int iplane=0;iplane<int(GraphZY[igeo].size());iplane++) {
      if(GraphZY[igeo][iplane].GetN() > 0) GraphZY[igeo][iplane].Draw("PEL");
    }
    lhZY->Draw("PEL");
    lvZY->Draw("PEL");
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hrefZX.Draw();
    for(int iplane=0;iplane<int(GraphZX[igeo].size());iplane++) {
      if(GraphZX[igeo][iplane].GetN() > 0) GraphZX[igeo][iplane].Draw("PEL");
    }
    lhZX->Draw("PEL");
    lvZX->Draw("PEL");
    c1->cd(3);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hrefXY.Draw();
    for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
      if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
    }
    lhXY->Draw("PEL");
    lvXY->Draw("PEL");
    c1->cd();
    PlotLogo(0.15,0.89,0.51,0.5);
    c1->Print(EPSName.Data());

#if 0
    TVector3  MyOrigin(0.0,0.0,0.0);
    std::vector<TVector3> MyMomentum;
    std::vector<TString>  MyParticle;
    std::vector<int>      MyColors;
    MyMomentum.clear();
    MyParticle.clear();
    MyColors.clear();
    double alpha;
    
    alpha = 60.0*global->GetUnit("deg");
    MyMomentum.push_back( (200*global->GetMomentumUnit("MeV/c"))*TVector3(TMath::Cos(alpha),TMath::Sin(alpha),0.0) );
    MyParticle.push_back(TString("mu+"));
    MyColors.push_back(kBlue);
    
    alpha = 120.0*global->GetUnit("deg");
    MyMomentum.push_back( (200*global->GetMomentumUnit("MeV/c"))*TVector3(TMath::Cos(alpha),TMath::Sin(alpha),0.0) );
    MyParticle.push_back(TString("mu-"));
    MyColors.push_back(kRed);
    
    alpha = 240.0*global->GetUnit("deg");
    MyMomentum.push_back( (200*global->GetMomentumUnit("MeV/c"))*TVector3(TMath::Cos(alpha),TMath::Sin(alpha),0.0) );
    MyParticle.push_back(TString("mu+"));
    MyColors.push_back(kGreen+2);
    
    alpha = -60.0*global->GetUnit("deg");
    MyMomentum.push_back( (200*global->GetMomentumUnit("MeV/c"))*TVector3(TMath::Cos(alpha),TMath::Sin(alpha),0.0) );
    MyParticle.push_back(TString("mu-"));
    MyColors.push_back(kMagenta);
    
    const int Nmoms(MyMomentum.size());
    TGraph* gMytrackXY[Nmoms];
    TGraph* gMytrackXY_points[Nmoms];
    for(int i=0;i<Nmoms;i++) {
      double MarkerSize  = 0.8;
      int    MarkerStyle = 20;
      
      HistName   = TString("gMytrackXY_points_mon") + long(i+1);
      HistTitle  = TString("Track geometry intersections projection in XY plane for momentum ") + long(i+1);
      gMytrackXY_points[i] = new TGraph();
      gMytrackXY_points[i]->SetName(HistName.Data());
      gMytrackXY_points[i]->SetTitle(HistTitle.Data());
      gMytrackXY_points[i]->SetMarkerColor(MyColors[i]);
      gMytrackXY_points[i]->SetMarkerSize(MarkerSize);
      gMytrackXY_points[i]->SetMarkerStyle(MarkerStyle);
      gMytrackXY_points[i]->SetLineColor(MyColors[i]);
      gMytrackXY_points[i]->SetLineWidth(2);
      
      HistName   = TString("gMytrackXY_mon") + long(i+1);
      HistTitle  = TString("Track geometry intersections projection in XY plane for momentum ") + long(i+1);
      gMytrackXY[i] = new TGraph();
      gMytrackXY[i]->SetName(HistName.Data());
      gMytrackXY[i]->SetTitle(HistTitle.Data());
      gMytrackXY[i]->SetMarkerColor(MyColors[i]);
      gMytrackXY[i]->SetLineColor(MyColors[i]);
      gMytrackXY[i]->SetLineWidth(2);
      
      std::vector<IntersectionHit_t> ItersectionHitList;
      ItersectionHitList.clear();
      TrackerList[igeo]->GetTrajectory()->SetParticle(MyParticle[i]);
      TrackerList[igeo]->SetTrajectoryInitConditions(MyOrigin,MyMomentum[i]);
      TrackerList[igeo]->GetIntersectionsWithGeometry(ItersectionHitList);
      int counter_hit = 0;
      for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) { //begin intersection hits loops
	double s;
	TVector3 PositionXYZ;
	
	//First intersection with geometry element
	s = ItersectionHitList[ihit].s_in;
	PositionXYZ = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	gMytrackXY_points[i]->SetPoint(counter_hit,PositionXYZ.X(),PositionXYZ.Y());
	counter_hit++;
	
	//Middle point of intersection inside geometry element
	s = ItersectionHitList[ihit].s;
	PositionXYZ = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	gMytrackXY_points[i]->SetPoint(counter_hit,PositionXYZ.X(),PositionXYZ.Y());
	counter_hit++;
	
	//Second intersection with geometry element
	s = ItersectionHitList[ihit].s_out;
	PositionXYZ = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	gMytrackXY_points[i]->SetPoint(counter_hit,PositionXYZ.X(),PositionXYZ.Y());
	counter_hit++;
      } //end intersection hits loops
      ItersectionHitList.clear();
      
      double MyRs[2];
      MyRs[0] = 0.0;
      MyRs[1] = 7*global->GetUnit("cm");
      int MyNbins_s = int((MyRs[1] - MyRs[0])/(1.0*global->GetUnit("mm")));
      for(int is=0;is<MyNbins_s+1;is++) {
	double s = MyRs[0] + is*(MyRs[1] - MyRs[0])/MyNbins_s;
	TVector3 point = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	
	gMytrackXY[i]->SetPoint(is,point.X(),point.Y());
      }
      
    }

    TLatex* latex = new TLatex();
    latex->SetTextAlign(12);
    latex->SetTextSize(0.08);
    latex->SetTextColor(kRed);
    
    double MyLimit = 7.0;
    TH2F MyhrefXY("MyhrefXY","",
		  100,-MyLimit,MyLimit,
		  100,-MyLimit,MyLimit);
    MyhrefXY.SetXTitle("X (cm)");
    MyhrefXY.GetXaxis()->CenterTitle(true);
    MyhrefXY.SetYTitle("Y (cm)");
    MyhrefXY.GetYaxis()->CenterTitle(true);
    MyhrefXY.SetLineColor(1);
    MyhrefXY.SetLineWidth(1);
    MyhrefXY.GetXaxis()->SetTitleSize(TheSize);
    MyhrefXY.GetXaxis()->SetLabelSize(TheSize);
    MyhrefXY.GetYaxis()->SetTitleSize(TheSize);
    MyhrefXY.GetYaxis()->SetLabelSize(TheSize);
    MyhrefXY.SetStats(false);
    
    c1->Clear();
    MyhrefXY.Draw();
    
    TImage *img = TImage::Open("frijoles_charros.jpg");
    img->SetConstRatio(kFALSE);
    double frac = 0.10;
    TVector3 pos_pad(0.531,0.515,0.0);
    TPad *l = new TPad("l","l",
		       pos_pad.X() - frac,
		       pos_pad.Y() - frac,
		       pos_pad.X() + frac,
		       pos_pad.Y() + frac);
    l->Draw();
    l->cd();
    img->Draw("xxx");
    
    c1->cd();
    for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
      if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
    }
    for(int i=0;i<Nmoms;i++) {
      if(gMytrackXY[i]->GetN() > 0)        gMytrackXY[i]->Draw("PEL");
      if(gMytrackXY_points[i]->GetN() > 0) gMytrackXY_points[i]->Draw("P");
    }
    
    latex->DrawLatex(-1.2,+5.5,"Guari");
    latex->DrawLatex(-1.7,-5.5,"guanchi");
    
    //lhXY->Draw("PEL"); lvXY->Draw("PEL");
#if 0
    TGraph* g_tmp = new TGraph();
    //double R = 1.6 - 300e-4;
    double R = 9.37127;
    int Nbins_tmp = 100;
    for(int i=0;i<Nbins_tmp+1;i++) {
      double alpha = 0.0 + i*(2*TMath::Pi())/Nbins_tmp;
      double x = R*TMath::Cos(alpha);
      double y = R*TMath::Sin(alpha);
      g_tmp->SetPoint(i,x,y);
    }
    g_tmp->SetLineColor(kRed);
    g_tmp->SetLineWidth(2);
    g_tmp->SetLineStyle(2);
    g_tmp->Draw("PEL");
#endif
    c1->Print(EPSName.Data());
#endif
    
    if(DoRZGeoRepresentation) {
      c1->Clear();
      hrefZR.Draw();
      for(int iplane=0;iplane<int(GraphZR[igeo].size());iplane++) {
        if(GraphZR[igeo][iplane].GetN() > 0) GraphZR[igeo][iplane].Draw("PEL");
      }
      lhZR->Draw("PEL");
      lvZR->Draw("PEL");
      
      c1->cd();
      PlotLogo(0.17,0.89,0.9,0.9);
      c1->Print(EPSName.Data());
    }
    
    for(int irange=0;irange<NgeoRanges;irange++) {
      HistName = GeometryList[igeo]->GetName() + TString(" geometry ZY (Range") + long(irange+1) + TString(")");
      hrefZY_geoRange[irange]->SetTitle(HistName.Data());
      HistName = GeometryList[igeo]->GetName() + TString(" geometry ZX (Range") + long(irange+1) + TString(")");
      hrefZX_geoRange[irange]->SetTitle(HistName.Data());
      HistName = GeometryList[igeo]->GetName() + TString(" geometry XY (Range") + long(irange+1) + TString(")");
      hrefXY_geoRange[irange]->SetTitle(HistName.Data());
      HistName = GeometryList[igeo]->GetName() + TString(" geometry ZR (Range") + long(irange+1) + TString(")");
      hrefZR_geoRange[irange]->SetTitle(HistName.Data());
    
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      hrefZY_geoRange[irange]->Draw();
      for(int iplane=0;iplane<int(GraphZY[igeo].size());iplane++) {
	if(GraphZY[igeo][iplane].GetN() > 0) GraphZY[igeo][iplane].Draw("PEL");
      }
      lhZY->Draw("PEL");
      lvZY->Draw("PEL");
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      hrefZX_geoRange[irange]->Draw();
      for(int iplane=0;iplane<int(GraphZX[igeo].size());iplane++) {
	if(GraphZX[igeo][iplane].GetN() > 0) GraphZX[igeo][iplane].Draw("PEL");
      }
      lhZX->Draw("PEL");
      lvZX->Draw("PEL");
      c1->cd(3);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      hrefXY_geoRange[irange]->Draw();
      for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
	if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
      }
      lhXY->Draw("PEL");
      lvXY->Draw("PEL");
      c1->cd();
      PlotLogo(0.17,0.89,0.5,0.5);
      c1->Print(EPSName.Data());
    
      if(DoRZGeoRepresentation) {
        c1->Clear();
        hrefZR_geoRange[irange]->Draw();
        for(int iplane=0;iplane<int(GraphZR[igeo].size());iplane++) {
	  if(GraphZR[igeo][iplane].GetN() > 0) GraphZR[igeo][iplane].Draw("PEL");
        }
        lhZR->Draw("PEL");
        lvZR->Draw("PEL");
	
	c1->cd();
	PlotLogo(0.17,0.89,0.9,0.9);
        c1->Print(EPSName.Data());
      }
      
    }
  } //End of loop over geometries
  
  c1->Print(EPSNameC.Data());

  
  if(SavePlots) {
    cout << endl;
    TString ROOTName = TheOutputFile + TString(".root");
    cout << "Saving geometry visualisations to " << ROOTName.Data() << " file" << endl;
    TFile file(ROOTName.Data(),"UPDATE");
    
    TCanvas* c_geo[NGeometries];
    TCanvas* c_geo2[NGeometries];
    for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
      HistName = GeometryList[igeo]->GetName() + TString(" geometry ZY");
      hrefZY.SetTitle(HistName.Data());
      HistName = GeometryList[igeo]->GetName() + TString(" geometry ZX");
      hrefZX.SetTitle(HistName.Data());
      HistName = GeometryList[igeo]->GetName() + TString(" geometry XY");
      hrefXY.SetTitle(HistName.Data());
      HistName = GeometryList[igeo]->GetName() + TString(" geometry ZR");
      hrefZR.SetTitle(HistName.Data());
    
      HistName  = TString("c_geo") + long(igeo+1);
      HistTitle = GeometryList[igeo]->GetName();
      c_geo[igeo] = new TCanvas(HistName.Data(),HistTitle.Data(),LargeCanvasX,LargeCanvasY);
      //c_geo[igeo] = new TCanvas(HistName.Data(),HistTitle.Data());
      
      HistName  = TString("c_geo2") + long(igeo+1);
      HistTitle = GeometryList[igeo]->GetName();
      c_geo2[igeo] = new TCanvas(HistName.Data(),HistTitle.Data(),LargeCanvasX,LargeCanvasY);
      //c_geo2[igeo] = new TCanvas(HistName.Data(),HistTitle.Data());
    
      c_geo[igeo]->Clear();
      c_geo[igeo]->Divide(2,2);
      c_geo[igeo]->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      hrefZY.Draw();
      for(int iplane=0;iplane<int(GraphZY[igeo].size());iplane++) {
	if(GraphZY[igeo][iplane].GetN() > 0) GraphZY[igeo][iplane].Draw("PEL");
      }
      lhZY->Draw("PEL");
      lvZY->Draw("PEL");
      c_geo[igeo]->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      hrefZX.Draw();
      for(int iplane=0;iplane<int(GraphZX[igeo].size());iplane++) {
	if(GraphZX[igeo][iplane].GetN() > 0) GraphZX[igeo][iplane].Draw("PEL");
      }
      lhZX->Draw("PEL");
      lvZX->Draw("PEL");
      c_geo[igeo]->cd(3);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      hrefXY.Draw();
      for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
	if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
      }
      lhXY->Draw("PEL");
      lvXY->Draw("PEL");
      
      c_geo[igeo]->cd();
      PlotLogo(0.15,0.89,0.51,0.5);
      c_geo[igeo]->Write();
      
      if(DoRZGeoRepresentation) {
        c_geo2[igeo]->Clear();
        hrefZR.Draw();
        for(int iplane=0;iplane<int(GraphZR[igeo].size());iplane++) {
	  if(GraphZR[igeo][iplane].GetN() > 0) GraphZR[igeo][iplane].Draw("PEL");
        }
        lhZR->Draw("PEL");
        lvZR->Draw("PEL");
	
	c_geo2[igeo]->cd();
	PlotLogo(0.17,0.89,0.9,0.9);
	c_geo2[igeo]->Write();
      }
    } //End of loop over geometries
    
    file.Close();
  }


  if(DoPlotSomeTracks) { //Begin if plot some tracks
    double phi;
    int Somebins = Nbins_mom_redu;
    int Thebins  = momArr.size();
    if(Thebins > Somebins) Thebins = Somebins;
    if(UseAllMomVals_GeoPrint)  Thebins = momArr.size();
    const int MyNbins_mom(Thebins);

    //int Nbins_s = 500;
    //double Rs[2];
    //Rs[0] = 0;
    //Rs[1] = 10*sqrt(pow(RX[1] - RX[0],2) + pow(RY[1] - RY[0],2) + pow(RZ[1] - RZ[0],2))*global->GetDistanceUnit("cm");
    
    TString* EPSName_trk  = new TString[NGeometries];
    TString* EPSNameO_trk = new TString[NGeometries];
    TString* EPSNameC_trk = new TString[NGeometries];
    TCanvas* c_geo_trk[NGeometries];
    for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
      GlobalFileCounter++;
      EPSName_trk[igeo]  = TheOutputFile + TString("_") + long(GlobalFileCounter) + TString(".eps");
      EPSNameO_trk[igeo] = EPSName_trk[igeo] + TString("[");
      EPSNameC_trk[igeo] = EPSName_trk[igeo] + TString("]");
      
      HistName = TString("c_geo_trk") + long(igeo+1);
      c_geo_trk[igeo] = new TCanvas(HistName.Data(),HistName.Data(),LargeCanvasX,LargeCanvasY);
      //c_geo_trk[igeo] = new TCanvas(HistName.Data(),HistName.Data());
      c_geo_trk[igeo]->SetFillColor(10);
      c_geo_trk[igeo]->SetFrameFillColor(10);
      c_geo_trk[igeo]->SetTickx(1);
      c_geo_trk[igeo]->SetTicky(1);
      c_geo_trk[igeo]->SetLeftMargin(0.15);
      c_geo_trk[igeo]->SetBottomMargin(0.15);      
      c_geo_trk[igeo]->Print(EPSNameO_trk[igeo].Data());
    }
    
    cout << endl;
    cout << "Start plot some tracks  ";
    global->fWatch.Print();
    global->fWatch.Continue();
    fPrintFreq = 50;
    long counter_bin = 0;
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //theta loop
      double theta = thetaArr[itheta];
      for(int iphi=0;iphi<NphiPoints;iphi++) { //phi loop
        phi = phiArr[iphi];
        int counter = 0;

	TGraph** gtrackXY[MyNbins_mom];
	TGraph** gtrackZY[MyNbins_mom];
	TGraph** gtrackZX[MyNbins_mom];
	TGraph** gtrackZR[MyNbins_mom];
	
	TGraph** gtrackXY_points[MyNbins_mom];
	TGraph** gtrackZY_points[MyNbins_mom];
	TGraph** gtrackZX_points[MyNbins_mom];
	TGraph** gtrackZR_points[MyNbins_mom];

	int**      Nbins_s_mom = new int*[MyNbins_mom];
	double***  Rs_mom      = new double**[MyNbins_mom];
        for(int ip=0;ip<MyNbins_mom;ip++) { //momentum loop
	  Rs_mom[ip]      = new double*[NGeometries];
	  Nbins_s_mom[ip] = new int[NGeometries];
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    Rs_mom[ip][igeo] = new double[2];
	    
	    Nbins_s_mom[ip][igeo] = 500;
	    Rs_mom[ip][igeo][0]   = 0.0;
	    Rs_mom[ip][igeo][1]   = 10*sqrt(pow(RX[1] - RX[0],2) + pow(RY[1] - RY[0],2) + pow(RZ[1] - RZ[0],2))*global->GetDistanceUnit("cm");
	  }
	  
	  double mom;
	  if(UseAllMomVals_GeoPrint) mom = momArr[ip];
	  else                       mom = momArr[0] + (ip+0.5)*(momArr[momArr.size()-1] - momArr[0])/MyNbins_mom;

	  TVector3 momentum = global->GetMomentum(global->GetMomentumFromMomVar(mom,particle,momVariable),theta,phi);

	  counter_bin++;
	  
	  int mycolor = (counter+1==10)?49:(counter+1);
          if(mycolor == 5 || mycolor == 7) {
            mycolor++;
            counter++;
          }
          if(mycolor == 3) mycolor = kGreen+2;
          counter++;

	  std::vector<IntersectionHit_t>* ItersectionHitList = new std::vector<IntersectionHit_t>[NGeometries];
	  for(int igeo=0;igeo<NGeometries;igeo++) { //geometries loop
	    TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);
	    TrackerList[igeo]->GetIntersectionsWithGeometry(ItersectionHitList[igeo]);
	  } //end of loop on geometries

	  double smax = -1.0e+20;
	  for(int igeo=0;igeo<NGeometries;igeo++) { //geometries loop
	    for(int ihit=0;ihit<int(ItersectionHitList[igeo].size());ihit++) { //loop on hits
	      //if(!ItersectionHitList[igeo][ihit].IsSensitivePoint) continue;
	      
	      double s_tmp = ItersectionHitList[igeo][ihit].s;
	      if(smax < s_tmp) smax = s_tmp;
	      
	      s_tmp = ItersectionHitList[igeo][ihit].s_out;
	      if(smax < s_tmp) smax = s_tmp;
	      
	    } // end of loop on hits
	    double s_world = TrackerList[igeo]->GetIntersectionWithWorldVolume();
	    if(s_world < 1.0e+4*global->GetUnit("m")) {
	      if(smax < s_world) smax = s_world;
	    }
	    else {
	      if(smax >= 0) smax = Rs_mom[ip][igeo][1];
	    }
	    
	    double Extension = 0.0;
	    double h = 1.0*global->GetDistanceUnit("mm");
	    if(smax > -1.0) {
	      Rs_mom[ip][igeo][1]   = smax*(1.0 + Extension);
	    }
	    Nbins_s_mom[ip][igeo] = int((Rs_mom[ip][igeo][1] - Rs_mom[ip][igeo][0])/h);
	  } //end of loop on geometries

	  gtrackXY_points[ip] = new TGraph*[NGeometries];
	  gtrackZY_points[ip] = new TGraph*[NGeometries];
	  gtrackZX_points[ip] = new TGraph*[NGeometries];
	  gtrackZR_points[ip] = new TGraph*[NGeometries];
	  for(int igeo=0;igeo<NGeometries;igeo++) { //geometries loop
	    char ccc1[100];
	    char ccc2[100];
	    char ccc3[100];
            sprintf(ccc1,"%.1f",mom/global->GetMomentumUnit("GeV/c"));
	    if(polarVariable == TString("theta"))         sprintf(ccc2,"#theta = %.1f deg",theta/global->GetAngleUnit("deg"));
	    else if(polarVariable == TString("costheta")) sprintf(ccc2,"cos(#theta) = %.3f",global->FromThetaToCosTheta(theta));
	    else if(polarVariable == TString("eta"))      sprintf(ccc2,"#eta = %.3f",global->FromThetaToEta(theta));
	    sprintf(ccc2,"%.1f",theta/global->GetAngleUnit("deg"));
	    sprintf(ccc3,"%.1f",phi/global->GetAngleUnit("deg"));
	    
	    //double MarkerSize  = 0.8;
	    double MarkerSize  = 1.2;
	    int    MarkerStyle = 20;
	    
	    HistName   = TString("gtrackXY_points_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track geometry intersections projection in XY plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackXY_points[ip][igeo] = new TGraph();
	    gtrackXY_points[ip][igeo]->SetName(HistName.Data());
            gtrackXY_points[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackXY_points[ip][igeo]->SetMarkerColor(mycolor);
	    gtrackXY_points[ip][igeo]->SetMarkerSize(MarkerSize);
	    gtrackXY_points[ip][igeo]->SetMarkerStyle(MarkerStyle);
            gtrackXY_points[ip][igeo]->SetLineColor(mycolor);
            gtrackXY_points[ip][igeo]->SetLineWidth(2);
	    
	    HistName   = TString("gtrackZY_points_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track geometry intersections projection in ZY plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackZY_points[ip][igeo] = new TGraph();
	    gtrackZY_points[ip][igeo]->SetName(HistName.Data());
            gtrackZY_points[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackZY_points[ip][igeo]->SetMarkerColor(mycolor);
	    gtrackZY_points[ip][igeo]->SetMarkerSize(MarkerSize);
	    gtrackZY_points[ip][igeo]->SetMarkerStyle(MarkerStyle);
            gtrackZY_points[ip][igeo]->SetLineColor(mycolor);
            gtrackZY_points[ip][igeo]->SetLineWidth(2);
	    
	    HistName   = TString("gtrackZX_points_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track geometry intersections projection in ZX plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackZX_points[ip][igeo] = new TGraph();
	    gtrackZX_points[ip][igeo]->SetName(HistName.Data());
            gtrackZX_points[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackZX_points[ip][igeo]->SetMarkerColor(mycolor);
	    gtrackZX_points[ip][igeo]->SetMarkerSize(MarkerSize);
	    gtrackZX_points[ip][igeo]->SetMarkerStyle(MarkerStyle);
            gtrackZX_points[ip][igeo]->SetLineColor(mycolor);
            gtrackZX_points[ip][igeo]->SetLineWidth(2);
	    
	    HistName   = TString("gtrackZR_points_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track geometry intersections projection in ZR plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackZR_points[ip][igeo] = new TGraph();
	    gtrackZR_points[ip][igeo]->SetName(HistName.Data());
            gtrackZR_points[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackZR_points[ip][igeo]->SetMarkerColor(mycolor);
	    gtrackZR_points[ip][igeo]->SetMarkerSize(MarkerSize);
	    gtrackZR_points[ip][igeo]->SetMarkerStyle(MarkerStyle);
            gtrackZR_points[ip][igeo]->SetLineColor(mycolor);
            gtrackZR_points[ip][igeo]->SetLineWidth(2);

	    int counter_hit = 0;
	    for(int ihit=0;ihit<int(ItersectionHitList[igeo].size());ihit++) { //begin intersection hits loops
	      double s;
	      TVector3 PositionXYZ;

	      //First intersection with geometry element
	      s = ItersectionHitList[igeo][ihit].s_in;
	      PositionXYZ = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	      gtrackXY_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.X(),PositionXYZ.Y());
	      gtrackZY_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),PositionXYZ.Y());
	      gtrackZX_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),PositionXYZ.X());
	      gtrackZR_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),sqrt(pow(PositionXYZ.X(),2) + pow(PositionXYZ.Y(),2)));
	      counter_hit++;
	      
	      //Middle point of intersection inside geometry element
	      s = ItersectionHitList[igeo][ihit].s;
	      PositionXYZ = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	      gtrackXY_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.X(),PositionXYZ.Y());
	      gtrackZY_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),PositionXYZ.Y());
	      gtrackZX_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),PositionXYZ.X());
	      gtrackZR_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),sqrt(pow(PositionXYZ.X(),2) + pow(PositionXYZ.Y(),2)));
	      counter_hit++;
	      
	      //Second intersection with geometry element
	      s = ItersectionHitList[igeo][ihit].s_out;
	      PositionXYZ = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	      gtrackXY_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.X(),PositionXYZ.Y());
	      gtrackZY_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),PositionXYZ.Y());
	      gtrackZX_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),PositionXYZ.X());
	      gtrackZR_points[ip][igeo]->SetPoint(counter_hit,PositionXYZ.Z(),sqrt(pow(PositionXYZ.X(),2) + pow(PositionXYZ.Y(),2)));
	      counter_hit++;
	      
	    } //end intersection hits loops
	    ItersectionHitList[igeo].clear();
	    
	  } // end of geometries loop

	  gtrackXY[ip] = new TGraph*[NGeometries];
	  gtrackZY[ip] = new TGraph*[NGeometries];
	  gtrackZX[ip] = new TGraph*[NGeometries];
	  gtrackZR[ip] = new TGraph*[NGeometries];
	  for(int igeo=0;igeo<NGeometries;igeo++) { //begin geometries loop
	    char ccc1[100];
	    char ccc2[100];
	    char ccc3[100];
            sprintf(ccc1,"%.1f",mom/global->GetMomentumUnit("GeV/c"));
	    if(polarVariable == TString("theta"))         sprintf(ccc2,"#theta = %.1f deg",theta/global->GetAngleUnit("deg"));
	    else if(polarVariable == TString("costheta")) sprintf(ccc2,"cos(#theta) = %.3f",global->FromThetaToCosTheta(theta));
	    else if(polarVariable == TString("eta"))      sprintf(ccc2,"#eta = %.3f",global->FromThetaToEta(theta));	    
	    sprintf(ccc3,"%.1f",phi/global->GetAngleUnit("deg"));

	    HistName   = TString("gtrackXY_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track projection in XY plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackXY[ip][igeo] = new TGraph(Nbins_s_mom[ip][igeo]);
	    gtrackXY[ip][igeo]->SetName(HistName.Data());
            gtrackXY[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackXY[ip][igeo]->SetMarkerColor(mycolor);
            gtrackXY[ip][igeo]->SetLineColor(mycolor);
            gtrackXY[ip][igeo]->SetLineWidth(2);

	    HistName   = TString("gtrackZY_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track projection in ZY plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackZY[ip][igeo] = new TGraph(Nbins_s_mom[ip][igeo]);
	    gtrackZY[ip][igeo]->SetName(HistName.Data());
            gtrackZY[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackZY[ip][igeo]->SetMarkerColor(mycolor);
            gtrackZY[ip][igeo]->SetLineColor(mycolor);
            gtrackZY[ip][igeo]->SetLineWidth(2);
	  
	    HistName   = TString("gtrackZX_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track projection in ZX plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackZX[ip][igeo] = new TGraph(Nbins_s_mom[ip][igeo]);
	    gtrackZX[ip][igeo]->SetName(HistName.Data());
            gtrackZX[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackZX[ip][igeo]->SetMarkerColor(mycolor);
            gtrackZX[ip][igeo]->SetLineColor(mycolor);
            gtrackZX[ip][igeo]->SetLineWidth(2);
	    
	    HistName   = TString("gtrackZR_mon") + long(ip+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_igeo") + long(igeo+1);
	    HistTitle  = TString("Track projection in ZR plane for momentum = ") + TString(ccc1) + 
	                 TString(" GeV/c, ") + TString(ccc2) + TString(" and #phi = ") + TString(ccc3) + TString(" deg") + TString(", geo-") + long(igeo+1);
	    gtrackZR[ip][igeo] = new TGraph(Nbins_s_mom[ip][igeo]);
	    gtrackZR[ip][igeo]->SetName(HistName.Data());
            gtrackZR[ip][igeo]->SetTitle(HistTitle.Data());
            gtrackZR[ip][igeo]->SetMarkerColor(mycolor);
            gtrackZR[ip][igeo]->SetLineColor(mycolor);
            gtrackZR[ip][igeo]->SetLineWidth(2);

	    for(int is=0;is<Nbins_s_mom[ip][igeo]+1;is++) {
	      double s = Rs_mom[ip][igeo][0] + is*(Rs_mom[ip][igeo][1] - Rs_mom[ip][igeo][0])/Nbins_s_mom[ip][igeo];
	      TVector3 point = (1.0/global->GetDistanceUnit("cm"))*TrackerList[igeo]->GetTrajectory()->GetTrueTrajectoryCoordinates(s);
	  
	      gtrackXY[ip][igeo]->SetPoint(is,point.X(),point.Y());
	      gtrackZY[ip][igeo]->SetPoint(is,point.Z(),point.Y());
	      gtrackZX[ip][igeo]->SetPoint(is,point.Z(),point.X());
	      gtrackZR[ip][igeo]->SetPoint(is,point.Z(),sqrt(pow(point.X(),2) + pow(point.Y(),2)));
	    }

	  } //end geometries loop

	  if(!((counter_bin)%fPrintFreq)) {
	    char ccc2[100];
	    if(polarVariable == TString("theta"))         sprintf(ccc2,"theta = %.1f deg",theta/global->GetAngleUnit("deg"));
	    else if(polarVariable == TString("costheta")) sprintf(ccc2,"cos(theta) = %.3f",global->FromThetaToCosTheta(theta));
	    else if(polarVariable == TString("eta"))      sprintf(ccc2,"eta = %.3f",global->FromThetaToEta(theta));
	    
	    cout << counter_bin << " momentum, theta and phi bins processed out of " << NthetaPoints*NphiPoints*MyNbins_mom << ". "
	         << " Mom-bin "   << ip+1     << ", momentum = " << mom/global->GetMomentumUnit("GeV/c") << " GeV/c; "
		 << " theta-bin " << itheta+1 << ", " << ccc2 << "; "
		 << " phi-bin "   << iphi+1   << ", phi = " << phi/global->GetAngleUnit("deg") << " deg !!! ";
	    global->fWatch.Print();
	    global->fWatch.Continue();
          }

          for(int igeo=0;igeo<NGeometries;igeo++) {
	    delete [] Rs_mom[ip][igeo];
	  }
          delete [] Rs_mom[ip];
	  delete [] Nbins_s_mom[ip];
          delete [] ItersectionHitList;
        } //end momentum loop
        delete [] Nbins_s_mom;
	delete [] Rs_mom;

        //Now plot the geometries and the tracks
        for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
	  TLegend* leg = new TLegend(0.1,0.1,0.9,0.9);
	  leg->SetFillColor(10);
	  
	  char thetheta[100];
	  if(polarVariable == TString("theta"))         sprintf(thetheta,"#theta = %.1f deg",theta/global->GetAngleUnit("deg"));
	  else if(polarVariable == TString("costheta")) sprintf(thetheta,"cos(#theta) = %.3f",global->FromThetaToCosTheta(theta));
	  else if(polarVariable == TString("eta"))      sprintf(thetheta,"#eta = %.3f",global->FromThetaToEta(theta));	    
	  
	  char thephi[100];
	  sprintf(thephi,"%.1f",phi/global->GetAngleUnit("deg"));

	  HistName = GeometryList[igeo]->GetName() + TString(" geometry ZY for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
          hrefZY.SetTitle(HistName.Data());
	  HistName = GeometryList[igeo]->GetName() + TString(" geometry ZX for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
          hrefZX.SetTitle(HistName.Data());
	  HistName = GeometryList[igeo]->GetName() + TString(" geometry XY for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
          hrefXY.SetTitle(HistName.Data());
	  HistName = GeometryList[igeo]->GetName() + TString(" geometry ZR for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
          hrefZR.SetTitle(HistName.Data());

	  c_geo_trk[igeo]->Clear();
          c_geo_trk[igeo]->Divide(2,2);
          c_geo_trk[igeo]->cd(1);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.15);
          hrefZY.Draw();
	  for(int iplane=0;iplane<int(GraphZY[igeo].size());iplane++) {
	    if(GraphZY[igeo][iplane].GetN() > 0) GraphZY[igeo][iplane].Draw("PEL");
	  }

	  for(int ip=0;ip<MyNbins_mom;ip++) { //momentum loop
	    double mom;
	    if(UseAllMomVals_GeoPrint) mom = momArr[ip];
	    else                       mom = momArr[0] + (ip+0.5)*(momArr[momArr.size()-1] - momArr[0])/MyNbins_mom;
	  
	    char themom[100];
	    TString MomUnit_tmp("GeV/c");
	    if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("GeV");
	    double mom_tmp = mom/global->GetUnit(MomUnit_tmp);
	    if(mom_tmp < 0.01) {
	      MomUnit_tmp = TString("MeV/c");
	      if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("MeV");
	      mom_tmp = mom/global->GetUnit(MomUnit_tmp);
	      
	      if(mom_tmp < 0.01) {
		MomUnit_tmp = TString("keV/c");
		if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("keV");
		mom_tmp = mom/global->GetUnit(MomUnit_tmp);
	      }
	    }
	    
	    sprintf(themom,"%.3f",mom_tmp);
	    TString MyMomVar("p");
	    if(momVariable == TString("E"))             MyMomVar = TString("E");
	    else if(momVariable == TString("Ekin"))     MyMomVar = TString("E_{kin}");
	    else if(momVariable == TString("EkinPerU")) MyMomVar = TString("E_{kin}/u");
	    
	    TString MyLegend = MyMomVar +  TString(" = ") + TString(themom) + TString(" ") + TString(MomUnit_tmp);
	    
	    leg->AddEntry(gtrackZY[ip][igeo],MyLegend.Data(),"l");

	    if(gtrackZY[ip][igeo]->GetN() > 0)        gtrackZY[ip][igeo]->Draw("PEL");
	    if(gtrackZY_points[ip][igeo]->GetN() > 0) gtrackZY_points[ip][igeo]->Draw("P");
	  }
          lhZY->Draw("PEL");
          lvZY->Draw("PEL");

          c_geo_trk[igeo]->cd(2);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.15);
          hrefZX.Draw();
	  for(int iplane=0;iplane<int(GraphZX[igeo].size());iplane++) {
	    if(GraphZX[igeo][iplane].GetN() > 0) GraphZX[igeo][iplane].Draw("PEL");
	  }
        
          for(int ip=0;ip<MyNbins_mom;ip++) {
	    if(gtrackZX[ip][igeo]->GetN() > 0)        gtrackZX[ip][igeo]->Draw("PEL");
	    if(gtrackZX_points[ip][igeo]->GetN() > 0) gtrackZX_points[ip][igeo]->Draw("P");
	  }
        
          lhZX->Draw("PEL");
          lvZX->Draw("PEL");

          c_geo_trk[igeo]->cd(3);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.15);
          hrefXY.Draw();
	  for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
	    if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
	  }
	
	  for(int ip=0;ip<MyNbins_mom;ip++) {
	    if(gtrackXY[ip][igeo]->GetN() > 0)        gtrackXY[ip][igeo]->Draw("PEL");
	    if(gtrackXY_points[ip][igeo]->GetN() > 0) gtrackXY_points[ip][igeo]->Draw("P");
	  }
          lhXY->Draw("PEL");
          lvXY->Draw("PEL");
	
	  c_geo_trk[igeo]->cd(4);
	  leg->Draw();
	  
	  c_geo_trk[igeo]->cd();
	  PlotLogo(0.15,0.89,0.51,0.5);
	  c_geo_trk[igeo]->Print(EPSName_trk[igeo].Data());
	  
	  if(DoRZGeoRepresentation) {
	    c_geo_trk[igeo]->Clear();
            hrefZR.Draw();
	    for(int iplane=0;iplane<int(GraphZR[igeo].size());iplane++) {
	      if(GraphZR[igeo][iplane].GetN() > 0) GraphZR[igeo][iplane].Draw("PEL");
	    }
	    for(int ip=0;ip<MyNbins_mom;ip++) {
	      if(gtrackZR[ip][igeo]->GetN() > 0)        gtrackZR[ip][igeo]->Draw("PEL");
	      if(gtrackZR_points[ip][igeo]->GetN() > 0) gtrackZR_points[ip][igeo]->Draw("P");
	    }
            lhZR->Draw("PEL");
            lvZR->Draw("PEL");
	    
	    c_geo_trk[igeo]->cd();
	    PlotLogo(0.17,0.89,0.9,0.9);
	    c_geo_trk[igeo]->Print(EPSName_trk[igeo].Data());
	  }
	  
	  for(int irange=0;irange<NgeoRanges;irange++) {
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZY for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefZY_geoRange[irange]->SetTitle(HistName.Data());
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZX for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefZX_geoRange[irange]->SetTitle(HistName.Data());
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry XY for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefXY_geoRange[irange]->SetTitle(HistName.Data());
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZR for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefZR_geoRange[irange]->SetTitle(HistName.Data());
	    
	    c_geo_trk[igeo]->Clear();
	    c_geo_trk[igeo]->Divide(2,2);
	    c_geo_trk[igeo]->cd(1);
	    gPad->SetFillColor(10);
	    gPad->SetFrameFillColor(10);
	    gPad->SetTickx(1);
	    gPad->SetTicky(1);
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    hrefZY_geoRange[irange]->Draw();
	    for(int iplane=0;iplane<int(GraphZY[igeo].size());iplane++) {
	      if(GraphZY[igeo][iplane].GetN() > 0) GraphZY[igeo][iplane].Draw("PEL");
	    }
	    for(int ip=0;ip<MyNbins_mom;ip++) {
	      if(gtrackZY[ip][igeo]->GetN() > 0)        gtrackZY[ip][igeo]->Draw("PEL");
	      if(gtrackZY_points[ip][igeo]->GetN() > 0) gtrackZY_points[ip][igeo]->Draw("P");
	    }
	    lhZY->Draw("PEL");
	    lvZY->Draw("PEL");
	    c_geo_trk[igeo]->cd(2);
	    gPad->SetFillColor(10);
	    gPad->SetFrameFillColor(10);
	    gPad->SetTickx(1);
	    gPad->SetTicky(1);
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    hrefZX_geoRange[irange]->Draw();
	    for(int iplane=0;iplane<int(GraphZX[igeo].size());iplane++) {
	      if(GraphZX[igeo][iplane].GetN() > 0) GraphZX[igeo][iplane].Draw("PEL");
	    }
	    for(int ip=0;ip<MyNbins_mom;ip++) {
	      if(gtrackZX[ip][igeo]->GetN() > 0)        gtrackZX[ip][igeo]->Draw("PEL");
	      if(gtrackZX_points[ip][igeo]->GetN() > 0) gtrackZX_points[ip][igeo]->Draw("P");
	    }
	    lhZX->Draw("PEL");
	    lvZX->Draw("PEL");
	    c_geo_trk[igeo]->cd(3);
	    gPad->SetFillColor(10);
	    gPad->SetFrameFillColor(10);
	    gPad->SetTickx(1);
	    gPad->SetTicky(1);
	    gPad->SetLeftMargin(0.15);
	    gPad->SetBottomMargin(0.15);
	    hrefXY_geoRange[irange]->Draw();
	    for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
	      if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
	    }
	    for(int ip=0;ip<MyNbins_mom;ip++) {
	      if(gtrackXY[ip][igeo]->GetN() > 0)        gtrackXY[ip][igeo]->Draw("PEL");
	      if(gtrackXY_points[ip][igeo]->GetN() > 0) gtrackXY_points[ip][igeo]->Draw("P");
	    }
	    lhXY->Draw("PEL");
	    lvXY->Draw("PEL");
	    c_geo_trk[igeo]->cd(4);
	    leg->Draw();
	    
	    c_geo_trk[igeo]->cd();
	    PlotLogo(0.15,0.89,0.51,0.5);
	    c_geo_trk[igeo]->Print(EPSName_trk[igeo].Data());
	    
	    if(DoRZGeoRepresentation) {
	      hrefZR_geoRange[irange]->Draw();
	      for(int iplane=0;iplane<int(GraphZR[igeo].size());iplane++) {
	        if(GraphZR[igeo][iplane].GetN() > 0) GraphZR[igeo][iplane].Draw("PEL");
	      }
	      for(int ip=0;ip<MyNbins_mom;ip++) {
	        if(gtrackZR[ip][igeo]->GetN() > 0)        gtrackZR[ip][igeo]->Draw("PEL");
	        if(gtrackZR_points[ip][igeo]->GetN() > 0) gtrackZR_points[ip][igeo]->Draw("P");
	      }
	      lhZR->Draw("PEL");
	      lvZR->Draw("PEL");
	      
	      c_geo_trk[igeo]->cd();
	      PlotLogo(0.17,0.89,0.9,0.9);
	      c_geo_trk[igeo]->Print(EPSName_trk[igeo].Data());
	    }
	    
	  }

	} //End of loop over geometries

        if(SavePlots) {
	  double MarkerSize  = 0.8;
	  
          TString ROOTName = TheOutputFile + TString(".root");
          TFile file(ROOTName.Data(),"UPDATE");

	  //Now plot the geometries and the tracks
          for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
	    TLegend* leg = new TLegend(0.1,0.1,0.9,0.9);
	    leg->SetFillColor(10);
	    
	    char thetheta[100];
	    if(polarVariable == TString("theta"))         sprintf(thetheta,"#theta = %.1f deg",theta/global->GetAngleUnit("deg"));
	    else if(polarVariable == TString("costheta")) sprintf(thetheta,"cos(#theta) = %.3f",global->FromThetaToCosTheta(theta));
	    else if(polarVariable == TString("eta"))      sprintf(thetheta,"#eta = %.3f",global->FromThetaToEta(theta));	    
	  
	    char thephi[100];
	    sprintf(thephi,"%.1f",phi/global->GetAngleUnit("deg"));

	    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZY for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefZY.SetTitle(HistName.Data());
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZX for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefZX.SetTitle(HistName.Data());
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry XY for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefXY.SetTitle(HistName.Data());
	    HistName = GeometryList[igeo]->GetName() + TString(" geometry ZR for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
            hrefZR.SetTitle(HistName.Data());

	    HistName  = TString("c_geo") + long(igeo+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1);
	    HistTitle = GeometryList[igeo]->GetName() + TString(" geometry for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
	    TCanvas* tmp_trk = new TCanvas(HistName.Data(),HistTitle.Data(),LargeCanvasX,LargeCanvasY);
	    //TCanvas* tmp_trk = new TCanvas(HistName.Data(),HistTitle.Data());
	    tmp_trk->SetFillColor(10);
	    tmp_trk->SetFrameFillColor(10);
	    tmp_trk->SetTickx(1);
	    tmp_trk->SetTicky(1);
	    tmp_trk->SetLeftMargin(0.15);
	    tmp_trk->SetBottomMargin(0.15);

	    tmp_trk->Clear();
	    tmp_trk->Divide(2,2);
	    tmp_trk->cd(1);
	    gPad->SetFillColor(10);
            gPad->SetFrameFillColor(10);
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            hrefZY.Draw();
	    for(int iplane=0;iplane<int(GraphZY[igeo].size());iplane++) {
	      if(GraphZY[igeo][iplane].GetN() > 0) GraphZY[igeo][iplane].Draw("PEL");
	    }

	    for(int ip=0;ip<MyNbins_mom;ip++) { //momentum loop
	      double mom;
	      if(UseAllMomVals_GeoPrint) mom = momArr[ip];
	      else                       mom = momArr[0] + (ip+0.5)*(momArr[momArr.size()-1] - momArr[0])/MyNbins_mom;
	  
	      char themom[100];
	      TString MomUnit_tmp("GeV/c");
	      if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("GeV");
	      double mom_tmp = mom/global->GetUnit(MomUnit_tmp);
	      if(mom_tmp < 0.01) {
	        MomUnit_tmp = TString("MeV/c");
	        if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("MeV");
	        mom_tmp = mom/global->GetUnit(MomUnit_tmp);
	      
	        if(mom_tmp < 0.01) {
		  MomUnit_tmp = TString("keV/c");
		  if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("keV");
		  mom_tmp = mom/global->GetUnit(MomUnit_tmp);
	        }
	      }
	      
	      sprintf(themom,"%.3f",mom_tmp);
	      TString MyMomVar("p");
	      if(momVariable == TString("E"))             MyMomVar = TString("E");
	      else if(momVariable == TString("Ekin"))     MyMomVar = TString("E_{kin}");
	      else if(momVariable == TString("EkinPerU")) MyMomVar = TString("E_{kin}/u");
	    
	      TString MyLegend = MyMomVar +  TString(" = ") + TString(themom) + TString(" ") + TString(MomUnit_tmp);
	      leg->AddEntry(gtrackZY[ip][igeo],MyLegend.Data(),"l");
	      
	      gtrackZY[ip][igeo]->SetMarkerSize(MarkerSize);
	      gtrackZY_points[ip][igeo]->SetMarkerSize(MarkerSize);
	      if(gtrackZY[ip][igeo]->GetN() > 0)        gtrackZY[ip][igeo]->Draw("PEL");
	      if(gtrackZY_points[ip][igeo]->GetN() > 0) gtrackZY_points[ip][igeo]->Draw("P");
	    }
            lhZY->Draw("PEL");
            lvZY->Draw("PEL");

            tmp_trk->cd(2);
            gPad->SetFillColor(10);
            gPad->SetFrameFillColor(10);
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            hrefZX.Draw();
	    for(int iplane=0;iplane<int(GraphZX[igeo].size());iplane++) {
	      if(GraphZX[igeo][iplane].GetN() > 0) GraphZX[igeo][iplane].Draw("PEL");
	    }
        
            for(int ip=0;ip<MyNbins_mom;ip++) {
	      gtrackZX[ip][igeo]->SetMarkerSize(MarkerSize);
	      gtrackZX_points[ip][igeo]->SetMarkerSize(MarkerSize);
	      if(gtrackZX[ip][igeo]->GetN() > 0)        gtrackZX[ip][igeo]->Draw("PEL");
	      if(gtrackZX_points[ip][igeo]->GetN() > 0) gtrackZX_points[ip][igeo]->Draw("P");
	    }
        
            lhZX->Draw("PEL");
            lvZX->Draw("PEL");

            tmp_trk->cd(3);
            gPad->SetFillColor(10);
            gPad->SetFrameFillColor(10);
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            hrefXY.Draw();
	    for(int iplane=0;iplane<int(GraphXY[igeo].size());iplane++) {
	      if(GraphXY[igeo][iplane].GetN() > 0) GraphXY[igeo][iplane].Draw("PEL");
	    }
	    for(int ip=0;ip<MyNbins_mom;ip++) {
	      gtrackXY[ip][igeo]->SetMarkerSize(MarkerSize);
	      gtrackXY_points[ip][igeo]->SetMarkerSize(MarkerSize);
	      if(gtrackXY[ip][igeo]->GetN() > 0)        gtrackXY[ip][igeo]->Draw("PEL");
	      if(gtrackXY_points[ip][igeo]->GetN() > 0) gtrackXY_points[ip][igeo]->Draw("P");
	    }
            lhXY->Draw("PEL");
            lvXY->Draw("PEL");
	
	    tmp_trk->cd(4);
	    leg->Draw();
	    
	    tmp_trk->cd();
	    PlotLogo(0.15,0.89,0.51,0.5);
	    tmp_trk->Write();
	  
	    if(DoRZGeoRepresentation) {
	      HistName  = TString("c_geo2") + long(igeo+1) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1);
	      HistTitle = GeometryList[igeo]->GetName() + TString(" geometry for ") + TString(thetheta) + TString(", #phi = ") + TString(thephi) + TString(" deg");
	      TCanvas* tmp_trk2 = new TCanvas(HistName.Data(),HistTitle.Data(),LargeCanvasX,LargeCanvasY);
	      //TCanvas* tmp_trk2 = new TCanvas(HistName.Data(),HistTitle.Data());
	      tmp_trk2->SetFillColor(10);
	      tmp_trk2->SetFrameFillColor(10);
	      tmp_trk2->SetTickx(1);
	      tmp_trk2->SetTicky(1);
	      tmp_trk2->SetLeftMargin(0.15);
	      tmp_trk2->SetBottomMargin(0.15);
	  
	      tmp_trk2->Clear();
              hrefZR.Draw();
	      for(int iplane=0;iplane<int(GraphZR[igeo].size());iplane++) {
	        if(GraphZR[igeo][iplane].GetN() > 0) GraphZR[igeo][iplane].Draw("PEL");
	      }
	      for(int ip=0;ip<MyNbins_mom;ip++) {
		gtrackZR[ip][igeo]->SetMarkerSize(MarkerSize);
	        gtrackZR_points[ip][igeo]->SetMarkerSize(MarkerSize);
	        if(gtrackZR[ip][igeo]->GetN() > 0)        gtrackZR[ip][igeo]->Draw("PEL");
	        if(gtrackZR_points[ip][igeo]->GetN() > 0) gtrackZR_points[ip][igeo]->Draw("P");
	      }
              lhZR->Draw("PEL");
              lvZR->Draw("PEL");
	      
	      tmp_trk2->cd();
	      PlotLogo(0.17,0.89,0.9,0.9);
	      tmp_trk2->Write();
	    }
	    
	  } //End of loop over geometries
	  
          file.Close();
	}

      } //end phi loop
    } //end theta loop
    
    for(int igeo=0;igeo<NGeometries;igeo++) {
      c_geo_trk[igeo]->Print(EPSNameC_trk[igeo].Data());
    }
    
    cout << "End plot some tracks  ";
    global->fWatch.Print();
    global->fWatch.Continue();
    cout << endl;
    
    delete [] EPSName_trk;
    delete [] EPSNameO_trk;
    delete [] EPSNameC_trk;
  } //end if Plots some tracks

  //Freeing the memory
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int i=0;i<int(GraphXY[igeo].size());i++) GraphXY[igeo].clear();
  }
  delete [] GraphXY;
  delete [] GraphZY;
  delete [] GraphZX;
  delete [] GraphZR;
  
  return;
  
}
//====================================================================


