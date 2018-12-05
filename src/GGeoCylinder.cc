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
#include "include/GGeoCylinder.h"

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
GGeoCylinder::GGeoCylinder(TString   aName,
			   int       aGeometryID,
			   int       aGeoObjectIdx,
			   TVector3  aPosition,
			   TMatrixD  aRot,
			   float     aThickness,
			   TString   aMaterial,
			   bool      aIsSensitive,
			   float     aLength,
			   float     aRadius,
			   TVector2  aVInsensitive,
			   GGlobalTools* aglobal,
			   float     aResolutionU,
			   float     aResolutionV,
			   float     aDetEffic,
			   float     aROtime,
			   float     aBkgRate) 
                           : GGeoObject(aName,
					aGeometryID,
					aGeoObjectIdx,
					aPosition,
					aRot,
					aThickness,
					aMaterial,
					aIsSensitive,
					aglobal,
					aResolutionU,aResolutionV,
					aDetEffic,aROtime,aBkgRate)
{
  
  GeoObjectType = TString("Cylinder");
  
  Length       = aLength;
  Radius       = aRadius;
  VInsensitive = aVInsensitive;
  
  CheckInputs();
  
  FillSurfaces();
  
  SetAllSurfacesGeoObjectType();
  
}
//====================================================================
GGeoCylinder::GGeoCylinder(const GGeoCylinder& other,TString aName)
                           : GGeoObject(aName,
				        other.GeometryID,
				        other.GeoObjectIdx,
				        other.Position,
				        other.Rot,
				        other.Thickness,
				        other.Material,
				        other.IsSensitive,
				        other.global,
				        other.ResolutionU,
				        other.ResolutionV,
				        other.DetEffic,
				        other.ROtime,
				        other.BkgRate)
{
  
  GeoObjectType = other.GeoObjectType;
  
  Length       = other.Length;
  Radius       = other.Radius;
  VInsensitive = other.VInsensitive;
  
  ResolutionModel   = other.ResolutionModel;
  ResolutionModelID = other.ResolutionModelID;
  
  CheckInputs();
  
  FillSurfaces();
  
}
//====================================================================
GGeoCylinder::~GGeoCylinder() 
{
  
}
//====================================================================
GGeoObject* GGeoCylinder::clone(TString aName) const
{
 
  return new GGeoCylinder(*this,aName);
  
}
//====================================================================
void  GGeoCylinder::CheckInputs()
{
  
  if(Length < 0.0) {
    cout << endl;
    cout << "Length of Cylinder geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius < 0.0) {
    cout << endl;
    cout << "Radius of Cylinder geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Cylinder geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Cylinder geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() + VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge V-insensive fraction of Cylinder geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
    
  return;
  
}
//====================================================================
void  GGeoCylinder::FillSurfaces()
{
  
  TVector3 MainUVector(1.0,0.0,0.0);
  TVector3 MainVVector(0.0,1.0,0.0);
  TVector3 MainWVector(0.0,0.0,1.0);
  
  TMatrixD  InvRot;
  InvRot.ResizeTo(3,3);
  InvRot = Rot;
  InvRot.Invert();
  global->RotateVector(InvRot,MainUVector);
  global->RotateVector(InvRot,MainVVector);
  global->RotateVector(InvRot,MainWVector);
  
  TString   NameTmp;
  TVector3  tmpUVector,tmpVVector,tmpWVector;
  TMatrixD  tmpRot;
  tmpRot.ResizeTo(3,3);
  TVector3  tmpPosition;
  double    tmpLength;
  double    tmpRadius;
  double    tmpRin;
  double    tmpRout;
  TVector2  tmpVInsensitive(0.0,0.0);
  
  //Creating the main surface object
  NameTmp = TString("Main cylinder surface for Cylinder geometry element ") + GeoObjectName;
  MainSurface = new GSurfaceCylinder(NameTmp,
				     Position,
				     Rot,
				     Length,
				     Radius,
				     VInsensitive,
				     IsSensitive,
				     global);
  
  //Filling the list of 
  BoundarySurfacesList.clear();
  
  //top-boundary
  NameTmp = TString("cylinder top-boundary surface for Cylinder geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength = Length;
  tmpRadius = Radius + Thickness/2.0;
  BoundarySurfacesList.push_back(new GSurfaceCylinder(NameTmp,
						      tmpPosition,
						      tmpRot,
						      tmpLength,
						      tmpRadius,
						      tmpVInsensitive,
						      false,
						      global));
  
  //bottom-boundary
  NameTmp = TString("cylinder bottom-boundary surface for Cylinder geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength = Length;
  tmpRadius = Radius - Thickness/2.0;
  if(tmpRadius > 1.0*global->GetDistanceUnit("um")) {
    BoundarySurfacesList.push_back(new GSurfaceCylinder(NameTmp,
						        tmpPosition,
						        tmpRot,
						        tmpLength,
						        tmpRadius,
						        tmpVInsensitive,
						        false,
						        global));
  }
  
  //fwd disk boundary
  NameTmp = TString("fwd-disk boundary surface for Cylinder geometry element ") + GeoObjectName;
  tmpPosition = Position + (Length/2)*MainWVector;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin    = Radius - (Thickness/2);
  tmpRout   = Radius + (Thickness/2);
  if(tmpRin < 0.0) tmpRin = 0.0;
  BoundarySurfacesList.push_back(new GSurfaceDisk(NameTmp,
						  tmpPosition,
						  tmpRot,
						  tmpRin,
						  tmpRout,
						  tmpVInsensitive,
						  false,
						  global));
  
  //bwd disk boundary
  NameTmp = TString("bwd-disk boundary surface for Cylinder geometry element ") + GeoObjectName;
  tmpPosition = Position - (Length/2)*MainWVector;
  tmpUVector  = (-1)*MainUVector;
  tmpVVector  = (-1)*MainVVector;
  tmpWVector  = (-1)*MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin    = Radius - (Thickness/2);
  tmpRout   = Radius + (Thickness/2);
  if(tmpRin < 0.0) tmpRin = 0.0;
  BoundarySurfacesList.push_back(new GSurfaceDisk(NameTmp,
						  tmpPosition,
						  tmpRot,
						  tmpRin,
						  tmpRout,
						  tmpVInsensitive,
						  false,
						  global));
  
  return;
  
}
//====================================================================
bool  GGeoCylinder::IsPointInsideGeometry(TVector3 PosXYZ)
{
  
  for(int kkk=0;kkk<int(BoundarySurfacesList.size());kkk++) {
    double W = (BoundarySurfacesList[kkk]->GetUVWFromXYZ(PosXYZ)).Z();
    if(BoundarySurfacesList[kkk]->GetName().BeginsWith("cylinder bottom-boundary surface")) W *= -1;
    
    if(W > 0.0) return false;
  }
  
  return true;
  
}
//====================================================================
TVector3  GGeoCylinder::GetBoundaryNormVector(int idx, TVector3 PosXYZ)
{
  
  if(idx < 0 || idx > int(BoundarySurfacesList.size()-1)) {
    cout << endl;
    cout << "GGeoCylinder::GetBoundaryNormVector:: index " << idx << " is out of boundaries for (0," << BoundarySurfacesList.size()-1 << ") for " 
         << GeoObjectName.Data() << " geometry objet. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  TVector3 NormVect = BoundarySurfacesList[idx]->GetWVector();
  if(BoundarySurfacesList[idx]->GetName().BeginsWith("cylinder bottom-boundary surface")) NormVect = -1*NormVect;
  
  return  NormVect;
  
}
//====================================================================
double  GGeoCylinder::GetVolume()
{
  
  return  2.0*TMath::Pi()*Length*Radius*Thickness;
  
}
//====================================================================
void  GGeoCylinder::Print()
{
  
  TVector3 tmpPosition = (1.0/global->GetDistanceUnit("cm"))*Position;
  TVector3 tmpAngles  = (1.0/global->GetAngleUnit("deg"))*global->GetRotationAngles(Rot);
  TVector3 UVector    = MainSurface->GetUVector();
  TVector3 VVector    = MainSurface->GetVVector();
  TVector3 WVector    = MainSurface->GetWVector();
  
  cout << "  Begin Cylinder Geometry Element" << endl;
  cout << "    Cylinder Name                  " << GeoObjectName.Data() << endl;
  if(LayerName  != TString("")) cout << "    Layer  Name                 " << LayerName.Data() << endl;
  else                          cout << "    Layer  Name                 None" << endl;
  if(SystemName != TString("")) cout << "    System Name                 " << SystemName.Data() << endl;
  else                          cout << "    System Name                 None" << endl;
  cout << "    Length                     " << Length/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << "    Radius                     " << Radius/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << "    Position   (X,Y,Z)         (" << tmpPosition.X() << "," << tmpPosition.Y() << "," << tmpPosition.Z() << ") cm"  << endl;
  cout << "    Rot-Angles (Z,Y,X)         (" << tmpAngles.Z()   << "," << tmpAngles.Y()   << "," << tmpAngles.X()   << ") deg" << endl;
  cout << "    U-vector   (X,Y,Z)         (" << UVector.X()     << "," << UVector.Y()     << "," << UVector.Z()     << ")"     << endl;
  cout << "    V-vector   (X,Y,Z)         (" << VVector.X()     << "," << VVector.Y()     << "," << VVector.Z()     << ")"     << endl;
  cout << "    W-vector   (X,Y,Z)         (" << WVector.X()     << "," << WVector.Y()     << "," << WVector.Z()     << ")"     << endl;
  cout << "    Material                   "  << Material.Data() << endl;
  cout << "    Thickness                  "  << Thickness/global->GetDistanceUnit("um") << " um" << endl;
  cout << "    X/X0                       "  << XOX0*100 << " %" << endl;
  if(IsSensitive) cout << "    IsSensitive                 Yes" << endl;
  else            cout << "    IsSensitive                 No"  << endl;
  if(IsSensitive) { 
    cout << "    V-insensitive (low,high)   (" << VInsensitive.X()*Length/global->GetDistanceUnit("um") << "," << VInsensitive.Y()*Length/global->GetDistanceUnit("um") << ") um" << endl;
    cout << "    Det-efficiency             "  << DetEffic*100 << " %" << endl;
    cout << "    RO-time                    "  << ROtime/global->GetTimeUnit("us") << " us" << endl;
    cout << "    Bkg-Rate                   "  << BkgRate/global->GetRateDensityUnit("MHz/cm2") << " MHz/cm2" << endl;
    cout << "    Resolution (U,V)           (" << ResolutionU/global->GetDistanceUnit("um") << "," << ResolutionV/global->GetDistanceUnit("um") << ") um" << endl;
    if(ResolutionModel != TString("")) {
      cout << "    Resolution Model           " << ResolutionModel.Data() << endl;
      cout << "    Resolution ModelID         " << ResolutionModelID << endl;
    }
  }
  cout << "  End   Cylinder Geometry Element" << endl;
  
  return;
  
}
//====================================================================


