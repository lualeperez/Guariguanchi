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
#include "include/GGeoDisk.h"

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
GGeoDisk::GGeoDisk(TString   aName,
		   int       aGeometryID,
		   int       aGeoObjectIdx,
		   TVector3  aPosition,
		   TMatrixD  aRot,
		   float     aThickness,
		   TString   aMaterial,
		   bool      aIsSensitive,
		   double    aRin,
		   double    aRout,
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
  
  GeoObjectType = TString("Disk");
  
  Rin          = aRin;
  Rout         = aRout;
  VInsensitive = aVInsensitive;
  
  CheckInputs();
  
  FillSurfaces();

  SetAllSurfacesGeoObjectType();
  
}
//====================================================================
GGeoDisk::GGeoDisk(const GGeoDisk& other,TString aName)
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
  
  Rin          = other.Rin;
  Rout         = other.Rout;
  VInsensitive = other.VInsensitive;
  
  ResolutionModel   = other.ResolutionModel;
  ResolutionModelID = other.ResolutionModelID;
  
  CheckInputs();
  
  FillSurfaces();
  
}
//====================================================================
GGeoDisk::~GGeoDisk() 
{
  
}
//====================================================================
GGeoObject* GGeoDisk::clone(TString aName) const
{
 
  return new GGeoDisk(*this,aName);
  
}
//====================================================================
void  GGeoDisk::CheckInputs()
{
  
  if(Rin < 0.0) {
    cout << endl;
    cout << "Inner Radius of Disk geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Rout < 0.0) {
    cout << endl;
    cout << "Outer Radius of Disk geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Rout <= Rin) {
    cout << endl;
    cout << "Outer Radius is smaller or equal to inner Radius of Disk geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Disk geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Disk geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() + VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge V-insensive fraction of Disk geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
    
  return;
  
}
//====================================================================
void  GGeoDisk::FillSurfaces()
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
  TVector2  tmpVInsensitive(0.0,0.0);
  double    tmpLength;
  double    tmpRadius;
  double    tmpRin;
  double    tmpRout;
  
  //Creating the main surface object
  NameTmp = TString("Main disk surface for Disk geometry element ") + GeoObjectName;
  MainSurface = new GSurfaceDisk(NameTmp,
				 Position,
				 Rot,
				 Rin,
				 Rout,
				 VInsensitive,
				 IsSensitive,
				 global);
  
  //Filling the list of 
  BoundarySurfacesList.clear();
  
  //top-positive boundary
  NameTmp = TString("top-boundary disk surface for Disk geometry element ") + GeoObjectName;
  tmpPosition = Position + (Thickness/2)*MainWVector;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin      = Rin;
  tmpRout     = Rout;
  BoundarySurfacesList.push_back(new GSurfaceDisk(NameTmp,
						  tmpPosition,
						  tmpRot,
						  tmpRin,
						  tmpRout,
						  tmpVInsensitive,
						  false,
						  global));
  
  //bottom-positive boundary
  NameTmp = TString("bottom-boundary disk surface for Disk geometry element ") + GeoObjectName;
  tmpPosition = Position - (Thickness/2)*MainWVector;
  tmpUVector  = (-1)*MainUVector;
  tmpVVector  = (-1)*MainVVector;
  tmpWVector  = (-1)*MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin      = Rin;
  tmpRout     = Rout;
  BoundarySurfacesList.push_back(new GSurfaceDisk(NameTmp,
						  tmpPosition,
						  tmpRot,
						  tmpRin,
						  tmpRout,
						  tmpVInsensitive,
						  false,
						  global));
  
  //Outer cylinder
  NameTmp = TString("outer cylinder boundary surface for Disk geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength   = Thickness;
  tmpRadius   = Rout;
  BoundarySurfacesList.push_back(new GSurfaceCylinder(NameTmp,
						      tmpPosition,
						      tmpRot,
						      tmpLength,
						      tmpRadius,
						      tmpVInsensitive,
						      false,
						      global));
  
  //Inner cylinder
  NameTmp = TString("inner cylinder boundary surface for Disk geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength   = Thickness;
  tmpRadius   = Rin;
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
  
  return;
  
}
//====================================================================
bool  GGeoDisk::IsPointInsideGeometry(TVector3 PosXYZ)
{
  
  for(int kkk=0;kkk<int(BoundarySurfacesList.size());kkk++) {
    double W = (BoundarySurfacesList[kkk]->GetUVWFromXYZ(PosXYZ)).Z();
    if(BoundarySurfacesList[kkk]->GetName().BeginsWith("inner cylinder boundary")) W *= -1;
    
    if(W > 0.0) return false;
  }
  
  return true;
  
}
//====================================================================
TVector3  GGeoDisk::GetBoundaryNormVector(int idx, TVector3 PosXYZ)
{
  
  if(idx < 0 || idx > int(BoundarySurfacesList.size()-1)) {
    cout << endl;
    cout << "GGeoDisk::GetBoundaryNormVector:: index " << idx << " is out of boundaries for (0," << BoundarySurfacesList.size()-1 << ") for " 
         << GeoObjectName.Data() << " geometry objet. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  TVector3 NormVect = BoundarySurfacesList[idx]->GetWVector();
  if(BoundarySurfacesList[idx]->GetName().BeginsWith("inner cylinder boundary")) NormVect = -1*NormVect;
  
  return  NormVect;
  
}
//====================================================================
double  GGeoDisk::GetVolume()
{
  
  return  TMath::Pi()*Thickness*(pow(Rout,2) - pow(Rin,2));
  
}
//====================================================================
void  GGeoDisk::Print()
{
    
  TVector3 tmpPosition = (1.0/global->GetDistanceUnit("cm"))*Position;
  TVector3 tmpAngles  = (1.0/global->GetAngleUnit("deg"))*global->GetRotationAngles(Rot);
  TVector3 UVector    = MainSurface->GetUVector();
  TVector3 VVector    = MainSurface->GetVVector();
  TVector3 WVector    = MainSurface->GetWVector();
  
  cout << "  Begin Disk Geometry Element" << endl;
  cout << "    Disk Name                  " << GeoObjectName.Data() << endl;
  if(LayerName  != TString("")) cout << "    Layer  Name                 " << LayerName.Data() << endl;
  else                          cout << "    Layer  Name                 None" << endl;
  if(SystemName != TString("")) cout << "    System Name                 " << SystemName.Data() << endl;
  else                          cout << "    System Name                 None" << endl;
  cout << "    Rin                        " << Rin/global->GetDistanceUnit("cm")  << " cm" << endl;
  cout << "    Rout                       " << Rout/global->GetDistanceUnit("cm") << " cm" << endl;
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
    cout << "    V-insensitive (low,high)   (" << VInsensitive.X()*(Rout - Rin)/global->GetDistanceUnit("um") << "," << VInsensitive.Y()*(Rout - Rin)/global->GetDistanceUnit("um") << ") um" << endl;
    cout << "    Det-efficiency             "  << DetEffic*100 << " %" << endl;
    cout << "    RO-time                    "  << ROtime/global->GetTimeUnit("us") << " us" << endl;
    cout << "    Bkg-Rate                   "  << BkgRate/global->GetRateDensityUnit("MHz/cm2") << " MHz/cm2" << endl;
    cout << "    Resolution (U,V)           (" << ResolutionU/global->GetDistanceUnit("um") << "," << ResolutionV/global->GetDistanceUnit("um") << ") um" << endl;
    if(ResolutionModel != TString("")) {
      cout << "    Resolution Model           " << ResolutionModel.Data() << endl;
      cout << "    Resolution ModelID         " << ResolutionModelID << endl;
    }
  }
  cout << "  End   Disk Geometry Element" << endl;
  
  return;
  
}
//====================================================================


