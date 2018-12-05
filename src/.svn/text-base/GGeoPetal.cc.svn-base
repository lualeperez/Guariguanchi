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
#include "include/GGeoPetal.h"

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
GGeoPetal::GGeoPetal(TString   aName,
		     int       aGeometryID,
		     int       aGeoObjectIdx,
		     TVector3  aPosition,
		     TMatrixD  aRot,
		     float     aThickness,
		     TString   aMaterial,
		     bool      aIsSensitive,
		     float     aWBase,
		     float     aWTop,
		     float     aHeight,
		     TVector2  aUInsensitive,
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
  
  GeoObjectType = TString("Petal");
  
  WBase        = aWBase;
  WTop         = aWTop;
  Height       = aHeight;
  UInsensitive = aUInsensitive;
  VInsensitive = aVInsensitive;
  
  CheckInputs();
  
  FillSurfaces();
  
  SetAllSurfacesGeoObjectType();
  
}
//====================================================================
GGeoPetal::GGeoPetal(const GGeoPetal& other,TString aName)
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
  
  WBase        = other.WBase;
  WTop         = other.WTop;
  Height       = other.Height;
  UInsensitive = other.UInsensitive;
  VInsensitive = other.VInsensitive;
  
  ResolutionModel   = other.ResolutionModel;
  ResolutionModelID = other.ResolutionModelID;
  
  CheckInputs();
  
  FillSurfaces();
  
}
//====================================================================
GGeoPetal::~GGeoPetal() 
{
  
}
//====================================================================
GGeoObject* GGeoPetal::clone(TString aName) const
{
 
  return new GGeoPetal(*this,aName);
  
}
//====================================================================
void  GGeoPetal::CheckInputs()
{
  
  if(WBase < 0.0) {
    cout << endl;
    cout << "Base width of Petal geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(WTop < 0.0) {
    cout << endl;
    cout << "Top width of Petal geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Height < 0.0) {
    cout << endl;
    cout << "Height width of Petal geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.X() < 0.0 || UInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge U-insensive fraction of Petal geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.Y() < 0.0 || UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge U-insensive fraction of Petal geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Petal geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Petal geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.X() + UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge U-insensive fraction of Petal geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() + VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge V-insensive fraction of Petal geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
    
  return;
  
}
//====================================================================
void  GGeoPetal::FillSurfaces()
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
  TVector2  tmpUVWidth;
  TVector2  tmpUInsensitive(0.0,0.0);
  TVector2  tmpVInsensitive(0.0,0.0);
  double    tmpWTop;
  double    tmpWBase;
  double    tmpHeight;
  double    alpha = TMath::ATan(0.5*(WTop - WBase)/Height);
  double    h     = sqrt(pow(0.5*(WTop - WBase),2) + pow(Height,2));
  
  //Creating the main surface object
  NameTmp = TString("Main petal surface for Petal geometry element ") + GeoObjectName;
  MainSurface = new GSurfacePetal(NameTmp,
				  Position,
				  Rot,
				  WBase,
				  WTop,
				  Height,
				  UInsensitive,
				  VInsensitive,
				  IsSensitive,
				  global);
  
  //Filling the list of 
  BoundarySurfacesList.clear();
  
  //W-positive boundary
  NameTmp = TString("W-positive petal surface for Petal geometry element ") + GeoObjectName;
  tmpPosition = Position + (Thickness/2)*MainWVector;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpWTop     = WTop;
  tmpWBase    = WBase;
  tmpHeight   = Height;
  BoundarySurfacesList.push_back(new GSurfacePetal(NameTmp,
						   tmpPosition,
						   tmpRot,
						   tmpWBase,
						   tmpWTop,
						   tmpHeight,
						   tmpUInsensitive,
						   tmpVInsensitive,
						   false,
						   global));
  
  //W-negative boundary
  NameTmp = TString("W-negative petal surface for Petal geometry element ") + GeoObjectName;
  tmpPosition = Position - (Thickness/2)*MainWVector;
  tmpUVector  = (-1)*MainUVector;
  tmpVVector  =      MainVVector;
  tmpWVector  = (-1)*MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpWTop     = WTop;
  tmpWBase    = WBase;
  tmpHeight   = Height;
  BoundarySurfacesList.push_back(new GSurfacePetal(NameTmp,
						   tmpPosition,
						   tmpRot,
						   tmpWBase,
						   tmpWTop,
						   tmpHeight,
						   tmpUInsensitive,
						   tmpVInsensitive,
						   false,
						   global));
  
  //U-positive boundary
  NameTmp = TString("U-positive plane surface for Petal geometry element ") + GeoObjectName;
  tmpPosition = Position + 0.5*(WBase + h*TMath::Sin(alpha))*MainUVector + 0.5*(-Height + h*TMath::Cos(alpha))*MainVVector;
  tmpUVector  = -TMath::Sin(alpha)*MainUVector - TMath::Cos(alpha)*MainVVector;
  tmpVVector  = -MainWVector;
  tmpWVector  = +TMath::Cos(alpha)*MainUVector - TMath::Sin(alpha)*MainVVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpUVWidth  = TVector2(h,Thickness);
  BoundarySurfacesList.push_back(new GSurfacePlane(NameTmp,
						   tmpPosition,
						   tmpRot,
						   tmpUVWidth,
						   tmpUInsensitive,
						   tmpVInsensitive,
						   false,
						   global));
  
  //U-negative boundary
  NameTmp = TString("U-negative plane surface for Petal geometry element ") + GeoObjectName;
  tmpPosition = Position - 0.5*(WBase + h*TMath::Sin(alpha))*MainUVector + 0.5*(-Height + h*TMath::Cos(alpha))*MainVVector;
  tmpUVector  = +TMath::Sin(alpha)*MainUVector - TMath::Cos(alpha)*MainVVector;
  tmpVVector  =  MainWVector;
  tmpWVector  = -TMath::Cos(alpha)*MainUVector - TMath::Sin(alpha)*MainVVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpUVWidth  = TVector2(h,Thickness);
  BoundarySurfacesList.push_back(new GSurfacePlane(NameTmp,
						   tmpPosition,
						   tmpRot,
						   tmpUVWidth,
						   tmpUInsensitive,
						   tmpVInsensitive,
						   false,
						   global));
    
  //V-positive boundary
  NameTmp = TString("V-positive plane surface for Petal geometry element ") + GeoObjectName;
  tmpPosition = Position + (Height/2)*MainVVector;
  tmpUVector  = MainWVector;
  tmpVVector  = MainUVector;
  tmpWVector  = MainVVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpUVWidth  = TVector2(Thickness,WTop);
  BoundarySurfacesList.push_back(new GSurfacePlane(NameTmp,
						   tmpPosition,
						   tmpRot,
						   tmpUVWidth,
						   tmpUInsensitive,
						   tmpVInsensitive,
						   false,
						   global));
  
  //V-negative boundary
  NameTmp = TString("V-negative plane surface for Petal geometry element ") + GeoObjectName;
  tmpPosition = Position - (Height/2)*MainVVector;
  tmpUVector  = (-1)*MainWVector;
  tmpVVector  = (-1)*MainUVector;
  tmpWVector  = (-1)*MainVVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpUVWidth  = TVector2(Thickness,WBase);
  BoundarySurfacesList.push_back(new GSurfacePlane(NameTmp,
						   tmpPosition,
						   tmpRot,
						   tmpUVWidth,
						   tmpUInsensitive,
						   tmpVInsensitive,
						   false,
						   global));
  
  return;
  
}
//====================================================================
bool  GGeoPetal::IsPointInsideGeometry(TVector3 PosXYZ)
{
  
  for(int kkk=0;kkk<int(BoundarySurfacesList.size());kkk++) {
    if((BoundarySurfacesList[kkk]->GetUVWFromXYZ(PosXYZ)).Z() > 0.0) return false;
  }
  
  return true;
  
}
//====================================================================
double  GGeoPetal::GetVolume()
{
  
  return  0.5*Height*(WBase + WTop)*Thickness;
  
}
//====================================================================
void  GGeoPetal::Print()
{
    
  TVector3 tmpPosition = (1.0/global->GetDistanceUnit("cm"))*Position;
  TVector3 tmpAngles  = (1.0/global->GetAngleUnit("deg"))*global->GetRotationAngles(Rot);
  TVector3 UVector    = MainSurface->GetUVector();
  TVector3 VVector    = MainSurface->GetVVector();
  TVector3 WVector    = MainSurface->GetWVector();
  
  cout << "  Begin Petal Geometry Element" << endl;
  cout << "    Petal Name                  " << GeoObjectName.Data() << endl;
  if(LayerName  != TString("")) cout << "    Layer  Name                 " << LayerName.Data() << endl;
  else                          cout << "    Layer  Name                 None" << endl;
  if(SystemName != TString("")) cout << "    System Name                 " << SystemName.Data() << endl;
  else                          cout << "    System Name                 None" << endl;
  cout << "    bottomWidth                " << WBase/global->GetDistanceUnit("cm")  << " cm" << endl;
  cout << "    topWidth                   " << WTop/global->GetDistanceUnit("cm")   << " cm" << endl;
  cout << "    Height                     " << Height/global->GetDistanceUnit("cm") << " cm" << endl;
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
    cout << "    U-insensitive (low,high)   (" << UInsensitive.X()*TMath::Max(WBase,WTop)/global->GetDistanceUnit("um") << "," << UInsensitive.Y()*TMath::Max(WBase,WTop)/global->GetDistanceUnit("um") << ") um" << endl;
    cout << "    V-insensitive (low,high)   (" << VInsensitive.X()*Height/global->GetDistanceUnit("um") << "," << VInsensitive.Y()*Height/global->GetDistanceUnit("um") << ") um" << endl;
    cout << "    Det-efficiency             "  << DetEffic*100 << " %" << endl;
    cout << "    RO-time                    "  << ROtime/global->GetTimeUnit("us") << " us" << endl;
    cout << "    Bkg-Rate                   "  << BkgRate/global->GetRateDensityUnit("MHz/cm2") << " MHz/cm2" << endl;
    cout << "    Resolution (U,V)           (" << ResolutionU/global->GetDistanceUnit("um") << "," << ResolutionV/global->GetDistanceUnit("um") << ") um" << endl;
    if(ResolutionModel != TString("")) {
      cout << "    Resolution Model           " << ResolutionModel.Data() << endl;
      cout << "    Resolution ModelID         " << ResolutionModelID << endl;
    }
  }
  cout << "  End   Petal Geometry Element" << endl;
  
  return;
  
}
//====================================================================


