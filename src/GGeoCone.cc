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
#include "include/GGeoCone.h"

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
GGeoCone::GGeoCone(TString   aName,
		   int       aGeometryID,
		   int       aGeoObjectIdx,
		   TVector3  aPosition,
		   TMatrixD  aRot,
		   float     aThickness1,
		   float     aThickness2,
		   TString   aMaterial,
		   bool      aIsSensitive,
		   double    aLength,
		   double    aRadius1,
		   double    aRadius2,
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
				aThickness1,
				aMaterial,
				aIsSensitive,
				aglobal,
				aResolutionU,aResolutionV,
				aDetEffic,aROtime,aBkgRate)
{
  
  GeoObjectType = TString("Cone");
  
  Length       = aLength;
  Radius1      = aRadius1;
  Radius2      = aRadius2;
  Thickness1   = aThickness1;
  Thickness2   = aThickness2;
  VInsensitive = aVInsensitive;
  
  CheckInputs();
  
  FillSurfaces();

  SetAllSurfacesGeoObjectType();
  
}
//====================================================================
GGeoCone::GGeoCone(const GGeoCone& other,TString aName)
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
  Radius1      = other.Radius1;
  Radius2      = other.Radius2;
  Thickness1   = other.Thickness1;
  Thickness2   = other.Thickness2;
  VInsensitive = other.VInsensitive;
  
  ResolutionModel   = other.ResolutionModel;
  ResolutionModelID = other.ResolutionModelID;
  
  CheckInputs();
  
  FillSurfaces();
  
}
//====================================================================
GGeoCone::~GGeoCone() 
{
  
}
//====================================================================
GGeoObject* GGeoCone::clone(TString aName) const
{
 
  return new GGeoCone(*this,aName);
  
}
//====================================================================
void  GGeoCone::CheckInputs()
{
  
  if(Length < 0.0) {
    cout << endl;
    cout << "Length of Cone geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius1 < 0.0) {
    cout << endl;
    cout << "Radius1 of Cone geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius2 < 0.0) {
    cout << endl;
    cout << "Radius2 of Cone geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Thickness1 < 0.0) {
    cout << endl;
    cout << "Thickness1 of Cone geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Thickness2 < 0.0) {
    cout << endl;
    cout << "Thickness2 of Cone geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Cone geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Cone geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() + VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge V-insensive fraction of Cone geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
    
  return;
  
}
//====================================================================
void  GGeoCone::FillSurfaces()
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
  double tmpLength;
  double tmpRadius1;
  double tmpRadius2;
  double tmpRin;
  double tmpRout;
  
  //Creating the main surface object
  NameTmp = TString("Main cone surface for Cone geometry element ") + GeoObjectName;
  MainSurface = new GSurfaceCone(NameTmp,
				 Position,
				 Rot,
				 Length,
				 Radius1,
				 Radius2,
				 VInsensitive,
				 IsSensitive,
				 global);
  
  //Filling the list of 
  BoundarySurfacesList.clear();
  
  //top-boundary
  NameTmp = TString("top-boundary cone surface for Cone geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength  = Length;
  tmpRadius1 = Radius1 + (Thickness1/2);
  tmpRadius2 = Radius2 + (Thickness2/2);
  BoundarySurfacesList.push_back(new GSurfaceCone(NameTmp,
						  tmpPosition,
						  tmpRot,
						  tmpLength,
						  tmpRadius1,
						  tmpRadius2,
						  tmpVInsensitive,
						  false,
						  global));
  
  //bottom-boundary
  NameTmp = TString("bottom-boundary cone surface for Cone geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength  = Length;
  tmpRadius1 = Radius1 - (Thickness1/2);
  tmpRadius2 = Radius2 - (Thickness2/2);
  if(tmpRadius1 > 1.0*global->GetDistanceUnit("um") && tmpRadius2 > 1.0*global->GetDistanceUnit("um")) {
    BoundarySurfacesList.push_back(new GSurfaceCone(NameTmp,
						    tmpPosition,
						    tmpRot,
						    tmpLength,
						    tmpRadius1,
						    tmpRadius2,
						    tmpVInsensitive,
						    false,
						    global));
  }
    
  //fwd disk boundary
  NameTmp = TString("fwd disk boundary surface for Cone geometry element ") + GeoObjectName;
  tmpPosition = Position + (Length/2)*MainWVector;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin    = Radius2 - (Thickness2/2);
  tmpRout   = Radius2 + (Thickness2/2);
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
  NameTmp = TString("bwd disk boundary surface for Cone geometry element ") + GeoObjectName;
  tmpPosition = Position - (Length/2)*MainWVector;
  tmpUVector  = (-1)*MainUVector;
  tmpVVector  = (-1)*MainVVector;
  tmpWVector  = (-1)*MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin    = Radius1 - (Thickness1/2);
  tmpRout   = Radius1 + (Thickness1/2);
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
bool  GGeoCone::IsPointInsideGeometry(TVector3 PosXYZ)
{
  
  for(int kkk=0;kkk<int(BoundarySurfacesList.size());kkk++) {
    double W = (BoundarySurfacesList[kkk]->GetUVWFromXYZ(PosXYZ)).Z();
    if(BoundarySurfacesList[kkk]->GetName().BeginsWith("bottom-boundary cone surface")) W *= -1;
    
    if(W > 0.0) return false;
  }
  
  return true;
  
}
//====================================================================
TVector3  GGeoCone::GetBoundaryNormVector(int idx, TVector3 PosXYZ)
{
  
  if(idx < 0 || idx > int(BoundarySurfacesList.size()-1)) {
    cout << endl;
    cout << "GGeoCone::GetBoundaryNormVector:: index " << idx << " is out of boundaries for (0," << BoundarySurfacesList.size()-1 << ") for " 
         << GeoObjectName.Data() << " geometry objet. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  TVector3 NormVect = BoundarySurfacesList[idx]->GetWVector();
  if(BoundarySurfacesList[idx]->GetName().BeginsWith("bottom-boundary cone surface")) NormVect = -1*NormVect;
  
  return  NormVect;
  
}
//====================================================================
double  GGeoCone::GetVolume()
{
  
  double r1,r2;
  
  //Vol_in
  r1 = Radius1 - 0.5*Thickness1;
  r2 = Radius2 - 0.5*Thickness2;
  double Vol_in = (1.0/3.0)*TMath::Pi()*Length*(pow(r1,2) + r1*r2 + pow(r2,2));
  
  //Vol_out
  r1 = Radius1 + 0.5*Thickness1;
  r2 = Radius2 + 0.5*Thickness2;
  double Vol_out = (1.0/3.0)*TMath::Pi()*Length*(pow(r1,2) + r1*r2 + pow(r2,2));
  
  return  Vol_out - Vol_in;
  
}
//====================================================================
void  GGeoCone::Print()
{
    
  TVector3 tmpPosition = (1.0/global->GetDistanceUnit("cm"))*Position;
  TVector3 tmpAngles  = (1.0/global->GetAngleUnit("deg"))*global->GetRotationAngles(Rot);
  TVector3 UVector    = MainSurface->GetUVector();
  TVector3 VVector    = MainSurface->GetVVector();
  TVector3 WVector    = MainSurface->GetWVector();
  
  double L = sqrt(pow(Radius2 - Radius1,2) + pow(Length,2));
  
  cout << "  Begin Cone Geometry Element" << endl;
  cout << "    Cone Name                  " << GeoObjectName.Data() << endl;
  if(LayerName  != TString("")) cout << "    Layer  Name                 " << LayerName.Data() << endl;
  else                          cout << "    Layer  Name                 None" << endl;
  if(SystemName != TString("")) cout << "    System Name                 " << SystemName.Data() << endl;
  else                          cout << "    System Name                 None" << endl;
  cout << "    Length                     " << Length/global->GetDistanceUnit("cm")  << " cm" << endl;
  cout << "    Radius1                    " << Radius1/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << "    Radius2                    " << Radius2/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << "    Position   (X,Y,Z)         (" << tmpPosition.X() << "," << tmpPosition.Y() << "," << tmpPosition.Z() << ") cm"  << endl;
  cout << "    Rot-Angles (Z,Y,X)         (" << tmpAngles.Z()   << "," << tmpAngles.Y()   << "," << tmpAngles.X()   << ") deg" << endl;
  cout << "    U-vector   (X,Y,Z)         (" << UVector.X()     << "," << UVector.Y()     << "," << UVector.Z()     << ")"     << endl;
  cout << "    V-vector   (X,Y,Z)         (" << VVector.X()     << "," << VVector.Y()     << "," << VVector.Z()     << ")"     << endl;
  cout << "    W-vector   (X,Y,Z)         (" << WVector.X()     << "," << WVector.Y()     << "," << WVector.Z()     << ")"     << endl;
  cout << "    Material                   "  << Material.Data() << endl;
  cout << "    Thickness1                 "  << Thickness1/global->GetDistanceUnit("um") << " um" << endl;
  cout << "    Thickness2                 "  << Thickness2/global->GetDistanceUnit("um") << " um" << endl;
  cout << "    X/X0                       "  << XOX0*100 << " %" << endl;
  if(IsSensitive) cout << "    IsSensitive                 Yes" << endl;
  else            cout << "    IsSensitive                 No"  << endl;
  if(IsSensitive) { 
    cout << "    V-insensitive (low,high)   (" << VInsensitive.X()*L/global->GetDistanceUnit("um") << "," << VInsensitive.Y()*L/global->GetDistanceUnit("um") << ") um" << endl;
    cout << "    Det-efficiency             "  << DetEffic*100 << " %" << endl;
    cout << "    RO-time                    "  << ROtime/global->GetTimeUnit("us") << " us" << endl;
    cout << "    Bkg-Rate                   "  << BkgRate/global->GetRateDensityUnit("MHz/cm2") << " MHz/cm2" << endl;
    cout << "    Resolution (U,V)           (" << ResolutionU/global->GetDistanceUnit("um") << "," << ResolutionV/global->GetDistanceUnit("um") << ") um" << endl;
    if(ResolutionModel != TString("")) {
      cout << "    Resolution Model           " << ResolutionModel.Data() << endl;
      cout << "    Resolution ModelID         " << ResolutionModelID << endl;
    }
  }
  cout << "  End   Cone Geometry Element" << endl;
  
  return;
  
}
//====================================================================


