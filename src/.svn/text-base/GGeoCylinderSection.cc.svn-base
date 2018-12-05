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
#include "include/GGeoCylinderSection.h"

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
GGeoCylinderSection::GGeoCylinderSection(TString   aName,
					 int       aGeometryID,
					 int       aGeoObjectIdx,
					 TVector3  aPosition,
					 TMatrixD  aRot,
					 float     aThickness,
					 TString   aMaterial,
					 bool      aIsSensitive,
					 float     aLength,
					 float     aRadius,
					 float     aDeltaPhi,
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
  
  GeoObjectType = TString("CylinderSection");
  
  Length       = aLength;
  Radius       = aRadius;
  DeltaPhi     = aDeltaPhi;
  UInsensitive = aUInsensitive;
  VInsensitive = aVInsensitive;
  
  CheckInputs();
  
  FillSurfaces();
  
  SetAllSurfacesGeoObjectType();
  
}
//====================================================================
GGeoCylinderSection::GGeoCylinderSection(const GGeoCylinderSection& other,TString aName)
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
GGeoCylinderSection::~GGeoCylinderSection() 
{
  
}
//====================================================================
GGeoObject* GGeoCylinderSection::clone(TString aName) const
{
 
  return new GGeoCylinderSection(*this,aName);
  
}
//====================================================================
void  GGeoCylinderSection::CheckInputs()
{
  
  if(Length < 0.0) {
    cout << endl;
    cout << "Length of Cylinder Section geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius < 0.0) {
    cout << endl;
    cout << "Radius of Cylinder Section geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(DeltaPhi < 0.0) {
    cout << endl;
    cout << "DeltaPhi of Cylinder Section geometry object " << GeoObjectName.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(DeltaPhi > 2.0*TMath::Pi()) {
    cout << endl;
    cout << "DeltaPhi of Cylinder Section geometry object " << GeoObjectName.Data() << " is higher than 2*pi. Setting it to 2*pi !!!" << endl;
    cout << endl;
    DeltaPhi = 2.0*TMath::Pi();
  }
  
  if(UInsensitive.X() < 0.0 || UInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge U-insensive fraction of Cylinder Section geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.Y() < 0.0 || UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge U-insensive fraction of Cylinder Section geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.X() + UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge U-insensive fraction of Cylinder Section geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Cylinder Section geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Cylinder Section geometry object " << GeoObjectName.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() + VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High + Low edge V-insensive fraction of Cylinder Section geometry object " << GeoObjectName.Data() << " is higher than 1. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
    
  return;
  
}
//====================================================================
void  GGeoCylinderSection::FillSurfaces()
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
  TVector2  tmpUInsensitive(0.0,0.0);
  TVector2  tmpVInsensitive(0.0,0.0);
  
  //Creating the main surface object
  NameTmp = TString("Main cylinder section surface for Cylinder Section geometry element ") + GeoObjectName;
  MainSurface = new GSurfaceCylinderSection(NameTmp,
					    Position,
					    Rot,
					    Length,
					    Radius,
					    DeltaPhi,
					    UInsensitive,
					    VInsensitive,
					    IsSensitive,
					    global);
  
  //Filling the list of 
  BoundarySurfacesList.clear();
  
  //top-boundary
  NameTmp = TString("cylinder section top-boundary surface for Cylinder Section geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength = Length;
  tmpRadius = Radius + Thickness/2.0;
  BoundarySurfacesList.push_back(new GSurfaceCylinderSection(NameTmp,
							     tmpPosition,
							     tmpRot,
							     tmpLength,
							     tmpRadius,
							     DeltaPhi,
							     tmpUInsensitive,
							     tmpVInsensitive,
							     false,
							     global));
  
  //bottom-boundary
  NameTmp = TString("cylinder section bottom-boundary surface for Cylinder Section geometry element ") + GeoObjectName;
  tmpPosition = Position;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpLength = Length;
  tmpRadius = Radius - Thickness/2.0;
  if(tmpRadius > 1.0*global->GetDistanceUnit("um")) {
    BoundarySurfacesList.push_back(new GSurfaceCylinderSection(NameTmp,
							       tmpPosition,
							       tmpRot,
							       tmpLength,
							       tmpRadius,
							       DeltaPhi,
							       tmpUInsensitive,
							       tmpVInsensitive,
							       false,
							       global));
  }
  
  //fwd disk boundary
  NameTmp = TString("fwd-disk section boundary surface for Cylinder Section geometry element ") + GeoObjectName;
  tmpPosition = Position + (Length/2)*MainWVector;
  tmpUVector  = MainUVector;
  tmpVVector  = MainVVector;
  tmpWVector  = MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin    = Radius - (Thickness/2);
  tmpRout   = Radius + (Thickness/2);
  if(tmpRin < 0.0) tmpRin = 0.0;
  BoundarySurfacesList.push_back(new GSurfaceDiskSection(NameTmp,
							 tmpPosition,
							 tmpRot,
							 tmpRin,
							 tmpRout,
							 DeltaPhi,
							 tmpUInsensitive,
							 tmpVInsensitive,
							 false,
							 global));
  
  //bwd disk boundary
  NameTmp = TString("bwd-disk section boundary surface for Cylinder Section geometry element ") + GeoObjectName;
  tmpPosition = Position - (Length/2)*MainWVector;
  tmpUVector  =      MainUVector;
  tmpVVector  = (-1)*MainVVector;
  tmpWVector  = (-1)*MainWVector;
  global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
  tmpRin    = Radius - (Thickness/2);
  tmpRout   = Radius + (Thickness/2);
  if(tmpRin < 0.0) tmpRin = 0.0;
  BoundarySurfacesList.push_back(new GSurfaceDiskSection(NameTmp,
							 tmpPosition,
							 tmpRot,
							 tmpRin,
							 tmpRout,
							 DeltaPhi,
							 tmpUInsensitive,
							 tmpVInsensitive,
							 false,
							 global));
  
  if(DeltaPhi < 2*TMath::Pi()) {
    //Positive DeltaPhi-end plane boundary
    NameTmp = TString("Positive DeltaPhi-end plane boundary surface for Cylinder section geometry element ") + GeoObjectName;
    tmpPosition = Position + (Radius*TMath::Cos(0.5*DeltaPhi))*MainUVector - (Radius*TMath::Sin(0.5*DeltaPhi))*MainVVector;
    tmpUVector  = -MainWVector;
    tmpVVector  =  TMath::Cos(0.5*DeltaPhi)*MainUVector - TMath::Sin(0.5*DeltaPhi)*MainVVector;
    tmpWVector  = -TMath::Sin(0.5*DeltaPhi)*MainUVector - TMath::Cos(0.5*DeltaPhi)*MainVVector;
    global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
    BoundarySurfacesList.push_back(new GSurfacePlane(NameTmp,
						     tmpPosition,
						     tmpRot,
						     TVector2(Length,Thickness),
						     tmpUInsensitive,
						     tmpVInsensitive,
						     false,
						     global));
    
    //Negative DeltaPhi-end plane boundary
    NameTmp = TString("Negative DeltaPhi-end plane boundary surface for Cylinder section geometry element ") + GeoObjectName;
    tmpPosition = Position + (Radius*TMath::Cos(0.5*DeltaPhi))*MainUVector + (Radius*TMath::Sin(0.5*DeltaPhi))*MainVVector;
    tmpUVector  =  MainWVector;
    tmpVVector  =  TMath::Cos(0.5*DeltaPhi)*MainUVector + TMath::Sin(0.5*DeltaPhi)*MainVVector;
    tmpWVector  = -TMath::Sin(0.5*DeltaPhi)*MainUVector + TMath::Cos(0.5*DeltaPhi)*MainVVector;
    global->GetRotationMatrixFromUnitVects(tmpUVector,tmpVVector,tmpWVector,tmpRot);
    BoundarySurfacesList.push_back(new GSurfacePlane(NameTmp,
						     tmpPosition,
						     tmpRot,
						     TVector2(Length,Thickness),
						     tmpUInsensitive,
						     tmpVInsensitive,
						     false,
						     global));
  }
  
  return;
  
}
//====================================================================
bool  GGeoCylinderSection::IsPointInsideGeometry(TVector3 PosXYZ)
{
  
  for(int kkk=0;kkk<int(BoundarySurfacesList.size());kkk++) {
    double W = (BoundarySurfacesList[kkk]->GetUVWFromXYZ(PosXYZ)).Z();
    if(BoundarySurfacesList[kkk]->GetName().BeginsWith("cylinder section bottom-boundary surface")) W *= -1;
    
    if(W > 0.0) return false;
  }
  
  return true;
  
}
//====================================================================
TVector3  GGeoCylinderSection::GetBoundaryNormVector(int idx, TVector3 PosXYZ)
{
  
  if(idx < 0 || idx > int(BoundarySurfacesList.size()-1)) {
    cout << endl;
    cout << "GGeoCylinderSection::GetBoundaryNormVector:: index " << idx << " is out of boundaries for (0," << BoundarySurfacesList.size()-1 << ") for " 
         << GeoObjectName.Data() << " geometry objet. Check your inputs. Exiting now !!!" 
	 << endl;
    cout << endl;
    assert(false);
  }
  
  TVector3 NormVect = BoundarySurfacesList[idx]->GetWVector();
  if(BoundarySurfacesList[idx]->GetName().BeginsWith("cylinder section bottom-boundary surface")) NormVect = -1*NormVect;
  
  return  NormVect;
  
}
//====================================================================
double  GGeoCylinderSection::GetVolume()
{
  
  return  2.0*TMath::Pi()*Length*Radius*Thickness*(DeltaPhi/(2.0*TMath::Pi()));
  
}
//====================================================================
void  GGeoCylinderSection::Print()
{
  
  TVector3 tmpPosition = (1.0/global->GetDistanceUnit("cm"))*Position;
  TVector3 tmpAngles  = (1.0/global->GetAngleUnit("deg"))*global->GetRotationAngles(Rot);
  TVector3 UVector    = MainSurface->GetUVector();
  TVector3 VVector    = MainSurface->GetVVector();
  TVector3 WVector    = MainSurface->GetWVector();
  
  cout << "  Begin Cylinder Section Geometry Element" << endl;
  cout << "    Cylinder Name                  " << GeoObjectName.Data() << endl;
  if(LayerName  != TString("")) cout << "    Layer  Name                 " << LayerName.Data() << endl;
  else                          cout << "    Layer  Name                 None" << endl;
  if(SystemName != TString("")) cout << "    System Name                 " << SystemName.Data() << endl;
  else                          cout << "    System Name                 None" << endl;
  cout << "    Length                     " << Length/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << "    Radius                     " << Radius/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << "    DeltaPhi                   " << DeltaPhi/global->GetAngleUnit("deg") << " deg" << endl;
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
  cout << "  End   Cylinder Section Geometry Element" << endl;
  
  return;
  
}
//====================================================================


