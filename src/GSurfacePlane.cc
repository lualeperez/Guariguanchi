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
#include "include/GSurfacePlane.h"

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
GSurfacePlane::GSurfacePlane(TString   aName,
			     TVector3  aPosition,
			     TMatrixD  aRot,
			     TVector2  aUVWidth,
			     TVector2  aUInsensitive,
			     TVector2  aVInsensitive,
			     bool      aIsTrackingLayer,
			     GGlobalTools* aglobal)
                             : GSurfaceObject(aName,
					      aPosition,
					      aRot,
					      aIsTrackingLayer,
					      aglobal)
{
   
  UVWidth         = aUVWidth;
  UInsensitive    = aUInsensitive;
  VInsensitive    = aVInsensitive;
  
  Type = TString("Plane");
  
  CheckInputs();
  
}
//====================================================================
GSurfacePlane::GSurfacePlane(const GSurfacePlane& other,TString Name)
                             : GSurfaceObject(Name,
					      other.Position,
					      other.Rot,
					      other.IsTrackingLayer,
					      other.global) 
{
  
  UVWidth         = other.UVWidth;
  UInsensitive    = other.UInsensitive;
  VInsensitive    = other.VInsensitive;
  Type            = other.Type;
  
  CheckInputs();
  
}
//====================================================================
GSurfacePlane::~GSurfacePlane() 
{
  
}
//====================================================================
GSurfaceObject* GSurfacePlane::clone(TString aName) const
{
 
  return new GSurfacePlane(*this,aName);
  
}
//====================================================================
void  GSurfacePlane::CheckInputs()
{
  
  if(UVWidth.X() < 0.0) {
    cout << endl;
    cout << "Uwidth of Plane surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UVWidth.Y() < 0.0) {
    cout << endl;
    cout << "Vwidth of Plane surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.X() < 0.0 || UInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge U-insensive fraction of Plane surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.Y() < 0.0 || UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge U-insensive fraction of Plane surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Plane surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Plane surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
}
//====================================================================
TVector3  GSurfacePlane::GetInitPositionForIntersection()
{
  
  return Position;
  
}
//====================================================================
TVector3 GSurfacePlane::GetUVWFromXYZ(TVector3 PositionXYZ)
{
  
  TVector3 DeltaPositionXYZ = PositionXYZ - Position;
  
  double U = DeltaPositionXYZ.Dot(UVector);
  double V = DeltaPositionXYZ.Dot(VVector);
  double W = DeltaPositionXYZ.Dot(WVector);
    
  return TVector3(U,V,W);
  
}
//====================================================================
TVector3 GSurfacePlane::GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ,int idx)
{

  //Return the derivative of the local (U,V,W) coordinates w.r.t to global XYZ coordinates
  
  CheckIndex(idx);
  
  return TVector3(UVector(idx),VVector(idx),WVector(idx));
  
}
//====================================================================
TVector3 GSurfacePlane::GetXYZFromUVW(TVector3 PositionUVW)
{
  
  //Return the glonal XYZ coordinates from the local UV coordinates
  
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  double W = PositionUVW.Z();
  
  return U*UVector + V*VVector + W*WVector + Position;
  
}
//====================================================================
bool  GSurfacePlane::IsInMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside plane limits
  
  bool IntersectedMaterial = false;
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  
  double widthU = UVWidth.X();
  double widthV = UVWidth.Y();
    
  if((U >= -0.5*widthU && U <= +0.5*widthU) && 
     (V >= -0.5*widthV && V <= +0.5*widthV)) IntersectedMaterial = true;
  
  return IntersectedMaterial;
  
}
//====================================================================
bool  GSurfacePlane::IsInSensitiveMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside the sensitive surface of a sentivie plane
  
  if(!IsTrackingLayer) return false;
  
  bool IntersectedSensMaterial = false;
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  
  double widthU           =        UVWidth.X();
  double delta_widthU_neg = widthU*UInsensitive.X();
  double delta_widthU_pos = widthU*UInsensitive.Y();
  double withUneg         = 0.5*widthU - delta_widthU_neg;
  double withUpos         = 0.5*widthU - delta_widthU_pos;
    
  double widthV           =        UVWidth.Y();
  double delta_widthV_neg = widthV*VInsensitive.X();
  double delta_widthV_pos = widthV*VInsensitive.Y();
  double withVneg         = 0.5*widthV - delta_widthV_neg;
  double withVpos         = 0.5*widthV - delta_widthV_pos;
    
  if((U >= -withUneg && U <= +withUpos) && 
     (V >= -withVneg && V <= +withVpos)) IntersectedSensMaterial = true;
  
  return IntersectedSensMaterial;
  
}
//====================================================================
TVector3  GSurfacePlane::GetNormVector(TVector3 PositionXYZ)
{
  
  //Returns a normal vector to the PlaneLayer_t object for a given Position in it
  
  return WVector;
  
}
//====================================================================
void  GSurfacePlane::FillSurfaceRepresentation(int color,
					       std::vector<TGraph>   &GraphXY,
					       std::vector<TGraph>   &GraphZY,
					       std::vector<TGraph>   &GraphZX,
					       std::vector<TGraph>   &GraphZR,
					       int linestyle)
{
  
  double WidthU    = UVWidth.X();
  double WidthV    = UVWidth.Y();
  
  TGraph gr_XY;
  gr_XY.SetLineColor(color);
  gr_XY.SetMarkerColor(color);
  gr_XY.SetLineWidth(line_width);
  gr_XY.SetLineStyle(linestyle);
  TGraph gr_ZY;
  gr_ZY.SetLineColor(color);
  gr_ZY.SetMarkerColor(color);
  gr_ZY.SetLineWidth(line_width);
  gr_ZY.SetLineStyle(linestyle);
  TGraph gr_ZX;
  gr_ZX.SetLineColor(color);
  gr_ZX.SetMarkerColor(color);
  gr_ZX.SetLineWidth(line_width);
  gr_ZX.SetLineStyle(linestyle);
  TGraph gr_ZR;
  gr_ZR.SetLineColor(color);
  gr_ZR.SetMarkerColor(color);
  gr_ZR.SetLineWidth(line_width);
  gr_ZR.SetLineStyle(linestyle);
    
  TGraph gr_ZR2;
  gr_ZR2.SetLineColor(color);
  gr_ZR2.SetMarkerColor(color);
  gr_ZR2.SetLineWidth(line_width);
  gr_ZR2.SetLineStyle(linestyle);
    
  TVector3 Pos[4];
  Pos[0] = Position + (-0.5*WidthU)*UVector + (-0.5*WidthV)*VVector;
  Pos[1] = Position + (+0.5*WidthU)*UVector + (-0.5*WidthV)*VVector;
  Pos[2] = Position + (+0.5*WidthU)*UVector + (+0.5*WidthV)*VVector;
  Pos[3] = Position + (-0.5*WidthU)*UVector + (+0.5*WidthV)*VVector;
  for(int i=0;i<4;i++) Pos[i] *= (1.0/global->GetDistanceUnit("cm"));
  
  int counter = 0;
  double x,y,z,r;
  for(int i=0;i<5;i++) {
    int idx = i;
    if(i == 4) idx = 0;
    x = Pos[idx].X();
    y = Pos[idx].Y();
    z = Pos[idx].Z();
    r = sqrt(pow(x,2) + pow(y,2));
    gr_XY.SetPoint(counter,x,y);
    gr_ZY.SetPoint(counter,z,y);
    gr_ZX.SetPoint(counter,z,x);
    
    //gr_ZR.SetPoint(counter,z,r);        
    //if(LadderType == TString("Spiral") || LadderType == TString("Alternating")) gr_ZR2.SetPoint(counter,z,-r);
      
    counter++;
  }
  
  int Npoints = 20;
  counter = 0;
  for(int i=0;i<4;i++) {
    int idx1 = i;
    int idx2 = i+1;
    if(i == 3) {
      idx1 = i;
      idx2 = 0;
    }
    
    TVector3  director = Pos[idx2] - Pos[idx1];
    double d = director.Mag();
    director = director.Unit();
    for(int kkk=0;kkk<Npoints+1;kkk++) {
      double t = (double)kkk/Npoints;
      TVector3 Pos_tmp = Pos[idx1] + (t*d)*director;
      x = Pos_tmp.X();
      y = Pos_tmp.Y();
      z = Pos_tmp.Z();
      r = sqrt(pow(x,2) + pow(y,2));
      gr_ZR.SetPoint(counter,z,r);  
      if(LadderType == TString("Spiral") || LadderType == TString("Alternating") || GeoObjectType == TString("Petal")) gr_ZR2.SetPoint(counter,z,-r);
      
      counter++;
    }
    
  }
    
  GraphXY.push_back(gr_XY);
  GraphZY.push_back(gr_ZY);
  GraphZX.push_back(gr_ZX);
  GraphZR.push_back(gr_ZR);
  if(LadderType == TString("Spiral") || LadderType == TString("Alternating") || GeoObjectType == TString("Petal"))  GraphZR.push_back(gr_ZR2);
  
  return;
  
}
//====================================================================


