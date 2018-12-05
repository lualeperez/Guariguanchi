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
#include "include/GSurfaceCylinder.h"

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
GSurfaceCylinder::GSurfaceCylinder(TString   aName,
				   TVector3  aPosition,
				   TMatrixD  aRot,
				   double    aLength,
				   double    aRadius,
				   TVector2  aVInsensitive,
				   bool      aIsTrackingLayer,
				   GGlobalTools* aglobal) 
                                  : GSurfaceObject(aName,
						   aPosition,
						   aRot,
						   aIsTrackingLayer,
						   aglobal)
{
  
  Length       = aLength;
  Radius       = aRadius;
  VInsensitive = aVInsensitive;
  
  Type = TString("Cylinder");
  
  CheckInputs();
  
}
//====================================================================
GSurfaceCylinder::GSurfaceCylinder(const GSurfaceCylinder& other,TString Name)
                                   : GSurfaceObject(Name,
						    other.Position,
						    other.Rot,
						    other.IsTrackingLayer,
						    other.global) 
{
  
  Length       = other.Length;
  Radius       = other.Radius;
  VInsensitive = other.VInsensitive;
  Type         = other.Type;
  
  CheckInputs();
  
}
//====================================================================
GSurfaceCylinder::~GSurfaceCylinder() 
{
  
}
//====================================================================
GSurfaceObject* GSurfaceCylinder::clone(TString aName) const
{
 
  return new GSurfaceCylinder(*this,aName);
  
}
//====================================================================
void  GSurfaceCylinder::CheckInputs()
{
  
  if(Length < 0.0) {
    cout << endl;
    cout << "Length of Cylinder surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius < 0.0) {
    cout << endl;
    cout << "Radius of Cylinder surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge Lenght-wise insensive fraction of Cylinder surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge Length-wise insensive fraction of Cylinder surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
}
//====================================================================
TVector3 GSurfaceCylinder::GetUVWFromXYZ(TVector3 PositionXYZ)
{
  
  TVector3 DeltaPositionXYZ      = PositionXYZ - Position;
  TVector3 DeltaPositionXYZ_perp = DeltaPositionXYZ - (DeltaPositionXYZ.Dot(WVector))*WVector;
    
  double x = DeltaPositionXYZ.Dot(UVector);
  double y = DeltaPositionXYZ.Dot(VVector);
  double angle = global->GetAngle(x,y);
  if(angle > TMath::Pi()) angle -= 2.0*TMath::Pi();
    
  double U = Radius*angle;
  double V = DeltaPositionXYZ.Dot(WVector);
  double W = DeltaPositionXYZ_perp.Mag() - Radius;

  return TVector3(U,V,W);
  
}
//====================================================================
TVector3  GSurfaceCylinder::GetInitPositionForIntersection()
{
  
  return GetXYZFromUVW(TVector3(0,0,0));
  
}
//====================================================================
TVector3 GSurfaceCylinder::GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ, int idx)
{

  //Return the derivative of the local (U,V,W) coordinates w.r.t to global XYZ coordinates
  
  CheckIndex(idx);
    
  TVector3 DeltaPositionXYZ      = PositionXYZ - Position;
  TVector3 DeltaPositionXYZ_perp = DeltaPositionXYZ - (DeltaPositionXYZ.Dot(WVector))*WVector;
  double x = DeltaPositionXYZ.Dot(UVector);
  double y = DeltaPositionXYZ.Dot(VVector);
  
  double ui        = UVector(idx);
  double vi        = VVector(idx);
  double ni        = WVector(idx);
  double DeltaXi   = DeltaPositionXYZ(idx);
  double Der_alpha = (-y*ui + x*vi)/(pow(x,2) + pow(y,2));
      
  double DerU = Radius*Der_alpha;
  double DerV = ni;
  double DerW = (DeltaXi - ni*(DeltaPositionXYZ.Dot(WVector)))/DeltaPositionXYZ_perp.Mag();
 
  return  TVector3(DerU,DerV,DerW);
  
}
//====================================================================
TVector3 GSurfaceCylinder::GetXYZFromUVW(TVector3 PositionUVW)
{
  
  //Return the glonal XYZ coordinates from the local UV coordinates
  
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  double W = PositionUVW.Z();
    
  return  (Radius + W)*TMath::Cos(U/Radius)*UVector + (Radius + W)*TMath::Sin(U/Radius)*VVector + V*WVector + Position;
  
}
//====================================================================
bool  GSurfaceCylinder::IsInMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside plane limits
  
  bool IntersectedMaterial = false;
  double V = PositionUVW.Y();
  
  if(V >= -0.5*Length && V <= +0.5*Length) IntersectedMaterial = true;
  
  return IntersectedMaterial;
  
}
//====================================================================
bool  GSurfaceCylinder::IsInSensitiveMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside the sensitive surface of a sentivie plane
  
  if(!IsTrackingLayer) return false;
  
  bool IntersectedSensMaterial = false;
  double V = PositionUVW.Y();
    
  double delta_Length_neg = Length*VInsensitive.X();
  double delta_Length_pos = Length*VInsensitive.Y();
  double Lengthneg        = 0.5*Length - delta_Length_neg;
  double Lengthpos        = 0.5*Length - delta_Length_pos;
    
  if(V >= -Lengthneg && V <= +Lengthpos) IntersectedSensMaterial = true;
  
  return IntersectedSensMaterial;
  
}
//====================================================================
// TVector3  GSurfaceCylinder::GetNormVector(TVector3 PositionXYZ)
// {
//   
//   //Returns a normal vector to the PlaneLayer_t object for a given Position in it
//   
//   TVector3 DeltaPositionXYZ      = PositionXYZ - Position;    
//   double x     = DeltaPositionXYZ.Dot(UVector);
//   double y     = DeltaPositionXYZ.Dot(VVector);
//   double angle = global->GetAngle(x,y);
//     
//   return TMath::Cos(angle)*UVector + TMath::Sin(angle)*VVector;
//   
// }
//====================================================================
void  GSurfaceCylinder::FillSurfaceRepresentation(int color,
					          std::vector<TGraph>   &GraphXY,
					          std::vector<TGraph>   &GraphZY,
					          std::vector<TGraph>   &GraphZX,
					          std::vector<TGraph>   &GraphZR,
						  int linestyle)
{
  
  TGraph gr_XY1;
  gr_XY1.SetLineColor(color);
  gr_XY1.SetMarkerColor(color);
  gr_XY1.SetLineWidth(line_width);
  gr_XY1.SetLineStyle(linestyle);
  TGraph gr_ZY1;
  gr_ZY1.SetLineColor(color);
  gr_ZY1.SetMarkerColor(color);
  gr_ZY1.SetLineWidth(line_width);
  gr_ZY1.SetLineStyle(linestyle);
  TGraph gr_ZX1;
  gr_ZX1.SetLineColor(color);
  gr_ZX1.SetMarkerColor(color);
  gr_ZX1.SetLineWidth(line_width);
  gr_ZX1.SetLineStyle(linestyle);
  TGraph gr_ZR1;
  gr_ZR1.SetLineColor(color);
  gr_ZR1.SetMarkerColor(color);
  gr_ZR1.SetLineWidth(line_width);
  gr_ZR1.SetLineStyle(linestyle);
    
  TGraph gr_XY2;
  gr_XY2.SetLineColor(color);
  gr_XY2.SetMarkerColor(color);
  gr_XY2.SetLineWidth(line_width);
  gr_XY2.SetLineStyle(linestyle);
  TGraph gr_ZY2;
  gr_ZY2.SetLineColor(color);
  gr_ZY2.SetMarkerColor(color);
  gr_ZY2.SetLineWidth(line_width);
  gr_ZY2.SetLineStyle(linestyle);
  TGraph gr_ZX2;
  gr_ZX2.SetLineColor(color);
  gr_ZX2.SetMarkerColor(color);
  gr_ZX2.SetLineWidth(line_width);
  gr_ZX2.SetLineStyle(linestyle);
  TGraph gr_ZR2;
  gr_ZR2.SetLineColor(color);
  gr_ZR2.SetMarkerColor(color);
  gr_ZR2.SetLineWidth(line_width);
  gr_ZR2.SetLineStyle(linestyle);
    
  int Npoints = 100;
  TVector3 Pos;
  for(int ipoint=0;ipoint<Npoints+1;ipoint++) {
    int idx = ipoint;
    if(ipoint == Npoints) idx = 0;
    double alpha = (idx+0.5)*2.0*TMath::Pi()/Npoints;
    double x,y,z,r;
    int sign = 1;
    if(alpha >= TMath::Pi()) sign = -1;
      
    //BWD
    Pos  = Position + Radius*TMath::Cos(alpha)*UVector + Radius*TMath::Sin(alpha)*VVector + (-0.5*Length)*WVector;
    Pos *= 1.0/global->GetDistanceUnit("cm");
    x = Pos.X();
    y = Pos.Y();
    z = Pos.Z();
    r = sign*sqrt(pow(x,2) + pow(y,2));
    gr_XY1.SetPoint(ipoint,x,y);
    gr_ZY1.SetPoint(ipoint,z,y);
    gr_ZX1.SetPoint(ipoint,z,x);
    gr_ZR1.SetPoint(ipoint,z,r);
      
    //FWD
    Pos  = Position + Radius*TMath::Cos(alpha)*UVector + Radius*TMath::Sin(alpha)*VVector + (+0.5*Length)*WVector;
    Pos *= 1.0/global->GetDistanceUnit("cm");
    x = Pos.X();
    y = Pos.Y();
    z = Pos.Z();
    r = sign*sqrt(pow(x,2) + pow(y,2));
    gr_XY2.SetPoint(ipoint,x,y);
    gr_ZY2.SetPoint(ipoint,z,y);
    gr_ZX2.SetPoint(ipoint,z,x);
    gr_ZR2.SetPoint(ipoint,z,r);
  }
    
  GraphXY.push_back(gr_XY1);
  GraphZY.push_back(gr_ZY1);
  GraphZX.push_back(gr_ZX1);
  //GraphZR.push_back(gr_ZR1);
    
  GraphXY.push_back(gr_XY2);
  GraphZY.push_back(gr_ZY2);
  GraphZX.push_back(gr_ZX2);
  //GraphZR.push_back(gr_ZR2);
    
  Npoints = 10;
  for(int ipoint=0;ipoint<Npoints;ipoint++) {
    int idx = ipoint;
    double alpha = (idx+0.5)*2.0*TMath::Pi()/Npoints;
    double x,y,z,r;
    int sign = 1;
    if(alpha >= TMath::Pi()) sign = -1;
      
    TGraph gr_XY_tmp;
    gr_XY_tmp.SetLineColor(color);
    gr_XY_tmp.SetMarkerColor(color);
    gr_XY_tmp.SetLineWidth(line_width);
    gr_XY_tmp.SetLineStyle(linestyle);
    TGraph gr_ZY_tmp;
    gr_ZY_tmp.SetLineColor(color);
    gr_ZY_tmp.SetMarkerColor(color);
    gr_ZY_tmp.SetLineWidth(line_width);
    gr_ZY_tmp.SetLineStyle(linestyle);
    TGraph gr_ZX_tmp;
    gr_ZX_tmp.SetLineColor(color);
    gr_ZX_tmp.SetMarkerColor(color);
    gr_ZX_tmp.SetLineWidth(line_width);
    gr_ZX_tmp.SetLineStyle(linestyle);
    TGraph gr_ZR_tmp;
    gr_ZR_tmp.SetLineColor(color);
    gr_ZR_tmp.SetMarkerColor(color);
    gr_ZR_tmp.SetLineWidth(line_width);
    gr_ZR_tmp.SetLineStyle(linestyle);
      
    //BWD
    Pos  = Position + Radius*TMath::Cos(alpha)*UVector + Radius*TMath::Sin(alpha)*VVector + (-0.5*Length)*WVector;
    Pos *= 1.0/global->GetDistanceUnit("cm");
    x = Pos.X();
    y = Pos.Y();
    z = Pos.Z();
    r = sign*sqrt(pow(x,2) + pow(y,2));
    gr_XY_tmp.SetPoint(0,x,y);
    gr_ZY_tmp.SetPoint(0,z,y);
    gr_ZX_tmp.SetPoint(0,z,x);
    gr_ZR_tmp.SetPoint(0,z,r);
      
    //FWD
    Pos  = Position + Radius*TMath::Cos(alpha)*UVector + Radius*TMath::Sin(alpha)*VVector + (+0.5*Length)*WVector;
    Pos *= 1.0/global->GetDistanceUnit("cm");
    x = Pos.X();
    y = Pos.Y();
    z = Pos.Z();
    r = sign*sqrt(pow(x,2) + pow(y,2));
    gr_XY_tmp.SetPoint(1,x,y);
    gr_ZY_tmp.SetPoint(1,z,y);
    gr_ZX_tmp.SetPoint(1,z,x);
    gr_ZR_tmp.SetPoint(1,z,r);
      
    GraphXY.push_back(gr_XY_tmp);
    GraphZY.push_back(gr_ZY_tmp);
    GraphZX.push_back(gr_ZX_tmp);
    GraphZR.push_back(gr_ZR_tmp);
  }
  
  return;
  
}
//====================================================================


