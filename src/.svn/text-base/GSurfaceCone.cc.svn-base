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
#include "include/GSurfaceCone.h"

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
GSurfaceCone::GSurfaceCone(TString   aName,
			   TVector3  aPosition,
			   TMatrixD  aRot,
			   double    aLength,
			   double    aRadius1,
			   double    aRadius2,
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
  Radius1      = aRadius1;
  Radius2      = aRadius2;
  VInsensitive = aVInsensitive;
  
  Type = TString("Cone");
  
  CheckInputs();
  
}
//====================================================================
GSurfaceCone::GSurfaceCone(const GSurfaceCone& other,TString Name)
                           : GSurfaceObject(Name,
					    other.Position,
					    other.Rot,
					    other.IsTrackingLayer,
					    other.global) 
{
  
  Length        = other.Length;
  Radius1       = other.Radius1;
  Radius2       = other.Radius2;
  VInsensitive  = other.VInsensitive;
  Type          = other.Type;
  
  CheckInputs();
  
}
//====================================================================
GSurfaceCone::~GSurfaceCone() 
{
  
}
//====================================================================
GSurfaceObject* GSurfaceCone::clone(TString aName) const
{
 
  return new GSurfaceCone(*this,aName);
  
}
//====================================================================
void  GSurfaceCone::CheckInputs()
{
  
  if(Length < 0.0) {
    cout << endl;
    cout << "Length of Cone surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius1 < 0.0) {
    cout << endl;
    cout << "Radius1 of Cone surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Radius2 < 0.0) {
    cout << endl;
    cout << "Radius2 of Cone surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Cone surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Cone surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
}
//====================================================================
TVector3  GSurfaceCone::GetInitPositionForIntersection()
{
  
  return  GetXYZFromUVW(TVector3(0,0,0));
  
}
//====================================================================
TVector3 GSurfaceCone::GetUVWFromXYZ(TVector3 PositionXYZ)
{
  
  TVector3 DeltaPositionXYZ = PositionXYZ - Position;
  
  double x = DeltaPositionXYZ.Dot(UVector);
  double y = DeltaPositionXYZ.Dot(VVector);
  double z = DeltaPositionXYZ.Dot(WVector);
    
  double beta = global->GetAngle(x,y);
  if(beta > TMath::Pi()) beta -= 2.0*TMath::Pi();
  double cosBeta  = TMath::Cos(beta);
  double sinBeta  = TMath::Sin(beta);
    
  double DeltaR   = Radius2 - Radius1;
  double cosAlpha = Length/sqrt(pow(Length,2) + pow(DeltaR,2));
  double sinAlpha = DeltaR/sqrt(pow(Length,2) + pow(DeltaR,2));
  double Rave     = 0.5*(Radius1 + Radius2);
  double a        = DeltaR/Length;
    
  double V = sinAlpha*(x*cosBeta + y*sinBeta - Rave) + cosAlpha*z;
  double W = cosAlpha*(x*cosBeta + y*sinBeta - Rave) - sinAlpha*z;
  double R = a*cosAlpha*V + Rave;
  double U = R*beta;
    
  return TVector3(U,V,W);
  
}
//====================================================================
TVector3 GSurfaceCone::GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ, int idx)
{

  //Return the derivative of the local (U,V,W) coordinates w.r.t to global XYZ coordinates
  
  CheckIndex(idx);
  
  TVector3 DeltaPositionXYZ = PositionXYZ - Position;
  double x = DeltaPositionXYZ.Dot(UVector);
  double y = DeltaPositionXYZ.Dot(VVector);
  double z = DeltaPositionXYZ.Dot(WVector);
    
  double beta     = global->GetAngle(x,y);
  double cosBeta  = TMath::Cos(beta);
  double sinBeta  = TMath::Sin(beta);
    
  double DeltaR   = Radius2 - Radius1;
  double cosAlpha = Length/sqrt(pow(Length,2) + pow(DeltaR,2));
  double sinAlpha = DeltaR/sqrt(pow(Length,2) + pow(DeltaR,2));
  double Rave     = 0.5*(Radius2 + Radius1);
  double a        = DeltaR/Length;
  double V        = sinAlpha*(x*cosBeta + y*sinBeta - Rave) + cosAlpha*z;
  double R        = a*cosAlpha*V + Rave;
    
  double ui      = UVector(idx);
  double vi      = VVector(idx);
  double ni      = WVector(idx);
  double DerBeta = (-y*ui + x*vi)/(pow(x,2) + pow(y,2));
  
  double DerV    = sinAlpha*(ui*cosBeta + vi*sinBeta + DerBeta*(-x*sinBeta + y*cosBeta)) + ni*cosAlpha;
  double DerW    = cosAlpha*(ui*cosBeta + vi*sinBeta + DerBeta*(-x*sinBeta + y*cosBeta)) - ni*sinAlpha;
  double DerR    = a*cosAlpha*DerV;
  double DerU    = R*DerBeta + DerR*beta;
  
  return  TVector3(DerU,DerV,DerW);
      
}
//====================================================================
TVector3 GSurfaceCone::GetXYZFromUVW(TVector3 PositionUVW)
{
  
  //Return the glonal XYZ coordinates from the local UV coordinates
  
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  double W = PositionUVW.Z();
  
  double DeltaR   = Radius2 - Radius1;
  double cosAlpha = Length/sqrt(pow(Length,2) + pow(DeltaR,2));
  double sinAlpha = DeltaR/sqrt(pow(Length,2) + pow(DeltaR,2));
  double a        = DeltaR/Length;
  double Rave     = 0.5*(Radius1 + Radius2);
    
  double z = V*cosAlpha - W*sinAlpha + Position.Z();
  double R = a*cosAlpha*V + Rave;
  double x = TMath::Abs(V*sinAlpha + W*cosAlpha + Rave)*TMath::Cos(U/R) + Position.X();
  double y = TMath::Abs(V*sinAlpha + W*cosAlpha + Rave)*TMath::Sin(U/R) + Position.Y();
    
  TVector3 PositionXYZ(x,y,z);
    
  return PositionXYZ;
  
}
//====================================================================
bool  GSurfaceCone::IsInMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside plane limits
  
  bool IntersectedMaterial = false;
  double V = PositionUVW.Y();
  
  double widthV = sqrt(pow(Length,2) + pow(Radius2 - Radius1,2));
  if(V >= -0.5*widthV && V <= +0.5*widthV) IntersectedMaterial = true;
  
  return IntersectedMaterial;
  
}
//====================================================================
bool  GSurfaceCone::IsInSensitiveMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside the sensitive surface of a sentivie plane
  
  if(!IsTrackingLayer) return false;
  
  bool IntersectedSensMaterial = false;
  double V = PositionUVW.Y();
  
  double widthV           = sqrt(pow(Length,2) + pow(Radius2 - Radius1,2));
  double delta_widthV_neg = widthV*VInsensitive.X();
  double delta_widthV_pos = widthV*VInsensitive.Y();
  double withVneg         = 0.5*widthV - delta_widthV_neg;
  double withVpos         = 0.5*widthV - delta_widthV_pos;
    
  if(V >= -withVneg && V <= +withVpos) IntersectedSensMaterial = true;
  
  return IntersectedSensMaterial;
  
}
//====================================================================
// TVector3  GSurfaceCone::GetNormVector(TVector3 PositionXYZ)
// {
//   
//   //Returns a normal vector to the PlaneLayer_t object for a given Position in it
//   
//   TVector3 DeltaPositionXYZ      = PositionXYZ - Position;    
//   double x     = DeltaPositionXYZ.Dot(UVector);
//   double y     = DeltaPositionXYZ.Dot(VVector);
//   double angle = global->GetAngle(x,y);
//     
//   double cosAlpha =              Length/sqrt(pow(Length,2) + pow(Radius2 - Radius1,2));
//   double sinAlpha = (Radius2 - Radius1)/sqrt(pow(Length,2) + pow(Radius2 - Radius1,2));
//     
//   TVector3 xVect(1.0,0.0,0.0);
//   TVector3 yVect(1.0,0.0,0.0);
//   TVector3 zVect(1.0,0.0,0.0);
//   
//   return TMath::Cos(angle)*cosAlpha*xVect + TMath::Sin(angle)*cosAlpha*yVect - sinAlpha*zVect;
//   
// }
//====================================================================
void  GSurfaceCone::FillSurfaceRepresentation(int color,
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
    Pos  = Position + Radius1*TMath::Cos(alpha)*UVector + Radius1*TMath::Sin(alpha)*VVector + (-0.5*Length)*WVector;
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
    Pos  = Position + Radius2*TMath::Cos(alpha)*UVector + Radius2*TMath::Sin(alpha)*VVector + (+0.5*Length)*WVector;
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
    Pos  = Position + Radius1*TMath::Cos(alpha)*UVector + Radius1*TMath::Sin(alpha)*VVector + (-0.5*Length)*WVector;
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
    Pos  = Position + Radius2*TMath::Cos(alpha)*UVector + Radius2*TMath::Sin(alpha)*VVector + (+0.5*Length)*WVector;
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
void  GSurfaceCone::GetVRange(double& Vmin, double& Vmax)
{
  
  double L = sqrt(pow(Length,2) + pow(Radius2 - Radius1,2));
  
  Vmin = -0.5*L;
  Vmax = +0.5*L;
  
  return;
  
}
//====================================================================
void  GSurfaceCone::GetURange(double V, double& Umin, double& Umax)
{
  
  double a   = (Radius2 - Radius1)/Length;
  double b   = 0.5*(Radius2 + Radius1);
  double cos = Length/sqrt(pow(Length,2) + pow(Radius2 - Radius1,2));
  double r   = a*V*cos + b;
  
  Umin = -TMath::Pi()*r;
  Umax = +TMath::Pi()*r;
  
  return;
  
}
//====================================================================

