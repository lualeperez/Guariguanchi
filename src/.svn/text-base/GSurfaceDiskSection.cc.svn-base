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
#include "include/GSurfaceDiskSection.h"

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
GSurfaceDiskSection::GSurfaceDiskSection(TString   aName,
					 TVector3  aPosition,
					 TMatrixD  aRot,
					 double    aRin,
					 double    aRout,
					 double    aDeltaPhi,
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
   
  Rin           = aRin;
  Rout          = aRout;
  DeltaPhi      = aDeltaPhi;
  UInsensitive  = aUInsensitive;
  VInsensitive  = aVInsensitive;
  
  Type = TString("DiskSection");
  
  CheckInputs();
  
}
//====================================================================
GSurfaceDiskSection::GSurfaceDiskSection(const GSurfaceDiskSection& other,TString Name)
                                         : GSurfaceObject(Name,
							  other.Position,
							  other.Rot,
							  other.IsTrackingLayer,
							  other.global) 
{
  
  Rin          = other.Rin;
  Rout         = other.Rout;
  DeltaPhi     = other.DeltaPhi;
  UInsensitive = other.UInsensitive;
  VInsensitive = other.VInsensitive;
  Type         = other.Type;
  
  CheckInputs();
  
}
//====================================================================
GSurfaceDiskSection::~GSurfaceDiskSection() 
{
  
}
//====================================================================
GSurfaceObject* GSurfaceDiskSection::clone(TString aName) const
{
 
  return new GSurfaceDiskSection(*this,aName);
  
}
//====================================================================
void  GSurfaceDiskSection::CheckInputs()
{
  
  if(Rin < 0.0) {
    cout << endl;
    cout << "Inner Radius of Disk Section surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Rout < 0.0) {
    cout << endl;
    cout << "Outer Radius of Disk Section surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Rin > Rout) {
    cout << endl;
    cout << "Inner radius if higher than Outer Radius of Disk Section surface object " << Name.Data() << ". Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(DeltaPhi < 0.0) {
    cout << endl;
    cout << "DeltaPhi of Disk Section surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(DeltaPhi > 2*TMath::Pi()) {
    cout << endl;
    cout << "DeltaPhi of Disk Section surface object " << Name.Data() << " is higher than 2*pi. Setting it to 2*pi !!!" << endl;
    cout << endl;
    DeltaPhi = 2*TMath::Pi();
  }
  if(UInsensitive.X() < 0.0 || UInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge angle-wise insensive fraction of Disk Section surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.Y() < 0.0 || UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge angle-wise insensive fraction of Disk Section surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge radius-wise insensive fraction of Disk Section surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge radius-wise insensive fraction of Disk Section surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
}
//====================================================================
TVector3  GSurfaceDiskSection::GetInitPositionForIntersection()
{
  
  return GetXYZFromUVW(TVector3(0,0,0));
  
}
//====================================================================
TVector3 GSurfaceDiskSection::GetUVWFromXYZ(TVector3 PositionXYZ)
{
  
  TVector3 DeltaPositionXYZ = PositionXYZ - Position;
    
  double x = DeltaPositionXYZ.Dot(UVector);
  double y = DeltaPositionXYZ.Dot(VVector);
    
  double angle = global->GetAngle(x,y);
  if(angle > TMath::Pi()) angle -= 2.0*TMath::Pi();
  double R     = sqrt(pow(x,2) + pow(y,2));
    
  double U = R*angle;
  double V = R - 0.5*(Rin + Rout);
  double W = DeltaPositionXYZ.Dot(WVector);
  
  return TVector3(U,V,W);
  
}
//====================================================================
TVector3 GSurfaceDiskSection::GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ, int idx)
{

  //Return the derivative of the local (U,V,W) coordinates w.r.t to global XYZ coordinates
  
  CheckIndex(idx);
  
  TVector3 DeltaPositionXYZ = PositionXYZ - Position;
    
  double x = DeltaPositionXYZ.Dot(UVector);
  double y = DeltaPositionXYZ.Dot(VVector);
  double angle = global->GetAngle(x,y);
  if(angle > TMath::Pi()) angle -= 2.0*TMath::Pi();
  double R     = sqrt(pow(x,2) + pow(y,2));
    
  TVector3 u_der = ((angle*x - y)/R)*UVector + ((angle*y + x)/R)*VVector;
  TVector3 v_der =             (x/R)*UVector +             (y/R)*VVector;
  TVector3 w_der =                   WVector;
  
  return  TVector3(u_der(idx),v_der(idx),w_der(idx));
  
}
//====================================================================
TVector3 GSurfaceDiskSection::GetXYZFromUVW(TVector3 PositionUVW)
{
  
  //Return the glonal XYZ coordinates from the local UV coordinates
  
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  double W = PositionUVW.Z();
  
  double Rave = 0.5*(Rin + Rout);
  
  TVector3 PositionXYZ = Position;
  PositionXYZ += (V + Rave)*TMath::Cos(U/(V + Rave))*UVector;
  PositionXYZ += (V + Rave)*TMath::Sin(U/(V + Rave))*VVector;
  PositionXYZ += W*WVector;
    
  return PositionXYZ;
  
}
//====================================================================
bool  GSurfaceDiskSection::IsInMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside plane limits
  
  bool IntersectedMaterial = false;
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  
  double Vmin,Vmax;
  GetVRange(Vmin,Vmax);
  double Umin,Umax;
  GetURange(V,Umin,Umax);
  
  if((V >= Vmin && V <= Vmax) && 
     (U >= Umin && U <= Umax)) IntersectedMaterial = true;
  
  return IntersectedMaterial;
  
}
//====================================================================
bool  GSurfaceDiskSection::IsInSensitiveMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside the sensitive surface of a sentivie plane
  
  if(!IsTrackingLayer) return false;
  
  bool IntersectedSensMaterial = false;
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  
  double Vmin,Vmax;
  GetVRange(Vmin,Vmax);
  double DeltaV = Vmax - Vmin;
  double Umin,Umax;
  GetURange(V,Umin,Umax);
  double DeltaU = Umax - Umin;
    
  Vmin += DeltaV*VInsensitive.X();
  Vmax += DeltaV*VInsensitive.Y();
  
  if(DeltaPhi < 2.0*TMath::Pi()) {
    Umin += DeltaU*UInsensitive.X();
    Umax += DeltaU*UInsensitive.Y();
  }
     
  if((V >= Vmin && V <= Vmax) && 
     (U >= Umin && U <= Umax)) IntersectedSensMaterial = true;
  
  return IntersectedSensMaterial;
  
}
//====================================================================
TVector3  GSurfaceDiskSection::GetNormVector(TVector3 PositionXYZ)
{
  
  //Returns a normal vector to the PlaneLayer_t object for a given Position in it
  
  return WVector;
  
}
//====================================================================
void  GSurfaceDiskSection::FillSurfaceRepresentation(int color,
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
    double alpha = -0.5*DeltaPhi + (idx+0.5)*DeltaPhi/Npoints;
    double x,y,z,r;
    int sign = 1;
      
    Pos  = Position + Rin*TMath::Cos(alpha)*UVector + Rin*TMath::Sin(alpha)*VVector;
    Pos *= 1.0/global->GetDistanceUnit("cm");
    x = Pos.X();
    y = Pos.Y();
    z = Pos.Z();
    r = sign*sqrt(pow(x,2) + pow(y,2));
    gr_XY1.SetPoint(ipoint,x,y);
    gr_ZY1.SetPoint(ipoint,z,y);
    gr_ZX1.SetPoint(ipoint,z,x);
    gr_ZR1.SetPoint(ipoint,z,r);
    
    Pos  = Position + Rout*TMath::Cos(alpha)*UVector + Rout*TMath::Sin(alpha)*VVector;
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
  GraphZR.push_back(gr_ZR1);
    
  GraphXY.push_back(gr_XY2);
  GraphZY.push_back(gr_ZY2);
  GraphZX.push_back(gr_ZX2);
  GraphZR.push_back(gr_ZR2);
    
  Npoints = 9;
  for(int ipoint=0;ipoint<Npoints+1;ipoint++) {
    int idx = ipoint;
    double alpha = -0.5*DeltaPhi + idx*DeltaPhi/Npoints;
    double x,y,z,r;
    int sign = 1;
      
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
      
    Pos  = Position + Rin*TMath::Cos(alpha)*UVector + Rin*TMath::Sin(alpha)*VVector;
    Pos *= 1.0/global->GetDistanceUnit("cm");
    x = Pos.X();
    y = Pos.Y();
    z = Pos.Z();
    r = sign*sqrt(pow(x,2) + pow(y,2));
    gr_XY_tmp.SetPoint(0,x,y);
    gr_ZY_tmp.SetPoint(0,z,y);
    gr_ZX_tmp.SetPoint(0,z,x);
    gr_ZR_tmp.SetPoint(0,z,r);
      
    Pos  = Position + Rout*TMath::Cos(alpha)*UVector + Rout*TMath::Sin(alpha)*VVector;
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



