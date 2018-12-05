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
#include "include/GSurfacePetal.h"

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
GSurfacePetal::GSurfacePetal(TString   aName,
			     TVector3  aPosition,
			     TMatrixD  aRot,
			     double    aWBase,
			     double    aWTop,
			     double    aHeight,
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
   
  WBase        = aWBase;
  WTop         = aWTop;
  Height       = aHeight;
  UInsensitive = aUInsensitive;
  VInsensitive = aVInsensitive;
  
  Type = TString("Petal");
  
  CheckInputs();
  
}
//====================================================================
GSurfacePetal::GSurfacePetal(const GSurfacePetal& other,TString Name)
                             : GSurfaceObject(Name,
					      other.Position,
					      other.Rot,
					      other.IsTrackingLayer,
					      other.global) 
{
  
  WBase         = other.WBase;
  WTop          = other.WTop;
  Height        = other.Height;
  UInsensitive  = other.UInsensitive;
  VInsensitive  = other.VInsensitive;
  Type          = other.Type;
  
  CheckInputs();
  
}
//====================================================================
GSurfacePetal::~GSurfacePetal() 
{
  
}
//====================================================================
GSurfaceObject* GSurfacePetal::clone(TString aName) const
{
 
  return new GSurfacePetal(*this,aName);
  
}
//====================================================================
void  GSurfacePetal::CheckInputs()
{
  
  if(WBase < 0.0) {
    cout << endl;
    cout << "Base width of Petal surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(WTop < 0.0) {
    cout << endl;
    cout << "Top width of Petal surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Height < 0.0) {
    cout << endl;
    cout << "Height of Petal surface object " << Name.Data() << " is negative. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.X() < 0.0 || UInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge U-insensive fraction of Petal surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(UInsensitive.Y() < 0.0 || UInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge U-insensive fraction of Petal surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.X() < 0.0 || VInsensitive.X() > 1.0) {
    cout << endl;
    cout << "Low edge V-insensive fraction of Petal surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(VInsensitive.Y() < 0.0 || VInsensitive.Y() > 1.0) {
    cout << endl;
    cout << "High edge V-insensive fraction of Petal surface object " << Name.Data() << " is not in [0,1] range. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
}
//====================================================================
TVector3  GSurfacePetal::GetInitPositionForIntersection()
{
  
  return Position;
  
}
//====================================================================
TVector3 GSurfacePetal::GetUVWFromXYZ(TVector3 PositionXYZ)
{
  
  TVector3 DeltaPositionXYZ = PositionXYZ - Position;
  
  double U = DeltaPositionXYZ.Dot(UVector);
  double V = DeltaPositionXYZ.Dot(VVector);
  double W = DeltaPositionXYZ.Dot(WVector);
    
  return TVector3(U,V,W);
  
}
//====================================================================
TVector3 GSurfacePetal::GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ,int idx)
{

  //Return the derivative of the local (U,V,W) coordinates w.r.t to global XYZ coordinates
  
  CheckIndex(idx);
  
  return TVector3(UVector(idx),VVector(idx),WVector(idx));
  
}
//====================================================================
TVector3 GSurfacePetal::GetXYZFromUVW(TVector3 PositionUVW)
{
  
  //Return the glonal XYZ coordinates from the local UV coordinates
  
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  double W = PositionUVW.Z();
  
  return U*UVector + V*VVector + W*WVector + Position;
  
}
//====================================================================
bool  GSurfacePetal::IsInMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside plane limits
  
  bool IntersectedMaterial = false;
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  
  double RV[2];
  double RU[2];
  GetVRange(RV[0],RV[1]);
  GetURange(V,RU[0],RU[1]);
  
  if((U >= RU[0] && U <= RU[1]) && (V >= RV[0] && V <= RV[1])) IntersectedMaterial = true;
  
#if 0
  cout << endl;
  cout << "(U,V) = (" << U/global->GetDistanceUnit("mm") << "," << V/global->GetDistanceUnit("mm") << ") mm, " 
       << "RangeU = (" << RU[0]/global->GetDistanceUnit("mm") << "," << RU[1]/global->GetDistanceUnit("mm") << ") mm, "
       << "RangeV = (" << RV[0]/global->GetDistanceUnit("mm") << "," << RV[1]/global->GetDistanceUnit("mm") << ") mm";
       
  if(U >= RU[0] && U <= RU[1]) cout << " (is Inside U-range) ";
  else                         cout << " (is not Inside U-range) ";
  
  if(V >= RV[0] && V <= RV[1]) cout << " (is Inside V-range) ";
  else                         cout << " (is not Inside V-range) ";
  
  if(IntersectedMaterial) cout << " (Is inside)";
  else                    cout << " (Is Not inside)";
  cout << endl;
  cout << endl;
#endif
  
  return IntersectedMaterial;
  
}
//====================================================================
bool  GSurfacePetal::IsInSensitiveMaterial(TVector3 PositionUVW)
{
  
  //Checks is a given local UV position is inside the sensitive surface of a sentivie plane
  
  if(!IsTrackingLayer) return false;
  
  bool IntersectedSensMaterial = false;
  double U = PositionUVW.X();
  double V = PositionUVW.Y();
  
  double RV[2];
  double RU[2];
  GetVRange(RV[0],RV[1]);
  GetURange(V,RU[0],RU[1]);
  
  double delta_widthU_neg = (RU[1] - RU[0])*UInsensitive.X();
  double delta_widthU_pos = (RU[1] - RU[0])*UInsensitive.Y();
  RU[0] += delta_widthU_neg;
  RU[1] -= delta_widthU_pos;
  
  double delta_widthV_neg = (RV[1] - RV[0])*VInsensitive.X();
  double delta_widthV_pos = (RV[1] - RV[0])*VInsensitive.Y();
  RV[0] += delta_widthV_neg;
  RV[1] -= delta_widthV_pos;
    
  if((U >= RU[0] && U <= RU[1]) && 
     (V >= RV[0] && V <= RV[1])) IntersectedSensMaterial = true;
  
  return IntersectedSensMaterial;
  
}
//====================================================================
TVector3  GSurfacePetal::GetNormVector(TVector3 PositionXYZ)
{
  
  //Returns a normal vector to the PlaneLayer_t object for a given Position in it
  
  return WVector;
  
}
//====================================================================
void  GSurfacePetal::FillSurfaceRepresentation(int color,
					       std::vector<TGraph>   &GraphXY,
					       std::vector<TGraph>   &GraphZY,
					       std::vector<TGraph>   &GraphZX,
					       std::vector<TGraph>   &GraphZR,
					       int linestyle)
{
  
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
  Pos[0] = Position + (-0.5*WBase)*UVector + (-0.5*Height)*VVector;
  Pos[1] = Position + (+0.5*WBase)*UVector + (-0.5*Height)*VVector;
  Pos[2] = Position + (+0.5*WTop)*UVector  + (+0.5*Height)*VVector;
  Pos[3] = Position + (-0.5*WTop)*UVector  + (+0.5*Height)*VVector;
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
    //gr_ZR2.SetPoint(counter,z,-r);
    
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
    
    //cout << "Pos["<<idx2<<"] - Pos["<<idx1<<"]" << endl;
    //cout << "Pos["<<idx1<<"] = " << Pos[idx1].X()/global->GetDistanceUnit("cm") << "," << Pos[idx1].Y()/global->GetDistanceUnit("cm") << "," << Pos[idx1].Z()/global->GetDistanceUnit("cm") << ") cm" << endl;
    //cout << "Pos["<<idx2<<"] = " << Pos[idx2].X()/global->GetDistanceUnit("cm") << "," << Pos[idx2].Y()/global->GetDistanceUnit("cm") << "," << Pos[idx2].Z()/global->GetDistanceUnit("cm") << ") cm" << endl;
    //cout << "d = " << d/global->GetDistanceUnit("cm") << " cm" << endl;
    //cout << "director  = " << director.X() << "," << director.Y() << "," << director.Z() << ")" << endl;
    
    for(int kkk=0;kkk<Npoints+1;kkk++) {
      double t = (double)kkk/Npoints;
      TVector3 Pos_tmp = Pos[idx1] + (t*d)*director;
      x = Pos_tmp.X();
      y = Pos_tmp.Y();
      z = Pos_tmp.Z();
      r = sqrt(pow(x,2) + pow(y,2));
      gr_ZR.SetPoint(counter,z,r);  
      gr_ZR2.SetPoint(counter,z,-r);
      
      //cout << "kkk = " << kkk << ", " 
      //     << "t = " << t << ", "
	//   << "t*d = " << t*d/global->GetDistanceUnit("cm") << " cm"
	//   << "Pos_tmp = (" << Pos_tmp.X()/global->GetDistanceUnit("cm") << "," << Pos_tmp.Y()/global->GetDistanceUnit("cm") << "," << Pos_tmp.Z()/global->GetDistanceUnit("cm") << ") cm, "
        //   << "r = " << r/global->GetDistanceUnit("cm") << " cm" 
	//   << endl;
      
      counter++;
    }
    
  }

  GraphXY.push_back(gr_XY);
  GraphZY.push_back(gr_ZY);
  GraphZX.push_back(gr_ZX);
  
  GraphZR.push_back(gr_ZR);  
  GraphZR.push_back(gr_ZR2);
  
  return;
  
}
//====================================================================
void  GSurfacePetal::GetURange(double V, double& Umin, double& Umax)
{
  
  double a = (WTop - WBase)/(2*Height);
  double b = 0.5*(WTop + WBase)/2;
  
  Umin = -a*V - b;
  Umax = +a*V + b;
  
  return;
  
}
//====================================================================


