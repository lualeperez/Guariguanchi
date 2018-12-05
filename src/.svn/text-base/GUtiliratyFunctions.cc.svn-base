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
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include "include/GUtiliratyFunctions.h"

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

//====================================================================
double    IntersectCoordinates(GSurfaceObject* aSurface,
			       GTrajectory* Trajectory,
			       GGlobalTools* global,
			       double sinit,
			       bool DoPositiveS,
			       double scut,
			       bool verbose)
{
  
  //Calculate the value of the trajectory s parameter for the intersection between the track and a Surface

  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  TVector3 Intersect(0.0,0.0,0.0);
  double limit = 1.0*global->GetDistanceUnit("nm");
  
  const int MaxIterations(MaxIterations_intersection);
  int counter = 0;
  double s0,s1;
  double delta0;
  double Derdelta0;
  double deltaODerdelta;
  
  s1 = 0.0;
  
  int Npush = 1;
  if(TMath::Abs(scut) > 1.0*global->GetUnit("um")) Npush = 10;

  counter = 0;
  double w = 0.0;
  if(DoPositiveS) {
    if(sinit == Dummy_value) s0 = (InitPosition - aSurface->GetInitPositionForIntersection()).Mag();
    else                     s0 = sinit;
  }
  else {
    if(sinit == Dummy_value) s0 = -(InitPosition - aSurface->GetInitPositionForIntersection()).Mag();
    else                     s0 = sinit;
  }
    
  TVector3 pos         = Trajectory->GetTrueTrajectoryCoordinates(s0);
  TVector3 Derpos      = Trajectory->GetTrueTrajectoryUnitMon(s0);
  TVector3 posUVW      = aSurface->GetUVWFromXYZ(pos);
  double DerposW = 0.0;
  for(int i=0;i<3;i++) DerposW += aSurface->GetDerUVWFromXYZ_Analytical(pos,i).Z()*Derpos(i);

  delta0    = w  - posUVW.Z();
  Derdelta0 =    - DerposW;

  if(TMath::Abs(Derdelta0) < 1.0e-8) deltaODerdelta = -1.0;
  else                                deltaODerdelta = delta0/Derdelta0;

  s1 = s0 - deltaODerdelta;
  
  if(DoPositiveS) {
    if(s1 < scut) s1 = scut + TMath::Abs(s1 - scut)*Npush;
  }
  else {
    if(s1 > scut) s1 = scut - TMath::Abs(s1 - scut)*Npush;
  }
  posUVW = aSurface->GetUVWFromXYZ(Trajectory->GetTrueTrajectoryCoordinates(s1));

  if(verbose) {
    cout << "s1 = " << s1/global->GetUnit("cm") <<  ", cm, delta s = " << (s0 - s1)/global->GetUnit("um") << " um, W = " << posUVW.Z()/global->GetUnit("um") << " um, interation " << counter << endl;
  }
  
  while((TMath::Abs(posUVW.Z()) > limit || TMath::Abs(s1 - s0) > limit) && counter <= MaxIterations) {
    if(verbose) {
      cout << "s1 = " << s1/global->GetUnit("cm") <<  ", cm, delta s = " << (s0 - s1)/global->GetUnit("um") << " um, W = " << posUVW.Z()/global->GetUnit("um") << " um, interation " << counter << endl;
    }
    
    s0 = s1;
    pos         = Trajectory->GetTrueTrajectoryCoordinates(s0);
    Derpos      = Trajectory->GetTrueTrajectoryUnitMon(s0);
    posUVW      = aSurface->GetUVWFromXYZ(pos);
    DerposW = 0.0;
    for(int i=0;i<3;i++) DerposW += aSurface->GetDerUVWFromXYZ_Analytical(pos,i).Z()*Derpos(i);

    delta0    = w  - posUVW.Z();
    Derdelta0 =    - DerposW;
      
    if(TMath::Abs(Derdelta0) < 1.0e-10) deltaODerdelta = -1.0;
    else                                deltaODerdelta = delta0/Derdelta0;

    s1 = s0 - deltaODerdelta;
    
    if(DoPositiveS) {
      if(s1 < scut) s1 = scut + TMath::Abs(s1 - scut)*Npush;
    }
    else {
      if(s1 > scut) s1 = scut - TMath::Abs(s1 - scut)*Npush;
    }
    posUVW = aSurface->GetUVWFromXYZ(Trajectory->GetTrueTrajectoryCoordinates(s1));

    counter++;
  }

  if(TMath::Abs(posUVW.Z()) > limit || TMath::Abs(s1 - s0) > limit) {
    if(verbose) {
      cout << "WARNNING in function GUtiliratyFunctions::IntersectCoordinates, iteration difference " << TMath::Abs(posUVW.Z())/global->GetDistanceUnit("cm") 
           << " cm is bigger than limit " << limit/global->GetDistanceUnit("um") << " um after " << counter << " iterations. Result not precise!." 
	   << endl;
    }
    s1 = Dummy_value;
  }

  if(verbose) {
    cout << "s1 = " << s1/global->GetDistanceUnit("cm") << " cm" << endl;
  }
    
  return s1;
  
}
//====================================================================
bool  DoGeometryElementsOverlap(GGeoObject* aGeoObject1,
				GGeoObject* aGeoObject2,
				GGlobalTools* global,
				double Precision)
{
  
  //This function check for overlaps between two geometry elements
  //It will make a scan on the W local coordinates from -w/2 to w/2,
  //where w is the thickness of the geometry element.
  //For each value of W it will built a new surface object located at 
  //the shifted position. Then checks if there is an intersection
  //between the shifted surface objects.
  
  double w1, w2;
  w1 = w2 = 0.0;
  aGeoObject2->GetThicknesses(w1,w2);
  
  int Nsteps_w = int(TMath::Max(w1,w2)/Precision);
  for(int iw=0;iw<Nsteps_w;iw++) {
    double wtmp1 = 0.0;
    double wtmp2 = 0.0;
    wtmp1 = -0.5*w1 + (iw+0.5)*w1/Nsteps_w;
    wtmp2 = -0.5*w2 + (iw+0.5)*w2/Nsteps_w;
    
    GSurfaceObject* surface_test_2 = aGeoObject2->GetMainSurface()->clone(TString("tmp2"));
    surface_test_2->VarySurface(wtmp1,wtmp2);
      
    if(DoSurfaceIntersectsObject(aGeoObject1,surface_test_2,global,Precision)) {
      delete  surface_test_2;
      return  true;
    }
      
    delete  surface_test_2;
  }
  
  return false;
  
}
//====================================================================
bool  DoSurfaceIntersectsObject(GGeoObject*     aGeoObject,
				GSurfaceObject* aSurfaceObject,
				GGlobalTools* global,
				double Precision)
{
  
  //This function scans the local coordinates of surface object aSurfaceObject1
  //transform to the global coordinate system
  //then transform to local coordinate system of surface object aSurfaceObject2
  //and checks if the point is within the limits of this one
  
  double RV[2];
  aSurfaceObject->GetVRange(RV[0],RV[1]);
  int Nsteps_v = int((RV[1] - RV[0])/Precision); //number of bins is defined by GeoCheckPrecision parameter
  for(int iv=0;iv<Nsteps_v;iv++) {
    double V = RV[0] + (iv+0.5)*(RV[1] - RV[0])/Nsteps_v;
    
    double RU[2];
    aSurfaceObject->GetURange(V,RU[0],RU[1]);
    
    int Nsteps_u = int((RU[1] - RU[0])/Precision); //number of bins is defined by GeoCheckPrecision parameter
    for(int iu=0;iu<Nsteps_u;iu++) {
      double U = RU[0] + (iu+0.5)*(RU[1] - RU[0])/Nsteps_u;
      
      TVector3 PosUVW1(U,V,0);
      TVector3 PosXYZ  = aSurfaceObject->GetXYZFromUVW(PosUVW1); //Transform local coordinates (U,V,0) of plane APlane1 to global reference frame
      
      if(aGeoObject->IsPointInsideGeometry(PosXYZ)) {
	cout << endl;
	cout << "Intersection found at point (x,y,z,r) = (" 
	     << PosXYZ.X()/global->GetDistanceUnit("cm") << "," 
	     << PosXYZ.Y()/global->GetDistanceUnit("cm") << "," 
	     << PosXYZ.Z()/global->GetDistanceUnit("cm") << ","
	     << sqrt(pow(PosXYZ.X(),2) + pow(PosXYZ.Y(),2))/global->GetDistanceUnit("cm") << ") cm"
	     << endl;
	return true;
      }
      
    }
  }
  
  return false;
  
}
//====================================================================
