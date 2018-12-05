/***************************************************************************//**
 * @brief:      
 * @Description: Class for cylinder surfaces
 *
 *
 * @createdby:  PEREZ PEREZ Alejandro <luis_alejandro.perez_perez@iphc.cnrs.fr> at 2017-10-01 14:07:38
 * @copyright:  (c)2017 IPHC - CNRS - Université de Strasbourg. All Rights Reserved.
 * 
 * @License: You are free to use this source files for your own development as long
 *           as it stays in a public research context. You are not allowed to use it
 *           for commercial purpose. You must put this header with laboratory and
 *           authors names in all development based on this library.
 *           When results obtained with this package (Guariguanchi) are communicated, 
 *           you should quote the package name.
 *
 * @lastchange: $Revision$
 *              $Author$
 *              $Date$
 *
 *******************************************************************************/

#ifndef GSurfaceCylinder_h
#define GSurfaceCylinder_h

#include "TMath.h"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TCutG.h>
#include "TStopwatch.h"
#include "include/GSurfaceObject.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GSurfaceCylinder : public GSurfaceObject {

private:
  
public:
  
  double    Length;              // cylinder length
  double    Radius;              // cylinder radius
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total Length
  
  // Constructor
  GSurfaceCylinder(TString   aName,
		   TVector3  aPosition,
		   TMatrixD  aRot,
		   double    Length,
		   double    Radius,
		   TVector2  aVInsensitive,
		   bool      aIsTrackingLayer,
		   GGlobalTools* aglobal);
  
  GSurfaceCylinder(const GSurfaceCylinder& other, TString Name = TString(""));

  // Destructor
  ~GSurfaceCylinder();
  
  //Set of generic functions
  double    GetCylinderLength()  const { return  Length; }
  double    GetCylinderRadius()  const { return  Radius; }
  TVector2  GetVInsensitive()    const { return  VInsensitive; }
  
  //Functions to access internal variables
  void     CheckInputs();
  TVector3  GetInitPositionForIntersection();
  
  void    SetCylinderLength(double aLength)    { Length       = aLength; }
  void    SetCylinderRadius(double aRadius)    { Radius       = aRadius; }
  void    SetVInsensitive(TVector2 vSensitive) { VInsensitive = vSensitive; }
  
  //Set of functions to be defined for the daughter classes
  TVector3  GetUVWFromXYZ(TVector3 PositionXYZ);
  TVector3  GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ, int idx);
  TVector3  GetXYZFromUVW(TVector3 PositionUVW);
  
  GSurfaceObject* clone(TString aName) const;
  bool      IsInMaterial(TVector3 PositionUVW);
  bool      IsInSensitiveMaterial(TVector3 PositionUVW);
  //TVector3  GetNormVector(TVector3 PositionXYZ);
  void      VarySurface(double w1, double w2) { Radius += w1; }
  void      GetVRange(double& Vmin, double& Vmax) {Vmin = -0.5*Length; Vmax = +0.5*Length; }
  void      GetURange(double V, double& Umin, double& Umax) {Umin = -TMath::Pi()*Radius; Umax = +TMath::Pi()*Radius; }
  void      FillSurfaceRepresentation(int color,
				      std::vector<TGraph>   &GraphXY,
				      std::vector<TGraph>   &GraphZY,
				      std::vector<TGraph>   &GraphZX,
				      std::vector<TGraph>   &GraphZR,
				      int linestyle = 1);
  
protected:
  
};

#endif //~ GSurfaceCylinder_h

