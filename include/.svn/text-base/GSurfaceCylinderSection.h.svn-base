/***************************************************************************//**
 * @brief:      
 * @Description: Class for cylinder section surfaces
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

#ifndef GSurfaceCylinderSection_h
#define GSurfaceCylinderSection_h

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

class GSurfaceCylinderSection : public GSurfaceObject {

private:
  
public:
  
  double    Length;              // cylinder section length
  double    Radius;              // cylinder section radius
  double    DeltaPhi;            // cylinder section angular extension
  TVector2  UInsensitive;        // Fraction of insensitive width w.r.t. to total angular extension
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total Length
  
  // Constructor
  GSurfaceCylinderSection(TString   aName,
			  TVector3  aPosition,
			  TMatrixD  aRot,
			  double    aLength,
			  double    aRadius,
			  double    aDeltaPhi,
			  TVector2  aUInsensitive,
			  TVector2  aVInsensitive,
			  bool      aIsTrackingLayer,
			  GGlobalTools* aglobal);
  
  GSurfaceCylinderSection(const GSurfaceCylinderSection& other, TString Name = TString(""));

  // Destructor
  ~GSurfaceCylinderSection();
  
  //Set of generic functions
  double    GetCylinderSectionLength()    const { return  Length; }
  double    GetCylinderSectionRadius()    const { return  Radius; }
  double    GetCylinderSectionDeltaPhi()  const { return  DeltaPhi; }
  TVector2  GetUInsensitive()             const { return  UInsensitive; }
  TVector2  GetVInsensitive()             const { return  VInsensitive; }
  
  //Functions to access internal variables
  void     CheckInputs();
  TVector3  GetInitPositionForIntersection();
  
  void    SetCylinderSectionLength(double aLength)     { Length       = aLength; }
  void    SetCylinderSectionRadius(double aRadius)     { Radius       = aRadius; }
  void    SetCylinderSectionDeltaPhi(double aDeltaPhi) { DeltaPhi     = aDeltaPhi; }
  void    SetUInsensitive(TVector2 uSensitive)         { UInsensitive = uSensitive; }
  void    SetVInsensitive(TVector2 vSensitive)         { VInsensitive = vSensitive; }
  
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
  void      GetURange(double V, double& Umin, double& Umax) {Umin = -0.5*DeltaPhi*Radius; Umax = +0.5*DeltaPhi*Radius; }
  void      FillSurfaceRepresentation(int color,
				      std::vector<TGraph>   &GraphXY,
				      std::vector<TGraph>   &GraphZY,
				      std::vector<TGraph>   &GraphZX,
				      std::vector<TGraph>   &GraphZR,
				      int linestyle = 1);
  
protected:
  
};

#endif //~ GSurfaceCylinderSection_h

