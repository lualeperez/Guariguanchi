/***************************************************************************//**
 * @brief:      
 * @Description: Class for cone section surfaces
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
 

#ifndef GSurfaceConeSection_h
#define GSurfaceConeSection_h

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

class GSurfaceConeSection : public GSurfaceObject {

private:
  
public:
  
  double    Length;              // Cone Section length
  double    Radius1;             // Cone Section R1
  double    Radius2;             // Cone Section R2
  double    DeltaPhi;            // Cone Section DeltaPhi
  TVector2  UInsensitive;        // Fraction of insensitive width w.r.t. to total U width
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total V width
  
  // Constructor
  GSurfaceConeSection(TString   aName,
	       TVector3  aPosition,
	       TMatrixD  aRot,
	       double    aLength,
	       double    aRadius1,
	       double    aRadius2,
	       double    aDeltaPhi,
	       TVector2  aUInsensitive,
	       TVector2  aVInsensitive,
	       bool      aIsTrackingLayer,
	       GGlobalTools* aglobal);
  
  GSurfaceConeSection(const GSurfaceConeSection& other, TString Name = TString(""));

  // Destructor
  ~GSurfaceConeSection();
  
  //Set of generic functions
  double    GetConeSectionLength()    const { return  Length;  }
  double    GetConeSectionRadius1()   const { return  Radius1; }
  double    GetConeSectionRadius2()   const { return  Radius2; }
  double    GetConeSectionDeltaPhi()  const { return  DeltaPhi; }
  TVector2  GetUInsensitive()         const { return  VInsensitive; }
  TVector2  GetVInsensitive()         const { return  VInsensitive; }
  
  //Functions to access internal variables
  void      SetConeSectionLength(double aLength)        { Length   = aLength;  }
  void      SetConeSectionRadius1(double aRadius1)      { Radius1  = aRadius1; }
  void      SetConeSectionRadius2(double aRadius2)      { Radius2  = aRadius2; }
  void      SetConeSectionDeltaPhi(double aDeltaPhi)    { DeltaPhi = aDeltaPhi; }
  void      SetUInsensitive(TVector2 uSensitive) { UInsensitive = uSensitive; }
  void      SetVInsensitive(TVector2 vSensitive) { VInsensitive = vSensitive; }
  
  //Set of functions to be defined for the daughter classes
  void      CheckInputs();
  TVector3  GetInitPositionForIntersection();
  
  TVector3  GetUVWFromXYZ(TVector3 PositionXYZ);
  TVector3  GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ, int idx);
  TVector3  GetXYZFromUVW(TVector3 PositionUVW);
  
  GSurfaceObject* clone(TString aName) const;
  bool      IsInMaterial(TVector3 PositionUVW);
  bool      IsInSensitiveMaterial(TVector3 PositionUVW);
  //TVector3  GetNormVector(TVector3 PositionXYZ);
  void      VarySurface(double w1, double w2) {Radius1 += w1; Radius2 += w2; }
  void      GetVRange(double& Vmin, double& Vmax);
  void      GetURange(double V, double& Umin, double& Umax);
  void      FillSurfaceRepresentation(int color,
				      std::vector<TGraph>   &GraphXY,
				      std::vector<TGraph>   &GraphZY,
				      std::vector<TGraph>   &GraphZX,
				      std::vector<TGraph>   &GraphZR,
				      int linestyle = 1);

protected:
  
};

#endif //~ GSurfaceConeSection_h

