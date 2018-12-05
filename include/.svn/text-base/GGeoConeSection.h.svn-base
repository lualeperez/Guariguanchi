/***************************************************************************//**
 * @brief:      
 * @Description: Class for cone section geometrical object
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


#ifndef GGeoConeSection_h
#define GGeoConeSection_h

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
#include "include/GSurfacePlane.h"
#include "include/GSurfacePetal.h"
#include "include/GSurfaceCylinder.h"
#include "include/GSurfaceCylinderSection.h"
#include "include/GSurfaceCone.h"
#include "include/GSurfaceConeSection.h"
#include "include/GSurfaceDisk.h"
#include "include/GSurfaceDiskSection.h"
#include "include/GGlobalTools.h"
#include "include/GGeoObject.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GGeoConeSection : public GGeoObject {
  
private:

public:
  
  double    Length;              // Cone Section length
  double    Radius1;             // Cone Section R1
  double    Radius2;             // Cone Section R2
  double    Thickness1;          // Thickness at R1
  double    Thickness2;          // Thickness at R2
  double    DeltaPhi;
  TVector2  UInsensitive;        // Fraction of insensitive width w.r.t. to total U width
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total V width
  
  // Constructor
  GGeoConeSection(TString   aName,
		  int       aGeometryID,
		  int       aGeoObjectIdx,
		  TVector3  aPosition,
		  TMatrixD  aRot,
		  float     aThickness1,
		  float     aThickness2,
		  TString   aMaterial,
		  bool      aIsSensitive,
		  double    aLength,
		  double    aRadius1,
		  double    aRadius2,
		  double    aDeltaPhi,
		  TVector2  aUInsensitive,
		  TVector2  aVInsensitive,
		  GGlobalTools* aglobal,
		  float     aResolutionU = -1.0,
		  float     aResolutionV = -1.0,
		  float     aDetEffic    = 0.0,
		  float     aROtime      = -1,
		  float     aBkgRate     = 0.0);
  
  GGeoConeSection(const GGeoConeSection& other, TString Name = TString(""));

  // Destructor
  ~GGeoConeSection();
  
  //Functions to access internal variables
  double    GetConeSectionLength()         const { return  Length;  }
  double    GetConeSectionRadius1()        const { return  Radius1; }
  double    GetConeSectionRadius2()        const { return  Radius2; }
  double    GetConeSectionThickness1()     const { return  Thickness1; }
  double    GetConeSectionThickness2()     const { return  Thickness2; }
  double    GetConeSectionDeltaPhi()       const { return  DeltaPhi; }
  TVector2  GetUInsensitive()       const { return  UInsensitive; }
  TVector2  GetVInsensitive()       const { return  VInsensitive; }
  
  //Runctions to set internal variables
  void      GetConeSectionLength(double aLength)        { Length     = aLength;  }
  void      GetConeSectionRadius1(double aRadius1)      { Radius1    = aRadius1; }
  void      GetConeSectionRadius2(double aRadius2)      { Radius2    = aRadius2; }
  void      GetConeSectionThickness1(double aThickness) { Thickness1 = aThickness;  }
  void      GetConeSectionThickness2(double aThickness) { Thickness2 = aThickness;  }
  void      GetConeSectionDeltaPhi(double aDeltaPhi)    { DeltaPhi   = aDeltaPhi;  }
  void      SetUInsensitive(TVector2 aUInsensitive)  { UInsensitive = aUInsensitive; }
  void      SetVInsensitive(TVector2 aVInsensitive)  { VInsensitive = aVInsensitive; }
  
  //Set of functions to be defined for the daughter classes
  void  GetThicknesses(double& w1, double& w2) {w1 = Thickness1; w2 = Thickness2; }
  GGeoObject* clone(TString aName) const;
  double    GetVolume();
  void      CheckInputs();
  void      FillSurfaces();
  bool      IsPointInsideGeometry(TVector3 PosXYZ);
  TVector3  GetBoundaryNormVector(int idx, TVector3 PosXYZ);
  void      Print();

protected:

};

#endif //~ GGeoConeSection_h

