/***************************************************************************//**
 * @brief:      
 * @Description: Class for cone geometrical object
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


#ifndef GGeoCone_h
#define GGeoCone_h

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
#include "include/GSurfaceCylinder.h"
#include "include/GSurfaceCone.h"
#include "include/GSurfaceDisk.h"
#include "include/GGlobalTools.h"
#include "include/GGeoObject.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GGeoCone : public GGeoObject {
  
private:

public:
  
  double    Length;              // Cone length
  double    Radius1;             // Cone R1
  double    Radius2;             // Cone R2
  double    Thickness1;          // Thickness at R1
  double    Thickness2;          // Thickness at R2
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total V width
  
  // Constructor
  GGeoCone(TString   aName,
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
	   TVector2  aVInsensitive,
	   GGlobalTools* aglobal,
	   float     aResolutionU = -1.0,
	   float     aResolutionV = -1.0,
	   float     aDetEffic    = 0.0,
	   float     aROtime      = -1,
	   float     aBkgRate     = 0.0);
  
  GGeoCone(const GGeoCone& other, TString Name = TString(""));

  // Destructor
  ~GGeoCone();
  
  //Functions to access internal variables
  double    GetConeLength()         const { return  Length;  }
  double    GetConeRadius1()        const { return  Radius1; }
  double    GetConeRadius2()        const { return  Radius2; }
  double    GetConeThickness1()     const { return  Thickness1; }
  double    GetConeThickness2()     const { return  Thickness2; }
  TVector2  GetVInsensitive()       const { return  VInsensitive; }
  
  //Runctions to set internal variables
  void      GetConeLength(double aLength)        { Length     = aLength;  }
  void      GetConeRadius1(double aRadius1)      { Radius1    = aRadius1; }
  void      GetConeRadius2(double aRadius2)      { Radius2    = aRadius2; }
  void      GetConeThickness1(double aThickness) { Thickness1 = aThickness;  }
  void      GetConeThickness2(double aThickness) { Thickness2 = aThickness;  }
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

#endif //~ GGeoCone_h

