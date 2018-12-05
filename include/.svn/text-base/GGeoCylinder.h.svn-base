/***************************************************************************//**
 * @brief:      
 * @Description: Class for cylinder geometrical object
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
 

#ifndef GGeoCylinder_h
#define GGeoCylinder_h

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

class GGeoCylinder : public GGeoObject {
  
private:

public:
  
  double    Length;              // cylinder length
  double    Radius;              // cylinder radius
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total V width
  
  // Constructor
  GGeoCylinder(TString   aName,
	       int       aGeometryID,
	       int       aGeoObjectIdx,
               TVector3  aPosition,
               TMatrixD  aRot,
	       float     aThickness,
	       TString   aMaterial,
	       bool      aIsSensitive,
	       float     aLength,
	       float     aRadius,
	       TVector2  aVInsensitive,
	       GGlobalTools* aglobal,
	       float     aResolutionU = -1.0,
	       float     aResolutionV = -1.0,
	       float     aDetEffic    = 0.0,
	       float     aROtime      = -1,
	       float     aBkgRate     = 0.0);

  GGeoCylinder(const GGeoCylinder& other, TString Name = TString(""));
  
  // Destructor
  ~GGeoCylinder();
  
  //Functions to access internal variables
  double    GetLength()        const { return  Length; }
  double    GetRadius()        const { return  Radius; }
  TVector2  GetVInsensitive()  const { return  VInsensitive; }
  
  //Runctions to set internal variables
  void      SetLength(double aLength)                { Length       = aLength; }
  void      SetRadius(double aRadius)                { Radius       = aRadius; }
  void      SetVInsensitive(TVector2 aVInsensitive)  { VInsensitive = aVInsensitive; }
  
  //Set of functions to be defined for the daughter classes
  //void  GetThicknesses(double& w1, double& w2) {w1 = Thickness; w2 = 0.0; }
  GGeoObject* clone(TString aName) const;
  double    GetVolume();
  void      CheckInputs();
  void      FillSurfaces();
  bool      IsPointInsideGeometry(TVector3 PosXYZ);
  TVector3  GetBoundaryNormVector(int idx, TVector3 PosXYZ);
  void      Print();

protected:

};

#endif //~ GGeoCylinder_h

