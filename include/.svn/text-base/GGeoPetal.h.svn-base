/***************************************************************************//**
 * @brief:      
 * @Description: Class for petal geometrical object
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
 

#ifndef GGeoPetal_h
#define GGeoPetal_h

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
#include "include/GGlobalTools.h"
#include "include/GGeoObject.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GGeoPetal : public GGeoObject {
  
private:

public:
  
  double    WBase;               // Width at base of petal
  double    WTop;                // Width at top  of petal
  double    Height;              // High of petal
  
  TVector2  UVWidth;             // Plane withds in local U and V coordinates.
  TVector2  UInsensitive;        // Fraction of insensitive width w.r.t. to total U width
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total V width
  
  // Constructor
  GGeoPetal(TString   aName,
	    int       aGeometryID,
	    int       aGeoObjectIdx,
            TVector3  aPosition,
            TMatrixD  aRot,
	    float     aThickness,
	    TString   aMaterial,
	    bool      aIsSensitive,
	    float     aWBase,
	    float     aWTop,
	    float     aHeight,
	    TVector2  aUInsensitive,
	    TVector2  aVInsensitive,
	    GGlobalTools* aglobal,
	    float     aResolutionU = -1.0,
	    float     aResolutionV = -1.0,
	    float     aDetEffic    = 0.0,
	    float     aROtime      = -1,
	    float     aBkgRate     = 0.0);
  
  GGeoPetal(const GGeoPetal& other, TString Name = TString(""));

  // Destructor
  ~GGeoPetal();
  
  //Functions to access internal variables
  double    GetPetalBaseWidth()     const { return  WBase;  }
  double    GetPetalTopWidth()      const { return  WTop;   }
  double    GetPetalHeight()        const { return  Height; }
  TVector2  GetUInsensitive()       const { return  UInsensitive; }
  TVector2  GetVInsensitive()       const { return  VInsensitive; }
  
  //Runctions to set internal variables
  void      SetPetalBaseWidth(double aWBase)         { WBase  = aWBase;  }
  void      SetPetalTopWidth(double  aWTop)          { WTop   = aWTop;   }
  void      SetPetalHeight(double    aHeight)        { Height = aHeight; }
  void      SetUInsensitive(TVector2 aUInsensitive)  { UInsensitive = aUInsensitive; }
  void      SetVInsensitive(TVector2 aVInsensitive)  { VInsensitive = aVInsensitive; }
  
  //Set of functions to be defined for the daughter classes
  //void  GetThicknesses(double& w1, double& w2) {w1 = Thickness; w2 = 0.0; }
  GGeoObject* clone(TString aName) const;
  double  GetVolume();
  void    CheckInputs();
  void    FillSurfaces();
  bool    IsPointInsideGeometry(TVector3 PosXYZ);
  void    Print();

protected:

};

#endif //~ GGeoPetal_h

