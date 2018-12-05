/***************************************************************************//**
 * @brief:      
 * @Description: Set of utilitary functions
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
 

#ifndef GUtiliratyFunctions_h
#define GUtiliratyFunctions_h

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
#include "include/GGeoObject.h"
#include "include/GGlobalTools.h"
#include "include/GTrajectory.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

double    IntersectCoordinates(GSurfaceObject* aSurface,
			       GTrajectory* Trajectory,
			       GGlobalTools* global,
			       double sinit = Dummy_value,
			       bool DoPositiveS = true,
			       double scut = 0,
			       bool verbose = false);


bool      DoGeometryElementsOverlap(GGeoObject*  aGeoObject1, GGeoObject* aGeoObject2, GGlobalTools* global, double Precision);

bool      DoSurfaceIntersectsObject(GGeoObject*  aGeoObject,  GSurfaceObject* aSurfaceObject, GGlobalTools* global, double Precision);


#endif //~ GUtiliratyFunctions_h

