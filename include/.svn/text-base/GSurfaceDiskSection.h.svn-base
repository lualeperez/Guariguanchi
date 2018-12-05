/***************************************************************************//**
 * @brief:      
 * @Description: Class for disk section surfaces
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

#ifndef GSurfaceDiskSection_h
#define GSurfaceDiskSection_h

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

class GSurfaceDiskSection : public GSurfaceObject {

private:
  
public:
  
  double    Rin;                 // Disk inner radius
  double    Rout;                // Disk outer radius
  double    DeltaPhi;            // Disk opening angle
  TVector2  UInsensitive;        // Fraction of insensitive width w.r.t. to total U width
  TVector2  VInsensitive;        // Fraction of insensitive width w.r.t. to total V width
  
  // Constructor
  GSurfaceDiskSection(TString   aName,
		      TVector3  aPosition,
		      TMatrixD  aRot,
		      double    aRin,
		      double    aDeltaPhi,
		      double    aRout,
		      TVector2  aUInsensitive,
		      TVector2  aVInsensitive,
		      bool      aIsTrackingLayer,
		      GGlobalTools* aglobal);
  
  GSurfaceDiskSection(const GSurfaceDiskSection& other, TString Name = TString(""));

  // Destructor
  ~GSurfaceDiskSection();
  
  //Set of generic functions
  double    GetDiskSectionInnerRadius() const { return  Rin;  }
  double    GetDiskSectionOuterRadius() const { return  Rout; }
  double    GetDiskSectionDeltaPhi()    const { return  DeltaPhi; }
  TVector2  GetUInsensitive()           const { return  UInsensitive; }
  TVector2  GetVInsensitive()           const { return  VInsensitive; }
  
  //Functions to access internal variables
  void      SetDiskSectionInnerRadius(double aRin)   { Rin      = aRin; }
  void      SetDiskSectionOuterRadius(double aRout)  { Rout     = aRout; }
  void      SetDiskSectionDeltaPhi(double aDeltaPhi) { DeltaPhi = aDeltaPhi; }
  void      SetUInsensitive(TVector2 uSensitive)     { UInsensitive = uSensitive; }
  void      SetVInsensitive(TVector2 vSensitive)     { VInsensitive = vSensitive; }
  
  //Set of functions to be defined for the daughter classes
  void      CheckInputs();
  TVector3  GetInitPositionForIntersection();
  
  TVector3  GetUVWFromXYZ(TVector3 PositionXYZ);
  TVector3  GetDerUVWFromXYZ_Analytical(TVector3 PositionXYZ, int idx);
  TVector3  GetXYZFromUVW(TVector3 PositionUVW);
  
  GSurfaceObject* clone(TString aName) const;
  bool      IsInMaterial(TVector3 PositionUVW);
  bool      IsInSensitiveMaterial(TVector3 PositionUVW);
  TVector3  GetNormVector(TVector3 PositionXYZ);
  void      VarySurface(double w1, double w2) { Position += w1*WVector; }
  void      GetVRange(double& Vmin, double& Vmax) {Vmin = -0.5*(Rout - Rin); Vmax = +0.5*(Rout - Rin); }
  void      GetURange(double V, double& Umin, double& Umax) {Umin = -0.5*DeltaPhi*(V + 0.5*(Rout + Rin)); Umax = +0.5*DeltaPhi*(V + 0.5*(Rout + Rin)); }
  void      FillSurfaceRepresentation(int color,
				      std::vector<TGraph>   &GraphXY,
				      std::vector<TGraph>   &GraphZY,
				      std::vector<TGraph>   &GraphZX,
				      std::vector<TGraph>   &GraphZR,
				      int linestyle = 1);

protected:
  
};

#endif //~ GSurfaceDiskSection_h

