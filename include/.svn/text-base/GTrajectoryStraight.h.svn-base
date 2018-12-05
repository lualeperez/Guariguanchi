/***************************************************************************//**
 * @brief:      
 * @Description: Class for straight trajectories
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
 

#ifndef GTrajectoryStraight_h
#define GTrajectoryStraight_h

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
#include "include/GGlobalTools.h"
#include "include/GBField.h"
#include "include/GTrajectory.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GTrajectoryStraight : public GTrajectory {
  
private:
  
public:
  
  // Constructor
  GTrajectoryStraight(TString   aName,
		      TString   aParticle,
		      TVector3  aPos0,
		      TVector3  aMom0,
		      TVector3  aRefPoint,
		      GBField*  aBfield,
		      GGlobalTools* aglobal);
  
  GTrajectoryStraight(const GTrajectoryStraight& other, TString Name = TString(""));

  // Destructor
  virtual ~GTrajectoryStraight();
  
  //Set of generic functions
  
  //Functions to access internal variables
  
  //Functions to set internal variables
  
  //Set of functions to be defined for the daughter classes
  GTrajectory*  clone(TString aName) const;
  
  int        GetTheDOCAParIndex()                { return   1;}
  double     GetPtAtDOCA()                       { return  -1.0; }
  double     GetOmega()                          { return 0.0;}
  double     GetSLimit(int nloops)               { return 1.0e+3*global->GetDistanceUnit("m"); }
  int        GetMaximumIntersections()           { return 2; }
  double     GetDummyParamEpsilon()              { return global->GetDistanceUnit("nm"); }
  TString    GetDummyUnits()                     { return TString("mm"); }
  double     DummyValueCorrection(double aDummy) { return aDummy; }
  double     GetInitValueParForDerivativeCalculation(int idx);
  int        GetMinHits() { return 2; }
  double     GetInitialSFor2ndIntersection(double s) { return s*10; }
  
  //Equations of motion functions
  void       CalculatePrimedReframe();
  TVector3   GetTrueTrajectoryCoordinates(double s);
  TVector3   GetTrueTrajectoryUnitMon(double s);
  
  //Parameterized track functions  
  void       FillParametersNames(TString PtParamFormat = TString(""));
  void       GetFitParsFromInitConds();
  TVector3   GetFitTrackCoordinates(double dummy, double sinit = 0);
  TVector3   GetFitTrackMomentum(double dummy, double sinit = 0);
  double     GetFitTrackDummyParFromS(double s);
  TVector3   GetFitTrackCoorDerWRTDummy(double dummy, double sinit = 0);
  void       PrintParameters();

protected:
  
};

#endif //~ GTrajectoryStraight_h

