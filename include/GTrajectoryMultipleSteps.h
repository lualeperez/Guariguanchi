/***************************************************************************//**
 * @brief:      
 * @Description: Class for multiple step trajectories, i.e. a trajectory within an Multiple Steps Bfield
 *               (see GBFieldMultipleSteps class)
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
 

#ifndef GTrajectoryMultipleSteps_h
#define GTrajectoryMultipleSteps_h

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
#include "include/GBFieldConstant.h"
#include "include/GBFieldMultipleSteps.h"
#include "include/GTrajectory.h"
#include "include/GTrajectoryStraight.h"
#include "include/GTrajectoryHelix.h"
#include "include/GGeoObject.h"
#include "include/GUtiliratyFunctions.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

struct  TrackOrder_t {
  bool IsInsideVolTrack;
  int  idx_InVolTrack;
  int  idx_OutVolTrack;
  double s_referece;
};

class GTrajectoryMultipleSteps : public GTrajectory {
  
private:
  
  double small_distance;
  double large_distance;
  int    Ninside_plus;
  int    Noutside_plus;
  
  double MaxOmega;
  
public:

  GBFieldMultipleSteps* TheBField;
  
  std::vector<GGeoObject*>        VolumeList;
  
  std::vector<TVector3>           SInterVals_true;
  std::vector<TrackOrder_t>       TrueTrackOrderList;
  std::vector<GTrajectory*>       TrueInsideVolTrajectoryList;
  std::vector<GTrajectory*>       TrueOutsideVolTrajectoryList;
  
  double z_ref;
  double s_ref;
  
  std::vector<TVector3>           SInterVals_fit;
  std::vector<TrackOrder_t>       FitTrackOrderList;
  std::vector<GTrajectory*>       FitInsideVolTrajectoryList;
  std::vector<GTrajectory*>       FitOutsideVolTrajectoryList;
  
  // Constructor
  GTrajectoryMultipleSteps(TString   aName,
			   TString   aParticle,
			   TVector3  aPos0,
			   TVector3  aMom0,
			   TVector3  aRefPoint,
			   GBField*  aBfield,
			   GGlobalTools* aglobal);
  
  GTrajectoryMultipleSteps(const GTrajectoryMultipleSteps& other, TString Name = TString(""));

  // Destructor
  virtual ~GTrajectoryMultipleSteps();
  
  //Set of generic functions
  
  //Functions to access internal variables
  TVector3  GetInBField(int idx)  const { return  (dynamic_cast<GBFieldMultipleSteps*>(Bfield))->GetInBField(idx);  }
  TVector3  GetOutBField()        const { return  (dynamic_cast<GBFieldMultipleSteps*>(Bfield))->GetOutBField(); }
  
  //Functions to set internal variables
  void  SetBfield(GBField* aBfield);
  void  FillTrajectories(void);
  
  void  SetTrueTrajectories(TVector3 x0, TVector3 p0);
  
  void  SetFitTrajectories(TVector3 x0, TVector3 p0);
  
  void  SetTrueTrajectories(void);
  
  void  SetFitTrajectories(void);
  
  double  GetSFromZ(double zval, double sinit = 0, bool FromTrueTrajectory = true);
  
  //Set of functions to be defined for the daughter classes
  void      SetAllParameters(std::vector<double> aFitParams);
  
  GTrajectory*  clone(TString aName) const;

  int        GetTheDOCAParIndex()                { return   0;}
  double     GetParamError(int idx, TMatrixD  FitCovMatrix);
  double     GetPtAtDOCA();
  double     GetOmega();
  double     GetSLimit(int nloops);
  int        GetMaximumIntersections()           { return 3; }
  double     GetDummyParamEpsilon()              { return global->GetAngleUnit("urad"); }
  TString    GetDummyUnits()                     { return TString("deg"); }
  double     DummyValueCorrection(double aDummy);
  double     GetInitValueParForDerivativeCalculation(int idx);
  int        GetMinHits() { return 3; }
  double     GetInitialSFor2ndIntersection(double s);
  
  //Equations of motion functions
  void       CalculatePrimedReframe();
  TVector3   GetTrueTrajectoryCoordinates(double s);
  TVector3   GetTrueTrajectoryUnitMon(double s);
  
  TVector3   GetFitTrajectoryCoordinates(double s);
  TVector3   GetFitTrajectoryMon(double s);
  
  //Parameterized track functions
  void       FillParametersNames(TString PtParamFormat = TString(""));
  void       GetFitParsFromInitConds();
  TVector3   GetFitTrackCoordinates(double dummy,double sinit = 0);
  TVector3   GetFitTrackMomentum(double dummy, double sinit = 0);
  double     GetFitTrackDummyParFromS(double s);
  TVector3   GetFitTrackCoorDerWRTDummy(double dummy, double sinit = 0);
  void       PrintParameters();

protected:

};

#endif //~ GTrajectoryMultipleSteps_h

