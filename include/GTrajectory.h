/***************************************************************************//**
 * @brief:      
 * @Description: Base class for trajectories
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

#ifndef GTrajectory_h
#define GTrajectory_h

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

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GTrajectory {
  
private:
  
public:
  
  TString   Name;                            // Trajectory name
  TString   Type;                            // Trajectory type: Straight, Helix, ...
  
  TString   Particle;                        // particle name
  int       charge;                          // particle charge
  double    mass;                            // particle mass
  TVector3  pos0;                            // particle initial position
  TVector3  mom0;                            // particle initial momentum
  
  TVector3  RefPoint;                        // Reference point for tracking
  
  std::vector<double>  FitParams;            // List of fit parameters
  
  std::vector<TString> FitParamNames;        // List of fit parameters names
  std::vector<TString> FitParamUnits;        // List of fit parameters units
  std::vector<TString> FitParamUnitsTitles;  // List of fit parameters units
  
  std::vector<TString> FitParamErrorNames;        // List of fit parameters titles
  std::vector<TString> FitParamErrorUnits;        // List of fit parameters units
  std::vector<TString> FitParamErrorUnitsTitles;  // List of fit parameters units
  
  GBField*  Bfield;                   // Bfield
  
  //Internal coordinate system
  TVector3 xprimeVect;
  TVector3 yprimeVect;
  TVector3 zprimeVect;
  
  GGlobalTools* global;
  
  // Constructor
  GTrajectory(TString   aName,
	      TString   aParticle,
	      TVector3  aPos0,
	      TVector3  aMom0,
	      TVector3  aRefPoint,
	      GBField*  aBfield,
	      GGlobalTools* aglobal);
  
  GTrajectory(const GTrajectory& other, TString Name = TString(""));

  // Destructor
  virtual ~GTrajectory();
  
  //Set of generic functions
  
  //Functions to access internal variables
  TString   GetName()             const { return  Name; }
  TString   GetType()             const { return  Type; }
  TString   GetParticle()         const { return  Particle; }
  TVector3  GetInitPosition()     const { return  pos0; }
  TVector3  GetInitMomentum()     const { return  mom0; }
  GBField*  GetBfield()           const { return  Bfield; }
  int       GetNParameters()      const { return  FitParams.size(); } 
  double    GetParameter(int idx) const;
  void      GetAllParameters(std::vector<double>& aFitParams);
  bool      GetVerbose() { return verbose; }
  
  void      GetInternalRF(TVector3& aXvect, TVector3& aYvect, TVector3& aZvect)  { aXvect = xprimeVect; aYvect = yprimeVect; aZvect = zprimeVect; };
  
  TString   GetParameterName(int idx) const;
  TString   GetParameterUnit(int idx) const;
  TString   GetParameterUnitTitle(int idx) const;
  TString   GetParameterErrorName(int idx) const;
  TString   GetParameterErrorUnit(int idx) const;
  TString   GetParameterErrorUnitTitle(int idx) const;
  
  //Functions to set internal variables
  void      SetName(TString aName)                  { Name = aName; }
  void      SetType(TString aType)                  { Type = aType; }
  void      SetParticle(TString aParticle);
  void      SetInitPosition(TVector3 aInitPosition);
  void      SetInitMomentum(TVector3 aInitMomentum);
  void      SetInitPositionAndMomentum(TVector3 aInitPosition, TVector3 aInitMomentum);
  virtual   void      SetAllParameters(std::vector<double> aFitParams);
  
  void      SetVerbose(bool aVerbose) { verbose = aVerbose; }
  
  double     GetParamsCorrelation(int idx1, int idx2, TMatrixD  FitCovMatrix);
  
  virtual   void  SetBfield(GBField* aBfield);
  
  TVector3  GetTrueTrajectoryMon(double s);
  
  //Set of functions to be defined for the daughter classes
  virtual  GTrajectory*  clone(TString aName) const;
  
  virtual  int        GetTheDOCAParIndex() { return -1;}
  virtual  double     GetParamError(int idx, TMatrixD  FitCovMatrix);
  virtual  double     GetPtAtDOCA()                       { return  -1.0; }
  virtual  double     GetOmega()                          { return   0.0; }
  virtual  double     GetSLimit(int nloops)               { return   0.0; }
  virtual  int        GetMaximumIntersections()           { return   1; }
  virtual  double     GetDummyParamEpsilon()              { return   global->GetDistanceUnit("mm"); }
  virtual  TString    GetDummyUnits()                     { return   TString("mm"); }
  virtual  double     DummyValueCorrection(double aDummy) { return   aDummy; }
  virtual  double     GetInitValueParForDerivativeCalculation(int idx) { return 1.0e-8; }
  virtual  int        GetMinHits() { return 0; }
  virtual  double     GetInitialSFor2ndIntersection(double s) { return s*10; }
  
  //Equations of motion functions
  virtual  void       CalculatePrimedReframe() {;}
  virtual  TVector3   GetTrueTrajectoryCoordinates(double s) { return TVector3(0.0,0.0,0.0); }
  virtual  TVector3   GetTrueTrajectoryUnitMon(double s) { return TVector3(0.0,0.0,0.0); }
  
  //Parameterized track functions
  virtual  void       FillParametersNames(TString PtParamFormat = TString("")) {;}
  virtual  void       GetFitParsFromInitConds() {;}
  virtual  TVector3   GetFitTrackCoordinates(double dummy, double sinit = 0) { return TVector3(0.0,0.0,0.0); }
  virtual  TVector3   GetFitTrackMomentum(double dummy, double sinit = 0) { return TVector3(0.0,0.0,0.0); }
  virtual  double     GetFitTrackDummyParFromS(double s) { return 0.0; }
  virtual  TVector3   GetFitTrackCoorDerWRTDummy(double dummy, double sinit = 0) { return TVector3(0.0,0.0,0.0); }
  virtual  void       PrintParameters() {;}
  

protected:

  bool verbose;
  
};

#endif //~ GTrajectory_h

