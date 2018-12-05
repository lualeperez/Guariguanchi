/***************************************************************************//**
 * @brief:      
 * @Description: Class for tracking related calculations
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
 

#ifndef GTracker_h
#define GTracker_h

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
#include "include/GGeoObject.h"
#include "include/GGeoPlane.h"
#include "include/GGeoCylinder.h"
#include "include/GGeoDisk.h"
#include "include/GGeoPetal.h"
#include "include/GGeoCone.h"
#include "include/GSurfaceObject.h"
#include "include/GSurfacePlane.h"
#include "include/GSurfaceCylinder.h"
#include "include/GSurfaceDisk.h"
#include "include/GSurfaceCone.h"
#include "include/GSurfacePetal.h"
#include "include/GGlobalTools.h"
#include "include/GBField.h"
#include "include/GBFieldConstant.h"
#include "include/GResolutionModel.h"
#include "include/GGeometry.h"
#include "include/GTrajectory.h"
#include "include/GTrajectoryStraight.h"
#include "include/GTrajectoryHelix.h"
#include "include/GTrajectoryMultipleSteps.h"

#include "include/GTrackFinderAlgo.h"
#include "include/GTrackFinderAlgoFPCCD.h"

#include "include/GUtiliratyFunctions.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GTracker {
  
private:
  
  //Set of temporary variables needed for calculation
  TMatrixD** F;
  TMatrixD** K;
  double**** Cov_deltaR_deltaR;
  double**   CovUU;
  double**   CovVV;
  double**   CovUV;
  
  TMatrixD* Fr_local;
  TMatrixD* Fp_local;
  TMatrixD* Kr_local;
  TMatrixD* Kp_local;
  double*   SigmaMS2_local;
  TVector3*  m1Unit_local;
  TVector3*  m2Unit_local;
  double*    SigmaEloss2_local;
  TVector3*  pUnit_local;
  double*** DerUVwrtPos_local;
  double**  GammaUDer_local;
  double**  GammaVDer_local;
  
  TMatrixD* Fr_global;
  TMatrixD* Fp_global;
  TMatrixD* Kr_global;
  TMatrixD* Kp_global;
  double*   SigmaMS2_global;
  TVector3*  m1Unit_global;
  TVector3*  m2Unit_global;
  double*    SigmaEloss2_global;
  TVector3*  pUnit_global;
  double*** DerUVwrtPos_global;
  double**  GammaUDer_global;
  double**  GammaVDer_global;
  double*** Der_Params_prime_wrt_thetaMS_global;
  double**  Der_Params_prime_wrt_Eloss_global;
  
  double   MinMSEffect; //Minimum multimple scatteting angle to include in cov matrix calculation

public:
  
  TString Name;  // Tracker Name
  
  GGlobalTools* global;
  
  GTrajectory* Trajectory;
  
  GGeometry*   Geometry;
  
  // Constructor
  GTracker(TString       aName,
	   GGeometry*    aGeometry,
	   TString       aParticle,
	   TVector3      aPos0,
	   TVector3      aMom0,
	   TVector3      aRefPoint,
	   GGlobalTools* aglobal);

  // Destructor
  virtual ~GTracker();
  
  void          InitializeTrajectory(TString   aParticle,
				     TVector3  aPos0,
				     TVector3  aMom0,
				     TVector3  aRefPoint);
  
  //Functions to access internal variables
  TString       GetName()       const { return  Name; }
  GTrajectory*  GetTrajectory() const { return  Trajectory; }
  GGeometry*    GetGeometry()   const { return  Geometry; }
  bool          GetVerbose() { return verbose; }
  
  //Functions to set internal variables
  void      SetName(TString aName)  { Name = aName; }
  void      SetVerbose(bool aVerbose) { verbose = aVerbose; }
  void      SetTrajectoryInitConditions(TVector3 pos0, TVector3 mom0) { Trajectory->SetInitPositionAndMomentum(pos0,mom0); }
  
  void      SetIncludeEloss(bool aIncludeEloss) { IncludeEloss = aIncludeEloss; }
  
  //Tracker Functions
  //Equations of motion functions
  double    GetIntersectCoordinates(GSurfaceObject* aSurface, double sinit = Dummy_value, double scut = 0.0, bool Myverbose = false);
  void      GetAllIntersectCoordinates(GSurfaceObject* aSurface, std::vector<double>& SListIntersections);
  double    GetIntersectionWithWorldVolume();
  void      GetIntersectionsWithGeometry(std::vector<IntersectionHit_t>& ItersectionHitList);
  double    GetIntersectionsWithGeoElement(double s0, GGeoObject* aGeoElement);
  double    GetIntersectionXOX0(IntersectionHit_t AHit, GGeoObject* aGeoElement);
  void      GetMaterialBudget(double& MatBudget, int& Nhits, double& dist1stPointToRef);
  void      GetMaterialBudget_FromPlaneList(std::vector<IntersectionHit_t> ItersectionHitList,
					    double& MatBudget, int& Nhits, double& dist1stPointToRef);
  void      GetGeometrySystemsMaterialBudget(std::vector<double>& SystemMatBudget, double& TotMatBudget);
  
  TVector3  GetTrueTrackInterCoordDerWRTPos_Numerical(double s0, GGeoObject* aGeoElement, int der_index);
  TVector3  GetTrueTrackInterCoordDerWRTMom_Numerical(double s0, GGeoObject* aGeoElement, int der_index);
  TVector3  GetTrueTrackInterMomDerWRTPos_Numerical(  double s0, GGeoObject* aGeoElement, int der_index);
  TVector3  GetTrueTrackInterMomDerWRTMom_Numerical(  double s0, GGeoObject* aGeoElement, int der_index);
  
  void      GetTrueTrackInterCoordAndMomDerWRTPos_Numerical(double s0, GGeoObject* aGeoElement, int der_index, TVector3& DerCoor, TVector3& DerMom);
  void      GetTrueTrackInterCoordAndMomDerWRTMom_Numerical(double s0, GGeoObject* aGeoElement, int der_index, TVector3& DerCoor, TVector3& DerMom);
  
  //Parameterized track functions
  double    GetFitTrackIntersectionCoordinates(GSurfaceObject* aSurface,double dummy_init, double sinit);
  TVector3  GetFitTrackIntersectionCoordinatesUVW(GSurfaceObject* aSurface,double dummy_init, double sinit);
  TVector3  GetFitTrackIntersectionCoordUVWDerWRTPar_Numerical(GSurfaceObject* aSurface, int  parameter, double dummy_init, double sinit);
  
  void      ConfigureTelescopePlanes(std::vector<IntersectionHit_t>&  ItersectionHitList);
  void      GetTelescopeResolutionAtDUTPlanes(std::vector<IntersectionHit_t> ItersectionHitList, 
					      TMatrixD  FitCovMatrix,
					      std::vector<TelResolAtDUT_t>&  TelResolAtDUTList);
  
  void      GetAllTrackParamDerWrtMS(std::vector<IntersectionHit_t>  ItersectionHitList,
				     double*                         SigmaMS2,
				     TVector3*                       m1Unit,
				     TVector3*                       m2Unit,
				     double***                       Der_Params_prime_wrt_thetaMS);
  
  void      GetAllTrackParamDerWrtEloss(std::vector<IntersectionHit_t>  ItersectionHitList,
					double*                         SigmaEloss2,
					TVector3*                       pUnit,
					double**                        Der_Params_prime_wrt_Eloss);
  
  void      GetNewParCovarianceMatrix(std::vector<IntersectionHit_t>  ItersectionHitList_global,
				      std::vector<IntersectionHit_t>  ItersectionHitList_current,
				      int                             nlayer_current,
				      double***                       Der_Params_prime_wrt_thetaMS,
				      double*                         SigmaMS2,
				      double**                        Der_Params_prime_wrt_Eloss,
				      double*                         SigmaEloss2,
				      TMatrixD                        FitCovMatrix_old,
				      TMatrixD&                       FitCovMatrix_new);
  
  void      GetCovMatrixUVOfTrackIntersection(std::vector<IntersectionHit_t>  ItersectionHitList_global,
					      std::vector<IntersectionHit_t>  ItersectionHitList_current,
					      int                             nlayer_current,
					      double**                        GammaUDer,
					      double**                        GammaVDer,
					      TMatrixD                        FitCovMatrix,
					      TMatrixD&                       CovUV);
  
  void      GetUVDerivatiesWrtPars(std::vector<IntersectionHit_t>  ItersectionHitList,
				   double**                        GammaUDer,
				   double**                        GammaVDer,
				   const int                       Npars);
  
  void      GetFsAndKsMS(std::vector<IntersectionHit_t>  ItersectionHitList,
			 TMatrixD*                       Fr,
			 TMatrixD*                       Fp,
			 TMatrixD*                       Kr,
			 TMatrixD*                       Kp);
  
  void      GetSigmaMSAndPerpVectors(std::vector<IntersectionHit_t> ItersectionHitList,
				     double*                        SigmaMS2,
				     TVector3*                      m1Unit,
				     TVector3*                      m2Unit);
  
  void      GetSigmaElossAndUnitVector(std::vector<IntersectionHit_t> ItersectionHitList,
				       double*                        SigmaEloss2,
				       TVector3*                      pUnit);
  
  void      GetUVDerWrtPosForAllLayers(std::vector<IntersectionHit_t>  ItersectionHitList,
				       double***                       DerUVwrtPos);
  
  void      GetHitUVCovMatrixWithMSFromGlobalCalc(std::vector<IntersectionHit_t> ItersectionHitList_global,
						  std::vector<IntersectionHit_t> ItersectionHitList_current,
						  TMatrixD*  Fr, TMatrixD* Fp, TMatrixD* Kr, TMatrixD* Kp,
						  double*    SigmaMS2, TVector3*  m1Unit, TVector3*  m2Unit,
						  double*    SigmaEloss2, TVector3*  pUnit,
						  double***  DerUVwrtPos,
						  TMatrixD&  HitUVCovMatrix);
  
  void      GetHitUVCovMatrixWithMS(std::vector<IntersectionHit_t> ItersectionHitList,
				    TMatrixD&            HitUVCovMatrix);
  
  bool      GetFitTrackParsCovMatrix(int                  NhitsMin,
				     TMatrixD&            FitCovMatrix);
  
  bool      doTrkResolAnalysis(int        NhitsMin,
			       bool       DoTelescopeAnalysis,
			       int&       Nhits, double& Material_budget, double &dist1stPointToRef,
			       TMatrixD&  FitCovMatrix,
			       std::vector<TelResolAtDUT_t>& TelResolAtDUTList);
  
  bool      GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc(std::vector<IntersectionHit_t>  ItersectionHitList_global,
								  std::vector<IntersectionHit_t>  ItersectionHitList_current,
								  double**                        GammaUDer,
								  double**                        GammaVDer,
								  TMatrixD                        MeasCovMatrix,
								  TMatrixD&                       FitCovMatrix);
  
  bool      GetFitTrackParsCovMatrix_FromPlaneList(std::vector<IntersectionHit_t> ItersectionHitList,
						   int                  NhitsMin,
						   TMatrixD&            FitCovMatrix);
  
  void      GetFitTrackPseudoEfficiency(Efficiencies_t&      Efficiencies,
					TMatrixD&            AveFitCovMatrix);
  
  void      GetFitTrackPseudoEfficiency_FPCCD(GTrackFinderAlgoFPCCD* aTrackFinderAlgo,
					      Efficiencies_t&      Efficiencies,
					      TMatrixD&            AveFitCovMatrix);
  
  
  
protected:

  bool IncludeEloss;
  
  bool verbose;
  
};

#endif //~ GTracker_h

