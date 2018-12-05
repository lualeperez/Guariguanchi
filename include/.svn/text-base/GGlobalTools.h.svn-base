/***************************************************************************//**
 * @brief:      
 * @Description: Class with global variables and general functions
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

#ifndef GGlobalTools_h
#define GGlobalTools_h

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
#include <TFile.h>
#include <TF1.h>
#include <TRandom.h>
#include "TStopwatch.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

//Some global vairables and structures

const int     NMC_SeedEffic = 100;

const TString  WorldVolMaterial("Vacuum");

const int MAX_CHARS_PER_LINE  = 512;
const int MAX_TOKENS_PER_LINE = 200;
const char* const DELIMITER   = " ";
const char* const COMMENT     = "//";

const int      line_width     = 1;
const double   Dummy_value    = -999.0;  // a dummy value
const TVector3 Dummy_vector(-999.9,-999.9,-999.9); // dummy vector
const TString  Dummy_name("dummy");
const int      Nbins_mom_redu = 10;      // number of bins for momentum with plotting tracks in geometry visualizations

const int      MaxNLayers = 500;         // Maximum number of layers traversed by a particle
const int      MaxNhits   = 500;         // Maximum number of hits (traversing of sensitive layers)
const int      MaxNpars   =   5;         // Maximum number of parameters in a track parametrization

const TString  TheAxis[3] = {"X","Y","Z"}; //Regular laboratory axes

const double   epsilon_derivatives        = 1.0e-4; // Small number reference for numerical calculations
const int      Nmax_iterations            = 500;   // Maximum number of iterations for numerical derivatives calculation
const int      MaxIterations_intersection = 200;    // Maximum number of iterations for calculating trajectory <-> geometry-element intersection point

// Structure with the tracjectory parameter values at entrance and exit of a geometry-element
struct  GeoElementInOut_t {
  double s_in;
  double s_out;
};

// Structure with all needed informacion of an intersections point
struct  IntersectionHit_t {
  double s_in;              // trajectory parameter at the entrance of tracjectory at geometry-element
  double s_out;             // trajectory parameter at the exit     of tracjectory at geometry-element
  double s;                 // tracjectory parameter at intersection (average between s_in and s_out)
  int    geoElement_idx;    // geometry element index
  bool   IsSensitivePoint;  // is the intersection with a sensitive geometry-element?
};

// Structure with seed layers configurations
struct  SeedLayers_t {
  std::vector<TString> SeedLayers;
};

// Structure with seed layers probs
struct  LayersProbs_t {
  std::vector<TVector3> Probs; //0 -> good hit, 1-> fake hit, 2 -> null hit
};

// Structure with probabilities of track-hit association: correct, fake and null
struct  HitsConfiguration_t {
  long                  HitConfig;
  std::vector<TVector3> Probs; //0<-> good, 1<-> fake, 2 <-> null
  std::vector<long>     CorrectAndFakeHitConfigList;
};

// Structure with a voxel variables
struct  Voxel_t {
  double  Rx[2];      // X-interval
  double  Ry[2];      // Y-interval
  double  Rz[2];      // Z-interval
  double  Rtheta[2];  // theta-interval
  double  Rphi[2];    // phi-interval
  double  Rr[2];      // 3d-radial distance 
};

// Structure with momentum voxel
struct  MomVoxel_t {
  double Rp[2];     //momentum magnitude interval
  double Rtheta[2]; //momentum theta     interval
  double Rphi[2];   //momentum phi       interval
};

struct  TrackFindingRegion_t {
  Voxel_t     posRange;
  MomVoxel_t  momRange;
};

// Structure with a set of ranges in X, Y and Z
struct ARange_t {
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
};

// Structure with all possible cuts on a track
struct Cut_t {
  double   pMin;     // Minimum value of momentum
  double   pMax;     // Maximum value of momentum
  double   phiMin;   // Minimum value of phi   at origin
  double   phiMax;   // Maximum value of phi   at origin
  double   thetaMin; // Minimum value of theta at origin
  double   thetaMax; // Maximum value of theta at origin
  double   NhitsMin; // Minimum number of hits on a given system
  double   NhitsMax; // Maximum number of hits on a given system
};

// Structure with a list of track cuts 
struct TrackCuts_t {
  TString  SystemName;         // system name
  std::vector<Cut_t> CutList;  // List of cuts
};

//Structure with the parameters for the track resolution analysis
struct TrkResolAnalysisPars_t {
  int   NhitsMin;                //Minimum number of hits. Default value is tne minimum number of hits for tracking depending on the trajectory, i.e. 2 (3) for straigth (helix) tracks
  
  //Set of flags to control plots
  bool  SameRange;               //flag to produce track resolution parameters plots with the same vertical ranges. Default value is true
  bool  UseAllMomVals;           //flag to use all the momentum values specified in datacard. Default value is false
  bool  PlotMaterialBudget;      //flag to plot material budget encounter by track vs momentum, theta, phi. Default value is false
  bool  PlotDOCAatHighMom;       //flag to plot track's distance of closes approch to pivot point vs momentum, theta, phi. Default value is zero.
  bool  PlotDOCAvsMonFit;        //flag to plot the doca vs momentum fit
  bool  PlotPhiAveraged;         //flag to plot phi-averged track resolution performances. Default is false. Supersided by parameter below
  bool  PlotOnlyPhiAveraged;     //flag to plot only phi-avergaged track resolution performances. Default value is false
  bool  PlotPerformancesVsTheta; //flag to plot phi-avergaed performances vs theta for the various specified momenta. Default value is false
  bool  UseLogYAxes;             //flag to use log-y-axis
};

//Structure with the parameters for the material budget analysis
struct  MatBudgetAnalysisPars_t {
  float  mom_min;  // minimum momentum value for material budget analysis
  float  mom_max;  // maximum momentum value for material budget analysis
};

//System tructure
struct  System_t  {
  TString Name;
  std::vector<TString> LayersList;
};

//Structure with the Telescope resolution at DUT
struct  TelResolAtDUT_t {
  int     intersection_idx;
  int     geoElement_idx;
  double  s;
  double  TelResolU;
  double  TelResolV;
  double  TelCorrUV;
  double  Shadow_area;
  double  Nbkg;
};

//Structure with the parameters for the track efficiency analysis
struct EfficAnalysisPars_t {
  //Set of flags to control plots
  bool  SameRange;               //flag to produce track resolution parameters plots with the same vertical ranges. Default value is true
  bool  UseAllMomVals;           //flag to use all the momentum values specified in datacard. Default value is false
  bool  PlotPhiAveraged;         //flag to plot phi-averged track resolution performances. Default is false. Supersided by parameter below
  bool  PlotOnlyPhiAveraged;     //flag to plot only phi-avergaged track resolution performances. Default value is false
  bool  PlotPerformancesVsTheta; //flag to plot phi-avergaed performances vs theta for the various specified momenta. Default value is false
  bool  UseLogYAxes;             //flag to use log-y-axis
};

//Structure with the different efficiencies
struct Efficiencies_t {
  double Effic_tot;
  double Effic_NoFakes;
  double Effic_1Fake;
  double Effic_2orMoreFakes;
};

//Structure with a particle attributes
struct Particle_t {
  int      charge;       // charge
  int      Aweight;      // atomic weight (only in case of heavy ions)
  double   mass;         // mass
  TString  AntiParticle; // anti-particle
};

//Structure with a material attributes
struct Material_t {
  double  X0;                   // Radiation length
  double  density;              // density
  double  density_effect_1st;   // density effect first  junction point
  double  density_effect_2nd;   // density effect second junction point
  double  mI;                   // mean excitation energy
  double  mZA;                  // mean Z/A
  double  mZ;                   // mean Z
  int     color;                // color
};

class GGlobalTools {
  
private:
  
public:
  /// constructor(s)
  GGlobalTools();

  /// destructor
  ~GGlobalTools();
  
  void      OrderListList(std::vector<double>& List);
  TVector3  ProbsForTrackClusterAssociation(double det_effic,double rate, double corr, double sigmaU, double sigmaV,double chi2Ondf = 3.0, int ndf = 2);
  double    MSAngle(TString Aparticle, TVector3 momentum, double xOverX0);
  void      FillTables(void);
  bool      CheckParticleName(TString particle);
  double    GetX0FromMaterial(TString material);
  double    GetDensityFromMaterial(TString material);
  int       GetMaterialColor(TString material);
  double    GetDensityEffect1stFromMaterial(TString material);
  double    GetDensityEffect2ndFromMaterial(TString material);
  double    GetmIFromMaterial(TString material);
  double    GetmAZFromMaterial(TString material);
  
  double    GetParticleMass(TString particle);
  int       GetParticleCharge(TString particle);
  int       GetParticleAweight(TString particle);
  
  //Get all units
  double    GetUnit(const char* myunit);
  double    GetUnit(TString myunit);
  //Get distance units
  double    GetDistanceUnit(const char* myunit);
  double    GetDistanceUnit(TString myunit);
  bool      IsDistanceUnit(const char* myunit);
  bool      IsDistanceUnit(TString myunit);
  //Get time units
  double    GetTimeUnit(const char* myunit);
  double    GetTimeUnit(TString myunit);
  bool      IsTimeUnit(const char* myunit);
  bool      IsTimeUnit(TString myunit);
  //Get angle units
  double    GetAngleUnit(const char* myunit);
  double    GetAngleUnit(TString myunit);
  bool      IsAngleUnit(const char* myunit);
  bool      IsAngleUnit(TString myunit);
  //Get energy units
  double    GetEnergyUnit(const char* myunit);
  double    GetEnergyUnit(TString myunit);
  bool      IsEnergyUnit(const char* myunit);
  bool      IsEnergyUnit(TString myunit);
  //Get momentum units
  double    GetMomentumUnit(const char* myunit);
  double    GetMomentumUnit(TString myunit);
  bool      IsMomentumUnit(const char* myunit);
  bool      IsMomentumUnit(TString myunit);
  //Get mass units
  double    GetMassUnit(const char* myunit);
  double    GetMassUnit(TString myunit);
  bool      IsMassUnit(const char* myunit);
  bool      IsMassUnit(TString myunit);
  //Get B-field units
  double    GetBfieldUnit(const char* myunit);
  double    GetBfieldUnit(TString myunit);
  bool      IsBfieldUnit(const char* myunit);
  bool      IsBfieldUnit(TString myunit);
  //Get Rate density units
  double    GetRateDensityUnit(const char* myunit);
  double    GetRateDensityUnit(TString myunit);  
  bool      IsRateDensityUnit(const char* myunit);
  bool      IsRateDensityUnit(TString myunit);
  
  double    GetAngle(double x,double y);
  
  int       delta_cronequer(int i,int j);
  bool      CheckIfGoodSquareMatrix(TMatrixT<double> I);
  bool      CheckUnitMatrix(TMatrixT<double> I);
  bool      CheckIfSymmetricMatrix(TMatrixT<double> I);
  bool      CheckIfGoodErrorMatrix(TMatrixT<double> I);
  void      GetGeometryImpacParameterParameters(TGraphErrors* gr,double theta, TString sigma_units, TString momentum_units, bool FitPowerForImpactParam, double &a,double &b, double& exp);
  
  void      OrderIntersectionHitList(std::vector<IntersectionHit_t>& ItersectionHitList);
  
  void      GetListOfLayersTypes(long Config, std::vector<int>& LayerTypes);
  long      GetConfigFromListOfLayersTypes(std::vector<int> LayerTypes);
  void      IncludeFakesOnConfigList(long config,int Nfakes,bool IsSeeding,std::vector<long>& NewConfigList);
  void      GetNGoodFakesAndNull(long config,int& Ngood, int& Nfake, int& Nnull);
  void      FillConfigurations(int depth,
			       std::vector<int> Position,
			       int  Nlevels,
			       int* Config,
			       int* HitKind,
			       std::vector<long>& ConfigList);
  void      FillSeedConfigsWithFakes(std::vector<long>& SeedConfigList, int NfakesSeedMax);
  
  void      GetRotationMatrix(double angle,TString axis,TMatrixD& Rot);
  void      GetGlobalRotationMatrix(TVector3 angles,TMatrixD& Rot);
  void      GetGlobalRotationMatrix_FromList(std::vector<TString> axis,std::vector<double>  angles,TMatrixD& Rot);
  TVector3  GetRotationAngles(TMatrixD RotMatrix, bool Myverbose = false);
  void      RotateVector(TMatrixD Rot,TVector3& vector);
  void      GetRotationMatrixFromUnitVects(TVector3 UVector,TVector3 VVector,TVector3 WVector,TMatrixD& Rot);
  
  double    GetMomentumFromMomVar(double mom, TString Aparticle, TString momVar);
  
  //Mosaic geometry parameter calculation
  void      GetMosaicGeoParams(TString GeoConfig,
			       double& R,double w, double& rf, int n_ladders,
			       double& length, double& shift,  double& alpha,
			       TString TheVarPar,bool ShiftFix);
  void      GetMosaicGeoParams_Spiral(double& R,double w,double& rf, int n_ladders,
				      double& length, double& shift, double& alpha,
				      TString TheVarPar,bool ShiftFix);
  void      GetMosaicGeoParams_Alternating(double& R,double w,double& rf, int n_ladders,
					   double& length, double& shift, double& alpha,
					   TString TheVarPar,bool ShiftFix);
  void      CheckResolutionModelType(TString ModelType);
  
  TVector3  GetMomentum(double p,double theta,double phi);
  
  void      GetSurfAndMomOrthVects(TVector3 momentum, TVector3& m1,TVector3& m2);
  
  bool      IsPointInVoxel(Voxel_t aVoxel,TVector3 Point);
  
  bool      IsInitCondsInTrackFindingRegion(TrackFindingRegion_t aRegion, TVector3 x0, TVector3 p0);
  
  TString   GetOutputDirectory(TString TheOutputFile);
  
  bool      SetBoolFromString(TString string);
  
  TString   GetStringFromBool(bool fff);
  
  double    FromCosThetaToTheta(double costheta);
  double    FromThetaToCosTheta(double theta);
  double    FromEtaToTheta(double eta);
  double    FromThetaToEta(double theta);
  
  double    GetEllipseArea(double sigmaU, double sigmaV, double corr);
  
  TString   GetAntiParticle(TString particle);
  
  void      RemoveElementFromList(int idx, std::vector<int> &List);
  
  double    BetheBlochGeant(double  momentum,
			    TString particle,
			    TString material);
  
  double    BetheBlochGeant(double bg,
			    int    charge,
			    double density,
			    double density_effect_1st,
			    double density_effect_2nd,
			    double mI,
			    double mZA);
  
  double    BetheBlochPDG(double  momentum,
			  TString particle,
			  TString material);
  
  void     GetBetheBlochParams(TString material,
			       double& density,
			       double& density_effect_1st,
			       double& density_effect_2nd,
			       double& mI,
			       double& mZA,
			       double& mZ);
  
  double   GetKappa(void) { return kappa; }
  void     SetKappa(double akappa);
  
  double   GetSigmaELoss(double  momentum,
			 TString particle,
			 TString material,
			 double dE);
  
  double   GetSigmaELossPDG(double  momentum,
			    TString particle,
			    TString material,
			    double  X);

  std::map<TString,Material_t>  MaterialMap;        // Map with material attributes
  std::map<TString,Particle_t>  ParticleMap;        // Map with particle attributes
  
  std::map<TString,double>  units;                  // Map with set of units used by code
  std::map<TString,double>  units_distance;         // Map with set of distance units
  std::map<TString,double>  units_time;             // Map with set of time     units
  std::map<TString,double>  units_angle;            // Map with set of angle    units
  std::map<TString,double>  units_Energy;           // Map with set of Energy   units
  std::map<TString,double>  units_Momentum;         // Map with set of Momentum units
  std::map<TString,double>  units_Mass;             // Map with set of Mass     units
  std::map<TString,double>  units_Bfield;           // Map with set of Bfield   units
  std::map<TString,double>  units_RateDensity;      // Map with set of Rate density units
  
  std::vector<TString>      GeoLadderTypes;         // List of implemented ladder types
  std::vector<TString>      GeoPlaneTypes;          // List of implemented plane  types
  std::vector<TString>      WorldVolumeTypes;       // List of implemented world volumes types
  std::vector<TString>      ResolutionModelTypes;   // List of implemented resolution models
  std::vector<TString>      PatternRecogAlgoTypes;  // List of implemented pattern recognition algorithms for pseudo-efficiency calculation
  
  double betaCurvRadius; // constan defining the ralationship between transverse momentum (Pt), B-field (B) and curvature radius (R): betaCurvRadius = Pt / (B * R)
    
  double InfiniteMon;    // high value of momentum

  double epsilon_Bfield;   // small reference value for a B-field
  double epsilon_distance; // small reference value for a distance

  bool GlobalDoRichardson; // bool to decide if using richarson method's for derivative numerical calculation
  
  TRandom* rand;
  
  TStopwatch fWatch;
  
protected:

  bool verbose; // bool for verbosity
  
  double K;     // Constant for Bethe-Bloch calculations
  
  double kappa; // Constant for the fluctuations on the E-loss
  
};

#endif //~ GGlobalTools_h

