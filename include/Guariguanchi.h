/***************************************************************************//**
 * @brief:      
 * @Description: Analytical tool for studiying tracking performances
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

#ifndef Guariguanchi_h
#define Guariguanchi_h

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
#include <TImage.h>
#include <TASImage.h>
#include <TAttImage.h>
#include <TPad.h>
#include "include/GBField.h"
#include "include/GBFieldConstant.h"
#include "include/GBFieldMultipleSteps.h"
#include "include/GGeometry.h"
#include "include/GGeoObject.h"
#include "include/GSurfaceObject.h"
#include "include/GTrajectory.h"
#include "include/GTracker.h"
#include "include/GGlobalTools.h"
#include "include/GResolutionModel.h"
#include "include/GResolutionModelTPC.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class Guariguanchi {
  
private:
  
  //Logo image
  TASImage* logo;
  
  std::vector<Voxel_t>    VoxelList;     //List of voxels
  std::vector<ARange_t>   GeoRanges;     //List of geometry ranges for plotting geometry
  std::vector<GGeometry*> GeometryList;  //List of geometies
  std::vector<GTracker*>  TrackerList;   //List of Trackers 
  
  std::vector<TrkResolAnalysisPars_t>  TrkResolAnalysisPars;  //Track Resolution Analysis Parameters
  std::vector<MatBudgetAnalysisPars_t> MatBudgetAnalysisPars; //Material budget analysis parameters 
  std::vector<EfficAnalysisPars_t>     EfficAnalysisPars;  //Track Resolution Analysis Parameters
  
  //TStopwatch fWatch;
  
  //Setup:
  void    ReadDataCard(void);
  void    FillGeometryFromDataCard(const char* datacard,int index);
  void    ReadGeometryDataCard(const char* datacard, GGeometry* aGeometry);
  void    ReadBField(TString FieldName,TString  BlockName,ifstream* fin, GBField* &aBField, const char* datacard, int &line_number);
  void    ReadWorldVolume(TString WorldVolName,TString  BlockName,ifstream* fin,GGeoObject* &aWorldolume, int geoID, const char* datacard, int &line_number);
  void    ReadResolutionModel(TString  BlockName,ifstream* fin,GResolutionModel* &aResolModel, const char* datacard, int &line_number);
  void    ReadGasDetector(ifstream* fin,GGeometry* aGeometry, const char* datacard, int &line_number);
  void    ReadLadderBlock(TString  BlockName,ifstream* fin,GGeometry* aGeometry, const char* datacard, int &line_number);
  void    ReadSystemsConfiguration(ifstream* fin,GGeometry* aGeometry,const char* datacard, int &line_number);
  void    ReadTrackFinderAlgo(TString  BlockName,ifstream* fin,GGeometry* aGeometry,const char* datacard, int &line_number);
  
  //GeometryHandler:
  void    FillGeometries(void);
  void    FillTrackers(void);
  void    PrintGeometries(void);
  void    PrintGeometryWeights(void);
  void    PlotGeometries(void);
  void    CheckGeometries(void);
  void    CheckGeometriesComparability(void);
  void    CheckTelescopeConfigurations(void);
  
  void    BuildSingleFile(void);
  
  void    PrintTrkResolAnalysisParamsBlock();
  void    PrintMatBudgetAnalysisParamsBlock();
  void    PrintEfficAnalysisParamsBlock();
  
  void    PlotLogo(double Logo_Height,
		   double aspectRatio,
		   double LogoX,
		   double LogoY);
    
public:
  /// constructor(s)
  Guariguanchi(const char* datacard = "DataCards/Datacard_file.txt");

  /// destructor
  ~Guariguanchi();
  
  //Analysis:
  void   DoAnalysis(void);
  
  //Do calculation of material budget vs momentum, theta and phi
  void   doMaterialBudgetAnalysis(void);
  
  //Track parameters resolutions vs momentum
  void   doTrkResolAnalysis(void);
  
  //Tracking pseudo-efficiency vs momentum
  void   doTrackingPseudoEfficVsMomentum(void);


protected:

  GGlobalTools* global;
  
  bool verbose;
  
  //The datacard
  TString DataCard;
  //List with all the geometry datacards specified
  std::vector<TString> GeometryDataCardList;
  
  //particle characteristics
  TString  particle;
  TVector3 ParticleOrigin;

  //Reference point for distance of closest approach for tracking
  TVector3 ReferencePoint;
  
  int   SinglePointMarkerSytle;
  float SinglePointMarkerSize;
  
  std::vector<double> momArr;
  std::vector<double> thetaArr;
  std::vector<double> phiArr;
  double InfiniteMon;

  TString  polarVariable;
  TString  momVariable;
  
  double TheSize;
  double TitleOffSet;
  
  TString TheOutputFile;
  
  int GlobalFileCounter;
    
  //Pseudo efficiency calculation parameters
  int    Nhits_min_sel;
  int    Nmin_Layers_track_seed_sel;
  int    Nfakes_max_sel;
  int    ndf_sel;
  double chi2Ondf_sel;
  bool   SameRange_sel;
  bool   SeedExternal_sel;
  
  double GeoCheckPrecision;
  
  //global B-field for all geometries
  GBField* GlobalBfield;
  
  bool SavePlots;
  
  bool IncludeEloss;
  
  bool DoGeometryCheck;
  bool DoPrintGeometry;
  bool DoPrintGeometryWeight;
  bool DoPlotGeometry;
  bool DoPlotWorldVolume;
  bool DoPlotStepBfieldVolume;
  bool DoRZGeoRepresentation;
  bool DoPlotSomeTracks;
  bool UseAllMomVals_GeoPrint;
  
  bool DoMaterialBugdetAnalysis;
  bool DoTelescopeAnalysis;
  bool DoTrkResolAnalysis;
  bool DoPseudoEfficVsMon;
  
  long fPrintFreq;
  
  bool  FitPowerForImpactParam;
  
  TString MonResolRepresentation;
  
  int LargeCanvasX;
  int LargeCanvasY;
  
  long MCSeed;
  
  double kappa;
  
};

#endif //~ Guariguanchi_h

