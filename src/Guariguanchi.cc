/***************************************************************************//**
 * SVN File:    $Id$
 *
 * Project:     ILC
 *
 * @brief:      
 *
 * Description: Analytical tool for studiying tracking performances
 *
 *
 * @createdby:  PEREZ PEREZ Alejandro <luis_alejandro.perez_perez@iphc.cnrs.fr> at 2016-06-01 14:07:38
 * @copyright:  (c)2016 IPHC - CNRS - Université de Strasbourg. All Rights Reserved.
 *
 * @lastchange: $Revision$
 *              $Author$
 *              $Date$
 *
 *******************************************************************************/
 
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2D.h>
#include <TText.h>
#include <TSystem.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TVector2.h>
#include <TVector3.h>
#include "include/Guariguanchi.h"

//C++, C
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <cstring>

using namespace std;

//====================================================================
Guariguanchi::Guariguanchi(const char* datacard)
{
 
  global = new GGlobalTools();

  SinglePointMarkerSytle = 20;
  SinglePointMarkerSize  = 0.8;
  
  GlobalFileCounter = 0;
  
  Nhits_min_sel     = 4;
  Nmin_Layers_track_seed_sel = 0;
  Nfakes_max_sel    = 0;
  ndf_sel           = 2;
  chi2Ondf_sel      = 3.0;
  SameRange_sel     = true;
  SeedExternal_sel  = false;
  
  VoxelList.clear();
  GeoRanges.clear();
  GeometryList.clear();
  TrackerList.clear();
  
  GlobalBfield = NULL;
    
  GeoCheckPrecision    = 1.0*global->GetDistanceUnit("mm");
  
  SavePlots            = false;
  
  IncludeEloss         = true;
  
  DoGeometryCheck        = false;
  DoPrintGeometry        = false;
  DoPrintGeometryWeight  = false;
  DoPlotGeometry         = false;
  DoPlotWorldVolume      = false;
  DoPlotStepBfieldVolume = false;
  DoRZGeoRepresentation  = false;
  DoPlotSomeTracks       = false;
  UseAllMomVals_GeoPrint = false;
  
  DoMaterialBugdetAnalysis = false;
  DoTelescopeAnalysis = false;
  DoTrkResolAnalysis  = false;
  DoPseudoEfficVsMon  = false;

  TheSize     = 0.06;
  TitleOffSet = 1.5;
    
  TrkResolAnalysisPars.clear();
  MatBudgetAnalysisPars.clear();
  EfficAnalysisPars.clear();
  
  momArr.clear();
  thetaArr.clear();
  phiArr.clear();
  
  particle       = Dummy_name;
  ParticleOrigin = Dummy_vector;
  ReferencePoint = Dummy_vector;
  
  FitPowerForImpactParam = false;
  
  MonResolRepresentation = TString("sigma(Pt)/Pt");
  
  polarVariable          = TString("theta");
  momVariable            = TString("p");
  
  TheOutputFile = Dummy_name;
  
  DataCard = TString(datacard);
  GeometryDataCardList.clear();

  verbose = false;
  
  //LargeCanvasX = 2000;
  //LargeCanvasY = 1600;
  LargeCanvasX = 1800;
  LargeCanvasY = 1400;
  
  MCSeed = Dummy_value;
  
  kappa  = Dummy_value;
  
  //TString logo_image = TString("Guariguanchi_logo.eps");
  TString logo_image = TString("Guariguanchi_logo.png");
  logo = new TASImage(logo_image.Data());
  logo->SetImageQuality(TAttImage::kImgBest);
  
  fPrintFreq = 10;
  global->fWatch.Start();
  
}
//====================================================================
Guariguanchi::~Guariguanchi()
{

  delete global;
  //Cleaning up the memory
  for(int igeo=0;igeo<int(GeometryList.size());igeo++)  {
    if(GeometryList[igeo] != NULL) delete GeometryList[igeo];
    if(TrackerList[igeo]  != NULL) delete TrackerList[igeo];
  }
  GeometryList.clear();
  TrackerList.clear();
  if(GlobalBfield != NULL) delete GlobalBfield;
      
}
//====================================================================
void Guariguanchi::BuildSingleFile()
{
  
  if(GlobalFileCounter <= 0) return;
  
  cout << endl;
  cout << "Start Building Single File  ";
  global->fWatch.Print();
  global->fWatch.Continue();
  
  TString PDFName      = TheOutputFile + TString(".pdf");
  TString FinalPDFName = PDFName;
  
  TString command;
  command = TString("gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=") + PDFName + TString(" -dBATCH ");
  for(int ifile=0;ifile<GlobalFileCounter;ifile++) {
    PDFName = TheOutputFile + TString("_") + long(ifile+1) + TString(".eps");
    command += PDFName + TString("  ");
  }
  cout << command.Data() << endl;
  gSystem->Exec(command.Data());
  
  command = TString("rm -rf ");
  for(int ifile=0;ifile<GlobalFileCounter;ifile++) {
    PDFName = TheOutputFile + TString("_") + long(ifile+1) + TString(".eps");
    command += PDFName + TString("  ");
  }
  gSystem->Exec(command.Data());
  
  command = TString("ls -tlr ") + FinalPDFName;
  cout << endl;
  gSystem->Exec(command.Data());
  cout << endl;
  
  cout << endl;
  cout << "End Building Single File  ";
  global->fWatch.Print();
  global->fWatch.Continue();
  
  return;
  
}
//====================================================================
void  Guariguanchi::PrintTrkResolAnalysisParamsBlock()
{
  
  if(TrkResolAnalysisPars.size() == 0) {
    cout << "BeginTrkResolAnalysisParams" << endl;
    cout << "  NhitsMin                 N     // N > 0 (Mandatory)"    << endl;
    cout << "  SameRange                bool  // (Optional) default is true " << endl;
    cout << "  UseAllMomVals            bool  // (Optional) default is false" << endl;
    cout << "  PlotMaterialBudget       bool  // (Optional) default is false" << endl;
    cout << "  PlotDOCAatHighMom        bool  // (Optional) default is false" << endl;
    cout << "  PlotDOCAvsMonFit         bool  // (Optional) default is false" << endl;
    cout << "  PlotPhiAveraged          bool  // (Optional) default is false" << endl;
    cout << "  PlotOnlyPhiAveraged      bool  // (Optional) default is false" << endl;
    cout << "  PlotPerformancesVsTheta  bool  // (Optional) default is false" << endl;
    cout << "  UseLogYAxes              bool  // (Optional) default is false" << endl;
    cout << "EndTrkResolAnalysisParams"   << endl;
  }
  else {
    cout << "BeginTrkResolAnalysisParams" << endl;
    cout << "  NhitsMin                 " <<  TrkResolAnalysisPars[0].NhitsMin << endl;
    cout << "  SameRange                " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].SameRange).Data() << endl;
    cout << "  UseAllMomVals            " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].UseAllMomVals).Data() << endl;
    cout << "  PlotMaterialBudget       " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].PlotMaterialBudget).Data() << endl;
    cout << "  PlotDOCAatHighMom        " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].PlotDOCAatHighMom).Data() << endl;
    cout << "  PlotDOCAvsMonFit         " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].PlotDOCAvsMonFit).Data() << endl;
    cout << "  PlotPhiAveraged          " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].PlotPhiAveraged).Data() << endl;
    cout << "  PlotOnlyPhiAveraged      " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].PlotOnlyPhiAveraged).Data() << endl;
    cout << "  PlotPerformancesVsTheta  " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].PlotPerformancesVsTheta).Data() << endl;
    cout << "  UseLogYAxes              " <<  global->GetStringFromBool(TrkResolAnalysisPars[0].UseLogYAxes).Data() << endl;
    cout << "EndTrkResolAnalysisParams"   << endl;    
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::PrintMatBudgetAnalysisParamsBlock()
{
  
  if(MatBudgetAnalysisPars.size() == 0) {
    cout << "BeginMatBudgetAnalysisParams" << endl;
    cout << "  mom_min    val  units    // mom_min > 0 (Mandatory)"    << endl;
    cout << "  mom_max    val  units    // mom_max > 0 (Mandatory)"    << endl;
    cout << "EndMatBudgetAnalysisParams"   << endl;
    cout << "// mom_min < mom_max" << endl;
  }
  else {
    cout << "BeginMatBudgetAnalysisParams" << endl;
    cout << "  mom_min    " << MatBudgetAnalysisPars[0].mom_min/global->GetUnit("GeV/c") << " GeV/c" << endl;
    cout << "  mom_max    " << MatBudgetAnalysisPars[0].mom_max/global->GetUnit("GeV/c") << " GeV/c" << endl;
    cout << "EndMatBudgetAnalysisParams"  << endl;    
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::PrintEfficAnalysisParamsBlock()
{
  
  if(EfficAnalysisPars.size() == 0) {
    cout << "BeginEfficAnalysisParams" << endl;
    cout << "  SameRange                bool  // (Optional) default is true " << endl;
    cout << "  UseAllMomVals            bool  // (Optional) default is false" << endl;
    cout << "  PlotPhiAveraged          bool  // (Optional) default is false" << endl;
    cout << "  PlotOnlyPhiAveraged      bool  // (Optional) default is false" << endl;
    cout << "  PlotPerformancesVsTheta  bool  // (Optional) default is false" << endl;
    cout << "  UseLogYAxes              bool  // (Optional) default is false" << endl;
    cout << "EndEfficAnalysisParams"   << endl;
  }
  else {
    cout << "BeginEfficAnalysisParams" << endl;
    cout << "  SameRange                " <<  global->GetStringFromBool(EfficAnalysisPars[0].SameRange).Data() << endl;
    cout << "  UseAllMomVals            " <<  global->GetStringFromBool(EfficAnalysisPars[0].UseAllMomVals).Data() << endl;
    cout << "  PlotPhiAveraged          " <<  global->GetStringFromBool(EfficAnalysisPars[0].PlotPhiAveraged).Data() << endl;
    cout << "  PlotOnlyPhiAveraged      " <<  global->GetStringFromBool(EfficAnalysisPars[0].PlotOnlyPhiAveraged).Data() << endl;
    cout << "  PlotPerformancesVsTheta  " <<  global->GetStringFromBool(EfficAnalysisPars[0].PlotPerformancesVsTheta).Data() << endl;
    cout << "  UseLogYAxes              " <<  global->GetStringFromBool(EfficAnalysisPars[0].UseLogYAxes).Data() << endl;
    cout << "EndEfficAnalysisParams"   << endl;    
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::PlotLogo(double Logo_Height,
			     double aspectRatio,
			     double LogoX,
			     double LogoY)
{
  
  double H    = Logo_Height;
  double W    = aspectRatio * H * gPad->GetWh() * gPad->GetHNDC() / (gPad->GetWNDC() * gPad->GetWw());
  double lowX = LogoX - 0.5*W;
  double lowY = LogoY - 0.5*H;
  
  TPad *logo_pad = new TPad("logo_pad","Logo pad",lowX,lowY,lowX + W,lowY+H);
  logo_pad->Draw();
  logo_pad->cd();
  logo->Draw("X");
  
  return;
  
}
//====================================================================



