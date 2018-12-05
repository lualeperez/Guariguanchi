#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2D.h>
#include <TText.h>
#include <TSystem.h>
#include <TLine.h>
#include <TFile.h>
#include <TEllipse.h>
#include <TVector2.h>
#include <TVector3.h>
#include "include/GTrajectoryStraight.h"

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
GTrajectoryStraight::GTrajectoryStraight(TString   aName,
					 TString   aParticle,
					 TVector3  aPos0,
					 TVector3  aMom0,
					 TVector3  aRefPoint,
					 GBField*  aBfield,
					 GGlobalTools* aglobal)
                                         : GTrajectory(aName,
						       aParticle,
						       aPos0,
						       aMom0,
						       aRefPoint,
						       aBfield,
						       aglobal)
{
 
  if(aBfield->GetType() != TString("Constant")) {
    cout << endl;
    cout << "ERROR in GTrajectoryStraight::GTrajectoryStraight:: GBField object is not GBFieldConstant. Check you inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Type   = TString("Straight");
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  FillParametersNames();
  
}
//====================================================================
GTrajectoryStraight::GTrajectoryStraight(const GTrajectoryStraight& other,TString aName)
                                         : GTrajectory(aName,
						       other.Particle,
						       other.pos0,
						       other.mom0,
						       other.RefPoint,
						       other.Bfield,
						       other.global)
{
  
  Type   = other.Type;
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  FillParametersNames();
  
}
//====================================================================
GTrajectoryStraight::~GTrajectoryStraight() 
{
  
}
//====================================================================
GTrajectory* GTrajectoryStraight::clone(TString aName) const
{
 
  return new GTrajectoryStraight(*this,aName);
  
}
//====================================================================
void  GTrajectoryStraight::FillParametersNames(TString PtParamFormat)
{
  
  FitParamNames.clear();
  FitParamUnits.clear();
  FitParamUnitsTitles.clear();
  
  FitParamErrorNames.clear();
  FitParamErrorUnits.clear();
  FitParamErrorUnitsTitles.clear();
  
  FitParamNames.push_back(TString("tanX"));
  FitParamUnits.push_back(TString("mrad"));
  FitParamUnitsTitles.push_back(TString("mrad"));
  FitParamErrorNames.push_back(TString("#sigma(tan(#alpha_{x}))"));
  FitParamErrorUnits.push_back(TString("urad"));
  FitParamErrorUnitsTitles.push_back(TString("#murad"));
  
  FitParamNames.push_back(TString("x_{0}"));
  FitParamUnits.push_back(TString("mm"));
  FitParamUnitsTitles.push_back(TString("mm"));
  FitParamErrorNames.push_back(TString("#sigma(x_{0})"));
  FitParamErrorUnits.push_back(TString("um"));
  FitParamErrorUnitsTitles.push_back(TString("#mum"));
  
  FitParamNames.push_back(TString("tanY"));
  FitParamUnits.push_back(TString("mrad"));
  FitParamUnitsTitles.push_back(TString("mrad"));
  FitParamErrorNames.push_back(TString("#sigma(tan(#alpha_{y}))"));
  FitParamErrorUnits.push_back(TString("urad"));
  FitParamErrorUnitsTitles.push_back(TString("#murad"));
  
  FitParamNames.push_back(TString("y_{0}"));
  FitParamUnits.push_back(TString("mm"));
  FitParamUnitsTitles.push_back(TString("mm"));
  FitParamErrorNames.push_back(TString("#sigma(y_{0})"));
  FitParamErrorUnits.push_back(TString("um"));
  FitParamErrorUnitsTitles.push_back(TString("#mum"));
  
  return;
  
}
//====================================================================
double  GTrajectoryStraight::GetInitValueParForDerivativeCalculation(int idx)
{
  
  if(idx < 0 || idx > int(FitParams.size())-1) {
    cout << endl;
    cout << "ERROR inside GTrajectoryStraight::GetInitValueParForDerivativeCalculation:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(idx == 0 || idx == 2) return  global->GetAngleUnit("mrad");
  else                     return  global->GetDistanceUnit("mm");
  
}
//====================================================================
void  GTrajectoryStraight::CalculatePrimedReframe()
{
  
  //Calculation primed reference frame, which depends on parcle's initial direction as well as on the magntic field (assumed constant)
  
  TVector3 InitMomentumUnit = mom0.Unit();
  double angleX = TMath::ACos(InitMomentumUnit.Dot(TVector3(1.0,0.0,0.0)));
  double angleY = TMath::ACos(InitMomentumUnit.Dot(TVector3(0.0,1.0,0.0)));
  double angleZ = TMath::ACos(InitMomentumUnit.Dot(TVector3(0.0,0.0,1.0)));
    
  double angle_min = 1.0e+20;
  if(angle_min > angleX) angle_min = angleX;
  if(angle_min > angleY) angle_min = angleY;
  if(angle_min > angleZ) angle_min = angleZ;
    
  if(TMath::Abs(angle_min - angleX) < 1.0e-8) {
    zprimeVect = TVector3(1.0,0.0,0.0);
    xprimeVect = TVector3(0.0,1.0,0.0);
    yprimeVect = TVector3(0.0,0.0,1.0);
  }
  else if(TMath::Abs(angle_min - angleY) < 1.0e-8) {
    zprimeVect = TVector3(0.0,1.0,0.0);
    xprimeVect = TVector3(0.0,0.0,1.0);
    yprimeVect = TVector3(1.0,0.0,0.0);
  }
  else if(TMath::Abs(angle_min - angleZ) < 1.0e-8) {
    zprimeVect = TVector3(0.0,0.0,1.0);
    xprimeVect = TVector3(1.0,0.0,0.0);
    yprimeVect = TVector3(0.0,1.0,0.0);
  }
  
  return;
  
}
//====================================================================
TVector3  GTrajectoryStraight::GetTrueTrajectoryCoordinates(double s)
{
  
  // Function returns the position vector of the particle's trajectory for a value of the "s" parameter and the set of particle initial parameters
  // - particle type
  // - particle's initial position
  // - particle's initial momentum
  // - magnetic field
  
  double xprime0 = pos0.Dot(xprimeVect);
  double yprime0 = pos0.Dot(yprimeVect);
  double zprime0 = pos0.Dot(zprimeVect);
  
  TVector3 InitMomentumUnit = mom0.Unit();
  double xprime = xprime0 + InitMomentumUnit.Dot(xprimeVect)*s;
  double yprime = yprime0 + InitMomentumUnit.Dot(yprimeVect)*s;
  double zprime = zprime0 + InitMomentumUnit.Dot(zprimeVect)*s;
  
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
TVector3  GTrajectoryStraight::GetTrueTrajectoryUnitMon(double s)
{
  
  // Function returns the momentum vector direction of the particles trajectory for a value of the "s" parameter and the set of particle initial parameters,
  // - particle type
  // - particle's initial position
  // - particle's initial momentum
  // - magnetic field
  
  TVector3 InitMomentumUnit = mom0.Unit();
  double  xprime = InitMomentumUnit.Dot(xprimeVect);
  double  yprime = InitMomentumUnit.Dot(yprimeVect);
  double  zprime = InitMomentumUnit.Dot(zprimeVect);
  
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
void  GTrajectoryStraight::GetFitParsFromInitConds()
{
  
  //Gets the track parameters from the initial position, momentum and magnetic field
  //In current case of zero Bfield  => stright-track: 4 parameters
  
  FitParams.clear();
  
  double xprime0    = pos0.Dot(xprimeVect);
  double yprime0    = pos0.Dot(yprimeVect);
  double zprime0    = pos0.Dot(zprimeVect);
  
  //double xprime_ref = RefPoint.Dot(xprimeVect);
  //double yprime_ref = RefPoint.Dot(yprimeVect);
  double zprime_ref = RefPoint.Dot(zprimeVect);
  
  //Straight line track: => 4 parameters
    
  // - x = (px/pz)*(z - z_ref) + x0;
  // - y = (py/pz)*(z - z_ref) + y0;
  //parameters are:
  //0 -> px/pz
  //1 -> x0
  //2 -> py/pz
  //3 -> y0
  double pxOpz = mom0.Dot(xprimeVect)/mom0.Dot(zprimeVect);
  double x0    = xprime0 + pxOpz*(zprime_ref - zprime0);
  double pyOpz = mom0.Dot(yprimeVect)/mom0.Dot(zprimeVect);
  double y0    = yprime0 + pyOpz*(zprime_ref - zprime0);
  FitParams.push_back(pxOpz);
  FitParams.push_back(x0);
  FitParams.push_back(pyOpz);
  FitParams.push_back(y0);
    
  if(verbose) PrintParameters();
  
  return;
  
}
//====================================================================
TVector3  GTrajectoryStraight::GetFitTrackCoordinates(double dummy, double sinit)
{
  
  //Straight line track
    
  // - x = (px/pz)*(z-zr) + x0;
  // - y = (py/pz)*(z-zr) + y0;
  //Dummy parameter is z
  //parameters are four
  //0 -> px/pz
  //1 -> x0
  //2 -> py/pz
  //3 -> y0
     
  double pxOpz  = FitParams[0];
  double x0     = FitParams[1];
  double pyOpz  = FitParams[2];
  double y0     = FitParams[3];
  double dummy0 = RefPoint.Dot(zprimeVect);
      
  double xprime = pxOpz*(dummy - dummy0) + x0;
  double yprime = pyOpz*(dummy - dummy0) + y0;
  double zprime = dummy;
    
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
TVector3  GTrajectoryStraight::GetFitTrackMomentum(double dummy, double sinit)
{
  
  //Straight line track
    
  // - pz = sign(pz)*p/sqrt(1 + (px/pz)**2 + (py/pz)**2)
  // - px = (px/pz)*pz;
  // - py = (py/pz)*pz;
  //parameters are four
  //0 -> px/pz
  //1 -> x0
  //2 -> py/pz
  //3 -> y0
      
  double pxOpz  = FitParams[0];
  double pyOpz  = FitParams[2];
    
  double PzSign = mom0.Dot(zprimeVect)/TMath::Abs(mom0.Dot(zprimeVect));
    
  double zprime = PzSign*mom0.Mag()/sqrt(1.0 + pow(pxOpz,2) + pow(pyOpz,2));
  double xprime = pxOpz*zprime;
  double yprime = pyOpz*zprime;
    
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
double  GTrajectoryStraight::GetFitTrackDummyParFromS(double s)
{
  
  //Get the dummy track transport parameter from the s parameter, track initial position and momentum, and track parameters
  //Take a look of the function GetFitTrackCoordinates to understand with the dummy parameter is
  
  TVector3 pos_true = GetTrueTrajectoryCoordinates(s);
  
  //Straight line track, four parameters
  
  TString units_ttt = TString("mm");
  double dummy = pos_true.Dot(zprimeVect);
  
  TVector3 pos_fit  = GetFitTrackCoordinates(dummy);
  TVector3 pos_diff = pos_true - pos_fit;
  if(pos_diff.Mag() > 1.0e-6) {
    cout << endl;
    cout << "ERROR in GTrajectoryStraight::GetFitTrackDummyParFromS:  track position from true function and fit function with parameters s = " << s/global->GetDistanceUnit("cm") 
         << " cm and dummy = " << dummy/global->GetUnit(units_ttt) << " " << units_ttt.Data() << " give different positions." << endl;
    cout << "True position = (" << pos_true.X()/global->GetDistanceUnit("cm") << "," << pos_true.Y()/global->GetDistanceUnit("cm") << "," << pos_true.Z()/global->GetDistanceUnit("cm") << ") cm" << endl;
    cout << "Fit  position = (" << pos_fit.X()/global->GetDistanceUnit("cm")  << "," << pos_fit.Y()/global->GetDistanceUnit("cm")  << "," << pos_fit.Z()/global->GetDistanceUnit("cm")  << ") cm" << endl;
    cout << "diff position = (" << pos_diff.X()/global->GetDistanceUnit("cm") << "," << pos_diff.Y()/global->GetDistanceUnit("cm") << "," << pos_diff.Z()/global->GetDistanceUnit("cm") << ") cm" << endl;
    cout << endl;
    assert(false);
  }

  return  dummy;
  
}
//====================================================================
TVector3  GTrajectoryStraight::GetFitTrackCoorDerWRTDummy(double dummy, double sinit)
{
  
  //Straight line track
    
  // - x = (px/pz)*(z-zr) + x0;
  // - y = (py/pz)*(z-zr) + y0;
  //Dummy parameter is z
  //parameters are four
  //0 -> px/pz
  //1 -> x0
  //2 -> py/pz
  //3 -> y0
      
  double pxOpz  = FitParams[0];
  double pyOpz  = FitParams[2];
      
  double xprime = pxOpz;
  double yprime = pyOpz;
  double zprime = 1.0;
  
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
    
}
//====================================================================
void  GTrajectoryStraight::PrintParameters()
{
  
  cout << endl;
  cout << "Parameters: " << endl;
  cout << " - x0   = " << FitParams[0]/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << " - tanX = " << FitParams[1]          << endl;
  cout << " - y0   = " << FitParams[2]/global->GetDistanceUnit("cm") << " cm" << endl;
  cout << " - tanY = " << FitParams[3]          << endl;
  cout << endl;
  
  return;
  
}
//====================================================================

