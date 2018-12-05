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
#include "include/GTrajectoryHelix.h"

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
GTrajectoryHelix::GTrajectoryHelix(TString   aName,
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
    cout << "ERROR in GTrajectoryHelix::GTrajectoryHelix:: GBField object is not GBFieldConstant. Check you inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Type   = TString("Helix");
  Bfield = NULL;
  SetBfield(aBfield);
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  FillParametersNames();
  
}
//====================================================================
GTrajectoryHelix::GTrajectoryHelix(const GTrajectoryHelix& other,TString aName)
                                   : GTrajectory(aName,
						 other.Particle,
						 other.pos0,
						 other.mom0,
						 other.RefPoint,
						 other.Bfield,
						 other.global)
{
  
  Type   = other.Type;
  SetBfield(other.Bfield);
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  FillParametersNames();
  
}
//====================================================================
GTrajectoryHelix::~GTrajectoryHelix() 
{
  
  //delete  Bfield;
  
}
//====================================================================
GTrajectory* GTrajectoryHelix::clone(TString aName) const
{
 
  return new GTrajectoryHelix(*this,aName);
  
}
//====================================================================
double  GTrajectoryHelix::GetPtAtDOCA()
{
  
  if(ConstBField.Mag() < 1.0e-6*global->GetBfieldUnit("T")) {
    cout << endl;
    cout << "ERROR in GTrajectoryHelix::GetPtAtDOCA:: ConstBField is zero. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  double alpha    = mom0.Dot(ConstBField)/pow(ConstBField.Mag(),2);
  TVector3 PtVect = alpha*ConstBField;
  PtVect          = mom0 - PtVect;
  
  return PtVect.Mag();
  
}
//====================================================================
double  GTrajectoryHelix::GetOmega()
{
  
  double pbar = mom0.Mag()/charge;
  //if(charge < 0.0) pbar *= -1.0;
  
  return (global->betaCurvRadius*ConstBField.Mag())/pbar;
  
}
//====================================================================
double  GTrajectoryHelix::GetInitialSFor2ndIntersection(double s)
{
  
  double omega = TMath::Abs(GetOmega());
  
  if(omega < 1.0/(1.0e+3*global->GetUnit("m"))) return s*10;
  else                                          return s + TMath::Pi()/omega;
  
}
//====================================================================
double  GTrajectoryHelix::GetSLimit(int nloops)
{
  
  return nloops*2.0*TMath::Pi()/TMath::Abs(GetOmega());
  
}
//====================================================================
double  GTrajectoryHelix::DummyValueCorrection(double aDummy)
{
  
  if(charge > 0) return  -TMath::Abs(aDummy);
  else           return  +TMath::Abs(aDummy);
  
}
//====================================================================
void  GTrajectoryHelix::SetBfield(GBField* aBfield)
{
  
  if(Bfield != NULL) return;
  
  Bfield = aBfield->clone(aBfield->GetName());
  ConstBField = Bfield->GetBFieldValue(TVector3(0,0,0));
 
  return;
  
}
//====================================================================
void  GTrajectoryHelix::FillParametersNames(TString PtParamFormat)
{
  
  FitParamNames.clear();
  FitParamUnits.clear();
  FitParamUnitsTitles.clear();
  
  FitParamErrorNames.clear();
  FitParamErrorUnits.clear();
  FitParamErrorUnitsTitles.clear();
  
  FitParamNames.push_back(TString("d_{#rho}"));
  FitParamUnits.push_back(TString("mm"));
  FitParamUnitsTitles.push_back(TString("mm"));
  FitParamErrorNames.push_back(TString("#sigma(d_{#rho})"));
  FitParamErrorUnits.push_back(TString("um"));
  FitParamErrorUnitsTitles.push_back(TString("#mum"));
  
  FitParamNames.push_back(TString("#phi_{0}"));
  FitParamUnits.push_back(TString("deg"));
  FitParamUnitsTitles.push_back(TString("deg"));
  FitParamErrorNames.push_back(TString("#sigma(#phi_{0})"));
  FitParamErrorUnits.push_back(TString("urad"));
  FitParamErrorUnitsTitles.push_back(TString("#murad"));

  FitParamNames.push_back(TString("d_{z}"));
  FitParamUnits.push_back(TString("mm"));
  FitParamUnitsTitles.push_back(TString("mm"));
  FitParamErrorNames.push_back(TString("#sigma(d_{z})"));
  FitParamErrorUnits.push_back(TString("um"));
  FitParamErrorUnitsTitles.push_back(TString("#mum"));
  
  FitParamNames.push_back(TString("tan#lambda"));
  FitParamUnits.push_back(TString("rad"));
  FitParamUnitsTitles.push_back(TString(""));
  FitParamErrorNames.push_back(TString("#sigma(tan#lambda)"));
  FitParamErrorUnits.push_back(TString("urad"));
  FitParamErrorUnitsTitles.push_back(TString("#murad"));
  
  FitParamNames.push_back(TString("p_{t}"));
  FitParamUnits.push_back(TString("GeV/c"));
  FitParamUnitsTitles.push_back(TString("GeV/c"));
  if(PtParamFormat == TString("sigma(Pt)/Pt")) {
    FitParamErrorNames.push_back(TString("#sigma(p_{t})/p_{t}"));
    FitParamErrorUnits.push_back(TString("rad"));
    FitParamErrorUnitsTitles.push_back(TString("%"));
  }
  else if(PtParamFormat == TString("sigma(1/Pt)")) {
    FitParamErrorNames.push_back(TString("#sigma(1/p_{t})"));
    FitParamErrorUnits.push_back(TString("1/(GeV/c)"));
    FitParamErrorUnitsTitles.push_back(TString("(GeV/c)^{-1}"));
  }
  else if(PtParamFormat == TString("sigma(Pt)")) {
    FitParamErrorNames.push_back(TString("#sigma(p_{t})"));
    FitParamErrorUnits.push_back(TString("MeV/c"));
    FitParamErrorUnitsTitles.push_back(TString("MeV/c"));
  }
  else {
    FitParamErrorNames.push_back(TString("#sigma(p_{t})/p_{t}"));
    FitParamErrorUnits.push_back(TString("rad"));
    FitParamErrorUnitsTitles.push_back(TString("%"));
  }
  
  return;
  
}
//====================================================================
double  GTrajectoryHelix::GetParamError(int idx, TMatrixD  FitCovMatrix)
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParamError:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(FitCovMatrix.GetNcols() != int(FitParams.size()) || FitCovMatrix.GetNrows() != int(FitParams.size())) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParamError:: FitCovMatrix is not of size Nparameters x Nparameters, with Nparameters = " << FitParams.size() << ". Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(idx <= 3)  return  sqrt(FitCovMatrix(idx,idx))/global->GetUnit(FitParamErrorUnits[idx]);
  else {
    if(FitParamErrorNames[idx] == TString("#sigma(p_{t})/p_{t}")) {
      return  100.0*sqrt(FitCovMatrix(idx,idx))/GetPtAtDOCA();
    }
    else if(FitParamErrorNames[idx] == TString("#sigma(1/p_{t})")) {
      double error  = sqrt(FitCovMatrix(idx,idx))/pow(GetPtAtDOCA(),2);
      error        /= global->GetUnit(FitParamErrorUnits[idx]);

      return  error;
    }
    else if(FitParamErrorNames[idx] == TString("#sigma(p_{t})")) {
      return  sqrt(FitCovMatrix(idx,idx))/global->GetUnit(FitParamErrorUnits[idx]);
    }
    else {
      return  100.0*sqrt(FitCovMatrix(idx,idx))/GetPtAtDOCA();
    }
  }
  
}
//====================================================================
double  GTrajectoryHelix::GetInitValueParForDerivativeCalculation(int idx)
{
  
  if(idx < 0 || idx > int(FitParams.size())-1) {
    cout << endl;
    cout << "ERROR inside GTrajectoryHelix::GetInitValueParForDerivativeCalculation:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  //0 -> drho
  //1 -> phi0
  //2 -> dz
  //3 -> tanLambda
  //4 -> Ptbar
    
  if(idx == 0 || idx == 2)       return  1.0*global->GetDistanceUnit("mm");
  else if(idx == 1 || idx == 3)  return  1.0*global->GetAngleUnit("mrad");
  else                           return  FitParams[4]*1.0e-2;
  
}
//====================================================================
//Equations of motion functions
void  GTrajectoryHelix::CalculatePrimedReframe()
{
  
  //Calculation primed reference frame, which depends on parcle's initial direction as well as on the magntic field (assumed constant)
  
  // - z-prime axis is along B-field,
  zprimeVect       = ConstBField.Unit();
  TVector3 UnitMon = mom0.Unit();
    
  TVector3 PerpMomentum = UnitMon - (UnitMon.Dot(zprimeVect))*zprimeVect;
    
  if(PerpMomentum.Mag() < 1.0e-8) {
    // If Initial momentum and B-field are parallel
      
    // - x-prime and y-prime are two unit vectors perpendicular to the B-field and between thenselves
    xprimeVect = zprimeVect.Orthogonal();
    xprimeVect = xprimeVect.Unit();
     
    yprimeVect = zprimeVect.Cross(xprimeVect);
  }
  else {
    // - x-prime axis along the perpendicular component of initial momentum to B-field
    xprimeVect = PerpMomentum.Unit();
    
    // - y-prime is the cross product of z-prime x x-prime director vectors
    yprimeVect = zprimeVect.Cross(xprimeVect);
  }
  
  return;
  
}
//====================================================================
double  GTrajectoryHelix::GetCurvRadius()  const
{
  
  double pbar = mom0.Mag()/charge;
  //if(charge < 0.0) pbar *= -1.0;
  double cosLambda  = mom0.Dot(xprimeVect)/mom0.Mag();
  double CurvRadius = (cosLambda*pbar)/(global->betaCurvRadius*ConstBField.Mag());
  
  return CurvRadius;
  
}
//====================================================================
double  GTrajectoryHelix::GetCosLambda()   const
{
  
  return  mom0.Dot(xprimeVect)/mom0.Mag();
  
}
//====================================================================
double  GTrajectoryHelix::GetSinLambda()   const
{
  
  return   mom0.Dot(zprimeVect)/mom0.Mag();
  
}
//====================================================================
TVector3  GTrajectoryHelix::GetTrueTrajectoryCoordinates(double s)
{
  
  // Function returns the position vector of the particle's trajectory for a value of the "s" parameter and the set of particle initial parameters
  // - particle type
  // - particle's initial position
  // - particle's initial momentum
  // - magnetic field
  
  double xprime0 = pos0.Dot(xprimeVect);
  double yprime0 = pos0.Dot(yprimeVect);
  double zprime0 = pos0.Dot(zprimeVect);
  
  double pbar = mom0.Mag()/charge;
  //if(charge < 0.0) pbar *= -1.0;
  double cosLambda  = mom0.Dot(xprimeVect)/mom0.Mag();
  double sinLambda  = mom0.Dot(zprimeVect)/mom0.Mag();
  double CurvRadius = (cosLambda*pbar)/(global->betaCurvRadius*ConstBField.Mag());
  double omega      = (global->betaCurvRadius*ConstBField.Mag())/pbar;
  
  double xprime = xprime0 + CurvRadius *  TMath::Sin(omega*s);
  double yprime = yprime0 + CurvRadius * (TMath::Cos(omega*s) - 1.0);
  double zprime = zprime0 + sinLambda*s;
  
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
TVector3  GTrajectoryHelix::GetTrueTrajectoryUnitMon(double s)
{
  
  // Function returns the momentum vector direction of the particles trajectory for a value of the "s" parameter and the set of particle initial parameters,
  // - particle type
  // - particle's initial position
  // - particle's initial momentum
  // - magnetic field
  
  double pbar = mom0.Mag()/charge;
  //if(charge < 0.0) pbar *= -1.0;
  double cosLambda  = mom0.Dot(xprimeVect)/mom0.Mag();
  double sinLambda  = mom0.Dot(zprimeVect)/mom0.Mag();
  double CurvRadius = (cosLambda*pbar)/(global->betaCurvRadius*ConstBField.Mag());
  double omega      = (global->betaCurvRadius*ConstBField.Mag())/pbar;
  
  double xprime =  omega*CurvRadius*TMath::Cos(omega*s);
  double yprime = -omega*CurvRadius*TMath::Sin(omega*s);
  double zprime =  sinLambda;
  
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
void  GTrajectoryHelix::GetFitParsFromInitConds()
{
  
  //Gets the track parameters from the initial position, momentum and magnetic field
  //In current case of constant non-zero Bfield  => helix-track: 5 parameters
  
  FitParams.clear();
  
  double xprime0    = pos0.Dot(xprimeVect);
  double yprime0    = pos0.Dot(yprimeVect);
  //double zprime0    = pos0.Dot(zprimeVect);
  
  double xprime_ref = RefPoint.Dot(xprimeVect);
  double yprime_ref = RefPoint.Dot(yprimeVect);
  double zprime_ref = RefPoint.Dot(zprimeVect);
  
  //Helix track
  // - x' = xr' + drho*cos(phi0) +           (Ptbar)/(global->betaCurvRadius*ConstBField)*(cos(phi0) - cos(phi0 + phi))
  // - y' = yr' + drho*sin(phi0) +           (Ptbar)/(global->betaCurvRadius*ConstBField)*(sin(phi0) - sin(phi0 + phi))
  // - z' = zr' + dz             - (tanLambda*Ptbar)/(global->betaCurvRadius*ConstBField)*phi
  //Parameterezation taken from http://www-glc.kek.jp/subg/offl/lib/docs/helix_manip/node3.html
  //parameters are 5
  //0 -> drho
  //1 -> phi0
  //2 -> dz
  //3 -> tanLambda
  //4 -> Ptbar

  double drho;
  double phi0;
  double dz;
  double tanLambda;
  double Ptbar;
    
  double pbar       = mom0.Mag()/charge;
  //if(charge < 0.0) pbar *= -1.0;
  double cosLambda  = mom0.Dot(xprimeVect)/mom0.Mag();
  double CurvRadius = (cosLambda*pbar)/(global->betaCurvRadius*ConstBField.Mag());
  double omega      = (global->betaCurvRadius*ConstBField.Mag())/pbar;
    
  tanLambda         = mom0.Dot(zprimeVect)/mom0.Dot(xprimeVect);
  Ptbar             = pbar*cosLambda;
    
  double xprime_c = xprime0;
  double yprime_c = yprime0 - CurvRadius;
    
  TVector3 diff_RefPoint_Center(xprime_ref - xprime_c,yprime_ref - yprime_c,0.0);
  drho = diff_RefPoint_Center.Mag() - TMath::Abs(CurvRadius);
  if(charge < 0.0) drho *= -1.0;
    
  if(charge < 0.0) phi0 = global->GetAngle( diff_RefPoint_Center.X(), diff_RefPoint_Center.Y());
  else             phi0 = global->GetAngle(-diff_RefPoint_Center.X(),-diff_RefPoint_Center.Y());
    
  TVector3 Min_aproch_point_XY_loc     = diff_RefPoint_Center.Unit()*TMath::Abs(CurvRadius);
  TVector3 ParticleOrigin_loc(xprime0 - xprime_c,yprime0 - yprime_c,0.0);
  TVector3 Unit_ParticleOrigin_loc      = ParticleOrigin_loc.Unit();
  TVector3 Unit_Perp_ParticleOrigin_loc = ParticleOrigin_loc.Unit();
  if(charge > 0.0) Unit_Perp_ParticleOrigin_loc.RotateZ(-0.5*TMath::Pi());
  else             Unit_Perp_ParticleOrigin_loc.RotateZ(+0.5*TMath::Pi());
    
  double X_local   = Min_aproch_point_XY_loc.Dot(Unit_ParticleOrigin_loc);
  double Y_local   = Min_aproch_point_XY_loc.Dot(Unit_Perp_ParticleOrigin_loc);
  double s_closest = global->GetAngle(X_local,Y_local)/TMath::Abs(omega);

  double xprime_closest = Min_aproch_point_XY_loc.X() + xprime_c;
  double yprime_closest = Min_aproch_point_XY_loc.Y() + yprime_c;
  TVector3 Closest_approach_point = GetTrueTrajectoryCoordinates(s_closest);

  double epsilon_distance = 10*global->GetUnit("nm");
  if(sqrt(pow(xprimeVect.Dot(Closest_approach_point) - xprime_closest,2) + pow(yprimeVect.Dot(Closest_approach_point) - yprime_closest,2)) > epsilon_distance) {
    cout << endl;
    cout << "ERROR in GTrajectoryHelix::GetFitParsFromInitConds: closest approach point to reference point in perpendicular plane are different." << endl;
    cout << " s_closest = " << s_closest/global->GetDistanceUnit("cm") << " cm"
         << " gives : (X,Y) = (" << xprimeVect.Dot(Closest_approach_point)/global->GetDistanceUnit("cm") << "," << yprimeVect.Dot(Closest_approach_point)/global->GetDistanceUnit("cm") << ") cm" << endl;
    cout << " calculated point gives               = (" << xprime_closest/global->GetDistanceUnit("cm") << "," << yprime_closest/global->GetDistanceUnit("cm") << ") cm" << endl;
    cout << endl;
    assert(false);
  }
  dz = zprimeVect.Dot(Closest_approach_point) - zprime_ref;
    
  double s_closest0 = s_closest;
  int Revolutions   = 20;
  for(int irev=0;irev<Revolutions;irev++) {
    double dz_tmp     = 0.0;
      
    s_closest = s_closest0 + (irev+1)*2.0*TMath::Pi()/TMath::Abs(omega);
    Closest_approach_point = GetTrueTrajectoryCoordinates(s_closest);
    dz_tmp = zprimeVect.Dot(Closest_approach_point) - zprime_ref;
    if(TMath::Abs(dz) > TMath::Abs(dz_tmp)) dz = dz_tmp;
      
    s_closest = s_closest0 - (irev+1)*2.0*TMath::Pi()/TMath::Abs(omega);
    Closest_approach_point = GetTrueTrajectoryCoordinates(s_closest);
    dz_tmp = zprimeVect.Dot(Closest_approach_point) - zprime_ref;
    if(TMath::Abs(dz) > TMath::Abs(dz_tmp)) dz = dz_tmp;
  }

  FitParams.push_back(drho);
  FitParams.push_back(phi0);
  FitParams.push_back(dz);
  FitParams.push_back(tanLambda);
  FitParams.push_back(Ptbar);
    
  if(verbose) PrintParameters();
  
  return;
  
}
//====================================================================
TVector3  GTrajectoryHelix::GetFitTrackCoordinates(double dummy, double sinit)
{
  
  //Helix track
  // - x' = xr' + drho*cos(phi0) +           (Ptbar)/(global->betaCurvRadius*MyBfield)*(cos(phi0) - cos(phi0 + phi))
  // - y' = yr' + drho*sin(phi0) +           (Ptbar)/(global->betaCurvRadius*MyBfield)*(sin(phi0) - sin(phi0 + phi))
  // - z' = zr' + dz             - (tanLambda*Ptbar)/(global->betaCurvRadius*MyBfield)*phi
  //Dummy parameter is phi
  //parameters are 5
  //0 -> drho
  //1 -> phi0
  //2 -> dz
  //3 -> tanLambda
  //4 -> Ptbar
    
  double drho      = FitParams[0];
  double phi0      = FitParams[1];
  double dz        = FitParams[2];
  double tanLambda = FitParams[3];
  double Ptbar     = FitParams[4];
    
  double xprime_ref = RefPoint.Dot(xprimeVect);
  double yprime_ref = RefPoint.Dot(yprimeVect);
  double zprime_ref = RefPoint.Dot(zprimeVect);
    
  double CurvRadius = Ptbar/(global->betaCurvRadius*ConstBField.Mag());
    
  double xprime = xprime_ref + drho*TMath::Cos(phi0) + CurvRadius*(TMath::Cos(phi0) - TMath::Cos(phi0 + dummy));
  double yprime = yprime_ref + drho*TMath::Sin(phi0) + CurvRadius*(TMath::Sin(phi0) - TMath::Sin(phi0 + dummy));
  double zprime = zprime_ref + dz                    - CurvRadius*tanLambda*dummy;
    
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
TVector3  GTrajectoryHelix::GetFitTrackMomentum(double dummy, double sinit)
{
  
  //Helix track
  // - px = -|Ptbar| * sin(phi0 + phi)
  // - py = +|Ptbar| * cos(phi0 + phi)
  // - pz = +|Ptbar| * tanLambda
  //Dummy parameter is phi
  //parameters are 5
  //0 -> drho
  //1 -> phi0
  //2 -> dz
  //3 -> tanLambda
  //4 -> Ptbar

  double phi0      = FitParams[1];
  double tanLambda = FitParams[3];
  double Ptbar     = FitParams[4];

  double xprime = -TMath::Abs(Ptbar)*TMath::Sin(phi0 + dummy);
  double yprime = +TMath::Abs(Ptbar)*TMath::Cos(phi0 + dummy);
  double zprime = +TMath::Abs(Ptbar)*tanLambda;
    
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
  
}
//====================================================================
double  GTrajectoryHelix::GetFitTrackDummyParFromS(double s)
{
  
  //Get the dummy track transport parameter from the s parameter, track initial position and momentum, and track parameters
  //Take a look of the function GetFitTrackCoordinates to understand with the dummy parameter is
  
  TVector3 pos_true = GetTrueTrajectoryCoordinates(s);
  
  //Helix track
  // - x' = xr' + drho*cos(phi0) -             (Ptbar)/(global->betaCurvRadius*MyBfield)*(cos(phi0) - cos(phi0 + phi))
  // - y' = yr' + drho*sos(phi0) -             (Ptbar)/(global->betaCurvRadius*MyBfield)*(sin(phi0) - sin(phi0 + phi))
  // - z' = zr' + dz             - (tan(lambda)*Ptbar)/(global->betaCurvRadius*MyBfield)*phi
  //Dummy parameter is phi
  //parameters are 5
  //0 -> drho
  //1 -> phi0
  //2 -> dz
  //3 -> tanLambda
  //4 -> Ptbar
  TString units_ttt = TString("rad");
  double dummy;
    
  double Ptbar      = FitParams[4];
  double CurvRadius = (Ptbar)/(global->betaCurvRadius*ConstBField.Mag());
    
  double xprime_true = pos_true.Dot(xprimeVect);
  double yprime_true = pos_true.Dot(yprimeVect);
    
  double xprime_ref  = RefPoint.Dot(xprimeVect);
  double yprime_ref  = RefPoint.Dot(yprimeVect);
  double xprime0     = pos0.Dot(xprimeVect);
  double yprime0     = pos0.Dot(yprimeVect);
  double xprime_c    = xprime0;
  double yprime_c    = yprime0 - CurvRadius;
    
  TVector3 Unit_X(xprime_ref - xprime_c,yprime_ref - yprime_c,0.0);
  Unit_X = Unit_X.Unit();
  if(charge > 0.0) Unit_X.RotateZ(TMath::Pi());
  TVector3 Unit_Y = Unit_X;
  if(charge > 0.0) Unit_Y.RotateZ(-0.5*TMath::Pi());
  else             Unit_Y.RotateZ(+0.5*TMath::Pi());
  Unit_Y = Unit_Y.Unit();
    
  TVector3 pos_true_locXY(xprime_true - xprime_c,yprime_true - yprime_c,0.0);
    
  double X = pos_true_locXY.Dot(Unit_X);
  double Y = pos_true_locXY.Dot(Unit_Y);
  if(charge > 0.0) {
    X *= -1.0;
    Y *= -1.0;
  }
    
  dummy = global->GetAngle(X,Y);
  if(charge > 0.0) dummy *= -1.0;
    
  TVector3 ppp = GetFitTrackCoordinates(dummy);
  double my_epsilon_ppp = 1.0e-6*global->GetDistanceUnit("mm");
    
  if((ppp - pos_true).Mag() > my_epsilon_ppp) {
    double dummy0   = dummy;
    int Sign = -1;
    if(charge > 0.0) Sign *= -1;
    int Revolutions = 100000;
    for(int irev=0;irev<Revolutions;irev++) {
      double dummy_tmp;
      
      dummy_tmp = dummy0 + Sign*(irev+1)*2.0*TMath::Pi();
      ppp  = GetFitTrackCoordinates(dummy_tmp);
      if((ppp - pos_true).Mag() < my_epsilon_ppp) {
	dummy = dummy_tmp;
	break;
      }
        
      dummy_tmp = dummy0 - Sign*(irev+1)*2.0*TMath::Pi();
      ppp  = GetFitTrackCoordinates(dummy_tmp);
      if((ppp - pos_true).Mag() < my_epsilon_ppp) {
	dummy = dummy_tmp;
	break;
      }
      
    }
  }
  
  TVector3 pos_fit  = GetFitTrackCoordinates(dummy);
  TVector3 pos_diff = pos_true - pos_fit;
  if(pos_diff.Mag() > 1.0e-6) {
    cout << endl;
    cout << "ERROR in GTrajectoryHelix::GetFitTrackDummyParFromS:  track position from true function and fit function with parameters s = " << s/global->GetDistanceUnit("cm") 
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
TVector3  GTrajectoryHelix::GetFitTrackCoorDerWRTDummy(double dummy, double sinit)
{
  
  //Helix track
  // - x' = xr' + drho*cos(phi0) +           (Ptbar)/(global->betaCurvRadius*MyBfield)*(cos(phi0) - cos(phi0 + phi))
  // - y' = yr' + drho*sin(phi0) +           (Ptbar)/(global->betaCurvRadius*MyBfield)*(sin(phi0) - sin(phi0 + phi))
  // - z' = zr' + dz             - (tanLambda*Ptbar)/(global->betaCurvRadius*MyBfield)*phi
  //Dummy parameter is phi
  //parameters are 5
  //0 -> drho
  //1 -> phi0
  //2 -> dz
  //3 -> tanLambda
  //4 -> Ptbar
    
  double phi0      = FitParams[1];
  double tanLambda = FitParams[3];
  double Ptbar     = FitParams[4];
    
  double CurvRadius = Ptbar/(global->betaCurvRadius*ConstBField.Mag());
    
  double xprime = +CurvRadius*TMath::Sin(phi0 + dummy);
  double yprime = -CurvRadius*TMath::Cos(phi0 + dummy);
  double zprime = -CurvRadius*tanLambda;
  
  return  xprime*xprimeVect + yprime*yprimeVect + zprime*zprimeVect;
    
}
//====================================================================
void  GTrajectoryHelix::PrintParameters()
{
  
  cout << endl;
  cout << "Parameters: " << endl;
  cout << " - drho       = " << FitParams[0]/global->GetDistanceUnit("cm")     << " cm"    << endl;
  cout << " - phi0       = " << FitParams[1]/global->GetAngleUnit("deg")       << " deg"   << endl;
  cout << " - dz         = " << FitParams[2]/global->GetDistanceUnit("cm")     << " cm"    << endl;
  cout << " - tanLambda  = " << FitParams[3]                                   << ""       << endl;
  cout << " - Ptbar      = " << FitParams[4]/global->GetMomentumUnit("GeV/c")  << " GeV/c" << endl;
  cout << endl;
  
  return;
  
}
//====================================================================
