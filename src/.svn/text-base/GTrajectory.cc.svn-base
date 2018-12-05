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
#include "include/GTrajectory.h"

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
GTrajectory::GTrajectory(TString   aName,
			 TString   aParticle,
			 TVector3  aPos0,
			 TVector3  aMom0,
			 TVector3  aRefPoint,
			 GBField*  aBfield,
			 GGlobalTools* aglobal)
{
  
  Bfield = NULL;
  
  if(aBfield == NULL) {
    cout << endl;
    cout << "ERROR in GTrajectory::GTrajectory:: Input GBField* pointer is null. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Name     = aName;
  global   = aglobal;
  Type     = TString("");
  pos0     = aPos0;
  mom0     = aMom0;
  RefPoint = aRefPoint;
  SetParticle(aParticle);
  
  verbose = false;
    
}
//====================================================================
GTrajectory::GTrajectory(const GTrajectory& other,TString aName)
{
  
  Name     = aName;
  global   = other.global;
  pos0     = other.pos0;
  mom0     = other.mom0;
  RefPoint = other.RefPoint;
  SetParticle(other.Particle);
  Type     = other.Type;
  verbose  = other.verbose;
    
}
//====================================================================
GTrajectory::~GTrajectory() 
{
  
  delete Bfield;
  
  FitParams.clear();
  FitParamNames.clear();
  
}
//====================================================================
GTrajectory* GTrajectory::clone(TString aName) const
{
 
  return new GTrajectory(*this,aName);
  
}
//====================================================================
void  GTrajectory::SetParticle(TString aParticle)
{
  
  Particle = aParticle;
  charge   = global->GetParticleCharge(Particle);
  mass     = global->GetParticleMass(Particle);
  
  return;
  
}
//====================================================================
void   GTrajectory::SetInitPosition(TVector3 aInitPosition)
{
  
  pos0 = aInitPosition;
  GetFitParsFromInitConds();
  
  return;
  
}
//====================================================================
void  GTrajectory::SetInitMomentum(TVector3 aInitMomentum)
{
  
  mom0 = aInitMomentum;
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  
  return;
  
}
//====================================================================
void  GTrajectory::SetInitPositionAndMomentum(TVector3 aInitPosition, TVector3 aInitMomentum)
{
  
  pos0 = aInitPosition;
  mom0 = aInitMomentum;
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  
}
//====================================================================
void  GTrajectory::SetBfield(GBField* aBfield)
{
  
  return;
  
}
//====================================================================
TVector3  GTrajectory::GetTrueTrajectoryMon(double s)
{
  
  return mom0.Mag()*GetTrueTrajectoryUnitMon(s);
  
}
//====================================================================
double  GTrajectory::GetParameter(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameter:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParams[idx]/global->GetUnit(FitParamUnits[idx]);
  
}
//====================================================================
double  GTrajectory::GetParamError(int idx, TMatrixD  FitCovMatrix)
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
    
  return  sqrt(FitCovMatrix(idx,idx))/global->GetUnit(FitParamErrorUnits[idx]);
  
}
//====================================================================
TString  GTrajectory::GetParameterName(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameterName:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParamNames[idx];
  
}
//====================================================================
TString  GTrajectory::GetParameterUnit(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameterUnit:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParamUnits[idx];
  
}
//====================================================================
TString  GTrajectory::GetParameterUnitTitle(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameterUnitTitle:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParamUnitsTitles[idx];
  
}
//====================================================================
TString  GTrajectory::GetParameterErrorName(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameterErrorName:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParamErrorNames[idx];
  
}
//====================================================================
TString  GTrajectory::GetParameterErrorUnit(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameterErrorUnit:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParamErrorUnits[idx];
  
}
//====================================================================
TString  GTrajectory::GetParameterErrorUnitTitle(int idx) const
{
  
  if(idx < 0 || idx > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParameterErrorUnitTitle:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return  FitParamErrorUnitsTitles[idx];
  
}
//====================================================================
void  GTrajectory::GetAllParameters(std::vector<double>& aFitParams)
{
  
  aFitParams.clear();
  for(int ipar=0;ipar<int(FitParams.size());ipar++) aFitParams.push_back(FitParams[ipar]);
  
  return;
  
}
//====================================================================
void  GTrajectory::SetAllParameters(std::vector<double> aFitParams)
{
  
  FitParams.clear();
  for(int ipar=0;ipar<int(aFitParams.size());ipar++) FitParams.push_back(aFitParams[ipar]);
  
  return;
  
}
//====================================================================
double  GTrajectory::GetParamsCorrelation(int idx1, int idx2, TMatrixD  FitCovMatrix)
{
  
  if(idx1 < 0 || idx1 > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParamsCorrelation:: Parameter index " << idx1 << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(idx2 < 0 || idx2 > int(FitParams.size()-1)) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParamsCorrelation:: Parameter index " << idx1 << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(FitCovMatrix.GetNcols() != int(FitParams.size()) || FitCovMatrix.GetNrows() != int(FitParams.size())) {
    cout << endl;
    cout << "ERROR inside GTrajectory::GetParamError:: FitCovMatrix is not of size Nparameters x Nparameters, with Nparameters = " << FitParams.size() << ". Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  
  return   FitCovMatrix(idx1,idx2)/sqrt(FitCovMatrix(idx1,idx1)*FitCovMatrix(idx2,idx2));
  
}
//====================================================================
