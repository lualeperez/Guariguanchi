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
#include "include/GTrajectoryMultipleSteps.h"

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
GTrajectoryMultipleSteps::GTrajectoryMultipleSteps(TString   aName,
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
  
  if(aBfield->GetType() != TString("MultipleSteps")) {
    cout << endl;
    cout << "ERROR in GTrajectoryMultipleSteps::GTrajectoryMultipleSteps:: GBField object is not GBFieldMultipleSteps. Check you inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  Type   = TString("MultipleSteps");
  Bfield = NULL;
  SetBfield(aBfield);
  FillTrajectories();

  CalculatePrimedReframe();
  FillParametersNames();  
  GetFitParsFromInitConds();
  
}
//====================================================================
GTrajectoryMultipleSteps::GTrajectoryMultipleSteps(const GTrajectoryMultipleSteps& other,TString aName)
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
  FillTrajectories();
  
  CalculatePrimedReframe();
  GetFitParsFromInitConds();
  FillParametersNames();
  
}
//====================================================================
GTrajectoryMultipleSteps::~GTrajectoryMultipleSteps() 
{
  
  //delete  Bfield;
  
  for(int i=0;i<int(TrueInsideVolTrajectoryList.size());i++) {
    if(TrueInsideVolTrajectoryList[i] != NULL) delete  TrueInsideVolTrajectoryList[i];
  }
  for(int i=0;i<int(TrueOutsideVolTrajectoryList.size());i++) {
    if(TrueOutsideVolTrajectoryList[i] != NULL) delete  TrueOutsideVolTrajectoryList[i];    
  }
  TrueInsideVolTrajectoryList.clear();
  TrueOutsideVolTrajectoryList.clear();
  for(int i=0;i<int(FitInsideVolTrajectoryList.size());i++) {
    if(FitInsideVolTrajectoryList[i] != NULL) delete  FitInsideVolTrajectoryList[i];
  }
  for(int i=0;i<int(FitOutsideVolTrajectoryList.size());i++) {
    if(FitOutsideVolTrajectoryList[i] != NULL) delete  FitOutsideVolTrajectoryList[i];
  }
  FitInsideVolTrajectoryList.clear();
  FitOutsideVolTrajectoryList.clear();
  
  VolumeList.clear();
  
}
//====================================================================
GTrajectory* GTrajectoryMultipleSteps::clone(TString aName) const
{
 
  return new GTrajectoryMultipleSteps(*this,aName);
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetPtAtDOCA()
{

  return mom0.Mag();
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetOmega()
{
  
  return  Dummy_value;
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetInitialSFor2ndIntersection(double s)
{
  
  //if(MaxOmega < 1.0/(1.0e+3*global->GetUnit("m"))) return s*10;
  //else                                             return s + TMath::Pi()/MaxOmega;
  
  return 10*s;
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetSLimit(int nloops)
{
    
  return  1.0e+3*global->GetDistanceUnit("m");
  
}
//====================================================================
double  GTrajectoryMultipleSteps::DummyValueCorrection(double aDummy)
{
  
  return  aDummy;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::SetBfield(GBField* aBfield)
{
  
  if(Bfield != NULL) return;
  Bfield = aBfield->clone(aBfield->GetName());
 
  return;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::FillTrajectories(void)
{
  
  VolumeList.clear();
  TrueInsideVolTrajectoryList.clear();
  TrueOutsideVolTrajectoryList.clear();
  FitInsideVolTrajectoryList.clear();
  FitOutsideVolTrajectoryList.clear();
  
  small_distance = 1.0*global->GetUnit("um");
  large_distance = 1.0e+3*global->GetUnit("m");
  
  double small_Bfield = 1.0e-6*global->GetUnit("T");
  
  TheBField = dynamic_cast<GBFieldMultipleSteps*>(Bfield);
  
  TString NameTmp;
  NameTmp = TString("Constant B-field outside volumes");
  GBFieldConstant  OutBField(NameTmp,TheBField->GetOutBField(),global);
  
  int counter_outside = 0;
  Noutside_plus = 4;
  Ninside_plus  = 3;
  for(int kkk=0;kkk<Ninside_plus;kkk++) {
    for(int ivol=0;ivol<TheBField->GetNVolumes();ivol++) {
      VolumeList.push_back(TheBField->GetVolume(ivol));
    
      NameTmp = TString("Constant B-field inside volume ") + long(ivol+1);
      GBFieldConstant  InBField(NameTmp,TheBField->GetInBField(ivol),global);
      
      if(TheBField->GetInBField(ivol).Mag() > small_Bfield) {
        NameTmp = TString("True trajectory inside volume ") + long(ivol+1);
        TrueInsideVolTrajectoryList.push_back(new GTrajectoryHelix(NameTmp,
							           Particle,
							           pos0,mom0,
							           RefPoint,
							           &InBField,
							           global));
    
        NameTmp = TString("Fit trajectory inside volume ") + long(ivol+1);
        FitInsideVolTrajectoryList.push_back(new GTrajectoryHelix(NameTmp,
							          Particle,
							          pos0,mom0,
							          RefPoint,
							          &InBField,
							          global));
      }
      else {
        NameTmp = TString("True trajectory inside volume ") + long(ivol+1);
        TrueInsideVolTrajectoryList.push_back(new GTrajectoryStraight(NameTmp,
							              Particle,
							              pos0,mom0,
							              RefPoint,
							              &InBField,
							              global));
    
        NameTmp = TString("Fit trajectory inside volume ") + long(ivol+1);
        FitInsideVolTrajectoryList.push_back(new GTrajectoryStraight(NameTmp,
							             Particle,
							             pos0,mom0,
							             RefPoint,
							             &InBField,
							             global));
      }
    }
  }
   
  for(int ivol=0;ivol<TheBField->GetNVolumes();ivol++) {
    for(int kkk=0;kkk<Noutside_plus;kkk++) {
      if(TheBField->GetOutBField().Mag() > small_Bfield) {
        NameTmp = TString("True trajectory outside volumes ") + long(counter_outside+1);
        TrueOutsideVolTrajectoryList.push_back(new GTrajectoryHelix(NameTmp,
								    Particle,
								    pos0,mom0,
								    RefPoint,
								    &OutBField,
								    global));
    
        NameTmp = TString("Fit trajectory outside volumes ") + long(counter_outside+1);
        FitOutsideVolTrajectoryList.push_back(new GTrajectoryHelix(NameTmp,
							           Particle,
							           pos0,mom0,
							           RefPoint,
							           &OutBField,
							           global));
      }
      else {
	NameTmp = TString("True trajectory outside volumes ") + long(counter_outside+1);
        TrueOutsideVolTrajectoryList.push_back(new GTrajectoryStraight(NameTmp,
								       Particle,
								       pos0,mom0,
								       RefPoint,
								       &OutBField,
								       global));
    
        NameTmp = TString("Fit trajectory outside volumes ") + long(counter_outside+1);
        FitOutsideVolTrajectoryList.push_back(new GTrajectoryStraight(NameTmp,
							              Particle,
							              pos0,mom0,
							              RefPoint,
							              &OutBField,
							              global));
      }
      counter_outside++;
    }
    
  }
  
  MaxOmega = -1.0e+10;
  for(int i=0;i<int(TrueInsideVolTrajectoryList.size());i++) {
    double omega = TMath::Abs(TrueInsideVolTrajectoryList[i]->GetOmega());
    if(MaxOmega < omega) MaxOmega = omega;
  }
  double omega = TMath::Abs(TrueOutsideVolTrajectoryList[0]->GetOmega());
  if(MaxOmega < omega) MaxOmega = omega;
    
  return;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::FillParametersNames(TString PtParamFormat)
{
  
  FitParamNames.clear();
  FitParamUnits.clear();
  FitParamUnitsTitles.clear();
  
  FitParamErrorNames.clear();
  FitParamErrorUnits.clear();
  FitParamErrorUnitsTitles.clear();
  
  FitParamNames.push_back(TString("x_{0}"));
  FitParamUnits.push_back(TString("mm"));
  FitParamUnitsTitles.push_back(TString("mm"));
  FitParamErrorNames.push_back(TString("#sigma(x_{0})"));
  FitParamErrorUnits.push_back(TString("um"));
  FitParamErrorUnitsTitles.push_back(TString("#mum"));
  
  FitParamNames.push_back(TString("y_{0}"));
  FitParamUnits.push_back(TString("mm"));
  FitParamUnitsTitles.push_back(TString("mm"));
  FitParamErrorNames.push_back(TString("#sigma(y_{0})"));
  FitParamErrorUnits.push_back(TString("um"));
  FitParamErrorUnitsTitles.push_back(TString("#mum"));

  FitParamNames.push_back(TString("t^{0}_{x}"));
  FitParamUnits.push_back(TString("rad"));
  FitParamUnitsTitles.push_back(TString(""));
  FitParamErrorNames.push_back(TString("#sigma(t^{0}_{x})"));
  FitParamErrorUnits.push_back(TString("urad"));
  FitParamErrorUnitsTitles.push_back(TString("#murad"));
  
  FitParamNames.push_back(TString("t^{0}_{y}"));
  FitParamUnits.push_back(TString("rad"));
  FitParamUnitsTitles.push_back(TString(""));
  FitParamErrorNames.push_back(TString("#sigma(t^{0}_{y})"));
  FitParamErrorUnits.push_back(TString("urad"));
  FitParamErrorUnitsTitles.push_back(TString("#murad"));
  
  FitParamNames.push_back(TString("p"));
  FitParamUnits.push_back(TString("GeV/c"));
  FitParamUnitsTitles.push_back(TString("GeV/c"));
  if(PtParamFormat == TString("sigma(Pt)/Pt")) {
    FitParamErrorNames.push_back(TString("#sigma(p)/p"));
    FitParamErrorUnits.push_back(TString("rad"));
    FitParamErrorUnitsTitles.push_back(TString("%"));
  }
  else if(PtParamFormat == TString("sigma(1/Pt)")) {
    FitParamErrorNames.push_back(TString("#sigma(1/p)"));
    FitParamErrorUnits.push_back(TString("1/(GeV/c)"));
    FitParamErrorUnitsTitles.push_back(TString("(GeV/c)^{-1}"));
  }
  else if(PtParamFormat == TString("sigma(Pt)")) {
    FitParamErrorNames.push_back(TString("#sigma(p)"));
    FitParamErrorUnits.push_back(TString("MeV/c"));
    FitParamErrorUnitsTitles.push_back(TString("MeV/c"));
  }
  else {
    FitParamErrorNames.push_back(TString("#sigma(p)/p"));
    FitParamErrorUnits.push_back(TString("rad"));
    FitParamErrorUnitsTitles.push_back(TString("%"));
  }
  
  return;
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetParamError(int idx, TMatrixD  FitCovMatrix)
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
      return  100.0*sqrt(FitCovMatrix(idx,idx))/mom0.Mag();
    }
    else if(FitParamErrorNames[idx] == TString("#sigma(1/p_{t})")) {
      double error  = sqrt(FitCovMatrix(idx,idx))/pow(mom0.Mag(),2);
      error        /= global->GetUnit(FitParamErrorUnits[idx]);

      return  error;
    }
    else if(FitParamErrorNames[idx] == TString("#sigma(p_{t})")) {
      return  sqrt(FitCovMatrix(idx,idx))/global->GetUnit(FitParamErrorUnits[idx]);
    }
    else {
      return  100.0*sqrt(FitCovMatrix(idx,idx))/mom0.Mag();
    }
  }
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetInitValueParForDerivativeCalculation(int idx)
{
  
  if(idx < 0 || idx > int(FitParams.size())-1) {
    cout << endl;
    cout << "ERROR inside GTrajectoryMultipleSteps::GetInitValueParForDerivativeCalculation:: Parameter index " << idx << " is out of limits (" << 0 << "," << FitParams.size()-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  //0 -> x0
  //1 -> y0
  //2 -> tx
  //3 -> ty
  //4 -> Pbar
    
  if(idx == 0 || idx == 1)       return  1.0e-2*global->GetDistanceUnit("mm");
  else if(idx == 2 || idx == 3)  return  1.0e-2*global->GetAngleUnit("mrad");
  else                           return  FitParams[4]*1.0e-3;
  
}
//====================================================================
//Equations of motion functions
void  GTrajectoryMultipleSteps::CalculatePrimedReframe()
{
  
  //Calculation primed reference frame, which depends on parcle's initial direction as well as on the magntic field (assumed constant)
  
  xprimeVect = TVector3(1.0,0.0,0.0);
  yprimeVect = TVector3(0.0,1.0,0.0);
  zprimeVect = TVector3(0.0,0.0,1.0);
  
  return;
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetTrueTrajectoryCoordinates(double s)
{
  
  // Function returns the position vector of the particle's trajectory for a value of the "s" parameter and the set of particle initial parameters
  // - particle type
  // - particle's initial position
  // - particle's initial momentum
  // - magnetic field
  
  for(int i=0;i<int(SInterVals_true.size());i++) {
    if(s >= SInterVals_true[i](0) && s <= SInterVals_true[i](1)) {
      if(TrueTrackOrderList[i].IsInsideVolTrack) return  TrueInsideVolTrajectoryList[TrueTrackOrderList[i].idx_InVolTrack]->GetTrueTrajectoryCoordinates(s - TrueTrackOrderList[i].s_referece);
      else                                       return  TrueOutsideVolTrajectoryList[TrueTrackOrderList[i].idx_OutVolTrack]->GetTrueTrajectoryCoordinates(s - TrueTrackOrderList[i].s_referece);
    }
  }
  
  return  Dummy_vector;
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetTrueTrajectoryUnitMon(double s)
{
  
  // Function returns the momentum vector direction of the particles trajectory for a value of the "s" parameter and the set of particle initial parameters,
  // - particle type
  // - particle's initial position
  // - particle's initial momentum
  // - magnetic field

  for(int i=0;i<int(SInterVals_true.size());i++) {
    if(s >= SInterVals_true[i](0) && s <= SInterVals_true[i](1)) {
      if(TrueTrackOrderList[i].IsInsideVolTrack) return  TrueInsideVolTrajectoryList[TrueTrackOrderList[i].idx_InVolTrack]->GetTrueTrajectoryUnitMon(s - TrueTrackOrderList[i].s_referece);
      else                                       return  TrueOutsideVolTrajectoryList[TrueTrackOrderList[i].idx_OutVolTrack]->GetTrueTrajectoryUnitMon(s - TrueTrackOrderList[i].s_referece);
    }
  }
  
  return  Dummy_vector;
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetFitTrajectoryCoordinates(double s)
{
  
  for(int i=0;i<int(SInterVals_fit.size());i++) {
    if(s >= SInterVals_fit[i](0) && s <= SInterVals_fit[i](1)) {
      if(FitTrackOrderList[i].IsInsideVolTrack) return  FitInsideVolTrajectoryList[FitTrackOrderList[i].idx_InVolTrack]->GetTrueTrajectoryCoordinates(s - FitTrackOrderList[i].s_referece);
      else                                      return  FitOutsideVolTrajectoryList[FitTrackOrderList[i].idx_OutVolTrack]->GetTrueTrajectoryCoordinates(s - FitTrackOrderList[i].s_referece);
    }
  }
  
  return  Dummy_vector;
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetFitTrajectoryMon(double s)
{

  for(int i=0;i<int(SInterVals_fit.size());i++) {
    if(s >= SInterVals_fit[i](0) && s <= SInterVals_fit[i](1)) {
      if(FitTrackOrderList[i].IsInsideVolTrack) return  FitInsideVolTrajectoryList[FitTrackOrderList[i].idx_InVolTrack]->GetTrueTrajectoryMon(s - FitTrackOrderList[i].s_referece);
      else                                      return  FitOutsideVolTrajectoryList[FitTrackOrderList[i].idx_OutVolTrack]->GetTrueTrajectoryMon(s - FitTrackOrderList[i].s_referece);
    }
  }
  
  return  Dummy_vector;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::SetTrueTrajectories(TVector3 x0, TVector3 p0)
{
  
  //Calculate the trajectory sections inside and outside the B-field volumes

  bool Test = false;
  //Test = true;
  
  TVector3  pos = x0;
  TVector3  mom = p0;
  
  std::vector<int> IdxVolumeList;
  IdxVolumeList.clear();
  for(int ivol=0;ivol<int(VolumeList.size());ivol++) IdxVolumeList.push_back(ivol);
  
  int Outside_Trk_counter = 0;
  
  TrackOrder_t  aTrackOrder;
  aTrackOrder.IsInsideVolTrack = false;
  aTrackOrder.idx_InVolTrack   = -1;
  aTrackOrder.idx_OutVolTrack  = -1;
  
  GTrajectory* aTrajectory = NULL;
  
  TVector3  SRange(0,0,0);
  double s = 0.0;
  
  SInterVals_true.clear();
  TrueTrackOrderList.clear();
  
  bool Intersection = false;
  do {
    aTrajectory = NULL;
    SRange(0) = s;
    
    int Vol_idx  = -1;
    for(int ivol=0;ivol<int(IdxVolumeList.size());ivol++) {
      int tmp_vol = IdxVolumeList[ivol];
      if(VolumeList[tmp_vol]->IsPointInsideGeometry(pos + small_distance*mom.Unit())) Vol_idx = tmp_vol;
    }
    
    if(Vol_idx < 0) {
      //Point outside all volumes
      aTrackOrder.IsInsideVolTrack = false;
      aTrackOrder.idx_InVolTrack   = -1;
      aTrackOrder.idx_OutVolTrack  = Outside_Trk_counter;
      Outside_Trk_counter++;
     
      if(Test) {
        cout << "Position = (" 
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << "Is outside all volumes!!!"
	     << " Outside_Trk_counter = " << Outside_Trk_counter-1 << ", max is = " << TrueOutsideVolTrajectoryList.size()
             << endl;
      }
      
      aTrajectory = TrueOutsideVolTrajectoryList[aTrackOrder.idx_OutVolTrack];
    }
    else {
      //Point inside a volume
      aTrackOrder.IsInsideVolTrack = true;
      aTrackOrder.idx_InVolTrack   = Vol_idx;
      aTrackOrder.idx_OutVolTrack  = -1;
      
      if(Test) {
        cout << "Position = (" 
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << "Is inside volume " << VolumeList[Vol_idx]->GetName().Data()
             << endl;
      }
      
      aTrajectory = TrueInsideVolTrajectoryList[Vol_idx];
    }
    
    aTrajectory->SetInitPositionAndMomentum(pos,mom);
      
    Intersection = false;
    double s_min = +1.0e+6*global->GetUnit("m");
    std::vector<TString>  VolNames;
    VolNames.clear();
    for(int ivol=0;ivol<int(IdxVolumeList.size());ivol++) { //begin loop over volumes
      GGeoObject* Volume = VolumeList[IdxVolumeList[ivol]];
      bool AlreadyUsed = false;
      for(int kkk=0;kkk<int(VolNames.size());kkk++) {
	if(VolNames[kkk] == Volume->GetName()) AlreadyUsed = true;
      }
      if(AlreadyUsed) continue;
      else            VolNames.push_back(Volume->GetName());
      
      if(Test) cout << "Looking at Volume = " << Volume->GetName().Data() << " (" << IdxVolumeList[ivol] << ")" << endl;
      
      int Nboundaries = Volume->GetNBoundarySurfaces();
      for(int iboundary=0;iboundary<Nboundaries;iboundary++) { // being loop over boundaries of geo element
        GSurfaceObject* aBoundary = Volume->GetBoundarySurface(iboundary);
       
	if(Test) cout << "  Looking at Boundary = " << aBoundary->GetName().Data() << endl;
	
        double s = IntersectCoordinates(aBoundary,aTrajectory,global,Dummy_value,true);
	if(Test) cout << "s = " << s/global->GetUnit("cm") << " cm; 1st try" << endl;
	if(s == Dummy_value) continue;
	if(std::isnan(s)) continue;
        if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
        if(s < small_distance) {
	  if(aTrajectory->GetType() == TString("Straight")) continue;
	  double omega   = aTrajectory->GetOmega();
	  s = TMath::Pi()/(TMath::Abs(omega));
	  s = IntersectCoordinates(aBoundary,aTrajectory,global,s,true);
	  if(Test) cout << "s = " << s/global->GetUnit("cm") << " cm; 2nd try" << endl;
	  if(s == Dummy_value) continue;
          if(std::isnan(s)) continue;
	  if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
          if(s < small_distance) continue;
	}

        TVector3 IntersectionPointXYZ = aTrajectory->GetTrueTrajectoryCoordinates(s);
        TVector3 IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);
	
	TVector3 IntersectionMom      = aTrajectory->GetTrueTrajectoryMon(s);
	bool IsminusInside = Volume->IsPointInsideGeometry(IntersectionPointXYZ - small_distance*IntersectionMom.Unit());
	bool IsplusInside  = Volume->IsPointInsideGeometry(IntersectionPointXYZ + small_distance*IntersectionMom.Unit());
	
	bool IsTraversing = false;
	if((IsminusInside && !IsplusInside) || (!IsminusInside && IsplusInside)) IsTraversing = true;

	if(Test) {
	  TString minus("minus outside");
	  TString plus("plus outside");
	  if(IsminusInside) minus = TString("minus inside");
	  if(IsplusInside)  plus  = TString("plus inside");
	  TVector3 Pos_minus = (1.0/global->GetUnit("cm"))*(IntersectionPointXYZ - small_distance*IntersectionMom.Unit());
	  TVector3 Pos_plus  = (1.0/global->GetUnit("cm"))*(IntersectionPointXYZ + small_distance*IntersectionMom.Unit());
	
	  cout << "  s = " << s/global->GetUnit("cm") << " cm, position = ("
	       << IntersectionPointXYZ(0)/global->GetUnit("cm") << ","
	       << IntersectionPointXYZ(1)/global->GetUnit("cm") << ","
	       << IntersectionPointXYZ(2)/global->GetUnit("cm") << ") cm, "
	       << minus.Data() << "  "
	       << "pos- = (" << Pos_minus(0) << "," << Pos_minus(1) << "," << Pos_minus(2) << ") cm; "
	       << plus.Data() << "   "
	       << "pos+ = (" << Pos_plus(0)  << "," << Pos_plus(1)  << "," << Pos_plus(2)  << ") cm; "
	       << endl;
	}
	
        if(aBoundary->IsInMaterial(IntersectionPointUVW) && IsTraversing) {
	  if(s_min > s) s_min = s;
	  Intersection     = true;
	  
	  if(Test) cout << "  Is inside matrial!!!" << endl;
        }
      } // end loop over boundaries of geo element
    } //end loop over volumes
      
    if(Intersection) {
      s        += s_min;
      SRange(1) = s;
      pos = aTrajectory->GetTrueTrajectoryCoordinates(s_min);
      mom = aTrajectory->GetTrueTrajectoryMon(s_min);
      
      if(Test) {
        cout << "Found intersection at for s = " << s/global->GetUnit("cm") << " cm (s_min = " << s_min/global->GetUnit("cm") << " cm) at pos = ("
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << endl;
      }

    }
    else SRange(1) = large_distance;
    
    aTrackOrder.s_referece = SRange(0);
    
    SInterVals_true.push_back(SRange);
    TrueTrackOrderList.push_back(aTrackOrder);
    
    if(Vol_idx >= 0) global->RemoveElementFromList(Vol_idx, IdxVolumeList);
  }
  while(Intersection);
  
  if(Test) {
    cout << endl;
    cout << "S true intervals:" << endl;
    for(int i=0;i<int(SInterVals_true.size());i++) {
      cout << "S-Intervals " << i+1 << " (" << SInterVals_true[i](0)/global->GetUnit("cm") << "," << SInterVals_true[i](1)/global->GetUnit("cm") << ") cm" << endl;
    }
    cout << endl;
  }

  return;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::SetFitTrajectories(TVector3 x0, TVector3 p0)
{

  //Calculate the trajectory sections inside and outside the B-field volumes

  bool Test = false;
  //Test = true;
  
  TVector3  pos = x0;
  TVector3  mom = p0;
  
  std::vector<int> IdxVolumeList_posi;
  std::vector<int> IdxVolumeList_nega;
  IdxVolumeList_posi.clear();
  IdxVolumeList_nega.clear();
  for(int ivol=0;ivol<int(VolumeList.size());ivol++) {
    IdxVolumeList_posi.push_back(ivol);
    IdxVolumeList_nega.push_back(ivol);
  }
  
  int Outside_Trk_counter = 0;
  
  TrackOrder_t  aTrackOrder;
  aTrackOrder.IsInsideVolTrack = false;
  aTrackOrder.idx_InVolTrack   = -1;
  aTrackOrder.idx_OutVolTrack  = -1;
  
  GTrajectory* aTrajectory = NULL;
  
  TVector3  SRange(0,0,0);
  double s = 0.0;
  bool Intersection = false;
  
  SInterVals_fit.clear();
  FitTrackOrderList.clear();
  
  std::vector<TVector3>      SInterVals_posi;
  std::vector<TrackOrder_t>  TrackOrderList_posi;
  SInterVals_posi.clear();
  TrackOrderList_posi.clear();
  
  pos = x0;
  mom = p0;
  s = 0.0;
  Intersection = false;
  do {
    aTrajectory = NULL;
    SRange(0) = s;
    
    int Vol_idx  = -1;
    for(int ivol=0;ivol<int(IdxVolumeList_posi.size());ivol++) {
      int tmp_vol = IdxVolumeList_posi[ivol];
      if(VolumeList[tmp_vol]->IsPointInsideGeometry(pos + small_distance*mom.Unit())) Vol_idx = tmp_vol;
    }
    
    if(Vol_idx < 0) {
      //Point outside all volumes
      aTrackOrder.IsInsideVolTrack = false;
      aTrackOrder.idx_InVolTrack   = -1;
      aTrackOrder.idx_OutVolTrack  = Outside_Trk_counter;
      Outside_Trk_counter++;
      

      if(Test) {
        cout << "Position = (" 
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << "Is outside all volumes!!!"
	     << " Outside_Trk_counter = " << Outside_Trk_counter-1 << ", max is = " << FitOutsideVolTrajectoryList.size()
             << endl;
      }

      aTrajectory = FitOutsideVolTrajectoryList[aTrackOrder.idx_OutVolTrack];
    }
    else {
      //Point inside a volume
      aTrackOrder.IsInsideVolTrack = true;
      aTrackOrder.idx_InVolTrack   = Vol_idx;
      aTrackOrder.idx_OutVolTrack  = -1;
      
      if(Test) {
        cout << "Position = (" 
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << "Is insie volume " << VolumeList[Vol_idx]->GetName().Data()
             << endl;
      }

      aTrajectory = FitInsideVolTrajectoryList[Vol_idx];
    }
    
    aTrajectory->SetInitPositionAndMomentum(pos,mom);
      
    Intersection = false;
    double s_min = +1.0e+6*global->GetUnit("m");
    std::vector<TString>  VolNames;
    VolNames.clear();
    for(int ivol=0;ivol<int(IdxVolumeList_posi.size());ivol++) { //begin loop over volumes
      GGeoObject* Volume = VolumeList[IdxVolumeList_posi[ivol]];
      bool AlreadyUsed = false;
      for(int kkk=0;kkk<int(VolNames.size());kkk++) {
	if(VolNames[kkk] == Volume->GetName()) AlreadyUsed = true;
      }
      if(AlreadyUsed) continue;
      else            VolNames.push_back(Volume->GetName());
      
      if(Test) cout << "Looking at Volume = " << Volume->GetName().Data() << " (" << IdxVolumeList_posi[ivol] << ")" << endl;
      
      int Nboundaries = Volume->GetNBoundarySurfaces();
      for(int iboundary=0;iboundary<Nboundaries;iboundary++) { // being loop over boundaries of geo element
        GSurfaceObject* aBoundary = Volume->GetBoundarySurface(iboundary);
       
        if(Test) cout << "  Looking at Boundary = " << aBoundary->GetName().Data() << endl;
       
        double s = IntersectCoordinates(aBoundary,aTrajectory,global,Dummy_value,true);
	if(s == Dummy_value) continue;
	if(std::isnan(s)) continue;
        if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
	if(s < small_distance) {
	  if(aTrajectory->GetType() == TString("Straight")) continue;
	  double omega   = aTrajectory->GetOmega();
	  s = TMath::Pi()/(TMath::Abs(omega));
	  s = IntersectCoordinates(aBoundary,aTrajectory,global,s,true);
	  if(s == Dummy_value) continue;
          if(std::isnan(s)) continue;
	  if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
          if(s < small_distance) continue;
	}

        TVector3 IntersectionPointXYZ = aTrajectory->GetTrueTrajectoryCoordinates(s);
        TVector3 IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);
	
	TVector3 IntersectionMom      = aTrajectory->GetTrueTrajectoryMon(s);
	bool IsminusInside = Volume->IsPointInsideGeometry(IntersectionPointXYZ - small_distance*IntersectionMom.Unit());
	bool IsplusInside  = Volume->IsPointInsideGeometry(IntersectionPointXYZ + small_distance*IntersectionMom.Unit());

	if(Test) {
	  TString minus("minus outside");
	  TString plus("plus outside");
	  if(IsminusInside) minus = TString("minus inside");
	  if(IsplusInside)  plus  = TString("plus inside");
	
	  cout << "  s = " << s/global->GetUnit("cm") << " cm, position = ("
	       << IntersectionPointXYZ(0)/global->GetUnit("cm") << ","
	       << IntersectionPointXYZ(1)/global->GetUnit("cm") << ","
	       << IntersectionPointXYZ(2)/global->GetUnit("cm") << ") cm, "
	       << minus.Data() << "  " << plus.Data()
	       << endl;
        }
	
	bool IsTraversing = false;
	if((IsminusInside && !IsplusInside) || (!IsminusInside && IsplusInside)) IsTraversing = true;
      
        if(aBoundary->IsInMaterial(IntersectionPointUVW) && IsTraversing) {
	  if(s_min > s) s_min = s;
	  Intersection     = true;
	  
	  if(Test) cout << "  Is inside matrial!!!" << endl;
        }
      } // end loop over boundaries of geo element
    } //end loop over volumes
      
    if(Intersection) {
      s        += s_min;
      SRange(1) = s;
      pos = aTrajectory->GetTrueTrajectoryCoordinates(s_min);
      mom = aTrajectory->GetTrueTrajectoryMon(s_min);
     
      if(Test) {
        cout << "found intersection for s = " << s/global->GetUnit("cm") << " cm ( s_min = " << s_min/global->GetUnit("cm") << ") at ("
             << pos(0)/global->GetUnit("cm") << ","
	     << pos(1)/global->GetUnit("cm") << ","
	     << pos(2)/global->GetUnit("cm") << ") cm"
             << endl;
      }
      
    }
    else SRange(1) = large_distance;
    
    aTrackOrder.s_referece = SRange(0);
    
    SInterVals_posi.push_back(SRange);
    TrackOrderList_posi.push_back(aTrackOrder);

    if(Test) cout << "s-Range = (" << SRange(0)/global->GetUnit("cm") << "," << SRange(1)/global->GetUnit("cm") << ") cm" << endl;
    
    if(Vol_idx >= 0) global->RemoveElementFromList(Vol_idx, IdxVolumeList_posi);
  }
  while(Intersection);


  std::vector<TVector3>      SInterVals_nega;
  std::vector<TrackOrder_t>  TrackOrderList_nega;
  SInterVals_nega.clear();
  TrackOrderList_nega.clear();
  
  pos = x0;
  mom = p0;
  s = 0.0;
  Intersection = false;
  do {
    aTrajectory = NULL;
    SRange(0) = s;
    
    int Vol_idx  = -1;
    for(int ivol=0;ivol<int(IdxVolumeList_nega.size());ivol++) {
      int tmp_vol = IdxVolumeList_nega[ivol];
      if(VolumeList[tmp_vol]->IsPointInsideGeometry(pos - small_distance*mom.Unit())) Vol_idx = tmp_vol;
    }
        
    if(Vol_idx < 0) {
      //Point outside all volumes
      aTrackOrder.IsInsideVolTrack = false;
      aTrackOrder.idx_InVolTrack   = -1;
      aTrackOrder.idx_OutVolTrack  = Outside_Trk_counter;
      Outside_Trk_counter++;

      if(Test) {
        cout << "Position = (" 
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << "Is outside all volumes!!!"
	     << " Outside_Trk_counter = " << Outside_Trk_counter-1 << ", max is = " << FitOutsideVolTrajectoryList.size()
             << endl;
      }
      
      aTrajectory = FitOutsideVolTrajectoryList[aTrackOrder.idx_OutVolTrack];
    }
    else {
      //Point inside a volume
      aTrackOrder.IsInsideVolTrack = true;
      aTrackOrder.idx_InVolTrack   = Vol_idx;
      aTrackOrder.idx_OutVolTrack  = -1;

      if(Test) {
        cout << "Position = (" 
             << pos(0)/global->GetUnit("cm") << "," 
	     << pos(1)/global->GetUnit("cm") << "," 
	     << pos(2)/global->GetUnit("cm") << ") cm, "
	     << " and momentum = ("
	     << mom(0)/global->GetUnit("GeV/c") << "," 
	     << mom(1)/global->GetUnit("GeV/c") << "," 
	     << mom(2)/global->GetUnit("GeV/c") << ") GeV/c, "
	     << "Is insie volume " << VolumeList[Vol_idx]->GetName().Data()
             << endl;
      }
      
      aTrajectory = FitInsideVolTrajectoryList[Vol_idx];
    }
    
    aTrajectory->SetInitPositionAndMomentum(pos,mom);
    
    Intersection = false;
    double s_min = +1.0e+6*global->GetUnit("m");
    std::vector<TString>  VolNames;
    VolNames.clear();
    for(int ivol=0;ivol<int(IdxVolumeList_nega.size());ivol++) { //begin loop over volumes
      GGeoObject* Volume = VolumeList[IdxVolumeList_nega[ivol]];
      bool AlreadyUsed = false;
      for(int kkk=0;kkk<int(VolNames.size());kkk++) {
	if(VolNames[kkk] == Volume->GetName()) AlreadyUsed = true;
      }
      if(AlreadyUsed) continue;
      else            VolNames.push_back(Volume->GetName());
      
      if(Test) cout << "Looking at Volume = " << Volume->GetName().Data() << " (" << IdxVolumeList_nega[ivol] << ")" << endl;
      
      int Nboundaries = Volume->GetNBoundarySurfaces();
      for(int iboundary=0;iboundary<Nboundaries;iboundary++) { // being loop over boundaries of geo element
       GSurfaceObject* aBoundary = Volume->GetBoundarySurface(iboundary);
       
        if(Test) cout << "  Looking at Boundary = " << aBoundary->GetName().Data() << endl;
       
        double s = IntersectCoordinates(aBoundary,aTrajectory,global,Dummy_value,false);
	if(s == Dummy_value) continue;
	if(std::isnan(s)) continue;
	s *= -1;
        if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
	if(s < small_distance) {
	  if(aTrajectory->GetType() == TString("Straight")) continue;
	  double omega   = aTrajectory->GetOmega();
	  s = -TMath::Pi()/(TMath::Abs(omega));
	  s = IntersectCoordinates(aBoundary,aTrajectory,global,s,false);
	  if(s == Dummy_value) continue;
          if(std::isnan(s)) continue;
	  s *= -1;
	  if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
          if(s < small_distance) continue;
	}

        TVector3 IntersectionPointXYZ = aTrajectory->GetTrueTrajectoryCoordinates(-s);
        TVector3 IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);
      
	TVector3 IntersectionMom      = aTrajectory->GetTrueTrajectoryMon(-s);
	bool IsminusInside = Volume->IsPointInsideGeometry(IntersectionPointXYZ - small_distance*IntersectionMom.Unit());
	bool IsplusInside  = Volume->IsPointInsideGeometry(IntersectionPointXYZ + small_distance*IntersectionMom.Unit());

	if(Test) {
	  TString minus("minus outside");
	  TString plus("plus outside");
	  if(IsminusInside) minus = TString("minus inside");
	  if(IsplusInside)  plus  = TString("plus inside");
	  
	  cout << "  s = " << s/global->GetUnit("cm") << " cm, position = ("
	       << IntersectionPointXYZ(0)/global->GetUnit("cm") << ","
	       << IntersectionPointXYZ(1)/global->GetUnit("cm") << ","
	       << IntersectionPointXYZ(2)/global->GetUnit("cm") << ") cm, "
	       << minus.Data() << "  " << plus.Data()
	       << endl;
	}
	
	bool IsTraversing = false;
	if((IsminusInside && !IsplusInside) || (!IsminusInside && IsplusInside)) IsTraversing = true;
	
        if(aBoundary->IsInMaterial(IntersectionPointUVW) && IsTraversing) {
	  if(s_min > s) s_min = s;
	  Intersection     = true;
	  
	  if(Test) cout << "  Is inside matrial!!!" << endl;
        }
      } // end loop over boundaries of geo element
    } //end loop over volumes
    
    if(Intersection) {
      s        += s_min;
      SRange(1) = s;
      pos = aTrajectory->GetTrueTrajectoryCoordinates(-s_min);
      mom = aTrajectory->GetTrueTrajectoryMon(-s_min);
     
      if(Test) {
        cout << "found intersection for s = " << s/global->GetUnit("cm") << " cm ( s_min = " << s_min/global->GetUnit("cm") << ") at ("
             << pos(0)/global->GetUnit("cm") << ","
	     << pos(1)/global->GetUnit("cm") << ","
	     << pos(2)/global->GetUnit("cm") << ") cm"
             << endl;
      }

    }
    else SRange(1) = large_distance;
    
    aTrackOrder.s_referece = -SRange(0);
    
    SInterVals_nega.push_back(SRange);
    TrackOrderList_nega.push_back(aTrackOrder);
    
    if(Test) cout << "s-Range = (" << SRange(0)/global->GetUnit("cm") << "," << SRange(1)/global->GetUnit("cm") << ") cm" << endl;

    if(Vol_idx >= 0) global->RemoveElementFromList(Vol_idx, IdxVolumeList_nega);
  }
  while(Intersection);

  for(int i=0;i<int(SInterVals_nega.size());i++) {
    int idx = SInterVals_nega.size() - 1 - i;
    SRange(0) = -SInterVals_nega[idx](1);
    SRange(1) = -SInterVals_nega[idx](0);
    
    SInterVals_fit.push_back(SRange);
    FitTrackOrderList.push_back(TrackOrderList_nega[idx]);
  }

  for(int i=0;i<int(SInterVals_posi.size());i++) {
    int idx = i;
    SRange(0) = SInterVals_posi[idx](0);
    SRange(1) = SInterVals_posi[idx](1);
    
    SInterVals_fit.push_back(SRange);
    FitTrackOrderList.push_back(TrackOrderList_posi[idx]);
  }
  

  if(Test) {
    cout << endl;
    cout << "S fit intervals:" << endl;
    for(int i=0;i<int(SInterVals_fit.size());i++) {
      cout << "S-Intervals " << i+1 << " (" << SInterVals_fit[i](0)/global->GetUnit("cm") << "," << SInterVals_fit[i](1)/global->GetUnit("cm") << ") cm" << endl;
    }
    cout << endl;
  }

  return;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::SetTrueTrajectories(void)
{

  SetTrueTrajectories(pos0,mom0);

  return;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::SetFitTrajectories(void)
{
  
  z_ref = RefPoint.Dot(zprimeVect);
  s_ref = GetSFromZ(z_ref,0,true);
    
  TVector3 pos = GetTrueTrajectoryCoordinates(s_ref);
  double x_ref = pos.Dot(xprimeVect);
  double y_ref = pos.Dot(yprimeVect);
  
  TVector3 mom = GetTrueTrajectoryMon(s_ref);
  
  double tx = mom.X()/mom.Z();
  double ty = mom.Y()/mom.Z();
  
  //double Pbar = global->GetParticleCharge(Particle)*mom.Mag();
  double Pbar = (mom.Z()/TMath::Abs(mom.Z()))*mom.Mag();
  
  FitParams.push_back(x_ref);
  FitParams.push_back(y_ref);
  FitParams.push_back(tx);
  FitParams.push_back(ty);
  FitParams.push_back(Pbar);
  
  SetFitTrajectories(pos,mom);
  
  return;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::SetAllParameters(std::vector<double> aFitParams)
{
  
  TVector3  pos;
  TVector3  mom;
  
  //pos = GetTrueTrajectoryCoordinates(s_ref);
  //mom = GetTrueTrajectoryMon(s_ref);
  //cout << endl;
  //cout << "Ref pos = (" << pos(0)/global->GetUnit("cm") << "," << pos(1)/global->GetUnit("cm") << "," << pos(2)/global->GetUnit("cm") << ") cm; "
  //     << "Ref mom = (" << mom(0)/global->GetUnit("GeV/c") << "," << mom(1)/global->GetUnit("GeV/c") << "," << mom(2)/global->GetUnit("GeV/c") << ") GeV/c; "
  //     << endl;
  
  FitParams.clear();
  for(int ipar=0;ipar<int(aFitParams.size());ipar++) FitParams.push_back(aFitParams[ipar]);
  
  double x_ref = FitParams[0];
  double y_ref = FitParams[1];
  
  double tx     = FitParams[2];
  double ty     = FitParams[3];
  double p      = TMath::Abs(FitParams[4]);
  double signPz = FitParams[4]/TMath::Abs(FitParams[4]);
  
  double den = sqrt(1 + pow(tx,2) + pow(ty,2));
  double pz  = signPz*p/den;
  double px  = tx*pz;
  double py  = ty*pz;
  
  pos = TVector3(x_ref,y_ref,z_ref);
  mom = TVector3(px,py,pz);
  
  //cout << "New pos = (" << pos(0)/global->GetUnit("cm") << "," << pos(1)/global->GetUnit("cm") << "," << pos(2)/global->GetUnit("cm") << ") cm; "
  //     << "New pos = (" << mom(0)/global->GetUnit("GeV/c") << "," << mom(1)/global->GetUnit("GeV/c") << "," << mom(2)/global->GetUnit("GeV/c") << ") GeV/c; "
  //     << endl;
  //cout << endl;
  
  SetFitTrajectories(pos,mom);
  
  return;
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetSFromZ(double zval,
					    double sinit,
					    bool FromTrueTrajectory)
{
  
  //Calculate the value of the trajectory s parameter for which f_z(s) = z
  
  double limit = 1.0*global->GetDistanceUnit("nm");
  
  const int MaxIterations(MaxIterations_intersection);
  int counter = 0;
  double s0,s1;
  double delta0;
  double Derdelta0;
  double deltaODerdelta;
  
  s1 = 0.0;
  s0 = sinit;
  
  counter = 0;
  double Z    = 0.0;
  double DerZ = 0.0;
  
  if(FromTrueTrajectory) {
    Z    = GetTrueTrajectoryCoordinates(s0).Z();
    DerZ = GetTrueTrajectoryUnitMon(s0).Z();
  }
  else {
    Z    = GetFitTrajectoryCoordinates(s0).Z();
    
    TVector3 mom_tmp = GetFitTrajectoryMon(s0);
    DerZ = mom_tmp.Z()/mom_tmp.Mag();
  }

  delta0    = zval - Z;
  Derdelta0 =      - DerZ;

  if(TMath::Abs(Derdelta0) < 1.0e-10) deltaODerdelta = -1.0;
  else                                deltaODerdelta = delta0/Derdelta0;
  s1 = s0 - deltaODerdelta;
  if(FromTrueTrajectory) delta0 = zval - GetTrueTrajectoryCoordinates(s1).Z();
  else                   delta0 = zval - GetFitTrajectoryCoordinates(s1).Z();

  while(TMath::Abs(delta0) > limit && counter <= MaxIterations) {
    s0 = s1;
    
    if(FromTrueTrajectory) {
      Z    = GetTrueTrajectoryCoordinates(s0).Z();
      DerZ = GetTrueTrajectoryUnitMon(s0).Z();
    }
    else {
      Z    = GetFitTrajectoryCoordinates(s0).Z();
      
      TVector3 mom_tmp = GetFitTrajectoryMon(s0);
      DerZ = mom_tmp.Z()/mom_tmp.Mag();
    }

    delta0    = zval - Z;
    Derdelta0 =      - DerZ;
      
    if(TMath::Abs(Derdelta0) < 1.0e-10) deltaODerdelta = -1.0;
    else                                deltaODerdelta = delta0/Derdelta0;
    
    s1 = s0 - deltaODerdelta;
    if(FromTrueTrajectory) delta0 = zval - GetTrueTrajectoryCoordinates(s1).Z();
    else                   delta0 = zval - GetFitTrajectoryCoordinates(s1).Z();

    counter++;
  }

  if(TMath::Abs(delta0) > limit) {
    if(verbose) {
      cout << "WARNNING in function GTrajectoryMultipleSteps::GetSFromZ, iteration difference " << TMath::Abs(delta0)/global->GetDistanceUnit("cm") 
           << " cm is bigger than limit " << limit/global->GetDistanceUnit("um") << " um after " << counter << " iterations. Result not precise!." 
	   << endl;
    }
    s1 = Dummy_value;
  }

  return s1;
  
}
//====================================================================
void  GTrajectoryMultipleSteps::GetFitParsFromInitConds()
{
  
  //Gets the track parameters from the initial position, momentum and magnetic field
  //In current case of constant non-zero Bfield  => helix-track: 5 parameters
  
  SetTrueTrajectories();
  
  FitParams.clear();
  SetFitTrajectories();
  
  if(verbose) PrintParameters();
  
  return;
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetFitTrackCoordinates(double dummy, double sinit)
{
  
  double s = GetSFromZ(dummy,sinit - s_ref,false);
  
  //cout << "s = " << s/global->GetUnit("cm") << " cm, sinit = " << (sinit - s_ref)/global->GetUnit("cm") << " cm" << endl;
    
  return  GetFitTrajectoryCoordinates(s);
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetFitTrackMomentum(double dummy, double sinit)
{
  
  double s = GetSFromZ(dummy,sinit - s_ref,false);
  
  return  GetFitTrajectoryMon(s);
  
}
//====================================================================
double  GTrajectoryMultipleSteps::GetFitTrackDummyParFromS(double s)
{
  
  //Get the dummy track transport parameter from the s parameter, track initial position and momentum, and track parameters
  //Take a look of the function GetFitTrackCoordinates to understand with the dummy parameter is

  //cout << "Inside GetFitTrackDummyParFromS" << endl;
  
  TVector3 pos_true = GetFitTrajectoryCoordinates(s - s_ref);
  double dummy      = pos_true.Z();
  
  TString units_ttt = TString("cm");
  TVector3 pos_fit  = GetFitTrackCoordinates(dummy,s);
  TVector3 pos_diff = pos_true - pos_fit;
  if(pos_diff.Mag() > 1.0*global->GetUnit("nm")) {
    cout << endl;
    cout << "ERROR in GTrajectoryMultipleSteps::GetFitTrackDummyParFromS:  track position from true function and fit function with parameters s = " << s/global->GetDistanceUnit("cm") 
         << " cm (s_ref = " << s_ref/global->GetDistanceUnit("cm") << " cm, and s - s_ref = " << (s - s_ref)/global->GetDistanceUnit("cm") <<  " cm), "
	 << "and dummy = " << dummy/global->GetUnit(units_ttt) << " " << units_ttt.Data() << " give different positions." << endl;
    cout << "True position = (" << pos_true.X()/global->GetDistanceUnit("cm") << "," << pos_true.Y()/global->GetDistanceUnit("cm") << "," << pos_true.Z()/global->GetDistanceUnit("cm") << ") cm" << endl;
    cout << "Fit  position = (" << pos_fit.X()/global->GetDistanceUnit("cm")  << "," << pos_fit.Y()/global->GetDistanceUnit("cm")  << "," << pos_fit.Z()/global->GetDistanceUnit("cm")  << ") cm" << endl;
    cout << "diff position = (" << pos_diff.X()/global->GetDistanceUnit("cm") << "," << pos_diff.Y()/global->GetDistanceUnit("cm") << "," << pos_diff.Z()/global->GetDistanceUnit("cm") << ") cm" << endl;
    cout << endl;
    assert(false);
  }

  return  dummy;
  
}
//====================================================================
TVector3  GTrajectoryMultipleSteps::GetFitTrackCoorDerWRTDummy(double dummy, double sinit)
{
  
  TVector3 mom = GetFitTrackMomentum(dummy,sinit);
  
  return  TVector3(mom.X()/mom.Z(),mom.Y()/mom.Z(),1.0);
  
}
//====================================================================
void  GTrajectoryMultipleSteps::PrintParameters()
{
  
  cout << endl;
  cout << "Parameters: " << endl;
  cout << " - x0    = " << FitParams[0]/global->GetDistanceUnit("cm")     << " cm"    << endl;
  cout << " - y0    = " << FitParams[1]/global->GetDistanceUnit("cm")     << " cm"    << endl;
  cout << " - tx    = " << FitParams[2]                                   << ""       << endl;
  cout << " - ty    = " << FitParams[3]                                   << ""       << endl;
  cout << " - Pbar  = " << FitParams[4]/global->GetMomentumUnit("GeV/c")  << " GeV/c" << endl;
  cout << endl;
  
  return;
  
}
//====================================================================
