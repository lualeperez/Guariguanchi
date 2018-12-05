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
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TMatrixDSymEigen.h>
#include <TMatrixTSym.h>
#include <TVectorD.h>
#include "include/GTracker.h"

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

//const bool UseMonOrigin = true;
const bool UseMonOrigin = false;

//====================================================================
GTracker::GTracker(TString       aName,
		   GGeometry*    aGeometry,
		   TString       aParticle,
		   TVector3      aPos0,
		   TVector3      aMom0,
		   TVector3      aRefPoint,
		   GGlobalTools* aglobal)
{
  
  verbose = false;
  
  IncludeEloss = false;
  
  Name     = aName;
  Geometry = aGeometry;
  global   = aglobal;
  
  InitializeTrajectory(aParticle,aPos0,aMom0,aRefPoint);
  
  MinMSEffect = 1.0e-2*global->GetUnit("nm");
  
  //Allocating memory for some arrays used in calculations
  F = new TMatrixD*[MaxNLayers];
  K = new TMatrixD*[MaxNLayers];
  for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
    F[nlayer] = new TMatrixD[MaxNLayers];
    K[nlayer] = new TMatrixD[MaxNLayers];
    for(int mlayer=0;mlayer<MaxNLayers;mlayer++) {
      F[nlayer][mlayer].ResizeTo(3,3);
      K[nlayer][mlayer].ResizeTo(3,3);
    }
  }
  
  Cov_deltaR_deltaR = new double***[3];
  for(int i=0;i<3;i++) {
    Cov_deltaR_deltaR[i] = new double**[3];
    for(int j=0;j<3;j++) {
      Cov_deltaR_deltaR[i][j] = new double*[MaxNLayers];
      for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
	Cov_deltaR_deltaR[i][j][nlayer] = new double[MaxNLayers];
      }
    }
  }
  
  Fr_local  = new  TMatrixD[MaxNLayers];
  Fp_local  = new  TMatrixD[MaxNLayers];
  Kr_local  = new  TMatrixD[MaxNLayers];
  Kp_local  = new  TMatrixD[MaxNLayers];
  Fr_global = new  TMatrixD[MaxNLayers];
  Fp_global = new  TMatrixD[MaxNLayers];
  Kr_global = new  TMatrixD[MaxNLayers];
  Kp_global = new  TMatrixD[MaxNLayers];
  for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
    Fr_local[nlayer].ResizeTo(3,3);
    Fp_local[nlayer].ResizeTo(3,3);
    Kr_local[nlayer].ResizeTo(3,3);
    Kp_local[nlayer].ResizeTo(3,3);
    
    Fr_global[nlayer].ResizeTo(3,3);
    Fp_global[nlayer].ResizeTo(3,3);
    Kr_global[nlayer].ResizeTo(3,3);
    Kp_global[nlayer].ResizeTo(3,3);
  }
  
  CovUU = new double*[MaxNhits];
  CovVV = new double*[MaxNhits];
  CovUV = new double*[MaxNhits];
  for(int nhit=0;nhit<MaxNhits;nhit++) { //1st loop on the hit local coordinates
    CovUU[nhit] = new double[MaxNhits];
    CovVV[nhit] = new double[MaxNhits];
    CovUV[nhit] = new double[MaxNhits];
  }
  
  SigmaMS2_local     = new double[MaxNLayers];
  m1Unit_local       = new TVector3[MaxNLayers];
  m2Unit_local       = new TVector3[MaxNLayers];
  SigmaEloss2_local  = new double[MaxNLayers];
  pUnit_local        = new TVector3[MaxNLayers];
  DerUVwrtPos_local  = new double**[MaxNLayers];
  SigmaMS2_global    = new double[MaxNLayers];
  m1Unit_global      = new TVector3[MaxNLayers];
  m2Unit_global      = new TVector3[MaxNLayers];
  SigmaEloss2_global = new double[MaxNLayers];
  pUnit_global       = new TVector3[MaxNLayers];
  DerUVwrtPos_global = new double**[MaxNLayers];
  for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
    DerUVwrtPos_local[nlayer]   = new double*[2];
    DerUVwrtPos_global[nlayer]  = new double*[2];
    for(int iuv=0;iuv<2;iuv++) {
      DerUVwrtPos_local[nlayer][iuv]  = new double[3];
      DerUVwrtPos_global[nlayer][iuv] = new double[3];
    }
  }
  
  Der_Params_prime_wrt_thetaMS_global = new double**[MaxNLayers];
  Der_Params_prime_wrt_Eloss_global   = new double*[MaxNLayers];
  for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
    Der_Params_prime_wrt_thetaMS_global[nlayer] = new double*[MaxNpars];
    Der_Params_prime_wrt_Eloss_global[nlayer]   = new double[MaxNpars];
    for(int ipar=0;ipar<MaxNpars;ipar++) {
      Der_Params_prime_wrt_thetaMS_global[nlayer][ipar] = new double[2];
    }
  }
  
  GammaUDer_local  = new double*[MaxNLayers];
  GammaVDer_local  = new double*[MaxNLayers];
  GammaUDer_global = new double*[MaxNLayers];
  GammaVDer_global = new double*[MaxNLayers];
  for(int nhit=0;nhit<MaxNLayers;nhit++) {
    GammaUDer_local[nhit]  = new double[MaxNpars];
    GammaVDer_local[nhit]  = new double[MaxNpars];
    GammaUDer_global[nhit] = new double[MaxNpars];
    GammaVDer_global[nhit] = new double[MaxNpars];
  }
  
}
//====================================================================
GTracker::~GTracker() 
{
 
  delete Trajectory;
  
  for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
    delete [] F[nlayer];
    delete [] K[nlayer];
  }
  delete [] F;
  delete [] K;
  
  //freeing the allocated memory
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
	delete [] Cov_deltaR_deltaR[i][j][nlayer];
      }
      delete [] Cov_deltaR_deltaR[i][j];
    }
    delete [] Cov_deltaR_deltaR[i];
  }
  delete [] Cov_deltaR_deltaR;

  
  delete [] Fr_local;
  delete [] Fp_local;
  delete [] Kr_local;
  delete [] Kp_local;
  delete [] Fr_global;
  delete [] Fp_global;
  delete [] Kr_global;
  delete [] Kp_global;
  
  for(int nlayer=0;nlayer<MaxNLayers;nlayer++) {
    for(int iuv=0;iuv<2;iuv++) {
      delete [] DerUVwrtPos_local[nlayer][iuv];
      delete [] DerUVwrtPos_global[nlayer][iuv];
    }
    
    for(int ipar=0;ipar<MaxNpars;ipar++) {
      delete [] Der_Params_prime_wrt_thetaMS_global[nlayer][ipar];
    }
    delete [] Der_Params_prime_wrt_thetaMS_global[nlayer];
    delete [] Der_Params_prime_wrt_Eloss_global[nlayer];
    
    delete [] DerUVwrtPos_local[nlayer];
    delete [] DerUVwrtPos_global[nlayer];
    
  }
  delete [] m1Unit_local;
  delete [] m2Unit_local;
  delete [] SigmaMS2_local;
  delete [] pUnit_local;
  delete [] SigmaEloss2_local;
  delete [] DerUVwrtPos_local;
  delete [] m1Unit_global;
  delete [] m2Unit_global;
  delete [] SigmaMS2_global;
  delete [] pUnit_global;
  delete [] SigmaEloss2_global;
  delete [] DerUVwrtPos_global;
  delete [] Der_Params_prime_wrt_thetaMS_global;
  delete [] Der_Params_prime_wrt_Eloss_global;
  
  for(int nhit=0;nhit<MaxNhits;nhit++) {
    delete [] GammaUDer_local[nhit];
    delete [] GammaVDer_local[nhit];
    delete [] GammaUDer_global[nhit];
    delete [] GammaVDer_global[nhit];
  }
  delete [] GammaUDer_local;
  delete [] GammaVDer_local;
  delete [] GammaUDer_global;
  delete [] GammaVDer_global;
  
  for(int nhit=0;nhit<MaxNhits;nhit++) {
    delete [] CovUU[nhit];
    delete [] CovVV[nhit];
    delete [] CovUV[nhit];
  }
  delete [] CovUU;
  delete [] CovVV;
  delete [] CovUV;
  
}
//====================================================================
void  GTracker::InitializeTrajectory(TString   aParticle,
				     TVector3  aPos0,
				     TVector3  aMom0,
				     TVector3  aRefPoint)
{
  
  //Initialization of the Trajectory

  //Get the geometry B-field
  GBField* aGeoBField = Geometry->GetBField();
  
  Trajectory = NULL;  
  if(Geometry->GetBField()->GetType() == TString("Constant")) {
    // Constant B field over all space    
    TVector3 ConstantBfield = (dynamic_cast<GBFieldConstant*>(aGeoBField))->GetConstantBField();

    if(ConstantBfield.Mag() < 1.0e-6*global->GetBfieldUnit("T")) {
      // B-field is zero => straight line trajectory
      TString Name = TString("Straight trajectory for geometry ") + Geometry->GetName();
      Trajectory = new GTrajectoryStraight(Name,
					   aParticle,
					   aPos0,
					   aMom0,
					   aRefPoint,
					   aGeoBField,
					   global);
    }
    else {
      // B-field is constant and different from zero => helix trajectory
      TString Name = TString("Helix trajectory for geometry ") + Geometry->GetName();
      Trajectory = new GTrajectoryHelix(Name,
					aParticle,
					aPos0,
					aMom0,
					aRefPoint,
					aGeoBField,
					global);
    }
    
  }
  else if(Geometry->GetBField()->GetType() == TString("MultipleSteps")) {
    TString Name = TString("Multiple Steps trajectory for geometry ") + Geometry->GetName();
    Trajectory = new GTrajectoryMultipleSteps(Name,
					      aParticle,
					      aPos0,
					      aMom0,
					      aRefPoint,
					      aGeoBField,
					      global);
  }

  if(Trajectory == NULL) {
    cout << endl;
    cout << "ERROR in GTracker Trajectory Initialization. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  return;
  
}
//====================================================================
double   GTracker::GetIntersectCoordinates(GSurfaceObject* aSurface, double sinit,double scut, bool Myverbose)
{
  
  //Calculate the value of the trajectory s parameter for the intersection between the track and a Surface
  
  return    IntersectCoordinates(aSurface,Trajectory,global,sinit,true,scut,Myverbose);
  
}
//====================================================================
void     GTracker::GetAllIntersectCoordinates(GSurfaceObject* aSurface, std::vector<double>& SListIntersections)
{
  
  SListIntersections.clear();
  
  double s = GetIntersectCoordinates(aSurface);
  if(s > Dummy_value) SListIntersections.push_back(s);
  
  return;
  
}
//====================================================================
double   GTracker::GetIntersectionWithWorldVolume()
{
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  //TVector3 HighMomentum = (1.0e+3*global->GetMomentumUnit("GeV/c"))*InitMomentum.Unit();
  
  double s_min = 1.0e+4*global->GetUnit("m");
  GGeoObject* WorldVolume = Geometry->GetWorldVolume();
  
  int Nboundaries = WorldVolume->GetNBoundarySurfaces();
  for(int iboundary=0;iboundary<Nboundaries;iboundary++) {
    GSurfaceObject* aBoundary = WorldVolume->GetBoundarySurface(iboundary);
    
    double si;
  
    //SetTrajectoryInitConditions(InitPosition,HighMomentum);
    //si = GetIntersectCoordinates(aBoundary,Dummy_value);
    
    //SetTrajectoryInitConditions(InitPosition,InitMomentum);
    //si = GetIntersectCoordinates(aBoundary,si);
    si = GetIntersectCoordinates(aBoundary,Dummy_value);
    
    if(std::isnan(si)) continue;
    if(si == Dummy_value) continue;
    if(si < 0.0) continue;
    
    TVector3 IntersectionPointXYZ = Trajectory->GetTrueTrajectoryCoordinates(si);
    TVector3 IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);
    
    if(!aBoundary->IsInMaterial(IntersectionPointUVW)) continue;
    
    if(s_min > si) s_min = si;
  }
  
  if(s_min < 1.0*global->GetDistanceUnit("nm")) s_min = 1.0e+10;
  
  //SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return s_min;
  
}
//====================================================================
void     GTracker::GetIntersectionsWithGeometry(std::vector<IntersectionHit_t>& ItersectionHitList)
{
  
  //Get the intersection points of a track with the objects of a geometry

  bool Myverbose = false;
  //Myverbose = true;
  
  double n_loops = 1.0;
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  //TVector3 HighMomentum = (1.0e+3*global->GetMomentumUnit("GeV/c"))*InitMomentum.Unit();
  
  double s_world = GetIntersectionWithWorldVolume();
  if(verbose || Myverbose) {
    cout << "s_world = " << s_world/global->GetDistanceUnit("cm") << " cm" << endl;
    TVector3 point = Trajectory->GetTrueTrajectoryCoordinates(s_world);
    point *= (1.0/global->GetDistanceUnit("cm"));
    cout << "Intersection point with World volume (x,y,z,r) =  (" << point.X() << "," << point.Y() << "," << point.Z() << "," << sqrt(pow(point.X(),2) + pow(point.Y(),2)) << ") cm" << endl;
  }

  //double omega   = Trajectory->GetOmega();
  double s_limit = Trajectory->GetSLimit(n_loops);
  
  ItersectionHitList.clear();

  if(verbose || Myverbose) {
    cout << endl;
    cout << "particle = " << Trajectory->GetParticle().Data() << ", " 
         << "pos0 = (" << InitPosition.X()/global->GetDistanceUnit("cm") << "," << InitPosition.Y()/global->GetDistanceUnit("cm") << "," << InitPosition.Z()/global->GetDistanceUnit("cm") << ") cm,  "
         << "mom0 = " << InitMomentum.Mag()/global->GetMomentumUnit("GeV/c") << " GeV/c in direction (" << InitMomentum.Unit().X() << "," << InitMomentum.Unit().Y() << "," << InitMomentum.Unit().Z() << "),  "
	 << "s_world = " << s_world/global->GetDistanceUnit("cm") << " cm, "
	 << "s_limit = " << s_limit/global->GetDistanceUnit("cm") << " cm, "
         << endl;
  }

  int Nmax_Intersections = Trajectory->GetMaximumIntersections();

  //double s_prev = -1*global->GetDistanceUnit("cm");
  //double s_prev = Dummy_value;
  for(int k=0;k<Geometry->GetNVoxelesGeoElements();k++) { // begin of loop on voxeled planes
    std::vector<double> BoundaryIntersectons;
    BoundaryIntersectons.clear();
    
    GGeoObject* aGeoElement = Geometry->GetVoxeledGeometryElement(k);
    if(verbose || Myverbose) {
      cout << endl;
      aGeoElement->Print();
      cout << endl;
    }
    
    int Nboundaries = aGeoElement->GetNBoundarySurfaces();
    for(int iboundary=0;iboundary<Nboundaries;iboundary++) { // being loop over boundaries of geo element
      GSurfaceObject* aBoundary = aGeoElement->GetBoundarySurface(iboundary);
      
      bool Intersection = false;
      int  counter_intersection = 0;
      double s;
      
      //Estimation of the initia value of s parameter by supposing a straight trajectory => very large momentum
      //SetTrajectoryInitConditions(InitPosition,HighMomentum);
      //s = GetIntersectCoordinates(aBoundary,s_prev);
      //if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
      //if(s >= 0) s_prev = s;

      //SetTrajectoryInitConditions(InitPosition,InitMomentum);
      //s = GetIntersectCoordinates(aBoundary,s_prev,0.0,(Myverbose || verbose));
      s = GetIntersectCoordinates(aBoundary,Dummy_value,0.0,(Myverbose || verbose));
      if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
      if(std::isnan(s)) continue;
      if(s == Dummy_value) continue;
      if(s > s_world) continue;
      if(s > s_limit) continue;
      //if(s >= 0) s_prev = s;
      //else       continue;

      TVector3 IntersectionPointXYZ = Trajectory->GetTrueTrajectoryCoordinates(s);
      TVector3 IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);
      
      if(verbose || Myverbose) {
	cout << endl;
	cout << "s_intersection (" << counter_intersection+1 << ") = " << s/global->GetDistanceUnit("cm") << " cm at point (x,y,z) = (" 
	     << IntersectionPointXYZ(0)/global->GetDistanceUnit("cm") << ","
	     << IntersectionPointXYZ(1)/global->GetDistanceUnit("cm") << ","
	     << IntersectionPointXYZ(2)/global->GetDistanceUnit("cm") << ") cm, "
	     << "for boundary " << aBoundary->GetName().Data() << ", Type = " << aBoundary->GetType().Data()
	     << endl;
      }
      
      if(aBoundary->IsInMaterial(IntersectionPointUVW)) {
	BoundaryIntersectons.push_back(s);
	Intersection = true;
	counter_intersection++;
	
	if(verbose || Myverbose) {
	  cout << "Intersection (" << counter_intersection << ") is within boundaries!!!" << endl;
        }
	
      }
      
      while(Intersection && counter_intersection < Nmax_Intersections) {
	if(verbose || Myverbose) {
	  cout << "Looking for intersection " << counter_intersection+1 << " on boundary surfave " << aBoundary->GetName().Data() << " of volume " << aGeoElement->GetName().Data() << endl;
	}
	Intersection = false;
	double scut   = s;
	double s_init = Trajectory->GetInitialSFor2ndIntersection(s);
	s = GetIntersectCoordinates(aBoundary,s_init,scut,(Myverbose || verbose));
        if(TMath::Abs(s) < 2.0*global->GetUnit("nm")) s = 0.0;
        if(std::isnan(s)) continue;
        if(s == Dummy_value) continue;
        if(s > s_world) continue;
        if(s > s_limit) continue;

        IntersectionPointXYZ = Trajectory->GetTrueTrajectoryCoordinates(s);
        IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);
      
        if(verbose || Myverbose) {
	  cout << endl;
	  cout << "s_intersection (" << counter_intersection+1 << ") = " << s/global->GetDistanceUnit("cm") << " cm at point (x,y,z) = (" 
	       << IntersectionPointXYZ(0)/global->GetDistanceUnit("cm") << ","
	       << IntersectionPointXYZ(1)/global->GetDistanceUnit("cm") << ","
	       << IntersectionPointXYZ(2)/global->GetDistanceUnit("cm") << ") cm, "
	       << "for boundary " << aBoundary->GetName().Data() << ", Type = " << aBoundary->GetType().Data()
	       << endl;
        }
      
        if(aBoundary->IsInMaterial(IntersectionPointUVW)) {
	  bool NotInside = true;
	  for(int kkk=0;kkk<int(BoundaryIntersectons.size());kkk++) {
	    if(BoundaryIntersectons[kkk] == s) {
	      NotInside = false;
	      break;
	    }
	  }
	  if(NotInside) {
	    BoundaryIntersectons.push_back(s);
	    Intersection = true;
	    counter_intersection++;
	
	    if(verbose || Myverbose) {
	     cout << "Intersection (" << counter_intersection << ") is within boundaries!!!" << endl;
            }
	  }
	
        }

      }
      
      if(verbose || Myverbose) cout << endl;

    } // end loop over boundaries of geo element
    
    if(BoundaryIntersectons.size() == 0) continue;

    global->OrderListList(BoundaryIntersectons); //Order the intersections

    std::vector<GeoElementInOut_t> SIntersectionWithBoundaries;
    SIntersectionWithBoundaries.clear();

    for(int inter=0;inter<int(BoundaryIntersectons.size() - 1);inter++) {
      double s1 = BoundaryIntersectons[inter];
      double s2 = BoundaryIntersectons[inter+1];
      double s  = 0.5*(s1 + s2);
      
      TVector3 PosXYZ;
      
      PosXYZ         = Trajectory->GetTrueTrajectoryCoordinates(s1 + 5.0*global->GetDistanceUnit("nm"));
      bool In1       = aGeoElement->IsPointInsideGeometry(PosXYZ);
      
      PosXYZ         = Trajectory->GetTrueTrajectoryCoordinates(s2 + 5.0*global->GetDistanceUnit("nm"));
      bool In2       = aGeoElement->IsPointInsideGeometry(PosXYZ);
      
      PosXYZ         = Trajectory->GetTrueTrajectoryCoordinates(s);
      bool In_middle = aGeoElement->IsPointInsideGeometry(PosXYZ);
      
      if(In1 && !In2 && In_middle) {
	GeoElementInOut_t AGeoElementInOut;
	AGeoElementInOut.s_in  = s1;
	AGeoElementInOut.s_out = s2;
	SIntersectionWithBoundaries.push_back(AGeoElementInOut);
      }
    }

    
    for(int imult=0;imult<int(SIntersectionWithBoundaries.size());imult++) {  //begin of loop over multiple intersections
      bool SensPoint             = false;
      bool IntersectsMainSurface = true;
      
      GSurfaceObject* MainSurface = aGeoElement->GetMainSurface();
      
      double s_intersection0 = 0.5*(SIntersectionWithBoundaries[imult].s_in + SIntersectionWithBoundaries[imult].s_out);
      double s_intersection  = GetIntersectCoordinates(MainSurface,s_intersection0);
      if(s_intersection < 0.0) {
	IntersectsMainSurface = false;
	s_intersection        = s_intersection0;
      }
      if(s_intersection > s_world) continue;
      if(s_intersection > s_limit) continue;
      
      TVector3 IntersectionPointXYZ = Trajectory->GetTrueTrajectoryCoordinates(s_intersection);
      TVector3 IntersectionPointUVW = MainSurface->GetUVWFromXYZ(IntersectionPointXYZ);
      
      if(verbose) {
        cout << "Intersection of particle with "
             << "geoElement = " << aGeoElement->GetIndex() << ",  "
	     << "GeoElement Name = " << aGeoElement->GetName().Data() << ", "
	     << "s = " << s_intersection/global->GetDistanceUnit("mm") << " mm, "
	     << "Itersection(X,Y,Z) = (" 
	     << IntersectionPointXYZ.X()/global->GetDistanceUnit("cm") << "," 
	     << IntersectionPointXYZ.Y()/global->GetDistanceUnit("cm") << "," 
	     << IntersectionPointXYZ.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	     << "Itersection(U,V,W) = (" 
	     << IntersectionPointUVW.X()/global->GetDistanceUnit("cm") << "," 
	     << IntersectionPointUVW.Y()/global->GetDistanceUnit("cm") << "," 
	     << IntersectionPointUVW.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	     << endl;
      }

      if(!MainSurface->IsInMaterial(IntersectionPointUVW)) {
        if(verbose) {
	  cout << "Intersection of particle with "
	       << "geoElement = " << aGeoElement->GetIndex() << ",  "
	       << "GeoElement Name = " << aGeoElement->GetName().Data() << ", "
	       << "s = " << s_intersection/global->GetDistanceUnit("mm") << " mm, "
	       << "Itersection(X,Y,Z) = (" 
	       << IntersectionPointXYZ.X()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointXYZ.Y()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointXYZ.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	       << "Itersection(U,V,W) = (" 
	       << IntersectionPointUVW.X()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointUVW.Y()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointUVW.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	       << "Not in material!!!"
	       << endl;
        }
        continue;
      }
      else {
        if(verbose) {
	  cout << "Intersection of particle with "
	       << "geoElement = " << aGeoElement->GetIndex() << ",  "
	       << "GeoElement Name = " << aGeoElement->GetName().Data() << ", "
	       << "s = " << s_intersection/global->GetDistanceUnit("mm") << " mm, "
	       << "Itersection(X,Y,Z) = (" 
	       << IntersectionPointXYZ.X()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointXYZ.Y()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointXYZ.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	       << "Itersection(U,V,W) = (" 
	       << IntersectionPointUVW.X()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointUVW.Y()/global->GetDistanceUnit("cm") << "," 
	       << IntersectionPointUVW.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	       << "In material!!!"
	       << endl;
        }
      }

      if(aGeoElement->GetIsSensitive() && IntersectsMainSurface) {
        if(MainSurface->IsInSensitiveMaterial(IntersectionPointUVW)) {
	  SensPoint = true;
	  if(verbose) {
	    cout << "geoElement = " << aGeoElement->GetIndex() << ",  "
	         << "GeoElement Name = " << aGeoElement->GetName().Data() << ", "
                 << "s = " << s_intersection/global->GetDistanceUnit("mm") << " mm, "
	         << "Itersection(X,Y,Z) = (" 
		 << IntersectionPointXYZ.X()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointXYZ.Y()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointXYZ.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	         << "Itersection(U,V,W) = (" 
		 << IntersectionPointUVW.X()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointUVW.Y()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointUVW.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	         << "In sensitive region of material!!!"
                 << endl;
	  }
        }
        else {
	  if(verbose) {
	    cout << "geoElement = " << aGeoElement->GetIndex() << ",  "
	         << "GeoElement Name = " << aGeoElement->GetName().Data() << ", "
                 << "s = " << s_intersection/global->GetDistanceUnit("mm") << " mm, "
	         << "Itersection(X,Y,Z) = (" 
		 << IntersectionPointXYZ.X()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointXYZ.Y()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointXYZ.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	         << "Itersection(U,V,W) = (" 
		 << IntersectionPointUVW.X()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointUVW.Y()/global->GetDistanceUnit("cm") << "," 
		 << IntersectionPointUVW.Z()/global->GetDistanceUnit("cm") << ") cm,  "
	         << "In non-sensitive region of material!!!"
                 << endl;
	  }
        }
      }

      IntersectionHit_t AHit;
      AHit.s                = s_intersection;
      AHit.geoElement_idx   = aGeoElement->GetIndex();
      AHit.IsSensitivePoint = SensPoint;
      AHit.s_in             = SIntersectionWithBoundaries[imult].s_in;
      AHit.s_out            = SIntersectionWithBoundaries[imult].s_out;
      ItersectionHitList.push_back(AHit);
    } //end of loop over multiple intersections

  } // end of loop on voxeled planes
  if(verbose) cout << endl;

  //Order hits wrt s parameter from low to high
  global->OrderIntersectionHitList(ItersectionHitList);
  
  //Cutting away hits on a given subsystem according to table of cuts on number of hits
  Geometry->ApplyTrackCuts(InitMomentum,ItersectionHitList);
  //Geometry->ApplyTrackCuts2(InitMomentum,ItersectionHitList);

  if(verbose) {
    for(int k=0;k<int(ItersectionHitList.size());k++) {
      TVector3 Pos = Trajectory->GetTrueTrajectoryCoordinates(ItersectionHitList[k].s);
      Pos *= 1.0/global->GetDistanceUnit("cm");
      cout << k+1 << "  ";
      cout << "s = " << ItersectionHitList[k].s/global->GetDistanceUnit("cm") << " cm,  ";
      cout << "geoElementID = " << ItersectionHitList[k].geoElement_idx+1 << " (" << Geometry->GetGeometryElement(ItersectionHitList[k].geoElement_idx)->GetName().Data() << "),  ";
      cout << "Intersection position = (" << Pos.X() << "," << Pos.Y() << "," << Pos.Z() << ") cm,  ";
      if(ItersectionHitList[k].IsSensitivePoint) cout << "Sensitive = yes";
      else                                       cout << "Sensitive = no";
      cout << endl;
    }
  }

  //SetTrajectoryInitConditions(InitPosition,InitMomentum);

  return;
  
}
//====================================================================
double   GTracker::GetIntersectionsWithGeoElement(double s0, GGeoObject* aGeoElement)
{

  double s_ave = 100.0*global->GetDistanceUnit("m");
#if 1
  std::vector<double> BoundaryIntersectons;
  BoundaryIntersectons.clear();
  int Nboundaries = aGeoElement->GetNBoundarySurfaces();
  for(int iboundary=0;iboundary<Nboundaries;iboundary++) { // being loop over boundaries of geo element
    double s;
    GSurfaceObject* aBoundary = aGeoElement->GetBoundarySurface(iboundary);
    s = GetIntersectCoordinates(aBoundary,s0);
    if(s < 0) continue;

    TVector3 IntersectionPointXYZ = Trajectory->GetTrueTrajectoryCoordinates(s);
    TVector3 IntersectionPointUVW = aBoundary->GetUVWFromXYZ(IntersectionPointXYZ);

    if(aBoundary->IsInMaterial(IntersectionPointUVW)) BoundaryIntersectons.push_back(s);
  } // end loop over boundaries of geo element
  global->OrderListList(BoundaryIntersectons); //Order the intersections

  if(BoundaryIntersectons.size() == 2) {
    s_ave = 0.5*(BoundaryIntersectons[0] + BoundaryIntersectons[1]);
  }
#endif

  //GSurfaceObject* aBoundary = aGeoElement->GetMainSurface();
  //s_ave = GetIntersectCoordinates(aBoundary,s0);
  
  return  s_ave;
  
}
//====================================================================
double   GTracker::GetIntersectionXOX0(IntersectionHit_t AHit, GGeoObject* aGeoElement)
{
  
  double XOX0  = (TMath::Abs(AHit.s_in - AHit.s_out)/aGeoElement->GetThickness());
  XOX0        *= aGeoElement->GetXOX0();
  
  return XOX0;
  
}
//====================================================================
void  GTracker::GetGeometrySystemsMaterialBudget(std::vector<double>& SystemMatBudget, double& TotMatBudget)
{
  
  //This function calculates the material budget, both the total and per system, that a particle of given momentum encounters
  
  //Obtain first the intersections of the track with the different materials
  std::vector<IntersectionHit_t> ItersectionHitList;
  ItersectionHitList.clear();
  GetIntersectionsWithGeometry(ItersectionHitList);
  
  TotMatBudget = 0.0;
  
  SystemMatBudget.clear();
  for(int isys=0;isys<Geometry->GetNSystemNames();isys++)  SystemMatBudget.push_back(0.0);
  
  int ihit_sens_max = -999;
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ItersectionHitList[ihit].IsSensitivePoint) {
      if(ihit_sens_max < ihit) ihit_sens_max = ihit;
    }
  }
  
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ihit > ihit_sens_max) continue;
    double s        = ItersectionHitList[ihit].s;
    int igeoElement = ItersectionHitList[ihit].geoElement_idx;
    
    GGeoObject* GeoElement = Geometry->GetGeometryElement(igeoElement);
    double XOX0 = GetIntersectionXOX0(ItersectionHitList[ihit],GeoElement);
    
    if(std::isnan(XOX0)) {
      cout << "X/X0 is nan for " 
           << "s = " << s/global->GetDistanceUnit("cm") << " cm, "
	   << "geometry element  " << igeoElement << " with name " << GeoElement->GetName().Data() << ","
	   << "plane material " << GeoElement->GetMaterial().Data() << ", "
	   << "X/X0 for ref thickness = " << GeoElement->GetXOX0()
           << endl;
    }
    
    TotMatBudget += XOX0;
    for(int isys=0;isys<Geometry->GetNSystemNames();isys++) {
      if(Geometry->GetGeometryElement(igeoElement)->GetSystemName() == Geometry->GetASystemName(isys)) SystemMatBudget[isys] += XOX0;
    }
  }
  
  return;
  
}
//====================================================================
void     GTracker::GetMaterialBudget(double& MatBudget, int& Nhits, double& dist1stPointToRef)
{
  
  //This function calculates the material budget that a particle of given momentum encounters for 
  //a geometry.
  //The function also calculates the number of "hits" that the particle produces, i.e., the number of intersections 
  //with sensitive material
  
  //Obtain first the intersections of the track with the different materials
  std::vector<IntersectionHit_t> ItersectionHitList;
  ItersectionHitList.clear();
  GetIntersectionsWithGeometry(ItersectionHitList);
  
  GetMaterialBudget_FromPlaneList(ItersectionHitList,MatBudget,Nhits,dist1stPointToRef);
  
  return;
  
}
//====================================================================
void  GTracker::GetMaterialBudget_FromPlaneList(std::vector<IntersectionHit_t> ItersectionHitList,
						double& MatBudget, int& Nhits, double& dist1stPointToRef)
{
  
  MatBudget          = 0.0;
  Nhits              = 0;
  dist1stPointToRef  = 0.0;
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  
  int highest_sens_hit = -1;
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ItersectionHitList[ihit].IsSensitivePoint) {
      Nhits++;
      if(Nhits == 1) dist1stPointToRef = (Trajectory->GetTrueTrajectoryCoordinates(ItersectionHitList[ihit].s) - InitPosition).Mag();
      if(highest_sens_hit < ihit) highest_sens_hit = ihit;
    }
  }
  
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ihit > highest_sens_hit) continue;
    
    double s        = ItersectionHitList[ihit].s;
    int igeoElement = ItersectionHitList[ihit].geoElement_idx;
    
    GGeoObject* GeoElement = Geometry->GetGeometryElement(igeoElement);
    double XOX0 = GetIntersectionXOX0(ItersectionHitList[ihit],GeoElement);
    
    if(std::isnan(XOX0)) {
      cout << "X/X0 is nan for " 
           << "s = " << s/global->GetDistanceUnit("cm") << " cm, "
	   << "geometry element  " << igeoElement << " with name " << GeoElement->GetName().Data() << ","
	   << "plane material " << GeoElement->GetMaterial().Data() << ", "
	   << "X/X0 for ref thickness = " << GeoElement->GetXOX0()
           << endl;
    }
    
    MatBudget += XOX0;
  }
  
  return;
  
}
//====================================================================
TVector3   GTracker::GetTrueTrackInterCoordDerWRTPos_Numerical(double s0, GGeoObject* aGeoElement, int der_index)
{
  
  //Return the derivative of the track position w.r.t. initial position with numerical method
  
  if(der_index <0 || der_index > 2) {
    cout << endl;
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTPos_Numerical:: ";
    if(der_index  > 2) cout << "derivative index with value " << der_index  << " is ouside range (0,2). Exting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }

  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
    
  TVector3 InitPosition_p;
  TVector3 InitPosition_m;
  TVector3 InitPosition_p2;
  TVector3 InitPosition_m2;  
  
  double h = 1.0e-2*global->GetDistanceUnit("mm");
  TVector3 Prev_der(-999.9,-999.9,-999.9);
  TVector3 Current_der(-999.9,-999.9,-999.9);
  TVector3 Current_der1,Current_der2;
  
  int counter = 0;
  TVector3 diff(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;

    if(counter == 1 && DoRichardson) {
      InitPosition_p = InitPosition;
      InitPosition_m = InitPosition;
      InitPosition_p(der_index) += 0.5*h1;
      InitPosition_m(der_index) -= 0.5*h1;
    }
    
    InitPosition_p2 = InitPosition;
    InitPosition_m2 = InitPosition;
    InitPosition_p2(der_index) += 0.5*h2;
    InitPosition_m2(der_index) -= 0.5*h2;
    
    if(DoRichardson) {
      if(counter == 1) {
	SetTrajectoryInitConditions(InitPosition_p,InitMomentum);
	double s_tmp_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1   = Trajectory->GetTrueTrajectoryCoordinates(s_tmp_p);
	
	SetTrajectoryInitConditions(InitPosition_m,InitMomentum);
	double s_tmp_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
        Current_der1  -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp_m);
	
        Current_der1  *= (1.0/h1);
      }
      else Current_der1 = Prev_der;
    }
    
    SetTrajectoryInitConditions(InitPosition_p2,InitMomentum);
    double s_tmp2_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2    = Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_p);
    
    SetTrajectoryInitConditions(InitPosition_m2,InitMomentum);
    double s_tmp2_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2   -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_m);
    
    Current_der2   *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der = (1.0/3.0)*(4.0*Current_der2 - Current_der1);
      if(counter == 1) Prev_der = Current_der2;
    }
    else Current_der = Current_der2;
    
    diff = Current_der - Prev_der;
    if(Current_der.Mag() > 1.0e-8) diff *= (1.0/Current_der.Mag());
    
    h /= 2.0;
    Prev_der = Current_der;
  }
  while(counter <= Nmax_iterations && diff.Mag() > epsilon_derivatives);
  
  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTPos_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << diff.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return Current_der;
  
}
//====================================================================
TVector3   GTracker::GetTrueTrackInterCoordDerWRTMom_Numerical(double s0, GGeoObject* aGeoElement, int der_index)
{
  
  //Return the derivative of the track position w.r.t. initial momentum with numerical method
  
  if(der_index <0 || der_index > 2) {
    cout << endl;
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTMom_Numerical:: ";
    if(der_index  > 2) cout << "derivative index with value " << der_index  << " is ouside range (0,2). Exting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }

  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  TVector3 InitMomentum_p;
  TVector3 InitMomentum_m;
  TVector3 InitMomentum_p2;
  TVector3 InitMomentum_m2;
  
  //double h = 1.0e-2*global->GetMomentumUnit("MeV/c");
  double h = 1.0e-3*InitMomentum.Mag();
  TVector3 Prev_der(-999.9,-999.9,-999.9);
  TVector3 Current_der(-999.9,-999.9,-999.9);
  TVector3 Current_der1,Current_der2;
  
  int counter = 0;
  TVector3 diff(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;
    
    if(counter == 1 && DoRichardson) {
      InitMomentum_p  = InitMomentum;
      InitMomentum_m  = InitMomentum;
      InitMomentum_p(der_index) += 0.5*h1;
      InitMomentum_m(der_index) -= 0.5*h1;
    }
    
    InitMomentum_p2  = InitMomentum;
    InitMomentum_m2  = InitMomentum;
    InitMomentum_p2(der_index) += 0.5*h2;
    InitMomentum_m2(der_index) -= 0.5*h2;
    
    if(DoRichardson) {
      if(counter == 1) {
	SetTrajectoryInitConditions(InitPosition,InitMomentum_p);
	double s_tmp_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1   = Trajectory->GetTrueTrajectoryCoordinates(s_tmp_p);
	
	SetTrajectoryInitConditions(InitPosition,InitMomentum_m);
	double s_tmp_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
        Current_der1  -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp_m);
	
        Current_der1  *= (1.0/h1);
      }
      else Current_der1 = Prev_der;
    }
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum_p2);
    double s_tmp2_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2    = Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_p);
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum_m2);
    double s_tmp2_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2   -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_m);
    
    Current_der2   *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der = (1.0/3.0)*(4.0*Current_der2 - Current_der1);
      if(counter == 1) Prev_der = Current_der2;
    }
    else Current_der = Current_der2;
    
    diff = Current_der - Prev_der;
    if(Current_der.Mag() > 1.0e-8) diff *= (1.0/Current_der.Mag());
    
    h /= 2.0;
    Prev_der = Current_der;
  }
  while(counter <= Nmax_iterations && diff.Mag() > epsilon_derivatives);
  
  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTMom_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << diff.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return  Current_der;
  
}
//====================================================================
TVector3   GTracker::GetTrueTrackInterMomDerWRTPos_Numerical(double s0, GGeoObject* aGeoElement, int der_index)
{
  
  
  
  //Return the derivative of the track momentum w.r.t. initial position with numerical method
  
  if(der_index <0 || der_index > 2) {
    cout << endl;
    cout << "WARNNIN inside GTracker::GetTrueTrackInterMomDerWRTPos_Numerical:: ";
    if(der_index  > 2) cout << "derivative index with value " << der_index  << " is ouside range (0,2). Exting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  TVector3 InitPosition_p;
  TVector3 InitPosition_m;
  TVector3 InitPosition_p2;
  TVector3 InitPosition_m2;
  
  double h = 1.0e-2*global->GetDistanceUnit("mm");
  TVector3 Prev_der(-999.9,-999.9,-999.9);
  TVector3 Current_der(-999.9,-999.9,-999.9);
  TVector3 Current_der1,Current_der2;
  
  int counter = 0;
  TVector3 diff(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;
    
    if(counter == 1 && DoRichardson) {
      InitPosition_p  = InitPosition;
      InitPosition_m  = InitPosition;
      InitPosition_p(der_index) += 0.5*h1;
      InitPosition_m(der_index) -= 0.5*h1;
    }
    InitPosition_p2  = InitPosition;
    InitPosition_m2  = InitPosition;
    InitPosition_p2(der_index) += 0.5*h2;
    InitPosition_m2(der_index) -= 0.5*h2;

    if(DoRichardson) {
      if(counter == 1) {
	SetTrajectoryInitConditions(InitPosition_p,InitMomentum);
	double s_tmp_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1   = Trajectory->GetTrueTrajectoryMon(s_tmp_p);
	
	SetTrajectoryInitConditions(InitPosition_m,InitMomentum);
	double s_tmp_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1  -= Trajectory->GetTrueTrajectoryMon(s_tmp_m);
	
        Current_der1  *= (1.0/h1);
      }
      else Current_der1 = Prev_der;
    }
    
    SetTrajectoryInitConditions(InitPosition_p2,InitMomentum);
    double s_tmp2_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2    = Trajectory->GetTrueTrajectoryMon(s_tmp2_p);
    
    SetTrajectoryInitConditions(InitPosition_m2,InitMomentum);
    double s_tmp2_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2   -= Trajectory->GetTrueTrajectoryMon(s_tmp2_m);
    
    Current_der2   *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der = (1.0/3.0)*(4.0*Current_der2 - Current_der1);
      if(counter == 1) Prev_der = Current_der2;
    }
    else Current_der = Current_der2;
    
    diff = Current_der - Prev_der;
    if(Current_der.Mag() > 1.0e-8) diff *= (1.0/Current_der.Mag());
    
    h /= 2.0;
    Prev_der = Current_der;
  }
  while(counter <= Nmax_iterations && diff.Mag() > epsilon_derivatives);

  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GTracker::GetTrueTrackInterMomDerWRTPos_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << diff.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return  Current_der;
  
}
//====================================================================
TVector3   GTracker::GetTrueTrackInterMomDerWRTMom_Numerical(double s0, GGeoObject* aGeoElement, int der_index)
{
  
  //Return the derivative of the track momentum w.r.t. initial momentum with numerical method
  
  if(der_index <0 || der_index > 2) {
    cout << endl;
    cout << "WARNNIN inside GTracker::GetTrueTrackInterMomDerWRTMom_Numerical:: ";
    if(der_index  > 2) cout << "derivative index with value " << der_index  << " is ouside range (0,2). Exting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  TVector3 InitMomentum_p;
  TVector3 InitMomentum_m;
  TVector3 InitMomentum_p2;
  TVector3 InitMomentum_m2;
  
  //double h = 1.0e-2*global->GetMomentumUnit("MeV/c");
  double h = 1.0e-3*InitMomentum.Mag();
  TVector3 Prev_der(-999.9,-999.9,-999.9);
  TVector3 Current_der(-999.9,-999.9,-999.9);
  TVector3 Current_der1,Current_der2;
  
  int counter = 0;
  TVector3 diff(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;

  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;
    
    if(counter == 1 && DoRichardson) {
      InitMomentum_p  = InitMomentum;
      InitMomentum_m  = InitMomentum;
      InitMomentum_p(der_index) += 0.5*h1;
      InitMomentum_m(der_index) -= 0.5*h1;
    }
    InitMomentum_p2  = InitMomentum;
    InitMomentum_m2  = InitMomentum;
    InitMomentum_p2(der_index) += 0.5*h2;
    InitMomentum_m2(der_index) -= 0.5*h2;

    if(DoRichardson) {
      if(counter == 1) {
	SetTrajectoryInitConditions(InitPosition,InitMomentum_p);
	double s_tmp_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1   = Trajectory->GetTrueTrajectoryMon(s_tmp_p);
	
	SetTrajectoryInitConditions(InitPosition,InitMomentum_m);
	double s_tmp_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1  -= Trajectory->GetTrueTrajectoryMon(s_tmp_m);
	
        Current_der1  *= (1.0/h1);
      }
      else Current_der1 = Prev_der;
    }
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum_p2);
    double s_tmp2_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2    = Trajectory->GetTrueTrajectoryMon(s_tmp2_p);
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum_m2);
    double s_tmp2_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2   -= Trajectory->GetTrueTrajectoryMon(s_tmp2_m);
    
    Current_der2   *= (1.0/h2);
 
    if(DoRichardson) {
      Current_der = (1.0/3.0)*(4.0*Current_der2 - Current_der1);
      if(counter == 1) Prev_der = Current_der2;
    }
    else Current_der = Current_der2;
    
    diff = Current_der - Prev_der;
    if(Current_der.Mag() > 1.0e-8) diff *= (1.0/Current_der.Mag());
    
    h /= 2.0;
    Prev_der = Current_der;
  }
  while(counter <= Nmax_iterations && diff.Mag() > epsilon_derivatives);
  
  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GTracker::GetTrueTrackInterMomDerWRTMom_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << diff.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return  Current_der;
  
}
//====================================================================
void      GTracker::GetTrueTrackInterCoordAndMomDerWRTPos_Numerical(double s0, GGeoObject* aGeoElement, int der_index, TVector3& DerCoor, TVector3& DerMom)
{
  
  //Return the derivative of the track position and track momentum w.r.t. initial position with numerical method
  
  if(der_index <0 || der_index > 2) {
    cout << endl;
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTPos_Numerical:: ";
    if(der_index  > 2) cout << "derivative index with value " << der_index  << " is ouside range (0,2). Exting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }

  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
    
  TVector3 InitPosition_p;
  TVector3 InitPosition_m;
  TVector3 InitPosition_p2;
  TVector3 InitPosition_m2;  
  
  double h = 1.0e-2*global->GetDistanceUnit("mm");
  
  TVector3 Prev_der_pos(-999.9,-999.9,-999.9);
  TVector3 Current_der_pos(-999.9,-999.9,-999.9);
  TVector3 Current_der1_pos,Current_der2_pos;
  
  TVector3 Prev_der_mom(-999.9,-999.9,-999.9);
  TVector3 Current_der_mom(-999.9,-999.9,-999.9);
  TVector3 Current_der1_mom,Current_der2_mom;
  
  int counter = 0;
  TVector3 diff_pos(0.0,0.0,0.0);
  TVector3 diff_mom(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;

    if(counter == 1 && DoRichardson) {
      InitPosition_p = InitPosition;
      InitPosition_m = InitPosition;
      InitPosition_p(der_index) += 0.5*h1;
      InitPosition_m(der_index) -= 0.5*h1;
    }
    
    InitPosition_p2 = InitPosition;
    InitPosition_m2 = InitPosition;
    InitPosition_p2(der_index) += 0.5*h2;
    InitPosition_m2(der_index) -= 0.5*h2;
    
    if(DoRichardson) {
      if(counter == 1) {
	SetTrajectoryInitConditions(InitPosition_p,InitMomentum);
	double s_tmp_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1_pos   = Trajectory->GetTrueTrajectoryCoordinates(s_tmp_p);
	Current_der1_mom   = Trajectory->GetTrueTrajectoryMon(s_tmp_p);
	
	SetTrajectoryInitConditions(InitPosition_m,InitMomentum);
	double s_tmp_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
        Current_der1_pos  -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp_m);
	Current_der1_mom  -= Trajectory->GetTrueTrajectoryMon(s_tmp_m);
	
        Current_der1_pos  *= (1.0/h1);
	Current_der1_mom  *= (1.0/h1);
      }
      else {
	Current_der1_pos = Prev_der_pos;
	Current_der1_mom = Prev_der_mom;
      }
    }
    
    SetTrajectoryInitConditions(InitPosition_p2,InitMomentum);
    double s_tmp2_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2_pos    = Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_p);
    Current_der2_mom    = Trajectory->GetTrueTrajectoryMon(s_tmp2_p);
    
    SetTrajectoryInitConditions(InitPosition_m2,InitMomentum);
    double s_tmp2_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2_pos   -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_m);
    Current_der2_mom   -= Trajectory->GetTrueTrajectoryMon(s_tmp2_m);
    
    Current_der2_pos   *= (1.0/h2);
    Current_der2_mom   *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der_pos = (1.0/3.0)*(4.0*Current_der2_pos - Current_der1_pos);
      Current_der_mom = (1.0/3.0)*(4.0*Current_der2_mom - Current_der1_mom);
      
      if(counter == 1) {
	Prev_der_pos = Current_der2_pos;
	Prev_der_mom = Current_der2_mom;
      }
    }
    else {
      Current_der_pos = Current_der2_pos;
      Current_der_mom = Current_der2_mom;
    }
    
    diff_pos = Current_der_pos - Prev_der_pos;
    if(Current_der_pos.Mag() > 1.0e-8) diff_pos *= (1.0/Current_der_pos.Mag());
    
    diff_mom = Current_der_mom - Prev_der_mom;
    if(Current_der_mom.Mag() > 1.0e-8) diff_mom *= (1.0/Current_der_mom.Mag());
    
    h /= 2.0;
    Prev_der_pos = Current_der_pos;
    Prev_der_mom = Current_der_mom;
  }
  while(counter <= Nmax_iterations && (diff_pos.Mag() > epsilon_derivatives || diff_mom.Mag() > epsilon_derivatives));
  
  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordAndMomDerWRTMom_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Position difference magnitude is " << diff_pos.Mag() << ", Momentum difference magnitude is " << diff_mom.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  DerCoor = Current_der_pos;
  DerMom  = Current_der_mom;
  
  return;
  
}
//====================================================================
void      GTracker::GetTrueTrackInterCoordAndMomDerWRTMom_Numerical(double s0, GGeoObject* aGeoElement, int der_index,
								    TVector3& DerCoor,
								    TVector3& DerMom)
{

  //Return the derivative of the track position and track momentum w.r.t. initial momentum with numerical method
  
  if(der_index <0 || der_index > 2) {
    cout << endl;
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTMom_Numerical:: ";
    if(der_index  > 2) cout << "derivative index with value " << der_index  << " is ouside range (0,2). Exting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }

  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  TVector3 InitMomentum_p;
  TVector3 InitMomentum_m;
  TVector3 InitMomentum_p2;
  TVector3 InitMomentum_m2;
  
  //double h = 1.0e-2*global->GetMomentumUnit("MeV/c");
  double h = 1.0e-3*InitMomentum.Mag();
  
  TVector3 Prev_der_pos(-999.9,-999.9,-999.9);
  TVector3 Current_der_pos(-999.9,-999.9,-999.9);
  TVector3 Current_der1_pos,Current_der2_pos;
  
  TVector3 Prev_der_mom(-999.9,-999.9,-999.9);
  TVector3 Current_der_mom(-999.9,-999.9,-999.9);
  TVector3 Current_der1_mom,Current_der2_mom;
  
  int counter = 0;
  TVector3 diff_pos(0.0,0.0,0.0);
  TVector3 diff_mom(0.0,0.0,0.0);
  
  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  do {
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;
    
    if(counter == 1 && DoRichardson) {
      InitMomentum_p  = InitMomentum;
      InitMomentum_m  = InitMomentum;
      InitMomentum_p(der_index) += 0.5*h1;
      InitMomentum_m(der_index) -= 0.5*h1;
    }
    
    InitMomentum_p2  = InitMomentum;
    InitMomentum_m2  = InitMomentum;
    InitMomentum_p2(der_index) += 0.5*h2;
    InitMomentum_m2(der_index) -= 0.5*h2;
    
    if(DoRichardson) {
      if(counter == 1) {
	SetTrajectoryInitConditions(InitPosition,InitMomentum_p);
	double s_tmp_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
	Current_der1_pos   = Trajectory->GetTrueTrajectoryCoordinates(s_tmp_p);
	Current_der1_mom   = Trajectory->GetTrueTrajectoryMon(s_tmp_p);
	
	SetTrajectoryInitConditions(InitPosition,InitMomentum_m);
	double s_tmp_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
        Current_der1_pos  -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp_m);
	Current_der1_mom  -= Trajectory->GetTrueTrajectoryMon(s_tmp_m);
	
        Current_der1_pos  *= (1.0/h1);
	Current_der1_mom  *= (1.0/h1);
      }
      else {
	Current_der1_pos = Prev_der_pos;
	Current_der1_mom = Prev_der_mom;
      }
    }
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum_p2);
    double s_tmp2_p = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2_pos    = Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_p);
    Current_der2_mom    = Trajectory->GetTrueTrajectoryMon(s_tmp2_p);
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum_m2);
    double s_tmp2_m = GetIntersectionsWithGeoElement(s0,aGeoElement);
    Current_der2_pos   -= Trajectory->GetTrueTrajectoryCoordinates(s_tmp2_m);
    Current_der2_mom   -= Trajectory->GetTrueTrajectoryMon(s_tmp2_m);
    
    Current_der2_pos   *= (1.0/h2);
    Current_der2_mom   *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der_pos = (1.0/3.0)*(4.0*Current_der2_pos - Current_der1_pos);
      Current_der_mom = (1.0/3.0)*(4.0*Current_der2_mom - Current_der1_mom);
      if(counter == 1) {
	Prev_der_pos = Current_der2_pos;
	Prev_der_mom = Current_der2_mom;
      }
    }
    else {
      Current_der_pos = Current_der2_pos;
      Current_der_mom = Current_der2_mom;
    }
    
    diff_pos = Current_der_pos - Prev_der_pos;
    if(Current_der_pos.Mag() > 1.0e-8) diff_pos *= (1.0/Current_der_pos.Mag());
    
    diff_mom = Current_der_mom - Prev_der_mom;
    if(Current_der_mom.Mag() > 1.0e-8) diff_mom *= (1.0/Current_der_mom.Mag());
    
    h /= 2.0;
    Prev_der_pos = Current_der_pos;
    Prev_der_mom = Current_der_mom;
  }
  while(counter <= Nmax_iterations && (diff_pos.Mag() > epsilon_derivatives));
  
  if(counter >= Nmax_iterations) {
    cout << "WARNNIN inside GTracker::GetTrueTrackInterCoordDerWRTMom_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Position difference magnitude is " << diff_pos.Mag() << ", Momentum difference magnitude is " << diff_mom.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  DerCoor = Current_der_pos;
  DerMom  = Current_der_mom;
  
  return;
  
}
//====================================================================
double    GTracker::GetFitTrackIntersectionCoordinates(GSurfaceObject* aSurface,double dummy_init, double sinit)
{
  
  //Calculate the value of the trajectory dummy parameter for the intersection between the track and a plane
  
  TVector3 Intersect(0.0,0.0,0.0);
  
  const int MaxIterations(MaxIterations_intersection);
  int counter = 0;
  double dummy0,dummy1;
  double delta0;
  double Derdelta0;
  double deltaODerdelta;
  
  double my_epsilon   = Trajectory->GetDummyParamEpsilon();
  TString dummy_units = Trajectory->GetDummyUnits();
 
  counter = 0;
  double w = 0.0;
  dummy0 = dummy_init;
    
  TVector3 pos         = Trajectory->GetFitTrackCoordinates(    dummy0, sinit);
  TVector3 Derpos      = Trajectory->GetFitTrackCoorDerWRTDummy(dummy0, sinit);
  TVector3 posUVW      = aSurface->GetUVWFromXYZ(pos);
  double UnitMon[3];
  UnitMon[0] = Derpos.X();
  UnitMon[1] = Derpos.Y();
  UnitMon[2] = Derpos.Z();
  double DerposW = 0.0;
  for(int i=0;i<3;i++) DerposW += aSurface->GetDerUVWFromXYZ_Analytical(pos,i).Z()*UnitMon[i];
  
  delta0    = w  - posUVW.Z();
  Derdelta0 =    - DerposW;
    
  if(TMath::Abs(Derdelta0) < 1.0e-10) deltaODerdelta = -1.0;
  else                                deltaODerdelta = delta0/Derdelta0;
  
  dummy1 = dummy0 - deltaODerdelta;
  dummy1 = Trajectory->DummyValueCorrection(dummy1);

  while(TMath::Abs(dummy0 - dummy1) > my_epsilon && counter <= MaxIterations) {
    dummy0 = dummy1;
    pos         = Trajectory->GetFitTrackCoordinates(    dummy0, sinit);
    Derpos      = Trajectory->GetFitTrackCoorDerWRTDummy(dummy0, sinit);
    posUVW      = aSurface->GetUVWFromXYZ(pos);
    UnitMon[0] = Derpos.X();
    UnitMon[1] = Derpos.Y();
    UnitMon[2] = Derpos.Z();
    DerposW = 0.0;
    for(int i=0;i<3;i++) DerposW += aSurface->GetDerUVWFromXYZ_Analytical(pos,i).Z()*UnitMon[i];
    
    delta0    = w  - posUVW.Z();
    Derdelta0 =    - DerposW;
 
    if(TMath::Abs(Derdelta0) < 1.0e-10) deltaODerdelta = -1.0;
    else                                deltaODerdelta = delta0/Derdelta0;
      
    dummy1 = dummy0 - deltaODerdelta;
    dummy1 = Trajectory->DummyValueCorrection(dummy1);
    
    counter++;
  }
  
  if(TMath::Abs(dummy0 - dummy1) > my_epsilon) {
    if(verbose) {
      cout << "WARNNING in function GTracker::GetFitTrackIntersectionCoordinates, iteration difference " << TMath::Abs(dummy0 - dummy1)/global->GetUnit(dummy_units) 
           << " " << dummy_units.Data() << " is bigger than limit " << my_epsilon/global->GetUnit(dummy_units) << " " << dummy_units.Data() << " after " << counter << " iterations. Result not precise!." 
	   << endl;
    }
    dummy1 = Dummy_value;
  }

  return dummy1;
  
}
//====================================================================
TVector3  GTracker::GetFitTrackIntersectionCoordinatesUVW(GSurfaceObject* aSurface,double dummy_init, double sinit)
{
  
  double dummy = GetFitTrackIntersectionCoordinates(aSurface,dummy_init,sinit);
  TVector3 IntersectionPosition = Trajectory->GetFitTrackCoordinates(dummy,sinit);
  
  return aSurface->GetUVWFromXYZ(IntersectionPosition);
  
}
//====================================================================
TVector3  GTracker::GetFitTrackIntersectionCoordUVWDerWRTPar_Numerical(GSurfaceObject* aSurface, int  parameter, double dummy_init, double sinit)
{
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  const int Npars(Trajectory->GetNParameters());
  
  if(parameter < 0 || parameter > Npars-1) {
    cout << endl;
    cout << "ERROR in GTracker::GetTrackDerivatives: parameter index " << parameter << " is out of limits (0," << Npars-1 << "). Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  TVector3 IntersectionUVW = aSurface->GetUVWFromXYZ(Trajectory->GetTrueTrajectoryCoordinates(sinit));
  
  double h = Trajectory->GetInitValueParForDerivativeCalculation(parameter);
  
  TVector3 Prev_der(-999.9,-999.9,-999.9);
  TVector3 Current_der(-999.9,-999.9,-999.9);
  TVector3 Current_der1,Current_der2;
  
  int counter = 0;
  TVector3 diff(0.0,0.0,0.0);

  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  std::vector<double> FitParams;
  Trajectory->GetAllParameters(FitParams);
  
  TVector3 tmpVector;
  
  do {
    SetTrajectoryInitConditions(InitPosition,InitMomentum);
    
    counter++;
    
    double h1 = h;
    double h2 = h/2.0;
    
    std::vector<double> FitParams_p;
    std::vector<double> FitParams_m;
    std::vector<double> FitParams_p2;
    std::vector<double> FitParams_m2;
    FitParams_p.clear();
    FitParams_m.clear();
    FitParams_p2.clear();
    FitParams_m2.clear();
    
    if(counter == 1 && DoRichardson) {
      Trajectory->GetAllParameters(FitParams_p);
      Trajectory->GetAllParameters(FitParams_m);
      FitParams_p[parameter]  += 0.5*h1;
      FitParams_m[parameter]  -= 0.5*h1;
    }
    Trajectory->GetAllParameters(FitParams_p2);
    Trajectory->GetAllParameters(FitParams_m2);
    FitParams_p2[parameter] += 0.5*h2;
    FitParams_m2[parameter] -= 0.5*h2;
    
    if(DoRichardson) {
      if(counter == 1) {
	Trajectory->SetAllParameters(FitParams_p);
	tmpVector = GetFitTrackIntersectionCoordinatesUVW(aSurface,dummy_init,sinit);
        Current_der1  = tmpVector;
	
	Trajectory->SetAllParameters(FitParams_m);
	tmpVector = GetFitTrackIntersectionCoordinatesUVW(aSurface,dummy_init,sinit);
        Current_der1 -= tmpVector;
	
        Current_der1 *= (1.0/h1);
      }
      else Current_der1 = Prev_der;
    }
    
    //cout << endl;
    //cout << "Interaction s = " << sinit/global->GetUnit("cm") << " cm, UVW-Derivative wrt parameter " << Trajectory->GetParameterName(parameter).Data() << " iteration " << counter << endl;
    //cout << "par val (" << Trajectory->GetParameterUnitTitle(parameter).Data() << ") = " << FitParams[parameter]/global->GetUnit(Trajectory->GetParameterUnit(parameter)) << ", "
    //     << "h2 = " << h2/global->GetUnit(Trajectory->GetParameterUnit(parameter)) << ", "
    //     << "positive/negative variations = " << FitParams_p2[parameter]/global->GetUnit(Trajectory->GetParameterUnit(parameter)) << "/" << FitParams_m2[parameter]/global->GetUnit(Trajectory->GetParameterUnit(parameter))
    //     << endl;
    //cout << "Nominal intersection coordinates (U,V,W)  = (" << IntersectionUVW(0)/global->GetUnit("cm") << "," << IntersectionUVW(1)/global->GetUnit("cm") << "," << IntersectionUVW(2)/global->GetUnit("cm") << ") cm" << endl;
    
    Trajectory->SetAllParameters(FitParams_p2);
    tmpVector = GetFitTrackIntersectionCoordinatesUVW(aSurface,dummy_init,sinit);
    Current_der2  = tmpVector;
    //cout << "Positive intersection coordinates (U,V,W) = (" << tmpVector(0)/global->GetUnit("cm") << "," << tmpVector(1)/global->GetUnit("cm") << "," << tmpVector(2)/global->GetUnit("cm") << ") cm" << endl;
    
    Trajectory->SetAllParameters(FitParams_m2);
    tmpVector = GetFitTrackIntersectionCoordinatesUVW(aSurface,dummy_init,sinit);
    Current_der2 -= tmpVector;
    //cout << "Positive intersection coordinates (U,V,W) = (" << tmpVector(0)/global->GetUnit("cm") << "," << tmpVector(1)/global->GetUnit("cm") << "," << tmpVector(2)/global->GetUnit("cm") << ") cm" << endl;
    
    Current_der2 *= (1.0/h2);
    
    if(DoRichardson) {
      Current_der = (1.0/3.0)*(4.0*Current_der2 - Current_der1);
      if(counter == 1) Prev_der = Current_der2;
    }
    else Current_der = Current_der2;
    
    //cout << "Derivative (U,V,W) = (" 
    //     << Current_der(0)*(global->GetUnit(Trajectory->GetParameterUnit(parameter))/global->GetUnit("mm")) << ","
	// << Current_der(1)*(global->GetUnit(Trajectory->GetParameterUnit(parameter))/global->GetUnit("mm")) << "," 
	// << Current_der(2)*(global->GetUnit(Trajectory->GetParameterUnit(parameter))/global->GetUnit("mm")) << ") " 
	// << "mm/" << Trajectory->GetParameterUnitTitle(parameter).Data() 
	// << endl;
    
    //cout << endl;
    
    diff = Current_der - Prev_der;
    if(Current_der.Mag() > 1.0e-8) diff *= (1.0/Current_der.Mag());
    
    h /= 2.0;
    Prev_der = Current_der;
  }
  while(counter <= Nmax_iterations && diff.Mag() > epsilon_derivatives);

  if(counter >= Nmax_iterations) {
    cout << "WARNNING in function GTracker::GetFitTrackIntersectionCoordUVDerWRTPar_Numerical:: ";
    cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << diff.Mag() << endl;
  }
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return Current_der;
  
}
//====================================================================
void  GTracker::GetAllTrackParamDerWrtMS(std::vector<IntersectionHit_t>  ItersectionHitList,
					 double*                         SigmaMS2,
					 TVector3*                       m1Unit,
					 TVector3*                       m2Unit,
					 double***                       Der_Params_prime_wrt_thetaMS)
{
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  const int Npars(Trajectory->GetNParameters());
  const int NLayers(ItersectionHitList.size()); // Number of intersection geometry layers by the particle in global geometrys

  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  TString MSthetaUnit("urad");
  TString* ParUnit = new TString[Npars];
  for(int ipar=0;ipar<Npars;ipar++) ParUnit[ipar] = Trajectory->GetParameterUnit(ipar);
  
  double Diff_Der_Magnitude;
  double Params_prime_wrt_thetaMS_prev[Npars];     //Previous derivative
  double Params_prime_wrt_thetaMS_current[Npars];  //Current  derivative
  double Params_prime_wrt_thetaMS_current1[Npars]; //Current  derivative
  double Params_prime_wrt_thetaMS_current2[Npars]; //Current  derivative
  
  for(int nlayer=0;nlayer<NLayers;nlayer++) { //layer loop
    TVector3 Mym[2];
    Mym[0] = m1Unit[nlayer];
    Mym[1] = m2Unit[nlayer];
    double MySigmaMS2 = SigmaMS2[nlayer];
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum);
    double s         = ItersectionHitList[nlayer].s;
    TVector3 r1      = Trajectory->GetTrueTrajectoryCoordinates(s);
    TVector3 p1      = Trajectory->GetTrueTrajectoryMon(s);
    TVector3 p1_unit = p1.Unit();
    
    double theta_p[2]; // positive variations of MS angles
    double theta_m[2]; // negative variations of MS angles
    for(int itheta_ms=0;itheta_ms<2;itheta_ms++) { //Begin of loop on MS angles
      //Current and previous derivatives Initialization to dummy values
      for(int ipar=0;ipar<Npars;ipar++) Params_prime_wrt_thetaMS_prev[ipar] = Params_prime_wrt_thetaMS_current[ipar] = -999.0;

      //step parameter Initialization for derivative calculation
      double h = sqrt(MySigmaMS2);
      int counter = 0;
      do {
        counter++;
      
        double h1 = h;
        double h2 = h/2.0;

        //Calculation of the new track parameters using the variated momentum
        std::vector<double>  FitParams_tmp_p1;
        std::vector<double>  FitParams_tmp_m1;
        std::vector<double>  FitParams_tmp_p2;
        std::vector<double>  FitParams_tmp_m2;
      
        if(counter == 1 && DoRichardson) {
	  //positive and negative variations of the MS angles
          for(int jtheta_ms=0;jtheta_ms<2;jtheta_ms++) theta_p[jtheta_ms] = theta_m[jtheta_ms] = 0.0;
          theta_p[itheta_ms] = +0.5*h1;
          theta_m[itheta_ms] = -0.5*h1;
      
          //Variated momentum vectors after passing the surface APlane
          TVector3 mom_p = p1_unit + Mym[0]*theta_p[0] + Mym[1]*theta_p[1];
          TVector3 mom_m = p1_unit + Mym[0]*theta_m[0] + Mym[1]*theta_m[1];
	  mom_p *= p1.Mag();
	  mom_m *= p1.Mag();
	
	  SetTrajectoryInitConditions(r1,mom_p);
          Trajectory->GetFitParsFromInitConds();
	  Trajectory->GetAllParameters(FitParams_tmp_p1);
	  
	  SetTrajectoryInitConditions(r1,mom_m);
          Trajectory->GetFitParsFromInitConds();
	  Trajectory->GetAllParameters(FitParams_tmp_m1);
        }
      
        //positive and negative variations of the MS angles
        for(int jtheta_ms=0;jtheta_ms<2;jtheta_ms++) theta_p[jtheta_ms] = theta_m[jtheta_ms] = 0.0;
        theta_p[itheta_ms] = +0.5*h2;
        theta_m[itheta_ms] = -0.5*h2;
      
        //Variated momentum vectors after passing the surface APlane
        TVector3 mom_p = p1_unit + Mym[0]*theta_p[0] + Mym[1]*theta_p[1];
        TVector3 mom_m = p1_unit + Mym[0]*theta_m[0] + Mym[1]*theta_m[1];
	mom_p *= p1.Mag();
	mom_m *= p1.Mag();
	
	SetTrajectoryInitConditions(r1,mom_p);
	Trajectory->GetFitParsFromInitConds();
	Trajectory->GetAllParameters(FitParams_tmp_p2);
	
	SetTrajectoryInitConditions(r1,mom_m);
	Trajectory->GetFitParsFromInitConds();
	Trajectory->GetAllParameters(FitParams_tmp_m2);

        //Numerical calculation of the new parameters derivatives for step parameter h
        Diff_Der_Magnitude = 0.0;
        for(int ipar=0;ipar<Npars;ipar++) {
	  if(DoRichardson) {
            if(counter == 1) {
              Params_prime_wrt_thetaMS_current1[ipar]  = FitParams_tmp_p1[ipar];
              Params_prime_wrt_thetaMS_current1[ipar] -= FitParams_tmp_m1[ipar];
              Params_prime_wrt_thetaMS_current1[ipar] /= h1;
            }
            else Params_prime_wrt_thetaMS_current1[ipar] = Params_prime_wrt_thetaMS_prev[ipar];
          }
        
          Params_prime_wrt_thetaMS_current2[ipar]  = FitParams_tmp_p2[ipar];
          Params_prime_wrt_thetaMS_current2[ipar] -= FitParams_tmp_m2[ipar];
          Params_prime_wrt_thetaMS_current2[ipar] /= h2;
	
	  if(DoRichardson) {
            Params_prime_wrt_thetaMS_current[ipar] = (1.0/3.0)*(4.0*Params_prime_wrt_thetaMS_current2[ipar] - Params_prime_wrt_thetaMS_current1[ipar]);
            if(counter == 1) Params_prime_wrt_thetaMS_prev[ipar] = Params_prime_wrt_thetaMS_current2[ipar];
          }
          else Params_prime_wrt_thetaMS_current[ipar] = Params_prime_wrt_thetaMS_current2[ipar];
	
	  double delta = Params_prime_wrt_thetaMS_current[ipar] - Params_prime_wrt_thetaMS_prev[ipar];
	  double unit_factor = global->GetUnit(MSthetaUnit)/global->GetUnit(ParUnit[ipar]);
	  if(Params_prime_wrt_thetaMS_current[ipar]*unit_factor > 1.0e-8) delta /= Params_prime_wrt_thetaMS_current[ipar];
	  Diff_Der_Magnitude = TMath::Max(TMath::Abs(delta),Diff_Der_Magnitude);
	
	  Params_prime_wrt_thetaMS_prev[ipar] = Params_prime_wrt_thetaMS_current[ipar];
        }
        
        //Re-sizing of step parameter
        h /= 2.0;
      }
      while(counter <= Nmax_iterations && Diff_Der_Magnitude > epsilon_derivatives);

      if(counter >= Nmax_iterations) {
        cout << "WARNNING in function GTracker::GetAllTrackParamDerWrtMS:: ";
        cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << Diff_Der_Magnitude << " for itheta_ms = " << itheta_ms << endl;
      }
    
      for(int ipar=0;ipar<Npars;ipar++) Der_Params_prime_wrt_thetaMS[nlayer][ipar][itheta_ms] = Params_prime_wrt_thetaMS_current[ipar];
    } //End of loop on MS angles
    SetTrajectoryInitConditions(InitPosition,InitMomentum);
  } //End of layer loop
  
  delete [] ParUnit;
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return;
  
}
//====================================================================
void  GTracker::GetAllTrackParamDerWrtEloss(std::vector<IntersectionHit_t>  ItersectionHitList,
					    double*                         SigmaEloss2,
					    TVector3*                       pUnit,
					    double**                        Der_Params_prime_wrt_Eloss)
{
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  //double   EOp2         = (InitMomentum.Mag2() + pow(global->GetParticleMass(Trajectory->GetParticle()),2))/InitMomentum.Mag2();
  
  const int Npars(Trajectory->GetNParameters());
  const int NLayers(ItersectionHitList.size()); // Number of intersection geometry layers by the particle in global geometrys

  bool DoRichardson = false;
  DoRichardson      = global->GlobalDoRichardson;
  
  TString particle = Trajectory->GetParticle();
  
  TString ElossUnit("MeV/c");
  TString* ParUnit = new TString[Npars];
  for(int ipar=0;ipar<Npars;ipar++) ParUnit[ipar] = Trajectory->GetParameterUnit(ipar);
  
  double Diff_Der_Magnitude;
  double Params_prime_wrt_Eloss_prev[Npars];     //Previous derivative
  double Params_prime_wrt_Eloss_current[Npars];  //Current  derivative
  double Params_prime_wrt_Eloss_current1[Npars]; //Current  derivative
  double Params_prime_wrt_Eloss_current2[Npars]; //Current  derivative

  for(int nlayer=0;nlayer<NLayers;nlayer++) { //layer loop
    //TVector3 Mym = pUnit[nlayer];
    //double MySigmaP2 = EOp2*SigmaEloss2[nlayer];

    SetTrajectoryInitConditions(InitPosition,InitMomentum);
    double s         = ItersectionHitList[nlayer].s;
    TVector3 r1      = Trajectory->GetTrueTrajectoryCoordinates(s);
    TVector3 p1      = Trajectory->GetTrueTrajectoryMon(s);
    TVector3 p1_unit = p1.Unit();
    
    //Current and previous derivatives Initialization to dummy values
    for(int ipar=0;ipar<Npars;ipar++) Params_prime_wrt_Eloss_prev[ipar] = Params_prime_wrt_Eloss_current[ipar] = -999.0;

    //step parameter Initialization for derivative calculation
    double h = 1.0e-3*InitMomentum.Mag();
    int counter = 0;
    do {
      counter++;
      
      double h1 = h;
      double h2 = h/2.0;

      //Calculation of the new track parameters using the variated momentum
      std::vector<double>  FitParams_tmp_p1;
      std::vector<double>  FitParams_tmp_m1;
      std::vector<double>  FitParams_tmp_p2;
      std::vector<double>  FitParams_tmp_m2;
      
      if(counter == 1 && DoRichardson) {
	//Variated momentum vectors after passing the surface APlane
        TVector3 mom_p = (p1.Mag() + 0.5*h1)*p1_unit;
        TVector3 mom_m = (p1.Mag() - 0.5*h1)*p1_unit;
	
	SetTrajectoryInitConditions(r1,mom_p);
        Trajectory->GetFitParsFromInitConds();
	Trajectory->GetAllParameters(FitParams_tmp_p1);
	  
	SetTrajectoryInitConditions(r1,mom_m);
        Trajectory->GetFitParsFromInitConds();
	Trajectory->GetAllParameters(FitParams_tmp_m1);
      }
      
      //Variated momentum vectors after passing the surface APlane
      TVector3 mom_p = (p1.Mag() + 0.5*h2)*p1_unit;
      TVector3 mom_m = (p1.Mag() - 0.5*h2)*p1_unit;
	
      SetTrajectoryInitConditions(r1,mom_p);
      Trajectory->GetFitParsFromInitConds();
      Trajectory->GetAllParameters(FitParams_tmp_p2);
	
      SetTrajectoryInitConditions(r1,mom_m);
      Trajectory->GetFitParsFromInitConds();
      Trajectory->GetAllParameters(FitParams_tmp_m2);

      //Numerical calculation of the new parameters derivatives for step parameter h
      Diff_Der_Magnitude = 0.0;
      for(int ipar=0;ipar<Npars;ipar++) {
	if(DoRichardson) {
          if(counter == 1) {
            Params_prime_wrt_Eloss_current1[ipar]  = FitParams_tmp_p1[ipar];
            Params_prime_wrt_Eloss_current1[ipar] -= FitParams_tmp_m1[ipar];
            Params_prime_wrt_Eloss_current1[ipar] /= h1;
          }
          else Params_prime_wrt_Eloss_current1[ipar] = Params_prime_wrt_Eloss_prev[ipar];
        }
        
        Params_prime_wrt_Eloss_current2[ipar]  = FitParams_tmp_p2[ipar];
        Params_prime_wrt_Eloss_current2[ipar] -= FitParams_tmp_m2[ipar];
        Params_prime_wrt_Eloss_current2[ipar] /= h2;
	
	if(DoRichardson) {
          Params_prime_wrt_Eloss_current[ipar] = (1.0/3.0)*(4.0*Params_prime_wrt_Eloss_current2[ipar] - Params_prime_wrt_Eloss_current1[ipar]);
          if(counter == 1) Params_prime_wrt_Eloss_prev[ipar] = Params_prime_wrt_Eloss_current2[ipar];
        }
        else Params_prime_wrt_Eloss_current[ipar] = Params_prime_wrt_Eloss_current2[ipar];
	
	double delta = Params_prime_wrt_Eloss_current[ipar] - Params_prime_wrt_Eloss_prev[ipar];
	double unit_factor = global->GetUnit(ElossUnit)/global->GetUnit(ParUnit[ipar]);
	if(Params_prime_wrt_Eloss_current[ipar]*unit_factor > 1.0e-8) delta /= Params_prime_wrt_Eloss_current[ipar];
	Diff_Der_Magnitude = TMath::Max(TMath::Abs(delta),Diff_Der_Magnitude);
	
	Params_prime_wrt_Eloss_prev[ipar] = Params_prime_wrt_Eloss_current[ipar];
      }
        
      //Re-sizing of step parameter
      h /= 2.0;
    }
    while(counter <= Nmax_iterations && Diff_Der_Magnitude > epsilon_derivatives);

    if(counter >= Nmax_iterations) {
      cout << "WARNNING in function GTracker::GetAllTrackParamDerWrtEloss:: ";
      cout << "iterations reached maximum " << Nmax_iterations << ". Difference magnitude is " << Diff_Der_Magnitude << endl;
    }
    
    for(int ipar=0;ipar<Npars;ipar++) Der_Params_prime_wrt_Eloss[nlayer][ipar] = Params_prime_wrt_Eloss_current[ipar];
    
    SetTrajectoryInitConditions(InitPosition,InitMomentum);
  } //End of layer loop
  
  delete [] ParUnit;
  
  SetTrajectoryInitConditions(InitPosition,InitMomentum);

  return;
  
}
//====================================================================
void  GTracker::GetNewParCovarianceMatrix(std::vector<IntersectionHit_t>  ItersectionHitList_global,
					  std::vector<IntersectionHit_t>  ItersectionHitList_current,
					  int                             nlayer_current,
					  double***                       Der_Params_prime_wrt_thetaMS,
					  double*                         SigmaMS2,
					  double**                        Der_Params_prime_wrt_Eloss,
					  double*                         SigmaEloss2,
					  TMatrixD                        FitCovMatrix_old,
					  TMatrixD&                       FitCovMatrix_new)
{
  
  int nlayer_global = -999;
  for(int klayer=0;klayer<int(ItersectionHitList_global.size());klayer++) {
    if(ItersectionHitList_global[klayer].geoElement_idx == ItersectionHitList_current[nlayer_current].geoElement_idx) {
      nlayer_global = klayer;
      break;
    }
  }
  if(nlayer_global == -999) {
    cout << endl;
    cout << "ERROR inside GTracker::GetNewParCovarianceMatrix:: global layer not found. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  //Now calculate the new track parameters covariance matrix including the new material MS
  const int Npars(Trajectory->GetNParameters());
  FitCovMatrix_new.ResizeTo(Npars,Npars);
  for(int ipar=0;ipar<Npars;ipar++) {
    for(int jpar=0;jpar<Npars;jpar++) {
      FitCovMatrix_new(ipar,jpar) = 0.0;
      
      //First include the multiple scattering component
      for(int itheta_ms=0;itheta_ms<2;itheta_ms++) {
	FitCovMatrix_new(ipar,jpar) += Der_Params_prime_wrt_thetaMS[nlayer_global][ipar][itheta_ms]*Der_Params_prime_wrt_thetaMS[nlayer_global][jpar][itheta_ms];
      }
      FitCovMatrix_new(ipar,jpar) *= SigmaMS2[nlayer_global];
    }
  }
  
  //Now calculate the new track parameters covariance matrix including the new material Eloss
  if(IncludeEloss) {
    TVector3 InitMomentum = Trajectory->GetInitMomentum();
    double   EOp2         = (InitMomentum.Mag2() + pow(global->GetParticleMass(Trajectory->GetParticle()),2))/InitMomentum.Mag2();
    
    for(int ipar=0;ipar<Npars;ipar++) {
      for(int jpar=0;jpar<Npars;jpar++) {
	FitCovMatrix_new(ipar,jpar) += Der_Params_prime_wrt_Eloss[nlayer_global][ipar]*Der_Params_prime_wrt_Eloss[nlayer_global][jpar]*EOp2*SigmaEloss2[nlayer_global];
      }
    }
  }
  
  //Add the previous parameters covariance matrix
  FitCovMatrix_new += FitCovMatrix_old;
  
  return;
  
}
//====================================================================
void  GTracker::GetCovMatrixUVOfTrackIntersection(std::vector<IntersectionHit_t>  ItersectionHitList_global,
						  std::vector<IntersectionHit_t>  ItersectionHitList_current,
						  int                             nlayer_current,
						  double**                        GammaUDer,
						  double**                        GammaVDer,
						  TMatrixD                        FitCovMatrix,
						  TMatrixD&                       CovUV)
{
  
  //Calculate the covariance matrix of the local coordinates of the track intersection with surface APlane
  //This cov matrix is the propagation of the track parameters uncertainties into the track intersection local coordinates at surface APlane
  
  int nlayer_global = -999;
  for(int klayer=0;klayer<int(ItersectionHitList_global.size());klayer++) {
    if(ItersectionHitList_global[klayer].geoElement_idx == ItersectionHitList_current[nlayer_current].geoElement_idx) {
      nlayer_global = klayer;
      break;
    }
  }
  if(nlayer_global == -999) {
    cout << endl;
    cout << "ERROR inside GTracker::GetCovMatrixUVOfTrackIntersection:: global layer not found. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  const int Npars(Trajectory->GetNParameters());
  
  TVector3 pos_n = Trajectory->GetTrueTrajectoryCoordinates(ItersectionHitList_global[nlayer_global].s);
  TVector3 mom_n = Trajectory->GetTrueTrajectoryMon(        ItersectionHitList_global[nlayer_global].s);
  if(UseMonOrigin) mom_n = Trajectory->GetInitPosition();
  
  double ResolutionU = 0.0;
  double ResolutionV = 0.0;
  //int geoElement_idx = ItersectionHitList_global[nlayer_global].geoElement_idx;
  //TVector2 Resolution = Geometry->GetResolutionUV(geoElement_idx,mom_n,pos_n);
  //ResolutionU = Resolution.X();
  //ResolutionV = Resolution.Y();
  
  //initialize local coordinates covariance matrix
  CovUV.ResizeTo(2,2);
  for(int i=0;i<2;i++) {
    for(int j=0;j<2;j++) {
      CovUV(i,j) = 0.0;
      if(i != j) continue;
      //Intrinsique layer point resolution
      
      //Still need to specify the momentum and coordinates of the hit to the function
      if(i == 0) CovUV(i,i) = pow(ResolutionU,2);
      if(i == 1) CovUV(i,i) = pow(ResolutionV,2);
    }
  }
  
  //Add the contribution from the track parameter uncertainty
  for(int ipar=0;ipar<Npars;ipar++) {
    for(int jpar=0;jpar<Npars;jpar++) {
      CovUV(0,0) += GammaUDer[nlayer_global][ipar]*FitCovMatrix(ipar,jpar)*GammaUDer[nlayer_global][jpar];
      CovUV(0,1) += GammaUDer[nlayer_global][ipar]*FitCovMatrix(ipar,jpar)*GammaVDer[nlayer_global][jpar];
      CovUV(1,0) += GammaVDer[nlayer_global][ipar]*FitCovMatrix(ipar,jpar)*GammaUDer[nlayer_global][jpar];
      CovUV(1,1) += GammaVDer[nlayer_global][ipar]*FitCovMatrix(ipar,jpar)*GammaVDer[nlayer_global][jpar];
    }
  }
  
  return;
  
}
//====================================================================
void  GTracker::GetUVDerivatiesWrtPars(std::vector<IntersectionHit_t>  ItersectionHitList,
				       double**                        GammaUDer,
				       double**                        GammaVDer,
				       const int                       Npars)
{
  
  for(int nhit=0;nhit<int(ItersectionHitList.size());nhit++) { //1st loop on the hit local coordinates
    //Get the hit coordiates in the lab frame for hit nhit
    int nhit_layer     = nhit;
    int geoElement_idx = ItersectionHitList[nhit_layer].geoElement_idx;
    double sn          = ItersectionHitList[nhit_layer].s;
    TVector3 pos_n     = Trajectory->GetTrueTrajectoryCoordinates(sn);
    double dummy       = Trajectory->GetFitTrackDummyParFromS(sn);

    GSurfaceObject*  aSurface = Geometry->GetGeometryElement(geoElement_idx)->GetMainSurface();
    
#if 0
    double dummy_tmp = GetFitTrackIntersectionCoordinates(aSurface,dummy,sn);
    TVector3 pos_dummy = (1.0/global->GetUnit("cm"))*Trajectory->GetFitTrackCoordinates(dummy_tmp,sn);
    TVector3 pos_s     = (1.0/global->GetUnit("cm"))*pos_n;
    cout << "nhit = " << nhit << ", "
         << "s = " << sn/global->GetUnit("cm") << " cm, "
	 << "dummy_init = " << dummy/global->GetUnit("cm") << " cm, "
	 << "dummy_tmp = "  << dummy_tmp/global->GetUnit("cm") << " cm, "
         << "pos_n = ("     << pos_s(0)     << "," << pos_s(1)     << "," << pos_s(2)     << ") cm; "
	 << "pos_dummy = (" << pos_dummy(0) << "," << pos_dummy(1) << "," << pos_dummy(2) << ") cm; "
	 << "diff = ("      << pos_s(0) - pos_dummy(0) << "," << pos_s(1) - pos_dummy(1) << "," << pos_s(2) - pos_dummy(2) << ") cm; "
	 << endl;
#endif
    
    for(int ipar=0;ipar<Npars;ipar++) {  //loop on track parameters      
      TVector3 GammaUVWDer   = GetFitTrackIntersectionCoordUVWDerWRTPar_Numerical(aSurface,ipar,dummy,sn);
      GammaUDer[nhit][ipar]  = GammaUVWDer(0);
      GammaVDer[nhit][ipar]  = GammaUVWDer(1);
    } //end of loop on track parameters
  } // end of 1st loop of hit local coordinates
  
  return;
  
}
//====================================================================
void  GTracker::GetFsAndKsMS(std::vector<IntersectionHit_t>  ItersectionHitList,
			     TMatrixD*                       Fr,
			     TMatrixD*                       Fp,
			     TMatrixD*                       Kr,
			     TMatrixD*                       Kp)
{
  
  TVector3 InitPosition = Trajectory->GetInitPosition();
  TVector3 InitMomentum = Trajectory->GetInitMomentum();
  
  //Filling up of the F and K arrays
  int nCoor   = 3;
  int NLayers = ItersectionHitList.size();
  for(int nlayer=0;nlayer<NLayers;nlayer++) { //layer loop
    if(nlayer == 0) continue;

    double   s_nm1     = ItersectionHitList[nlayer-1].s;
    double   s_n       = ItersectionHitList[nlayer].s;
    TVector3 r_nm1     = Trajectory->GetTrueTrajectoryCoordinates(s_nm1);
    TVector3 p_nm1     = Trajectory->GetTrueTrajectoryMon(s_nm1);
    int geoElement_idx = ItersectionHitList[nlayer].geoElement_idx;
    
    GGeoObject* aGeoElement = Geometry->GetGeometryElement(geoElement_idx);

    SetTrajectoryInitConditions(r_nm1,p_nm1);
    for(int j=0;j<nCoor;j++) { //derivative loop
      TVector3  TrueTrackInterCoordDerWRTPos;
      TVector3  TrueTrackInterCoordDerWRTMom;
      TVector3  TrueTrackInterMomDerWRTPos;
      TVector3  TrueTrackInterMomDerWRTMom;
      
      GetTrueTrackInterCoordAndMomDerWRTPos_Numerical(s_n - s_nm1,aGeoElement,j,TrueTrackInterCoordDerWRTPos,TrueTrackInterMomDerWRTPos);
      GetTrueTrackInterCoordAndMomDerWRTMom_Numerical(s_n - s_nm1,aGeoElement,j,TrueTrackInterCoordDerWRTMom,TrueTrackInterMomDerWRTMom);
      
      for(int i=0;i<nCoor;i++) {  //coordinate loop
	Fr[nlayer](i,j) = TrueTrackInterCoordDerWRTPos(i);
	Fp[nlayer](i,j) = TrueTrackInterCoordDerWRTMom(i);
	Kr[nlayer](i,j) = TrueTrackInterMomDerWRTPos(i);
	Kp[nlayer](i,j) = TrueTrackInterMomDerWRTMom(i);
      }
    }
    SetTrajectoryInitConditions(InitPosition,InitMomentum);
  } // end of layer loop

  SetTrajectoryInitConditions(InitPosition,InitMomentum);
  
  return;
  
}
//====================================================================
void  GTracker::GetSigmaMSAndPerpVectors(std::vector<IntersectionHit_t> ItersectionHitList,
					 double*                        SigmaMS2,
					 TVector3*                      m1Unit,
					 TVector3*                      m2Unit)
{
  
  TString particle = Trajectory->GetParticle();
  
  int NLayers = ItersectionHitList.size();
  for(int nlayer=0;nlayer<NLayers;nlayer++) { //layer loop
    double s                 = ItersectionHitList[nlayer].s;
    int geoElement_idx       = ItersectionHitList[nlayer].geoElement_idx;
    TVector3 MonIntersection = Trajectory->GetTrueTrajectoryMon(s);
    
    TVector3 Mym1,Mym2;
    global->GetSurfAndMomOrthVects(MonIntersection,Mym1,Mym2);
    
    //m1 unit vector
    m1Unit[nlayer] = Mym1;
    
    //m2 unit vector
    m2Unit[nlayer] = Mym2;
    
    GGeoObject* aGeoElement = Geometry->GetGeometryElement(geoElement_idx);
    double XOX0 = GetIntersectionXOX0(ItersectionHitList[nlayer],aGeoElement);
    
    SigmaMS2[nlayer] = pow(global->MSAngle(particle,MonIntersection,XOX0),2);
    if(std::isnan(XOX0)) {
      TVector3 InitMomentum = Trajectory->GetInitMomentum();
      cout << "X/X0 is nan inside GetSigmaMSAndPerpVectors for " 
           << "s = " << s/global->GetDistanceUnit("cm") << " cm, "
	   << "geoElement_idx = " << geoElement_idx << " with name " << aGeoElement->GetName().Data() << ", "
	   << "plane material " << aGeoElement->GetMaterial().Data() << ", "
	   << "X/X0 for ref thickness = " << aGeoElement->GetXOX0() << ", " 
	   << "mom(p,theta,phi) = (" << InitMomentum.Mag()/global->GetMomentumUnit("GeV/c") << " GeV/c," 
	   << atan2(sqrt(pow(InitMomentum.X(),2) + pow(InitMomentum.Y(),2)),InitMomentum.Z())*180.0/TMath::Pi() << " deg,"
	   << atan2(InitMomentum.Y(),InitMomentum.X())*180.0/TMath::Pi() << " deg)"
           << endl;
    }
  }
  
  return;
  
}
//====================================================================
void  GTracker::GetSigmaElossAndUnitVector(std::vector<IntersectionHit_t> ItersectionHitList,
					   double*                        SigmaEloss2,
					   TVector3*                      pUnit)
{
  
  bool Myverbose = false;
  //Myverbose = true;
  
  TString particle = Trajectory->GetParticle();
  
  int NLayers = ItersectionHitList.size();
  if(Myverbose) cout << endl;
  for(int nlayer=0;nlayer<NLayers;nlayer++) { //layer loop
    double s                 = ItersectionHitList[nlayer].s;
    int geoElement_idx       = ItersectionHitList[nlayer].geoElement_idx;
    TVector3 MonIntersection = Trajectory->GetTrueTrajectoryMon(s);
    
    //unit vector along the momentum at the intersection
    pUnit[nlayer]    = MonIntersection.Unit();
    
    TString material = Geometry->GetGeometryElement(geoElement_idx)->GetMaterial();
    double X         = TMath::Abs(ItersectionHitList[nlayer].s_in - ItersectionHitList[nlayer].s_out);    
    double dEOdx     = global->BetheBlochPDG(MonIntersection.Mag(),particle,material);
    
    double sigma_Ivanov = global->GetSigmaELoss(MonIntersection.Mag(),particle,material,dEOdx*X);
    double sigma_pdg    = global->GetSigmaELossPDG(MonIntersection.Mag(),particle,material,X);
    
    SigmaEloss2[nlayer] = pow(sigma_Ivanov,2);
    //SigmaEloss2[nlayer] = pow(sigma_pdg,2);
    
    if(Myverbose) {
      double E = sqrt(MonIntersection.Mag2() + pow(global->GetParticleMass(particle),2));
      TVector3 InitMomentum = Trajectory->GetInitMomentum();
      cout << "s = " << s/global->GetDistanceUnit("cm") << " cm, "
	   << "geoElement_idx = " << geoElement_idx << " with name " << Geometry->GetGeometryElement(geoElement_idx)->GetName().Data() << ", "
	   << "plane material " << material.Data() << ", "
	   << "X = " << X/global->GetDistanceUnit("um") << " um; "
	   << "dE/dx_geant = " << global->BetheBlochGeant(MonIntersection.Mag(),particle,material)*(global->GetUnit("um")/global->GetUnit("MeV")) << " MeV/um; "
	   << "dE/dx_pdg   = " << global->BetheBlochPDG(MonIntersection.Mag(),particle,material)*(global->GetUnit("um")/global->GetUnit("MeV")) << " MeV/um; "
	   << "dE/dx = " << dEOdx*(global->GetUnit("um")/global->GetUnit("MeV")) << " MeV/um; "
	   << "Eloss = " << dEOdx*X/global->GetUnit("MeV") << " MeV; "
	   << "sigma(Eloss) Ivanov = " << sigma_Ivanov/global->GetUnit("MeV") << " MeV; "
           << "sigma(Eloss) pdg = " << sigma_pdg/global->GetUnit("MeV") << " MeV; "
	   << "E = " << E/global->GetUnit("GeV") << " GeV; "
	   << "Eloss/E = " << 100*dEOdx*X/E << " %; "
	   << "sigmaEloss_Ivanov/E = " << 100*sigma_Ivanov/E << "%; "
	   << "sigmaEloss_pdg/E = " << 100*sigma_pdg/E << "%; "
	   << "pUnit = (" << pUnit[nlayer](0) << "," << pUnit[nlayer](1) << "," << pUnit[nlayer](2) << "); "
	   << "mom(p,theta,phi) = (" << InitMomentum.Mag()/global->GetMomentumUnit("GeV/c") << " GeV/c," 
	   << atan2(sqrt(pow(InitMomentum.X(),2) + pow(InitMomentum.Y(),2)),InitMomentum.Z())*180.0/TMath::Pi() << " deg,"
	   << atan2(InitMomentum.Y(),InitMomentum.X())*180.0/TMath::Pi() << " deg)"
          << endl;
    }
  }
  if(Myverbose) cout << endl;
  
  return;
  
}
//====================================================================
void  GTracker::GetUVDerWrtPosForAllLayers(std::vector<IntersectionHit_t>  ItersectionHitList,
					   double***                       DerUVwrtPos)
{
  
  const int NLayers(ItersectionHitList.size());
  for(int nlayer=0;nlayer<NLayers;nlayer++) {
    //Get the hit coordiates in the lab frame for hit nhit
    int geoElement_idx = ItersectionHitList[nlayer].geoElement_idx;
    double sn          = ItersectionHitList[nlayer].s;
    TVector3 pos_n     = Trajectory->GetTrueTrajectoryCoordinates(sn);
    
    GSurfaceObject*  aSurface = Geometry->GetGeometryElement(geoElement_idx)->GetMainSurface();
    
    for(int i=0;i<3;i++) {
      TVector3 DerGUVW_nhit_i = aSurface->GetDerUVWFromXYZ_Analytical(pos_n,i);
      
      for(int k=0;k<2;k++) DerUVwrtPos[nlayer][k][i] = DerGUVW_nhit_i(k);
    }
  }
  
  return;
  
}
//====================================================================
void  GTracker::GetHitUVCovMatrixWithMSFromGlobalCalc(std::vector<IntersectionHit_t> ItersectionHitList_global,
						      std::vector<IntersectionHit_t> ItersectionHitList_current,
						      TMatrixD*  Fr, TMatrixD* Fp, TMatrixD* Kr, TMatrixD* Kp,
						      double*    SigmaMS2, TVector3*  m1Unit, TVector3*  m2Unit,
						      double*    SigmaEloss2, TVector3*  pUnit,
						      double***  DerUVwrtPos,
						      TMatrixD&  HitUVCovMatrix)
{
  
  double    InitMomMag2 = Trajectory->GetInitMomentum().Mag2();
  double    EOp2        = (InitMomMag2 + pow(global->GetParticleMass(Trajectory->GetParticle()),2))/InitMomMag2;
  
  const int NLayers_global(ItersectionHitList_global.size()); // Number of intersection geometry layers by the particle in global geometry
  const int NLayers_current(ItersectionHitList_current.size()); // Number of intersection geometry layers by the particle in sub-geomtry
  
  std::vector<int> SensLayersList_current;
  SensLayersList_current.clear();
  for(int ihit=0;ihit<NLayers_current;ihit++) {
    if(ItersectionHitList_current[ihit].IsSensitivePoint) SensLayersList_current.push_back(ihit);
  }
  const int Nhits(SensLayersList_current.size()); //number of hits
  HitUVCovMatrix.ResizeTo(2*Nhits,2*Nhits);

  std::vector<int> global_layer;
  global_layer.clear();
  for(int nlayer=0;nlayer<NLayers_current;nlayer++) { //layer loop
    for(int klayer=0;klayer<NLayers_global;klayer++) {
      if(ItersectionHitList_global[klayer].geoElement_idx == ItersectionHitList_current[nlayer].geoElement_idx && 
	 ItersectionHitList_global[klayer].s              == ItersectionHitList_current[nlayer].s) {
	global_layer.push_back(klayer);
	break;
      }
    }
  }

  //Calculate matrix elements, F and K arrays, which contain geometry information as well as particle initial position and momentum.
  //There are later used to calculate the correlations between the hit coordinates induced by multiple scattering.
  //double F[3][3][NLayers_current][NLayers_current];    //^ik_nl
  //double K[3][3][NLayers_current][NLayers_current];    //^ik_nl
  if(NLayers_current > MaxNLayers) {
    cout << endl;
    cout << "GTracker::GetHitUVCovMatrixWithMSFromGlobalCalc:: Number of layers " << NLayers_current << " is higher than maximum allowed number " << MaxNLayers << endl;
    cout << endl;
    assert(false);
  }
  if(Nhits > MaxNhits) {
    cout << endl;
    cout << "GTracker::GetHitUVCovMatrixWithMSFromGlobalCalc:: Number of layers " << Nhits << " is higher than maximum allowed number " << MaxNhits << endl;
    cout << endl;
    assert(false);
  }
  
  for(int nlayer=0;nlayer<NLayers_current;nlayer++) {
    for(int llayer=0;llayer<NLayers_current;llayer++) {
      F[nlayer][llayer].Zero();
      K[nlayer][llayer].Zero();
    }
  }

  double Fmax[3][NLayers_current];
  for(int nlayer=0;nlayer<NLayers_current;nlayer++) {
    for(int i=0;i<3;i++) Fmax[i][nlayer] = -1.0e+20;
  }
  
  TMatrixD tmp_matrix;
  tmp_matrix.ResizeTo(3,3);
  
#if 0
  for(int ihit=0;ihit<int(ItersectionHitList_global.size());ihit++) {
    double   s_tmp   = ItersectionHitList_global[ihit].s;
    TVector3 pos_tmp = Trajectory->GetTrueTrajectoryCoordinates(s_tmp);
    TVector3 mom_tmp = Trajectory->GetTrueTrajectoryMon(s_tmp);
    s_tmp   /= global->GetUnit("cm");
    pos_tmp *= (1.0/global->GetUnit("cm"));
    mom_tmp *= (1.0/global->GetUnit("GeV/c"));
    cout << "layer = " << ihit << ", s = " << s_tmp << " cm, pos = (" 
         << pos_tmp(0) << "," << pos_tmp(1) << "," << pos_tmp(2) << ") cm; mom = ("
	 << mom_tmp(0) << "," << mom_tmp(1) << "," << mom_tmp(2) << ") GeV/c; "
         << endl;
	 
    cout << "Fr matrix:" << endl;
    tmp_matrix = Fr[ihit];
    tmp_matrix.Print();
    cout << endl;
    
    cout << "Fp matrix (cm/GeV/c):" << endl;
    tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*Fp[ihit];
    tmp_matrix.Print();
    cout << endl;
    
    cout << "Kr matrix (GeV/c/cm):" << endl;
    tmp_matrix = (global->GetUnit("cm")/global->GetUnit("GeV/c"))*Kr[ihit];
    tmp_matrix.Print();
    cout << endl;
    
    cout << "Kp matrix:" << endl;
    tmp_matrix = Kp[ihit];
    tmp_matrix.Print();
    cout << endl;
  }
#endif
  
  //Filling up of the F and K arrays
  for(int nlayer=0;nlayer<NLayers_current;nlayer++) { //layer loop
    if(nlayer == 0) continue;
    
    int nlayer_global = global_layer[nlayer];
    
    for(int llayer=0;llayer<nlayer;llayer++) { //loop on layers previous to this layer
      //cout << "Begin nlayer = " << nlayer << ", llayer = " << llayer << endl;
      if(llayer < nlayer-1) {
	F[nlayer][llayer] = Fr[nlayer_global]*F[nlayer-1][llayer] + Fp[nlayer_global]*K[nlayer-1][llayer];
	K[nlayer][llayer] = Kr[nlayer_global]*F[nlayer-1][llayer] + Kp[nlayer_global]*K[nlayer-1][llayer];
	
#if 0
	cout << endl;
	cout << "F matrix:" << endl;
	cout << "Fr[" << nlayer_global << "] matrix:" << endl;
	tmp_matrix = Fr[nlayer_global];
        tmp_matrix.Print();
	cout << "F[" << nlayer-1 << "][" << llayer << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*F[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "Fr[" << nlayer_global << "] x F[" << nlayer-1 << "][" << llayer << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*Fr[nlayer_global]*F[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "Fp[" << nlayer_global << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*Fp[nlayer_global];
        tmp_matrix.Print();
	cout << "K[" << nlayer-1 << "][" << llayer << "]:" << endl;
	tmp_matrix = K[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "Fp[" << nlayer_global << "] x K[" << nlayer-1 << "][" << llayer << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*Fp[nlayer_global]*K[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "F[" << nlayer << "][" << llayer << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*F[nlayer][llayer];
        tmp_matrix.Print();
	
	cout << endl;
	
	cout << "K matrix:" << endl;
	cout << "Kr[" << nlayer_global << "] matrix:" << endl;
	tmp_matrix = (global->GetUnit("cm")/global->GetUnit("GeV/c"))*Kr[nlayer_global];
        tmp_matrix.Print();
	cout << "F[" << nlayer-1 << "][" << llayer << "]:" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*F[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "Kr[" << nlayer_global << "] x F[" << nlayer-1 << "][" << llayer << "] matrix:" << endl;
	tmp_matrix = Kr[nlayer_global]*F[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "Kp[" << nlayer_global << "] matrix:" << endl;
	tmp_matrix = Kp[nlayer_global];
        tmp_matrix.Print();
	cout << "K[" << nlayer-1 << "][" << llayer << "]:" << endl;
	tmp_matrix = K[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "Kp[" << nlayer_global << "] x K[" << nlayer-1 << "][" << llayer << "] matrix:" << endl;
	tmp_matrix = Kp[nlayer_global]*K[nlayer-1][llayer];
        tmp_matrix.Print();
	cout << "K[" << nlayer << "][" << llayer << "] matrix:" << endl;
	tmp_matrix = K[nlayer][llayer];
        tmp_matrix.Print();
	cout << endl;
#endif
      }
      else {
	F[nlayer][llayer] = Fp[nlayer_global];
	K[nlayer][llayer] = Kp[nlayer_global];
#if 0
	cout << endl;
	cout << "F matrix:" << endl;
	cout << "Fp[" << nlayer_global << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*Fp[nlayer_global];
        tmp_matrix.Print();
	cout << "F[" << nlayer << "][" << llayer << "] matrix (cm/GeV/c):" << endl;
	tmp_matrix = (global->GetUnit("GeV/c")/global->GetUnit("cm"))*F[nlayer][llayer];
        tmp_matrix.Print();
	
	cout << endl;
	
	cout << "K matrix:" << endl;
	cout << "Kp[" << nlayer_global << "] matrix:" << endl;
	tmp_matrix = Kp[nlayer_global];
        tmp_matrix.Print();
	cout << "K[" << nlayer << "][" << llayer << "] matrix:" << endl;
	tmp_matrix = K[nlayer][llayer];
        tmp_matrix.Print();
	cout << endl;
#endif
      }
      for(int i=0;i<3;i++) {  //coordinate loop
	for(int j=0;j<3;j++) { //derivative loop
	  if(Fmax[i][nlayer] < TMath::Abs(F[nlayer][llayer](i,j))) Fmax[i][nlayer] = TMath::Abs(F[nlayer][llayer](i,j));
	}
      }
      //cout << "End   nlayer = " << nlayer << ", llayer = " << llayer << endl;
      //cout << endl;
    } //end of loop on layers previous to this layer
  } // end of layer loop
  
  //Now it is calculated the covariance matrix of the intersections of the 
  //track with the geometry
  //Now the covariance matrix is filled-up
  //double Cov_deltaR_deltaR[3][3][NLayers_current][NLayers_current];  
  
  for(int nhit=0;nhit<Nhits;nhit++) { //1st loop on layers
    int nlayer = SensLayersList_current[nhit];
    
    //Get the hit coordiates in the lab frame for hit nhit
    int geoElement_idx = ItersectionHitList_current[nlayer].geoElement_idx;
    double sn          = ItersectionHitList_current[nlayer].s;
    TVector3 pos_n     = Trajectory->GetTrueTrajectoryCoordinates(sn);
    TVector3 mom_n     = Trajectory->GetTrueTrajectoryMon(sn);
    if(UseMonOrigin) mom_n = Trajectory->GetInitMomentum();
    
    int nlayer_global  = global_layer[nlayer];
    
    //Get the sensitive plane U and V spatial resolutions    
    //Still need to specify the momentum and coordinates of the hit to the function
    TVector2 Resolution = Geometry->GetResolutionUV(geoElement_idx,mom_n,pos_n);
    double Sigma2_sp_Un = pow(Resolution.X(),2);
    double Sigma2_sp_Vn = pow(Resolution.Y(),2);
    
    for(int mhit=0;mhit<Nhits;mhit++) { //1st loop on layers
      int mlayer = SensLayersList_current[mhit];
      
      if(nhit > mhit) continue;
      
      int mlayer_global = global_layer[mlayer];
      
      CovUU[nhit][mhit] = global->delta_cronequer(nhit,mhit)*Sigma2_sp_Un;
      CovVV[nhit][mhit] = global->delta_cronequer(nhit,mhit)*Sigma2_sp_Vn;
      CovUV[nhit][mhit] = 0.0;
      if(nhit != mhit) CovUV[mhit][nhit] = 0.0;
      
      for(int i=0;i<3;i++) {  //1st loop on coordinate
	for(int j=0;j<3;j++) { //2nd loop on coordinate
	  
	  Cov_deltaR_deltaR[i][j][nlayer][mlayer] = 0.0;

	  for(int llayer=0;llayer<nlayer;llayer++) {
	    int llayer_global = global_layer[llayer];
	    
	    //Calculate only those contributions to MS which are significant, i.e. higher than MinMSEffect
	    double MyminEffect = Fmax[i][nlayer]*Fmax[j][mlayer]*SigmaMS2[llayer_global]*InitMomMag2;
	    if(IncludeEloss) MyminEffect += Fmax[i][nlayer]*Fmax[j][mlayer]*SigmaEloss2[llayer_global]*EOp2;
	    
	    if(sqrt(MyminEffect) < MinMSEffect) continue;
	    
	    for(int alayer=0;alayer<mlayer;alayer++) {
	      int alayer_global = global_layer[alayer];
	    
	      if(llayer_global != alayer_global) continue;
	      
	      for(int k=0;k<3;k++) {
                for(int b=0;b<3;b++) {
		  //The covariance matrix is made of a sum of elements with geometry factors
		  double quantity  = F[nlayer][llayer](i,k)*F[mlayer][alayer](j,b);
		  
		  // ... and multiple by multiple scattering factors
		  double quantity_MS  = (m1Unit[llayer_global](k)*m1Unit[alayer_global](b) + m2Unit[llayer_global](k)*m2Unit[alayer_global](b))*SigmaMS2[llayer_global];
		  quantity_MS        *= InitMomMag2;
		  
		  double quantity_Eloss = 0.0;
		  if(IncludeEloss) {
		    quantity_Eloss  = (pUnit[llayer_global](k)*pUnit[alayer_global](b))*SigmaEloss2[llayer_global];
		    quantity_Eloss *= EOp2;
		  }
		  
		  Cov_deltaR_deltaR[i][j][nlayer][mlayer] += quantity*quantity_MS;
		  if(IncludeEloss)  Cov_deltaR_deltaR[i][j][nlayer][mlayer] += quantity*quantity_Eloss;
		  
	        }
	      }
	    }

	  }
	}  // end of 2nd loop on coordinate
      } // end of 1st loop on coordinate
      
      for(int i=0;i<3;i++) {  //1st loop on coordinate
        for(int j=0;j<3;j++) { //2nd loop on coordinate
	  CovUU[nhit][mhit] += DerUVwrtPos[nlayer_global][0][i]*DerUVwrtPos[mlayer_global][0][j]*Cov_deltaR_deltaR[i][j][nlayer][mlayer];
	  CovVV[nhit][mhit] += DerUVwrtPos[nlayer_global][1][i]*DerUVwrtPos[mlayer_global][1][j]*Cov_deltaR_deltaR[i][j][nlayer][mlayer];
	  
	  
	  CovUV[nhit][mhit] += DerUVwrtPos[nlayer_global][0][i]*DerUVwrtPos[mlayer_global][1][j]*Cov_deltaR_deltaR[i][j][nlayer][mlayer];
	  if(nhit != mhit) CovUV[mhit][nhit] += DerUVwrtPos[mlayer_global][0][i]*DerUVwrtPos[nlayer_global][1][j]*Cov_deltaR_deltaR[j][i][nlayer][mlayer];
	} //end of 2nd loop on coordinate
      } // end of 1st loop on coordinate
     
      CovUU[mhit][nhit] = CovUU[nhit][mhit];
      CovVV[mhit][nhit] = CovVV[nhit][mhit];
     
    } //end of end loop on layers
  } //end of 1st loop on layers
  
  //Now put the covariances between U and V in a single matrix
  for(int i=0;i<2*Nhits;i++) {
    for(int j=0;j<2*Nhits;j++) {
      HitUVCovMatrix(i,j) = 0.0;
      
      if(i%2 == 0 && j%2 == 0)      HitUVCovMatrix(i,j) = CovUU[i/2][j/2];
      else if(i%2 != 0 && j%2 != 0) HitUVCovMatrix(i,j) = CovVV[(i-1)/2][(j-1)/2];
      else if(i%2 == 0 && j%2 != 0) HitUVCovMatrix(i,j) = CovUV[i/2][(j-1)/2];
      else if(i%2 != 0 && j%2 == 0) HitUVCovMatrix(i,j) = CovUV[j/2][(i-1)/2];
    }
  }

  return;
  
}
//====================================================================
void  GTracker::GetHitUVCovMatrixWithMS(std::vector<IntersectionHit_t> ItersectionHitList,
					TMatrixD&                      HitUVCovMatrix)
{
  
  const int NLayers(ItersectionHitList.size()); // Number of intersection geometry layers by the particle

  if(NLayers > MaxNLayers) {
    cout << endl;
    cout << "GTracker::GetHitUVCovMatrixWithMS:: Number of layers " << NLayers << " is higher than maximum allowed number " << MaxNLayers << endl;
    cout << endl;
    assert(false);
  }
  
  //initialize some arrays
  for(int nlayer=0;nlayer<NLayers;nlayer++) {
    SigmaMS2_local[nlayer] = 0.0;
    m1Unit_local[nlayer] = TVector3(0,0,0);
    m2Unit_local[nlayer] = TVector3(0,0,0);
    
    SigmaEloss2_local[nlayer] = 0.0;
    pUnit_local[nlayer]       = TVector3(0,0,0);
    
    Fr_local[nlayer].Zero();
    Fp_local[nlayer].Zero();
    Kr_local[nlayer].Zero();
    Kp_local[nlayer].Zero();
    
    for(int i=0;i<3;i++) { 
      for(int iuv=0;iuv<2;iuv++) {
        DerUVwrtPos_local[nlayer][iuv][i] = 0.0;
      }
      
    }
  }

  //Getting the values of the Fr,Fp,Kr,Kp arrays which contain the derivatives of position 
  //and momentum at layer nLayer w.r.t. position and momentum at previous layer
  GetFsAndKsMS(ItersectionHitList,Fr_local,Fp_local,Kr_local,Kp_local);
  
  //At this stage it is calculated the multiple scattering angle at each geometry intersection.
  //It is also calculated a couple of othogonal unit vectors normal to the particles direction at 
  //the given geometry intersection.
  GetSigmaMSAndPerpVectors(ItersectionHitList,SigmaMS2_local,m1Unit_local,m2Unit_local);
  
  //At this stage it is calculated the Eloss resolution at each geometry intersection.
  //It is also calculated the unit vector pointing in the momentum direction  at the given geometry intersection.
  if(IncludeEloss) GetSigmaElossAndUnitVector(ItersectionHitList,SigmaEloss2_local,pUnit_local);

  //Calculate the derivatives of U and V w.r.t position X,Y,Z for all the layers
  GetUVDerWrtPosForAllLayers(ItersectionHitList,DerUVwrtPos_local);

  //Now calculate the covariace matrix between the U and V coordiates of the sensitive planes crossed by the particle
  GetHitUVCovMatrixWithMSFromGlobalCalc(ItersectionHitList,ItersectionHitList,
					Fr_local,Fp_local,Kr_local,Kp_local,
					SigmaMS2_local,m1Unit_local,m2Unit_local,
					SigmaEloss2_local,pUnit_local,
					DerUVwrtPos_local,
					HitUVCovMatrix);
  return;
  
}
//====================================================================
void  GTracker::ConfigureTelescopePlanes(std::vector<IntersectionHit_t>&  ItersectionHitList)
{
  
  if(Geometry->GetNTelescopePlanes() == 0) return;
  
  //Modify intersection list to set as sensitive only the geo-elements declared as telescope planes
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(!ItersectionHitList[ihit].IsSensitivePoint) continue;
    
    bool IsInTelescopeList = false;
    for(int itel=0;itel<Geometry->GetNTelescopePlanes();itel++) {
      if(Geometry->GetTelescopePlaneIndexInGeometry(itel) == ItersectionHitList[ihit].geoElement_idx) {
	IsInTelescopeList = true;
	break;
      }
    }
    if(!IsInTelescopeList) ItersectionHitList[ihit].IsSensitivePoint = false;
    
  }
  
  return;
  
}
//====================================================================
void  GTracker::GetTelescopeResolutionAtDUTPlanes(std::vector<IntersectionHit_t> ItersectionHitList,
						  TMatrixD  FitCovMatrix,
						  std::vector<TelResolAtDUT_t>&  TelResolAtDUTList)
{
  
  //get telescope resolution at DUT planes
  
  TelResolAtDUTList.clear();
  
  if(Geometry->GetNDUTPlanes() == 0) return;
  
  //s-range of the telescope planes
  double Rs_Tel[2];
  Rs_Tel[0] = +1.0e+20;
  Rs_Tel[1] = -1.0e+20;
  
  //Get 1st list of DUT planes
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    //Get s-range of telescope planes
    bool IsTelescope = false;
    for(int itel=0;itel<Geometry->GetNTelescopePlanes();itel++) {
      if(Geometry->GetTelescopePlaneIndexInGeometry(itel) == ItersectionHitList[ihit].geoElement_idx) {
	if(Rs_Tel[0] > ItersectionHitList[ihit].s)  Rs_Tel[0] = ItersectionHitList[ihit].s;
	if(Rs_Tel[1] < ItersectionHitList[ihit].s)  Rs_Tel[1] = ItersectionHitList[ihit].s;
	IsTelescope = true;
	break;
      }
    }
    if(IsTelescope) continue;
    
    bool IsDUT = false;
    for(int itel=0;itel<Geometry->GetNDUTPlanes();itel++) {
      if(Geometry->GetDUTPlaneIndexInGeometry(itel) == ItersectionHitList[ihit].geoElement_idx) {
	IsDUT = true;
	break;
      }
    }
    if(IsDUT) {
      TelResolAtDUT_t aTelResolDUT;
      aTelResolDUT.intersection_idx = ihit;
      aTelResolDUT.geoElement_idx   = ItersectionHitList[ihit].geoElement_idx;
      aTelResolDUT.s                = ItersectionHitList[ihit].s;      
      aTelResolDUT.TelResolU        = -999.9;
      aTelResolDUT.TelResolV        = -999.9;
      aTelResolDUT.TelCorrUV        = -999.9;
      aTelResolDUT.Shadow_area      = -999.9;
      aTelResolDUT.Nbkg             = -999.9;
    
      TelResolAtDUTList.push_back(aTelResolDUT);
    }
    
  }
  
  if(TelResolAtDUTList.size() == 0) return;
  
  const int Npars(Trajectory->GetNParameters());
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    for(int ipar=0;ipar<Npars;ipar++) GammaUDer_local[ihit][ipar] = GammaVDer_local[ihit][ipar] = 0.0;
  }
  GetUVDerivatiesWrtPars(ItersectionHitList,GammaUDer_local,GammaVDer_local,Npars);
  
  //initialize some arrays
  for(int nlayer=0;nlayer<int(ItersectionHitList.size());nlayer++) {
    SigmaMS2_local[nlayer] = 0.0;
    m1Unit_local[nlayer] = TVector3(0,0,0);
    m2Unit_local[nlayer] = TVector3(0,0,0);
  }
  GetSigmaMSAndPerpVectors(ItersectionHitList,SigmaMS2_local,m1Unit_local,m2Unit_local);
  if(IncludeEloss) GetSigmaElossAndUnitVector(ItersectionHitList,SigmaEloss2_local,pUnit_local);
  
  std::vector<int>  DUTPlanes_below_Tel;
  std::vector<int>  DUTPlanes_within_Tel;
  std::vector<int>  DUTPlanes_above_Tel;
  DUTPlanes_below_Tel.clear();
  DUTPlanes_within_Tel.clear();
  DUTPlanes_above_Tel.clear();
  for(int idut=0;idut<int(TelResolAtDUTList.size());idut++) {
    if(TelResolAtDUTList[idut].s < Rs_Tel[0])      DUTPlanes_below_Tel.push_back(idut);
    else if(TelResolAtDUTList[idut].s > Rs_Tel[1]) DUTPlanes_above_Tel.push_back(idut);
    else                                           DUTPlanes_within_Tel.push_back(idut);
  }
  
  //Fill 1st the easiest ones => the DUT resolutions within the telescope
  for(int idut=0;idut<int(DUTPlanes_within_Tel.size());idut++) {
    int dut_idx = DUTPlanes_within_Tel[idut];
    int ihit    = TelResolAtDUTList[dut_idx].intersection_idx;
    TMatrixD CovUV;
    GetCovMatrixUVOfTrackIntersection(ItersectionHitList,ItersectionHitList,
				      ihit,
				      GammaUDer_local,GammaVDer_local,FitCovMatrix,CovUV);
    
    TelResolAtDUTList[dut_idx].TelResolU   = sqrt(CovUV(0,0));
    TelResolAtDUTList[dut_idx].TelResolV   = sqrt(CovUV(1,1));
    TelResolAtDUTList[dut_idx].TelCorrUV   = CovUV(0,1)/sqrt(CovUV(0,0)*CovUV(1,1));
    TelResolAtDUTList[dut_idx].Shadow_area = global->GetEllipseArea(TelResolAtDUTList[dut_idx].TelResolU,
								    TelResolAtDUTList[dut_idx].TelResolV,
								    TelResolAtDUTList[dut_idx].TelCorrUV);
    
    double Nbkg = 0.0;
    if(Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetIsSensitive()) {
      Nbkg  = Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetBkgRate();
      Nbkg *= Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetROtime();
      Nbkg *= TelResolAtDUTList[dut_idx].Shadow_area;
    }
    TelResolAtDUTList[dut_idx].Nbkg = Nbkg;
    
  }

  
  if(DUTPlanes_below_Tel.size() + DUTPlanes_above_Tel.size() == 0) return;
  GetAllTrackParamDerWrtMS(ItersectionHitList,SigmaMS2_local,m1Unit_local,m2Unit_local,Der_Params_prime_wrt_thetaMS_global);
  if(IncludeEloss) GetAllTrackParamDerWrtEloss(ItersectionHitList,SigmaEloss2_local,pUnit_local,Der_Params_prime_wrt_Eloss_global);
  
  int counter_dut = 0;

  //Fill the DUT resolutions below the telescope
  counter_dut = 0;
  TMatrixD  FitCovMatrix_below(FitCovMatrix);
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    int hit_idx = ItersectionHitList.size() - ihit - 1;
    if(ItersectionHitList[hit_idx].s > Rs_Tel[0]) continue;
    
    bool IsDUTbelow = false;
    int  dut_idx = -1;
    for(int idut=0;idut<int(DUTPlanes_below_Tel.size());idut++) {
      if(TelResolAtDUTList[DUTPlanes_below_Tel[idut]].intersection_idx == hit_idx) {
	IsDUTbelow = true;
	dut_idx = DUTPlanes_below_Tel[idut];
	break;
      }
    }
    
    if(IsDUTbelow) {
      counter_dut++;
      TMatrixD CovUV;
      GetCovMatrixUVOfTrackIntersection(ItersectionHitList,ItersectionHitList,
					hit_idx,
					GammaUDer_local,GammaVDer_local,FitCovMatrix_below,CovUV);
    
      TelResolAtDUTList[dut_idx].TelResolU = sqrt(CovUV(0,0));
      TelResolAtDUTList[dut_idx].TelResolV = sqrt(CovUV(1,1));
      TelResolAtDUTList[dut_idx].TelCorrUV = CovUV(0,1)/sqrt(CovUV(0,0)*CovUV(1,1));
      TelResolAtDUTList[dut_idx].Shadow_area = global->GetEllipseArea(TelResolAtDUTList[dut_idx].TelResolU,
								      TelResolAtDUTList[dut_idx].TelResolV,
								      TelResolAtDUTList[dut_idx].TelCorrUV);
    
      double Nbkg = 0.0;
      if(Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetIsSensitive()) {
        Nbkg  = Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetBkgRate();
        Nbkg *= Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetROtime();
        Nbkg *= TelResolAtDUTList[dut_idx].Shadow_area;
      }
      TelResolAtDUTList[dut_idx].Nbkg = Nbkg;
    }
    
    if(counter_dut == int(DUTPlanes_below_Tel.size())) continue;
    
    TMatrixD FitCovMatrix_current(FitCovMatrix_below);
    GetNewParCovarianceMatrix(ItersectionHitList,ItersectionHitList,
			      hit_idx,
			      Der_Params_prime_wrt_thetaMS_global,SigmaMS2_local,
			      Der_Params_prime_wrt_Eloss_global,SigmaEloss2_local,
			      FitCovMatrix_below,FitCovMatrix_current);
    FitCovMatrix_below = FitCovMatrix_current;
  }
  
  //Fill the DUT resolutions above the telescope
  counter_dut = 0;
  TMatrixD  FitCovMatrix_above(FitCovMatrix);
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    int hit_idx = ihit;
    if(ItersectionHitList[hit_idx].s < Rs_Tel[1]) continue;
    
    bool IsDUTabove = false;
    int  dut_idx = -1;
    for(int idut=0;idut<int(DUTPlanes_above_Tel.size());idut++) {
      if(TelResolAtDUTList[DUTPlanes_above_Tel[idut]].intersection_idx == hit_idx) {
	IsDUTabove = true;
	dut_idx = DUTPlanes_above_Tel[idut];
	break;
      }
    }
    
    if(IsDUTabove) {
      counter_dut++;
      TMatrixD CovUV;
      GetCovMatrixUVOfTrackIntersection(ItersectionHitList,ItersectionHitList,
					hit_idx,
					GammaUDer_local,GammaVDer_local,FitCovMatrix_above,CovUV);
    
      TelResolAtDUTList[dut_idx].TelResolU = sqrt(CovUV(0,0));
      TelResolAtDUTList[dut_idx].TelResolV = sqrt(CovUV(1,1));
      TelResolAtDUTList[dut_idx].TelCorrUV = CovUV(0,1)/sqrt(CovUV(0,0)*CovUV(1,1));
      TelResolAtDUTList[dut_idx].Shadow_area = global->GetEllipseArea(TelResolAtDUTList[dut_idx].TelResolU,
								      TelResolAtDUTList[dut_idx].TelResolV,
								      TelResolAtDUTList[dut_idx].TelCorrUV);
      
      double Nbkg = 0.0;
      if(Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetIsSensitive()) {
        Nbkg  = Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetBkgRate();
        Nbkg *= Geometry->GetGeometryElement(TelResolAtDUTList[dut_idx].geoElement_idx)->GetROtime();
        Nbkg *= TelResolAtDUTList[dut_idx].Shadow_area;
      }
      TelResolAtDUTList[dut_idx].Nbkg = Nbkg;
    }
    
    if(counter_dut == int(DUTPlanes_above_Tel.size())) continue;
    
    TMatrixD FitCovMatrix_current(FitCovMatrix_above);
    GetNewParCovarianceMatrix(ItersectionHitList,ItersectionHitList,
			      hit_idx,
			      Der_Params_prime_wrt_thetaMS_global,SigmaMS2_local,
			      Der_Params_prime_wrt_Eloss_global,SigmaEloss2_local,
			      FitCovMatrix_above,FitCovMatrix_current);
    FitCovMatrix_above = FitCovMatrix_current;
  }
  
  return;
  
}
//====================================================================
bool  GTracker::GetFitTrackParsCovMatrix(int        NhitsMin,
					 TMatrixD&  FitCovMatrix)
{
  
  //This function returns the track parameters and the corresponding covariance matrix
  //The inputs are
  // - A geometry
  // - A particle: nature, initial position and momentum
  // - The magntic field
  // - A reference point for trajectory parameterezation

  //Obtain first the intersections of the track with the different materials
  std::vector<IntersectionHit_t> ItersectionHitList;
  ItersectionHitList.clear();
  GetIntersectionsWithGeometry(ItersectionHitList);

  //Then return the track parameters and covariance matrix
  bool Status = true;
  Status = GetFitTrackParsCovMatrix_FromPlaneList(ItersectionHitList,NhitsMin,FitCovMatrix);
  
  return  Status;
    
}
//====================================================================
bool  GTracker::GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc(std::vector<IntersectionHit_t>  ItersectionHitList_global,
								      std::vector<IntersectionHit_t>  ItersectionHitList_current,
								      double**                        GammaUDer,
								      double**                        GammaVDer,
								      TMatrixD                        MeasCovMatrix,
								      TMatrixD&                       FitCovMatrix)
{
  
  const int NLayers_global(ItersectionHitList_global.size()); // Number of intersection geometry layers by the particle in global geometry
  const int NLayers_current(ItersectionHitList_current.size()); // Number of intersection geometry layers by the particle in sub-geomtry

  std::vector<int> SensLayersList_current;
  SensLayersList_current.clear();
  for(int ihit=0;ihit<NLayers_current;ihit++) {
    if(ItersectionHitList_current[ihit].IsSensitivePoint) SensLayersList_current.push_back(ihit);
  }
  std::vector<int> SensLayersList_global;
  SensLayersList_global.clear();
  for(int ihit=0;ihit<NLayers_global;ihit++) {
    if(ItersectionHitList_global[ihit].IsSensitivePoint) SensLayersList_global.push_back(ihit);
  }
  const int Nhits(SensLayersList_current.size()); //number of hits
  const int Npars(Trajectory->GetNParameters());  //number of track parameters
  
  std::vector<int> sens_global_layer;
  sens_global_layer.clear();
  for(int k_current=0;k_current<int(SensLayersList_current.size());k_current++) {
    for(int ihit_global=0;ihit_global<NLayers_global;ihit_global++) {
      if(ItersectionHitList_current[SensLayersList_current[k_current]].geoElement_idx == ItersectionHitList_global[ihit_global].geoElement_idx && 
	 ItersectionHitList_current[SensLayersList_current[k_current]].s              == ItersectionHitList_global[ihit_global].s) {
	sens_global_layer.push_back(ihit_global);
	break;
      }
    }
  }

#if 0
  for(int ihit=0;ihit<NLayers_global;ihit++) {
    cout << endl;
    for(int ipar=0;ipar<Npars;ipar++) {
      cout << "Gamma[layer = " << ihit << "](U,V)DerPar-" << Trajectory->GetParameterName(ipar) << " = (" 
           << GammaUDer[ihit][ipar]*(global->GetUnit(Trajectory->GetParameterUnit(ipar))/global->GetUnit("mm")) << "," 
	   << GammaVDer[ihit][ipar]*(global->GetUnit(Trajectory->GetParameterUnit(ipar))/global->GetUnit("mm")) << ")" 
	   << " mm/" << Trajectory->GetParameterUnitTitle(ipar).Data()
           << endl;
    }
    cout << endl;
  }
#endif
  
  double GammaDer[2*Nhits][Npars];
  for(int i=0;i<2*Nhits;i++) {
    for(int ipar=0;ipar<Npars;ipar++) {  //loop on track parameters
      if(i%2 == 0) {
	int index = i/2;
	int index_global = sens_global_layer[index];
	
	GammaDer[i][ipar] = GammaUDer[index_global][ipar];
      }
      else if(i%2 != 0) {
	int index = (i-1)/2;
	int index_global = sens_global_layer[index];

	GammaDer[i][ipar] = GammaVDer[index_global][ipar];
      }
    } //end of loop on track parameters
  }

  //Check that the measurement covariance matrix is symmetric
  for(int i=0;i<2*Nhits;i++) {
    for(int j=0;j<2*Nhits;j++) {
      if(i >= j) continue;
      if(TMath::Abs(MeasCovMatrix(j,i) - MeasCovMatrix(i,j))/pow(global->GetDistanceUnit("mm"),2) > 1.0e-8) {
	cout << "WARNNING in GTracker::GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc MeasCovMatrix matrix element (" << i << "," << j << ") is not equal to its transposed one. Difference is " 
	     << MeasCovMatrix(j,i) - MeasCovMatrix(i,j) 
	     << endl;
      }
      //MeasCovMatrix(j,i) = MeasCovMatrix(i,j);
    }
  }

  double epsilon_tol = 1.0e-20;
  
  //Invert the measurement covariane matrix
  TMatrixD InvMeasCovMatrix;
  InvMeasCovMatrix.ResizeTo(2*Nhits,2*Nhits);
  TDecompLU lu_meas(MeasCovMatrix);
  lu_meas.SetTol(epsilon_tol);
  int nr  = 0;
  while(!lu_meas.Invert(InvMeasCovMatrix) && nr < 10) {
    lu_meas.SetMatrix(MeasCovMatrix);
    lu_meas.SetTol(0.1*lu_meas.GetTol());
    nr++;
  }
  if(nr > 0) cout << "nr(meas) = " << nr << endl;
  
  if(verbose) {
    TMatrixD MeasCovUnitMatrix;
    MeasCovUnitMatrix.ResizeTo(2*Nhits,2*Nhits);
    for(int i=0;i<2*Nhits;i++) {
      for(int j=0;j<2*Nhits;j++) {
        MeasCovUnitMatrix(i,j) = 0.0;
        for(int m=0;m<2*Nhits;m++) MeasCovUnitMatrix(i,j) += MeasCovMatrix(i,m)*InvMeasCovMatrix(m,j);
      }
    }
  
    for(int i=0;i<Nhits;i++) cout << "Hit uncertainties(U,V) for hit " << i << " = (" 
                                  << sqrt(MeasCovMatrix(2*i,2*i))/global->GetDistanceUnit("um") << "," 
				  << sqrt(MeasCovMatrix(2*i+1,2*i+1))/global->GetDistanceUnit("um") << ") um" << endl;
    cout << endl;
    cout << "Measurements Cov matrix:" << endl;
    MeasCovMatrix.Print();
    cout << "Inverted Measurements Cov matrix:" << endl;
    InvMeasCovMatrix.Print();
    cout << endl;
    cout << "Unit matrix Measurements Cov matrix:" << endl;
    MeasCovUnitMatrix.Print();
    cout << endl;
  }

  //Now calculate the inverse cov matrix of the fit parameters
  TMatrixD InvFitCovMatrix;
  InvFitCovMatrix.ResizeTo(Npars,Npars);
  for(int ipar=0;ipar<Npars;ipar++) {  //1st loop on track parameters
    for(int jpar=0;jpar<Npars;jpar++) {  //2nd loop on track parameters
      InvFitCovMatrix(ipar,jpar) = 0.0;
      
      for(int nhit=0;nhit<2*Nhits;nhit++) {
	for(int mhit=0;mhit<2*Nhits;mhit++) {
	  InvFitCovMatrix(ipar,jpar) += GammaDer[nhit][ipar]*InvMeasCovMatrix(nhit,mhit)*GammaDer[mhit][jpar];
	}
      }
      
    }  //end of 2nd loop on track parameters
  }  //end of 1st loop on track parameters
  
  //Check that the inverse of the track parameters covariance matrix is symmetric
  for(int ipar=0;ipar<Npars;ipar++) {
    for(int jpar=0;jpar<Npars;jpar++) {
      if(ipar >= jpar) continue;
      //if(TMath::Abs(InvFitCovMatrix(jpar,ipar) - InvFitCovMatrix(ipar,jpar)) > 1.0e-8) {
      //  cout << "WARNNING:: InvFitCovMatrix matrix element (" << ipar << "," << jpar << ") is not equal to its transposed one. Difference is " << InvFitCovMatrix(jpar,ipar) - InvFitCovMatrix(ipar,jpar) << endl;
      //}
      InvFitCovMatrix(jpar,ipar) = InvFitCovMatrix(ipar,jpar);
    }
  }
  
  for(int ipar=0;ipar<Npars;ipar++) {
    for(int jpar=0;jpar<Npars;jpar++) {
      if(TMath::Abs(InvFitCovMatrix(ipar,jpar)) < 1.0e-10) InvFitCovMatrix(ipar,jpar) = 0.0;
    }
  }
  
  //Obtain the fit track parameters covariance matrix by inverting its inverse
  FitCovMatrix.ResizeTo(Npars,Npars);
  TDecompLU lu_fit(InvFitCovMatrix);
  lu_fit.SetTol(epsilon_tol);
  nr = 0;
  while(!lu_fit.Invert(FitCovMatrix) && nr < 10) {
    lu_fit.SetMatrix(InvFitCovMatrix);
    lu_fit.SetTol(0.1*lu_fit.GetTol());
    nr++;
  }
  if(nr > 0) cout << "nr(fit) = " << nr << endl;
  
  //if(verbose || true) {
  if(verbose) {
    TMatrixD FitCovUnitMatrix;
    FitCovUnitMatrix.ResizeTo(Npars,Npars);
    for(int i=0;i<Npars;i++) {
      for(int j=0;j<Npars;j++) {
        FitCovUnitMatrix(i,j) = 0.0;
        for(int m=0;m<Npars;m++) FitCovUnitMatrix(i,j) += FitCovMatrix(i,m)*InvFitCovMatrix(m,j);
      }
    }
  
    cout << "Fit parameters inverse covariance matrix:" << endl;
    InvFitCovMatrix.Print();
    
    cout << "Fit parameters covariance matrix times inverse:" << endl;
    FitCovUnitMatrix.Print();
    
    cout << "Fit parameters covariance matrix:" << endl;
    FitCovMatrix.Print();
    
    for(int ipar=0;ipar<Npars;ipar++) {
      cout << "Resolution parameter " << ipar << " = " << sqrt(FitCovMatrix(ipar,ipar)) << endl;
    }
    cout << endl;
  }

  
  return true;
  
}
//====================================================================
bool  GTracker::GetFitTrackParsCovMatrix_FromPlaneList(std::vector<IntersectionHit_t> ItersectionHitList,
						       int                  NhitsMin,
						       TMatrixD&            FitCovMatrix)
{
  
  //This function returns the track parameters and the corresponding covariance matrix
  //The inputs are
  // - A geometry
  // - A set of insertion points with the geometry
  // - A particle: nature, initial position and momentum
  // - The magntic field
  // - A reference point for trajectory parameterezation
  
  bool Status = true;

  //Get the sub-list of sensitive layers intersected by the particle
  std::vector<int> SensLayersList;
  SensLayersList.clear();
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    if(ItersectionHitList[ihit].IsSensitivePoint) SensLayersList.push_back(ihit);
  }
  const int Nhits(SensLayersList.size()); //number of hits
  
  int MyNhitsMin = NhitsMin;
  if(MyNhitsMin < 0) MyNhitsMin = Trajectory->GetMinHits();
  else if(MyNhitsMin < Trajectory->GetMinHits()) MyNhitsMin = Trajectory->GetMinHits();
  if(Nhits < MyNhitsMin) return false;

  //Get the particle track fitting parameters
  Trajectory->GetFitParsFromInitConds();
  const int Npars(Trajectory->GetNParameters()); //This is the number of track parameters: 4 for a stright-line and 5 for a helix
    
  if(Npars > MaxNpars) {
    cout << endl;
    cout << "GTracker::GetFitTrackParsCovMatrix_FromPlaneList:: Number of Npars " << Npars << " is higher than maximum allowed number " << MaxNpars << endl;
    cout << endl;
    assert(false);
  }

  //This matrices are the derivatives of the local hit coordinates (1st index) w.r.t the track parameters (second index)
  //Now filling-up the "hits" convariance matrix
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
    for(int ipar=0;ipar<Npars;ipar++) GammaUDer_local[ihit][ipar] = GammaVDer_local[ihit][ipar] = 0.0;
  }
  GetUVDerivatiesWrtPars(ItersectionHitList,GammaUDer_local,GammaVDer_local,Npars);
  
  //Calculate the U,V hits covariance matrix, including MS and layer single-point resolution
  TMatrixD MeasCovMatrix;
  GetHitUVCovMatrixWithMS(ItersectionHitList,MeasCovMatrix);

  Status = GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc(ItersectionHitList,ItersectionHitList,GammaUDer_local,GammaVDer_local,MeasCovMatrix,FitCovMatrix);

  return Status;
  
}
//====================================================================
bool  GTracker::doTrkResolAnalysis(int        NhitsMin,
				   bool       DoTelescopeAnalysis,
				   int&       Nhits, double& Material_budget, double &dist1stPointToRef,
				   TMatrixD&  FitCovMatrix,
				   std::vector<TelResolAtDUT_t>& TelResolAtDUTList)
{
  
  //This functions does the track resoltion analysis for given track initial conditions, x0 and p0
  
  bool Status = true;

  //Get Intersections with geometry
  std::vector<IntersectionHit_t> ItersectionHitList;
  ItersectionHitList.clear();
  GetIntersectionsWithGeometry(ItersectionHitList);

  std::vector<IntersectionHit_t> ItersectionHitList_all;
  ItersectionHitList_all.clear();
  for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) ItersectionHitList_all.push_back(ItersectionHitList[ihit]);

  if(DoTelescopeAnalysis) {
    //Configure Telescope planes in list of geo-elements intersections
    ConfigureTelescopePlanes(ItersectionHitList);
  }
  
  GetMaterialBudget_FromPlaneList(ItersectionHitList,Material_budget,Nhits,dist1stPointToRef);

  if(DoTelescopeAnalysis) {
    double Rs[2];
    Rs[0] = +1.0e+20;
    Rs[1] = -1.0e+20;
    for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
      if(!ItersectionHitList[ihit].IsSensitivePoint) continue;
      if(Rs[0] > ItersectionHitList[ihit].s) Rs[0] = ItersectionHitList[ihit].s;
      if(Rs[1] < ItersectionHitList[ihit].s) Rs[1] = ItersectionHitList[ihit].s;
    }
    std::vector<IntersectionHit_t> ItersectionHitList_withinTracker;
    ItersectionHitList_withinTracker.clear();
    for(int ihit=0;ihit<int(ItersectionHitList.size());ihit++) {
      if(ItersectionHitList[ihit].s < Rs[0] || ItersectionHitList[ihit].s > Rs[1]) continue;
      ItersectionHitList_withinTracker.push_back(ItersectionHitList[ihit]);
    }
    ItersectionHitList.clear();
    for(int ihit=0;ihit<int(ItersectionHitList_withinTracker.size());ihit++) ItersectionHitList.push_back(ItersectionHitList_withinTracker[ihit]);
    ItersectionHitList_withinTracker.clear();
  }
  
  Status = GetFitTrackParsCovMatrix_FromPlaneList(ItersectionHitList,NhitsMin,FitCovMatrix);
  
  TelResolAtDUTList.clear();
  if(DoTelescopeAnalysis) {
    //Get Telescope resolution at DUT's
    GetTelescopeResolutionAtDUTPlanes(ItersectionHitList_all,FitCovMatrix,TelResolAtDUTList);
  }

  return  Status;
  
}
//====================================================================
void  GTracker::GetFitTrackPseudoEfficiency(Efficiencies_t&      Efficiencies,
					    TMatrixD&            AveFitCovMatrix)
{
  
  //This functions calculates the track pseudo-efficiency for given track initial conditions, x0 and p0
  
  bool Myverbose = false;
  Myverbose = verbose;
  //Myverbose = true;
  
  //Initialization  of the efficiencies
  Efficiencies.Effic_tot          = 0.0;
  Efficiencies.Effic_NoFakes      = 0.0;
  Efficiencies.Effic_1Fake        = 0.0;
  Efficiencies.Effic_2orMoreFakes = 0.0;
  
  //Get track initial conditions
  TVector3 x0 = Trajectory->GetInitPosition();
  TVector3 p0 = Trajectory->GetInitMomentum();
  
  //Identify the track finder algorithm to use from initial conditions
  int trackfinder_index = Geometry->FindTrackFinderAlgo(x0,p0);
  //Return zero efficiency if no track finder is found
  if(trackfinder_index < 0) return;
  
  //Get the tracl finder algorithm
  GTrackFinderAlgo* aTrackFinderAlgo = Geometry->GetTrackFinderAlgo(trackfinder_index);
  
  //Call the corresponding function depending on the requested track finder algorithm
  if(aTrackFinderAlgo->GetType() == TString("FPCCDTrackFinder")) {
    //This algorithm only works for constant B-field over all space => Helix track
    if(Trajectory->GetType() != TString("Helix")) {
      cout << endl;
      cout << "Error in GTracker::GetFitTrackPseudoEfficiency: Geometry " << Geometry->GetName().Data() << "; for the FPCCDTrackFinder the trajectory has to be a Helix. Check your inputs. Exiting now!!!" << endl;
      cout << endl;
      assert(false);
    }
    GetFitTrackPseudoEfficiency_FPCCD(dynamic_cast<GTrackFinderAlgoFPCCD*>(aTrackFinderAlgo),Efficiencies,AveFitCovMatrix);
  }
  
  if(Myverbose) {
    cout << "Total       efficiency = " << Efficiencies.Effic_tot       << endl;
    cout << "NoFakes     efficiency = " << Efficiencies.Effic_NoFakes   << endl;
    cout << "1-Fake      efficiency = " << Efficiencies.Effic_1Fake << endl;
    cout << ">= 2 Fakes  efficiency = " << Efficiencies.Effic_2orMoreFakes << endl;
    cout << endl;
    cout << endl;
  }
  
  return;
  
}
//====================================================================
void  GTracker::GetFitTrackPseudoEfficiency_FPCCD(GTrackFinderAlgoFPCCD* aTrackFinderAlgo,
						  Efficiencies_t&      Efficiencies,
						  TMatrixD&            AveFitCovMatrix)

{

  //This functions calculates the track pseudo-efficiency using the so-called FPCCD track finder algorithm  for given track initial conditions, x0 and p0
  //See guariguanchi documentation for more details
  
  bool Myverbose = false;
  Myverbose = verbose;
  //Myverbose = true;
  
  double MinSeedProb = 1.0e-6;
  
  //Parameters for the minimum Pt cut applied to the track seeding
  TVector3 CenterPosition = aTrackFinderAlgo->GetCenterPosition();                               // Position of center 
  double PtCutMin         = aTrackFinderAlgo->GetPtMinCut();                                     // Minimum Pt
  TVector3 Bfield         = (dynamic_cast<GTrajectoryHelix*>(Trajectory))->GetConstantBField();  // B-field
  TVector3 UnitBfield     = Bfield.Unit();
  double RadiusMin        = PtCutMin/(global->betaCurvRadius*Bfield.Mag());                      // Minimum curvature radius
  
  double Chi2OndfSeedCut = aTrackFinderAlgo->GetChi2OndfSeedCut();
  
  if(Myverbose) {
    cout << endl;
    cout << "Start GTracker::GetFitTrackPseudoEfficiency_FPCCD  ";
    global->fWatch.Print();
    global->fWatch.Continue();
    
    TVector3 x0   = (1.0/global->GetDistanceUnit("cm"))   *Trajectory->GetInitPosition();
    TVector3 p0   = (1.0/global->GetMomentumUnit("GeV/c"))*Trajectory->GetInitMomentum();
    double   pt   = (p0 - (p0.Dot(UnitBfield))*UnitBfield).Mag();
    double   Rtrk = (1.0/global->GetDistanceUnit("cm"))*(global->GetMomentumUnit("GeV/c"))*pt/(global->betaCurvRadius*Bfield.Mag());
    
    cout << endl;
    cout << "Geometry " << Geometry->GetName().Data() << ", "
         << "x0 = (" << x0(0) << "," << x0(1) << "," << x0(2) << ") cm and p0 = (" << p0(0) << "," << p0(1) << "," << p0(2) << ") GeV/c; "
	 << "p = " << p0.Mag() << " GeV/c; "
	 << "pt = " << pt << " GeV/c; Rtrk = " << Rtrk << " cm"
         << endl;
    
    TVector3 CP    = (1.0/global->GetDistanceUnit("cm"))*CenterPosition;
    double   ptcut = (1.0/global->GetMomentumUnit("GeV/c"))*PtCutMin;
    double   r     = (1.0/global->GetDistanceUnit("cm"))*RadiusMin;
    cout << "CenterPosition = (" << CP(0) << "," << CP(1) << "," << CP(2) << ") cm; "
         << "Pt min cut = " << ptcut << " GeV/c" << "; "
         << "Bfield_unit = (" << UnitBfield(0) << "," << UnitBfield(1) << "," << UnitBfield(2) << "); "
	 << "Bfield_mag = " << Bfield.Mag()/global->GetUnit("T") << " T; "
	 << "RminCut    = " << r << " cm"
         << endl;
  }
  
  //Get the list of both sensitive and non-sensitive layers of the geometry intersected by the particle
  std::vector<IntersectionHit_t> ItersectionHitList_init; //List with the track intersections
  ItersectionHitList_init.clear();
  GetIntersectionsWithGeometry(ItersectionHitList_init);
  const int NLayers_init(ItersectionHitList_init.size());
  
  //Get the sub-list of sensitive layers intersected by the particle
  std::vector<int> SensLayersList_init;
  SensLayersList_init.clear();
  for(int ihit=0;ihit<NLayers_init;ihit++) {
    if(ItersectionHitList_init[ihit].IsSensitivePoint) SensLayersList_init.push_back(ihit);
  }
  const int Nhits_init(SensLayersList_init.size()); //number of hits, i.e. sensitive layers crossed by the particle
  
  if(Myverbose) {
    cout << endl;
    cout << "Original intersections (" << NLayers_init << "):" << endl;
    for(int ihit=0;ihit<NLayers_init;ihit++) {
      TVector3 pos = (1.0/global->GetDistanceUnit("cm"))*Trajectory->GetTrueTrajectoryCoordinates(ItersectionHitList_init[ihit].s);
      cout << "  ihit = " << ihit+1 << ": s = " << ItersectionHitList_init[ihit].s/global->GetDistanceUnit("cm") << " cm at "
           << "geo-element " <<  ItersectionHitList_init[ihit].geoElement_idx << ", " 
	   << "Name = " << Geometry->GetGeometryElement(ItersectionHitList_init[ihit].geoElement_idx)->GetName().Data() << ", "
	   << "LayerName = " << Geometry->GetGeometryElement(ItersectionHitList_init[ihit].geoElement_idx)->GetLayerName().Data() << ", "
	   << "SystemName = " << Geometry->GetGeometryElement(ItersectionHitList_init[ihit].geoElement_idx)->GetSystemName().Data() << "  ";
      if(ItersectionHitList_init[ihit].IsSensitivePoint) cout << "(Sensitive)";
      else                                               cout << "(Insensitive)";
      cout << "; pos(x,y,z,r) = (" << pos(0) << "," << pos(1) << "," << pos(2) << "," << sqrt(pow(pos(0),2) + pow(pos(1),2)) << ") cm";
      cout << endl;
    }
  }
  
  //Get the list of layers used in the track-finding algorithm
  std::vector<TString>  FullSystemList = aTrackFinderAlgo->GetFullSystemList();
  if(Myverbose) {
    cout << endl;
    cout << "List of systems for track finding:" << endl;
    for(int isys=0;isys<int(FullSystemList.size());isys++) {
      cout << "  isys = " << isys+1 << ": "
	   << "Name = " << FullSystemList[isys].Data()
           << endl;
    }
  }
  
  //Get the maximum and minimum path-lengths for the layers used in the track-finding algorithm
  double smin = +1.0e+3*global->GetUnit("m");
  double smax = -1.0e+3*global->GetUnit("m");
  for(int ihit=0;ihit<NLayers_init;ihit++) { //begin loop over all the intersected layers
    for(int ilayer=0;ilayer<int(FullSystemList.size());ilayer++) { //begin loop over all the layers used for the track-finding
      if(Geometry->GetGeometryElement(ItersectionHitList_init[ihit].geoElement_idx)->GetSystemName() == FullSystemList[ilayer]) {
	double s = ItersectionHitList_init[ihit].s;
	if(smin > s) smin = s;
	if(smax < s) smax = s;
	break;
      }
    } //end loop over all the layers used for the track-finding
  } //end loop over all the intersected layers
  
  if(Myverbose) {
    cout << endl;
    cout << "tracking systems elements limited between smin = " << smin/global->GetDistanceUnit("cm") << " cm and smax = " << smax/global->GetDistanceUnit("cm") << " cm" << endl;
  }
  
  //Now only consider the layers intersected by the track which are used for the track-finder algorithm 
  std::vector<IntersectionHit_t> ItersectionHitList;
  ItersectionHitList.clear();
  for(int ihit=0;ihit<int(ItersectionHitList_init.size());ihit++) { //begin loop over all the intersected layers
    double s = ItersectionHitList_init[ihit].s;
    if(s < smin || s > smax) continue;
    
    ItersectionHitList.push_back(ItersectionHitList_init[ihit]);
  } //end loop over all the intersected layers
  //ItersectionHitList_init.clear();
  
  //Number of layers used in track-finder algorithm crossed by the particle
  const int NLayers(ItersectionHitList.size());

  if(Myverbose) {
    cout << endl;
    cout << "Reduced intersections (" << NLayers << "):" << endl;
    for(int ihit=0;ihit<NLayers;ihit++) {
      TVector3 pos = (1.0/global->GetDistanceUnit("cm"))*Trajectory->GetTrueTrajectoryCoordinates(ItersectionHitList_init[ihit].s);
      cout << "  ihit = " << ihit+1 << ": s = " << ItersectionHitList[ihit].s/global->GetDistanceUnit("cm") << " cm at "
           << "geo-element " <<  ItersectionHitList[ihit].geoElement_idx << ", " 
	   << "Name = " << Geometry->GetGeometryElement(ItersectionHitList[ihit].geoElement_idx)->GetName().Data() << ", "
	   << "LayerName = " << Geometry->GetGeometryElement(ItersectionHitList[ihit].geoElement_idx)->GetLayerName().Data() << ", "
	   << "SystemName = " << Geometry->GetGeometryElement(ItersectionHitList[ihit].geoElement_idx)->GetSystemName().Data() << "  ";
      if(ItersectionHitList[ihit].IsSensitivePoint) cout << "(Sensitive)";
      else                                          cout << "(Insensitive)";
      cout << "; pos(x,y,z,r) = (" << pos(0) << "," << pos(1) << "," << pos(2) << "," << sqrt(pow(pos(0),2) + pow(pos(1),2)) << ") cm";
      cout << endl;
    }
  }
  
  //Get the sub-list of sensitive layers intersected by the particle
  std::vector<int> SensLayersList;
  SensLayersList.clear();
  for(int ihit=0;ihit<NLayers;ihit++) {
    if(ItersectionHitList[ihit].IsSensitivePoint) SensLayersList.push_back(ihit);
  }
  const int Nhits(SensLayersList.size()); //number of hits, i.e. sensitive layers crossed by the particle
  
  //Some cut parameters for the track-finding algorithm
  int Nhits_min           = aTrackFinderAlgo->GetNhitsMin();         // minimum number of hits for track reco
  double PurityMin_cut    = aTrackFinderAlgo->GetPurityMinCut();     // minimum purity => (# real track-hits)/(# all track-hits)
  int Nfakes_max          = int((1.0 - PurityMin_cut)*Nhits);        // maximum number of fake hits allowed
  int NfakesSeed_max      = aTrackFinderAlgo->GetNfakesMaxSeeding(); // maximum number of fakes in seed
  if(NfakesSeed_max > 1)          NfakesSeed_max = 1;
  if(NfakesSeed_max > Nfakes_max) NfakesSeed_max = 0;
  
  //Get the track number of parameters
  const int Npars(Trajectory->GetNParameters());
  AveFitCovMatrix.ResizeTo(Npars,Npars);
  
  //Maximum number of sensitive layers without a hit
  int Nmax_null_Layers = Nhits - Nhits_min;
  if(Myverbose) {
    cout << endl;
    cout << "Geometry Name = " << Geometry->GetName().Data() << endl;
    cout << "Tracking algorithm type = " << aTrackFinderAlgo->GetType().Data() << " and name = " << aTrackFinderAlgo->GetName().Data() << endl;
    cout << "Nhits = " << Nhits << ", Nhits_min = " << Nhits_min << ", Nfakes_max = " << Nfakes_max << ", NfakesSeed_max = " << NfakesSeed_max << ", "
         << "Nmax_null_Layers = " << Nmax_null_Layers << endl;
    cout << "Npars = " << Npars << endl;
  }
  
  if(Nhits < Nhits_min)                 return; // do nothing if number of hits is smaller than min number of requested hits for tracking
  if(Nhits < Trajectory->GetMinHits())  return; // do nothing if number of hits is smaller than min number of hits to reconstruct track
  
  if(NLayers > MaxNLayers) {
    cout << endl;
    cout << "GTracker::GetFitTrackPseudoEfficiency:: Number of layers " << NLayers << " is higher than maximum allowed number " << MaxNLayers << endl;
    cout << endl;
    assert(false);
  }
  if(Nhits > MaxNhits) {
    cout << endl;
    cout << "GTracker::GetFitTrackPseudoEfficiency:: Number of hits " << Nhits << " is higher than maximum allowed number " << MaxNhits << endl;
    cout << endl;
    assert(false);
  }
  if(Npars > MaxNpars) {
    cout << endl;
    cout << "GTracker::GetFitTrackPseudoEfficiency:: Number of hits " << Npars << " is higher than maximum allowed number " << MaxNpars << endl;
    cout << endl;
    assert(false);
  }

  //Get the list of seed layers configurations 
  std::vector<SeedLayers_t>  SeedLayerConfigs = aTrackFinderAlgo->GetSeedLayerConfigs();
  
  if(Myverbose) {
    cout << endl;
    cout << "Get seed configs actually touched by the particle   ";
    global->fWatch.Print();
    global->fWatch.Continue();
  }
  //Now get the seed configurations which the particle actually touches
  std::vector<int> GoodSeedConfigs;
  GoodSeedConfigs.clear();
  for(int iseedconf=0;iseedconf<int(SeedLayerConfigs.size());iseedconf++) {
    //Check if this seed config has at least the minimum requered number of hits for track resconstruction
    if(int(SeedLayerConfigs[iseedconf].SeedLayers.size()) < Trajectory->GetMinHits()) continue;
    
    bool MyGoodSeedConfig = true;
    
    std::vector<bool> GoodLayers;
    GoodLayers.clear();
    for(int ilayer=0;ilayer<int(SeedLayerConfigs[iseedconf].SeedLayers.size());ilayer++) GoodLayers.push_back(false);
    
    for(int ihit=0;ihit<Nhits;ihit++) { //begin loop over sensitive layers
      TString HitLayerName = Geometry->GetGeometryElement(ItersectionHitList[SensLayersList[ihit]].geoElement_idx)->GetLayerName();
      for(int ilayer=0;ilayer<int(SeedLayerConfigs[iseedconf].SeedLayers.size());ilayer++) {
	if(SeedLayerConfigs[iseedconf].SeedLayers[ilayer] == HitLayerName) {
	  GoodLayers[ilayer] = true;
	  break;
	}
      }
    } //end loop over sensitive layers
    
    for(int ilayer=0;ilayer<int(SeedLayerConfigs[iseedconf].SeedLayers.size());ilayer++) MyGoodSeedConfig &= GoodLayers[ilayer];
    
    if(MyGoodSeedConfig) {
      if(Myverbose) {
        cout << "  Seed config: ";
        for(int ilayer=0;ilayer<int(SeedLayerConfigs[iseedconf].SeedLayers.size());ilayer++) {
	  cout << SeedLayerConfigs[iseedconf].SeedLayers[ilayer].Data() << "  ";
        }
        cout << endl;
      }
      
      GoodSeedConfigs.push_back(iseedconf);
    }
  }
  if(GoodSeedConfigs.size() == 0) return;
  
  if(Myverbose) {
    cout << endl;
    cout << "List of good seed layers configuration (" << GoodSeedConfigs.size() << "):   ";
    global->fWatch.Print();
    global->fWatch.Continue();
    for(int i=0;i<int(GoodSeedConfigs.size());i++) {
      cout << "    ";
      for(int j=0;j<int(SeedLayerConfigs[GoodSeedConfigs[i]].SeedLayers.size());j++) {
	cout << SeedLayerConfigs[GoodSeedConfigs[i]].SeedLayers[j].Data();
	if(j+1 < int(SeedLayerConfigs[GoodSeedConfigs[i]].SeedLayers.size())) cout << "   ";
      }
      cout << endl;
    }
  }

  //Now get the minimum and maximum path-lengths for the seeding system
  double smin_seed = +1.0e+3*global->GetUnit("m");
  double smax_seed = -1.0e+3*global->GetUnit("m");
  for(int ihit=0;ihit<Nhits;ihit++) { //begin loop over sensitive layers
    double s = ItersectionHitList[SensLayersList[ihit]].s;
    TString HitLayerName = Geometry->GetGeometryElement(ItersectionHitList[SensLayersList[ihit]].geoElement_idx)->GetLayerName();
    for(int iseedconf=0;iseedconf<int(GoodSeedConfigs.size());iseedconf++) {
      for(int ilayer=0;ilayer<int(SeedLayerConfigs[GoodSeedConfigs[iseedconf]].SeedLayers.size());ilayer++) {
	if(SeedLayerConfigs[GoodSeedConfigs[iseedconf]].SeedLayers[ilayer] == HitLayerName) {
	  if(smin_seed > s) smin_seed = s;
	  if(smax_seed < s) smax_seed = s;
	}
      }
    }
  } //end loop over sensitive layers
  
  if(Myverbose) {
    cout << endl;
    cout << "Seeding system max and min path lengths, smin = " << smin_seed/global->GetUnit("cm") << " cm and smax = " << smax_seed/global->GetUnit("cm") << " cm" << endl;
  }
  
  if(Myverbose) {
    cout << endl;
    cout << "List of seed configurations :  ";
    global->fWatch.Print();
    global->fWatch.Continue();
  }
  
  //Now get the seed configurations and the corresponding base probabilities
  std::vector<long> SeedConfigList;
  SeedConfigList.clear();
  std::vector<LayersProbs_t> SeedConfigIntrinsicProbsList;
  SeedConfigIntrinsicProbsList.clear();
  std::vector<double> SeedConfigPtCutProbsList;
  SeedConfigPtCutProbsList.clear();
  std::vector<double> SeedConfigChi2CutProbsList;
  SeedConfigChi2CutProbsList.clear();
  std::vector<double> SeedConfigTotProbsList;
  SeedConfigTotProbsList.clear();
  for(int iseedconf=0;iseedconf<int(GoodSeedConfigs.size());iseedconf++) { //begin loop over seed configurations
    //First get the seed configuration in a long format
    long config = 0;
    int counter = 0;
    for(int ihit=0;ihit<Nhits;ihit++) { //begin loop over sensitive layers
      int index = Nhits - ihit - 1;
      double s    = ItersectionHitList[SensLayersList[index]].s;
      if(s < smin_seed || s > smax_seed) continue;
      
      int igeoElm = ItersectionHitList[SensLayersList[index]].geoElement_idx;
      TString HitLayerName = Geometry->GetGeometryElement(igeoElm)->GetLayerName();
      
      int ON = 0;
      for(int ilayer=0;ilayer<int(SeedLayerConfigs[GoodSeedConfigs[iseedconf]].SeedLayers.size());ilayer++) {
	if(SeedLayerConfigs[GoodSeedConfigs[iseedconf]].SeedLayers[ilayer] == HitLayerName) {
	  ON = 1;
	  break;
	}
      }
      
      config += ON*pow(10,counter);
      
      counter++;
    }
    SeedConfigList.push_back(config);
    
    if(Myverbose) {
      cout << "  config " << config << " : ";
      for(int ilayer=0;ilayer<int(SeedLayerConfigs[GoodSeedConfigs[iseedconf]].SeedLayers.size());ilayer++) {
	cout << SeedLayerConfigs[GoodSeedConfigs[iseedconf]].SeedLayers[ilayer].Data() << "   ";
      }
      global->fWatch.Print();
      global->fWatch.Continue();
    }
    
    //Then calculate the base probability using the sensor detector efficiency
    LayersProbs_t ListOfProbs;
    ListOfProbs.Probs.clear();
    std::vector<int>  SeedLayersTypeList;
    SeedLayersTypeList.clear();
    global->GetListOfLayersTypes(config,SeedLayersTypeList);
    ListOfProbs.Probs.resize(SeedLayersTypeList.size());
    for(int ihit=0;ihit<int(SeedLayersTypeList.size());ihit++) { //begin loop seed layers
      int index     = Nhits - ihit - 1;
      int the_layer = SeedLayersTypeList.size() - ihit - 1;
      
      double s    = ItersectionHitList[SensLayersList[index]].s;
      int igeoElm = ItersectionHitList[SensLayersList[index]].geoElement_idx;
      TVector3 pos = Trajectory->GetTrueTrajectoryCoordinates(s);
      TVector3 mom = Trajectory->GetTrueTrajectoryMon(s);
      double DetEffic = Geometry->GetEfficiency(igeoElm,mom,pos,Trajectory->GetParticle());
      
      ListOfProbs.Probs[the_layer] = TVector3(DetEffic,0.0,0.0);
    } //end loop over seed layers
    
    //Store the different components of seed probabilities
    SeedConfigIntrinsicProbsList.push_back(ListOfProbs);
    SeedConfigPtCutProbsList.push_back(0.0);
    SeedConfigChi2CutProbsList.push_back(0.0);
    SeedConfigTotProbsList.push_back(0.0);
  } //end loop over seed configurations
    
  if(Myverbose) {
    cout << endl;
    cout << "List of seed configurations (" << SeedConfigList.size() << "):   ";
    global->fWatch.Print();
    global->fWatch.Continue();
  }
  
  //Now calculate some global quantities to be used later
  for(int nlayer=0;nlayer<NLayers_init;nlayer++) {
    SigmaMS2_global[nlayer] = 0.0;
    m1Unit_global[nlayer] = TVector3(0,0,0);
    m2Unit_global[nlayer] = TVector3(0,0,0);
    
    SigmaEloss2_global[nlayer] = 0.0;
    pUnit_global[nlayer] = TVector3(0,0,0);
    
    Fr_global[nlayer].Zero();
    Fp_global[nlayer].Zero();
    Kr_global[nlayer].Zero();
    Kp_global[nlayer].Zero();
    
    for(int ipar=0;ipar<Npars;ipar++) {
      Der_Params_prime_wrt_thetaMS_global[nlayer][ipar][0] = 0.0;
      Der_Params_prime_wrt_thetaMS_global[nlayer][ipar][1] = 0.0;
      
      Der_Params_prime_wrt_Eloss_global[nlayer][ipar]      = 0.0;
    }
    
    for(int i=0;i<3;i++) {
      for(int iuv=0;iuv<2;iuv++) {
        DerUVwrtPos_global[nlayer][iuv][i] = 0.0;
      }
      
    }
  }

  //Getting the values of the Fr,Fp,Kr,Kp arrays which contain the derivatives of position 
  //and momentum at layer nLayer w.r.t. position and momentum at previous layer
  GetFsAndKsMS(ItersectionHitList_init,Fr_global,Fp_global,Kr_global,Kp_global);

  //At this stage it is calculated the multiple scattering angle at each geometry intersection.
  //It is also calculated a couple of othogonal unit vectors normal to the particles direction at 
  //a given geometry intersection.
  GetSigmaMSAndPerpVectors(ItersectionHitList_init,SigmaMS2_global,m1Unit_global,m2Unit_global);
  
  //At this stage it is calculated the Eloss resolution at each geometry intersection.
  //It is also calculated the unit vector pointing in the momentum direction  at the given geometry intersection.
  if(IncludeEloss) GetSigmaElossAndUnitVector(ItersectionHitList_init,SigmaEloss2_global,pUnit_global);
  
  //Calculate the derivatives of U and V w.r.t position X,Y,Z for all the layers
  GetUVDerWrtPosForAllLayers(ItersectionHitList_init,DerUVwrtPos_global);

  //These matrices are the derivatives of the local hit coordinates (1st index) w.r.t the track parameters (second index)
  //Now filling-up the "hits" convariance matrix
  for(int ihit=0;ihit<Nhits;ihit++) {
    for(int ipar=0;ipar<Npars;ipar++) GammaUDer_global[ihit][ipar] = GammaVDer_global[ihit][ipar] = 0.0;
  }
  GetUVDerivatiesWrtPars(ItersectionHitList_init,GammaUDer_global,GammaVDer_global,Npars);

  //Get the derivatives of the track parameters w.r.t the MS angles at the different planes intersected by the particle
  //This will be used later to update the track parameters covariance matrix by the effect of MS on a given layer
  GetAllTrackParamDerWrtMS(ItersectionHitList_init,SigmaMS2_global,m1Unit_global,m2Unit_global,Der_Params_prime_wrt_thetaMS_global);
  
  //Get the derivatives of the track parameters w.r.t the Eloss at the different planes intersected by the particle
  //This will be used later to update the track parameters covariance matrix by the effect of Eloss on a given layer
  if(IncludeEloss) GetAllTrackParamDerWrtEloss(ItersectionHitList,SigmaEloss2_global,pUnit_global,Der_Params_prime_wrt_Eloss_global);
  
  //Get the global covariance matrix to be used to calculate the probability of track seeding
  TMatrixD GlobalMeasCovMatrix;
  GetHitUVCovMatrixWithMSFromGlobalCalc(ItersectionHitList_init,ItersectionHitList_init,
					Fr_global,Fp_global,Kr_global,Kp_global,
					SigmaMS2_global,m1Unit_global,m2Unit_global,
					SigmaEloss2_global,pUnit_global,
					DerUVwrtPos_global,
					GlobalMeasCovMatrix);
  
  if(Myverbose) {
    cout << endl;
    cout << "full list of hits resolutions and correlations (" << SensLayersList_init.size() << "):  ";
    global->fWatch.Print();
    global->fWatch.Continue();
    for(int ihit=0;ihit<Nhits_init;ihit++) {
      int index = Nhits_init - ihit - 1;
      int indexU = 2*index;
      int indexV = 2*index+1;
      
      double resU = sqrt(GlobalMeasCovMatrix(indexU,indexU));
      double resV = sqrt(GlobalMeasCovMatrix(indexV,indexV));
      double corr = GlobalMeasCovMatrix(indexU,indexV)/(resU*resV);
      resU /= global->GetUnit("um");
      resV /= global->GetUnit("um");
      
      cout << "  hit " << index << " res(U,V) = (" << resU << "," << resV << ") um; corr = " << corr*100 << " %" << endl;
    }
  }
  
  if(Myverbose) {
    cout << endl;
    cout << "List of seed configurations (" << SeedConfigList.size() << "):  ";
    global->fWatch.Print();
    global->fWatch.Continue();
  }
  
  //Now calculate the seeding probability applying the track finder cuts (chi2Ondf and PtMin)
  for(int iseed=0;iseed<int(SeedConfigList.size());iseed++) { //begin loop over seed configurations
    
    //First decode the seed configuration from the single long number
    std::vector<int>  SeedLayersTypeList;
    SeedLayersTypeList.clear();
    global->GetListOfLayersTypes(SeedConfigList[iseed],SeedLayersTypeList);
    
    if(Myverbose) {
      cout << endl << "  config " << SeedConfigList[iseed] << " (" << iseed+1 << " out of " << SeedConfigList.size() << "):   ";
      global->fWatch.Print();
      global->fWatch.Continue();
    }
    
    //Include seed configurations with fake hits
    std::vector<long> SeedConfigListWithFakes;
    SeedConfigListWithFakes.clear();
    if(NfakesSeed_max > 0) {
      if(Myverbose) {
        cout << "Begin calculation of fakes in seed   ";
        global->fWatch.Print();
        global->fWatch.Continue();
      }
      
      global->IncludeFakesOnConfigList(SeedConfigList[iseed],NfakesSeed_max,true,SeedConfigListWithFakes);
      
      if(Myverbose) {
	cout << "End calculation of fakes in seed   ";
        global->fWatch.Print();
        global->fWatch.Continue();
	cout << "  Seed configs with Nfakes <= NfakesSeed_max = " << NfakesSeed_max << endl;
	for(int kkk=0;kkk<int(SeedConfigListWithFakes.size());kkk++) {
	  cout << "    Config: " << SeedConfigListWithFakes[kkk] << endl;
	}
      }
    }
    
    //In the next loop perform the following calculations
    // - seed intrinsic prob
    // - for each seeding layer, the probabilities for good, fake and bull hit association
    // - number of null and good hits for this seed config
    // - maximum and minimum path length of seeding layers
    
    double Rs[2]; // max and min of seeding layers path lengths
    Rs[0] = +1.0e+20;
    Rs[1] = -1.0e+20;
    double Prob_intrinsic = 1.0; //Intrinsic seeding prob
    int    Nnull_seed     = 0;   //Number of null associations in track seeding
    int    Ngood_seed     = 0;   //Number of good associations in track seeding
    std::vector<int> Insensitive_seed_layers;
    Insensitive_seed_layers.clear();
    for(int ilayer=0;ilayer<int(SeedLayersTypeList.size());ilayer++) { //begin loop over seed config layers
      int index     = Nhits - ilayer - 1;
      int the_layer = SeedLayersTypeList.size() - ilayer - 1;
      
      double sn = ItersectionHitList[SensLayersList[index]].s;
      if(Rs[0] > sn) Rs[0] = sn;
      if(Rs[1] < sn) Rs[1] = sn;
      
      if(SeedLayersTypeList[ilayer] == 0) {
	Insensitive_seed_layers.push_back(SensLayersList[index]);
	Nnull_seed++;
      }
      if(SeedLayersTypeList[ilayer] == 1) Ngood_seed++;
      
      int TheIndex = -1;
      for(int kkk=0;kkk<Nhits_init;kkk++) {
	if(ItersectionHitList[SensLayersList[index]].s == ItersectionHitList_init[SensLayersList_init[kkk]].s) {
	  TheIndex = kkk;
	  break;
	}
      }
      int indexU  = 2*TheIndex;
      int indexV  = 2*TheIndex+1;
      double resU = sqrt(GlobalMeasCovMatrix(indexU,indexU));
      double resV = sqrt(GlobalMeasCovMatrix(indexV,indexV));
      double corr = GlobalMeasCovMatrix(indexU,indexV)/(resU*resV);
      
      double ROtime   = Geometry->GetGeometryElement(ItersectionHitList[SensLayersList[index]].geoElement_idx)->GetROtime();
      double BkgRate  = Geometry->GetGeometryElement(ItersectionHitList[SensLayersList[index]].geoElement_idx)->GetBkgRate();
      double Surface  = sqrt(Chi2OndfSeedCut)*global->GetEllipseArea(resU,resV,corr);
      double Nbkg     = Surface*ROtime*BkgRate;
      double P_noFake = TMath::PoissonI(0,Nbkg);
      
      double Edet      = SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](0);
      SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](0) = Edet;                         // Prob good hit
      SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](1) = (1 - Edet)*(1.0 - P_noFake);  // Prob fake hit
      SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](2) = (1 - Edet)*P_noFake;          // Prob null hit
      
      if(Myverbose) {
	cout << "  ihit = " << index << ", ilayer = " << the_layer << "   ";
	global->fWatch.Print();
        global->fWatch.Continue();
	cout << "    ROtime   = " << ROtime/global->GetUnit("us") << " us, " << endl;
	cout << "    BkgRate  = " << BkgRate/global->GetUnit("MHz/um2") << " MHz/um2, " << endl;
	cout << "    Surface  = " << Surface/global->GetUnit("um2") << " um2, " << endl;
	cout << "    Nbkg     = " << Nbkg << ", " << endl;
	cout << "    P_noFake = " << P_noFake*100 << " %," << endl;
	cout << "    Phit     = " << SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](0)*100 << " %" << endl;
	cout << "    Pfake    = " << SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](1)*100 << " %" << endl;
	cout << "    Pnull    = " << SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](2)*100 << " %" << endl;
      }
      
      if(SeedLayersTypeList[the_layer] == 1)      Prob_intrinsic *= SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](0);
      else if(SeedLayersTypeList[the_layer] == 2) Prob_intrinsic *= SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](1);
      else if(SeedLayersTypeList[the_layer] == 0) Prob_intrinsic *= SeedConfigIntrinsicProbsList[iseed].Probs[the_layer](2);
    } //end loop over seed config layers
    
    int Nmax_null_nonSeed = Nmax_null_Layers - Nnull_seed; //Maximum number of null associations in non-Seed layers
    
    //Fill up list of layers used in track seeding
    int Nhits_seed_p = 0;
    std::vector<IntersectionHit_t> ItersectionHitList_Seed;
    ItersectionHitList_Seed.clear();
    for(int ihit=0;ihit<NLayers;ihit++) {
      double s = ItersectionHitList[ihit].s;
      if(s < Rs[0] || s > Rs[1]) continue;
      
      IntersectionHit_t AHit;
      AHit.s                = s;
      AHit.geoElement_idx   = ItersectionHitList[ihit].geoElement_idx;
      AHit.IsSensitivePoint = ItersectionHitList[ihit].IsSensitivePoint;
      
      bool IsInsensitive = false;
      if(AHit.IsSensitivePoint) {
	for(int kkk=0;kkk<int(Insensitive_seed_layers.size());kkk++) {
	  if(ihit == Insensitive_seed_layers[kkk]) {
	    IsInsensitive = true;
	    break;
	  }
	}
      }
      if(IsInsensitive) AHit.IsSensitivePoint = false;
      
      ItersectionHitList_Seed.push_back(AHit);
      if(ItersectionHitList[ihit].IsSensitivePoint) Nhits_seed_p++;
    }
    global->OrderIntersectionHitList(ItersectionHitList_Seed);
    
    //Get only those seed layers with a hit
    std::vector<int> SeedSensLayers;
    SeedSensLayers.clear();
    std::vector<TVector3> SeedSensUVCoords;
    SeedSensUVCoords.clear();
    for(int ilayer=0;ilayer<int(SeedLayersTypeList.size());ilayer++) { //begin loop over seed config layers
      int the_layer = SeedLayersTypeList.size() - ilayer - 1;
      int index = Nhits - ilayer - 1;
      if(SeedLayersTypeList[the_layer] == 0) continue;
      
      SeedSensLayers.push_back(index);
      double s     = ItersectionHitList[SensLayersList[index]].s;
      TVector3 pos = Trajectory->GetTrueTrajectoryCoordinates(s);
      TVector3 UVW = Geometry->GetGeometryElement(ItersectionHitList[SensLayersList[index]].geoElement_idx)->GetMainSurface()->GetUVWFromXYZ(pos);
      SeedSensUVCoords.push_back(UVW);
    } //end  loop over seed config layers
    const int Nseed_layers(SeedSensLayers.size());
    
    //Now get the seed layers covariance matrix
    TMatrixTSym<double>  SeedMeasCovMatrix;
    SeedMeasCovMatrix.ResizeTo(TMatrixTSym<double>(2*Nseed_layers));
    for(int ihit1=0;ihit1<Nseed_layers;ihit1++) { //begin 1st loop over seed layers
      int index1 = SeedSensLayers[ihit1];
      int TheIndex1 = -1;
      //Get global index of sensitive layer ihit1
      for(int kkk=0;kkk<Nhits_init;kkk++) {
	if(ItersectionHitList[SensLayersList[index1]].s == ItersectionHitList_init[SensLayersList_init[kkk]].s) {
	  TheIndex1 = kkk;
	  break;
	}
      }
      
      int indexU1      = 2*TheIndex1;
      int indexV1      = 2*TheIndex1+1;
      int indexU1_seed = 2*ihit1;
      int indexV1_seed = 2*ihit1+1;
      SeedMeasCovMatrix(indexU1_seed,indexU1_seed) = GlobalMeasCovMatrix(indexU1,indexU1);
      SeedMeasCovMatrix(indexV1_seed,indexV1_seed) = GlobalMeasCovMatrix(indexV1,indexV1);
      SeedMeasCovMatrix(indexU1_seed,indexV1_seed) = GlobalMeasCovMatrix(indexU1,indexV1);
      SeedMeasCovMatrix(indexV1_seed,indexU1_seed) = GlobalMeasCovMatrix(indexV1,indexU1);
            
      for(int ihit2=0;ihit2<Nseed_layers;ihit2++) {  //begin 2nd loop over seed layers
	if(ihit1 >= ihit2) continue;
	int index2 = SeedSensLayers[ihit2];
        int TheIndex2 = -1;
	//Get global index of sensitive layer ihit2
        for(int kkk=0;kkk<Nhits_init;kkk++) {
	  if(ItersectionHitList[SensLayersList[index2]].s == ItersectionHitList_init[SensLayersList_init[kkk]].s) {
	    TheIndex2 = kkk;
	    break;
	  }
        }
        
        int indexU2      = 2*TheIndex2;
        int indexV2      = 2*TheIndex2+1;
        int indexU2_seed = 2*ihit2;
        int indexV2_seed = 2*ihit2+1;
	SeedMeasCovMatrix(indexU2_seed,indexU2_seed) = GlobalMeasCovMatrix(indexU2,indexU2);
        SeedMeasCovMatrix(indexV2_seed,indexV2_seed) = GlobalMeasCovMatrix(indexV2,indexV2);
        SeedMeasCovMatrix(indexU2_seed,indexV2_seed) = GlobalMeasCovMatrix(indexU2,indexV2);
        SeedMeasCovMatrix(indexV2_seed,indexU2_seed) = GlobalMeasCovMatrix(indexV2,indexU2);
	
	SeedMeasCovMatrix(indexU1_seed,indexU2_seed) = GlobalMeasCovMatrix(indexU1,indexU2);
	SeedMeasCovMatrix(indexU2_seed,indexU1_seed) = GlobalMeasCovMatrix(indexU2,indexU1);
	SeedMeasCovMatrix(indexU1_seed,indexV2_seed) = GlobalMeasCovMatrix(indexU1,indexV2);
	SeedMeasCovMatrix(indexV2_seed,indexU1_seed) = GlobalMeasCovMatrix(indexV2,indexU1);
	SeedMeasCovMatrix(indexV1_seed,indexU2_seed) = GlobalMeasCovMatrix(indexV1,indexU2);
	SeedMeasCovMatrix(indexU2_seed,indexV1_seed) = GlobalMeasCovMatrix(indexU2,indexV1);
	SeedMeasCovMatrix(indexV1_seed,indexV2_seed) = GlobalMeasCovMatrix(indexV1,indexV2);
	SeedMeasCovMatrix(indexV2_seed,indexV1_seed) = GlobalMeasCovMatrix(indexV2,indexV1);
      }  //end 2nd loop over seed layers
    }  //end 1st loop over seed layers
    
    if(Myverbose) {
      cout << endl;
      cout << "List of resolutions and correlations of seed config " << SeedConfigList[iseed] << "(" << Nseed_layers << "):  ";
      global->fWatch.Print();
      global->fWatch.Continue();
      cout << "Nmax_null_nonSeed = " << Nmax_null_nonSeed << endl;
      
      for(int ihit=0;ihit<Nseed_layers;ihit++) {
	int index = SeedSensLayers[ihit];
	
	int indexU = 2*ihit;
        int indexV = 2*ihit+1;
	
	double resU = sqrt(SeedMeasCovMatrix(indexU,indexU));
        double resV = sqrt(SeedMeasCovMatrix(indexV,indexV));
        double corr = SeedMeasCovMatrix(indexU,indexV)/(resU*resV);
        resU /= global->GetUnit("um");
        resV /= global->GetUnit("um");
      
        cout << "  hit " << index << " res(U,V) = (" << resU << "," << resV << ") um; corr = " << corr*100 << " %" << endl;
      }
    }
    
    //Calculate the Chi2/ndf cut for seed tracks
    int    ndf            = 2*Nseed_layers - Npars;           // get ndf
    double Chi2_cut       = Chi2OndfSeedCut*ndf;              // get Chi2 from Chi2/ndf and ndf
    double Prob_Chi2Cut   = 1.0 - TMath::Prob(Chi2_cut,ndf);  // calculate Chi2/ndf cut prob
    
    //Now calculate the PtMin efficiency by a MC method
    //First diagonalize the covariance matrix for the MC generation
    // - transform to a system of uncorrelated variables
    // - in this system generate 2*Nseed_layers independent gaussian distributed random numbers
    // - transform back to the original system
    TMatrixDSymEigen AMatrix_eigen(SeedMeasCovMatrix);
    TVectorD  Eigen_values      =  AMatrix_eigen.GetEigenValues();
    TMatrixD  Eigen_vectors     =  AMatrix_eigen.GetEigenVectors();
    TMatrixD  Eigen_vectors_Inv(2*Nseed_layers,2*Nseed_layers);
    Eigen_vectors_Inv.Transpose(Eigen_vectors);
    TMatrixD Diag;
    Diag.ResizeTo(2*Nseed_layers,2*Nseed_layers);
    Diag = (Eigen_vectors_Inv*SeedMeasCovMatrix)*Eigen_vectors;
    
    bool MCGenverbose = false;
    //MCGenverbose = true;
    
    double Effic_PtCut     = 0;
    double Effic_PtCut_err = 0;
    int Nmc = aTrackFinderAlgo->GetNmcSeedEffic();
    for(int idx_gen=0;idx_gen<Nmc;idx_gen++) { //begin loop over MC trails
      TVectorD Cvals_p(2*Nseed_layers);
      TVectorD Gens(2*Nseed_layers);
      for(int i=0;i<2*Nseed_layers;i++) {
	//Generate 2*Nseed_layers independent gaussian distributed random numbers
	Cvals_p(i) = global->rand->Gaus(0.0,sqrt(Diag(i,i)));
      }
      Gens = Eigen_vectors*Cvals_p; //Transform back to the original system
      
      //For each generation check if the set of seed points pass the minimum Pt cut
      bool Passed_PtCuts = true;
      TVector3 Xc1;
      TVector3 Xc2;
      for(int ihit=0;ihit<Nseed_layers;ihit++) { //begin loop over seed config sensitive layers
	int index  = SeedSensLayers[ihit];
	
	int indexU = 2*ihit;
	int indexV = 2*ihit+1;
	
	TVector3 PosUVW = SeedSensUVCoords[ihit]; //Coordinates of the actual hit in the local (U,V) frame
	//Smering
	PosUVW(0) += Gens(indexU);
	PosUVW(1) += Gens(indexV);
	//Transform from local (U,V) to global coordinate system
	TVector3 PosXYZ = Geometry->GetGeometryElement(ItersectionHitList[SensLayersList[index]].geoElement_idx)->GetMainSurface()->GetXYZFromUVW(PosUVW);
	
	if(ihit == 0) {
	  //Use the outermost hit to define the circle passing by this point, the CenterPosition and with radius R = Ptcut/(0.3*Bfield)
	  TVector3 X1_perp = PosXYZ         - (PosXYZ.Dot(UnitBfield))        *UnitBfield; // outermost seed point projection in plane perpendicular to Bfield
	  TVector3 X0_perp = CenterPosition - (CenterPosition.Dot(UnitBfield))*UnitBfield; // CenterPosition       projection in plane perpendicular to Bfield
	  TVector3 uUnit   = X1_perp - X0_perp;
	  double d         = 0.5*uUnit.Mag();
	  uUnit            = uUnit.Unit();
	  TVector3 vUnit   = uUnit.Cross(UnitBfield);
	  TVector3 Xave    = 0.5*(X1_perp + X0_perp);
	  double  m        = pow(RadiusMin,2) - pow(d,2);
	  if(m < 0) {
	    Passed_PtCuts = true;
	    break;
	  }
	  else {
	    m = sqrt(m);
	    //Possible circle center poisitions
	    Xc1 = Xave + m*vUnit;
	    Xc2 = Xave - m*vUnit;
	  }
	  
	  if(Myverbose && MCGenverbose) {
	    cout << endl;
	    cout << "idx_gen = " << idx_gen+1 << ", index = " << index << ":" << endl;
	    cout << "X1      = (" << PosXYZ(0)/global->GetUnit("cm")  << "," << PosXYZ(1)/global->GetUnit("cm")  << "," << PosXYZ(2)/global->GetUnit("cm")  << ") cm" << endl;
	    cout << "CP      = (" << CenterPosition(0)/global->GetUnit("cm") << "," << CenterPosition(1)/global->GetUnit("cm") << "," << CenterPosition(2)/global->GetUnit("cm") << ") cm" << endl;
	    cout << "X1_perp = (" << X1_perp(0)/global->GetUnit("cm") << "," << X1_perp(1)/global->GetUnit("cm") << "," << X1_perp(2)/global->GetUnit("cm") << ") cm" << endl;
	    cout << "X0_perp = (" << X0_perp(0)/global->GetUnit("cm") << "," << X0_perp(1)/global->GetUnit("cm") << "," << X0_perp(2)/global->GetUnit("cm") << ") cm" << endl;
	    cout << "Xave    = (" << Xave(0)/global->GetUnit("cm")    << "," << Xave(1)/global->GetUnit("cm")    << "," << Xave(2)/global->GetUnit("cm")    << ") cm" << endl;
	    cout << "BUnit   = (" << UnitBfield(0) << "," << UnitBfield(1) << "," << UnitBfield(2) << ")" << endl;
	    cout << "uUnit   = (" << uUnit(0) << "," << uUnit(1) << "," << uUnit(2) << ")" << endl;
	    cout << "vUnit   = (" << vUnit(0) << "," << vUnit(1) << "," << vUnit(2) << ")" << endl;
	    cout << "d       = " << d/global->GetUnit("cm") << " cm" << endl;
	    cout << "R       = " << RadiusMin/global->GetUnit("cm") << " cm" << endl;
	    cout << "m       = " << m/global->GetUnit("cm") << " cm" << endl;
	    cout << "Xc1     = (" << Xc1(0)/global->GetUnit("cm") << "," << Xc1(1)/global->GetUnit("cm") << "," << Xc1(2)/global->GetUnit("cm") << ") cm" << endl;
	    cout << "Xc2     = (" << Xc2(0)/global->GetUnit("cm") << "," << Xc2(1)/global->GetUnit("cm") << "," << Xc2(2)/global->GetUnit("cm") << ") cm" << endl;
	  }
	  
	  continue;
	}
	
	//Apply the condition for minimum Pt cut
	TVector3 X1_perp = PosXYZ - (PosXYZ.Dot(UnitBfield))*UnitBfield;
	double r1 = (X1_perp - Xc1).Mag();
	double r2 = (X1_perp - Xc2).Mag();
	if(r1 > RadiusMin || r2 > RadiusMin) Passed_PtCuts = false;
	
	if(Myverbose && MCGenverbose) {
	  cout << endl;
	  cout << "idx_gen = " << idx_gen+1 << ", index = " << index << ":" << endl;
	  cout << "X1      = (" << PosXYZ(0)/global->GetUnit("cm")  << "," << PosXYZ(1)/global->GetUnit("cm")  << "," << PosXYZ(2)/global->GetUnit("cm")  << ") cm" << endl;
	  cout << "X1_perp = (" << X1_perp(0)/global->GetUnit("cm") << "," << X1_perp(1)/global->GetUnit("cm") << "," << X1_perp(2)/global->GetUnit("cm") << ") cm" << endl;
	  cout << "r1      = " << r1/global->GetUnit("cm") << " cm" << endl;
	  cout << "r2      = " << r1/global->GetUnit("cm") << " cm" << endl;
	  cout << "R       = " << RadiusMin/global->GetUnit("cm") << " cm" << endl;
	}
	
	if(!Passed_PtCuts) break;
      } //end loop over seed config sensitive layers
      
      //Check if this trail passed min pt cut
      if(Passed_PtCuts) Effic_PtCut += 1.0;
      
    } //end loop over MC trails
    
    //Caculate minimum pt cut efficiency and error
    Effic_PtCut /= Nmc;
    Effic_PtCut_err = sqrt(Effic_PtCut*(1.0 - Effic_PtCut)/Nmc);
    
    SeedConfigPtCutProbsList[iseed]   = Effic_PtCut;                              // Min Pt cut prob
    SeedConfigChi2CutProbsList[iseed] = Prob_Chi2Cut;                             // Seed track Maximum Chi2/ndf cut
    SeedConfigTotProbsList[iseed]     = Prob_intrinsic*Prob_Chi2Cut*Effic_PtCut;  // Total seeding prob
    
    if(Myverbose) {
      TVector3 p0 = (1.0/global->GetMomentumUnit("GeV/c"))*Trajectory->GetInitMomentum();
      double   pt = (p0 - (p0.Dot(UnitBfield))*UnitBfield).Mag();
      double   ptcut = (1.0/global->GetMomentumUnit("GeV/c"))*PtCutMin;
    
      cout << "  " << SeedConfigList[iseed] << " : Prob_intrinsic = " << Prob_intrinsic*100 << " %, Chi2/ndf_cut prob = " << Prob_Chi2Cut*100 << " %," 
           << " PtMin_cut prob = (" << Effic_PtCut*100 << " +/- " << Effic_PtCut_err*100 << ") %; full prob = " << SeedConfigTotProbsList[iseed]*100 << " %; "
	   << "pt = " << pt << " GeV/c; Pt_cut = " << ptcut << " GeV/c    ";
      global->fWatch.Print();
      global->fWatch.Continue();
    }
    
    if(SeedConfigTotProbsList[iseed] < MinSeedProb) continue;
    
    //Track parameters covariance parameters by fitting the track seeding hits
    TMatrixD             FitCovMatrix;
    TMatrixD             MeasCovMatrix;
    GetHitUVCovMatrixWithMSFromGlobalCalc(ItersectionHitList,ItersectionHitList_Seed,
					  Fr_global,Fp_global,Kr_global,Kp_global,
					  SigmaMS2_global,m1Unit_global,m2Unit_global,
					  SigmaEloss2_global,pUnit_global,
					  DerUVwrtPos_global,
					  MeasCovMatrix);
    bool IsGoodFit = GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc(ItersectionHitList,ItersectionHitList_Seed,GammaUDer_global,GammaVDer_global,MeasCovMatrix,FitCovMatrix);

    if(Myverbose) {
      cout << endl;
      if(IsGoodFit) {
	cout << "  Seed par cov matrix:" << endl;
	cout << "    Is good fit" << endl;
        FitCovMatrix.Print();
        cout << endl;
      }
      else {
	cout << "  Seed par cov matrix is not good!" << endl;
      }
    }
    if(!IsGoodFit) continue;
    
    if(Myverbose) {
      cout << endl;
      cout << "  Ngood_seed = " << Ngood_seed << ", Nnull_seed = " << Nnull_seed << endl;
    }
    
    //Get the layers not used for track seeding
    std::vector<IntersectionHit_t> ItersectionHitList_NonSeed;
    ItersectionHitList_NonSeed.clear();
    for(int ilayer=0;ilayer<NLayers;ilayer++) {
      if(ItersectionHitList[ilayer].s >= Rs[0]) continue;
      ItersectionHitList_NonSeed.push_back(ItersectionHitList[ilayer]);
    }
    global->OrderIntersectionHitList(ItersectionHitList_NonSeed);
    
    //Now count the configurations for the non-seed layers
    std::vector<int>  Position_nonSeed;
    std::vector<long> NonSeedLayersConfigList;
    NonSeedLayersConfigList.clear();
    for(int inull=0;inull<Nmax_null_nonSeed+1;inull++) { //begin of loop over number of null associations on non-seed layers
      //The corresponding number of layers involved in track seedingNmax_null_Layers
      int Nnull_nonseed = inull;
      int Nhits_nonseed = Nhits - (Ngood_seed + Nnull_seed);
      if(Nhits_nonseed <= 0) continue;
      
      int Ngood_nonseed = Nhits_nonseed - Nnull_nonseed;
      if(Ngood_nonseed < 0) continue;
      
      if(Myverbose) cout << "  Nnull_layers_nonseed = " << Nnull_nonseed << ", Nhits_nonseed = " << Nhits_nonseed << ", and Ngood_nonseed = " << Ngood_nonseed  << endl;
      
      const int tmp_Nhits_nonseed(Nhits_nonseed);
      int Config[tmp_Nhits_nonseed];
      int HitKind[tmp_Nhits_nonseed];
      for(int k=0;k<Nhits_nonseed;k++) Config[k] = -1;
      
      int counter_kind = 0;
      for(int jgood=0;jgood<Ngood_nonseed;jgood++) {
	HitKind[counter_kind] = 1;
	counter_kind++;
      }
      for(int jnull=0;jnull<Nnull_nonseed;jnull++) {
	HitKind[counter_kind] = 3;
	counter_kind++;
      }

      //Fill the configurations for the combinations of good and null associations for the track seeding
      global->FillConfigurations(Nhits_nonseed,Position_nonSeed,Nhits_nonseed,Config,HitKind,NonSeedLayersConfigList);
      
    } //end of loop over number of null associations on non-seed layers
    
    std::vector<HitsConfiguration_t> NonSeedHitConfiguration;
    NonSeedHitConfiguration.clear();
    for(int iconf_nonseed=0;iconf_nonseed<int(NonSeedLayersConfigList.size());iconf_nonseed++) { //Begin of loop over non-seed configurations
      std::vector<int>  NonSeedLayersTypeList;
      NonSeedLayersTypeList.clear();
      global->GetListOfLayersTypes(NonSeedLayersConfigList[iconf_nonseed],NonSeedLayersTypeList);
      const int NonSeed_Layers(NonSeedLayersTypeList.size());
      
      //Cut on the number of hits of this config
      std::vector<int>  NonSeedLayersTypeListNew;
      NonSeedLayersTypeListNew.clear();
      for(int inonseed=0;inonseed<NonSeed_Layers;inonseed++) {
	if(NonSeedLayersTypeList[inonseed] == 3) NonSeedLayersTypeListNew.push_back(NonSeedLayersTypeList[inonseed]);
	else                                     NonSeedLayersTypeListNew.push_back(1);
      }
      long New_config = global->GetConfigFromListOfLayersTypes(NonSeedLayersTypeListNew);
      
      //bool IsIn = false;
      //for(int k=0;k<int(NonSeedHitConfiguration.size());k++) {
	//if(NonSeedHitConfiguration[k].HitConfig == New_config) {
	  //NonSeedHitConfiguration[k].CorrectAndFakeHitConfigList.push_back(NonSeedLayersConfigList[iconf_nonseed]);
	  //IsIn = true;
	//}
      //}
      //if(!IsIn) {
	HitsConfiguration_t AConfig;
        AConfig.HitConfig = New_config;
	AConfig.Probs.clear();
	AConfig.CorrectAndFakeHitConfigList.clear();
	AConfig.CorrectAndFakeHitConfigList.push_back(NonSeedLayersConfigList[iconf_nonseed]);
	NonSeedHitConfiguration.push_back(AConfig);
      //}
    }
    
    if(Myverbose) { 
      cout << endl;
      cout <<  "Non-seed configurations   ";
      global->fWatch.Print();
      global->fWatch.Continue();
      int counter_config = 0;
      for(int k=0;k<int(NonSeedHitConfiguration.size());k++) {
        cout << "config = " << NonSeedHitConfiguration[k].HitConfig << endl;
        for(int j=0;j<int(NonSeedHitConfiguration[k].CorrectAndFakeHitConfigList.size());j++) {
	  //cout << " - Sub-config = " << NonSeedHitConfiguration[k].CorrectAndFakeHitConfigList[j] << endl;
	  counter_config++;
        }
      }
      cout << "total N-configs = " << counter_config << ", initial total configs = " << NonSeedLayersConfigList.size() << endl;
      cout << endl;
    }
    
   for(int iconf_nonseed=0;iconf_nonseed<int(NonSeedHitConfiguration.size());iconf_nonseed++) { //begin of loop over non-seed configurations
      std::vector<int>  NonSeedLayersTypeList;
      NonSeedLayersTypeList.clear();
      global->GetListOfLayersTypes(NonSeedHitConfiguration[iconf_nonseed].HitConfig,NonSeedLayersTypeList);
      const int NonSeed_Layers(NonSeedLayersTypeList.size());

      if(Myverbose) {
	cout << endl << "Non-Seed layers config, for seed-config " << iseed+1 << "(" << SeedConfigList[iseed] << ", out of " << SeedConfigList.size()
	     << ") and non-seed config " << iconf_nonseed+1 << " (out of " << NonSeedLayersConfigList.size() << ") : ";
        for(int inonseed=0;inonseed<NonSeed_Layers;inonseed++) {
	  if(NonSeedLayersTypeList[inonseed] == 3) cout << 0 << "  ";
	  else                                     cout << NonSeedLayersTypeList[inonseed] << "  ";
	}
	global->fWatch.Print();
        global->fWatch.Continue();
	cout << endl << endl;
      }

      //Repository of the already navigated layers
      std::vector<IntersectionHit_t> ItersectionHitList_Navigated;
      ItersectionHitList_Navigated.clear();
      //Start by adding the layers used for track seeding
      for(int ilayer=0;ilayer<int(ItersectionHitList_Seed.size());ilayer++) ItersectionHitList_Navigated.push_back(ItersectionHitList_Seed[ilayer]);
      global->OrderIntersectionHitList(ItersectionHitList_Navigated);

      if(Myverbose) {
        cout << endl;
        cout << "Start navigation of nonseed-config " << iconf_nonseed+1 << " of seed-config " << iseed+1 << "   ";
	global->fWatch.Print();
        global->fWatch.Continue();
      }

      //Now start the proceed to inward navigation from the track-seeding layers
      TMatrixD FitCovMatrix_prev;
      TMatrixD FitCovMatrix_current;
      FitCovMatrix_prev.ResizeTo(Npars,Npars);
      FitCovMatrix_current.ResizeTo(Npars,Npars);
      FitCovMatrix_prev = FitCovMatrix;
      int counter_sens_layer = NonSeed_Layers;
      for(int inonseed_layer=0;inonseed_layer<int(ItersectionHitList_NonSeed.size());inonseed_layer++) { //begin of loop inward navigation on non-seed layers
        int layer_idx      = ItersectionHitList_NonSeed.size() - inonseed_layer - 1;
	int geoElement_idx = ItersectionHitList_NonSeed[layer_idx].geoElement_idx;

	if(Myverbose) {
	  cout << " - geoElement_idx = " << geoElement_idx << "   ";
	  global->fWatch.Print();
          global->fWatch.Continue();
	}

        //Add current layer to the list navigated layers
	IntersectionHit_t AHit;
        AHit.s                = ItersectionHitList_NonSeed[layer_idx].s;
	AHit.geoElement_idx   = ItersectionHitList_NonSeed[layer_idx].geoElement_idx;
        AHit.IsSensitivePoint = ItersectionHitList_NonSeed[layer_idx].IsSensitivePoint;
	if(AHit.IsSensitivePoint) {
	  counter_sens_layer -= 1;
	  if(NonSeedLayersTypeList[counter_sens_layer] == 3) AHit.IsSensitivePoint = false;
	}
        ItersectionHitList_Navigated.push_back(AHit);
        global->OrderIntersectionHitList(ItersectionHitList_Navigated);

        if(!ItersectionHitList_NonSeed[layer_idx].IsSensitivePoint) {
	  if(Myverbose) {
	    cout << "    * Layer was not initially sensitive!!! Perform correction of par cov matrix due to MS." << endl;
	    Geometry->GetGeometryElement(geoElement_idx)->Print();
	  }

	  //If layer was initially non-sensitive, then recalculate the parameters covariance matrix by adding the layer MS contribution
	  GetNewParCovarianceMatrix(ItersectionHitList,ItersectionHitList_Navigated,0,
			            Der_Params_prime_wrt_thetaMS_global,SigmaMS2_global,
			            Der_Params_prime_wrt_Eloss_global,SigmaEloss2_global,
			            FitCovMatrix_prev,FitCovMatrix_current);
	  FitCovMatrix_prev = FitCovMatrix_current;
	  
	  if(Myverbose) {
	    cout << "Interated par cov matrix from MS:" << endl;
	    FitCovMatrix_current.Print();
	    cout << endl;
	  }

	}
        else {
	  //If the layer was initially sensitive calculate the track extrapolation cone into the layer surface
	  TMatrixD CovUV;
	  GetCovMatrixUVOfTrackIntersection(ItersectionHitList,ItersectionHitList_Navigated,0,
				            GammaUDer_global,GammaVDer_global,FitCovMatrix_prev,CovUV);

	  //Use this information to calculate probabilities of correct, fake and null associantion
	  GGeoObject* aGeoElement = Geometry->GetGeometryElement(geoElement_idx);
	  
	  //Get the geo-element intrinsic resoltion
	  double sn_hit       = ItersectionHitList_NonSeed[layer_idx].s;
	  TVector3 pos_n      = Trajectory->GetTrueTrajectoryCoordinates(sn_hit);
	  TVector3 mom_n      = Trajectory->GetTrueTrajectoryMon(        sn_hit);
	  TVector2 Resolution = Geometry->GetResolutionUV(geoElement_idx,mom_n,pos_n);
	  
	  int    ndf_add      = 2;
	  double sigmaU       = sqrt(CovUV(0,0) + pow(Resolution.X(),2));
	  double sigmaV       = sqrt(CovUV(1,1) + pow(Resolution.Y(),2));
	  double corr         = CovUV(0,1)/(sigmaU*sigmaV);
	  double det_effic    = Geometry->GetEfficiency(geoElement_idx,mom_n,pos_n,Trajectory->GetParticle());
	  double rate         = aGeoElement->GetBkgRate()*aGeoElement->GetROtime();
	  TVector3  TheProbs  = global->ProbsForTrackClusterAssociation(det_effic,rate,corr,sigmaU,sigmaV,aTrackFinderAlgo->GetChi2OndfAddCut(),ndf_add);
	  NonSeedHitConfiguration[iconf_nonseed].Probs.push_back(TheProbs);

	  if(Myverbose) {
	    cout << "    * Layer was initially sensitive!!!" << endl;
	    cout << "      - (sigmaU,sigmaV,corr) = (" << sigmaU/global->GetDistanceUnit("um") << "," << sigmaV/global->GetDistanceUnit("um") << "," << corr*100 << ") (um,um,%)"<< endl;
	    cout << "      - det_effic            =  " << det_effic*100 << " %" << endl;
	    aGeoElement->Print();
	    cout << "      - ROtime               =  " << aGeoElement->GetROtime()/global->GetUnit("us") << " us" << endl;
	    cout << "      - bkgRate              =  " << aGeoElement->GetBkgRate()/global->GetUnit("MHz/um2") << " Mhz/um2" << endl;
	    cout << "      - rate                 =  " << rate/(1.0/global->GetUnit("um2")) << " 1/um2" << endl;
	    cout << "      - area                 =  " << TMath::Pi()*(1.0 - pow(corr,2))*sigmaU*sigmaV/global->GetUnit("um2") << " um2" << endl;
	    cout << "      - S                    =  " << TMath::Pi()*(1.0 - pow(corr,2))*sigmaU*sigmaV*rate << " hits" << endl;
	    cout << "      - chi2/ndf add         =  " << aTrackFinderAlgo->GetChi2OndfAddCut() << endl;
	    cout << "      - ndf add              =  " << ndf_add << endl;
	    cout << "      - chi2 add             =  " << aTrackFinderAlgo->GetChi2OndfAddCut()*ndf_add << endl;
	    cout << "      - Pcorr                =  " << TheProbs(0)*100 << " %" << endl;
	    cout << "      - Pfake                =  " << TheProbs(1)*100 << " %" << endl;
	    cout << "      - Pnull                =  " << TheProbs(2)*100 << " %" << endl;
	  }

	  if(NonSeedLayersTypeList[counter_sens_layer] == 1) {
	    if(Myverbose) cout << "    * Layer was initially sensitive, and with hit (either correct of fake association). Perform track fit including this point!!!" << endl;
	    //Recalculate the parameter covariance matrix including this point
	    TMatrixD  MeasCovMatrix_tmp;
            GetHitUVCovMatrixWithMSFromGlobalCalc(ItersectionHitList,ItersectionHitList_Navigated,
						  Fr_global,Fp_global,Kr_global,Kp_global,
						  SigmaMS2_global,m1Unit_global,m2Unit_global,
						  SigmaEloss2_global,pUnit_global,
						  DerUVwrtPos_global,
						  MeasCovMatrix_tmp);
	    IsGoodFit = GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc(ItersectionHitList,ItersectionHitList_Navigated,GammaUDer_global,GammaVDer_global,MeasCovMatrix_tmp,FitCovMatrix_current);
	    FitCovMatrix_prev = FitCovMatrix_current;
	    
	    if(Myverbose) {
	      cout << "Interated par cov matrix, adding new hit to chi2 fit:" << endl;
	      FitCovMatrix_current.Print();
	      cout << endl;
	    }
          }
	  else {
	    if(Myverbose) cout << "    * Layer was initially sensitive, but missed hit.  Perform correction of par cov matrix due to MS." << endl;

	    GetNewParCovarianceMatrix(ItersectionHitList,ItersectionHitList_Navigated,0,
			              Der_Params_prime_wrt_thetaMS_global,SigmaMS2_global,
			              Der_Params_prime_wrt_Eloss_global,SigmaEloss2_global,
			              FitCovMatrix_prev,FitCovMatrix_current);
	    FitCovMatrix_prev = FitCovMatrix_current;
	    
	    if(Myverbose) {
	      cout << "Interated par cov matrix from MS:" << endl;
	      FitCovMatrix_current.Print();
	      cout << endl;
	    }
	  }

        }

        if(layer_idx == 0) {
	  //Got the final layer, either sensitive or not. 
	  //Perform fit with all the layers traversed y the particle to get the track parameters uncertainties
	  TMatrixD  MeasCovMatrix_tmp;
          GetHitUVCovMatrixWithMSFromGlobalCalc(ItersectionHitList,ItersectionHitList_Navigated,
						Fr_global,Fp_global,Kr_global,Kp_global,
						SigmaMS2_global,m1Unit_global,m2Unit_global,
						SigmaEloss2_global,pUnit_global,
						DerUVwrtPos_global,
						MeasCovMatrix_tmp);
	  IsGoodFit = GetFitTrackParsCovMatrix_FromPlaneList_FromGlobalCalc(ItersectionHitList,ItersectionHitList_Navigated,GammaUDer_global,GammaVDer_global,MeasCovMatrix_tmp,FitCovMatrix_current);
	  FitCovMatrix_prev = FitCovMatrix_current;
	  
	  if(Myverbose) {
	    cout << "Adding last intersected layer. Recalculate trk-params cov matrix with all layers:   ";
	    global->fWatch.Print();
            global->fWatch.Continue();
	    FitCovMatrix_current.Print();
	    cout << endl;
	  }
	  
	}

      } //end of loop inward navigation on non-seed layers
      if(Myverbose) {
        cout << endl;
        cout << "End   navigation of nonseed-config " << iconf_nonseed+1 << " of seed-config " << iseed+1 << "   ";
	global->fWatch.Print();
        global->fWatch.Continue();
      }
      
      if(Myverbose) {
	cout << endl;
	cout << "    Tracking prob for the seed config " << SeedConfigList[iseed] << " and non seeding config " << NonSeedHitConfiguration[iconf_nonseed].HitConfig << " :    ";
	global->fWatch.Print();
        global->fWatch.Continue();
        cout << "        Prob_seed = " << SeedConfigTotProbsList[iseed] << endl;
      }
      
      //Calculate the tracking probability for this seed and non-seed configuration
      double Prob_track_config = SeedConfigTotProbsList[iseed]; //Start with the seeding probability
      for(int kkk=0;kkk<NonSeed_Layers;kkk++) {
	int index_kkk = NonSeed_Layers - kkk - 1;
	
	double TheProb = 1.0;
	if(NonSeedLayersTypeList[index_kkk] == 1)      TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](0);
	else if(NonSeedLayersTypeList[index_kkk] == 3) TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](2);
	
	if(Myverbose) {
	    cout << "        " << NonSeedLayersTypeList[index_kkk] << ", prob   = " << TheProb << "  "
	         << "(" << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](0) 
		 << "," << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](1) 
		 << "," << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](2)
		 << ")"
	    << endl;
	  }
	
	//Incude the probabilities for including or excluding a hit as the case might be
	Prob_track_config *= TheProb;
      }
      
      //Tracking efficiency contributions with no fakes
      Efficiencies.Effic_tot     += Prob_track_config;
      Efficiencies.Effic_NoFakes += Prob_track_config;
      
      //Adding up all the fit param covariance matrix to weighted wrt the track config efficiency
      AveFitCovMatrix += FitCovMatrix_prev*Prob_track_config;
      
      if(Myverbose) {
	cout << "        Prob_tot  = " << Prob_track_config << endl;
	cout << endl;
	cout << "        Prob_tot_add          = " << Efficiencies.Effic_tot << endl;
	cout << "        Prob_NoFakes_add      = " << Efficiencies.Effic_NoFakes << endl;
	cout << "        Prob_1Fake_add        = " << Efficiencies.Effic_1Fake << endl;
	cout << "        Prob_2orMoreFakes_add = " << Efficiencies.Effic_2orMoreFakes << endl;
	cout << endl;
      }
      
      if(Nfakes_max == 0) continue;

      //Now lets include those configurations with fakes on the seeding and non-seeding layers respecting the Purity cut

      //Start with configurations with no fakes in seeding
      //Get the non-seeding configurations with fakes
      if(Myverbose) {
        cout << "Begin calculation of fakes in non-seed (fakesless seed)   ";
        global->fWatch.Print();
        global->fWatch.Continue();
      }
      std::vector<long> NonSeedConfigListWithFakes;
      NonSeedConfigListWithFakes.clear();
      global->IncludeFakesOnConfigList(NonSeedHitConfiguration[iconf_nonseed].HitConfig,Nfakes_max,false,NonSeedConfigListWithFakes);
      if(Myverbose) {
        cout << "End   calculation of fakes in non-seed (fakesless seed)   ";
        global->fWatch.Print();
        global->fWatch.Continue();
      }
      
      if(Myverbose) { 
        cout << endl;
	cout << "  Seed config " << SeedConfigList[iseed] << " and non-seed config " << NonSeedHitConfiguration[iconf_nonseed].HitConfig << ". List of non-seed configs with fakes:    ";
	global->fWatch.Print();
        global->fWatch.Continue();
	cout << "  NonSeed configs with Nfakes <= Nfakes_max = " << Nfakes_max << endl;
	for(int kkk=0;kkk<int(NonSeedConfigListWithFakes.size());kkk++) {
	  cout << "    Config: " << NonSeedConfigListWithFakes[kkk] << endl;
	}
      }

      //Now loop over then and calculate tracking probability
      for(int iconf_fakes=0;iconf_fakes<int(NonSeedConfigListWithFakes.size());iconf_fakes++) { //begin loop over non-seed configs with fakes
	//Start with the seeding probability with no fakes
	Prob_track_config = SeedConfigTotProbsList[iseed];
	
	std::vector<int>  NonSeedLayersTypeList_WidthFakes;
	NonSeedLayersTypeList_WidthFakes.clear();
	global->GetListOfLayersTypes(NonSeedConfigListWithFakes[iconf_fakes],NonSeedLayersTypeList_WidthFakes);
	
	int MyNgood = 0;
	int MyNfake = 0;
	int MyNnull = 0;
	global->GetNGoodFakesAndNull(NonSeedConfigListWithFakes[iconf_fakes],MyNgood,MyNfake,MyNnull);
	
	if(Myverbose) {
	  cout << endl;
	  cout << "    Tracking prob for the seed config " << SeedConfigList[iseed] << " and non seeding config " << NonSeedConfigListWithFakes[iconf_fakes] << " : "
	       << "Nfakes = " << MyNfake << " (Nfakes_max = " << Nfakes_max << ")     ";
	  global->fWatch.Print();
          global->fWatch.Continue();
	  cout << "        Prob_seed = " << SeedConfigTotProbsList[iseed] << endl;
	}
	
	for(int kkk=0;kkk<int(NonSeedLayersTypeList_WidthFakes.size());kkk++) {
	  int index_kkk = NonSeedLayersTypeList_WidthFakes.size() - kkk - 1;
	  
	  double TheProb = 1.0;
	  if(NonSeedLayersTypeList_WidthFakes[index_kkk] == 1)      TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](0);
	  if(NonSeedLayersTypeList_WidthFakes[index_kkk] == 2)      TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](1);
	  else if(NonSeedLayersTypeList_WidthFakes[index_kkk] == 3) TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](2);
	
	  if(Myverbose) {
	    cout << "        " << NonSeedLayersTypeList_WidthFakes[index_kkk] << ", prob   = " << TheProb << "  "
	         << "(" << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](0) 
		 << "," << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](1) 
		 << "," << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](2)
		 << ")"
	         << endl;
	  }
	  
	  //Incude the probabilities for including or excluding a hit as the case might be
	  Prob_track_config *= TheProb;
	}
	
	//Tracking efficiency contributions with no fakes
        Efficiencies.Effic_tot          += Prob_track_config;
        if(MyNfake == 0)      Efficiencies.Effic_NoFakes      += Prob_track_config;
	else if(MyNfake == 1) Efficiencies.Effic_1Fake        += Prob_track_config;
	else if(MyNfake >  1) Efficiencies.Effic_2orMoreFakes += Prob_track_config;
      
        //Adding up all the fit param covariance matrix to weighted wrt the track config efficiency
        AveFitCovMatrix += FitCovMatrix_prev*Prob_track_config;
	
	if(Myverbose) {
	  cout << "        Prob_tot  = " << Prob_track_config << endl;
	  cout << endl;
	  cout << "        Prob_tot_add          = " << Efficiencies.Effic_tot << endl;
	  cout << "        Prob_NoFakes_add      = " << Efficiencies.Effic_NoFakes << endl;
	  cout << "        Prob_1Fake_add        = " << Efficiencies.Effic_1Fake << endl;
	  cout << "        Prob_2orMoreFakes_add = " << Efficiencies.Effic_2orMoreFakes << endl;
	  cout << endl;
	}
	
      } //end loop over non-seed configs with fakes

      //Now continue with configurations with fakes in seeding
      if(Myverbose) {
	cout << endl;
	cout << "  Seed configs with Nfakes <= NfakesSeed_max = " << NfakesSeed_max << " from no-fakes config " << SeedConfigList[iseed] << "   ";
	global->fWatch.Print();
        global->fWatch.Continue();
      }
      for(int iconf_seedfakes=0;iconf_seedfakes<int(SeedConfigListWithFakes.size());iconf_seedfakes++) { //begin loop over seed configs with fakes
	//continue;
	int MyNgood_seed = 0;
	int MyNfake_seed = 0;
	int MyNnull_seed = 0;
	global->GetNGoodFakesAndNull(SeedConfigListWithFakes[iconf_seedfakes],MyNgood_seed,MyNfake_seed,MyNnull_seed);
	
	//Decode the seed configuration with fakes from the single long number
        std::vector<int>  SeedLayersTypeList_WithFakes;
        SeedLayersTypeList_WithFakes.clear();
        global->GetListOfLayersTypes(SeedConfigListWithFakes[iconf_seedfakes],SeedLayersTypeList_WithFakes);
	
	//Calculate the seeding prob with fakes
	if(Myverbose) {
	  cout << "    Config: " << SeedConfigListWithFakes[iconf_seedfakes] << " (" << SeedLayersTypeList_WithFakes.size() << ")    ";
	  global->fWatch.Print();
          global->fWatch.Continue();
	}
	double Prob_seed_fakes = 1.0;
	for(int kkk=0;kkk<int(SeedLayersTypeList_WithFakes.size());kkk++) {
	  int index_kkk = SeedLayersTypeList_WithFakes.size() - kkk - 1;
	  
	  double TheProb = 1.0;
	  if(SeedLayersTypeList_WithFakes[index_kkk] == 1)      TheProb = SeedConfigIntrinsicProbsList[iseed].Probs[index_kkk](0);
	  else if(SeedLayersTypeList_WithFakes[index_kkk] == 2) TheProb = SeedConfigIntrinsicProbsList[iseed].Probs[index_kkk](1);
	  else if(SeedLayersTypeList_WithFakes[index_kkk] == 0) TheProb = SeedConfigIntrinsicProbsList[iseed].Probs[index_kkk](2);
	  
	  if(Myverbose) {
	    cout << "        " << SeedLayersTypeList_WithFakes[index_kkk] << ", prob   = " << TheProb << "  "
	         << "(" << SeedConfigIntrinsicProbsList[iseed].Probs[index_kkk](0) 
		 << "," << SeedConfigIntrinsicProbsList[iseed].Probs[index_kkk](1)
		 << "," << SeedConfigIntrinsicProbsList[iseed].Probs[index_kkk](2)
		 << ")"
	         << endl;
	  }
	  
	  Prob_seed_fakes *= TheProb;
	  
	}
	//Assume same pt cut efficiency and chi2/ndf prob the same as for the seed configs with no fakes
	if(Myverbose) {
	  cout << "        Prob_ptCut   = " << SeedConfigPtCutProbsList[iseed] << endl;
	  cout << "        Prob_chi2Cut = " << SeedConfigChi2CutProbsList[iseed] << endl;
	}
	Prob_seed_fakes *= SeedConfigPtCutProbsList[iseed];
	Prob_seed_fakes *= SeedConfigChi2CutProbsList[iseed];
	
	if(Myverbose) {
	  cout << "    seed config no-fakes " << SeedConfigList[iseed] << ", Prob = " << SeedConfigTotProbsList[iseed]*100 << " %; "
	       << " seed config with fakes " << SeedConfigListWithFakes[iconf_seedfakes] << ", Prob = " << Prob_seed_fakes*100 << " %; "
	       << endl;
	}
	
	if(Prob_seed_fakes < MinSeedProb) continue;
	
	int Nfakes_max_tmp = Nfakes_max - MyNfake_seed;
	if(Nfakes_max_tmp <= 0) continue;
	
	//Get the non-seeding configurations with fakes
	if(Myverbose) {
          cout << "Begin calculation of fakes in non-seed (seed with fakes)   ";
          global->fWatch.Print();
          global->fWatch.Continue();
        }
        NonSeedConfigListWithFakes.clear();
        global->IncludeFakesOnConfigList(NonSeedHitConfiguration[iconf_nonseed].HitConfig,Nfakes_max_tmp,false,NonSeedConfigListWithFakes);
	if(Myverbose) {
          cout << "End   calculation of fakes in non-seed (seed with fakes)   ";
          global->fWatch.Print();
          global->fWatch.Continue();
        }
	
	if(Myverbose) { 
          cout << endl;
	  cout << "  Seed config " << SeedConfigListWithFakes[iseed] << " and non-seed config " << NonSeedHitConfiguration[iconf_nonseed].HitConfig << ". List of non-seed configs with fakes:    ";
	  global->fWatch.Print();
          global->fWatch.Continue();
	  cout << "  NonSeed configs with Nfakes <= Nfakes_max = " << Nfakes_max << endl;
	  for(int kkk=0;kkk<int(NonSeedConfigListWithFakes.size());kkk++) {
	    cout << "    Config: " << NonSeedConfigListWithFakes[kkk] << endl;
	  }
        }
	
	//Now loop over then and calculated tracking probability
	for(int iconf_fakes=0;iconf_fakes<int(NonSeedConfigListWithFakes.size());iconf_fakes++) { //begin loop over non-seed configs with fakes
	  //Start with the seeding probability with fakes
	  Prob_track_config = Prob_seed_fakes;
	
	  std::vector<int>  NonSeedLayersTypeList_WidthFakes;
	  NonSeedLayersTypeList_WidthFakes.clear();
	  global->GetListOfLayersTypes(NonSeedConfigListWithFakes[iconf_fakes],NonSeedLayersTypeList_WidthFakes);

	  for(int kkk=0;kkk<int(NonSeedLayersTypeList_WidthFakes.size());kkk++) {
	    int index_kkk = NonSeedLayersTypeList_WidthFakes.size() - kkk - 1;
	    
	    double TheProb = 1.0;
	    if(NonSeedLayersTypeList_WidthFakes[index_kkk] == 1)      TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](0);
	    if(NonSeedLayersTypeList_WidthFakes[index_kkk] == 2)      TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](1);
	    else if(NonSeedLayersTypeList_WidthFakes[index_kkk] == 3) TheProb = NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](2);
	
	    if(Myverbose) {
	      cout << "        " << NonSeedLayersTypeList_WidthFakes[index_kkk] << ", prob   = " << TheProb << "  "
	           << "(" << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](0) 
		   << "," << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](1) 
		   << "," << NonSeedHitConfiguration[iconf_nonseed].Probs[index_kkk](2)
		   << ")"
	           << endl;
	    }
	    
	    //Incude the probabilities for including or excluding a hit as the case might be
	    Prob_track_config *= TheProb;
	  }
	
	  int MyNgood = 0;
	  int MyNfake = 0;
	  int MyNnull = 0;
	  global->GetNGoodFakesAndNull(NonSeedConfigListWithFakes[iconf_fakes],MyNgood,MyNfake,MyNnull);
	
	  //Tracking efficiency contributions with no fakes
          Efficiencies.Effic_tot          += Prob_track_config;
          if(MyNfake == 0)      Efficiencies.Effic_NoFakes      += Prob_track_config;
	  else if(MyNfake == 1) Efficiencies.Effic_1Fake        += Prob_track_config;
	  else if(MyNfake >  1) Efficiencies.Effic_2orMoreFakes += Prob_track_config;
      
          //Adding up all the fit param covariance matrix to weighted wrt the track config efficiency
          AveFitCovMatrix += FitCovMatrix_prev*Prob_track_config;
	
	  if(Myverbose) {
	    cout << endl;
	    cout << "    Tracking prob for the seed config " << SeedConfigListWithFakes[iconf_seedfakes] << " and non seeding config " << NonSeedConfigListWithFakes[iconf_fakes] << " : "
	         << "Nfakes = " << MyNfake << " (Nfakes_max = " << Nfakes_max << "); "
	         << "Prob = " << Prob_track_config*100 << " %    ";
            global->fWatch.Print();
            global->fWatch.Continue();
	    cout << "        Prob_tot_add          = " << Efficiencies.Effic_tot << endl;
	    cout << "        Prob_NoFakes_add      = " << Efficiencies.Effic_NoFakes << endl;
	    cout << "        Prob_1Fake_add        = " << Efficiencies.Effic_1Fake << endl;
	    cout << "        Prob_2orMoreFakes_add = " << Efficiencies.Effic_2orMoreFakes << endl; 
	    cout << endl;
	  }

        } //end loop over non-seed configs with fakes
      } //end loop over seed configs with fakes
      
    } //end of loop over non-seed configurations
    
  } //end loop over seed configurations
  

  if(Efficiencies.Effic_tot*100 > 1.0e-6) {
    AveFitCovMatrix *= 1.0/Efficiencies.Effic_tot;
  }

  if(Myverbose) {
    cout << "Total              efficiency = " << Efficiencies.Effic_tot*100          << " %" << endl;
    cout << "NoFakes            efficiency = " << Efficiencies.Effic_NoFakes*100      << " %" << endl;
    cout << "At least one fake  efficiency = " << Efficiencies.Effic_1Fake*100        << " %" << endl;
    cout << "At least two fakes efficiency = " << Efficiencies.Effic_2orMoreFakes*100 << " %" << endl;
    cout << endl;
    cout << endl;
    
    cout << endl;
    cout << "End   GTracker::GetFitTrackPseudoEfficiency_FPCCD  ";
    global->fWatch.Print();
    global->fWatch.Continue();
    cout << endl;
  }

  return;
  
}
//====================================================================



