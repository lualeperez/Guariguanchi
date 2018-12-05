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
#include "include/GResolutionModelTPC.h"

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
GResolutionModelTPC::GResolutionModelTPC(TString   aName,
					 double  aa_rphi, double  ab_rphi, double  ac_rphi,
					 double  aa_z,    double  ab_z,
					 GGlobalTools* aglobal)
                                         : GResolutionModel(aName,aglobal)
{
  
  Type   = TString("TPC");
  a_rphi = aa_rphi;
  b_rphi = ab_rphi;
  c_rphi = ac_rphi;
  a_z    = aa_z;
  b_z    = ab_z;
  
  CheckInputs();
  
}
//====================================================================
GResolutionModelTPC::GResolutionModelTPC(const GResolutionModelTPC& other,TString aName)
                                         : GResolutionModel(aName,other.global)
{
  
  Type   = other.Type;
  a_rphi = other.a_rphi;
  b_rphi = other.b_rphi;
  c_rphi = other.c_rphi;
  a_z    = other.a_z;
  b_z    = other.b_z;
  
  CheckInputs();
  
}
//====================================================================
GResolutionModelTPC::~GResolutionModelTPC() 
{
  
}
//====================================================================
void  GResolutionModelTPC::CheckInputs()
{
  
  if(a_rphi < 0) {
    cout << endl;
    cout << "ERROR in GResolutionModelTPC::CheckInputs: a_rphi parameter is negtive. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  else if(b_rphi < 0) {
    cout << endl;
    cout << "ERROR in GResolutionModelTPC::CheckInputs: b_rphi parameter is negtive. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  else if(c_rphi < 0) {
    cout << endl;
    cout << "ERROR in GResolutionModelTPC::CheckInputs: c_rphi parameter is negtive. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  else if(a_z < 0) {
    cout << endl;
    cout << "ERROR in GResolutionModelTPC::CheckInputs: a_z parameter is negtive. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(b_z < 0) {
    cout << endl;
    cout << "ERROR in GResolutionModelTPC::CheckInputs: b_z parameter is negtive. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
GResolutionModel* GResolutionModelTPC::clone(TString aName) const
{
 
  return new GResolutionModelTPC(*this,aName);
  
}
//====================================================================
TVector2  GResolutionModelTPC::GetResolution(TVector2 ResolutionUV, TVector3 mom, TVector3 pos)
{
  
  //sigma^2_rphi (U) = a_rphi^2 + b_rphi^2 * sin^2(phi) + c_rphi * sin(theta) * z
  //sigma^2_z    (V) = a_z^2    + b_z * z
  
  double sinPhi     = mom.Y()/sqrt(pow(mom.X(),2) + pow(mom.Y(),2));
  double sinTheta   = sqrt(pow(mom.X(),2) + pow(mom.Y(),2))/mom.Mag();
  double Z          = TMath::Abs(pos.Z())/global->GetUnit("cm");
      
  double ResU = sqrt(pow(a_rphi,2) + pow(b_rphi*sinPhi,2) + c_rphi*sinTheta*Z);
  double ResV = sqrt(pow(a_z,2) + b_z*Z);
  
  return TVector2(ResU,ResV);
  
}
//====================================================================
void  GResolutionModelTPC::Print()
{
  
  cout << "  Begin TPC Resolution Model" << endl;
  cout << "    ModelType  " << Type.Data() << endl;
  cout << "    ModelName  " << Name.Data() << endl;
  cout << "    a_rphi     " << a_rphi/global->GetDistanceUnit("um") << " um" << endl;
  cout << "    b_rphi     " << b_rphi/global->GetDistanceUnit("um") << " um" << endl;
  cout << "    c_rphi     " << c_rphi                               << ""    << endl;    
  cout << "    a_z        " << a_z/global->GetDistanceUnit("um")    << " um" << endl;
  cout << "    b_z        " << b_z                                  << ""    << endl;
  cout << "  End   TPC Resolution Model" << endl;
  
  return;
  
}
//====================================================================

