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
#include "include/GBFieldMultipleSteps.h"

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
GBFieldMultipleSteps::GBFieldMultipleSteps(TString       aName,
					   std::vector<TVector3>     aInBFieldList,
					   std::vector<GGeoObject*>  aVolumeList,
					   TVector3  aOutBField,
					   GGlobalTools* aglobal)
                                           : GBField(aName,aglobal)
{
  
  SetBFields(aInBFieldList,aOutBField);
  
  SetVolumes(aVolumeList);
  
  Type      = TString("MultipleSteps");
  
  CheckInputs();
  
}
//====================================================================
GBFieldMultipleSteps::GBFieldMultipleSteps(const GBFieldMultipleSteps& other,TString aName)
                                           : GBField(aName,other.global)
{
  
  SetBFields(other.InBFieldList,other.OutBField);
  SetVolumes(other.VolumeList);
  Type      = other.Type;
  
  CheckInputs();
  
}
//====================================================================
GBFieldMultipleSteps::~GBFieldMultipleSteps() 
{
  
  for(int i=0;i<int(VolumeList.size());i++) {
    if(VolumeList[i] != NULL) delete VolumeList[i];
  }
  VolumeList.clear();
  InBFieldList.clear();
  
}
//====================================================================
GBField* GBFieldMultipleSteps::clone(TString aName) const
{
 
  return new GBFieldMultipleSteps(*this,aName);
  
}
//====================================================================
void   GBFieldMultipleSteps::SetBFields(std::vector<TVector3> aInBFieldList, TVector3 aOutBField)
{
  
  //Set the internal fields, both inside the volumes and outside.
  
  OutBField = aOutBField;
  
  InBFieldList.clear();
  for(int i=0;i<int(aInBFieldList.size());i++) InBFieldList.push_back(aInBFieldList[i]);
  
  return;
  
}
//====================================================================
void   GBFieldMultipleSteps::SetVolumes(std::vector<GGeoObject*> aVolumeList)
{
  
  //Set the volumes
  
  VolumeList.clear();
  for(int i=0;i<int(aVolumeList.size());i++) VolumeList.push_back(aVolumeList[i]->clone(aVolumeList[i]->GetName()));
  
  return;
  
}
//====================================================================
void  GBFieldMultipleSteps::CheckInputs(void)
{
  
  if(VolumeList.size() == 0 || InBFieldList.size() == 0) {
    cout << endl;
    cout << "ERROR inside GBFieldMultipleSteps::CheckInputs:: Volume or inBfield list has zero size. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  if(InBFieldList.size() != VolumeList.size()) {
    cout << endl;
    cout << "ERROR inside GBFieldMultipleSteps::CheckInputs:: Volume and inBfield list have different sizes. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  double min_bfield = 1.0*global->GetUnit("uT");
  
  bool AllZerofields = true;
  for(int i=0;i<int(InBFieldList.size());i++) {
    if(InBFieldList[i].Mag() > min_bfield) {
      AllZerofields = false;
      break;
    }
  }
  if(AllZerofields) {
    if(OutBField.Mag() > min_bfield) AllZerofields = false;
  }
  
  if(AllZerofields) {
    cout << endl;
    cout << "ERROR inside GBFieldMultipleSteps::CheckInputs:: all specified fields are zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  CheckVolumes();
  
  return;
  
}
//====================================================================
void  GBFieldMultipleSteps::CheckVolumes(void)
{
  
  //Check for overlaps among the volumes
  
  double Precision = 1.0*global->GetUnit("mm");
  
  if(VolumeList.size() == 1) return;
  
  int Ngeo_checks = VolumeList.size()*(VolumeList.size() - 1)/2;

  if(verbose) {  
    cout << endl;
    cout << "Begin of volumes overlap check for GBFieldMultipleSteps " << Name.Data() << endl;
    cout << "  N volumes = " << VolumeList.size() << ". Will perform " << Ngeo_checks << " overlap checks!" << endl;
  }
  
  int counter = 0;
  for(int i=0;i<int(VolumeList.size());i++) {
    for(int j=0;j<int(VolumeList.size());j++) {
      if(i >= j) continue;
      
      counter++;
      
      if(verbose) {
        cout << "  Checking overlap between Volumes " << i+1 << " (" << VolumeList[i]->GetName().Data() << ") and "
             << j+1 << "(" << VolumeList[j]->GetName().Data() << "). Check number " << counter << " out of " << Ngeo_checks << "."
             << endl;
      }

      if(DoGeometryElementsOverlap(VolumeList[i],VolumeList[j],global,Precision)) {
	cout << endl;
	cout << "    ERROR inside GBFieldMultipleSteps::CheckVolumes:: Overlap found between volumes " << i+1 <<  "(" << VolumeList[i]->GetName().Data() << ") and " 
	     << j+1 << " (" << VolumeList[j]->GetName().Data() << "). "
	     << endl;
	cout << endl;
	cout << "Volume " << i+1 << " : " << endl;
	VolumeList[i]->Print();
	cout << endl;
	cout << "Volume " << j+1 << " : " << endl;
	VolumeList[j]->Print();
	cout << endl;
	cout << "Check your inputs. Exiting now!!!" << endl;
	cout << endl;
	assert(false);
      }
      
    }
  }
  
  if(verbose) {
    cout << "  No overlaps found!!!" << endl;
    cout << "End of overlap check for GBFieldMultipleSteps " << Name.Data() << endl;
    cout << endl;
  }
  
  return;
  
}
//====================================================================
TVector3  GBFieldMultipleSteps::GetInBField(int idx)  const
{
  
  //Get the B-field inside volume of index idx
  
  if(idx < 0 || idx > int(InBFieldList.size()-1)) {
    cout << endl;
    cout << " ERROR inside GBFieldMultipleSteps::GetInBField:: index " << idx << " is outside allowed limits [0," << InBFieldList.size()-1 << "]" << endl;
    cout << endl;
    assert(false);
  }
  
  return  InBFieldList[idx];
  
}
//====================================================================
GGeoObject*  GBFieldMultipleSteps::GetVolume(int idx) const
{
  
  //Get Volume with index idx in volume list
  
  if(idx < 0 || idx > int(VolumeList.size()-1)) {
    cout << endl;
    cout << " ERROR inside GBFieldMultipleSteps::GetVolume:: index " << idx << " is outside allowed limits [0," << VolumeList.size()-1 << "]" << endl;
    cout << endl;
    assert(false);
  }
  
  return  VolumeList[idx];
  
}
//====================================================================
TVector3  GBFieldMultipleSteps::GetBFieldValue(TVector3 PositionXYZ)
{
  
  //Get the B-field at position PositionXYZ
  
  for(int i=0;i<int(InBFieldList.size());i++) {
    if(VolumeList[i]->IsPointInsideGeometry(PositionXYZ)) return  InBFieldList[i];
  }
  
  return OutBField;
  
}
//====================================================================
void  GBFieldMultipleSteps::Print()
{
  
  TVector3 unitVect;
  double   Mag;
  
  cout << "  Begin Multiple Steps B Field" << endl;
  for(int i=0;i<int(InBFieldList.size());i++) {
    unitVect  = InBFieldList[i].Unit();
    Mag       = InBFieldList[i].Mag();
    cout << "    Begin Volume " << i+1 << endl;
    cout << "      B-field Magnitude = "  << Mag/global->GetUnit("T") << " Tesla" << endl;
    cout << "      B-field Direction = (" << unitVect.X() << "," << unitVect.Y() << "," << unitVect.Z() << ")" << endl;
    cout << "      Volume:" << endl;
    VolumeList[i]->Print();
    cout << "    End   Volume " << i+1 << endl;
  }
  unitVect  = OutBField.Unit();
  Mag       = OutBField.Mag();
  cout << "    Begin Outside B-field" << endl;
  cout << "      B-field Magnitude = "  << Mag/global->GetUnit("T") << " Tesla" << endl;
  cout << "      B-field Direction = (" << unitVect.X() << "," << unitVect.Y() << "," << unitVect.Z() << ")" << endl;
  cout << "    End   Outside B-field" << endl;
  cout << "  End   Multiple Steps B Field" << endl;
  
}
//====================================================================


