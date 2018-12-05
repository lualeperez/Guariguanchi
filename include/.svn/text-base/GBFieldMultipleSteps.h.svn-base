/***************************************************************************//**
 * @brief:      
 * @Description: Class for multiple steps B-fields, i.e. a B-field which has different constants values
 *               inside a list of volumes, and another constant value outside the set of volumes
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


#ifndef GBFieldMultipleSteps_h
#define GBFieldMultipleSteps_h

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
#include "include/GGlobalTools.h"
#include "include/GBField.h"
#include "include/GGeoObject.h"
#include "include/GUtiliratyFunctions.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GBFieldMultipleSteps : public GBField {
  
private:
  
public:
  
  std::vector<TVector3> InBFieldList;  // List constant B-fields inside volumes
  std::vector<GGeoObject*> VolumeList; // List of volumes
  TVector3 OutBField;                  // Constant B-field outside volumes
  
  // Constructor
  GBFieldMultipleSteps(TString   aName,
	               std::vector<TVector3> aInBFieldList,
		       std::vector<GGeoObject*>  aVolumeList,
	               TVector3  aOutBField,
		       GGlobalTools* aglobal);
  
  GBFieldMultipleSteps(const GBFieldMultipleSteps& other, TString Name = TString(""));

  // Destructor
  ~GBFieldMultipleSteps();
  
  //Set of generic functions
  
  //Functions to access internal variables
  int         GetNVolumes()  const  { return  VolumeList.size(); }
  TVector3    GetInBField(int idx)  const;
  GGeoObject* GetVolume(int idx) const;
  TVector3    GetOutBField() const { return  OutBField; }
  
  //Functions to set internal variables
  void      SetBFields(std::vector<TVector3> aInBFieldlist, TVector3 aOutBField);
  void      SetVolumes(std::vector<GGeoObject*> aVolumeList);
  
  void      CheckInputs(void);
  
  void      CheckVolumes(void);
  
  //Set of functions to be defined for the daughter classes
  GBField*  clone(TString aName) const;
  TVector3  GetBFieldValue(TVector3 PositionXYZ);
  void      Print();

protected:

};

#endif //~ GBFieldMultipleSteps_h

