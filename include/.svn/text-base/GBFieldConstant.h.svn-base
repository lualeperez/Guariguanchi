/***************************************************************************//**
 * @brief:      
 * @Description: Class for constant B-fields
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


#ifndef GBFieldConstant_h
#define GBFieldConstant_h

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

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GBFieldConstant : public GBField {
  
private:
  
public:
  
  TVector3 GlobalBField;  // Constant B-field value
  
  // Constructor
  GBFieldConstant(TString   aName,
		  TVector3  aBField,
		  GGlobalTools* aglobal);
  
  GBFieldConstant(TString   aName,
		  TVector3  aBFieldDirection,
		  float     aBFieldMagnitude,
		  GGlobalTools* aglobal);
  
  GBFieldConstant(TString   aName,
		  float     aBFieldXComponent,
		  float     aBFieldYComponent,
		  float     aBFieldZComponent,
		  GGlobalTools* aglobal);
  
  GBFieldConstant(const GBFieldConstant& other, TString Name = TString(""));

  // Destructor
  virtual ~GBFieldConstant();
  
  //Set of generic functions
  
  //Functions to access internal variables
  TVector3  GetConstantBField() const { return GlobalBField; }
  
  //Functions to set internal variables
  void      SetConstantBField(TVector3 aBField) { GlobalBField = aBField; }
  
  //Set of functions to be defined for the daughter classes
  GBField*  clone(TString aName) const;
  TVector3  GetBFieldValue(TVector3 PositionXYZ) { return GlobalBField; }
  void      Print();

protected:

};

#endif //~ GBFieldConstant_h

