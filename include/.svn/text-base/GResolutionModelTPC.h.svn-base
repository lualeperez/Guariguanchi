/***************************************************************************//**
 * @brief:      
 * @Description: Class for TPC Resolution models
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
 

#ifndef GResolutionModelTPC_h
#define GResolutionModelTPC_h

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
#include "include/GResolutionModel.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GResolutionModelTPC : public GResolutionModel {
  
private:
  
public:
  
  double    a_rphi;
  double    b_rphi;
  double    c_rphi;
  double    a_z;
  double    b_z;
  
  // Constructor
  GResolutionModelTPC(TString  aName,
		      double  aa_rphi, double  ab_rphi, double  ac_rphi,
		      double  aa_z,    double  ab_z,
		      GGlobalTools* aglobal);
  
  GResolutionModelTPC(const GResolutionModelTPC& other, TString Name = TString(""));

  // Destructor
  virtual ~GResolutionModelTPC();
  
  //Set of generic functions
  
  //Functions to access internal variables
  
  //Functions to set internal variables
  
  //Set of functions to be defined for the daughter classes
  GResolutionModel*  clone(TString aName) const;
  void      CheckInputs();
  TVector2  GetResolution(TVector2 ResolutionUV, TVector3 mom, TVector3 pos);
  void      Print();

protected:

};

#endif //~ GResolutionModel_h

