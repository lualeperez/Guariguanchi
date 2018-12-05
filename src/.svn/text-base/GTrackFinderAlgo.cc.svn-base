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
#include "include/GTrackFinderAlgo.h"

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
GTrackFinderAlgo::GTrackFinderAlgo(TString   aName,
				   TrackFindingRegion_t aRegion,
				   GGlobalTools* aglobal)
{
  
  Name    = aName;
  global  = aglobal;
  Type    = TString("");
  verbose = false;
  SetRegion(aRegion);
  
}
//====================================================================
GTrackFinderAlgo::GTrackFinderAlgo(const GTrackFinderAlgo& other,TString aName)
{
  
  Name    = aName;
  global  = other.global;
  Type    = other.Type;
  verbose = other.verbose;
  SetRegion(other.Region);
  
}
//====================================================================
GTrackFinderAlgo::~GTrackFinderAlgo() 
{
  
}
//====================================================================
GTrackFinderAlgo* GTrackFinderAlgo::clone(TString aName) const
{
 
  return new GTrackFinderAlgo(*this,aName);
  
}
//====================================================================
void  GTrackFinderAlgo::SetRegion(TrackFindingRegion_t aRegion)
{
  
  Region.posRange.Rx[0]     = aRegion.posRange.Rx[0];      Region.posRange.Rx[1]     = aRegion.posRange.Rx[1];
  Region.posRange.Ry[0]     = aRegion.posRange.Ry[0];      Region.posRange.Ry[1]     = aRegion.posRange.Ry[1];
  Region.posRange.Rz[0]     = aRegion.posRange.Rz[0];      Region.posRange.Rz[1]     = aRegion.posRange.Rz[1];
  Region.posRange.Rr[0]     = aRegion.posRange.Rr[0];      Region.posRange.Rr[1]     = aRegion.posRange.Rr[1];
  Region.posRange.Rtheta[0] = aRegion.posRange.Rtheta[0];  Region.posRange.Rtheta[1] = aRegion.posRange.Rtheta[1];
  Region.posRange.Rphi[0]   = aRegion.posRange.Rphi[0];    Region.posRange.Rphi[1]   = aRegion.posRange.Rphi[1];
  
  Region.momRange.Rp[0]     = aRegion.momRange.Rp[0];      Region.momRange.Rp[1]     = aRegion.momRange.Rp[1];
  Region.momRange.Rtheta[0] = aRegion.momRange.Rtheta[0];  Region.momRange.Rtheta[1] = aRegion.momRange.Rtheta[1];
  Region.momRange.Rphi[0]   = aRegion.momRange.Rphi[0];    Region.momRange.Rphi[1]   = aRegion.momRange.Rphi[1];
  
  return;
  
}
//====================================================================

