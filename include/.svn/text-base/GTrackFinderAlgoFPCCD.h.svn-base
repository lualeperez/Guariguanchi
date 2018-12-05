/***************************************************************************//**
 * @brief:      
 * @Description: Base class for FPCCD track finder algorithm
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

#ifndef GTrackFinderAlgoFPCCD_h
#define GTrackFinderAlgoFPCCD_h

#include "TMath.h"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TCutG.h>
#include <TRandom.h>
#include "TStopwatch.h"
#include "include/GGlobalTools.h"
#include "include/GTrackFinderAlgo.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GTrackFinderAlgoFPCCD : public GTrackFinderAlgo {
  
private:
  
  std::vector<TString>       FullSystemList;    //Full list of layers used for track finding
  std::vector<SeedLayers_t>  SeedLayerConfigs;  //List seeding configurations
  
  int       Nhits_min;           //Minimum number of hits for tracking
  double    PtMin_cut;           //Minimum Pt cut for track seeding
  double    PurityMin_cut;       //Cut on minimum track purity, purity = (# track good hits) / (# track good hits  +  track fake hits)
  double    Chi2Ondf_Seed_cut;   //Cut on maximum Chi2/ndf at seeding
  double    Chi2Ondf_Add_cut;    //Cut on maximum Chi2/ndf at adding hit to track
  bool      InwardTracking;      //bool to specify if doing inward or outward tracking
  int       NfakesMax_seeding;   //Maximum number of fakes during seeding
  TVector3  CenterPosition;     //Coordinates of the center from whicn the Pt cut is applied
  
  int       Nmc_seed_effic;      //Number of trails for seeding efficiency calculation
  
public:
  
  // Constructor
  GTrackFinderAlgoFPCCD(TString   aName,
			std::vector<TString>       aFullSystemList,
			std::vector<SeedLayers_t>  aSeedLayerConfigs,
			TrackFindingRegion_t       aRegion,
			GGlobalTools* aglobal,
			int       aNhits_min         = -1,
			double    aPtMin_cut         = Dummy_value,
			double    aPurityMin_cut     = 1.0,
			double    aChi2Ondf_Seed_cut = 3.0,
			double    aChi2Ondf_Add_cut  = 3.0,
			bool      aInwardTracking    = true,
			int       aNfakesMax_seeding = 0,
			TVector3  aCenterPosition    = TVector3(0,0,0));
  
  GTrackFinderAlgoFPCCD(const GTrackFinderAlgoFPCCD& other, TString Name = TString(""));

  // Destructor
  virtual ~GTrackFinderAlgoFPCCD();
  
  //Set of generic functions
  
  //Functions to access internal variables
  int       GetNhitsMin()          const { return  Nhits_min; }
  double    GetPtMinCut()          const { return  PtMin_cut; }
  double    GetPurityMinCut()      const { return  PurityMin_cut; }
  double    GetChi2OndfSeedCut()   const { return  Chi2Ondf_Seed_cut; }
  double    GetChi2OndfAddCut()    const { return  Chi2Ondf_Add_cut; }
  bool      GetInwardTracking()    const { return  InwardTracking; }
  int       GetNmcSeedEffic()      const { return  Nmc_seed_effic; }
  int       GetNfakesMaxSeeding()  const { return  NfakesMax_seeding; }
  TVector3  GetCenterPosition()    const { return  CenterPosition; }
  
  std::vector<TString>       GetFullSystemList()   const { return  FullSystemList; }
  std::vector<SeedLayers_t>  GetSeedLayerConfigs() const { return  SeedLayerConfigs; }
  
  //Functions to set internal variables
  void      SetNhitsMin(int aNhits_min);
  void      SetPtMinCut(double aPtMin_cut);
  void      SetPurityMinCut(double aPurityMin_cut);
  void      SetChi2OndfSeedCut(double aChi2Ondf_Seed_cut);
  void      SetChi2OndfAddCut(double aChi2Ondf_Add_cut);
  void      SetInwardTracking(bool aInwardTracking);
  void      SetNmcSeedEffic(int aNmc_seed_effic);
  void      SetNfakesMaxSeeding(int aNfakesMax_seeding);
  void      SetCenterPosition(TVector3 aCenterPosition)  { CenterPosition = aCenterPosition; }
  
  void      FillSystemList(std::vector<TString>  aFullSystemList);
  void      FillSeedConfigs(std::vector<SeedLayers_t>  aSeedLayerConfigs);
  
  //Set of functions to be defined for the daughter classes
  GTrackFinderAlgoFPCCD*  clone(TString aName) const;
  void      CheckInputs();
  void      Print();

protected:
  
  bool verbose;
  
};

#endif //~ GTrackFinderAlgoFPCCD_h

