///////////////////////////////////////////////////////////////////
// In this files are implemented several tool functions used for //
// the analysitical calculations of tracking performances        //
///////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCutG.h"
#include "TMath.h"
#include "TString.h"
#include "TColor.h"
#include "include/GGlobalTools.h"

#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

//====================================================================
GGlobalTools::GGlobalTools()
{
  
  //Class constructor
  
  MaterialMap.clear();

  ParticleMap.clear();
  
  units.clear();
  units_distance.clear();
  units_time.clear();
  units_angle.clear();
  units_Energy.clear();
  units_Momentum.clear();
  units_Mass.clear();
  units_Bfield.clear();
  units_RateDensity.clear();
  
  betaCurvRadius = -999.9;
  
  FillTables(); //Fill up the list of particles and materials
  
  InfiniteMon = 30.0*units["GeV/c"];
  
  //Used for E-loss calculation
  K           = 0.307*units["MeV"]*units["cm2"]/units["gr"]; //K = 4*pi*Na*re^2*mec^2
  
  //Used for the fluctuations on the E-loss, sigma(Eloss) = kappa x sqrt(Eloss)
  //kappa       = 0.07*units["GeV"];
  kappa       = 0.015*units["GeV"]; //Tunned by comparing with ILD fullsim
  
  rand = new TRandom(8028300);
    
}
//====================================================================
GGlobalTools::~GGlobalTools()
{
  
  //Class destructor
  delete rand;
  
}
//====================================================================
TVector3  GGlobalTools::GetMomentum(double p,double theta,double phi)
{
  
  //returns the three-momentum vector from: module, polar and azymutal angles
  
  double px = p*TMath::Sin(theta)*TMath::Cos(phi);
  double py = p*TMath::Sin(theta)*TMath::Sin(phi);
  double pz = p*TMath::Cos(theta);
  
  return  TVector3(px,py,pz);
  
}
//====================================================================
void  GGlobalTools::GetSurfAndMomOrthVects(TVector3  momentum,
					   TVector3& m1,
					   TVector3& m2)
{
  
  //Returns a set of two vectors othogonal to the momentum
  
  m1 = momentum.Orthogonal();
  m1 = m1.Unit();
  m2 = momentum.Cross(m1);
  m2 = m2.Unit();
  
  return;
  
}
//====================================================================
TVector3  GGlobalTools::ProbsForTrackClusterAssociation(double det_effic,
							double rate,
							double corr,
							double sigmaU,
							double sigmaV,
							double chi2Ondf,
							int    ndf)
{
  
  //Calculation of the probabilities of the different possibilities of the track-hit association,
  // - good association
  // - fake association
  // - no   association
  
  double gamma   = TMath::Prob(chi2Ondf*ndf,ndf);
  int NtimesCone = sqrt(chi2Ondf);
  double Surface = NtimesCone*GetEllipseArea(sigmaU,sigmaV,corr); //Cone shadow on surface
  double S       = rate*Surface; //Number of background particles inside shadows
  
  //Probability calculation for correct, null and fake hit associations
  double Prob_good_association = det_effic*(1.0 - pow(gamma,1.0 + S))/(1.0 + S);
  double Prob_null_association = (1.0 - det_effic + det_effic*gamma)*pow(gamma,S);
  double Prob_fake_association = 1.0 - Prob_good_association - Prob_null_association;
  
  //Saving the probabilities inside a 3-vector
  TVector3 Prob_associations(Prob_good_association,Prob_fake_association,Prob_null_association);
  
  return Prob_associations;
  
}
//====================================================================
double  GGlobalTools::MSAngle(TString Aparticle,
			      TVector3 momentum,
			      double xOverX0)
{
  
  //Calculation of the RMS of multiple scattering angle as a function of particle's
  // mass, momentum and the thickness (X/X0) of the traversed material
  
  double p0    = 13.6*units["MeV/c"];
  
  double p     = momentum.Mag();
  double mass  = GetParticleMass(Aparticle);
  double beta  = p/sqrt(pow(p,2) + pow(mass,2));
  int    z     = GetParticleCharge(Aparticle);
  
  double theta2 = pow(p0/(beta*p),2)*TMath::Abs(xOverX0)*pow(z,2);
  double lt = 1.0 + 0.038*TMath::Log(TMath::Abs(xOverX0));
  if(lt > 0) theta2 *= pow(lt,2);
  
  return sqrt(theta2)*units["rad"];
  
}
//====================================================================
void GGlobalTools::FillTables(void)
{
  
  //This function fills up several tables and maps
  
  TString Reference_dist_unit;
  Reference_dist_unit = TString("mm");
  
  TString Reference_PEM_unit;
  Reference_PEM_unit = TString("MeV");
  
  //Fill units table
  if(Reference_dist_unit == TString("um")) {
    // reference distance unit = um
    units[TString("nm")] = 1.0e-3;
    units[TString("um")] = 1.0;
    units[TString("mm")] = 1.0e+3;
    units[TString("cm")] = 1.0e+4;
    units[TString("m")]  = 1.0e+6;
  }
  else if(Reference_dist_unit == TString("mm")) {
    // reference distance unit = mm
    units[TString("nm")] = 1.0e-6;
    units[TString("um")] = 1.0e-3;
    units[TString("mm")] = 1.0;
    units[TString("cm")] = 1.0e+1;
    units[TString("m")]  = 1.0e+3;
  }
  else if(Reference_dist_unit == TString("cm")) {
    // refence distance unit = cm
    units[TString("nm")] = 1.0e-7;
    units[TString("um")] = 1.0e-4;
    units[TString("mm")] = 1.0e-1;
    units[TString("cm")] = 1.0;
    units[TString("m")]  = 1.0e+2;
  }
  else {
    // default reference for distance unit is mm
    units[TString("nm")] = 1.0e-6;
    units[TString("um")] = 1.0e-3;
    units[TString("mm")] = 1.0;
    units[TString("cm")] = 1.0e+1;
    units[TString("m")]  = 1.0e+3;
  }
  units_distance[TString("nm")] = units[TString("nm")];
  units_distance[TString("um")] = units[TString("um")];
  units_distance[TString("mm")] = units[TString("mm")];
  units_distance[TString("cm")] = units[TString("cm")];
  units_distance[TString("m")]  = units[TString("m")];
  
  if(Reference_PEM_unit == TString("GeV")) {
    //Reference Momentum, Mass and Energy Units: GeV/c, GeV/c2 and GeV
    //Momentum Units:
    units[TString("GeV/c")] = 1.0;
    units[TString("TeV/c")] = 1.0e+3*units[TString("GeV/c")];
    units[TString("MeV/c")] = 1.0e-3*units[TString("GeV/c")];
    units[TString("keV/c")] = 1.0e-6*units[TString("GeV/c")];
    units[TString("eV/c")]  = 1.0e-9*units[TString("GeV/c")];
  }
  else if(Reference_PEM_unit == TString("MeV")) {
    //Reference Momentum, Mass and Energy Units: MeV/c, MeV/c2 and MeV
    //Momentum Units:
    units[TString("MeV/c")] = 1.0;
    units[TString("TeV/c")] = 1.0e+6*units[TString("MeV/c")];
    units[TString("GeV/c")] = 1.0e+3*units[TString("MeV/c")];
    units[TString("keV/c")] = 1.0e-3*units[TString("MeV/c")];
    units[TString("eV/c")]  = 1.0e-6*units[TString("MeV/c")];
  }
  else if(Reference_PEM_unit == TString("keV")) {
    //Reference Momentum, Mass and Energy Units: keV/c, keV/c2 and keV
    //Momentum Units:
    units[TString("keV/c")] = 1.0;
    units[TString("eV/c")]  = 1.0e-3*units[TString("keV/c")];
    units[TString("MeV/c")] = 1.0e+3*units[TString("keV/c")];
    units[TString("GeV/c")] = 1.0e+6*units[TString("keV/c")];
    units[TString("TeV/c")] = 1.0e+9*units[TString("keV/c")];
  }
  else {
    //Default Momentum, Mass and Energy Units: GeV/c, GeV/c2 and GeV
    //Momentum Units:
    units[TString("GeV/c")] = 1.0;
    units[TString("TeV/c")] = 1.0e+3*units[TString("GeV/c")];
    units[TString("MeV/c")] = 1.0e-3*units[TString("GeV/c")];
    units[TString("keV/c")] = 1.0e-6*units[TString("GeV/c")];
    units[TString("eV/c")]  = 1.0e-9*units[TString("GeV/c")];
  }
  
  //Mass Units:
  units[TString("eV/c2")]  = units[TString("eV/c")];
  units[TString("keV/c2")] = units[TString("keV/c")];
  units[TString("MeV/c2")] = units[TString("MeV/c")];
  units[TString("GeV/c2")] = units[TString("GeV/c")];
  units[TString("TeV/c2")] = units[TString("TeV/c")];
  //Energy Units:
  units[TString("eV")]     = units[TString("eV/c")];
  units[TString("keV")]    = units[TString("keV/c")];
  units[TString("MeV")]    = units[TString("MeV/c")];
  units[TString("GeV")]    = units[TString("GeV/c")];
  units[TString("TeV")]    = units[TString("TeV/c")];
  
  units_Momentum[TString("eV/c")]  = units[TString("eV/c")];
  units_Momentum[TString("keV/c")] = units[TString("keV/c")];
  units_Momentum[TString("MeV/c")] = units[TString("MeV/c")];
  units_Momentum[TString("GeV/c")] = units[TString("GeV/c")];
  units_Momentum[TString("TeV/c")] = units[TString("TeV/c")];
  units_Energy[TString("eV")]      = units[TString("eV")];
  units_Energy[TString("keV")]     = units[TString("keV")];
  units_Energy[TString("MeV")]     = units[TString("MeV")];
  units_Energy[TString("GeV")]     = units[TString("GeV")];
  units_Energy[TString("TeV")]     = units[TString("TeV")];
  units_Mass[TString("eV/c2")]     = units[TString("eV/c2")];
  units_Mass[TString("keV/c2")]    = units[TString("keV/c2")];
  units_Mass[TString("MeV/c2")]    = units[TString("MeV/c2")];
  units_Mass[TString("GeV/c2")]    = units[TString("GeV/c2")];
  units_Mass[TString("TeV/c2")]    = units[TString("TeV/c2")];
  
  //Regular mass units
  units[TString("gr")] = 1.0;
  units[TString("Kg")] = 1.0e+3*units[TString("gr")];
  
  //Now other derived units from the ones above
  
  //1/momentums untis
  units[TString("1/(eV/c)")]  = 1.0/units[TString("eV/c")];
  units[TString("1/(keV/c)")] = 1.0/units[TString("keV/c")];
  units[TString("1/(MeV/c)")] = 1.0/units[TString("MeV/c")];
  units[TString("1/(GeV/c)")] = 1.0/units[TString("GeV/c")];
  units[TString("1/(TeV/c)")] = 1.0/units[TString("TeV/c")];
  
  //Angle Units: reference is rad
  units[TString("rad")]        = 1.0;
  units[TString("mrad")]       = 1.0e-3*units[TString("rad")];
  units[TString("urad")]       = 1.0e-6*units[TString("rad")];
  units[TString("deg")]        = (TMath::Pi()/180.0)*units[TString("rad")];
  units_angle[TString("rad")]  = units[TString("rad")];
  units_angle[TString("mrad")] = units[TString("mrad")];
  units_angle[TString("urad")] = units[TString("urad")];
  units_angle[TString("deg")]  = units[TString("deg")];
  
  //Time units: reference is s (second)
  units[TString("s")]         = 1.0;
  units[TString("us")]        = 1.0e-6*units[TString("s")];
  units[TString("ns")]        = 1.0e-9*units[TString("s")];
  units[TString("ms")]        = 1.0e-3*units[TString("s")];
  units[TString("mim")]       = 6.0e+1*units[TString("s")];
  units[TString("hour")]      = 3.6e+3*units[TString("s")];
  units_time[TString("us")]   = units[TString("us")];
  units_time[TString("ns")]   = units[TString("ns")];
  units_time[TString("ms")]   = units[TString("ms")];
  units_time[TString("s")]    = units[TString("s")];
  units_time[TString("mim")]  = units[TString("mim")];
  units_time[TString("hour")] = units[TString("hour")];
  
  //Rate units: reference is Hz
  units[TString("Hz")]   = 1.0e+0/units["s"];
  units[TString("kHz")]  = 1.0e+3/units["s"];
  units[TString("MHz")]  = 1.0e+6/units["s"];
  units[TString("GHz")]  = 1.0e+9/units["s"];
  
  //Surface units: derived from distance units
  units[TString("nm2")] = pow(units["nm"],2);
  units[TString("um2")] = pow(units["um"],2);
  units[TString("mm2")] = pow(units["mm"],2);
  units[TString("cm2")] = pow(units["cm"],2);
  units[TString("m2")]  = pow(units["m"],2);
  
  //Volume units: derived from distance units
  units[TString("nm3")] = pow(units["nm"],3);
  units[TString("um4")] = pow(units["um"],3);
  units[TString("mm3")] = pow(units["mm"],3);
  units[TString("cm3")] = pow(units["cm"],3);
  units[TString("m3")]  = pow(units["m"],3);
  
  //Density units
  units[TString("gr/cm3")] = units["gr"]/units["cm3"];
  units[TString("Kg/cm3")] = units["Kg"]/units["cm3"];
  units[TString("Kg/m3")]  = units["Kg"]/units["m3"];
  
  //Rate density units: derived from the ones above
  units[TString("MHz/um2")]              = units["MHz"]/units["um2"];
  units[TString("MHz/mm2")]              = units["MHz"]/units["mm2"];
  units[TString("MHz/cm2")]              = units["MHz"]/units["cm2"];
  units[TString("GHz/um2")]              = units["GHz"]/units["um2"];
  units[TString("GHz/mm2")]              = units["GHz"]/units["mm2"];
  units[TString("GHz/cm2")]              = units["GHz"]/units["cm2"];
  units[TString("kHz/um2")]              = units["kHz"]/units["um2"];
  units[TString("kHz/mm2")]              = units["kHz"]/units["mm2"];
  units[TString("kHz/cm2")]              = units["kHz"]/units["cm2"];
  units[TString("Hz/um2")]               = units["Hz"] /units["um2"];
  units[TString("Hz/mm2")]               = units["Hz"] /units["mm2"];
  units[TString("Hz/cm2")]               = units["Hz"] /units["cm2"];
  units_RateDensity[TString("MHz/um2")]  = units[TString("MHz/um2")];
  units_RateDensity[TString("MHz/mm2")]  = units[TString("MHz/mm2")];
  units_RateDensity[TString("MHz/cm2")]  = units[TString("MHz/cm2")];
  units_RateDensity[TString("GHz/um2")]  = units[TString("GHz/um2")];
  units_RateDensity[TString("GHz/mm2")]  = units[TString("GHz/mm2")];
  units_RateDensity[TString("GHz/cm2")]  = units[TString("GHz/cm2")];
  units_RateDensity[TString("kHz/um2")]  = units[TString("kHz/um2")];
  units_RateDensity[TString("kHz/mm2")]  = units[TString("kHz/mm2")];
  units_RateDensity[TString("kHz/cm2")]  = units[TString("kHz/cm2")];
  units_RateDensity[TString("Hz/um2")]   = units[TString("Hz/um2")];
  units_RateDensity[TString("Hz/mm2")]   = units[TString("Hz/mm2")];
  units_RateDensity[TString("Hz/cm2")]   = units[TString("Hz/cm2")];
  
  //Magnetic-Field units: reference is T (Tesla)
  units[TString("T")]           = 1.0;
  units[TString("mT")]          = 1.0e-3*units[TString("T")];
  units[TString("uT")]          = 1.0e-6*units[TString("T")];
  units[TString("kT")]          = 1.0e+3*units[TString("T")];
  units[TString("MT")]          = 1.0e+6*units[TString("T")];
  units[TString("G")]           = 1.0e-4*units[TString("T")];
  units[TString("uG")]          = 1.0e-6*units[TString("G")];
  units[TString("mG")]          = 1.0e-3*units[TString("G")];
  units[TString("kG")]          = 1.0e+3*units[TString("G")];
  units[TString("MG")]          = 1.0e+6*units[TString("G")];
  units_Bfield[TString("T")]    = units[TString("T")];
  units_Bfield[TString("mT")]   = units[TString("mT")];
  units_Bfield[TString("uT")]   = units[TString("uT")];
  units_Bfield[TString("kT")]   = units[TString("kT")];
  units_Bfield[TString("MT")]   = units[TString("MT")];
  units_Bfield[TString("G")]    = units[TString("G")];
  units_Bfield[TString("uG")]   = units[TString("uG")];
  units_Bfield[TString("mG")]   = units[TString("mG")];
  units_Bfield[TString("kG")]   = units[TString("kG")];
  units_Bfield[TString("MG")]   = units[TString("MG")];

  
  //////////////////////////////////
  //// Particle properties block ///
  //////////////////////////////////
  Particle_t  aParticle;
  //e-
  aParticle.AntiParticle      = TString("e+");
  aParticle.charge            = -1;
  aParticle.Aweight           =  1;
  aParticle.mass              =  0.510999*units["MeV/c2"];
  ParticleMap[TString("e-")]  = aParticle;
  //e+
  aParticle.AntiParticle      =  TString("e-");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("e+")]  = aParticle;
  //mu+
  aParticle.AntiParticle      = TString("mu-");
  aParticle.charge            = +1;
  aParticle.Aweight           =  1;
  aParticle.mass              =  105.6584*units["MeV/c2"];
  ParticleMap[TString("mu+")] = aParticle;
  //mu-
  aParticle.AntiParticle      =  TString("mu+");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("mu-")] = aParticle;
  //pi+
  aParticle.AntiParticle      = TString("pi-");
  aParticle.charge            = +1;
  aParticle.Aweight           =  1;
  aParticle.mass              =  139.570*units["MeV/c2"];
  ParticleMap[TString("pi+")] = aParticle;
  //pi-
  aParticle.AntiParticle      =  TString("pi+");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("pi-")] = aParticle;
  //K+
  aParticle.AntiParticle      = TString("K-");
  aParticle.charge            = +1;
  aParticle.Aweight           =  1;
  aParticle.mass              =  493.68*units["MeV/c2"];
  ParticleMap[TString("K+")] = aParticle;
  //K-
  aParticle.AntiParticle      =  TString("K+");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("K-")] = aParticle;
  //p+
  aParticle.AntiParticle      = TString("p-");
  aParticle.charge            = +1;
  aParticle.Aweight           =  1;
  aParticle.mass              =  938.2723*units["MeV/c2"];
  ParticleMap[TString("p+")] = aParticle;
  //p-
  aParticle.AntiParticle      =  TString("p+");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("p-")] = aParticle;
  //4He+2
  aParticle.AntiParticle      = TString("4He-2");
  aParticle.charge            = +2;
  aParticle.Aweight           =  4;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("4He+2")] = aParticle;
  //4He--
  aParticle.AntiParticle      =  TString("4He+2");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("4He-2")] = aParticle;
  //7Li+3
  aParticle.AntiParticle      = TString("7Li-3");
  aParticle.charge            = +3;
  aParticle.Aweight           =  7;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("7Li+3")] = aParticle;
  //7Li-3
  aParticle.AntiParticle      =  TString("7Li+3");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("7Li-3")] = aParticle;
  //9Be+4
  aParticle.AntiParticle      = TString("9Be-4");
  aParticle.charge            = +4;
  aParticle.Aweight           =  9;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("9Be+4")] = aParticle;
  //9Be-4
  aParticle.AntiParticle      =  TString("9Be+4");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("9Be-4")] = aParticle;
  //10B+5
  aParticle.AntiParticle      = TString("10B-5");
  aParticle.charge            = +5;
  aParticle.Aweight           =  10;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("10B+5")] = aParticle;
  //10B-5
  aParticle.AntiParticle      =  TString("10B+5");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("10B-5")] = aParticle;
  //12C+6
  aParticle.AntiParticle      = TString("12C-6");
  aParticle.charge            = +6;
  aParticle.Aweight           =  12;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("12C+6")] = aParticle;
  //12C-6
  aParticle.AntiParticle      =  TString("12C+6");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("12C-6")] = aParticle;
  //14N+7
  aParticle.AntiParticle      = TString("14N-7");
  aParticle.charge            = +7;
  aParticle.Aweight           =  14;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("14N+7")] = aParticle;
  //14N-7
  aParticle.AntiParticle      =  TString("14N+7");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("14N-7")] = aParticle;
  //16O+8
  aParticle.AntiParticle      = TString("16O-8");
  aParticle.charge            = +8;
  aParticle.Aweight           =  16;
  aParticle.mass              =  aParticle.Aweight*ParticleMap[TString("p+")].mass;
  ParticleMap[TString("16O+8")] = aParticle;
  //16O-8
  aParticle.AntiParticle      =  TString("16O+8");
  aParticle.charge            = -ParticleMap[aParticle.AntiParticle].charge;
  aParticle.Aweight           =  ParticleMap[aParticle.AntiParticle].Aweight;
  aParticle.mass              =  ParticleMap[aParticle.AntiParticle].mass;
  ParticleMap[TString("16O-8")] = aParticle;
  
  //////////////////////////////////
  //// Material properties block ///
  //////////////////////////////////
  Material_t aMaterial;
  TString aMaterialName;
  
  //copper
  aMaterialName                = TString("copper");
  aMaterial.X0                 = 1.44*units["cm"];
  aMaterial.density            = 8.96*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 29;
  aMaterial.mZA                = aMaterial.mZ/63.55;
  aMaterial.mI                 = 322.0*units["eV"];
  aMaterial.color              = kOrange+10;
  MaterialMap[aMaterialName]   = aMaterial;  
  //silicon
  aMaterialName                = TString("silicon");
  aMaterial.X0                 = 9.36*units["cm"];
  aMaterial.density            = 2.3290*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 14;
  aMaterial.mZA                = aMaterial.mZ/28.09;
  aMaterial.mI                 = 173.0*units["eV"];
  aMaterial.color              = kOrange-3;
  MaterialMap[aMaterialName]   = aMaterial;
  //aluminum
  aMaterialName                = TString("aluminum");
  aMaterial.X0                 = 8.89*units["cm"];
  aMaterial.density            = 2.70*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 13;
  aMaterial.mZA                = aMaterial.mZ/26.98;
  aMaterial.mI                 = 166.0*units["eV"];
  aMaterial.color              = kBlue+2;
  MaterialMap[aMaterialName]   = aMaterial;
  //beryllium
  aMaterialName                = TString("beryllium");
  aMaterial.X0                 = 35.28*units["cm"];
  aMaterial.density            = 1.85*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 4;
  aMaterial.mZA                = aMaterial.mZ/9.01;
  aMaterial.mI                 = 63.7*units["eV"];
  aMaterial.color              = kGreen+1;
  MaterialMap[aMaterialName]   = aMaterial;
  //epoxy
  //Use molecular formula: C_21 H_25 Cl O_5
  aMaterialName                = TString("epoxy");
  aMaterial.X0                 = 35.70*units["cm"];
  aMaterial.density            = 1.16*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 6.10;
  aMaterial.mZA                = 0.5;
  aMaterial.mI                 = aMaterial.mZ*16.0*units["eV"];
  aMaterial.color              = kCyan-10;
  MaterialMap[aMaterialName]   = aMaterial;
  //kapton
  //Use the following composition: f_H = 0.026362; f_C = 0.691133; f_N = 0.073270; f_O = 0.209235
  aMaterialName                = TString("kapton");
  aMaterial.X0                 = 28.57*units["cm"];
  aMaterial.density            = 1.42*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 5.02;
  aMaterial.mZA                = 0.51264;
  aMaterial.mI                 = 79.6*units["eV"];
  aMaterial.color              = kOrange+3;
  MaterialMap[aMaterialName]   = aMaterial;
  //styropor
  //Use the following composition: f_H = 0.077418; f_C = 0.922582;
  aMaterialName                = TString("styropor");
  aMaterial.X0                 = 175.183*units["cm"];
  aMaterial.density            = 1.06*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 3.5;
  aMaterial.mZA                = 0.53768;
  aMaterial.mI                 = 68.7*units["eV"];
  aMaterial.color              = kAzure-9;
  MaterialMap[aMaterialName]   = aMaterial;
  //glass
  //Use the following composition: f_B = 0.040064; f_O = 0.539562; f_N = 0.028191; f_Al = 0.011644; f_Si = 0.377220; f_K = 0.003321;
  aMaterialName                = TString("glass");
  aMaterial.X0                 = 10.69*units["cm"];
  aMaterial.density            = 2.40*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 10.0;
  aMaterial.mZA                = 0.49731;
  aMaterial.mI                 = 145.4*units["eV"];
  aMaterial.color              = kCyan-8;
  MaterialMap[aMaterialName]   = aMaterial;
  //nai
  //Use the following composition: f_Na = 0.153373; f_I = 0.846627;
  aMaterialName                = TString("nai");
  aMaterial.X0                 = 2.59*units["cm"];
  aMaterial.density            = 3.67*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 32;
  aMaterial.mZA                = 0.42697;
  aMaterial.mI                 = 452.0*units["eV"];
  aMaterial.color              = kTeal+9;
  MaterialMap[aMaterialName]   = aMaterial;
  //csi
  //Use the following composition: f_Cs = 0.511549; f_I = 0.488451;
  aMaterialName                = TString("csi");
  aMaterial.X0                 = 1.86*units["cm"];
  aMaterial.density            = 4.51*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 54;
  aMaterial.mZA                = 0.41569;
  aMaterial.mI                 = 553.1*units["eV"];
  aMaterial.color              = kGreen+3;
  MaterialMap[aMaterialName]   = aMaterial;
  //tungsten
  aMaterialName                = TString("tungsten");
  aMaterial.X0                 = 0.35*units["cm"];
  aMaterial.density            = 19.30*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 74;
  aMaterial.mZA                = aMaterial.mZ/183.84;
  aMaterial.mI                 = 727.0*units["eV"];
  aMaterial.color              = kPink+8;
  MaterialMap[aMaterialName]   = aMaterial;
  //iron
  aMaterialName                = TString("iron");
  aMaterial.X0                 = 1.76*units["cm"];
  aMaterial.density            = 7.874*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 26;
  aMaterial.mZA                = aMaterial.mZ/55.85;
  aMaterial.mI                 = 286.0*units["eV"];
  aMaterial.color              = kGray+1;
  MaterialMap[aMaterialName]   = aMaterial;
  //water
  //Use the following composition: f_H = 0.111894; f_O = 0.888106;
  aMaterialName                = TString("water");
  aMaterial.X0                 = 36.08*units["cm"];
  aMaterial.density            = 1.00*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 3.33;
  aMaterial.mZA                = 0.55509;
  aMaterial.mI                 = 79.7*units["eV"];
  aMaterial.color              = kAzure-4;
  MaterialMap[aMaterialName]   = aMaterial;
  //mylar
  //Use the following composition: f_H = 0.143711; f_C = 0.856289;
  aMaterialName                = TString("mylar");
  aMaterial.X0                 = 28.54*units["cm"];
  aMaterial.density            = 1.40*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 4.55;
  aMaterial.mZA                = 0.52037;
  aMaterial.mI                 = 57.4*units["eV"];
  aMaterial.color              = kGray;
  MaterialMap[aMaterialName]   = aMaterial;
  //CarbonFiber
  aMaterialName                = TString("CarbonFiber");
  aMaterial.X0                 = 18.85*units["cm"];
  aMaterial.density            = 2.27*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 6;
  aMaterial.mZA                = aMaterial.mZ/12.01;
  aMaterial.mI                 = 78.0*units["eV"];
  aMaterial.color              = kGray+3;
  MaterialMap[aMaterialName]   = aMaterial;  
  //BeampipeBeCableMix
  //Look for it
  aMaterialName                = TString("BeampipeBeCableMix");
  aMaterial.X0                 = 9.87281*units["cm"];
  aMaterial.density            = 3.69*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 8.85;
  aMaterial.mZA                = 0.5;
  aMaterial.mI                 = 115.05*units["eV"];
  aMaterial.color              = kGreen+3;
  MaterialMap[aMaterialName]   = aMaterial;
  //DryAir
  //Use the following composition: f_C = 0.000124; f_N = 0.755267; f_O = 0.231781; f_Ar = 0.012827;
  aMaterialName                = TString("DryAir");
  aMaterial.X0                 = 303.90*units["m"];
  aMaterial.density            = 1.225e-3*units["gr/cm3"];
  aMaterial.density_effect_1st = 2.0;
  aMaterial.density_effect_2nd = 4.0;
  aMaterial.mZ                 = 7.34;
  aMaterial.mZA                = 0.49919;
  aMaterial.mI                 = 79.7*units["eV"];
  aMaterial.color              = kAzure-4;
  MaterialMap[aMaterialName]   = aMaterial;
  //CarbonFiber
  aMaterialName                = TString("CarbonFiber");
  aMaterial.X0                 = 18.85*units["cm"];
  aMaterial.density            = 2.27*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 6;
  aMaterial.mZA                = aMaterial.mZ/12.01;
  aMaterial.mI                 = 78.0*units["eV"];
  aMaterial.color              = kGray+3;
  MaterialMap[aMaterialName]   = aMaterial;  
  //BeampipeBeCableMix
  //Look for it
  aMaterialName                = TString("BeampipeBeCableMix");
  aMaterial.X0                 = 9.87281*units["cm"];
  aMaterial.density            = 3.69*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = 8.85;
  aMaterial.mZA                = 0.5;
  aMaterial.mI                 = 115.05*units["eV"];
  aMaterial.color              = kGreen+3;
  MaterialMap[aMaterialName]   = aMaterial;
  //DryAir
  //Use the following composition: f_C = 0.000124; f_N = 0.755267; f_O = 0.231781; f_Ar = 0.012827;
  aMaterialName                = TString("DryAir");
  aMaterial.X0                 = 303.90*units["m"];
  aMaterial.density            = 1.225e-3*units["gr/cm3"];
  aMaterial.density_effect_1st = 2.0;
  aMaterial.density_effect_2nd = 4.0;
  aMaterial.mZ                 = 7.34;
  aMaterial.mZA                = 0.49919;
  aMaterial.mI                 = 85.7*units["eV"];
  aMaterial.color              = kWhite;
  MaterialMap[aMaterialName]   = aMaterial;
  //TPCGAS
  //Use the following composition: Ar:CF_4:Iso-C_4H_10 => 95:3:2 %
  //Use the following composition: f_Ar = 0.95; f_C = 0.0206426; f_F = 0.0259091; f_H = 0.00344828;
  aMaterialName                = TString("TPCGAS");
  aMaterial.X0                 = 115.528*units["m"];
  aMaterial.density            = 1.731e-3*units["gr/cm3"];
  aMaterial.density_effect_1st = 2.0;
  aMaterial.density_effect_2nd = 4.0;
  aMaterial.mZ                 = 17.5;
  aMaterial.mZA                = 0.45526;
  aMaterial.mI                 = 177.4*units["eV"];
  aMaterial.color              = kAzure+1;
  MaterialMap[aMaterialName]   = aMaterial;
  //g10
  //Use the following composition: f_Cl = 0.080; f_O = 0.461; f_Si = 0.361; f_H = 0.019; f_C = 0.079;
  aMaterialName                = TString("g10");
  aMaterial.X0                 = 46.82*units["cm"];
  aMaterial.density            = 1.7*units["gr/cm3"];
  aMaterial.density_effect_1st = 0.2;
  aMaterial.density_effect_2nd = 3.0;
  aMaterial.mZ                 = (17*0.080 + 8*0.461 + 14*0.361 + 1*0.019 + 6*0.079)/(0.080 + 0.461 + 0.361 + 0.019 + 0.079);
  aMaterial.mZA                = ((17/35.45)*0.080 + (8/16.0)*0.461 + (14/28.09)*0.361 + (1/1.008)*0.019 + (6/12.0)*0.079)/(0.080 + 0.461 + 0.361 + 0.019 + 0.079);
  aMaterial.mI                 = 114.4*units["eV"];
  aMaterial.color              = kOrange-1;
  MaterialMap[aMaterialName]   = aMaterial;
  //Vacuum
  //Use the following composition (as Dry air): f_C = 0.000124; f_N = 0.755267; f_O = 0.231781; f_Ar = 0.012827;
  aMaterialName                = TString("Vacuum");
  aMaterial.X0                 = 1.00e+14*units["m"];
  aMaterial.density            = 1.00e-14*units["gr/cm3"];
  aMaterial.density_effect_1st = 2.0;
  aMaterial.density_effect_2nd = 4.0;
  aMaterial.mZ                 = 7.34;
  aMaterial.mZA                = 0.49919;
  aMaterial.mI                 = 85.7*units["eV"];
  aMaterial.color              = kWhite;
  MaterialMap[aMaterialName]   = aMaterial;
  
  //World Volume Types
  WorldVolumeTypes.clear();
  WorldVolumeTypes.push_back("Cylinder");
  WorldVolumeTypes.push_back("Box");
  
  //Geometry elements types table
  GeoPlaneTypes.clear();
  GeoPlaneTypes.push_back("Plane");
  GeoPlaneTypes.push_back("Cylinder");
  GeoPlaneTypes.push_back("Disk");
  GeoPlaneTypes.push_back("Cone");
  GeoPlaneTypes.push_back("Petal");
  GeoLadderTypes.clear();
  GeoLadderTypes.push_back("Plane");
  GeoLadderTypes.push_back("Cylinder");
  GeoLadderTypes.push_back("Disk");
  GeoLadderTypes.push_back("Spiral");
  GeoLadderTypes.push_back("Alternating");
  GeoLadderTypes.push_back("Cone");
  GeoLadderTypes.push_back("Petal");
  ResolutionModelTypes.clear();
  ResolutionModelTypes.push_back("TPC");
  
  //Pattern recognition algorithms list
  PatternRecogAlgoTypes.clear();
  PatternRecogAlgoTypes.push_back(TString("FPCCDTrkFinder"));
  PatternRecogAlgoTypes.push_back(TString("DBDILDTrkFinder"));
  
  //Factor for the calculation of the curvature radius: Pt(GeV/c) = betaCurvRadius*B(T)*R(mm)
  betaCurvRadius = 0.3*(GetUnit("GeV/c")/(GetUnit("T")*GetUnit("m")));
  cout << endl;
  cout << "betaCurvRadius = " << betaCurvRadius << " " << Reference_PEM_unit.Data() << "/(" << Reference_dist_unit.Data() << "*T)" << endl;
  cout << endl;
  
  epsilon_Bfield = 1.0e-6*GetUnit("T");
  
  epsilon_distance = 1.0*GetUnit("nm");
    
  return;
  
}
//====================================================================
bool GGlobalTools::CheckParticleName(TString particle)
{
  
  //Checking if particle's name is in internal table
#if 0
  if(mass.count(particle) < 1) {
    cout << "particle type " << particle.Data() << " is unknown!" << endl;
    cout << "particle should be: " << endl;
    for(map<TString,double>::iterator i=mass.begin(); i!=mass.end(); ++i) cout << " - " << ((*i).first).Data() << " " << endl;
    
    return false;
  }
#endif

  bool IsInList = false;
  for(map<TString,Particle_t>::iterator i=ParticleMap.begin(); i!=ParticleMap.end(); ++i) {
    if(particle == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  return IsInList;
  
}
//====================================================================
double GGlobalTools::GetX0FromMaterial(TString material)
{
  
  //Returns the material X0 from material name
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].X0;
  
}
//====================================================================
double GGlobalTools::GetDensityFromMaterial(TString material)
{
  
  //Returns the material density from material name
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].density;
  
}
//====================================================================
int  GGlobalTools::GetMaterialColor(TString material)
{
  
  //Returns the material color from material name
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].color;
  
}
//====================================================================
double  GGlobalTools::GetDensityEffect1stFromMaterial(TString material)
{
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].density_effect_1st;
  
}
//====================================================================
double  GGlobalTools::GetDensityEffect2ndFromMaterial(TString material)
{
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].density_effect_2nd;
  
}
//====================================================================
double  GGlobalTools::GetmIFromMaterial(TString material)
{
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].mI;
  
}
//====================================================================
double  GGlobalTools::GetmAZFromMaterial(TString material)
{
  
  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return MaterialMap[material].mZA;
  
}
//====================================================================
double GGlobalTools::GetParticleMass(TString particle)
{
  
  //Returns the particle's mass from name
  
  bool IsInList = false;
  for(map<TString,Particle_t>::iterator i=ParticleMap.begin(); i!=ParticleMap.end(); ++i) {
    if(particle == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "particle type " << particle.Data() << " is unknown!" << endl;
    cout << "particle should be: " << endl;
    for(map<TString,Particle_t>::iterator ii=ParticleMap.begin(); ii!=ParticleMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return  ParticleMap[particle].mass;
  
}
//====================================================================
int GGlobalTools::GetParticleCharge(TString particle)
{
  
  //Returns the particle's charge (in electron charge units) from name
  
  bool IsInList = false;
  for(map<TString,Particle_t>::iterator i=ParticleMap.begin(); i!=ParticleMap.end(); ++i) {
    if(particle == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "particle type " << particle.Data() << " is unknown!" << endl;
    cout << "particle should be: " << endl;
    for(map<TString,Particle_t>::iterator ii=ParticleMap.begin(); ii!=ParticleMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return  ParticleMap[particle].charge;
  
}
//====================================================================
int GGlobalTools::GetParticleAweight(TString particle)
{
  
  //Returns the particle's atomic weight from name. Mainly useful in cases of heavy ions

  bool IsInList = false;
  for(map<TString,Particle_t>::iterator i=ParticleMap.begin(); i!=ParticleMap.end(); ++i) {
    if(particle == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "particle type " << particle.Data() << " is unknown!" << endl;
    cout << "particle should be: " << endl;
    for(map<TString,Particle_t>::iterator ii=ParticleMap.begin(); ii!=ParticleMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return  ParticleMap[particle].Aweight;
  
}
//====================================================================
TString  GGlobalTools::GetAntiParticle(TString particle)
{
  
  //Returns the particle's antiparticle name from particle's name

  bool IsInList = false;
  for(map<TString,Particle_t>::iterator i=ParticleMap.begin(); i!=ParticleMap.end(); ++i) {
    if(particle == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "particle type " << particle.Data() << " is unknown!" << endl;
    cout << "particle should be: " << endl;
    for(map<TString,Particle_t>::iterator ii=ParticleMap.begin(); ii!=ParticleMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  return  ParticleMap[particle].AntiParticle;
  
}
//====================================================================
double GGlobalTools::GetUnit(const char* myunit)
{
  
  //Returns the value of a given mynuit in terms of the reference units
  
  if(units.count(TString(myunit)) < 1) {
    cout << "unit type " << myunit << " is unknown!" << endl;
    cout << "unit should be: " << endl;

    for(map<TString,double>::iterator i=units.begin(); i!=units.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetUnit(TString myunit)
{
  
  //Returns the value of a given mynuit in terms of the reference units
  
  if(units.count(myunit) < 1) {
    cout << "unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "unit should be: " << endl;

    for(map<TString,double>::iterator i=units.begin(); i!=units.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units[myunit];
  
}
//====================================================================
double GGlobalTools::GetDistanceUnit(const char* myunit)
{
  
  //Returns the value of a given distance mynuit in terms of the reference units
  
  if(units_distance.count(TString(myunit)) < 1) {
    cout << "distance unit type " << myunit << " is unknown!" << endl;
    cout << "distance unit should be: " << endl;

    for(map<TString,double>::iterator i=units_distance.begin(); i!=units_distance.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_distance[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetDistanceUnit(TString myunit)
{
  
  //Returns the value of a given distance mynuit in terms of the reference units
  
  if(units_distance.count(myunit) < 1) {
    cout << "distance unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "distance unit should be: " << endl;

    for(map<TString,double>::iterator i=units_distance.begin(); i!=units_distance.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_distance[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsDistanceUnit(const char* myunit)
{
  
  //Checks if myunit is a distance unit
  
  if(units_distance.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsDistanceUnit(TString myunit)
{
  
  //Checks if myunit is a distance unit
  
  if(units_distance.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetTimeUnit(const char* myunit)
{
  
  //Returns the value of a given time mynuit in terms of the reference units
  
  if(units_time.count(TString(myunit)) < 1) {
    cout << "time unit type " << myunit << " is unknown!" << endl;
    cout << "time unit should be: " << endl;

    for(map<TString,double>::iterator i=units_time.begin(); i!=units_time.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_time[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetTimeUnit(TString myunit)
{
  
  //Returns the value of a given time mynuit in terms of the reference units
  
  if(units_time.count(myunit) < 1) {
    cout << "time unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "time unit should be: " << endl;

    for(map<TString,double>::iterator i=units_time.begin(); i!=units_time.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_time[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsTimeUnit(const char* myunit)
{
  
  //Checks if myunit is a time unit
  
  if(units_time.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsTimeUnit(TString myunit)
{
  
  //Checks if myunit is a time unit
  
  if(units_time.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetAngleUnit(const char* myunit)
{
  
  //Returns the value of a given angle mynuit in terms of the reference units
  
  if(units_angle.count(TString(myunit)) < 1) {
    cout << "angle unit type " << myunit << " is unknown!" << endl;
    cout << "angle unit should be: " << endl;

    for(map<TString,double>::iterator i=units_angle.begin(); i!=units_angle.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_angle[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetAngleUnit(TString myunit)
{
  
  //Returns the value of a given angle mynuit in terms of the reference units
  
  if(units_angle.count(myunit) < 1) {
    cout << "angle unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "angle unit should be: " << endl;

    for(map<TString,double>::iterator i=units_angle.begin(); i!=units_angle.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_angle[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsAngleUnit(const char* myunit)
{
  
  //Checks if myunit is a angle unit
  
  if(units_angle.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsAngleUnit(TString myunit)
{
  
  //Checks if myunit is a angle unit
  
  if(units_angle.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetEnergyUnit(const char* myunit)
{
  
  //Returns the value of a given Energy mynuit in terms of the reference units
  
  if(units_Energy.count(TString(myunit)) < 1) {
    cout << "Energy unit type " << myunit << " is unknown!" << endl;
    cout << "Energy unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Energy.begin(); i!=units_Energy.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Energy[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetEnergyUnit(TString myunit)
{
  
  //Returns the value of a given Energy mynuit in terms of the reference units
  
  if(units_Energy.count(myunit) < 1) {
    cout << "Energy unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "Energy unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Energy.begin(); i!=units_Energy.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Energy[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsEnergyUnit(const char* myunit)
{
  
  //Checks if myunit is a Energy unit
  
  if(units_Energy.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsEnergyUnit(TString myunit)
{
  
  //Checks if myunit is a Energy unit
  
  if(units_Energy.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetMomentumUnit(const char* myunit)
{
  
  //Returns the value of a given Momentum mynuit in terms of the reference units
  
  if(units_Momentum.count(TString(myunit)) < 1) {
    cout << "Momentum unit type " << myunit << " is unknown!" << endl;
    cout << "Momentum unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Momentum.begin(); i!=units_Momentum.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Momentum[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetMomentumUnit(TString myunit)
{
  
  //Returns the value of a given Momentum mynuit in terms of the reference units
  
  if(units_Momentum.count(myunit) < 1) {
    cout << "Momentum unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "Momentum unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Momentum.begin(); i!=units_Momentum.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Momentum[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsMomentumUnit(const char* myunit)
{
  
  //Checks if myunit is a Momentum unit
  
  if(units_Momentum.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsMomentumUnit(TString myunit)
{
  
  //Checks if myunit is a Momentum unit
  
  if(units_Momentum.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetMassUnit(const char* myunit)
{
  
  //Returns the value of a given Mass mynuit in terms of the reference units
  
  if(units_Mass.count(TString(myunit)) < 1) {
    cout << "Mass unit type " << myunit << " is unknown!" << endl;
    cout << "Mass unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Mass.begin(); i!=units_Mass.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Mass[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetMassUnit(TString myunit)
{
  
  //Returns the value of a given Mass mynuit in terms of the reference units
  
  if(units_Mass.count(myunit) < 1) {
    cout << "Mass unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "Mass unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Mass.begin(); i!=units_Mass.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Mass[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsMassUnit(const char* myunit)
{
  
  //Checks if myunit is a Mass unit
  
  if(units_Mass.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsMassUnit(TString myunit)
{
  
  //Checks if myunit is a Mass unit
  
  if(units_Mass.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetBfieldUnit(const char* myunit)
{
  
  //Returns the value of a given Bfield mynuit in terms of the reference units
  
  if(units_Bfield.count(TString(myunit)) < 1) {
    cout << "Bfield unit type " << myunit << " is unknown!" << endl;
    cout << "Bfield unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Bfield.begin(); i!=units_Bfield.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Bfield[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetBfieldUnit(TString myunit)
{
  
  //Returns the value of a given Bfield mynuit in terms of the reference units
  
  if(units_Bfield.count(myunit) < 1) {
    cout << "Bfield unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "Bfield unit should be: " << endl;

    for(map<TString,double>::iterator i=units_Bfield.begin(); i!=units_Bfield.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_Bfield[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsBfieldUnit(const char* myunit)
{
  
  //Checks if myunit is a Bfield unit
  
  if(units_Bfield.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsBfieldUnit(TString myunit)
{
  
  //Checks if myunit is a Bfield unit
  
  if(units_Bfield.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetRateDensityUnit(const char* myunit)
{
  
  //Returns the value of a given RateDensity mynuit in terms of the reference units
  
  if(units_RateDensity.count(TString(myunit)) < 1) {
    cout << "RateDensity unit type " << myunit << " is unknown!" << endl;
    cout << "RateDensity unit should be: " << endl;

    for(map<TString,double>::iterator i=units_RateDensity.begin(); i!=units_RateDensity.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_RateDensity[TString(myunit)];
  
}
//====================================================================
double GGlobalTools::GetRateDensityUnit(TString myunit)
{
  
  //Returns the value of a given RateDensity mynuit in terms of the reference units
  
  if(units_RateDensity.count(myunit) < 1) {
    cout << "RateDensity unit type " << myunit.Data() << " is unknown!" << endl;
    cout << "RateDensity unit should be: " << endl;

    for(map<TString,double>::iterator i=units_RateDensity.begin(); i!=units_RateDensity.end(); ++i) {
      cout << " - " << ((*i).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);

  }
  
  return units_RateDensity[myunit];
  
}
//====================================================================
bool  GGlobalTools::IsRateDensityUnit(const char* myunit)
{
  
  //Checks if myunit is a RateDensity unit
  
  if(units_RateDensity.count(TString(myunit)) < 1) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsRateDensityUnit(TString myunit)
{
  
  //Checks if myunit is a RateDensity unit
  
  if(units_RateDensity.count(myunit) < 1) return false;
  
  return true;
  
}
//====================================================================
double GGlobalTools::GetAngle(double x,double y)
{
  
  //Returns polar angle in (0.0,2*Pi) range for a given set of x and y coordinates
  
  double angle = 0.0;
  
  double epsilon = 1.0e-10;
  if(TMath::Abs(x) < epsilon) {
    if(TMath::Abs(y) < epsilon) angle = 0.0;
    else if(y >  epsilon)       angle = +0.5*TMath::Pi();
    else if(y < -epsilon)       angle = +1.5*TMath::Pi();
  }
  else if(x > epsilon) {
    if(TMath::Abs(y) < epsilon) angle = 0.0;
    else if(y >  epsilon)       angle = +TMath::ATan(TMath::Abs(y/x));
    else if(y < -epsilon)       angle = 2.0*TMath::Pi() - TMath::ATan(TMath::Abs(y/x));
  }
  else if(x < -epsilon) {
    if(TMath::Abs(y) < epsilon) angle = +TMath::Pi();
    else if(y >  epsilon)       angle = +TMath::Pi() - TMath::ATan(TMath::Abs(y/x));
    else if(y < -epsilon)       angle = +TMath::Pi() + TMath::ATan(TMath::Abs(y/x));
  }
  
  return angle;
  
}
//====================================================================
int  GGlobalTools::delta_cronequer(int i,int j)
{
  
  //Delta of cronequer function
  
  if(i==j) return 1;
  else     return 0;
  
}
//====================================================================
bool GGlobalTools::CheckIfGoodSquareMatrix(TMatrixT<double> M)
{
  
  //Checks if matrix M is a squarre matrix
  
  int Ncols = M.GetNcols();
  int Nrows = M.GetNrows();
  
  if(Ncols != Nrows) return false;
  if(Ncols <= 0)     return false;
  
  return true;
  
}
//====================================================================
bool GGlobalTools::CheckUnitMatrix(TMatrixT<double> M)
{
  
  //Checks if matrix M is unit matrix
  
  if(!CheckIfGoodSquareMatrix(M)) return false;
  
  int Ncols = M.GetNcols();
  for(int icol=0;icol<Ncols;icol++) {
    for(int irow=0;irow<Ncols;irow++) {
      if(icol == irow) {
	if(TMath::Abs(M(icol,irow) - 1.0) > 1.0e-8) return false;
      }
      else {
	if(TMath::Abs(M(icol,irow) - 0.0) > 1.0e-8) return false;
      }
    }
  }
  
  return true;
  
}
//====================================================================
bool GGlobalTools::CheckIfSymmetricMatrix(TMatrixT<double> M)
{
 
  //Checks if matrix M is symmetric
  
  if(!CheckIfGoodSquareMatrix(M)) return false;
  
  int Ncols = M.GetNcols();
  for(int icol=0;icol<Ncols;icol++) {
    for(int irow=0;irow<Ncols;irow++) {
      if(icol >= irow) continue;
      if(TMath::Abs(M(icol,irow) - M(irow,icol)) > 1.0e-8) {
	return false;
      }
    }
  }
  
  return true;
  
}
//====================================================================
bool GGlobalTools::CheckIfGoodErrorMatrix(TMatrixT<double> M)
{
  
  //Checks if matrix M is a good covariance matrix by checking if the correlation coefficients are inside [-1,1] interval
  
  if(!CheckIfGoodSquareMatrix(M)) return false;
  if(!CheckIfSymmetricMatrix(M))  return false;
  
  int Ncols = M.GetNcols();
  for(int icol=0;icol<Ncols;icol++) {
    if(M(icol,icol) < 0.0) return false;
  }
  
  for(int icol=0;icol<Ncols;icol++) {
    for(int irow=0;irow<Ncols;irow++) {
      if(icol >= irow) continue;
      double AbsCorr = TMath::Abs(M(icol,irow)/sqrt(M(icol,icol)*M(irow,irow)));
      if(AbsCorr - 1.0 > 1.0e-8) return false;
    }
  }
  
  return true;
  
}
//====================================================================
void GGlobalTools::GetGeometryImpacParameterParameters(TGraphErrors* gr,
						       double  theta,
						       TString sigma_units,
						       TString momentum_units,
						       bool    FitPowerForImpactParam,
						       double &a,
						       double &b,
						       double &exp)
{
  
  // Get a and b parameters from sigma_drho vs momentum for a given theta angle from the parameterization
  // sigma_ip^2 = a^2 + (b/p*beta*pow(sin(theta),3/2)^2,
  // a and b parameters are obtained by fitting the graph gr
  
  double a_init,b_init;
  const int Npoints(gr->GetN());
  
  double p;
  double sigma;
  double Rp[2];
  double Rsigma[2];
  
  //Getting the 1st point in the curve
  std::vector<int> idx_list;
  idx_list.clear();
  Rsigma[1] = -1.0e+20;
  for(int i=0;i<gr->GetN();i++) {
    if(i+1 < gr->GetN()-2) {
      double p1,p2,p3;
      double sigma1,sigma2,sigma3;
      gr->GetPoint(i,  p1,sigma1);
      gr->GetPoint(i+1,p2,sigma2);
      gr->GetPoint(i+2,p3,sigma3);
      if(sigma2 > sigma1 && sigma2 > sigma3) {
	idx_list.push_back(i+1);
      }
    }
    
    gr->GetPoint(i,p,sigma);
    if(Rsigma[1] < sigma) {
      Rsigma[1] = sigma;
      Rp[0]     = p;
    }
  }
  if(idx_list.size() > 0) {
    Rsigma[1] = -1.0e+20;
    for(int i=0;i<int(idx_list.size());i++) {
      gr->GetPoint(idx_list[i],p,sigma);
      //cout << "Local minuma found at p = " << p << " with value " << sigma << endl;
      if(Rsigma[1] < sigma) {
        Rsigma[1] = sigma;
        Rp[0]     = p;
      }
    }
  }
  p     = Rp[0];
  sigma = Rsigma[1];
  //Initial estimate b using 1st point
  b_init = sigma*p*pow(TMath::Sin(theta),1.5);
  
  //Getting the last point in the curve
  gr->GetPoint(Npoints-1, p, sigma);
  Rp[1]     = p;
  Rsigma[0] = sigma;
  //Initial estimate a using last point
  a_init = sigma;
  
  if(Rsigma[0] > Rsigma[1]) {
    double tmp_sigma = Rsigma[0];
    Rsigma[0] = Rsigma[1];
    Rsigma[1] = tmp_sigma;
  }
  
  //Get the point where the function decreases by a given porcent w.r.t to the 1st point
  double porcent = 0.10;
  //porcent = 0.10;
  double Rp_max = -999.0;
  for(int i=0;i<Npoints;i++) {
    if(i==0) continue;
    gr->GetPoint(i, p, sigma);
    
    if(sigma < Rsigma[1]*porcent) {
      Rp_max = p;
      break;
    }
  }
  if(Npoints == 1)     Rp_max = Rp[1];
  if(Rp_max == -999.0) Rp_max = Rp[1];
  
  //Write the fit functions
  TF1* f_impactPar = new TF1("f_impactPar","sqrt( pow([0],2) + pow([1]/(x*pow(TMath::Sin([3]),[2])),2) )",Rp[0],Rp_max);
  f_impactPar->SetParameter(0,a_init);
  f_impactPar->SetParameter(1,b_init);
  f_impactPar->FixParameter(3,theta);
  f_impactPar->SetParLimits(0,0.0,a_init*20);
  f_impactPar->SetParLimits(1,0.0,b_init*20);
  
  if(FitPowerForImpactParam) {
    f_impactPar->SetParameter(2,3.0/2.0);
    f_impactPar->SetParLimits(2,0.0,2.0);
  }
  else f_impactPar->FixParameter(2,3.0/2.0);
  
  //Fitting the TGraph
  gr->Fit(f_impactPar,"WQRN0");
  //gr->Fit(f_impactPar,"WRN0");
  
  //Obtain the values of the a and b parameters
  a   = f_impactPar->GetParameter(0)*units[sigma_units];
  b   = f_impactPar->GetParameter(1)*units[sigma_units]*units[momentum_units];
  exp = f_impactPar->GetParameter(2);
  
  return;
  
}
//====================================================================
void  GGlobalTools::OrderIntersectionHitList(std::vector<IntersectionHit_t>& ItersectionHitList)
{
  
  //Order intersection hit elements in list wrt s parameter from low to high
  
  if(int(ItersectionHitList.size()) <= 1) return;
  
  for(int iii=2;iii<=int(ItersectionHitList.size());iii++) {
    for(int jjj=0;jjj<=int(ItersectionHitList.size())-iii;jjj++) {
      double s_jjj   = ItersectionHitList[jjj].s;
      double s_jjjp1 = ItersectionHitList[jjj+1].s;
      
      if(s_jjj > s_jjjp1) {
	IntersectionHit_t AHit;
        AHit.s                                     = ItersectionHitList[jjj].s;
	AHit.geoElement_idx                        = ItersectionHitList[jjj].geoElement_idx;
        AHit.IsSensitivePoint                      = ItersectionHitList[jjj].IsSensitivePoint;
	AHit.s_in                                  = ItersectionHitList[jjj].s_in;
	AHit.s_out                                 = ItersectionHitList[jjj].s_out;
	
	ItersectionHitList[jjj].s                  = ItersectionHitList[jjj+1].s;
	ItersectionHitList[jjj].geoElement_idx     = ItersectionHitList[jjj+1].geoElement_idx;
        ItersectionHitList[jjj].IsSensitivePoint   = ItersectionHitList[jjj+1].IsSensitivePoint;
	ItersectionHitList[jjj].s_in               = ItersectionHitList[jjj+1].s_in;
	ItersectionHitList[jjj].s_out              = ItersectionHitList[jjj+1].s_out;
	
	ItersectionHitList[jjj+1].s                = AHit.s;
	ItersectionHitList[jjj+1].geoElement_idx   = AHit.geoElement_idx;
        ItersectionHitList[jjj+1].IsSensitivePoint = AHit.IsSensitivePoint;
	ItersectionHitList[jjj+1].s_in             = AHit.s_in;
	ItersectionHitList[jjj+1].s_out            = AHit.s_out;
	
      }
    }
  }
  
  return;
  
}
//====================================================================
void  GGlobalTools::OrderListList(std::vector<double>& List)
{
  
  //Order list from low to high
  
  if(int(List.size()) <= 1) return;
  
  for(int iii=2;iii<=int(List.size());iii++) {
    for(int jjj=0;jjj<=int(List.size())-iii;jjj++) {
      double s_jjj   = List[jjj];
      double s_jjjp1 = List[jjj+1];
      
      if(s_jjj > s_jjjp1) {
	double aux  = List[jjj];
	List[jjj]   = List[jjj+1];
	List[jjj+1] = aux;
      }
    }
  }
  
  return;
  
}
//====================================================================
void  GGlobalTools::GetListOfLayersTypes(long Config, std::vector<int>& LayerTypes)
{
  
  // Get list of layers types encoded in a single long number (layer types in terms of hit type: missed hit, good association and fake association)
  
  if(Config <= 0) return;
  
  LayerTypes.clear();
  
  int order = 0;
  double ratio;
  
  do {
    order++;
    ratio  = Config;
    ratio /= pow(10,order);
  }
  while(ratio > 1.0);
  order -= 1;
  
  long number = Config;
  for(int i=0;i<order+1;i++) {
    int order_tmp = order - i;
    int digit  = int(number/pow(10,order_tmp));
    
    number -= long(pow(10,order_tmp))*digit;
    
    LayerTypes.push_back(digit);
  }
  
  return;
  
}
//====================================================================
long  GGlobalTools::GetConfigFromListOfLayersTypes(std::vector<int> LayerTypes) 
{
  
  //Convert list of layers type in a single long number
  
  long config = 0;
  
  int Ntypes = LayerTypes.size();
  for(int i=0;i<Ntypes;i++) {
    int order = Ntypes - i - 1;
    config += long(LayerTypes[i]*pow(10,order));
  }
  
  return config;
  
}
//====================================================================
void  GGlobalTools::IncludeFakesOnConfigList(long config,int Nfakes,bool IsSeeding,std::vector<long>& NewConfigList)
{
  
  NewConfigList.clear();
  
  if(Nfakes <= 0) return;
  
  int Nnull     = 0;
  int Nfake_ppp = 0;
  int Ngood     = 0;
  std::vector<int>  TheLayersTypeList;
  TheLayersTypeList.clear();
  GetListOfLayersTypes(config,TheLayersTypeList);
  
  std::vector<int>  NullPositions;
  NullPositions.clear();
  for(int ilayer=0;ilayer<int(TheLayersTypeList.size());ilayer++) {
    if(TheLayersTypeList[ilayer] == 1)      Ngood++;
    else if(TheLayersTypeList[ilayer] == 2) Nfake_ppp++;
    else {
      Nnull++;
      NullPositions.push_back(ilayer);
    }
  }
  
  int Nhits = Ngood + Nfake_ppp;
  if(Nhits <= 0) return;
  
  std::vector<long> TmpConfigList;
  TmpConfigList.clear();
  std::vector<int>  Position;
  for(int ifake=1;ifake<=Nfakes;ifake++) {
    int Nfake_tmp = ifake;
    int Ngood_hit = Nhits - Nfake_tmp;
    if(Ngood_hit < 0) continue;
    
    const int tmp_Nhits_nonseed(Nhits);
    int Config[tmp_Nhits_nonseed];
    int HitKind[tmp_Nhits_nonseed];
    for(int k=0;k<Nhits;k++) Config[k] = -1;
    int counter_kind = 0;
    for(int jgood=0;jgood<Ngood_hit;jgood++) {
      HitKind[counter_kind] = 1;
      counter_kind++;
    }
    for(int jfake=0;jfake<Nfake_tmp;jfake++) {
      HitKind[counter_kind] = 2;
      counter_kind++;
    }

    FillConfigurations(Nhits,Position,Nhits,Config,HitKind,TmpConfigList);
  }
  
  if(Nnull == 0) {
    for(int i=0;i<int(TmpConfigList.size());i++) NewConfigList.push_back(TmpConfigList[i]);
    
    return;
  }
  
  for(int i=0;i<int(TmpConfigList.size());i++) {
    std::vector<int>  LayersTypeList;
    LayersTypeList.clear();
    GetListOfLayersTypes(TmpConfigList[i],LayersTypeList);
    
    std::vector<int>  TmpLayersTypeList;
    TmpLayersTypeList.clear();
    int counter = 0;
    for(int ilayer=0;ilayer<int(TheLayersTypeList.size());ilayer++) {
      bool IsNullPosition = false;
      for(int kkk=0;kkk<int(NullPositions.size());kkk++) {
	if(ilayer == NullPositions[kkk]) {
	  IsNullPosition = true;
	  break;
	}
      }
      if(IsNullPosition) {
	if(IsSeeding) TmpLayersTypeList.push_back(0);
	else          TmpLayersTypeList.push_back(3);
	continue;
      }
      
      TmpLayersTypeList.push_back(LayersTypeList[counter]);
      counter++;
    }
    long  config_tmp = GetConfigFromListOfLayersTypes(TmpLayersTypeList);
    NewConfigList.push_back(config_tmp);
  }
  TmpConfigList.clear();
  
  return;
  
}
//====================================================================
void  GGlobalTools::GetNGoodFakesAndNull(long config,int& Ngood,int& Nfake,int& Nnull)
{
  
  Nnull = 0;
  Nfake = 0;
  Ngood = 0;
  std::vector<int>  LayersTypeList;
  LayersTypeList.clear();
  GetListOfLayersTypes(config,LayersTypeList);
  
  for(int ilayer=0;ilayer<int(LayersTypeList.size());ilayer++) {
    if(LayersTypeList[ilayer] == 1)      Ngood++;
    else if(LayersTypeList[ilayer] == 2) Nfake++;
    else                                 Nnull++;
  }
  
  return;
  
}
//====================================================================
void GGlobalTools::FillConfigurations(int depth,
				      std::vector<int> Position,
				      int  Nlevels,
				      int* Config,
				      int* HitKind,
				      std::vector<long>& ConfigList)
{
  
  // Fill all possible configurations of layer types
  
  if(depth > 0) {
    for(int i=0;i<Nlevels;i++) {
      std::vector<int> MyPos;
      MyPos.clear();
      
      if(depth < Nlevels) {
        bool IsAlreadyIn = false;
        for(int ipos=0;ipos<int(Position.size());ipos++) {
	  if(Position[ipos] == i) {
	    IsAlreadyIn = true;
	    break;
	  }
        }
        if(IsAlreadyIn) continue;
      }
      
      if(depth < Nlevels) {
	for(int ipos=0;ipos<int(Position.size());ipos++) MyPos.push_back(Position[ipos]);
      }
      MyPos.push_back(i);
      Config[i] = HitKind[Nlevels-depth];
      
      if(depth > 1) FillConfigurations(depth-1,MyPos,Nlevels,Config,HitKind,ConfigList);
      
      if(depth == 1) {
	long config = 0;
        for(int k=0;k<Nlevels;k++) {
	  int order = Nlevels - k - 1;
	  config += long(pow(10,order)*Config[k]);
        }
      
        bool IsAlreadyIn = false;
        for(int iconf=0;iconf<int(ConfigList.size());iconf++) {
	  if(ConfigList[iconf] == config) {
	    IsAlreadyIn = true;
	    break;
	  }
	}
	if(!IsAlreadyIn) ConfigList.push_back(config);
      }
      
    }
  }
  
  return;
  
}
//====================================================================
void  GGlobalTools::FillSeedConfigsWithFakes(std::vector<long>& SeedConfigList, int NfakesSeedMax)
{

  if(NfakesSeedMax < 1) return;
  
  return;
  
}
//====================================================================
void GGlobalTools::RotateVector(TMatrixD Rot,
				TVector3& vector)
{ 
  
  // Rotate vector with Rotation matrix Tor
  
  double vect_old[3];
  vect_old[0] = vector.X();
  vect_old[1] = vector.Y();
  vect_old[2] = vector.Z();
  double vect_new[3];
  
  for(int i=0;i<3;i++) {
    vect_new[i] = 0.0;
    for(int j=0;j<3;j++) vect_new[i] += Rot(i,j)*vect_old[j];
  }
  
  vector = TVector3(vect_new[0],vect_new[1],vect_new[2]);
  
  return;
  
}
//====================================================================
void GGlobalTools::GetGlobalRotationMatrix(TVector3 angles,TMatrixD& Rot)
{
  
  // Get a rotation matrix from the rotation angles in the order Z, Y, and X
  
  std::vector<TString> axisList;
  std::vector<double>  angleList;
  axisList.clear();
  angleList.clear();
  
  axisList.push_back(TString("Z")); angleList.push_back(angles.Z());
  axisList.push_back(TString("Y")); angleList.push_back(angles.Y());
  axisList.push_back(TString("X")); angleList.push_back(angles.X());
  
  GetGlobalRotationMatrix_FromList(axisList,angleList,Rot);

  return;
  
}
//====================================================================
void GGlobalTools::GetGlobalRotationMatrix_FromList(std::vector<TString> axis,
						    std::vector<double>  angles,
						    TMatrixD& Rot)
{
  
  // Get a total rotation matrix from a list of ration angles
  
  TMatrixD Rot_prev;
  Rot_prev.ResizeTo(3,3);
  Rot.ResizeTo(3,3);
  GetRotationMatrix(0.0,TString("X"),Rot_prev);
  for(int i=0;i<int(axis.size());i++) {
    TMatrixD Rot_tmp;
    Rot_tmp.ResizeTo(3,3);
    GetRotationMatrix(angles[i],axis[i],Rot_tmp);
    
    Rot      = Rot_tmp*Rot_prev;
    Rot_prev = Rot;
  }
  
  return;
  
}
//====================================================================
void GGlobalTools::GetRotationMatrix(double angle,
				     TString axis,
				     TMatrixD& Rot)
{
  
  // Get rotation matrix from a rotation angle w.r.t one of the X, Y or Z axes
  
  Rot.ResizeTo(3,3);
  
  double cos = TMath::Cos(angle);
  double sin = TMath::Sin(angle);
  
  if(axis == TString("X")) {
    Rot(0,0) =  1.0; Rot(0,1) =  0.0; Rot(0,2) =  0.0;
    Rot(1,0) =  0.0; Rot(1,1) =  cos; Rot(1,2) =  sin;
    Rot(2,0) =  0.0; Rot(2,1) = -sin; Rot(2,2) =  cos;
  }
  else if(axis == TString("Y")) {
    Rot(0,0) =  cos; Rot(0,1) =  0.0; Rot(0,2) = -sin;
    Rot(1,0) =  0.0; Rot(1,1) =  1.0; Rot(1,2) =  0.0;
    Rot(2,0) =  sin; Rot(2,1) =  0.0; Rot(2,2) =  cos;
  }
  else if(axis == TString("Z")) {
    Rot(0,0) =  cos; Rot(0,1) =  sin; Rot(0,2) =  0.0;
    Rot(1,0) = -sin; Rot(1,1) =  cos; Rot(1,2) =  0.0;
    Rot(2,0) =  0.0; Rot(2,1) =  0.0; Rot(2,2) =  1.0;
  }
  else {
    cout << endl;
    cout << "ERROR in GGlobalTools::GetRotationMatrix:: Specified axis " << axis.Data() << " is not either X, Y or Z. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
TVector3 GGlobalTools::GetRotationAngles(TMatrixD RotMatrix, bool Myverbose)
{
  
  // Decomposes the rotation matrix, that is:
  // compute the values of the rotation angle around z, y, x
  // from the values of the rotation matrix
  
  double thetaX,thetaY,thetaZ; //0=x, 1=y, 2=z
  
  double R33 = RotMatrix(2,2); //cosX cosy
  double R23 = RotMatrix(1,2); //sinX cosy
  
  double R11 = RotMatrix(0,0); //cosy cosz
  double R12 = RotMatrix(0,1); //cosy sinz
  
  double R13 = RotMatrix(0,2); //-siny
  
  thetaX = atan2(R23,R33);
  thetaY = atan2(-R13,sqrt(pow(R11,2) + pow(R12,2)));
  thetaZ = atan2(R12,R11);
  
  if(Myverbose) {
    cout << endl;
    cout << "AngleX = " << thetaX/units["deg"] << " deg" << endl;
    cout << "AngleY = " << thetaY/units["deg"] << " deg" << endl;
    cout << "AngleZ = " << thetaZ/units["deg"] << " deg" << endl;
    cout << endl;
  }
  
  return  TVector3(thetaX,thetaY,thetaZ);
  
}
//====================================================================
void  GGlobalTools::GetRotationMatrixFromUnitVects(TVector3 UVector,TVector3 VVector,TVector3 WVector,TMatrixD& Rot)
{
  
  //Get rotation matrix from rotated frame unitary orthogonal vectors
  
  Rot.ResizeTo(3,3);
  Rot(0,0) = UVector.X(); Rot(0,1) = UVector.Y(); Rot(0,2) = UVector.Z();
  Rot(1,0) = VVector.X(); Rot(1,1) = VVector.Y(); Rot(1,2) = VVector.Z();
  Rot(2,0) = WVector.X(); Rot(2,1) = WVector.Y(); Rot(2,2) = WVector.Z();
  
}
//====================================================================
void  GGlobalTools::GetMosaicGeoParams(TString GeoConfig,
				       double& R,double w, double& rf, int n_ladders,
				       double& length, double& shift,  double& alpha,
				       TString TheVarPar,bool ShiftFix)
{
  
  // Calculate the parameters of mosaic gemetries
  
  if(GeoConfig == TString("Spiral"))            GetMosaicGeoParams_Spiral(     R,w,rf,n_ladders,length,shift,alpha,TheVarPar,ShiftFix);
  else if(GeoConfig == TString("Alternating"))  GetMosaicGeoParams_Alternating(R,w,rf,n_ladders,length,shift,alpha,TheVarPar,ShiftFix);
  
  return;  
  
}
//====================================================================
void GGlobalTools::GetMosaicGeoParams_Spiral(double& R,double w,double& rf, int n_ladders,
					     double& length, double& shift, double& alpha,
					     TString TheVarPar,bool ShiftFix)
{
  
  // Calculate parameter of spiral mosaic geometry
  
  alpha = 2.0*TMath::Pi()/n_ladders;
  
  if(TheVarPar == TString("Width")) {
    length  = 1.0/(1 - rf);
    length *= (1.0 - TMath::Cos(alpha))/TMath::Sin(alpha);
    length *= w + 2.0*R;
    
    if(ShiftFix) {
      cout << "ERROR: ShiftFix parameter can only be true when the varied parameter is Overlap. Check your inputs. Exiting now!!!" << endl;
      assert(false);
    }
    
  }
  else if(TheVarPar == TString("Overlap")) {
    if(!ShiftFix) {
      rf  = 1.0/length;
      rf *= (1.0 - TMath::Cos(alpha))/TMath::Sin(alpha);
      rf *= w + 2.0*R;
      rf  = 1.0 - rf;
    }
    else {
      rf  = shift*(1.0 - TMath::Cos(alpha))*(w + 2.0*R);
      rf /= R*(1.0 - TMath::Cos(alpha)) - w*TMath::Cos(alpha);
      rf  = 1.0 - rf;
    }
  }
  else if(TheVarPar == TString("Radius")) {
    R  = 0.5*((length*(1.0 - rf)*TMath::Sin(alpha)/(1.0 - TMath::Cos(alpha))) - w);
    
    if(ShiftFix) {
      cout << "ERROR: ShiftFix parameter can only be true when the varied parameter is Overlap. Check your inputs. Exiting now!!!" << endl;
      assert(false);
    }
  }
  else {
    cout << "WARNNING in GetGeoParams_Mosaic_Spiral:: varying parameter " << TheVarPar.Data() << " doesn't exist. Varying the overlap." << endl;
    
    rf  = 1.0/length;
    rf *= (1.0 - TMath::Cos(alpha))/TMath::Sin(alpha);
    rf *= w + 2.0*R;
    rf  = 1.0 - rf;
  }
  
  if(!ShiftFix) {
    shift   = 1.0 - rf;
    shift  *= R*(1.0 - TMath::Cos(alpha)) - w*TMath::Cos(alpha);
    shift  /= 1.0 - TMath::Cos(alpha);
    shift  /= w + 2.0*R;
  }
  
  return;
  
}
//====================================================================
void GGlobalTools::GetMosaicGeoParams_Alternating(double& R,double w,double& rf, int n_ladders,
						  double& length, double& shift, double& alpha,
						  TString TheVarPar,bool ShiftFix)
{
  
  // Calculate parameter of alternating mosaic geometry
  
  alpha = TMath::Pi()/(n_ladders/2);
  
  if(TheVarPar == TString("Width")) {
    length  = 2.0*(w + R)*TMath::Sin(alpha);
    length /= 1.0 - 2.0*rf + TMath::Cos(alpha);
  }
  else if(TheVarPar == TString("Overlap")) {
    rf = (length*(1.0 + TMath::Cos(alpha)) - 2.0*(w + R)*TMath::Sin(alpha))/(2.0*length);
  }
  else if(TheVarPar == TString("Radius")) {
    R = ((length*(1.0 - 2.0*rf + TMath::Cos(alpha)))/(2.0*TMath::Sin(alpha))) - w;
  }
  else {
    cout << "WARNNING in GetGeoParams_Mosaic_Spiral:: varying parameter " << TheVarPar.Data() << " doesn't exist. Varying the overlap." << endl;
    
    rf = (length*(1.0 + TMath::Cos(alpha)) - 2.0*(w + R)*TMath::Sin(alpha))/(2.0*length);
  }
  
  shift   = 1.0 + w/R;
  shift  *= 1.0 + (1.0 - 2.0*rf)*TMath::Cos(alpha);
  shift  /= 1.0 - 2.0*rf + TMath::Cos(alpha);
  
  return;
  
}
//====================================================================
void GGlobalTools::CheckResolutionModelType(TString ModelType)
{
  
  //This function checks if the resolution model type specified in datacard is in predefined list of possible resolution model types
  
  bool IsInList = false;
  for(int i=0;i<int(ResolutionModelTypes.size());i++) {
    if(ModelType == ResolutionModelTypes[i]) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Resolution model type type " << ModelType.Data() << " is unknown!" << endl;
    cout << "Resolution model type should be: " << endl;
    for(int i=0;i<int(ResolutionModelTypes.size());i++) {
      cout << " - " << ResolutionModelTypes[i].Data() << " " << endl;
    }
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
bool  GGlobalTools::IsPointInVoxel(Voxel_t aVoxel,TVector3 Point)
{
  
  //Checks if a given point is inside voxel range
  
  double x      = Point.X();
  double y      = Point.Y();
  double z      = Point.Z();
  double r      = Point.Mag();
  double theta  = Point.Theta();
  double phi    = Point.Phi();
  double phi_v2 = phi;
  if(phi_v2 < 0.0) phi_v2 += 2*TMath::Pi();
    
  if(     x     < aVoxel.Rx[0]     || x     > aVoxel.Rx[1])     return false;
  else if(y     < aVoxel.Ry[0]     || y     > aVoxel.Ry[1])     return false;
  else if(z     < aVoxel.Rz[0]     || z     > aVoxel.Rz[1])     return false;
  else if(r     < aVoxel.Rr[0]     || r     > aVoxel.Rr[1])     return false;
  else if(theta < aVoxel.Rtheta[0] || theta > aVoxel.Rtheta[1]) return false;
  else if(!(phi >= aVoxel.Rphi[0] && phi <= aVoxel.Rphi[1]) && !(phi_v2 >= aVoxel.Rphi[0] && phi_v2 <= aVoxel.Rphi[1])) return false;
  
  return true;
  
}
//====================================================================
bool  GGlobalTools::IsInitCondsInTrackFindingRegion(TrackFindingRegion_t aRegion,
						    TVector3 x0, TVector3 p0)
{
  
  
  bool IsInVoxel = IsPointInVoxel(aRegion.posRange,x0);
  
  bool IsInMomRange = true;
  
  double p = p0.Mag();
  double theta  = p0.Theta();
  double phi    = p0.Phi();
  double phi_v2 = phi;
  if(phi_v2 < 0.0) phi_v2 += 2*TMath::Pi();
  
  if(     p     < aRegion.momRange.Rp[0]     || p     > aRegion.momRange.Rp[1])     IsInMomRange = false;
  else if(theta < aRegion.momRange.Rtheta[0] || theta > aRegion.momRange.Rtheta[1]) IsInMomRange = false;
  else if(!(phi >= aRegion.momRange.Rphi[0] && phi <= aRegion.momRange.Rphi[1]) && !(phi_v2 >= aRegion.momRange.Rphi[0] && phi_v2 <= aRegion.momRange.Rphi[1])) IsInMomRange = false;
  
  return  (IsInVoxel && IsInMomRange);
  
}
//====================================================================
TString  GGlobalTools::GetOutputDirectory(TString TheOutputFile)
{
  
  //Get the output directory from output file generic name
  
  TString Dir("");
  TString Slash("/");
  
  int Last_slash_position = -999;
  for(int i=0;i<TheOutputFile.Length();i++) {
    int idx = TheOutputFile.Length() - i - 1;
    if(TString(TheOutputFile[idx]) == Slash) {
      Last_slash_position = idx;
      break;
    }
  }
  if(Last_slash_position >= 0) {
    Dir = TheOutputFile;
    Dir.Remove(Last_slash_position+1);
  }
  cout << Dir.Data() << endl;
  
  return Dir;
  
}
//====================================================================
bool  GGlobalTools::SetBoolFromString(TString string)
{
  
  if(     string == TString("true")  || string == TString("TRUE")  || string == TString("True"))  return  true;
  else if(string == TString("false") || string == TString("FALSE") || string == TString("False")) return  false;
  else {
    cout << endl;
    cout << "WARNING:: \"" << string.Data() << "\" bool parameter can only have the values true/TRUE/True or false/FALSE/False. Setting it to its defaul value \"false\"." << endl;
    cout << endl;
  }
  
  return  false;
  
}
//====================================================================
TString  GGlobalTools::GetStringFromBool(bool fff) {
  
  if(fff)  return TString("true");
  else     return TString("false");
  
}
//====================================================================
double  GGlobalTools::FromCosThetaToTheta(double costheta)
{
  
  //Get cos(theta) from theta
  
  return  TMath::ACos(costheta);
  
}
//====================================================================
double  GGlobalTools::FromThetaToCosTheta(double theta)
{
  
  //Get theta from cos(theta)
  
  return TMath::Cos(theta);
  
}
//====================================================================
double  GGlobalTools::FromEtaToTheta(double eta)
{
  
  //Get eta from theta
  
  return  2.0*TMath::ATan(TMath::Exp(-eta));
  
}
//====================================================================
double  GGlobalTools::FromThetaToEta(double theta)
{
  
  //Get theta from eta
  
  return  -TMath::Log(TMath::Tan(theta/2.0));
  
}
//====================================================================
double  GGlobalTools::GetEllipseArea(double sigmaU, double sigmaV, double corr)
{
  
  //Get ellipse area
  
  return  TMath::Pi()*sqrt(1.0 - pow(corr,2))*sigmaU*sigmaV;
  
}
//====================================================================
void  GGlobalTools::RemoveElementFromList(int idx, std::vector<int> &List)
{
  
  //Remove element idx from list
  
  if(List.size() == 0) return;
  
  for(int i=0;i<int(List.size());i++) {
    if(idx == List[i]) {
      List.erase(List.begin()+i);
      break;
    }
  }
  
  return;
  
}
//====================================================================
double  GGlobalTools::GetMomentumFromMomVar(double mom, TString Aparticle, TString momVar)
{
  
  double mass   = GetParticleMass(Aparticle);
  int    weight = GetParticleAweight(Aparticle);
  
  if(momVar == TString("p"))             return mom;
  else if(momVar == TString("E")) {
    if(mom < mass) {
      cout << endl;
      cout << "ERROR in GGlobalTools::GetMomentumFromMomVar:: Energy value " << mom/GetUnit("GeV") << " is smaller than particle " << Aparticle.Data() << " mass " << mass/GetUnit("GeV/c2") 
           << ". Check your inputs. Exiting now!!!" 
	   << endl;
      cout << endl;
    }
    return sqrt(pow(mom,2) - pow(mass,2));
  }
  else if(momVar == TString("Ekin"))     return sqrt(pow(mom + mass,2) - pow(mass,2));
  else if(momVar == TString("EkinPerU")) return sqrt(pow(mom*weight + mass,2) - pow(mass,2));
  else                                   return mom;
    
}
//====================================================================
double  GGlobalTools::BetheBlochGeant(double  momentum,
				      TString particle,
				      TString material)
{
 
  double bg = momentum/GetParticleMass(particle);
  
  double density,density_effect_1st,density_effect_2nd,mI,mZA,mZ;
  GetBetheBlochParams(material,
		      density,
		      density_effect_1st,
		      density_effect_2nd,
		      mI,
		      mZA,
		      mZ);
  
  return  BetheBlochGeant(bg,GetParticleCharge(particle),density,density_effect_1st,density_effect_2nd,mI,mZA);
  
}
//====================================================================
double  GGlobalTools::BetheBlochGeant(double bg,
				      int    charge,
				      double density,
				      double density_effect_1st,
				      double density_effect_2nd,
				      double mI,
				      double mZA) {

  // This is the parameterization of the Bethe-Bloch formula inspired by Geant and used by the ALICE collaboration
  //
  // bg                  - beta*gamma
  // charge              - particle's charge
  // density             - material density
  // density_effect_1st  - density effect first  junction point
  // density_effect_2nd  - density effect second junction point
  // mI                  - mean excitation energy [GeV]
  // mZA                 - mean Z/A
  //
  // The returned value is in [dE/distance].
  
  const double mK     = 0.307075e-3*units["GeV"]*units["cm2"]/units["gr"];
  const double me     = 0.511e-3   *units["GeV/c2"];
  const double x0     = density_effect_1st*2.303;   //no units
  const double x1     = density_effect_2nd*2.303;   //no units
  const double bg2    = pow(bg,2);
  const double maxT   = 2*me*bg2;    // neglecting the electron mass
  const double factor = 28.816e-9*units["GeV"]/sqrt(units["gr"]/pow(units["cm"],3));
  
  //*** Density effect
  double       d2   = 0.0; 
  const double x    = TMath::Log(bg);
  const double lhwI = TMath::Log(factor*TMath::Sqrt(density*mZA)/mI);
  if(x > x1) {
    d2 = lhwI + x - 0.5;
  }
  else if(x > x0) {
    const double r = (x1-x)/(x1-x0);
    
    d2 = lhwI + x - 0.5 + (0.5 - lhwI - x0)*pow(r,3);
  }

  double dEOdx  = pow(charge,2)*mK*mZA*(1+bg2)/bg2;
  dEOdx        *= 0.5*TMath::Log(2*me*bg2*maxT/pow(mI,2)) - bg2/(1+bg2) - d2;
  dEOdx        *= density;
  
  return  dEOdx;
  
}
//====================================================================
double  GGlobalTools::BetheBlochPDG(double  momentum,
				    TString particle,
				    TString material) {

  // This is the parameterization of the Bethe-Bloch formula inspired by the PDG
  // The returned value is in [dE/distance].
  
  int    charge     = GetParticleCharge(particle);
  double mass       = GetParticleMass(particle);
  double E          = sqrt(pow(momentum,2) + pow(mass,2));
  double beta       = momentum/E;
  double gamma      = 1.0/sqrt(1.0 - pow(beta,2));
  double beta_gamma = beta*gamma;

  double density,density_effect_1st,density_effect_2nd,mI,mZA,mZ;
  GetBetheBlochParams(material,
		      density,
		      density_effect_1st,
		      density_effect_2nd,
		      mI,
		      mZA,
		      mZ);

  double Factor = 0.5*K*density*mZA*pow(charge,2);
  double LogFactor;
  double F_tau;
  if(particle == TString("e-") || particle == TString("e+")) {
    double Te  = E - mass;
    double tau = Te/mass;

    LogFactor  = pow(tau,2)*(tau + 2.0);
    LogFactor /= (2.0*pow(mI/mass,2));
    LogFactor  = TMath::Log(LogFactor);
    LogFactor /= pow(beta,2);
     
    if(particle == TString("e+")) F_tau  = 2.0*TMath::Log(2.0) - (1.0/12.0)*pow(beta,2)*(23.0 + (14.0/(tau + 2.0)) + (10./pow(tau + 2.0,2)) + (4.0/pow(tau + 2.0,3)));
    else                          F_tau  = 1 - pow(beta,2) + ((1.0/8.0)*pow(tau,2) - (2.0*tau+1.0)*TMath::Log(2))/pow(tau + 1.0,2);
    F_tau     /= pow(beta,2);
  }
  else {
    double mass_e = GetParticleMass(TString("e-"));
    double mass_r = mass_e/mass;
    double Wmax   = 2.0*mass_e*pow(beta_gamma,2)/(1.0 + 2.0*gamma*mass_r + pow(mass_r,2));

    LogFactor  = (2.0*mass_e*pow(beta_gamma,2)*Wmax)/pow(mI,2);
    LogFactor  = TMath::Log(LogFactor);
    LogFactor /= pow(beta,2);

    F_tau      = -2.0;
  }
  
  //Density effect corrections
  //Currently for silicon
  double C_delta = -4.4351;
  double m_delta =  3.2546;
  double a_delta =  0.14921;
  double b_delta =  4.6052;

  double delta = 0.0;
  double X = TMath::Log10(beta_gamma);
  if(X >= density_effect_1st && X < density_effect_2nd) delta = b_delta*X + C_delta + a_delta*TMath::Power(density_effect_2nd - X,m_delta);
  else if(X >= density_effect_2nd)                      delta = b_delta*X + C_delta;
  delta /= pow(beta,2);

  //Shell corrections
  double Shell   = 0.0;
  double Shell_1 = 0.0;
  double Shell_2 = 0.0;

  Shell_1 +=  0.42237700*TMath::Power(beta_gamma,-2);
  Shell_1 +=  0.03040430*TMath::Power(beta_gamma,-4);
  Shell_1 += -0.00038106*TMath::Power(beta_gamma,-6);
  Shell_1 *= 1.0e-6*TMath::Power(mI/units["eV"],2);
  Shell   += Shell_1;

  Shell_2 +=  3.85019000*TMath::Power(beta_gamma,-2);
  Shell_2 +=  0.16679890*TMath::Power(beta_gamma,-4);
  Shell_2 += -0.00157955*TMath::Power(beta_gamma,-6);
  Shell_2 *= 1.0e-9*TMath::Power(mI/units["eV"],3);
  Shell   += Shell_2;
  
  Shell   *= 2.0/mZ;
  Shell   /= pow(beta,2);
  
  double dEOdx  = LogFactor + F_tau - delta - Shell;
  dEOdx        *= Factor;
  
  return  dEOdx;
  
}
//====================================================================
void  GGlobalTools::GetBetheBlochParams(TString material,
					double& density,
					double& density_effect_1st,
					double& density_effect_2nd,
					double& mI,
					double& mZA,
					double& mZ)
{

  bool IsInList = false;
  for(map<TString,Material_t>::iterator i=MaterialMap.begin(); i!=MaterialMap.end(); ++i) {
    if(material == (*i).first) {
      IsInList = true;
      break;
    }
  }
  
  if(!IsInList) {
    cout << "Material type " << material.Data() << " is unknown!" << endl;
    cout << "Material should be: " << endl;
    for(map<TString,Material_t>::iterator ii=MaterialMap.begin(); ii!=MaterialMap.end(); ++ii) {
      cout << " - " << ((*ii).first).Data() << " " << endl;
    }
    cout << endl;

    assert(false);
  }
  
  density            = MaterialMap[material].density;
  density_effect_1st = MaterialMap[material].density_effect_1st;
  density_effect_2nd = MaterialMap[material].density_effect_2nd;
  mI                 = MaterialMap[material].mI;
  mZA                = MaterialMap[material].mZA;
  mZ                 = MaterialMap[material].mZ;
  
  return;
  
}
//====================================================================
void  GGlobalTools::SetKappa(double akappa)
{
  
  if(akappa < 1.0e-6*units["GeV"]) {
    cout << endl;
    cout << "ERROR in function GGlobalTools::SetKappa:: kappa should only have positive values" << endl;
    cout << endl;
    
    assert(false);
  }
  
  kappa = akappa;
  
  return;
  
}
//====================================================================
double  GGlobalTools::GetSigmaELoss(double  momentum,
				    TString particle,
				    TString material,
				    double dE)
{

  //Approximate energy loss fluctuation (M.Ivanov)
  // - kappa to be tuned.
  double Eloss   = TMath::Abs(dE/units["GeV"]);
  double sigma   = kappa*sqrt(Eloss);
  
  return  sigma;
  
}
//====================================================================
double  GGlobalTools::GetSigmaELossPDG(double  momentum,
			               TString particle,
			               TString material,
				       double X)
{
  
  //Use equation for the Eloss fluctuations width at half maximum
  // PDG note 27: Passage of particles through matter, equation in text 
  // just after equation 27.10, w = 4ξ
  
  double beta    = momentum/sqrt(pow(momentum,2) + pow(GetParticleMass(particle),2));
  double density = MaterialMap[material].density;
  double mZA     = MaterialMap[material].mZA;
  
  double sigma = 2.0*K*mZA*X*density/pow(beta,2);
  
  return  sigma;
  
}
//====================================================================

