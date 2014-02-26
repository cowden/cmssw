#ifndef SimG4CMS_CFCShowerLibrary_h
#define SimG4CMS_CFCShowerLibrary_h 1
///////////////////////////////////////////////////////////////////////////////
// File: CFCShowerLibrary.h
// Description: Gets information from a shower library
///////////////////////////////////////////////////////////////////////////////

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimG4CMS/Calo/interface/HFFibre.h"

#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
 
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <memory>

struct photon {
  float wavelength;
//  float energy;
//  float na;
  float x;
  float y;
  float z;
  float t;
  int type; // 0: cherenkov   1: scintilation
//  float tprop;
  
  inline photon(double w, double gx, double gy,double gz, double gt, int tp)
    :wavelength(w),x(gx),y(gy),z(gz),t(gt),type(tp)
    { }
};

class G4Step;

class CFCShowerLibrary {
  
public:
  
  //Constructor and Destructor
  CFCShowerLibrary(edm::ParameterSet const & p, std::vector<double>& gpar);
  ~CFCShowerLibrary();

public:

  struct Hit {
    Hit() {}
    G4ThreeVector     position;  // local coordinate
    int               type;
    double            lambda;
    double            time;
  };

  void                initRun(G4ParticleTable * theParticleTable);
  std::vector<Hit>    getHits(G4Step * aStep, bool &ok);

protected:

  void                getRecord(int, std::string, int);
  void                loadEventInfo(TTree *);
  void                interpolate(int, double);
  void                extrapolate(int, double);
  void                storePhoton(int j, int type);

private:

  std::vector<double> gpar;
  TFile*              hfile;
  std::vector<TTree*> egamma_subTrees, hadron_subTrees;

  int                 emPDG, epPDG, gammaPDG, mumPDG, mupPDG;
  int                 pi0PDG, etaPDG, nuePDG, numuPDG, nutauPDG;
  int                 anuePDG, anumuPDG, anutauPDG, geantinoPDG;
/*
  int                 nMomBin, totEvents, evtPerBin;
  std::vector<double> pmom;
  std::vector<std::string> dirs;
*/
  int egamma_nMomBin, hadron_nMomBin;
  std::vector<int>    *egamma_evtPerBinVec, *hadron_evtPerBinVec;
  std::vector<double> *egamma_pmomVec, *hadron_pmomVec;
  std::vector<std::string> *egamma_dirsVec, *hadron_dirsVec;

  int                 scin_npe, opt_npe, npe;
  std::vector<photon> scin_pe, opt_pe, pe;
  
  std::vector<float> *scin_fx, *scin_fy, *scin_fz, *scin_t, *scin_wavelength;
  std::vector<float> *opt_fx, *opt_fy, *opt_fz, *opt_t, *opt_wavelength;

};
#endif
