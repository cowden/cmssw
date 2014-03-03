///////////////////////////////////////////////////////////////////////////////
// File: CFCShowerLibrary.cc
// Description: Shower library for Combined Forward Calorimeter
///////////////////////////////////////////////////////////////////////////////

#include "SimG4CMS/Calo/interface/CFCShowerLibrary.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "G4VPhysicalVolume.hh"
#include "G4NavigationHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

//#define DebugLog

CFCShowerLibrary::CFCShowerLibrary(edm::ParameterSet const & p,
				   std::vector<double> & gp) : gpar(gp) {
  

  edm::ParameterSet m_CFC = p.getParameter<edm::ParameterSet>("CFCSD");
  edm::FileInPath fp      = m_CFC.getParameter<edm::FileInPath>("FileName");
  std::string pTreeName   = fp.fullPath();
  if (pTreeName.find(".") == 0) pTreeName.erase(0,2);
  const char* nTree       = pTreeName.c_str();
  hfile                   = TFile::Open(nTree);

  if (!hfile->IsOpen()) { 
    edm::LogError("CFCShower") << "CFCShowerLibrary: opening " << nTree 
			       << " failed";
    throw cms::Exception("Unknown", "CFCShowerLibrary") 
      << "Opening of " << pTreeName << " fails\n";
  } else {
    edm::LogInfo("CFCShower") << "CFCShowerLibrary: opening " << nTree 
			      << " successfully"; 
  }

  TTree * infoTree = (TTree*) hfile->Get("info/infoTree");
  if( infoTree ){
     loadEventInfo(infoTree);
  }else{
    edm::LogError("CFCShower") << "CFCShowerLibrary: info Tree does not "
                              << "exist";
    throw cms::Exception("Unknown", "CFCShowerLibrary")
      << "Info tree absent\n";
  }

  for(unsigned int id=0; id<egamma_dirsVec->size(); id++){
     std::string treeName = "egamma/" + egamma_dirsVec->at(id) + "/event";
     TTree * subTree = (TTree*) hfile->Get(treeName.c_str());
     egamma_subTrees.push_back(subTree);
  }
  for(unsigned int id=0; id<hadron_dirsVec->size(); id++){
     std::string treeName = "hadron/" + hadron_dirsVec->at(id) + "/event";
     TTree * subTree = (TTree*) hfile->Get(treeName.c_str());
     hadron_subTrees.push_back(subTree);
  }

  scin_fx = new std::vector<float>(); scin_fy = new std::vector<float>(); scin_fz = new std::vector<float>(); scin_t = new std::vector<float>(); scin_wavelength = new std::vector<float>();
  opt_fx = new std::vector<float>(); opt_fy = new std::vector<float>(); opt_fz = new std::vector<float>(); opt_t = new std::vector<float>(); opt_wavelength = new std::vector<float>();

  emPDG = epPDG = gammaPDG = 0;
  pi0PDG = etaPDG = nuePDG = numuPDG = nutauPDG= 0;
  anuePDG= anumuPDG = anutauPDG = geantinoPDG = 0;
}

CFCShowerLibrary::~CFCShowerLibrary() {
  if (hfile)  hfile->Close();
}

void CFCShowerLibrary::initRun(G4ParticleTable * theParticleTable) {

  emPDG      = theParticleTable->FindParticle("e-")->GetPDGEncoding();
  epPDG      = theParticleTable->FindParticle("e+")->GetPDGEncoding();
  gammaPDG   = theParticleTable->FindParticle("gamma")->GetPDGEncoding();
  mumPDG     = theParticleTable->FindParticle("mu-")->GetPDGEncoding();
  mupPDG     = theParticleTable->FindParticle("mu+")->GetPDGEncoding();
  pi0PDG     = theParticleTable->FindParticle("pi0")->GetPDGEncoding();
  etaPDG     = theParticleTable->FindParticle("eta")->GetPDGEncoding();
  nuePDG     = theParticleTable->FindParticle("nu_e")->GetPDGEncoding();
  numuPDG    = theParticleTable->FindParticle("nu_mu")->GetPDGEncoding();
  nutauPDG   = theParticleTable->FindParticle("nu_tau")->GetPDGEncoding();
  anuePDG    = theParticleTable->FindParticle("anti_nu_e")->GetPDGEncoding();
  anumuPDG   = theParticleTable->FindParticle("anti_nu_mu")->GetPDGEncoding();
  anutauPDG  = theParticleTable->FindParticle("anti_nu_tau")->GetPDGEncoding();
  geantinoPDG= theParticleTable->FindParticle("geantino")->GetPDGEncoding();
#ifdef DebugLog
  edm::LogInfo("CFCShower") << "CFCShowerLibrary: Particle codes for e- = " 
			    << emPDG << ", e+ = " << epPDG << ", gamma = " 
			    << gammaPDG << ", pi0 = " << pi0PDG << ", eta = " 
			    << etaPDG << ", geantino = " << geantinoPDG 
			    << "\n        nu_e = " << nuePDG << ", nu_mu = " 
			    << numuPDG << ", nu_tau = " << nutauPDG 
			    << ", anti_nu_e = " << anuePDG << ", anti_nu_mu = " 
			    << anumuPDG << ", anti_nu_tau = " << anutauPDG;
#endif
}


std::vector<CFCShowerLibrary::Hit> CFCShowerLibrary::getHits(G4Step * aStep,
							     bool & ok) {

  G4StepPoint * preStepPoint  = aStep->GetPreStepPoint();
  G4StepPoint * postStepPoint = aStep->GetPostStepPoint();
  G4Track *     track    = aStep->GetTrack();
  // Get Z-direction
  const G4DynamicParticle *aParticle = track->GetDynamicParticle();
  G4ThreeVector momDir = aParticle->GetMomentumDirection();
//  double mom = aParticle->GetTotalMomentum();

  G4ThreeVector hitPoint = preStepPoint->GetPosition();
  G4String      partType = track->GetDefinition()->GetParticleName();
  int           parCode  = track->GetDefinition()->GetPDGEncoding();

  std::vector<CFCShowerLibrary::Hit> hit;
  ok = false;
  if (parCode == pi0PDG || parCode == etaPDG || parCode == nuePDG ||
      parCode == numuPDG || parCode == nutauPDG || parCode == anuePDG ||
      parCode == anumuPDG || parCode == anutauPDG || parCode == geantinoPDG)
    return hit;
  ok = true;

  double tSlice = (postStepPoint->GetGlobalTime())/nanosecond;
//  double pin    = preStepPoint->GetTotalEnergy();
  double pin    = preStepPoint->GetTotalEnergy()/GeV;
//  double pz     = momDir.z();
//  double zint   = hitPoint.z();

//  std::cout<<"partType : "<<partType<<"  hitPoint : "<<hitPoint<<"  mom : "<<mom/GeV<<"  pin : "<<pin<<std::endl;

  double sphi   = sin(momDir.phi());
  double cphi   = cos(momDir.phi());
  double ctheta = cos(momDir.theta());
  double stheta = sin(momDir.theta());

  if (parCode == emPDG || parCode == epPDG || parCode == gammaPDG ) {
    if (pin<egamma_pmomVec->at(egamma_nMomBin-1)) {
//      interpolate(0, pin); 
      interpolate2(0, pin); 
    } else {
//      extrapolate(0, pin); 
      extrapolate2(0, pin); 
    }
  } else { 
    if (pin<hadron_pmomVec->at(hadron_nMomBin-1)) {
//      interpolate(1, pin);
      interpolate2(1, pin);
    } else {
//      extrapolate(1, pin);
      extrapolate2(1, pin);
    } 
  }

  int nHit = 0;
  CFCShowerLibrary::Hit oneHit;
  for (int i = 0; i < npe; i++) {
    double zv = std::abs(pe[i].z); // abs local z
    // Updated coordinate transformation from local
    //  back to global using two Euler angles: phi and theta
    double pex = pe[i].x;
    double pey = pe[i].y;

    double xx = pex*ctheta*cphi - pey*sphi + zv*stheta*cphi;
    double yy = pex*ctheta*sphi + pey*cphi + zv*stheta*sphi;
    double zz = -pex*stheta + zv*ctheta;

    G4ThreeVector pos  = hitPoint + G4ThreeVector(xx,yy,zz);

    oneHit.position = pos;
//    oneHit.depth    = depth;
    oneHit.time     = tSlice+pe[i].t;
    oneHit.type     = pe[i].type;
    oneHit.weight   = pe[i].weight;
    hit.push_back(oneHit);
    nHit++;
  }

  return hit;
}


void CFCShowerLibrary::loadEventInfo(TTree* tree) {

  egamma_evtPerBinVec = new std::vector<int>(); hadron_evtPerBinVec = new std::vector<int>();
  egamma_pmomVec = new std::vector<double>(); hadron_pmomVec = new std::vector<double>();
  egamma_dirsVec = new std::vector<std::string>(); hadron_dirsVec = new std::vector<std::string>();

  tree->SetBranchAddress("egamma_nMomBin", &egamma_nMomBin);
  tree->SetBranchAddress("egamma_evtPerBinVec", &egamma_evtPerBinVec);
  tree->SetBranchAddress("egamma_pmomVec", &egamma_pmomVec);
  tree->SetBranchAddress("egamma_dirsVec", &egamma_dirsVec);
  tree->SetBranchAddress("hadron_nMomBin", &hadron_nMomBin);
  tree->SetBranchAddress("hadron_evtPerBinVec", &hadron_evtPerBinVec);
  tree->SetBranchAddress("hadron_pmomVec", &hadron_pmomVec);
  tree->SetBranchAddress("hadron_dirsVec", &hadron_dirsVec);
  tree->GetEntry(0);
//  for (unsigned int i=0; i<pmom.size(); i++) pmom[i] *= GeV;
}


void CFCShowerLibrary::interpolate2(int type, double pin) {

  int irc[2];
  double w = 0.;
  double r = G4UniformRand();

  int binIdx[2];

  std::vector<double> pmom; int nMomBin;
  std::vector<int> evtPerBin;
  std::vector<std::string> dirs;
  if( type == 0 ){
    pmom = (*egamma_pmomVec); nMomBin = egamma_nMomBin; evtPerBin = (*egamma_evtPerBinVec); dirs = (*egamma_dirsVec);
  }else{
    pmom = (*hadron_pmomVec); nMomBin = hadron_nMomBin; evtPerBin = (*hadron_evtPerBinVec); dirs = (*hadron_dirsVec);
  }

  irc[0] = -1; irc[1] = -1;
  binIdx[0] = -1; binIdx[1] = -1;

  if (pin<pmom[0]) {
    w = pin/pmom[0];
    irc[1] = int(evtPerBin[0]*r);
    irc[0] = 0;
    binIdx[0] = -1; binIdx[1] = 0;
  } else {
    for (int j=0; j<nMomBin-1; j++) {
      if (pin >= pmom[j] && pin < pmom[j+1]) {
        w = (pin-pmom[j])/(pmom[j+1]-pmom[j]);
        if (j == nMomBin-2) {
          irc[1] = int(evtPerBin[j+1]*0.5*r);
        } else {
          irc[1] = int(evtPerBin[j+1]*r);
        }
        r = G4UniformRand();
        irc[0] = int(evtPerBin[j]*r);
        binIdx[0] = j; binIdx[1] = j+1;
      }
    }
  }
  pe.clear();
  npe       = 0;

  double newWeight_opt = 0, newWeight_scin =0;
  if( binIdx[0] == -1 ){
     std::string dirName = dirs[binIdx[1]];
     getRecord (type, dirName, irc[1]);

     newWeight_opt = w*(float)opt_num/opt_t->size();
     newWeight_scin = w*(float)scin_num/scin_t->size();

     int noptPhoton = opt_t->size();
     int nscinPhoton = scin_t->size();
     for (int j=0; j<noptPhoton; j++){
        storePhoton (j, 0, newWeight_opt);
     }
     for (int j=0; j<nscinPhoton; j++){
        storePhoton (j, 1, newWeight_scin);
     }
  }else if( w < 0.5 ){
     std::string dirName = dirs[binIdx[1]];
     getRecord (type, dirName, irc[1]);
     newWeight_opt = w*(float)opt_num;
     newWeight_scin = w*(float)scin_num;

     dirName = dirs[binIdx[0]];
     getRecord (type, dirName, irc[0]);
     newWeight_opt += (1-w)*(float)opt_num;
     newWeight_scin += (1-w)*(float)scin_num;

     newWeight_opt /= opt_t->size();
     newWeight_scin /= scin_t->size();

     int noptPhoton = opt_t->size();
     int nscinPhoton = scin_t->size();
     for (int j=0; j<noptPhoton; j++){
        storePhoton (j, 0, newWeight_opt);
     }
     for (int j=0; j<nscinPhoton; j++){
        storePhoton (j, 1, newWeight_scin);
     }
  }else{
     std::string dirName = dirs[binIdx[0]];
     getRecord (type, dirName, irc[0]);
     newWeight_opt = (1-w)*(float)opt_num;
     newWeight_scin = (1-w)*(float)scin_num;

     dirName = dirs[binIdx[1]];
     getRecord (type, dirName, irc[1]);
     newWeight_opt += w*(float)opt_num;
     newWeight_scin += w*(float)scin_num;

     newWeight_opt /= opt_t->size();
     newWeight_scin /= scin_t->size();

     int noptPhoton = opt_t->size();
     int nscinPhoton = scin_t->size();
     for (int j=0; j<noptPhoton; j++){
        storePhoton (j, 0, newWeight_opt);
     }
     for (int j=0; j<nscinPhoton; j++){
        storePhoton (j, 1, newWeight_scin);
     }
  }
}


void CFCShowerLibrary::interpolate(int type, double pin) {

  int irc[2];
  double w = 0.;
  double r = G4UniformRand();

  int binIdx[2];

  std::vector<double> pmom; int nMomBin;
  std::vector<int> evtPerBin;
  std::vector<std::string> dirs;
  if( type == 0 ){
    pmom = (*egamma_pmomVec); nMomBin = egamma_nMomBin; evtPerBin = (*egamma_evtPerBinVec); dirs = (*egamma_dirsVec);
  }else{
    pmom = (*hadron_pmomVec); nMomBin = hadron_nMomBin; evtPerBin = (*hadron_evtPerBinVec); dirs = (*hadron_dirsVec);
  }

  if (pin<pmom[0]) {
    w = pin/pmom[0];
    irc[1] = int(evtPerBin[0]*r);
    irc[0] = 0;
    binIdx[0] = -1; binIdx[1] = 0;
  } else {
    for (int j=0; j<nMomBin-1; j++) {
      if (pin >= pmom[j] && pin < pmom[j+1]) {
        w = (pin-pmom[j])/(pmom[j+1]-pmom[j]);
        if (j == nMomBin-2) {
          irc[1] = int(evtPerBin[j+1]*0.5*r);
        } else {
          irc[1] = int(evtPerBin[j+1]*r);
        }
        r = G4UniformRand();
        irc[0] = int(evtPerBin[j]*r);
        binIdx[0] = j; binIdx[1] = j+1;
      }
    }
  }
  pe.clear();
  npe       = 0;
  for (int ir=0; ir < 2; ir++) {
    if (irc[ir]>0) {
      std::string dirName = dirs[binIdx[ir]];
      getRecord (type, dirName, irc[ir]);
      int noptPhoton = opt_t->size();
      int nscinPhoton = scin_t->size();
      for (int j=0; j<noptPhoton; j++){
        r = G4UniformRand();
        if ((ir==0 && r > w) || (ir > 0 && r < w)) {
          storePhoton (j, 0);
        }
      }
      for (int j=0; j<nscinPhoton; j++){
        r = G4UniformRand();
        if ((ir==0 && r > w) || (ir > 0 && r < w)) {
          storePhoton (j, 1);
        }
      }
    }
  }
}


void CFCShowerLibrary::extrapolate2(int type, double pin) {

  std::vector<double> pmom; int nMomBin = -1;
  std::vector<int> evtPerBin;
  std::vector<std::string> dirs;
  if( nMomBin == -1 ){}
  if( type == 0 ){
    pmom = (*egamma_pmomVec); nMomBin = egamma_nMomBin; evtPerBin = (*egamma_evtPerBinVec); dirs = (*egamma_dirsVec);
  }else{
    pmom = (*hadron_pmomVec); nMomBin = hadron_nMomBin; evtPerBin = (*hadron_evtPerBinVec); dirs = (*hadron_dirsVec);
  }

  double r = G4UniformRand();
  double w = pin/pmom.back();
  int pickedEvtIdx = int(evtPerBin.back()*r);

  pe.clear();
  npe       = 0;

  std::string dirName = dirs.back();
  getRecord (type, dirName, pickedEvtIdx);

  double newWeight_opt = w*(float)opt_num/opt_t->size();
  double newWeight_scin = w*(float)scin_num/scin_t->size();

  int noptPhoton = opt_t->size();
  int nscinPhoton = scin_t->size();
  for (int j=0; j<noptPhoton; j++){
     storePhoton (j, 0, newWeight_opt);
  }
  for (int j=0; j<nscinPhoton; j++){
     storePhoton (j, 1, newWeight_scin);
  }

}


void CFCShowerLibrary::extrapolate(int type, double pin) {

  std::vector<double> pmom; int nMomBin;
  std::vector<int> evtPerBin;
  std::vector<std::string> dirs;
  if( type == 0 ){
    pmom = (*egamma_pmomVec); nMomBin = egamma_nMomBin; evtPerBin = (*egamma_evtPerBinVec); dirs = (*egamma_dirsVec);
  }else{
    pmom = (*hadron_pmomVec); nMomBin = hadron_nMomBin; evtPerBin = (*hadron_evtPerBinVec); dirs = (*hadron_dirsVec);
  }

  int nrec   = int(pin/pmom[nMomBin-1]);
  double w   = (pin - pmom[nMomBin-1]*nrec)/pmom[nMomBin-1];
  nrec++;

  std::vector<int> irc(nrec);
  std::vector<int> binIdx(nrec);

  for (int ir=0; ir<nrec; ir++) {
    double r = G4UniformRand();
//    irc[ir] = int(evtPerBin.back()*0.5*r);
    irc[ir] = int(evtPerBin.back()*r);
    binIdx[ir] = nMomBin-1;
    if (irc[ir]<1) {
      edm::LogWarning("CFCShower") << "CFCShowerLibrary:: Illegal irc[" << ir
                                  << "] = " << irc[ir] << " now set to 1";
      irc[ir] = 1;
    } else if (irc[ir] > evtPerBin.back()) {
      edm::LogWarning("CFCShower") << "CFCShowerLibrary:: Illegal irc[" << ir
                                  << "] = " << irc[ir] << " now set to "
                                  << evtPerBin.back();
      irc[ir] = evtPerBin.back();
    }
  }

  pe.clear();
  npe       = 0;
  for (int ir=0; ir<nrec; ir++) {
    if (irc[ir]>=0) {
      std::string dirName = dirs[binIdx[ir]];
      getRecord (type, dirName, irc[ir]);
      int noptPhoton = opt_t->size();
      int nscinPhoton = scin_t->size();
      for (int j=0; j<noptPhoton; j++){
        double r = G4UniformRand();
        if( ir != nrec-1 || r < w ){
          storePhoton (j, 0);
        }
      }
      for (int j=0; j<nscinPhoton; j++){
        double r = G4UniformRand();
        if( ir != nrec-1 || r < w ){
          storePhoton (j, 1);
        }
      }
    }
  }
}


void CFCShowerLibrary::storePhoton(int j, int type){

   float wavelength, x, y, z, t, wt;
   if( type == 0 ){
      wavelength = opt_wavelength->at(j); x = opt_fx->at(j); y = opt_fy->at(j); z = opt_fz->at(j); t = opt_t->at(j); wt = (float)opt_num/opt_t->size();
   }else if(type ==1 ){
      wavelength = scin_wavelength->at(j); x = scin_fx->at(j); y = scin_fy->at(j); z = scin_fz->at(j); t = scin_t->at(j); wt = (float)scin_num/scin_t->size();
   }else{
      wavelength = -1; x = -1; y = -1; z = -1; t = -1; wt = -1;
   }
   photon ph(wavelength, x, y, z, t, type, wt);

   pe.push_back(ph);

   npe++;
   
}


void CFCShowerLibrary::storePhoton(int j, int type, double exWeight){

   float wavelength, x, y, z, t;
   if( type == 0 ){
      wavelength = opt_wavelength->at(j); x = opt_fx->at(j); y = opt_fy->at(j); z = opt_fz->at(j); t = opt_t->at(j);
   }else if(type ==1 ){
      wavelength = scin_wavelength->at(j); x = scin_fx->at(j); y = scin_fy->at(j); z = scin_fz->at(j); t = scin_t->at(j);
   }else{
      wavelength = -1; x = -1; y = -1; z = -1; t = -1;
   }
   photon ph(wavelength, x, y, z, t, type, exWeight);

   pe.push_back(ph);

   npe++;
   
}


void CFCShowerLibrary::getRecord(int type, std::string dirName, int record){

   std::vector<std::string> dirs;
   if( type == 0 ){
      dirs = (*egamma_dirsVec);
   }else{
      dirs = (*hadron_dirsVec);
   }
   for(unsigned int id=0; id<dirs.size(); id++){
      if(dirs[id] == dirName && type == 0 ){
         egamma_subTrees[id]->SetBranchAddress("scin_fx", &scin_fx);
         egamma_subTrees[id]->SetBranchAddress("scin_fy", &scin_fy);
         egamma_subTrees[id]->SetBranchAddress("scin_fz", &scin_fz);
         egamma_subTrees[id]->SetBranchAddress("scin_t", &scin_t);
         egamma_subTrees[id]->SetBranchAddress("scin_wavelength", &scin_wavelength);
         egamma_subTrees[id]->SetBranchAddress("scin_num", &scin_num);

         egamma_subTrees[id]->SetBranchAddress("opt_fx", &opt_fx);
         egamma_subTrees[id]->SetBranchAddress("opt_fy", &opt_fy);
         egamma_subTrees[id]->SetBranchAddress("opt_fz", &opt_fz);
         egamma_subTrees[id]->SetBranchAddress("opt_t", &opt_t);
         egamma_subTrees[id]->SetBranchAddress("opt_wavelength", &opt_wavelength);
         egamma_subTrees[id]->SetBranchAddress("opt_num", &opt_num);
         
         egamma_subTrees[id]->GetEntry(record);
      }
      if(dirs[id] == dirName && type == 1 ){
         hadron_subTrees[id]->SetBranchAddress("scin_fx", &scin_fx);
         hadron_subTrees[id]->SetBranchAddress("scin_fy", &scin_fy);
         hadron_subTrees[id]->SetBranchAddress("scin_fz", &scin_fz);
         hadron_subTrees[id]->SetBranchAddress("scin_t", &scin_t);
         hadron_subTrees[id]->SetBranchAddress("scin_wavelength", &scin_wavelength);
         hadron_subTrees[id]->SetBranchAddress("scin_num", &scin_num);

         hadron_subTrees[id]->SetBranchAddress("opt_fx", &opt_fx);
         hadron_subTrees[id]->SetBranchAddress("opt_fy", &opt_fy);
         hadron_subTrees[id]->SetBranchAddress("opt_fz", &opt_fz);
         hadron_subTrees[id]->SetBranchAddress("opt_t", &opt_t);
         hadron_subTrees[id]->SetBranchAddress("opt_wavelength", &opt_wavelength);
         hadron_subTrees[id]->SetBranchAddress("opt_num", &opt_num);
         
         hadron_subTrees[id]->GetEntry(record);
      }
   }
}
