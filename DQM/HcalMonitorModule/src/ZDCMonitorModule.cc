#include "DQM/HcalMonitorModule/interface/ZDCMonitorModule.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


// Geometry
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"


#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

//
// -------------------------------------- Constructor --------------------------------------------
//
ZDCMonitorModule::ZDCMonitorModule(const edm::ParameterSet& ps)
{
  edm::LogInfo("ZDCMonitorModule") <<  "Constructor  ZDCMonitorModule::ZDCMonitorModule " << std::endl;
  
 
}

//
// -- Destructor
//
ZDCMonitorModule::~ZDCMonitorModule()
{
  edm::LogInfo("ZDCMonitorModule") <<  "Destructor ZDCMonitorModule::~ZDCMonitorModule " << std::endl;
}

//
// -------------------------------------- beginRun --------------------------------------------
//
void ZDCMonitorModule::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("ZDCMonitorModule") <<  "ZDCMonitorModule::beginRun" << std::endl;
}
//
// -------------------------------------- bookHistos --------------------------------------------
//
void ZDCMonitorModule::bookHistograms(DQMStore::IBooker & ibooker_, edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("ZDCMonitorModule") <<  "ZDCMonitorModule::bookHistograms" << std::endl;

  ibooker_.setCurrentFolder("Hcal/ZDCMonitor");
  h_evNum = ibooker_.book1D("evNumerator","Numerator",50,0.,100.);
  h_evDen = ibooker_.book1D("evDenomenator","Denomenator",50,0.,100.);
  
}
//
// -------------------------------------- beginLuminosityBlock --------------------------------------------
//
void ZDCMonitorModule::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
  edm::LogInfo("ZDCMonitorModule") <<  "ZDCMonitorModule::beginLuminosityBlock" << std::endl;
}


//
// -------------------------------------- Analyze --------------------------------------------
//
void ZDCMonitorModule::analyze(edm::Event const& e, edm::EventSetup const& eSetup)
{
  edm::LogInfo("ZDCMonitorModule") <<  "ZDCMonitorModule::analyze" << std::endl;


}
//
// -------------------------------------- endLuminosityBlock --------------------------------------------
//
void ZDCMonitorModule::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
  edm::LogInfo("ZDCMonitorModule") <<  "ZDCMonitorModule::endLuminosityBlock" << std::endl;
}


//
// -------------------------------------- endRun --------------------------------------------
//
void ZDCMonitorModule::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
  edm::LogInfo("ZDCMonitorModule") <<  "ZDCMonitorModule::endRun" << std::endl;
}



