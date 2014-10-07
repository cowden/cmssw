#include "DQM/HcalMonitorClient/interface/ZDCMonitorClient.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

//
// -------------------------------------- Constructor --------------------------------------------
//
ZDCMonitorClient::ZDCMonitorClient(const edm::ParameterSet& ps)
{
  edm::LogInfo("ZDCMonitorClient") <<  "Constructor  ZDCMonitorClient::ZDCMonitorClient " << std::endl;

}

//
// -- Destructor
//
ZDCMonitorClient::~ZDCMonitorClient()
{
  edm::LogInfo("ZDCMonitorClient") <<  "Destructor ZDCMonitorClient::~ZDCMonitorClient " << std::endl;
}

//
// -------------------------------------- beginJob --------------------------------------------
//
void ZDCMonitorClient::beginJob()
{
  edm::LogInfo("ZDCMonitorClient") <<  "ZDCMonitorClient::beginJob " << std::endl;
}
//
// -------------------------------------- get and book in the endJob --------------------------------------------
//
void ZDCMonitorClient::dqmEndJob(DQMStore::IBooker& ibooker_, DQMStore::IGetter& igetter_)
{
  // create and cd into new folder
  ibooker_.setCurrentFolder("Hcal/ZDCMonitor");

  //get available histograms
  MonitorElement* numerator = igetter_.get("evNumerator");
  MonitorElement* denominator = igetter_.get("evDenomenator");

  if (!numerator || !denominator)
    {
      edm::LogError("ZDCMonitorClient") <<  "MEs not found!" << std::endl;
      return;
    }


  //book new histogram
  h_eventBased = ibooker_.book1D("ptRatio","pt ratio pf matched objects",50,0.,100.);
  h_eventBased->setAxisTitle("pt [GeV]");

  for (int iBin=1; iBin<numerator->getNbinsX(); ++iBin)
    {
      if(denominator->getBinContent(iBin) == 0)
	h_eventBased->setBinContent(iBin, 0.);
      else
	h_eventBased->setBinContent(iBin, numerator->getBinContent(iBin) / denominator->getBinContent(iBin));
    }
}

//
// -------------------------------------- get in the endLumi if needed --------------------------------------------
//
void ZDCMonitorClient::dqmEndLuminosityBlock(DQMStore::IBooker & ibooker_, DQMStore::IGetter & igetter_, edm::LuminosityBlock const & iLumi, edm::EventSetup const& iSetup) 
{
  edm::LogInfo("ZDCMonitorClient") <<  "ZDCMonitorClient::endLumi " << std::endl;

}
