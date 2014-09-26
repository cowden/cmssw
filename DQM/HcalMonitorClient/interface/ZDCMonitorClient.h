#ifndef ZDCMonitorClient_H
#define ZDCMonitorClient_H


#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/CPUTimer.h"

class TH2F;
class TH1F;
class TFile;

class ZDCMonitorClient : public DQMEDHarvester {
  
public:
  
  /// Constructor
  ZDCMonitorClient(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~ZDCMonitorClient();

protected:

  void beginJob();
  void dqmEndLuminosityBlock(DQMStore::IBooker &, DQMStore::IGetter &, edm::LuminosityBlock const &, edm::EventSetup const &);
  void dqmEndJob(DQMStore::IBooker &, DQMStore::IGetter &) override;

private:

  MonitorElement * h_eventBased;
  MonitorElement * h_LBbased;
  

};

#endif
