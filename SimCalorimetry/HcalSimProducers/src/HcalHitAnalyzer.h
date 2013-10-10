#ifndef HcalSimProducers_HcalHitAnalyzer_h
#define HcalSimProducers_HcalHitAnalyzer_h

#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HBHEHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HOHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HFHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/ZDCHitFilter.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include <string>

/** Compares HCAL RecHits to SimHit

  \Author Rick Wilkinson, Caltech
*/


class HcalHitAnalyzer : public edm::EDAnalyzer
{
public:

  explicit HcalHitAnalyzer(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& c);


private:
  HcalSimParameterMap simParameterMap_;
  HBHEHitFilter hbheFilter_;
  HOHitFilter hoFilter_;
  HFHitFilter hfFilter_;
  ZDCHitFilter zdcFilter_;
  CaloHitAnalyzer hbheAnalyzer_;
  CaloHitAnalyzer hoAnalyzer_;
  CaloHitAnalyzer hfAnalyzer_;
  CaloHitAnalyzer zdcAnalyzer_;

  edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HORecHitCollection> tok_ho_;
  edm::EDGetTokenT<HFRecHitCollection> tok_hf_;
  edm::EDGetTokenT<CrossingFrame<PCaloHit> > tok_cf_;
};

#endif
