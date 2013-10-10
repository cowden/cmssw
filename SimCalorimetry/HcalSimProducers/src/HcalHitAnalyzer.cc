#include "SimCalorimetry/HcalSimProducers/src/HcalHitAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include<iostream>


HcalHitAnalyzer::HcalHitAnalyzer(edm::ParameterSet const& conf) 
: simParameterMap_(conf),
  hbheFilter_(),
  hoFilter_(),
  hfFilter_(true),
  zdcFilter_(),
  hbheAnalyzer_("HBHE", 1., &simParameterMap_, &hbheFilter_),
  hoAnalyzer_("HO", 1., &simParameterMap_, &hoFilter_),
  hfAnalyzer_("HF", 1., &simParameterMap_, &hfFilter_),
  zdcAnalyzer_("ZDC", 1., &simParameterMap_, &zdcFilter_)
{
  tok_hbhe_ = consumes<HBHERecHitCollection>(conf.getParameter<edm::InputTag>("hbheRecHitCollectionTag"));
  tok_ho_ = consumes<HORecHitCollection>(conf.getParameter<edm::InputTag>("hoRecHitCollectionTag"));
  tok_hf_ = consumes<HFRecHitCollection>(conf.getParameter<edm::InputTag>("hfRecHitCollectionTag"));
  tok_cf_ = consumes<CrossingFrame<PCaloHit> >(edm::InputTag("mix", "g4SimHitsHcalHits"));
}


namespace HcalHitAnalyzerImpl {
  template<class Collection>
  void analyze(edm::Event const& e, CaloHitAnalyzer & analyzer, edm::EDGetTokenT<Collection>& tag) {
    edm::Handle<Collection> recHits;
    e.getByToken(tag, recHits);
    for(unsigned i = 0 ; i < recHits->size(); ++i) {
      analyzer.analyze((*recHits)[i].id().rawId(), (*recHits)[i].energy());
    }
  }
}


void HcalHitAnalyzer::analyze(edm::Event const& e, edm::EventSetup const& c) {
   // Step A: Get Inputs
  edm::Handle<CrossingFrame<PCaloHit> > cf, zdccf;
  e.getByToken(tok_cf_,cf);
  //e.getByLabel("mix", "ZDCHits", zdccf);

  // test access to SimHits for HcalHits and ZDC hits
  std::auto_ptr<MixCollection<PCaloHit> > hits(new MixCollection<PCaloHit>(cf.product()));
  //std::auto_ptr<MixCollection<PCaloHit> > zdcHits(new MixCollection<PCaloHit>(zdccf.product()));
  hbheAnalyzer_.fillHits(*hits); 
  //hoAnalyzer_.fillHits(*hits);
  //hfAnalyzer_.fillHits(*hits);
  //zdcAnalyzer_.fillHits(*hits);
  HcalHitAnalyzerImpl::analyze<HBHERecHitCollection>(e, hbheAnalyzer_, tok_hbhe_);
  HcalHitAnalyzerImpl::analyze<HORecHitCollection>(e, hoAnalyzer_, tok_ho_);
  HcalHitAnalyzerImpl::analyze<HFRecHitCollection>(e, hfAnalyzer_, tok_hf_);
  //HcalHitAnalyzerImpl::analyze<ZDCRecHitCollection>(e, zdcAnalyzer_);
}
