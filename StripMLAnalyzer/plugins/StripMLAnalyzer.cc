// -*- C++ -*-
//
// Package:    StripMLStudies/StripMLAnalyzer
// Class:      StripMLAnalyzer
//
/**\class StripMLAnalyzer StripMLAnalyzer.cc StripMLStudies/StripMLAnalyzer/plugins/StripMLAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vinaya Krishnan Nair
//         Created:  Wed, 25 Jun 2025 14:47:42 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class StripMLAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit StripMLAnalyzer(const edm::ParameterSet&);
  ~StripMLAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> clusterCollectionToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_; 
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
StripMLAnalyzer::StripMLAnalyzer(const edm::ParameterSet& iConfig) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
   tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
   clusterCollectionToken_  = consumes<edmNew::DetSetVector<SiStripCluster>>(iConfig.getParameter<edm::InputTag>("siStripClustersTag"));

   tkGeomToken_ = esConsumes();

}

StripMLAnalyzer::~StripMLAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void StripMLAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  const auto& tkGeom = &iSetup.getData(tkGeomToken_);
  const auto tkDets = tkGeom->dets();

  std::map<DetId,SiStripRecHit1D::ClusterRef> MapRecoHits;

  for (const auto& track: tracks){
    s
    for (auto const &hit : track.recHits()){
      DetId detid = hit->geographicalId();
      int subDet = detid.subdetId();

      bool hitInStrip = (subDet == SiStripDetId::TIB) || (subDet == SiStripDetId::TID) || (subDet == SiStripDetId::TOB) || (subDet == SiStripDetId::TEC);

      if (!hitInStrip) continue;
      
      const std::type_info &type = typeid(*hit);
      if (type == typeid(SiStripRecHit1D)){
	const SiStripRecHit1D *striphit = dynamic_cast<const SiStripRecHit1D *>(hit);
  	  if(striphit != nullptr){

	    SiStripRecHit1D::ClusterRef stripclust(striphit->cluster());

	    MapRecoHits[detid] = stripclust;
	  }
      }
    }
  }

  for (const auto& detSiStripClusters : *clusterCollection) {
    DetId detId = detSiStripClusters.detId();
    int subDet = detId.subdetId();
    for (const auto& stripCluster : detSiStripClusters) {
      GlobalPoint tar_gp = (tkGeom->idToDet(detId))->surface().toGlobal(p.localPosition((float) stripCluster.barycenter()));
      GlobalPoint ref_gp = (tkGeom->idToDet(detId))->surface().toGlobal(p.localPosition((float) MapRecoHits[detId]->barycenter()));

      
      
      
    }
  }
  
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void StripMLAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void StripMLAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void StripMLAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(StripMLAnalyzer);
