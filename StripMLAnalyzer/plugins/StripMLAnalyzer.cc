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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"
#include "TVector3.h"

using reco::TrackCollection;

class StripMLAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit StripMLAnalyzer(const edm::ParameterSet&);
  ~StripMLAnalyzer() override;

private:

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> clusterCollectionToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_; 
  edm::InputTag inputTagClusters;


  TTree* clusterTree;

  edm::Service<TFileService> fs;
  edm::EventNumber_t eventN;
  int runN;
  int lumi;

  uint32_t    detId_;
  uint16_t    size;
  
  int         charge;

  const static int nMax = 800000;
  uint16_t    adc[nMax];
  uint16_t    sig_adc[nMax];
  uint16_t    bkg_adc[nMax];
  
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
  
  //now do what ever initialization is needed
  inputTagClusters  = iConfig.getParameter<edm::InputTag>("siStripClustersTag");
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  clusterCollectionToken_  = consumes<edmNew::DetSetVector<SiStripCluster>>(inputTagClusters);

  tkGeomToken_ = esConsumes();

  usesResource("TFileService");

  clusterTree = fs->make<TTree>("clusterTree","clusterTree");
  clusterTree->Branch("event", &eventN, "event/i");
  clusterTree->Branch("run",   &runN, "run/I");
  clusterTree->Branch("lumi",  &lumi, "lumi/I");

  clusterTree->Branch("detId", &detId_, "detId/i");
  clusterTree->Branch("charge", &charge, "charge/I");
  clusterTree->Branch("size", &size, "size/s");
  clusterTree->Branch("adc", adc, "adc[size]/s");

  clusterTree->Branch("sig_adc", sig_adc, "sig_adc[size]/s");
  clusterTree->Branch("bkg_adc", bkg_adc, "bkg_adc[size]/s");
   
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

  //  edm::Handle<edmNew::DetSetVector<SiStripCluster>> clusterCollection = iEvent.getHandle(clusterCollectionToken_);

  const auto& clusterCollection = iEvent.getHandle(clusterCollectionToken_);
  const auto& tracksHandle = iEvent.getHandle(tracksToken_);

  if (!clusterCollection.isValid()){
    edm::LogError("BadRefCore") << "Cluster handle is invalid!";
  }

  std::map<DetId,SiStripRecHit1D::ClusterRef> MapRecoHits;

  for (const auto& track: *tracksHandle){
    
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

    eventN = iEvent.id().event();
    runN   = (int) iEvent.id().run();
    lumi   = (int) iEvent.id().luminosityBlock();
    
    DetId detId = detSiStripClusters.detId();

    int subDet = detId.subdetId();
    for (const auto& stripCluster : detSiStripClusters) {

      detId_ = detId;
      charge = stripCluster.charge();
      size   = stripCluster.size();
      
      
      const auto& _detId = detId; // for the capture clause in the lambda function
      auto det = std::find_if(tkDets.begin(), tkDets.end(), [_detId](auto& elem) -> bool {
	return (elem->geographicalId().rawId() == _detId);
      });
      const StripTopology& p = dynamic_cast<const StripGeomDetUnit*>(*det)->specificTopology();

      auto mapid = MapRecoHits.find(detId);

      if(mapid == MapRecoHits.end()) continue;
      
      LocalPoint tar_lp = p.localPosition((float)  stripCluster.barycenter());
      LocalPoint ref_lp = p.localPosition((float)  MapRecoHits[detId]->barycenter());
      
      GlobalPoint tar_gp = (tkGeom->idToDet(detId))->surface().toGlobal(tar_lp);
      GlobalPoint ref_gp = (tkGeom->idToDet(detId))->surface().toGlobal(ref_lp);

      //std::cout<<tar_gp.x()<<"  "<<ref_gp.x()<<std::endl;
      for (int strip = stripCluster.firstStrip(); strip < stripCluster.endStrip(); ++strip){

	adc [strip - stripCluster.firstStrip()] = stripCluster[strip - stripCluster.firstStrip()];

      }

      bool matched = false; 
      if(std::fabs(tar_lp.x() - ref_lp.x()) < 0.01 &&  (subDet == SiStripDetId::TIB || subDet == SiStripDetId::TOB)){
	matched = true;
      }
      else if(std::fabs(tar_gp.z() - ref_gp.z()) <= 0.0001 && (subDet == SiStripDetId::TID || subDet == SiStripDetId::TEC)){
	matched = true;
      }
      else{
	matched = false;
      }
      if (matched){
	for (int strip = stripCluster.firstStrip(); strip < stripCluster.endStrip(); ++strip){
	  sig_adc [strip - stripCluster.firstStrip()] = stripCluster[strip - stripCluster.firstStrip()];
	}
      }
      else{
	for (int strip = stripCluster.firstStrip(); strip < stripCluster.endStrip(); ++strip){
	  bkg_adc [strip - stripCluster.firstStrip()] = stripCluster[strip - stripCluster.firstStrip()];
	}
      }
      clusterTree->Fill();
    }
  }
  
  

}

//define this as a plug-in
DEFINE_FWK_MODULE(StripMLAnalyzer);
