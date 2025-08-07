// -*- C++ -*-
//
// Package:    ApproxDnnCluster/TrackEfficiencyAnalyzer
// Class:      TrackEfficiencyAnalyzer
//
/**\class TrackEfficiencyAnalyzer TrackEfficiencyAnalyzer.cc ApproxDnnCluster/TrackEfficiencyAnalyzer/plugins/TrackEfficiencyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vinaya Krishnan Nair
//         Created:  Mon, 04 Aug 2025 20:04:16 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
//#include "SimTracker/TrackAssociation/interface/trackAssociationByHits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class TrackEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackEfficiencyAnalyzer(const edm::ParameterSet&);
  ~TrackEfficiencyAnalyzer() override;



private:

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  void endJob() override;
  // ----------member data ---------------------------

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif

  edm::EDGetTokenT<reco::TrackCollection> recoTracksToken_;
  edm::EDGetTokenT<TrackingParticleCollection> simTracksToken_;
  double maxDeltaR_;
  double minSimPt_;

  TH1F* h_simPt;
  TH1F* h_matchedSimPt;
  TH1F* h_efficiency;
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
TrackEfficiencyAnalyzer::TrackEfficiencyAnalyzer(const edm::ParameterSet& iConfig){
  //now do what ever initialization is needed
  recoTracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("recoTracks"));
  simTracksToken_ = consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("simTracks"));
  maxDeltaR_ = iConfig.getParameter<double>("maxDeltaR");
  minSimPt_ = iConfig.getParameter<double>("minSimPt");

   edm::Service<TFileService> fs;

    h_simPt         = fs->make<TH1F>("h_simPt", "Simulated Tracks pT; pT [GeV]; Entries", 10, 0, 10);
    h_matchedSimPt  = fs->make<TH1F>("h_matchedSimPt", "Matched Simulated Tracks pT; pT [GeV]; Entries", 10, 0, 10);
    h_efficiency    = fs->make<TH1F>("h_efficiency", "Track Efficiency; pT [GeV]; Efficiency", 10, 0, 10);
}

TrackEfficiencyAnalyzer::~TrackEfficiencyAnalyzer() {

}

//
// member functions
//

// ------------ method called for each event  ------------
void TrackEfficiencyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  
  edm::Handle<reco::TrackCollection> recoTracks;
  iEvent.getByToken(recoTracksToken_, recoTracks);

  edm::Handle<TrackingParticleCollection> simTracks;
  iEvent.getByToken(simTracksToken_, simTracks);


  if (!simTracks.isValid()){
    edm::LogError("TrackEfficiency") << "Sim track collection not found!";
    return;
  }

  if (!recoTracks.isValid()){
    edm::LogError("TrackEfficiency") << "Reco track collection not found!";
    return;
  }


  // // Convert recoTracks to RefVector
  // edm::RefToBaseVector<reco::Track> trackRefs;
  // for (size_t i = 0; i < recoTracks->size(); ++i) {
  //   edm::Ref<reco::TrackCollection> ref(recoTracks, i);
  //   edm::RefToBase<reco::Track> baseRef(ref);
  //   trackRefs.push_back(baseRef);
  // }

  // // Convert simTracks to RefVector
  // edm::RefVector<TrackingParticleCollection> simTrackRefs;
  // for (size_t i = 0; i < simTracks->size(); ++i) {
  //   simTrackRefs.push_back(edm::Ref<TrackingParticleCollection>(simTracks, i));
  // }

  // // Association
  // reco::SimToRecoCollection simToReco = associator.associateSimToReco(trackRefs, simTrackRefs);


  int totalSim = 0;
  int matchedSim = 0;

  for (const auto& sim : *simTracks) {
    if (sim.charge() == 0) continue;  // skip neutrals
    if (sim.pt() < minSimPt_) continue;

    totalSim++;
    h_simPt->Fill(sim.pt());
    // Simulated track info
    float simEta = sim.eta();
    float simPhi = sim.phi();
    float simPt = sim.pt();

    bool matched = false;

    for (const auto& reco : *recoTracks) {
      if (reco.charge() == 0) continue;

      double dR = deltaR(simEta, simPhi, reco.eta(), reco.phi());

      if (dR < maxDeltaR_) {
	matched = true;
	break;
      }
    }

    if (matched){
      h_matchedSimPt->Fill(simPt);
      matchedSim++;
    }
  }

  double efficiency = (totalSim > 0) ? (double)matchedSim / totalSim : 0.0;

  std::cout<< "Simulated tracks: " << totalSim << ", Matched: " << matchedSim << ", Efficiency: " << efficiency<<std::endl;
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}
void TrackEfficiencyAnalyzer::endJob()
{
    // Compute efficiency histogram from matched / total
    h_efficiency->Divide(h_matchedSimPt, h_simPt, 1.0, 1.0, "B"); // binomial errors
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(TrackEfficiencyAnalyzer);
