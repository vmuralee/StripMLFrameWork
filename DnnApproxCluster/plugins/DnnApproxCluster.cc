// -*- C++ -*-
//
// Package:    StripMLStudies/DnnApproxCluster
// Class:      DnnApproxCluster
//
/**\class DnnApproxCluster DnnApproxCluster.cc StripMLStudies/DnnApproxCluster/plugins/DnnApproxCluster.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vinaya Krishnan Nair
//         Created:  Thu, 17 Jul 2025 14:55:09 GMT
//
//

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/FileInPath.h"

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripApproximateCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripApproximateClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "tensorflow/core/framework/tensor.h"

class DnnApproxCluster : public edm::stream::EDProducer<> {
public:
  explicit DnnApproxCluster(const edm::ParameterSet&);
  void produce(edm::Event&, const edm::EventSetup&) override;
  ~DnnApproxCluster() override;


private:
  std::unique_ptr<tensorflow::Session> session_;
  tensorflow::GraphDef* graphDef_;

   // Event Setup Data
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;

  uint16_t    firstStrip;
  uint16_t    endStrip;

  const static int nMax = 8000000;
  float       hitX[nMax];
  float       hitY[nMax];
  float       hitZ[nMax];
  
  uint16_t    adc[nMax];

  uint32_t  detId;
  uint16_t  size;

  float dnn_score;
  
  edm::EDGetTokenT<SiStripApproximateClusterCollection> approxClusterToken;

  std::string inputTensorName_;
  std::string outputTensorName_;

  double nn_threshold;
};


DnnApproxCluster::DnnApproxCluster(const edm::ParameterSet& iConfig) {
  auto modelPath = iConfig.getParameter<std::string>("modelPath");
  approxClusterToken  = consumes<SiStripApproximateClusterCollection>(iConfig.getParameter<edm::InputTag>("approxSiStripClustersTag"));

  graphDef_ = tensorflow::loadGraphDef(modelPath);
  session_.reset(tensorflow::createSession(graphDef_));

  tkGeomToken_ = esConsumes();

  inputTensorName_ = iConfig.getParameter<std::string>("inputTensorName");
  outputTensorName_ = iConfig.getParameter<std::string>("outputTensorName");

  nn_threshold = iConfig.getParameter<double>("NN_threshold");

  //produces<edmNew::DetSetVector<SiStripApproximateCluster>>();

  produces<SiStripApproximateClusterCollection>();
}

DnnApproxCluster::~DnnApproxCluster() {
  session_.reset();
  delete graphDef_;
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void DnnApproxCluster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Handle<SiStripApproximateClusterCollection>  approxClusterCollection = iEvent.getHandle(approxClusterToken);
  //auto result = std::make_unique<edmNew::DetSetVector<SiStripApproximateCluster>>();
  auto result = std::make_unique<SiStripApproximateClusterCollection>();
  

  const auto& tkGeom = &iSetup.getData(tkGeomToken_);
  const auto tkDets = tkGeom->dets();

  float mean_adc, std_adc;
  float hitx,hity,hitz;
  for (const auto& detApproxClusters : *approxClusterCollection) {
    //edmNew::DetSetVector<SiStripApproximateCluster>::FastFiller ff{*result, detApproxClusters.id()};
    auto ff = result->beginDet(detApproxClusters.id());
    for (const auto& cluster : detApproxClusters) {
      detId = detApproxClusters.id();

      uint16_t nStrips{0};
      
      const auto& _detId = detApproxClusters.id(); // for the capture clause in the lambda function
      auto det = std::find_if(tkDets.begin(), tkDets.end(), [_detId](auto& elem) -> bool {
	  return (elem->geographicalId().rawId() == _detId);
	});
      const StripTopology& p = dynamic_cast<const StripGeomDetUnit*>(*det)->specificTopology();
      nStrips = p.nstrips() - 1;
      const auto convertedCluster = SiStripCluster(cluster, nStrips);
      
      firstStrip = convertedCluster.firstStrip();
      endStrip   = convertedCluster.endStrip();

      size_t strip_size = 0;
      size = convertedCluster.size();
      for (int strip = firstStrip; strip < endStrip+1; ++strip)
	{
	  GlobalPoint gp = (tkGeom->idToDet(detApproxClusters.id()))->surface().toGlobal(p.localPosition((float) strip));
	  hitX   [strip - firstStrip] = gp.x();
	  hitY   [strip - firstStrip] = gp.y();
	  hitZ   [strip - firstStrip] = gp.z();
	  adc    [strip - firstStrip] = convertedCluster[strip - firstStrip];
	  ++strip_size;
	}


      auto meanADC = [](uint16_t* data, size_t size) {
	return static_cast<double>(std::accumulate(data,data + size,0))/size;
      };

      auto stddevADC = [meanADC](uint16_t * data, int size) {
        double m = meanADC(data, size);
        double sum_sq_diff = 0.0;
        for (int i = 0; i < size; ++i)
	  sum_sq_diff += std::pow(data[i] - m, 2);
        return std::sqrt(sum_sq_diff / size); // Population stddev
      };
      mean_adc = meanADC(adc,strip_size);
      std_adc  = stddevADC(adc,strip_size);

      auto mean = [](float* data, size_t size) {
	return static_cast<double>(std::accumulate(data,data + size,0))/size;
      };

      hitx = mean(hitX, strip_size);      
      hity = mean(hitY,strip_size);
      hitz = mean(hitZ,strip_size);

      tensorflow::Tensor input(tensorflow::DT_FLOAT, {1, 7});

      auto inputMap = input.tensor<float, 2>();
      inputMap(0,0) = mean_adc;
      inputMap(0,1) = std_adc;
      inputMap(0,2) = size;
      inputMap(0,3) = hitx;
      inputMap(0,4) = hity;
      inputMap(0,5) = detId/10000000;
      inputMap(0,6) = hitz;

      std::vector<tensorflow::Tensor> outputs;

      tensorflow::run(session_.get(),
		      {{inputTensorName_, input}},
		      {outputTensorName_},
		      &outputs);
      
      dnn_score = outputs[0].matrix<float>()(0, 0);
      
      if (dnn_score < nn_threshold) continue;
      // float barycenter = cluster.barycenter();
      // uint8_t width = cluster.width();
      // float avgCharge = cluster.avgCharge();
      // bool filter = cluster.filter();
      // bool isSaturated = cluster.isSaturated();
      
      // ff.push_back(SiStripApproximateCluster(barycenter, width, avgCharge, filter, isSaturated));

      ff.push_back(cluster);
    }
    
  }  

  iEvent.put(std::move(result));
  
  // Prepare input tensor
  
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(DnnApproxCluster);
