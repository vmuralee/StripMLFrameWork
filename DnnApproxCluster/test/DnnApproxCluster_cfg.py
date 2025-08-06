import FWCore.ParameterSet.Config as cms

process = cms.Process("nnHLT")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:/eos/home-v/vmuralee/PREanalysis/clusterMLstudies/RelValSinglePionNoPURawPrime1.root"))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'140X_mcRun3_2024_realistic_v21' , '')#'140X_dataRun3_Prompt_v3'


process.hltSiStripClusters2ApproxClusters = cms.EDProducer("DnnApproxCluster",
                                      modelPath = cms.string("/afs/cern.ch/work/v/vmuralee/private/PREanalysis/StripClusterMLStudies/CMSSW_14_1_0_pre6/src/ApproxDnnCluster/DnnApproxCluster/saved_model/frozen_graph_v1.pb"),
                                      approxSiStripClustersTag = cms.InputTag("hltSiStripClusters2ApproxClustersold"),

                                      inputTensorName = cms.string("input_layer:0"),
                                      outputTensorName = cms.string("StatefulPartitionedCall/sequential_1/dense_3_1/Sigmoid:0"),

                                      NN_threshold = cms.double(0.4)

                                      )



process.output_step = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('ZLIB'),
    compressionLevel = cms.untracked.int32(7),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('file:/eos/home-v/vmuralee/PREanalysis/clusterMLstudies/SinglePionApproxDnn.root'),
    # outputCommands = cms.untracked.vstring(
    #     'drop *',
    #     'keep *_*generalTracks*_*_*',
    #     'keep DetIds_hltSiStripRawToDigi_*_HLTX',
    #     'keep FEDRawDataCollection_*_*_*',
    #     'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_HLTX',
    #     'keep *_*_*_nnHLT',
    #     'keep triggerTriggerEvent_*_*_HLTX'
    # )
    outputCommands = cms.untracked.vstring('drop *_hltSiStripClusters2ApproxClustersold_*_HLTX',
                                           'keep *_*_*_*')
    # 'keep int_*_*_*',
    #   'keep *_*_*_nnHLT',
    #   'keep DetIds_hltSiStripRawToDigi_*_HLTX',
    #   'keep FEDRawDataCollection_*_*_*',
    #   'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_HLTX',
    #   'keep edmTriggerResults_*_*_HLTX',
    #   'keep triggerTriggerEvent_*_*_HLTX')
)

#process.schedule = cms.Schedule( *(process.dnnAproxclusProducer, process.output_step))

process.p = cms.Path(process.hltSiStripClusters2ApproxClusters)


process.e = cms.EndPath(process.output_step)
