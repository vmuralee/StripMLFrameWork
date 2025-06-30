import FWCore.ParameterSet.Config as cms

process = cms.Process("tupleProducer")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_14_1_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v21_STD_Recycled_2024_noPU-v1/2580000/a65a9d3e-c48b-4316-9947-64c1906fc136.root'

    ))

# Other statements
#process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'140X_mcRun3_2024_realistic_v21' , '')#'140X_dataRun3_Prompt_v3'

process.stripmlAnalyzer = cms.EDAnalyzer("StripMLAnalyzer",
                                tracks = cms.InputTag('generalTracks','','RECO'),
                                siStripClustersTag = cms.InputTag('siStripClusters','','RECO'),
                                )


process.TFileService = cms.Service("TFileService", fileName=cms.string('Test.root'))
process.p = cms.Path(process.stripmlAnalyzer)
