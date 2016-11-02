import FWCore.ParameterSet.Config as cms

process = cms.Process("SimpleMuonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/afs/cern.ch/user/g/goni/public/ForMihee/step3_25ns_RAW2DIGI_L1Reco_RECO_1.root'
    )
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(
        'file:./result_pA_EPOS.root'
    )
)

# Global Tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_pA_v1', '')


process.hionia = cms.EDAnalyzer('MuonAnalyzer',
    muonType         = cms.string("GLB"),
    oniaPDG          = cms.int32(443),
    isMC             = cms.bool(True),
    isHI             = cms.bool(False),
    checkIDCuts      = cms.bool(True),
    checkAcceptance  = cms.bool(True),
    srcMuon          = cms.InputTag("muons"),
    genParticles     = cms.InputTag("genParticles"),
    primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
    CentralitySrc    = cms.InputTag("hiCentrality"),
    CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
)


process.p = cms.Path(process.hionia)
