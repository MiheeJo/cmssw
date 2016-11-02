import FWCore.ParameterSet.Config as cms

process = cms.Process("SimpleMuonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #aodsim official pbpb sample(HINPbPbWinter16DR/Pythia8_JpsiMM_ptJpsi_12_15_Hydjet_MB/AODSIM)
      'file:/afs/cern.ch/work/m/miheejo/public/HitByHitMatching/Test/002033E1-69E8-E511-843B-001E67E34165.root'
#        'file:/afs/cern.ch/work/m/miheejo/public/HitByHitMatching/Test/HIN-HINPbPbWinter16DR-00028_reco.root'
    )
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(
        'file:./result.root'
    )
)

# Global Tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v13', '')


process.hionia = cms.EDAnalyzer('MuonAnalyzer',
    muonType         = cms.string("GLB"),
    oniaPDG          = cms.int32(443),
    isMC             = cms.bool(True),
    isHI             = cms.bool(True),
    checkIDCuts      = cms.bool(True),
    checkAcceptance  = cms.bool(True),
    srcMuon          = cms.InputTag("muons"),
    genParticles     = cms.InputTag("genParticles"),
    primaryVertexTag = cms.InputTag("hiSelectedVertex"),
    CentralitySrc    = cms.InputTag("hiCentrality"),
    CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
)


process.p = cms.Path(process.hionia)
