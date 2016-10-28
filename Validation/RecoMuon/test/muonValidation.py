import FWCore.ParameterSet.Config as cms

processName = "MuonSuite"
_isPbPb = False
sampleweight = 1 # 1 for pp, numbers vary for pbpb.

'''
pT03 = 8.50977e-6;
pT36 = 7.69819e-6;
pT69 = 1.24059e-6;
pT912 = 2.62048e-6;
pT1215 = 7.82144e-6;
pT1530 = 3.46790e-6;
'''

###############################################################################
process = cms.Process(processName)
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.ValidationHeavyIons_cff')
process.load('DQMOffline.Configuration.DQMOfflineHeavyIonsMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff")

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames =
    cms.untracked.vstring(
      "file:./HIN-HINppWinter16DR-00111_99.root",
#      "file:/afs/cern.ch/work/m/miheejo/private/cmssw758patch4/src/reco/Pythia8_JpsiMM_ptJpsi_00_03_Hydjet_MB/HIN-HINPbPbWinter16DR-00028_reco.root"
      )
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories = ['TrackAssociator', 'TrackValidator']
process.MessageLogger.debugModules = ['*']
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
    TrackAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
    TrackValidator = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    )
)
process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)


#---- For centrality ----#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v13', '')
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")

##---- Validation stuffs ----##
## Default validation modules
from Validation.RecoHI.muonValidationHeavyIons_cff import *
if (_isPbPb):
  process.cutsRecoTrkMuons.vertexTag = cms.InputTag("hiSelectedVertex","","")
  process.cutsRecoTracks.vertexTag = cms.InputTag("hiSelectedVertex","","")
else:
  process.cutsRecoTrkMuons.vertexTag = cms.InputTag("offlinePrimaryVertices","","")
  process.cutsRecoTracks.vertexTag = cms.InputTag("offlinePrimaryVertices","","")
  process.cutsRecoTracks.src = cms.InputTag("generalTracks","","")
  process.cutsRecoTrkMuons.src = cms.InputTag("generalTracks","","")

process.cutsTpMuons.ptMin = cms.double(0)
process.cutsRecoTrkMuons.ptMin = cms.double(0)
process.cutsRecoTracks.ptMin = cms.double(0)
process.staUpdMuonTrackVMuonAssoc.isPbPb = cms.bool(_isPbPb);
process.staUpdMuonTrackVMuonAssoc.outputFile = cms.string("file:./validation.root")
process.staUpdMuonTrackVMuonAssoc.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers");
process.staUpdMuonTrackVMuonAssoc.centralityRanges = cms.vdouble(0,5,10,15,20,30,40,50,60,100);
process.staUpdMuonTrackVMuonAssoc.nint = cms.int32(25)
process.staUpdMuonTrackVMuonAssoc.min = cms.double(0)
process.staUpdMuonTrackVMuonAssoc.useFabsEta = cms.bool(True)
process.staUpdMuonTrackVMuonAssoc.rapArr = cms.vdouble(0, 1.2, 2.1, 2.4)
process.staUpdMuonTrackVMuonAssoc.rapArrRes = cms.vdouble(0, 1.2, 2.1, 2.4)
process.staUpdMuonTrackVMuonAssoc.pTWeight = cms.double(sampleweight);
process.staUpdMuonTrackVMuonAssoc.ptArr = cms.vdouble(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20)
process.staUpdMuonTrackVMuonAssoc.ptArrRes = cms.vdouble(0, 1.5, 3, 4.5, 6, 7.5, 9, 10, 12, 14, 20)

process.prevalidation_ = cms.Sequence(process.cutsRecoTracks+process.cutsRecoTrkMuons+process.cutsTpMuons+process.tpToStaUpdMuonAssociation);
process.validation_ = cms.Sequence(process.staUpdMuonTrackVMuonAssoc);

process.prevalidation_step = cms.Path(process.prevalidation_)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.validation_step = cms.EndPath(process.validation_)

process.schedule = cms.Schedule(
    process.prevalidation_step,
    process.validation_step,
    process.endjob_step)

