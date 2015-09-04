import FWCore.ParameterSet.Config as cms

processName = "MuonSuite"
process = cms.Process(processName)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring("file:/data_CMS/cms/miheejo/ROOTFILES/PyquenEvtGen_jpsiMuMu_JPsiPt912/GEN-SIM-RECODEBUG/004E94BB-10FE-E211-B319-008CFA008C94.root")
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

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ValidationHeavyIons_cff')
process.load('Configuration.StandardSequences.ReMixingSeeds_cff')
process.load('DQMOffline.Configuration.DQMOfflineHeavyIonsMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff")

process.GlobalTag.globaltag = 'STARTHI44_V12::All'

#---- For centrality ----#
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

#---- Validation stuffs ----#
## Default validation modules
from Validation.RecoHI.muonValidationHeavyIons_cff import *

process.tpToStaUpdMuonAssociation.simtracksTag = cms.InputTag("hiSignalG4SimHits")
process.tpToStaUpdMuonAssociation.CSCsimHitsTag = cms.InputTag("hiSignalG4SimHits","MuonCSCHits")
process.tpToStaUpdMuonAssociation.DTsimhitsTag = cms.InputTag("hiSignalG4SimHits","MuonDTHits")
process.tpToStaUpdMuonAssociation.RPCsimhitsTag = cms.InputTag("hiSignalG4SimHits","MuonRPCHits")
process.cutsTpMuons.ptMin = cms.double(0)
process.cutsRecoTrkMuons.ptMin = cms.double(0)
process.cutsRecoTracks.ptMin = cms.double(0)
process.staUpdMuonTrackVMuonAssoc.outputFile = cms.string("file:./validation.root")
process.staUpdMuonTrackVMuonAssoc.srcCentrality = cms.InputTag("hiCentrality");
process.staUpdMuonTrackVMuonAssoc.centralityRanges = cms.vdouble(0,5,10,15,20,30,40,50,60,100);
process.staUpdMuonTrackVMuonAssoc.nint = cms.int32(25)
process.staUpdMuonTrackVMuonAssoc.min = cms.double(0)
process.staUpdMuonTrackVMuonAssoc.useFabsEta = cms.bool(True)
process.staUpdMuonTrackVMuonAssoc.rapArr = cms.vdouble(0, 0.8, 1.6, 2.4)
process.staUpdMuonTrackVMuonAssoc.rapArrRes = cms.vdouble(0, 0.8, 1.6, 2.4)
process.staUpdMuonTrackVMuonAssoc.ptArrRes = cms.vdouble(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20)
# if pTRanges is defined, nintpT, minpT and maxpT will not be used
process.staUpdMuonTrackVMuonAssoc.minpT = cms.double(0)
process.staUpdMuonTrackVMuonAssoc.maxpT = cms.double(20)
process.staUpdMuonTrackVMuonAssoc.nintpT = cms.int32(200)
# if pTRanges is defined, nintpT, minpT and maxpT will not be used
process.staUpdMuonTrackVMuonAssoc.pTRanges = cms.vdouble(0,1.5,3,4.5,6,7.5,9,10,12,14,20);
#process.staUpdMuonTrackVMuonAssoc.pTRanges = cms.vdouble(0,1.5,3.5,5.5,7.5,10,14,20);
process.staUpdMuonTrackVMuonAssoc.pTWeight = cms.double(2.62048e-6);
'''
pT03 = 8.50977e-6;
pT36 = 7.69819e-6;
pT69 = 1.24059e-6;
pT912 = 2.62048e-6;
pT1215 = 7.82144e-6;
pT1530 = 3.46790e-6;
'''

process.prevalidation_ = cms.Sequence(process.cutsRecoTracks+process.cutsRecoTrkMuons+process.cutsTpMuons+process.tpToStaUpdMuonAssociation);
process.validation_ = cms.Sequence(process.staUpdMuonTrackVMuonAssoc);
process.DQMOfflineHeavyIons_ = cms.Sequence(process.muonAnalyzer);
process.muonValidation_step = cms.Path(cms.SequencePlaceholder("mix")+process.recoMuonValidation)

process.prevalidation_step = cms.Path(process.prevalidation_)
process.dqmoffline_step = cms.Path(process.DQMOfflineHeavyIons_)
process.validation_step = cms.EndPath(process.validation_)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.prevalidation_step,
    process.validation_step,
#    process.dqmoffline_step,
#    process.muonValidation_step,
    process.endjob_step)
#    process.endjob_step,process.DQMoutput_step)

