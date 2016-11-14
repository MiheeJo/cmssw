import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


#----------------------------------------------------------------------------

# Setup Settings for ONIA TREE:

isMC           = False    # if input is MONTECARLO: True or if it's DATA: False
applyMuonCuts  = False    # Apply muon ID quality cuts
muonSelection  = "Trk"    # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker) are availale

#----------------------------------------------------------------------------

# Print Onia Skim settings:
print( " " )
print( "[INFO] Settings used for ONIA TREE DATA: " )
print( "[INFO] isMC          = " + ("True" if isMC else "False") )
print( "[INFO] applyMuonCuts = " + ("True" if applyMuonCuts else "False") )
print( "[INFO] muonSelection = " + muonSelection )
print( " " )


# set up process
process = cms.Process("HIOnia")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "OniaTree.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = 'file:/afs/cern.ch/user/a/anstahll/UPDATE/ONIA/CMSSW_8_0_23/src/HiSkim/HiOnia2MuMu/test/60F56674-D8A4-E611-A91A-FA163EB6F0FA.root'
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.extend(["HiOnia2MuMuPAT_muonLessSizeORpvTrkSize"])
process.MessageLogger.cerr.HiOnia2MuMuPAT_muonLessSizeORpvTrkSize = cms.untracked.PSet( limit = cms.untracked.int32(5) )

# load the Geometry and Magnetic Field for the TransientTrackBuilder
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Express_v15', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(
    record = cms.string("HeavyIonRcd"),
    tag = cms.string("CentralityTable_HFtowersPlusTrunc200_EPOS5TeV_v80x01_mc"),
    connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
    label = cms.untracked.string("HFtowersPlusTruncEpos")
    ),
  cms.PSet(
    record = cms.string('L1TUtmTriggerMenuRcd'),
    tag = cms.string("L1Menu_HeavyIons2016_v2_m2_xml"),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
    ),
  cms.PSet(
    record = cms.string('L1TGlobalPrescalesVetosRcd'),
    tag = cms.string("L1TGlobalPrescalesVetos_Stage2v0_hlt"),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
    )
  )

###################### Onia Skim Producer #################################################


# HLT Dimuon Triggers
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# HLT pPb MENU:  /users/anstahll/PA2016/PAMuon2016Full/V3
process.hltOniaHI.HLTPaths =  [
  "HLT_PAL1DoubleMuOpen_v1",
  "HLT_PAL1DoubleMuOpen_OS_v1",
  "HLT_PAL1DoubleMuOpen_SS_v1",
  "HLT_PAL1DoubleMu0_v1",
  "HLT_PAL1DoubleMu0_MGT1_v1",
  "HLT_PAL1DoubleMu0_HighQ_v1",
  "HLT_PAL2DoubleMu0_v1",
  "HLT_PAL3DoubleMu0_v1",
  "HLT_PAL3DoubleMu0_HIon_v1",
  "HLT_PAL1DoubleMu10_v1",
  "HLT_PAL2DoubleMu10_v1",
  "HLT_PAL3DoubleMu10_v1",
  "HLT_PAL2Mu12_v1",
  "HLT_PAL2Mu15_v1",
  "HLT_PAL3Mu3_v1",
  "HLT_PAL3Mu5_v1",
  "HLT_PAL3Mu7_v1",
  "HLT_PAL3Mu12_v1",
  "HLT_PAL3Mu15_v1"
  ]

process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True
process.hltOniaHI.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=isMC, HLT="HLT", Filter=True, useL1Stage2=True)

### Temporal fix for the PAT Trigger prescale warnings.
process.patTriggerFull.l1GtReadoutRecordInputTag = cms.InputTag("gtDigis","","RECO")
process.patTriggerFull.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
process.patTriggerFull.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
process.patTriggerFull.getPrescales      = cms.untracked.bool(True)
###

##### Onia2MuMuPAT input collections/options
process.onia2MuMuPatGlbGlb.dimuonSelection          = cms.string("mass > 0")
process.onia2MuMuPatGlbGlb.resolvePileUpAmbiguity   = True
process.onia2MuMuPatGlbGlb.srcTracks                = cms.InputTag("generalTracks")
process.onia2MuMuPatGlbGlb.primaryVertexTag         = cms.InputTag("offlinePrimaryVertices")
process.patMuonsWithoutTrigger.pvSrc                = cms.InputTag("offlinePrimaryVertices")
# Adding muonLessPV gives you lifetime values wrt. muonLessPV only
process.onia2MuMuPatGlbGlb.addMuonlessPrimaryVertex = True
if isMC:
  process.genMuons.src = "genParticles"
  process.onia2MuMuPatGlbGlb.genParticles = "genParticles"

##### Dimuon pair selection
commonP1 = "|| (innerTrack.isNonnull && genParticleRef(0).isNonnull)"
commonP2 = " && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35"
if muonSelection == "Glb":
  highP = "isGlobalMuon"; # At least one muon must pass this selection
  process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("("+highP+commonP1+")"+commonP2)
  lowP = "isGlobalMuon"; # BOTH muons must pass this selection
  process.onia2MuMuPatGlbGlb.lowerPuritySelection = cms.string("("+lowP+commonP1+")"+commonP2)
elif muonSelection == "GlbTrk":
  highP = "(isGlobalMuon && isTrackerMuon)";
  process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("("+highP+commonP1+")"+commonP2)
  lowP = "(isGlobalMuon && isTrackerMuon)";
  process.onia2MuMuPatGlbGlb.lowerPuritySelection = cms.string("("+lowP+commonP1+")"+commonP2)
elif muonSelection == "Trk":
  highP = "isTrackerMuon";
  process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("("+highP+commonP1+")"+commonP2)
  lowP = "isTrackerMuon";
  process.onia2MuMuPatGlbGlb.lowerPuritySelection = cms.string("("+lowP+commonP1+")"+commonP2)
else:
  print "ERROR: Incorrect muon selection " + muonSelection + " . Valid options are: Glb, Trk, GlbTrk"

##### Remove few paths for MC
if isMC:
  process.patMuonSequence.remove(process.hltOniaHI)


###################### HiOnia Analyzer #################################################

### For Centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("pACentrality")
process.centralityBin.centralityVariable = cms.string("HFtowersPlusTrunc")
process.centralityBin.nonDefaultGlauberModel = cms.string("Epos")
process.EventAna_step = cms.Path( process.centralityBin )

process.hionia = cms.EDAnalyzer('HiOniaAnalyzer',
                                #-- Collections
                                srcMuon             = cms.InputTag("patMuonsWithTrigger"),     # Name of PAT Muon Collection
                                srcMuonNoTrig       = cms.InputTag("patMuonsWithoutTrigger"),  # Name of PAT Muon Without Trigger Collection
                                src                 = cms.InputTag("onia2MuMuPatGlbGlb"),      # Name of Onia Skim Collection
                                EvtPlane            = cms.InputTag("hiEvtPlane",""),           # Name of Event Plane Collection. For RECO use: hiEventPlane,recoLevel

                                triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"), # Label of Trigger Results

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                                useRapidity = cms.bool(True),

                                #--
                                maxAbsZ = cms.double(24.0),

                                pTBinRanges      = cms.vdouble(0.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 40.0),
                                etaBinRanges     = cms.vdouble(0.0, 2.5),
                                centralityRanges = cms.vdouble(20,40,100),

                                onlyTheBest        = cms.bool(False),
                                applyCuts          = cms.bool(applyMuonCuts),
                                selTightGlobalMuon = cms.bool(False),
                                storeEfficiency    = cms.bool(False),

                                removeSignalEvents = cms.untracked.bool(False),  # Remove/Keep signal events
                                removeTrueMuons    = cms.untracked.bool(False),  # Remove/Keep gen Muons
                                storeSameSign      = cms.untracked.bool(True),   # Store/Drop same sign dimuons

                                #-- Gen Details
                                oniaPDG = cms.int32(443),
                                muonSel = cms.string(muonSelection),
                                isHI = cms.untracked.bool(False),
                                isPA = cms.untracked.bool(True),
                                isMC = cms.untracked.bool(isMC),
                                isPromptMC = cms.untracked.bool(True),
                                useEvtPlane = cms.untracked.bool(False),
                                useGeTracks = cms.untracked.bool(False),
                                runVersionChange = cms.untracked.uint32(182133),

                                #-- Histogram configuration
                                combineCategories = cms.bool(False),
                                fillRooDataSet    = cms.bool(False),
                                fillTree          = cms.bool(True),
                                fillHistos        = cms.bool(False),
                                minimumFlag       = cms.bool(False),
                                fillSingleMuons   = cms.bool(True),
                                fillRecoTracks    = cms.bool(False),
                                histFileName      = cms.string(options.outputFile),
                                dataSetName       = cms.string(options.secondaryOutputFile),

                                # HLT pPb MENU:  /users/anstahll/PA2016/PAMuon2016Full/V3

                                dblTriggerPathNames = cms.vstring("HLT_PAL1DoubleMuOpen_v1",
                                                                  "HLT_PAL1DoubleMuOpen_OS_v1",
                                                                  "HLT_PAL1DoubleMuOpen_SS_v1",
                                                                  "HLT_PAL1DoubleMu0_v1",
                                                                  "HLT_PAL1DoubleMu0_MGT1_v1",
                                                                  "HLT_PAL1DoubleMu0_HighQ_v1",
                                                                  "HLT_PAL2DoubleMu0_v1",
                                                                  "HLT_PAL3DoubleMu0_v1",
                                                                  "HLT_PAL3DoubleMu0_HIon_v1",
                                                                  "HLT_PAL1DoubleMu10_v1",
                                                                  "HLT_PAL2DoubleMu10_v1",
                                                                  "HLT_PAL3DoubleMu10_v1"),

                                dblTriggerFilterNames = cms.vstring("hltL1fL1sDoubleMuOpenBptxANDL1Filtered0",
                                                                    "hltL1fL1sDoubleMuOpenOSBptxANDL1Filtered0",
                                                                    "hltL1fL1sDoubleMuOpenSSBptxANDL1Filtered0",
                                                                    "hltL1fL1sDoubleMu0BptxANDL1Filtered0",
                                                                    "hltL1fL1sDoubleMu0MassGT1BptxANDL1Filtered0",
                                                                    "hltL1fL1sDoubleMu0BptxANDL1HighQFiltered0",
                                                                    "hltL2fL1sDoubleMuOpenBptxANDL1f0L2Filtered0",
                                                                    "hltL3fL1sDoubleMuOpenBptxANDL1f0L2f0L3Filtered0",
                                                                    "hltHIL3fL1sDoubleMuOpenBptxANDL1f0L2f0L3Filtered0",
                                                                    "hltL1fL1sDoubleMu10BptxANDL1Filtered0",
                                                                    "hltL2fL1sDoubleMu10BptxANDL1f0L2Filtered10",
                                                                    "hltL3fL1sDoubleMu10BptxANDL1f0L2f10L3Filtered10"),

                                sglTriggerPathNames = cms.vstring("HLT_PAL2Mu12_v1",
                                                                  "HLT_PAL2Mu15_v1",
                                                                  "HLT_PAL3Mu3_v1",
                                                                  "HLT_PAL3Mu5_v1",
                                                                  "HLT_PAL3Mu7_v1",
                                                                  "HLT_PAL3Mu12_v1",
                                                                  "HLT_PAL3Mu15_v1"),

                                sglTriggerFilterNames = cms.vstring("hltL2fL1sSingleMu7BptxANDL1f0L2Filtered12",
                                                                    "hltL2fL1sSingleMu7BptxANDL1f0L2Filtered15",
                                                                    "hltL3fL1sSingleMu3BptxANDL1f0L2f0L3Filtered3",
                                                                    "hltL3fL1sSingleMu5BptxANDL1f0L2f0L3Filtered5",
                                                                    "hltL3fL1sSingleMu5BptxANDL1f0L2f0L3Filtered7",
                                                                    "hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12",
                                                                    "hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered15")
                                )

process.hionia.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")
process.hionia.genParticles     = cms.InputTag("genParticles")
process.hionia.muonLessPV       = cms.bool(True)
process.hionia.EvtPlane         = cms.InputTag("hiEvtPlaneFlat","")
process.hionia.CentralitySrc    = cms.InputTag("pACentrality")
process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowersPlusTrunc")
process.hionia.srcTracks        = cms.InputTag("generalTracks")

process.Onia2MuMuPAT = cms.Path(
  process.patMuonSequence *
  process.onia2MuMuPatGlbGlb *
  process.onia2MuMuPatGlbGlbFilter *
  process.hionia
  )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source.fileNames      = cms.untracked.vstring(options.inputFiles)
process.maxEvents             = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.schedule              = cms.Schedule( process.EventAna_step, process.Onia2MuMuPAT )
