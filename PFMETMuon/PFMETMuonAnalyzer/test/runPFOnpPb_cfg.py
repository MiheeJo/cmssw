### -------------------------------------------------------------------
### VarParsing allows one to specify certain parameters in the command line
### e.g.
### cmsRun testElectronSequence_cfg.py print maxEvents=10
### -------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

_doMC = False;
_isHI = False;
_isPA = True;
_doJets = False;

options = VarParsing.VarParsing('analysis')
options.outputFile = "WwithPF.root"
#options.inputFiles = "/store/hidata/PARun2016B/PAEGJet1/AOD/PromptReco-v1/000/285/244/00000/FCBC1BFC-9EAC-E611-AC24-02163E01356A.root"
options.inputFiles = "/store/hidata/PARun2016B/PADoubleMuon/AOD/PromptReco-v1/000/285/244/00000/2E9A949F-A5AC-E611-9368-02163E0142A2.root"
options.maxEvents = -1

options.parseArguments()

##################################################################################
process = cms.Process("PFRECO")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("CommonTools.UtilAlgos.TFileService_cfi")

process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#                                        ignoreTotal=cms.untracked.int32(0),
#                                        oncePerEventMode = cms.untracked.bool(False)
#                                        )
#process.Timing = cms.Service("Timing")

### Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.2013collisionEventSelection_cff")

### pile up rejection


### Centrality for pPb
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("pACentrality")
process.centralityBin.centralityVariable = cms.string("HFtowersPlusTrunc")
process.centralityBin.nonDefaultGlauberModel = cms.string("Epos")
process.EventAna_step = cms.Path( process.centralityBin )

### Global tags for conditions data: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if _doMC:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_PIon', '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v15', '')
  process.GlobalTag.toGet = cms.VPSet(
      cms.PSet(
        record = cms.string("HeavyIonRcd"),
        tag =
        cms.string("CentralityTable_HFtowersPlusTrunc200_EPOS5TeV_v80x01_mc"),
        connect =
        cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowersPlusTruncEpos")
        ),
      cms.PSet(
        record = cms.string('L1TUtmTriggerMenuRcd'),
        tag = cms.string("L1Menu_HeavyIons2016_v2_m2_xml"),
        connect =
        cms.string('frontier://FrontierProd/CMS_CONDITIONS')
        ),
      cms.PSet(
        record = cms.string('L1TGlobalPrescalesVetosRcd'),
        tag = cms.string("L1TGlobalPrescalesVetos_Stage2v0_hlt"),
        connect =
        cms.string('frontier://FrontierProd/CMS_CONDITIONS')
        )
      )
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

### HLT muon path selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths =  [
  "HLT_PAL2Mu12_v1",
  "HLT_PAL2Mu15_v1",
  "HLT_PAL3Mu12_v1",
  "HLT_PAL3Mu15_v1"
  ]
process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True
process.hltOniaHI.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

### PAT muon production
from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=_doMC, HLT="HLT", Filter=False, useL1Stage2=True)

### Temporal fix for the PAT Trigger prescale warnings.
process.patTriggerFull.l1GtReadoutRecordInputTag = cms.InputTag("gtDigis","","RECO")
process.patTriggerFull.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
process.patTriggerFull.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
if _doMC:
  process.patTriggerFull.getPrescales      = cms.untracked.bool(False)
else:
  process.patTriggerFull.getPrescales      = cms.untracked.bool(True)
###

process.patMuons = cms.Sequence(
    process.PAcollisionEventSelection *
#    process.pileupVertexFilterCutGplus * 
#    process.pACentrality_step *
    process.patMuonSequence
    )

if _doMC:
  process.patMuons.remove(process.hltOniaHI)


### PAT muon selection
process.goodPatMuons = cms.EDFilter("PATMuonSelector",
                                     src = cms.InputTag("patMuonsWithTrigger"),
                                     cut = cms.string("pt>15."),
                                     filter = cms.bool(True)
                                   )
process.source       = cms.Source("PoolSource",
                                   fileNames = cms.untracked.vstring( options.inputFiles )
                                 )
process.TFileService = cms.Service("TFileService", 
                                    fileName = cms.string( options.outputFile )
                                  )
process.maxEvents    = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options      = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

### PFCandAnalyzer
import PFMETMuon.PFMETMuonAnalyzer.pfcandAnalyzer_pp_cfi
process.load("PFMETMuon.PFMETMuonAnalyzer.pfcandAnalyzer_pA_cfi")
process.pfcandAnalyzer.doMC = cms.untracked.bool(_doMC)
process.pfcandAnalyzer.isHI = cms.untracked.bool(_isHI)
process.pfcandAnalyzer.isPA = cms.untracked.bool(_isPA)
process.pfcandAnalyzer.doJets = cms.untracked.bool(_doJets)


process.ntuples = cms.Path(
    process.patMuons*
    process.goodPatMuons*
    process.pfcandAnalyzer
    #process.ak5PFJetAnalyzer*
    #process.akPu3PFJetAnalyzer
    #process.pftowerAna
    #process.anaMET
    )

process.schedule = cms.Schedule(process.ntuples)
