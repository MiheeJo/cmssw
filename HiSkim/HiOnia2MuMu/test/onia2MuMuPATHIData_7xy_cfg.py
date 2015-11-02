
# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples


import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# set up process
process = cms.Process("Onia2MuMuPAT")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.inputFiles = 'file:/afs/cern.ch/user/t/tuos/work/public/reco2AOD/round2April28/DIMUON/step2_RAW2DIGI_L1Reco_DIMUONskim_AOD.root'#step2_reRECO_740_100_1_lRV.root'
options.outputFile = '/tmp/camelia/onia2MuMuPAT_740_AOD.root'

options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

# event plane stuff
#process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
# tag for running on 2011 data in 7xy
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')


# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
    )

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths = ["HLT_HIL1DoubleMu0_HighQ_v*",
                              "HLT_HIL2Mu3_NHitQ_v*",
                              "HLT_HIL2Mu7_v*","HLT_HIL2Mu15_v*",
                              "HLT_HIL2DoubleMu0_NHitQ_v*",
                              "HLT_HIL2DoubleMu3_v*",
                              "HLT_HIL3Mu3_v*",
                              "HLT_HIL3DoubleMuOpen_v*","HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v*"
                              ]
process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True
process.hltOniaHI.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=False)
process.onia2MuMuPatGlbGlb.primaryVertexTag         = cms.InputTag("hiSelectedVertex")
process.onia2MuMuPatGlbGlb.dimuonSelection          = cms.string("mass > 0")
process.onia2MuMuPatGlbGlb.addMuonlessPrimaryVertex = False
process.onia2MuMuPatGlbGlb.resolvePileUpAmbiguity   = False

process.source.fileNames      = cms.untracked.vstring(options.inputFiles)
process.maxEvents             = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.outOnia2MuMu.fileName = cms.untracked.string( options.outputFile )


process.patMuonSequence = cms.Sequence(process.patMuonsWithTriggerSequence)
process.e               = cms.EndPath(process.outOnia2MuMu+process.outTnP)

process.schedule = cms.Schedule(process.Onia2MuMuPAT,
                               process.TagAndProbeTrig,
                               process.TagAndProbeSta, 
                               process.TagAndProbeMuID,
process.e)

from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process)
