import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIOnia")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = "Jpsi_Histos.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = 'file:/tmp/camelia/onia2MuMuPAT_740_AOD.root'
options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# tag for running on 2011 data in 7xy
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')

# centrality part
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality             = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable     = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")

process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
      tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v740x01_mc"),
      connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
      label = cms.untracked.string("HFtowersHydjetDrum5")
   ),
])


# event plane
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
'''
process.GlobalTag.toGet.extend([
        cms.PSet(record = cms.string("HeavyIonRPRcd"),
                 tag = cms.string('/afs/cern.ch/user/m/mironov/scratch0/CMSSW_7_4_0/src/HiAnalysis/HiOnia/test/HeavyIonRPRcd_Hydjet_74x_v02_mc.db'),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_PAT_000")
                 )
        ])
'''


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    )
)

process.hltDblMuOpen = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_HIL1DoubleMu0_HighQ_v*"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)

process.hionia = cms.EDAnalyzer('HiOniaAnalyzer',
                                srcMuon             = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig       = cms.InputTag("patMuonsWithoutTrigger"),
                                src                 = cms.InputTag("onia2MuMuPatGlbGlb"),
                                srcTracks           = cms.InputTag("hiGeneralTracks"),
                                genParticles        = cms.InputTag("hiGenParticles"),
                                primaryVertexTag    = cms.InputTag("hiSelectedVertex"),
                                triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
                             
                                CentralitySrc    = cms.InputTag("hiCentrality"),
                                CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
                                EvtPlane         = cms.InputTag("hiEvtPlane","recoLevel"),
                                EvtPlaneFlat     = cms.InputTag("hiEvtPlaneFlat",""),   

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                                useRapidity = cms.bool(True),
                                
                                #--
                                maxAbsZ = cms.double(24.0),
                                
                                pTBinRanges      = cms.vdouble(0.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 40.0),
                                etaBinRanges     = cms.vdouble(0.0, 2.5),
                                centralityRanges = cms.vdouble(20,40,100),

                                onlyTheBest     = cms.bool(False),		
                                applyCuts       = cms.bool(True),
                                storeEfficiency = cms.bool(False),
                      
                                removeSignalEvents = cms.untracked.bool(False),
                                removeTrueMuons    = cms.untracked.bool(False),
                                storeSameSign      = cms.untracked.bool(True),
                                muonLessPV         = cms.bool(False), # this has to be FALSE for HI reco algo
                                
                                #-- Gen Details
                                oniaPDG = cms.int32(443),
                                isHI = cms.untracked.bool(True),
                                isPA = cms.untracked.bool(False),
                                isMC = cms.untracked.bool(False),
                                isPromptMC = cms.untracked.bool(False),
                                useEvtPlane = cms.untracked.bool(False),
                                runVersionChange = cms.untracked.uint32(182133),

                                #-- Histogram configuration
                                combineCategories = cms.bool(False),
                                fillRooDataSet    = cms.bool(False),
                                fillTree          = cms.bool(True),
                                fillHistos        = cms.bool(False),
                                minimumFlag       = cms.bool(True),
                                fillSingleMuons   = cms.bool(True),
                                fillRecoTracks    = cms.bool(False),
                                histFileName      = cms.string(options.outputFile),		
                                dataSetName       = cms.string(options.secondaryOutputFile),
                                
                                #--
                                # NumberOfTriggers = cms.uint32(8),
                                dblTriggerPathNames = cms.vstring("HLT_HIL1DoubleMu0_HighQ_v*",
                                                                  "HLT_HIL2DoubleMu3_v*",
                                                                  "HLT_HIL3DoubleMuOpen_v*",
                                                                  "HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v*"),
                                dblTriggerFilterNames = cms.vstring("hltHIDoubleMuLevel1PathL1HighQFiltered",
                                                                    "hltHIL2DoubleMu3L2Filtered",
                                                                    "hltHIDimuonL3FilteredOpen",
                                                                    "hltHIDimuonL3FilteredMg2OSnoCowboy"),
                                sglTriggerPathNames = cms.vstring("HLT_HIL2Mu3_NHitQ_v*",
                                                                  "HLT_HIL2Mu7_v*",
                                                                  "HLT_HIL2Mu15_v*",
                                                                  "HLT_HIL3Mu3_v*"),
                                sglTriggerFilterNames = cms.vstring("hltHIL2Mu3NHitL2Filtered",
                                                                    "hltHIL2Mu7L2Filtered",
                                                                    "hltHIL2Mu15L2Filtered",
                                                                    "hltHISingleMu3L3Filtered")


                                )


#process.p = cms.Path(process.hltDblMuOpen*process.hionia)
process.p = cms.Path(process.centralityBin*process.hionia)
