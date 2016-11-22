import FWCore.ParameterSet.Config as cms

pfcandAnalyzer = cms.EDAnalyzer('PFMETMuonAnalyzer',
                                pfCandidateLabel = cms.InputTag("particleFlow"),
                                jetLabel = cms.InputTag("ak5patJets"),
                                genLabel = cms.InputTag("genParticles"),
                                CentralitySrc = cms.InputTag("pACentrality"),
                                CentralityBinSrc = cms.InputTag("centralityBin","HFtowersPlusTrunc"),
                                pfMETLabel = cms.InputTag("pfMet"),
                                trackLabel = cms.InputTag("generalTracks"),
                                trackQuality = cms.string("highPurity"),
                                muonLabel = cms.InputTag("goodPatMuons"),
                                vtxLabel = cms.InputTag("offlinePrimaryVertices"),

                                pfPtMin = cms.double(0.),
                                genPtMin = cms.double(0.5),
                                jetPtMin = cms.double(20.0),                                
                                etaBins = cms.int32(15),
                                fourierOrder = cms.int32(5),

                                verbosity = cms.untracked.int32(0),
                                skipCharged = cms.untracked.bool(False),
                                doVS = cms.untracked.bool(False),
                                bkg = cms.InputTag("voronoiBackgroundPF"),
                                doUEraw_ = cms.untracked.bool(False),

                                doJets = cms.untracked.bool(True),
                                doMC = cms.untracked.bool(False),
                                isHI = cms.untracked.bool(False),
                                isPA = cms.untracked.bool(True),

                                #Trigger information
                                triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
                                triggerPathNames = cms.vstring("HLT_PAL2Mu12_v1",
                                                               "HLT_PAL2Mu15_v1",
                                                               "HLT_PAL3Mu3_v1",
                                                               "HLT_PAL3Mu5_v1",
                                                               "HLT_PAL3Mu7_v1",
                                                               "HLT_PAL3Mu12_v1",
                                                               "HLT_PAL3Mu15_v1"),
                                triggerFilterNames = cms.vstring("hltL2fL1sSingleMu7BptxANDL1f0L2Filtered12",
                                                                 "hltL2fL1sSingleMu7BptxANDL1f0L2Filtered15",
                                                                 "hltL3fL1sSingleMu3BptxANDL1f0L2f0L3Filtered3",
                                                                 "hltL3fL1sSingleMu5BptxANDL1f0L2f0L3Filtered5",
                                                                 "hltL3fL1sSingleMu5BptxANDL1f0L2f0L3Filtered7",
                                                                 "hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12",
                                                                 "hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered15")
)
