from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'muonAnalyzerTest'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonAnalyzer_pbpb.py'
config.JobType.pyCfgParams = ['noprint']
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.inputDataset = '/Pythia8_JpsiMM_ptJpsi_12_15_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 10
config.Data.outLFNDirBase = '/store/user/miheejo/%s/Pythia8_JpsiMM_ptJpsi1215' % (config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.publishDBS = 'phys03'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_FR_CCIN2P3'
#config.Site.storageSite = 'T2_FR_GRIF_LLR'
