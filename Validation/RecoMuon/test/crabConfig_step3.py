from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'muonTrackValidator_pp'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonValidation.py'
config.JobType.outputFiles = ['validation.root']
config.JobType.pyCfgParams = ['noprint']
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.inputDataset = '/JpsiMM_5p02TeV_TuneCUETP8M1/miheejo-RAW2DIGI_RECO_LLRv2-0564587735dfa98972125c928a8975ef/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 10
config.Data.outLFNDirBase = '/store/user/miheejo/%s/Pythia8_Jpsi/' % (config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
#config.Data.publishDBS = 'global'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
