from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'muonTrackValidator'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonValidation.py'
config.JobType.pyCfgParams = ['noprint']
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.inputDataset =''
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/miheejo/%s/Pythia8_JpsiMM_ptJpsi0003/' % (config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
#config.Data.publishDBS = 'global'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
