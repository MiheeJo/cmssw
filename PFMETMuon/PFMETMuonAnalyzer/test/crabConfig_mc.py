from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Pyquen_WToMuNu_test1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runPFOnpPb_mc_cfg.py'
config.JobType.pyCfgParams = ['noprint']
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.inputDataset = '/MinBias/echapon-Pyquen_WToMuNu_step3_20161109_2-adec03854c5299734a1a1f4fec70dfbb/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/miheejo'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.publishDBS = 'phys03'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_FR_CCIN2P3'
#config.Site.storageSite = 'T2_FR_GRIF_LLR'
