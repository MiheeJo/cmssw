from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'ExpressPhysicsPA_test1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runPFOnpPb_cfg.py'
config.JobType.pyCfgParams = ['noprint']
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.inputDataset = '/ExpressPhysicsPA/PARun2016C-Express-v1/FEVT'
#config.Data.inputDataset = '/PASingleMuon/PARun2016C-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 200
#config.Data.totalUnits = 1
config.Data.runRange = '285090,285368,285480,285505,285517,285530,285537,285538,285539,285549'
#config.Data.outLFNDirBase = '/store/user/miheejo/%s' % (config.General.requestName)
config.Data.outLFNDirBase = '/store/user/miheejo'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.publishDBS = 'phys03'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_FR_CCIN2P3'
#config.Site.storageSite = 'T2_FR_GRIF_LLR'
