from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os
config = config()

config.General.workArea = '/nfs/dust/cms/user/jschindl/crab_recHit_study'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ConfFile_cfg.py'
#config.JobType.inputFiles = ['data']

#config.General.requestName = 'SingleMuon_Run2017A-PromptReco-v2' #no LS in golden json
#config.Data.inputDataset = '/SingleMuon/Run2017A-PromptReco-v2/MINIAOD' #no LS in golden json
#config.General.requestName = 'SingleMuon_Run2017A-PromptReco-v3' #no LS in golden json
#config.Data.inputDataset = '/SingleMuon/Run2017A-PromptReco-v3/MINIAOD' #no LS in golden json
#config.General.requestName = 'SingleMuon_Run2017B-PromptReco-v1'
#config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v1/MINIAOD'
#config.General.requestName = 'SingleMuon_Run2017B-PromptReco-v2'
#config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v2/MINIAOD'
config.General.requestName = 'DT_RecHit_study_WplusH_HToSSTobbbb_ms55_pl10000_v4'
config.Data.inputDataset = "/WplusH_HToSSTobbbb_ms55_pl10000_ev150000/sixie-crab_CMSSW_9_4_12_PrivateProduction_Fall17_DR_step1_WplusH_HToSSTobbbb_ms55_pl10000_v2_DR_CaltechT2-fb4de91b3672ad6012188656f7233fe2/USER"

config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'LumiBased'#if Data
config.Data.splitting = 'FileBased'#if MC

config.JobType.maxMemoryMB = 4000 #more memory
config.JobType.numCores = 8
config.JobType.maxJobRuntimeMin = 600

#if DATA:
#config.Data.lumiMask = 'data/JSON/Cert_294927-305364_13TeV_PromptReco_Collisions17_JSON.txt'
#'Cert_294927-302654_13TeV_PromptReco_Collisions17_JSON.txt'#17.85 fbinv
#'data/JSON/Cert_294927-302343_13TeV_PromptReco_Collisions17_JSON.txt'#'data/JSON/Cert_294927-301567_13TeV_PromptReco_Collisions17_JSON.txt'#golden json
#config.Data.lumiMask = 'data/JSON/json_DCSONLY.txt' #DCS only

#older
#config.Data.lumiMask = 'data/JSON/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt'
#'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'#proper way to link a json file from DQM certification
config.Data.unitsPerJob = 15

#config.Data.totalUnits = 10

#If storing in T2
config.Data.outLFNDirBase = '/store/user/jschindl/recstudy_ntuple'
config.Data.publication = False
#config.Data.outputDatasetTag = 'MET_trigger'

config.User.voGroup='dcms'
#If willing to publish on dbs
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'

config.Site.storageSite = 'T2_DE_DESY'#modify with a T2 where you have writing access
config.Site.blacklist   = ['T2_FR_IPHC']
