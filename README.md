# DT_Analyzer


cmsrel CMSSW_9_4_7
cd CMSSW_9_4_7/src
git clone git@github.com:JoergSchindler/DT_Analyzer.git Analyzer/Demo
scram b
cmsenv
# Run the ntuples
cmsRun python/ConfFile_cfg.py
