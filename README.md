# DT_Analyzer

```bash
# Done once to setup environment
cmsrel CMSSW_9_4_7
cd CMSSW_9_4_7/src
git clone git@github.com:JoergSchindler/DT_Analyzer.git Demo/DemoAnalyzer
scram b
cmsenv
# Run the ntuples
cmsRun python/ConfFile_cfg.py
```
