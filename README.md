# Strip Cluster MLFrameWork
The aim of this frame-work to differentiate cluster formed from the noise over the actual hits using the classification neural network. The classification model futher used to reduce the size of the strip cluster. The new strip cluster collection validated by reconstructed tracking performes. The framework details as follows. 

## Setup
   ```
    cmsrel CMSSW_14_1_0_pre3
    cd CMSSW_14_1_0_pre3/src/
    git clone https://github.com/vmuralee/StripMLFrameWork.git
    cd StripClusterMLStudies
    scram b -j 32
  ```

The details for each package shown below,

First, prepare a config file to produce `SiStripApproxClusterCollection`,
   
  ```
  hltGetConfiguration /users/vmuralee/PREmenu/V11
  --globaltag 140X_mcRun3_2024_realistic_v21 --mc --unprescale
  --max-events 100 --eras Run3
  --input /store/relval/CMSSW_14_1_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-RECO/140X_mcRun3_2024_realistic_v21_STD_Recycled_2024_noPU-v1/2580000/a65a9d3e-c48b-4316-9947-64c1906fc136.root > prehlt_mc.py
  ```

  
  The secondaryFileNames are ,
  ```
  '/store/relval/CMSSW_14_1_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v21_STD_Recycled_2024_noPU-v1/2580000/734b9a99-a98a-4c81-8e3e-21121f074b09.root',
  '/store/relval/CMSSW_14_1_0/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/140X_mcRun3_2024_realistic_v21_STD_Recycled_2024_noPU-v1/2580000/b0834551-1e2c-4d39-b51d-d8e554fa3f43.root'
  ```
 A modified version of config file added in `DnnApproxCluster/test/hltRelVal_cfg.py` 

  - #### StripMLAnalyzer
`cd StripMLAnalyzer/`
    It is a `EDAnalyzer` to prepare the input dataset for neural netowork training. The strip cluster matched to the recohits in the strip deteector are the signal and rest are the background. The input to the `EDAnalyzer` can obtain by rerun the HLT.
 
  
  - #### MLP
`cd MLP/`
    The classifier network is Architecture discussed in this folder. The trained model saved at `DnnApproxCluster/saved_model`


  - #### DnnApproxCluster
    `cd DnnApproxCluster/`
    
    It is a `EDProducer` which store the approximated cluster after applying `NN score` in each event. The configuration `test/DnnApproxCluster_cfg.py` has the parameter called `NN_threshold` which can vary the cut on the DNN score. It adds new reduced approximate cluster `*_hltSiStripClusters2ApproxClustersDNN{NN_threshold}_*_nnHLT`. While running the code add the collection of `NN_threshold` ranges {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9}.
  - #### Track reconstruction

```
cmsDriver.py step_reco --conditions 140X_mcRun3_2024_realistic_v21 -s RAW2DIGI,L1Reco,RECO  --datatier RECO --eventcontent RECO --mc --process reRECO --scenario pp --era Run3_pp_on_PbPb_approxSiStripClusters --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3 --hltProcess HLTX -n 100 --filein file:/eos/home-v/vmuralee/PREanalysis/clusterMLstudies/SinglePionApproxDnn.root --repacked --no_exec
```

A modified config file after `edmDumpConfig` is given in `DnnApproxCluster/test/reco_config.py`. Which rerun the **RECO** for the tracks. 

 - #### TrackEfficiencyAnalyzer
   `cd TrackEfficiencyAnalyzer/` is analyzer to perform the track efficiency. 

   
     
 
    
