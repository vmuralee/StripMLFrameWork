# Strip Cluster MLFrameWork
The aim of this frame-work to differentiate cluster formed from the noise over the actual hits using the classification neural network. The classification model futher used to reduce the size of the strip cluster. The new strip cluster collection validated by reconstructed tracking performes. The framework details as follows. 

## Setup
   ```
    cmsrel CMSSW_14_1_0_pre3
    cd CMSSW_14_1_0_pre3/src/
    git clone https://github.com/vmuralee/StripMLFrameWork.git
    cd StripClusterMLStudies
  ```

The details for each package shown below,

  - #### StripMLAnalyzer
`cd StripMLAnalyzer/`
    It is a `EDAnalyzer` to prepare the input dataset for neural netowork training. The strip cluster matched to the recohits in the strip deteector are the signal and rest are the background. The input to the `EDAnalyzer` can obtain by rerun the HLT.
    
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
  and update the `process.hltSiStripClusters2ApproxClusters` to `process.hltSiStripClusters2ApproxClustersold`. The output file after running  `cmsRun  python/prehlt_mc.py`, used for the making the flat ntuples by `cmsRun python/rawprimetupleproducer_cfg.py`. 
  
  - #### MLP
`cd MLP/`
    The classifier network is Architecture discussed in this folder. The trained model saved at `DnnApproxCluster/saved_model`


  - #### DnnApproxCluster
    `cd DnnApproxCluster/`
    
    It is a `EDProducer` which store the approximated cluster after applying `NN score` in each event. The configuration `test/DnnApproxCluster_cfg.py` has the parameter called `NN_threshold` which can vary the cut on the DNN score. It adds new reduced approximate cluster `*_hltSiStripClusters2ApproxClusters_*_nnHLT` and drop the `*_hltSiStripClusters2ApproxClustersold_*_HLTX` collection. 

  - #### Track reconstruction

```
cmsDriver.py step_reco --conditions 140X_mcRun3_2024_realistic_v21 -s RAW2DIGI,L1Reco,RECO  --datatier RECO --eventcontent RECO --mc --process reRECO --scenario pp --era Run3_pp_on_PbPb_approxSiStripClusters --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3 --hltProcess HLTX -n 100 --filein file:/eos/home-v/vmuralee/PREanalysis/clusterMLstudies/SinglePionApproxDnn.root --repacked --no_exec
```

Due to avoid the error I have command out following lines,
```
 from Configuration.Applications.ConfigBuilder import MassReplaceInputTag                                                                                                            MassReplaceInputTag(process, new="rawDataMapperByLabel", old="rawDataCollector")
```
The output file contains the reconstructed tracks from our new sistrip collection. 

 - #### TrackEfficiencyAnalyzer

   
     
 
    
