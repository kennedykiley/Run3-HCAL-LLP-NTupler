#!/bin/bash
root_dir=/mnt/hadoop/store/group/phys_exotica/privateProduction/DR/step1/RunIIFall17/GENSIM/WplusH_HToSSTobbbb_ms55_pl10000/v2/WplusH_HToSSTobbbb_ms55_pl10000_ev150000/crab_CMSSW_9_4_12_PrivateProduction_Fall17_DR_step1_WplusH_HToSSTobbbb_ms55_pl10000_v2_DR_CaltechT2/191014_015936/0000/
#list_dir=$CMSSW_BASE/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p7/MC_Summer16/v3/sixie/
list_dir=$CMSSW_BASE/src/cms_lpc_llp/llp_ntupler/lists/
echo $list_dir
mkdir -p $list_dir
for sample in \
WplusH_HToSSTobbbb_ms55_pl10000_ev150000
do
        echo "${list_dir}${sample}.txt"
        rm -f ${list_dir}${sample}.txt
        find ${root_dir} -name "*.root" >> ${list_dir}${sample}.txt
        sed -i '/failed/d' ${list_dir}${sample}.txt
        echo "input list created for $sample"

done
