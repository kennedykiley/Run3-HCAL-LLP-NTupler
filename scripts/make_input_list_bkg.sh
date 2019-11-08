#!/bin/bash
version=/privateProduction/DR/step1/RunIISummer16/GENSIM/WminusH_HToSSTobbbb_CSCDecayFilter_ms55_pl10000/v3/
root_dir=/mnt/hadoop/store/group/phys_exotica/${version}
#list_dir=$CMSSW_BASE/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p7/MC_Summer16/v3/sixie/
list_dir=${CMSSW_BASE}/src/cms_lpc_llp/llp_ntupler/lists/${version}
echo $list_dir
mkdir -p $list_dir
for sample in \
WminusH_HToSSTobbbb_ms55_pl10000_ev150000
do
        echo "${list_dir}${sample}.txt"
        rm -f ${list_dir}${sample}.txt
	mkdir -p ${list_dir}
        find ${root_dir} -name "*.root" >> ${list_dir}${sample}.txt
        sed -i '/failed/d' ${list_dir}${sample}.txt
        echo "input list created for $sample"

done
